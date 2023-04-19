setwd("/mnt/Others/Dataset/LONGCOVID_CC/")
library(data.table)
library(readxl)
library(survey)

PHASE <- "_postacute"

#cohort <- readRDS("cohort_acute_jiayi.RDS")
cohort <- readRDS("cohort_postacute_jiayi.RDS")

# Remove index date - latest vac > 180 days
for(i in 1:4){
  cohort[nchar(vacc_status)==i, time.since.last_dose := as.numeric(index.date - get(paste0("Date of vaccination.", c("1st","2nd","3rd","4th")[i])))]
}
cohort <- cohort[group=="X" | time.since.last_dose <= 180]   # n=2926865
cohort_180 <- copy(cohort)
saveRDS(cohort_180, "cohort_postacute_180_jiayi.RDS")


## Weighting ---------------------------------------------
cohort_full <- readRDS("cohort_postacute_180_jiayi.RDS")
# SELECT comparison
COMP <- "INF_VAC_dose_May"
cohort <- cohort_full[!is.na(group) & infected==1]
cohort$group <- cohort$group_dose
PHASE <- "_postacute"

cohort[, vaccinated := 0]
cohort[group_dose!="X", vaccinated := 1]
cohort_ <-cohort[vaccinated==1]
cohort_[, group := factor(sub("^$","X",vacc_status_new), levels=c("1d","2d","3d"))]

library(tableone)
#tb1 <- list(acute=table1(cohort), vaccinated=table1(cohort_))

# 2. Propensity score ----

library(cobalt)

PS_def <- setDT(read_excel("codes_new.xlsx", sheet="PS"))[!is.na(Name)]
PS_formula <- function(exp,lvl) paste0(exp, " ~ ", paste0(PS_def[Set<=lvl, Name], collapse=" + "))
PS_assess <- function(out) {
  print(summary(out))
  print(bal.tab(out))
  bal.plot(out, which="both")
}

library(WeightIt)
cohort[, postacute := 1]
W.glm <- weightit(as.formula(PS_formula("group",4)), data=cohort, estimand="ATE", method="ps")
PS_assess(W.glm)
cohort[, wts_postacute := W.glm$weights]

cohort[, vaccinated := 0]
cohort[group_dose!="X", vaccinated := 1]
cohort_ <-cohort[vaccinated==1]
cohort_[, group := factor(sub("^$","X",vacc_status_new), levels=c("1d","2d","3d"))]

library(survey)
#tb1.iptw <- list(postacute=table1(cohort, "wts_postacute"), vaccinated=table1(cohort_, "wts_postacute"))

#save(W.glm, tb1, tb1.iptw, file=paste0("COMP_", COMP, PHASE, "_W.glm.RData"))
#saveRDS(cohort, paste0("COMP_", COMP, PHASE, "_2.cohort_iptw_jiayi.RDS"))


# 3. Censoring ----
COMP <- "INF_VAC_dose_May_180"
#cohort<-readRDS(paste0("LONGCOVID_CC/COMP_", COMP, PHASE, "_2.cohort_iptw.RDS"))

OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv.acs" ,"cv_death", "death")
cohort[, `:=`(hx.cv_death = 0, hx.death = 0)]
cohort[, postacute.index.date := index.date]

for(oc in OUTCOMES) {
  cohort[,paste0("outcome.",oc,".date"):=(as.character(get(paste0("outcome.",oc,".date"))))]
  #cohort[,paste0("outcome.",oc,".date"):=(as.numeric(get(paste0("outcome.",oc,".date"))))]
  cohort[,paste0("outcome.",oc,".date"):=(as.Date(get(paste0("outcome.",oc,".date")),origin='1970-01-01'))]
  
}

# Postacute phase (after 28 days)
for(oc in OUTCOMES) {
  cohort[, c(paste0("postacute.hx.",oc)) := as.numeric(get(paste0("hx.",oc)) | (!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("outcome.",oc,".date")) < postacute.index.date))]
  cohort[, c(paste0("postacute.censor.date.", oc)) := pmin(get(paste0("outcome.",oc,".date")), death_date_ymd, as.Date("2023-01-23"), na.rm=T)]
  cohort[, c(paste0("postacute.time.to.censor.", oc)) := as.numeric(get(paste0("postacute.censor.date.", oc)) - postacute.index.date) + 1]
  cohort[, c(paste0("postacute.outcome.", oc)) := as.numeric(!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("postacute.censor.date.", oc))==get(paste0("outcome.",oc,".date")) & get(paste0("outcome.",oc,".date"))>=postacute.index.date)]
}
rm(oc)



# 4. Survival analysis ----

library(survival)

run_fit <- function(data, phase, outcomes) {
  
  res <- list()
  for(oc in outcomes) {
    cat(phase, oc, "\n")
    dt <- data[get(paste0(phase,".hx.",oc))==0]
    assertthat::assert_that(all(dt[, get(paste0(phase,".time.to.censor.",oc))]>=0))
    res[[oc]] <- list()
    res[[oc]]$events <- dt[, .(.N, Events=sum(get(paste0(phase,".outcome.",oc))), FU=sum(get(paste0(phase,".time.to.censor.",oc)))), keyby=group]
    if(sum(res[[oc]]$events$Events)==0) next
    res[[oc]]$crude <- coxph(as.formula(paste0("Surv(", phase, ".time.to.censor.", oc, ", ", phase, ".outcome.", oc, ") ~ group")), data = dt)
    design <- svydesign(ids = ~ 1, data = dt, weights = as.formula(paste0("~ wts_", phase)))
    res[[oc]]$adj <- svycoxph(as.formula(paste0("Surv(", phase, ".time.to.censor.", oc, ", ", phase, ".outcome.", oc, ") ~ group")), design=design)
  }
  return(res)
  gc()
}

run_fit_IR <- function(data, phase, outcomes) {
  res_IR <- list()
  for(oc in outcomes) {
    cat(phase, oc, "\n")
    dt <- data[get(paste0(phase,".hx.",oc))==0]
    events_fu <-dt[, lapply(.SD, function(x, w) sum(x*w), w=get(paste0("wts_", phase))), .SDcols = c(paste0(phase,".outcome.",oc),paste0(phase,".time.to.censor.",oc)),keyby=group_dose]
    res_IR[[oc]] <-setDT(pois.exact(events_fu[,2], events_fu[,3]))[,`:=`(group=events_fu$group_dose,N=x)][,c('group','N','pt','rate','lower','upper','conf.level')]
  }
  return(res_IR)
}

res_t2 <- sapply(c("postacute"), function(p) run_fit(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)
res_t2_IR <- sapply(c("postacute"), function(p) run_fit_IR(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)

sink(paste0("output/T2_postacute_28d_exl180_jiayi.txt"), append=F, split=T)
print("Table 2")
for(oc in OUTCOMES) {
  cat(oc, "\n")
  print(cbind(setnames(cbind(as.data.frame(res_t2_IR$postacute[[oc]])[,1:6],res_t2$postacute[[oc]]$events[,3]),c('group','wEvents','pt','rate','lower','upper','Events')),
              rbind("-",exp(cbind(coef(res_t2$postacute[[oc]]$adj),confint(res_t2$postacute[[oc]]$adj))))))
}
sink()
