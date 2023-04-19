setwd("/mnt/Others/Dataset/LONGCOVID_CC/")
library(data.table)
library(readxl)
library(survey)

# 0. Helpers ----

table1 <- function(cohort, wts_col=NULL, strata="group") {
  cohort_ <- copy(cohort)
  cohort_[, (grep("^(dx|px|ip|rx|hx|mm)",names(cohort_),value=T)) := lapply(.SD, as.logical), .SDcols=grep("^(dx|px|ip|rx|hx|mm)",names(cohort_),value=T)]
  
  if(!is.null(wts_col)) {
    cohort_ <- svydesign(ids = ~ 1, data = cohort_, weights = as.formula(paste0("~ ", wts_col)))
    tblfun <- svyCreateTableOne
  }
  else {
    tblfun <- CreateTableOne
  }
  
  tb1_def <- setDT(read_excel("codes_new.xlsx", sheet="table1"))
  t1 <- tblfun(tb1_def[!is.na(Name), Name], c(strata), data=cohort_, test=F)
  print(t1, smd=T, test=F)
  return(t1)
}

as.data.table.TableOne <- function(t) {
  tb1_def <- setDT(read_excel("codes_new.xlsx", sheet="table1"))[!is.na(Name), .(Name, Description)]
  tb1_def <- rbind(list(Name="n", Description="n"), tb1_def)
  t <- as.data.frame(print(t, test=F, smd=T, dropEqual=T, noSpaces=T))
  varlabels <- rownames(t)
  t$Name = sub("^([a-zA-Z0-9._]+).*$", "\\1", varlabels)
  t <- merge(as.data.table(t), tb1_def, by="Name", all.x=T, sort=F)
  t$Description = paste0(t$Description, sub("^([a-zA-Z0-9._]+)", "", varlabels))
  t[!is.na(Description), Name:=Description][, Description:=NULL]
  return(t)
}



# 1. Load cohort ----
cohort <- readRDS("cohort_acute_jiayi.RDS")
cohort[, postacute:=as.numeric(is.na(death_date_ymd) | death_date_ymd >= (index.date+90))]
cohort_postacute90d <- copy(cohort[postacute==1])
cohort_postacute90d[, index.date := (index.date+90)]
saveRDS(cohort_postacute90d, "cohort_postacute90d_jiayi.RDS")
rm(cohort, cohort_postacute90d)

## Weighting ---------------------------------------------
cohort_full <- readRDS("cohort_postacute90d_jiayi.RDS")
# SELECT comparison
COMP <- "INF_VAC_dose_May"
cohort <- cohort_full[!is.na(group) & infected==1]
cohort$group <- cohort$group_dose
PHASE <- "_postacute90d"

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
cohort[, wts_postacute90d := W.glm$weights]

saveRDS(cohort, paste0("COMP_", COMP, PHASE, "_2.cohort_iptw_jiayi.RDS"))

rm(cohort, cohort_full, W.glm)


library(epitools)
PHASE <- "_postacute90d"


# 3. Censoring ---- Table 2
COMP <- "INF_VAC_dose_May"
#cohort<-readRDS(paste0("COMP_", COMP, PHASE, "_2.cohort_iptw_jiayi.RDS"))

OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death", "death")
cohort[, `:=`(hx.cv_death = 0, hx.death = 0)]
cohort[, postacute90d.index.date := index.date]

cohort <- cohort[postacute90d.index.date <= "2023-01-23"]

# for(oc in OUTCOMES) {
#   cohort[,paste0("outcome.",oc,".date"):=(as.character(get(paste0("outcome.",oc,".date"))))]
#   #cohort[,paste0("outcome.",oc,".date"):=(as.numeric(get(paste0("outcome.",oc,".date"))))]
#   cohort[,paste0("outcome.",oc,".date"):=(as.Date(get(paste0("outcome.",oc,".date")),origin='1970-01-01'))]
#   
# }

# Postacute phase (after 90 days)
for(oc in OUTCOMES) {
  cohort[, c(paste0("postacute90d.hx.",oc)) := as.numeric(get(paste0("hx.",oc)) | (!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("outcome.",oc,".date")) < postacute90d.index.date))]
  cohort[, c(paste0("postacute90d.censor.date.", oc)) := pmin(get(paste0("outcome.",oc,".date")), death_date_ymd, as.Date("2023-01-23"), na.rm=T)]
  cohort[, c(paste0("postacute90d.time.to.censor.", oc)) := as.numeric(get(paste0("postacute90d.censor.date.", oc)) - postacute90d.index.date) + 1]
  cohort[, c(paste0("postacute90d.outcome.", oc)) := as.numeric(!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("postacute90d.censor.date.", oc))==get(paste0("outcome.",oc,".date")) & get(paste0("outcome.",oc,".date"))>=postacute90d.index.date)]
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

res_t2 <- sapply(c("postacute90d"), function(p) run_fit(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)
res_t2_IR <- sapply(c("postacute90d"), function(p) run_fit_IR(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)

sink(paste0("output/T2_postacute90d_log_jiayi.txt"), append=F, split=T)
for(oc in OUTCOMES) {
  cat(oc, "\n")
  print(cbind(setnames(cbind(as.data.frame(res_t2_IR$postacute[[oc]])[,1:6],res_t2$postacute[[oc]]$events[,3]),c('group','wEvents','pt','rate','lower','upper','Events')),
              rbind("-",exp(cbind(coef(res_t2$postacute[[oc]]$adj),confint(res_t2$postacute[[oc]]$adj))))))
}
sink()
