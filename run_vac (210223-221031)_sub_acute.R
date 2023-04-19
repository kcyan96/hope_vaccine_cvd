setwd("/mnt/Others/Dataset/LONGCOVID_CC")
library(data.table)
library(survey)
library(epitools)


# Subgroup
for (sub in(c("age.l65","age.ge65","sex.male","sex.female", "cci.l4", "cci.ge4", "sev1", "sev0", "bnt", "sino"))){
  
  SA <- sub
  cohort_full <- readRDS("cohort_acute_jiayi.RDS")
  cohort <- cohort_full[!is.na(group) & infected==1]
  cohort$group <- cohort$group_dose
  
  # subgroups
  if(SA=="age.l65") cohort <- cohort[Age<65]
  if(SA=="age.ge65") cohort <- cohort[Age>=65]
  if(SA=="sex.male") cohort <- cohort[sex=="M"]
  if(SA=="sex.female") cohort <- cohort[sex=="F"]
  if(SA=="cci.l4") cohort <- cohort[score.cci.age<4]
  if(SA=="cci.ge4") cohort <- cohort[score.cci.age>=4]
  if(SA=="sev1") cohort <- cohort[severity==1]
  if(SA=="sev0") cohort <- cohort[severity==0]
  if(SA=="bnt") cohort <- cohort[group=="X"|BNT162b2==1]
  if(SA=="sino") cohort <- cohort[group=="X"|CoronaVac==1]
  
  # SELECT comparison
  COMP <- "INF_VAC_dose_May"
  PHASE <- "_acute"
  
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
  cohort[, acute := 1]
  W.glm <- weightit(as.formula(PS_formula("group",4)), data=cohort, estimand="ATE", method="ps")
  PS_assess(W.glm)
  cohort[, wts_acute := W.glm$weights]
  
  
  saveRDS(cohort, paste0("COMP_", SA, "_", COMP,  PHASE, "_2.cohort_iptw_jiayi.RDS"))
  gc()
}


setwd("/mnt/Others/Dataset/LONGCOVID_CC")

for (sub in(c("age.l65","age.ge65","sex.male","sex.female", "cci.l4", "cci.ge4", "sev1", "sev0", "bnt", "sino"))){
  
  SA <- sub
  
  PHASE <- "_acute"
  
  # 3. Censoring ----
  COMP <- "INF_VAC_dose_May"
  cohort<-readRDS(paste0("COMP_", SA, "_", COMP,  PHASE, "_2.cohort_iptw_jiayi.RDS"))
  
  OUTCOMES <- c("cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death")  
  cohort[, `:=`(hx.cv_death = 0, hx.death = 0)]
  cohort[, postacute.index.date := index.date+28]
  
  # Acute phase (within 28 days)
  for(oc in OUTCOMES) {
    cohort[, c(paste0("acute.hx.",oc)) := as.numeric(get(paste0("hx.",oc)))]
    cohort[, c(paste0("acute.censor.date.", oc)) := pmin(get(paste0("outcome.",oc,".date")), death_date_ymd, postacute.index.date-1, na.rm=T)]
    cohort[, c(paste0("acute.time.to.censor.", oc)) := as.numeric(get(paste0("acute.censor.date.", oc)) - index.date) + 1]
    cohort[, c(paste0("acute.outcome.", oc)) := as.numeric(!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("acute.censor.date.", oc))==get(paste0("outcome.",oc,".date")))]
  }
  rm(oc)
  
  # 4. Survival analysis ----
  
  library(survival)
  
  run_fit <- function(data, phase, outcomes) {
    
    res <- list()
    for(oc in outcomes) {
      cat(phase, oc, "\n")
      dt <- data[get(phase)==1 & get(paste0(phase,".hx.",oc))==0]
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
      dt <- data[get(phase)==1 & get(paste0(phase,".hx.",oc))==0]
      events_fu <-dt[, lapply(.SD, function(x, w) sum(x*w), w=get(paste0("wts_", phase))), .SDcols = c(paste0(phase,".outcome.",oc),paste0(phase,".time.to.censor.",oc)),keyby=group_dose]
      res_IR[[oc]] <-setDT(pois.exact(events_fu[,2], events_fu[,3]))[,`:=`(group=events_fu$group_dose,N=x)][,c('group','N','pt','rate','lower','upper','conf.level')]
    }
    return(res_IR)
  }
  
  
  res_t2 <- sapply(c("acute"), function(p) run_fit(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)
  res_t2_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)
  
  
  sink(paste0("/mnt/Others/Dataset/LONGCOVID_CC/output/", SA, PHASE, "_jiayi_graph_cvds.txt"), append=F, split=T)
  print(SA)
  for(oc in OUTCOMES) {
  cat(oc, "\n")
  print(cbind(setnames(cbind(as.data.frame(res_t2_IR$acute[[oc]])[,1:6],res_t2$acute[[oc]]$events[,3]),c('group','wEvents','pt','rate','lower','upper','Events')),
              rbind("-",exp(cbind(coef(res_t2$acute[[oc]]$adj),confint(res_t2$acute[[oc]]$adj)))),
              setnames(rbind("-", as.data.frame(coef(res_t2$acute[[oc]]$adj))),"estimate"),
              setnames(rbind("-",as.data.frame(coef(summary(res_t2$acute[[oc]]$adj))[,3])),"se")))  
  }
  sink() 
}

