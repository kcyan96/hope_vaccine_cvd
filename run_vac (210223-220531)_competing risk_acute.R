setwd("/mnt/Others/Dataset/LONGCOVID_CC/")
library(data.table)
library(readxl)
library(survey)
library(dplyr)
library(survival)
library(epitools)

PHASE <- "_acute"

# 3. Censoring ----
COMP <- "INF_VAC_dose_May"
cohort<-readRDS(paste0("COMP_", COMP, PHASE, "_2.cohort_iptw_jiayi.RDS"))

#OUTCOMES <- c( "death") 
OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death", "death")
#c("cv.major" ,"cv.stroke" ,"cv.tia" ,"cv.af" ,"cv.afl" ,"cv.pericarditis" ,"cv.myocarditis" ,"cv.cad" ,"cv.acs" ,"cv.mi" ,"cv.stangina" ,"cv.unstangina" ,"cv.hf" ,"cv.cardarrest" ,"cv.cardshock" ,"cv.pe" ,"cv.dvt", "cv_death", "death")
#OUTCOMES <- c("cv.major", "cv.cad", "cv.stroke", "cv.hf", "cv.af", "cv.mi", "cv_death", "death")        # old c("cvd_comp","cad_comp","stroke_comp","hf_comp","dm_neuro","dm_retino","dm_nephro","esrd","cv_death","death")

#cohort[, `:=`(hx.death = 0)]
cohort[, `:=`(hx.cv_death = 0, hx.death = 0)]
cohort[, postacute.index.date := index.date+28]


# Acute phase (within 28 days)
OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death", "death")

for(oc in OUTCOMES) {
  cohort[, c(paste0("acute.hx.",oc)) := as.numeric(get(paste0("hx.",oc)))]
  cohort[, c(paste0("acute.censor.date.", oc)) := pmin(get(paste0("outcome.",oc,".date")), death_date_ymd, postacute.index.date-1, na.rm=T)]
  cohort[, c(paste0("acute.time.to.censor.", oc)) := as.numeric(get(paste0("acute.censor.date.", oc)) - index.date) + 1]
  cohort[, c(paste0("acute.outcome.", oc)) := as.numeric(!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("acute.censor.date.", oc))==get(paste0("outcome.",oc,".date")))]
}
rm(oc)

## Competing risk using Fine-gray model
OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death")

for(oc in OUTCOMES) {
  cohort[acute.outcome.death==1& get(paste0("acute.outcome.", oc))==0, (paste0("acute.outcome.", oc)) := 2]
  cohort[, paste0("acute.outcome.", oc, "1") := as.factor(recode(get(paste0("acute.outcome.", oc)), '0'='censor', '1'=oc, '2'='death'))]
  }

# 4. Survival analysis ----

run_fit_fg <- function(data, phase, outcomes) {
  
  res_fg <- list()
  for(oc in outcomes) {
    cat(phase, oc, "\n")
    dt <- data[get(phase)==1 & get(paste0(phase,".hx.",oc))==0]
    assertthat::assert_that(all(dt[, get(paste0(phase,".time.to.censor.",oc))]>=0))
    res_fg[[oc]] <- list()
    res_fg[[oc]] <- finegray(as.formula(paste0("Surv(acute.time.to.censor.",oc,", acute.outcome.",oc,"1) ~ group")), weights=wts_acute, data=dt)
    res_fg[[oc]] <- coxph(Surv(fgstart,fgstop,fgstatus)~group,weight = fgwt,data=res_fg[[oc]])
  }
  return(res_fg)
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
  gc()
}


res_cr <- sapply(c("acute"), function(p) run_fit_fg(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)
#res_cr_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)

sink(paste0("output/T2_acute_28d_competing_jiayi.txt"), append=F, split=T)
for(oc in OUTCOMES) {
  cat(oc, "\n")
  print(rbind("-",exp(cbind(coef(res_cr$acute[[oc]]),confint(res_cr$acute[[oc]])))))
}
sink()