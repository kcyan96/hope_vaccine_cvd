setwd("/mnt/Others/Dataset/")
library(data.table)
library(readxl)
library(survey)
library(dplyr)
library(survival)
library(epitools)

PHASE <- "_postacute"

# 3. Censoring ---- 
COMP <- "INF_VAC_dose_May"
cohort<-readRDS(paste0("LONGCOVID_CC/COMP_", COMP, PHASE, "_2.cohort_iptw.RDS"))

#OUTCOMES <- c( "death") 
OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death", "death")
#c("cv.major" ,"cv.stroke" ,"cv.tia" ,"cv.af" ,"cv.afl" ,"cv.pericarditis" ,"cv.myocarditis" ,"cv.cad" ,"cv.acs" ,"cv.mi" ,"cv.stangina" ,"cv.unstangina" ,"cv.hf" ,"cv.cardarrest" ,"cv.cardshock" ,"cv.pe" ,"cv.dvt", "cv_death", "death")
#OUTCOMES <- c("cv.major", "cv.cad", "cv.stroke", "cv.hf", "cv.af", "cv.mi", "cv_death", "death")        # old c("cvd_comp","cad_comp","stroke_comp","hf_comp","dm_neuro","dm_retino","dm_nephro","esrd","cv_death","death")

#cohort[, `:=`(hx.death = 0)]
cohort[, `:=`(hx.cv_death = 0, hx.death = 0)]
cohort[, postacute.index.date := index.date]


for(oc in OUTCOMES) {
  cohort[,paste0("outcome.",oc,".date"):=(as.character(get(paste0("outcome.",oc,".date"))))]
  #cohort[,paste0("outcome.",oc,".date"):=(as.numeric(get(paste0("outcome.",oc,".date"))))]
  cohort[,paste0("outcome.",oc,".date"):=(as.Date(get(paste0("outcome.",oc,".date")),origin='1970-01-01'))]
  
}

# Postacute phase (after 21 days)
for(oc in OUTCOMES) {
  cohort[, c(paste0("postacute.hx.",oc)) := as.numeric(get(paste0("hx.",oc)) | (!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("outcome.",oc,".date")) < postacute.index.date))]
  cohort[, c(paste0("postacute.censor.date.", oc)) := pmin(get(paste0("outcome.",oc,".date")), death_date_ymd, as.Date("2022-08-15"), na.rm=T)]
  cohort[, c(paste0("postacute.time.to.censor.", oc)) := as.numeric(get(paste0("postacute.censor.date.", oc)) - postacute.index.date) + 1]
  cohort[, c(paste0("postacute.outcome.", oc)) := as.numeric(!is.na(get(paste0("outcome.",oc,".date"))) & get(paste0("postacute.censor.date.", oc))==get(paste0("outcome.",oc,".date")) & get(paste0("outcome.",oc,".date"))>=postacute.index.date)]
}
rm(oc)

## Competing risk using Fine-gray model
OUTCOMES <- c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv_death")

for(oc in OUTCOMES) {
  cohort[postacute.outcome.death==1& get(paste0("postacute.outcome.", oc))==0, (paste0("postacute.outcome.", oc, "1")) := 2]
  cohort[, paste0("postacute.outcome.", oc, "1") := as.factor(recode(get(paste0("postacute.outcome.", oc)), '0'='censor', '1'=oc, '2'='death'))]
}

# 4. Survival analysis ----


run_fit_fg <- function(data, phase, outcomes) {
  
 res_fg <- list()
  for(oc in outcomes) {
    cat(phase, oc, "\n")
    dt <- data[get(phase)==1 & get(paste0(phase,".hx.",oc))==0]
    assertthat::assert_that(all(dt[, get(paste0(phase,".time.to.censor.",oc))]>=0))
    res_fg[[oc]] <- list()
    res_fg[[oc]] <- finegray(as.formula(paste0("Surv(postacute.time.to.censor.",oc,", postacute.outcome.",oc,"1) ~ group")), weights=wts_postacute, data=dt)
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

res_cr <- sapply(c("postacute"), function(p) run_fit_fg(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)
res_cr_IR <- sapply(c("postacute"), function(p) run_fit_IR(cohort, p, OUTCOMES), simplify=F, USE.NAMES=T)

sink(paste0("LONGCOVID_CC/output/CR_", COMP, PHASE, "_log.txt"), append=F, split=T)
print("Competing risk")
res_cr
print("cr IR")
res_cr_IR

gc()
#unlink(paste0("LONGCOVID_CC/output/CR_", COMP, PHASE, "_log.txt"))
