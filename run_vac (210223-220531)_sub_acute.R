setwd("/mnt/Others/Dataset/")
library(data.table)
library(survey)
library(epitools)

PHASE <- "_acute"

# 3. Censoring ----
COMP <- "INF_VAC_dose_May"
cohort<-readRDS(paste0("LONGCOVID_CC/COMP_", COMP, PHASE, "_2.cohort_iptw.RDS"))

OUTCOMES <- 
c("cv.major" ,"cv.stroke" ,"cv.cad" ,"cv.hf" ,"cv.mi" ,"cv_death", "death")
#c("cv.major" ,"cv.stroke" ,"cv.tia" ,"cv.af" ,"cv.afl" ,"cv.pericarditis" ,"cv.myocarditis" ,"cv.cad" ,"cv.acs" ,"cv.mi" ,"cv.stangina" ,"cv.unstangina" ,"cv.hf" ,"cv.cardarrest" ,"cv.cardshock" ,"cv.pe" ,"cv.dvt", "cv_death", "death")
#OUTCOMES <- c("cv.major", "cv.cad", "cv.stroke", "cv.hf", "cv.af", "cv.mi", "cv_death", "death")        # old c("cvd_comp","cad_comp","stroke_comp","hf_comp","dm_neuro","dm_retino","dm_nephro","esrd","cv_death","death")
cohort[, `:=`(hx.cv_death = 0, hx.death = 0)]
cohort[, postacute.index.date := index.date+21]

# Acute phase (within 21 days)
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

res_65 <- sapply(c("acute"), function(p) run_fit(cohort[Age<65], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_66 <- sapply(c("acute"), function(p) run_fit(cohort[Age>=65], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_m <- sapply(c("acute"), function(p) run_fit(cohort[sex=="M"], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_f <- sapply(c("acute"), function(p) run_fit(cohort[sex=="F"], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_cci3 <- sapply(c("acute"), function(p) run_fit(cohort[score.cci.age<4], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_cci4 <- sapply(c("acute"), function(p) run_fit(cohort[score.cci.age>=4], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_sev1 <- sapply(c("acute"), function(p) run_fit(cohort[severity==1], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_sev0 <- sapply(c("acute"), function(p) run_fit(cohort[severity==0], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_bnt <- sapply(c("acute"), function(p) run_fit(cohort[BNT162b2==1], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_sino <- sapply(c("acute"), function(p) run_fit(cohort[CoronaVac==1], p, OUTCOMES), simplify=F, USE.NAMES=T)


sink(paste0("LONGCOVID_CC/output/Sub_COMP_", COMP, PHASE, "_log.txt"), append=F, split=T)
print("Age<65")
res_65
print("Age>=65")
res_66
print("Sex=M")
res_m
print("Sex=F")
res_f
print("CCI<4")
res_cci3
print("CCI>=4")
res_cci4
print("severity==1")
res_sev1
print("severity==0")
res_sev0
print("BNT162b2==1")
res_bnt
print("CoronaVac==1")
res_sino
gc()

res_65_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[Age<65], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_66_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[Age>=65], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_m_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[sex=="M"], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_f_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[sex=="F"], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_cci3_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[score.cci.age<4], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_cci4_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[score.cci.age>=4], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_sev1_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[severity==1], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_sev0_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[severity==0], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_bnt_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[group=="X"|BNT162b2==1], p, OUTCOMES), simplify=F, USE.NAMES=T)
res_sino_IR <- sapply(c("acute"), function(p) run_fit_IR(cohort[group=="X"|CoronaVac==1], p, OUTCOMES), simplify=F, USE.NAMES=T)

sink(paste0("LONGCOVID_CC/output/Sub_COMP_IR_", COMP, PHASE, "_log.txt"), append=F, split=T)
print("Age<65")
res_65_IR
print("Age>=65")
res_66_IR
print("Sex=M")
res_m_IR
print("Sex=F")
res_f_IR
print("CCI<4")
res_cci3_IR
print("CCI>=4")
res_cci4_IR
print("severity==1")
res_sev1_IR
print("severity==0")
res_sev0_IR
print("BNT162b2==1")
res_bnt_IR
print("CoronaVac==1")
res_sino_IR

gc()
