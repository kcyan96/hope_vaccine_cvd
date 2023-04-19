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
#cohort <- readRDS("longcovid.cohort.matched_200401.221031_overall.RDS") [index.date>="2021-02-23"]# n=3787149
cohort <- readRDS("longcovid.cohort.matched_all_210223.221031_jiayi.RDS") [index.date>="2021-02-23"]


# Common variables
cohort[, Age := as.numeric(substr(index.date,1,4)) - dob_y]
cohort[, death_date_ymd := as.Date(death_date_ymd)]
cohort[, `:=`(`Date of vaccination.1st`=as.Date(`Date of vaccination.1st`), `Date of vaccination.2nd`=as.Date(`Date of vaccination.2nd`), `Date of vaccination.3rd`=as.Date(`Date of vaccination.3rd`), `Date of vaccination.4th`=as.Date(`Date of vaccination.4th`))]
cohort[, score.cci := (dx.mi+dx.chf+dx.pvd+dx.cbd+dx.copd+dx.dementia+dx.paralysis+(dx.dm_com0&!dx.dm_com1)+dx.dm_com1*2+dx.crf*2+(dx.liver_mild&!dx.liver_modsev)+dx.liver_modsev*3+dx.ulcers+dx.ra+dx.aids*6+dx.cancer*2+dx.cancer_mets*6)]
cohort[, score.cci.age := (score.cci + pmax(pmin(ceiling((Age-49)/10),4),0))]
cohort[, severity := 0]
cohort[icu.adm_date>=COVID.date & icu.adm_date<COVID.date+7 | vent.date>=COVID.date & vent.date<COVID.date+7, severity := 1]
cohort[,severity:= as.factor(severity)]

# Remove died before index date
cohort <- cohort[is.na(death_date_ymd) | death_date_ymd >= index.date]   # n=4490681

# Remove CVD before index date
cohort <- cohort[hx.cv.major==0]   # n=4173616

# Restrict age
cohort <- cohort[Age>=18]   # n=3784326

# Vaccination status before index date
cohort[, vacc_status := do.call(paste0, lapply(c("1st","2nd","3rd","4th"), function(c) ifelse(!is.na(get(paste0("Date of vaccination.",c))) & get(paste0("Date of vaccination.",c)) < index.date, substr(get(paste0("Vaccine Brand.",c)),1,1), "")))]
cohort[vacc_status!="", time.since.last_dose := as.numeric(index.date - get(paste0("Date of vaccination.", c("1st","2nd","3rd","4th")[nchar(vacc_status)])))]

# Define groups
# group by brand dose
cohort[, group := factor(sub("^$","X",vacc_status), levels=c("X","B","S","BB","SS","BBB","SSS"))]
# group by dose
library(dplyr)
cohort$vacc_status_new <- recode(cohort$vacc_status, 'B'='1d', 'S'='1d', 'BB'='2d','SS'='2d','BBB'='3d','SSS'='3d')
cohort[, group_dose := factor(sub("^$","X",vacc_status_new), levels=c("X","1d","2d","3d"))]


# Define vaccine type
cohort[, BNT162b2 := 0]
cohort[vacc_status=="B"|vacc_status=="BB"|vacc_status=="BBB", BNT162b2 := 1]
cohort$BNT162b2 <- as.factor(cohort$BNT162b2)
cohort[, CoronaVac := 0]
cohort[vacc_status=="S"|vacc_status=="SS"|vacc_status=="SSS", CoronaVac := 1]
cohort$CoronaVac <- as.factor(cohort$CoronaVac)
cohort[, postacute:=as.numeric(is.na(death_date_ymd) | death_date_ymd >= (index.date+28))]

cohort_full <- copy(cohort)

saveRDS(cohort_full, "cohort_acute_jiayi.RDS")

cohort_postacute <- copy(cohort[postacute==1])

cohort_postacute[, index.date := (index.date+28)]

saveRDS(cohort_postacute, "cohort_postacute_jiayi.RDS")

# # Remove index date - latest vac > 180 days
# for(i in 1:4){
#   cohort[nchar(vacc_status)==i, time.since.last_dose := as.numeric(index.date - get(paste0("Date of vaccination.", c("1st","2nd","3rd","4th")[i])))]
# }
# cohort <- cohort[group=="X" | time.since.last_dose <= 180]   # n=2926865
# 
# cohort_180 <- copy(cohort)
# saveRDS(cohort_180, "cohort_acute_180.RDS")
# 
# # Keep PCR+ only
# cohort <- cohort_full[COVID.test=="HA_241"|COVID.test=="HA_244"|COVID.test=="HA_245"|COVID.test=="PCR"]   # n=487190
# 
# cohort_pcr <- copy(cohort)
# saveRDS(cohort_pcr, "cohort_acute_pcr.RDS")

rm(cohort, cohort_full, cohort_180, cohort_pcr, cohort_postacute)

gc()


## Weighting ---------------------------------------------
cohort_full <- readRDS("cohort_postacute_jiayi.RDS")
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
#tb1 <- list(postacute=table1(cohort), vaccinated=table1(cohort_))

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
summary(cohort$wts_postacute)
cohort[, vaccinated := 0]
cohort[group_dose!="X", vaccinated := 1]
cohort_ <-cohort[vaccinated==1]
cohort_[, group := factor(sub("^$","X",vacc_status_new), levels=c("1d","2d","3d"))]

library(survey)
tb1.iptw <- list(postacute=table1(cohort, "wts_postacute"), vaccinated=table1(cohort_, "wts_postacute"))

#save(W.glm, tb1, tb1.iptw, file=paste0("COMP_", COMP, PHASE, "_W.glm.RData"))
saveRDS(cohort, paste0("COMP_", COMP, PHASE, "_2.cohort_iptw_jiayi.RDS"))

library(writexl)
xl_tabs <- c(lapply(c(tb1.iptw,tb1), function(t) as.data.table(t)))
write_xlsx(xl_tabs, paste0("output/tab_", COMP, PHASE, "jiayi.xlsx"))

rm(cohort, tb1, tb1.iptw, W.glm, W.glm.acute, xl_tabs)

gc()


# Subgroup

for (sub in(c("age.l65","age.ge65","sex.male","sex.female", "cci.l4", "cci.ge4", "sev1", "sev0", "bnt", "sino"))){
  
  SA <- sub
  cohort_full <- readRDS("cohort_acute.RDS")
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
  
  
  saveRDS(cohort, paste0("COMP_", SA, "_", COMP,  PHASE, "_2.cohort_iptw.RDS"))
  gc()
}




# # Sen 180
# 
# cohort_180 <- readRDS("cohort_acute_180.RDS")
# # SELECT comparison
# COMP <- "INF_VAC_dose_May_180"
# cohort <- cohort_180[!is.na(group) & infected==1]
# cohort$group <- cohort$group_dose
# PHASE <- "_acute"
# 
# # 2. Propensity score ----
# 
# library(cobalt)
# 
# PS_def <- setDT(read_excel("codes.xlsx", sheet="PS"))[!is.na(Name)]
# PS_formula <- function(exp,lvl) paste0(exp, " ~ ", paste0(PS_def[Set<=lvl, Name], collapse=" + "))
# PS_assess <- function(out) {
#   print(summary(out))
#   print(bal.tab(out))
#   bal.plot(out, which="both")
# }
# 
# library(WeightIt)
# cohort[, acute := 1]
# W.glm <- weightit(as.formula(PS_formula("group",4)), data=cohort, estimand="ATE", method="ps")
# PS_assess(W.glm)
# cohort[, wts_acute := W.glm$weights]
# 
# saveRDS(cohort, paste0("COMP_", COMP, PHASE, "_2.cohort_iptw.RDS"))
# gc()
# 
# 
# 
# # Sen PCR+
# 
# cohort_pcr <- readRDS("cohort_acute_pcr.RDS")
# # SELECT comparison
# COMP <- "INF_VAC_dose_May_pcr"
# cohort <- cohort_pcr[!is.na(group) & infected==1]
# cohort$group <- cohort$group_dose
# PHASE <- "_acute"
# 
# # 2. Propensity score ----
# 
# library(cobalt)
# 
# PS_def <- setDT(read_excel("codes.xlsx", sheet="PS"))[!is.na(Name)]
# PS_formula <- function(exp,lvl) paste0(exp, " ~ ", paste0(PS_def[Set<=lvl, Name], collapse=" + "))
# PS_assess <- function(out) {
#   print(summary(out))
#   print(bal.tab(out))
#   bal.plot(out, which="both")
# }
# 
# library(WeightIt)
# cohort[, acute := 1]
# W.glm <- weightit(as.formula(PS_formula("group",4)), data=cohort, estimand="ATE", method="ps")
# PS_assess(W.glm)
# cohort[, wts_acute := W.glm$weights]
# 
# 
# saveRDS(cohort, paste0("COMP_", COMP, PHASE, "_2.cohort_iptw.RDS"))
# gc()
# 