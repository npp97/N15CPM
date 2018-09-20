library(compiler)
library(doParallel)
library(doBy)
library(FME)
library(hydromad)

#
niter_1=100000;
niter_2=100000;
burnin.length=50000;
maxit.sce0=50000;
maxit.sce1=100000;
lfactor=200;
ufactor=200;
lfactorf=500;
ufactorf=200;

#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://N15CS//run"
setwd(path0)

if(file.exists('.//out')){file.rename(".//out",paste(".//out",substr(file.info("out")$ctime,1,10),sep='_'))};
dir.create(".//out")

#obs_ao<-read.csv('.//input//nh4no3_fcbs.csv')
#ini_val<-read.csv('.//input//Ini_value_fcbs.csv')

obs_ao<-read.csv('.//input//nh4no3_ne.csv')
ini_val<-read.csv('.//input//Ini_value_ne.csv')


#obs_ao<-read.csv('.//input//nh4no3_nbfb.csv')
#ini_val<-read.csv('.//input//ini_value_nbfb.csv')


#obs_ao<-read.csv('.//input//nh4no3_cbssuc_90raw.csv')
#ini_val<-read.csv('.//input//Ini_value_cbssuc_90raw.csv')


lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15')]
obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

#rndgrp<-c(-1234,747,1979,470,3000,6000,9000,11310)
rndgrp<-floor(rnorm(90)*10000)
rndgrp<-floor(runif(15,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)

no_cores <- detectCores() - 0
ncl<-max(min(no_cores,lgrpPal),2)

#TODO::prepare y0_label and obs_label



loadcmp("subs_all.Rc")
# 
#pars0 =c(C_NRec_2_NH4 = 0.13118,
#           C_NH4_2_NLab = 0.3318,
#           K_NLab_2_NH4 = 0.04396,
#           K_NO3_2_NRec = 0.000014,
#           C_NRec_2_NO3 = 0.04312,
#           K_NH4ads_2_NH4 = 0.12586,
#           C_NO3_2_NH4 =  0.0037,
#           K_NH4_2_NO3 = 0.52122,
#           K_NO3sto_2_NO3 = 0.000000001
#  );


#loadcmp("muller2003_batch_CS1.Rc")

#loadcmp("zhu2016_batch_p1.Rc")
#pars0 for N15EXPNE model:zhang2013a
# pars0 = c(
  # C_NRec_2_NH4 = 0.11170245,
  # K_NH4_2_NRec = 0.01774379,
  # C_NH4_2_NLab = 0.190611353,
  # K_NLab_2_NH4 = 0.0072679,
  # K_NO3_2_NRec = 0.000896081,
  # C_NRec_2_NO3 = 0.122305291,
  # K_NH4ads_2_NH4 = 0.021196383,
  # K_NH4_2_NH4ads = 0.000012044348,
  # C_NO3_2_NH4 =  0.012601542,
  # K_NH4_2_NO3 =0.021858307,
  # K_NO3_2_NO3sto = 0.0005222072,
  # K_NO3sto_2_NO3 = 0.004020091
# )
# pars0 = c(
# C_NRec_2_NH4 = 0.0535627,
# K_NH4_2_NRec = 3.12869E-10,
# C_NH4_2_NLab = 0.30305959,
# K_NLab_2_NH4 = 0.001473668,
# K_NO3_2_NRec = 0.000115873,
# C_NRec_2_NO3 = 0.273685075,
# K_NH4ads_2_NH4 = 0.028805643,
# K_NH4_2_NH4ads = 0.001361153,
# C_NO3_2_NH4 = 0.019215892,
# K_NH4_2_NO3 = 0.025135268,
# K_NO3_2_NO3sto = 0.00746159,
# K_NO3sto_2_NO3 = 0.026350967
# )

# loadcmp("zhang2013a_batch_p2nc.Rc")
#FOR zhang2013am
 pars0=c(
   C_NRec_2_NH4 = 0.0491716,
   Vmax_NH4_2_NRec = 0.0006520959,
   Km_NH4_2_NRec = 63.70075,
   Vmax_NH4_2_NLab = 0.01235169,
   Km_NH4_2_NLab = 49.1142,
   Vmax_NLab_2_NH4 = 0.01572523,
   Km_NLab_2_NH4 = 27.75582,
   Vmax_NO3_2_NRec = 0.007118786,
   Km_NO3_2_NRec = 129.8529,
   C_NRec_2_NO3 = 0.1631907,
   K_NH4ads_2_NH4 = 0.001804257,
   K_NH4_2_NH4ads = 0.01721704,
   Vmax_NO3_2_NH4 = 0.0003459241,
   Km_NO3_2_NH4 = 18.1045,
   Vmax_NH4_2_NO3 = 0.009894824,
   Km_NH4_2_NO3 = 19.2277,
   K_NO3_2_NO3sto = 0.0005222072,
   K_NO3sto_2_NO3 = 0.004020091
 )
 loadcmp("zhang2013am_batch_CS1.Rc")
#source("D://N15CS//srcs//batch-p//zhang2013am_batch_CS1.R")
#For Muller2003m
# pars0=c(
#   C_NRec_2_NH4 = 0.0491716,
#   Vmax_NH4_2_NLab = 0.01235169,
#   Km_NH4_2_NLab = 491.142,
#   Vmax_NLab_2_NH4 = 0.1572523,
#   Km_NLab_2_NH4 = 27.5582,
#   Vmax_NO3_2_NRec = 0.007118786,
#   Km_NO3_2_NRec = 1298.529,
#   C_NRec_2_NO3 = 0.1631907,
#   K_NH4ads_2_NH4 = 0.1804257,
#   Vmax_NO3_2_NH4 = 0.00459241,
#   Km_NO3_2_NH4 = 181.045,
#   Vmax_NH4_2_NO3 = 0.9894824,
#   Km_NH4_2_NO3 = 19.2277,
#   K_NO3sto_2_NO3 = 0.00004020091
# )
# loadcmp("muller2003m_batch_CS.Rc")
# 
# # 
# # #For zhang2013am11
# pars0=c(
#   C_NRec_2_NH4 = 0.0491716,
#   Vmax_NH4_2_NRec = 0.0006520959,
#   Km_NH4_2_NRec = 637.0075,
#   Vmax_NH4_2_NLab = 0.01235169,
#   Km_NH4_2_NLab = 491.142,
#   Vmax_NLab_2_NH4 = 0.01572523,
#   Km_NLab_2_NH4 = 277.5582,
#   Vmax_NO3_2_NRec = 0.007118786,
#   Km_NO3_2_NRec = 1298.529,
#   C_NRec_2_NO3 = 0.1631907,
#   K_NH4ads_2_NH4 = 0.001804257,
#   K_NH4_2_NH4ads = 0.01721704,
#   Vmax_NO3_2_NH4 = 0.0003459241,
#   Km_NO3_2_NH4 = 181.045,
#   Vmax_NH4_2_NO3 = 0.009894824,
#   Km_NH4_2_NO3 = 193.2277,
#   K_NO3sto_2_NO3 = 0.004020091
# )

##pars0 for N15EXP-CBS,zhang2013a
# pars0 = c(
#   C_NRec_2_NH4 = 0.031923477,
#   K_NH4_2_NRec = 0.0008489,
#   C_NH4_2_NLab = 0.028623952,
#   K_NLab_2_NH4 = 0.000239385,
#   K_NO3_2_NRec = 0.000131491,
#   C_NRec_2_NO3 = 0.056204438,
#   K_NH4ads_2_NH4 = 0.002753537,
#   K_NH4_2_NH4ads = 0.000915349,
#   C_NO3_2_NH4 =  0.001863671,
#   K_NH4_2_NO3 = 0.003059998
# )
  #loadcmp("zhu2016_batch_p2.Rc")
#loadcmp("zhang2013am_batch_p2nc.Rc")
#loadcmp("zhang2013am11_batch_p2nc.Rc")
  #loadcmp("zhang2013e_batch_p2.Rc") 
#loadcmp("muller2003m_batch_p2nc.Rc")



#try(loadcmp("muller2003_batch_p1.Rc"))
#try(loadcmp("muller2003_45_batch_p1.Rc"),silent = F);
# try(loadcmp("muller2003a_batch_p1.Rc"),silent = F);
# try(loadcmp("muller2003b_batch_p1.Rc"),silent = F);
# try(loadcmp("muller2003c_batch_p1.Rc"),silent = F);
# try(loadcmp("zhang2013a_batch_p1.Rc"),silent = F);
# try(loadcmp("zhang2013_batch_p1.Rc"),silent = F);
gc()
##################

# flst<-dir(path="..//out",pattern = ".RData")
# 
# nfile<-length(flst)
# setwd("..//out")
# 
# sink('compare_obs_sim.txt')
# for (i in 1:nfile){
#   
#   load(flst[i])
#   sss.nh4<-as.data.frame(sss.nh4)
#   sss.no3<-as.data.frame(sss.no3)
#   sss.nh4no3<-as.data.frame(sss.nh4no3)
#   corcott<-c(cor(obs_nh4$NH4,sss.nh4$NH4),cor.test(obs_nh4$NH4,sss.nh4$NH4)$p.value,
#              cor(obs_nh4$NO3,sss.nh4$NO3),cor.test(obs_nh4$NO3,sss.nh4$NO3)$p.value,
#              cor(obs_nh4$NH4_15,sss.nh4$NH4_15),cor.test(obs_nh4$NH4_15,sss.nh4$NH4_15)$p.value,
#              cor(obs_nh4$NO3_15,sss.nh4$NO3_15),cor.test(obs_nh4$NO3_15,sss.nh4$NO3_15)$p.value,
#              cor(obs_no3$NH4,sss.no3$NH4),cor.test(obs_no3$NH4,sss.no3$NH4)$p.value,
#              cor(obs_no3$NO3,sss.no3$NO3),cor.test(obs_no3$NO3,sss.no3$NO3)$p.value,
#              cor(obs_no3$NH4_15,sss.no3$NH4_15),cor.test(obs_no3$NH4_15,sss.no3$NH4_15)$p.value,
#              cor(obs_no3$NO3_15,sss.no3$NO3_15),cor.test(obs_no3$NO3_15,sss.no3$NO3_15)$p.value,
#              cor(obs_nh4no3$NH4,sss.nh4no3$NH4),cor.test(obs_nh4no3$NH4,sss.nh4no3$NH4)$p.value,
#              cor(obs_nh4no3$NO3,sss.nh4no3$NO3),cor.test(obs_nh4no3$NO3,sss.nh4no3$NO3)$p.value,
#              cor(obs_nh4no3$NH4_15,sss.nh4no3$NH4_15),cor.test(obs_nh4no3$NH4_15,sss.nh4no3$NH4_15)$p.value,
#              cor(obs_nh4no3$NO3_15,sss.nh4no3$NO3_15),cor.test(obs_nh4no3$NO3_15,sss.nh4no3$NO3_15)$p.value)
#   print(paste(flst[i],sfsmisc::vcat(format(corcott,digits=6),sep=','),sep=','))
#   
# }
# sink()
# 
# print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
# for (i in 1:nfile){
#   # rm(list=ls());gc();
#   load(flst[i])
#   rates<-c(NRec_2_NH4 = as.numeric(MCMC$bestpar['C_NRec_2_NH4']),
#            NH4_2_NRec = mean(MCMC$bestpar['C_NH4_2_NRec']),
#            NH4_2_NLab =  mean(MCMC$bestpar['C_NH4_2_NLab']),
#            NLab_2_NH4 = mean(MCMC$bestpar['K_NLab_2_NH4']*ss$NLab),
#            NO3_2_NRec = mean(MCMC$bestpar['C_NO3_2_NRec']),
#            NRec_2_NO3 = as.numeric(MCMC$bestpar['C_NRec_2_NO3']),
#            NO3sto_2_NO3 = as.numeric(MCMC$bestpar['C_NO3sto_2_NO3']),
#            NO3_2_NO3sto =mean(MCMC$bestpar['C_NO3_2_NO3sto']),
#            NH4ads_2_NH4 = mean(MCMC$bestpar['K_NH4ads_2_NH4']*ss$NH4ads),
#            NH4_2_NH4ads = mean(MCMC$bestpar['K_NH4_2_NH4ads']*ss$NH4),
#            NO3_2_NH4 = mean(MCMC$bestpar['K_NO3_2_NH4']*ss$NO3),
#            NH4_2_NO3 = mean(MCMC$bestpar['K_NH4_2_NO3']*ss$NH4),
#            NRec_2_N2O = as.numeric(MCMC$bestpar['C_NRec_2_N2O']),
#            NO3_2_N2O = mean(MCMC$bestpar['K_NO3_2_N2O']*ss$NO3),
#            NO3_2_N2 = mean(MCMC$bestpar['K_NO3_2_N2']*ss$NO3),
#            NRec_2_N2=as.numeric(MCMC$bestpar['C_NRec_2_N2']),
#            NH4_2_NH3=mean(MCMC$bestpar['K_NH4_2_NH3']*ss$NH4)
#   )
#   print(paste(flst[i],sfsmisc::vcat(format(rates,digits=3),sep=','),sep=','))
# }
# 
# print(names(rates))
# sink()
# 
# setwd("..//")
# 
file.rename(".//out",paste(".//out_ne_",substr(file.info("out")$ctime,1,10),sep='_'));

if(file.exists('.//out_cbs')){file.rename(".//out_cbs",paste(".//out_cbs",substr(file.info("out_cbs")$ctime,1,10),sep='_'))};
dir.create(".//out_cbs")

obs_ao<-read.csv('.//input//nh4no3_fcbs.csv')
ini_val<-read.csv('.//input//Ini_value_fcbs.csv')

#obs_ao<-read.csv('.//input//nh4no3_nbfb.csv')
#ini_val<-read.csv('.//input//ini_value_nbfb.csv')


#obs_ao<-read.csv('.//input//nh4no3_cbssuc_90raw.csv')
#ini_val<-read.csv('.//input//Ini_value_cbssuc_90raw.csv')


lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15')]
obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

#rndgrp<-c(-1234,747,1979,470,3000,6000,9000,11310)
rndgrp<-floor(rnorm(90)*10000)
rndgrp<-floor(runif(15,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)
pars0=c(
   C_NRec_2_NH4 = 0.0491716,
   Vmax_NH4_2_NRec = 0.0006520959,
   Km_NH4_2_NRec = 63.70075,
   Vmax_NH4_2_NLab = 0.01235169,
   Km_NH4_2_NLab = 49.1142,
   Vmax_NLab_2_NH4 = 0.01572523,
   Km_NLab_2_NH4 = 27.75582,
   Vmax_NO3_2_NRec = 0.007118786,
   Km_NO3_2_NRec = 129.8529,
   C_NRec_2_NO3 = 0.1631907,
   K_NH4ads_2_NH4 = 0.001804257,
   K_NH4_2_NH4ads = 0.01721704,
   Vmax_NO3_2_NH4 = 0.0003459241,
   Km_NO3_2_NH4 = 18.1045,
   Vmax_NH4_2_NO3 = 0.009894824,
   Km_NH4_2_NO3 = 19.2277,
   K_NO3_2_NO3sto = 0.0005222072,
   K_NO3sto_2_NO3 = 0.004020091
 )
 
loadcmp("subs_all.Rc")
loadcmp("zhang2013am_batch_CS1.Rc")

file.rename(".//out",paste(".//out_cbs",substr(file.info("out")$ctime,1,10),sep='_'));