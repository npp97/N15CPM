library(compiler)
library(doParallel)
library(doBy)
library(FME)
library(hydromad)

findpks<-function(parmset){
	icol<-ncol(parmset)
	parr<-parmset[1,];
	for (i in 1:ncol(parmset)){
		dp<-density(parmset[,i])
		ii<-which(dp$y == max(dp$y))
		parr[i]<-dp$x[ii]
	}
	return(parr)
}
#
lib_nm<-'zhang2013a.dll'
mod_type='Z13a'
Lb3<-c('NH4','NO3','Urea')

niter_1=10000;
niter_2=10000;
burnin.length=5000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=1000;
ufactor=1000;
lfactorf=1000;
ufactorf=1000;
out_dir='.//out_pre//SUSF-6090AVG-ZHANG2013a'
#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://N15CS//run"
setwd(path0)

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

#obs_ao<-read.csv('.//input//nh4no3_fcbs.csv')
#ini_val<-read.csv('.//input//Ini_value_fcbs.csv')

#obs_ao<-read.csv('.//input//nh4no3_nbfb.csv')
#ini_val<-read.csv('.//input//ini_value_nbfb.csv')

# obs_ao<-read.csv('.//input//nh4no3_cbssuc_90raw.csv')
# ini_val<-read.csv('.//input//Ini_value_cbssuc_90raw.csv')

# obs_ao<-read.csv('.//input//nh4no3_cbssuc_60raw.csv')
# ini_val<-read.csv('.//input//Ini_value_cbssuc_60raw.csv')

obs_ao<-read.csv('.//input//nh4no3_cbssuc_stgavg.csv')
ini_val<-read.csv('.//input//Ini_value_cbssuc_stgavg.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15'))]
obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

#rndgrp<-c(-1234,747,1979,470,3000,6000,9000,11310)
rndgrp<-floor(rnorm(90)*10000)
rndgrp<-floor(runif(45,-20000,20000))
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

##pars0 for N15EXP-CBS,zhang2013a
 pars0 = c(
   C_NRec_2_NH4 = 0.031923477,
   K_NH4_2_NRec = 0.0008489,
   C_NH4_2_NLab = 0.028623952,
   K_NLab_2_NH4 = 0.000239385,
   K_NO3_2_NRec = 0.000131491,
   C_NRec_2_NO3 = 0.056204438,
   K_NH4ads_2_NH4 = 0.002753537,
   K_NH4_2_NH4ads = 0.000915349,
   C_NO3_2_NH4 =  0.001863671,
   K_NH4_2_NO3 = 0.003059998
 )
 
 
 #define C_NRec_2_NH4 parms[0] //3 M_NREc
 #define K_NH4_2_NRec parms[1] // I_NH4
 #define C_NH4_2_NLab parms[2] //2 I_NH4
 #define K_NLab_2_NH4 parms[3] //1 M_Nlab
 #define K_NO3_2_NRec parms[4] //4 I_NO3
 #define C_NRec_2_NO3 parms[5] //5 O_NRec
 #define K_NH4ads_2_NH4 parms[6] //8 R_NH4a
 #define K_NH4_2_NH4ads parms[7] //A_NH4
 #define C_NO3_2_NH4 parms[8] //7 D_NO3
 #define K_NH4_2_NO3 parms[9]  //6 O_NH4 
 
loadcmp("Z13a_batch_CS1.Rc")


flst<-dir(pattern='.RData',path=out_dir,full.names=T)

load(flst[1])
pars.Z13a<-rbind(MCMC$bestpar,MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar)

for (i in 2:length(flst)){
load(flst[i])
if((min(sign(diag(cor_ss_nh4_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NH4$bestpar)
if((min(sign(diag(cor_ss_no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_no3_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NO3$bestpar)
if((min(sign(diag(cor_ss_nh4no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4no3_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NH4NO3$bestpar)

if((min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))),na.rm=T)>0) & (!is.na(min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))))))) pars.Z13a<-rbind(pars.Z13a,MCMC$bestpar)
}

pars0<-apply(pars.Z13a,2,median)


niter_1=150000;
niter_2=150000;
burnin.length=50000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=200;
ufactor=200;
lfactorf=500;
ufactorf=200;

out_dir='.//out//SUSF-6090AVG-ZHANG2013a'

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

rndgrp<-floor(runif(30,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)

no_cores <- detectCores() - 0
ncl<-max(min(no_cores,lgrpPal),2)

loadcmp("Z13a_batch_CS1.Rc")

########################################################################################################################
#########################################################################################################################
lib_nm<-'zhang2013a.dll'
mod_type='Z13a'
Lb3<-c('NH4','NO3','Urea')

niter_1=10000;
niter_2=10000;
burnin.length=5000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=1000;
ufactor=1000;
lfactorf=1000;
ufactorf=1000;
out_dir='.//out_pre//SUSF-60-ZHANG2013a'
#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://N15CS//run"
setwd(path0)

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

#obs_ao<-read.csv('.//input//nh4no3_fcbs.csv')
#ini_val<-read.csv('.//input//Ini_value_fcbs.csv')

#obs_ao<-read.csv('.//input//nh4no3_nbfb.csv')
#ini_val<-read.csv('.//input//ini_value_nbfb.csv')

# obs_ao<-read.csv('.//input//nh4no3_cbssuc_90raw.csv')
# ini_val<-read.csv('.//input//Ini_value_cbssuc_90raw.csv')

obs_ao<-read.csv('.//input//nh4no3_cbssuc_60raw.csv')
ini_val<-read.csv('.//input//Ini_value_cbssuc_60raw.csv')

#obs_ao<-read.csv('.//input//nh4no3_cbssuc_stgavg.csv')
#ini_val<-read.csv('.//input//Ini_value_cbssuc_stgavg.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15'))]
obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

#rndgrp<-c(-1234,747,1979,470,3000,6000,9000,11310)
rndgrp<-floor(rnorm(90)*10000)
rndgrp<-floor(runif(45,-20000,20000))
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

##pars0 for N15EXP-CBS,zhang2013a
 pars0 = c(
   C_NRec_2_NH4 = 0.031923477,
   K_NH4_2_NRec = 0.0008489,
   C_NH4_2_NLab = 0.028623952,
   K_NLab_2_NH4 = 0.000239385,
   K_NO3_2_NRec = 0.000131491,
   C_NRec_2_NO3 = 0.056204438,
   K_NH4ads_2_NH4 = 0.002753537,
   K_NH4_2_NH4ads = 0.000915349,
   C_NO3_2_NH4 =  0.001863671,
   K_NH4_2_NO3 = 0.003059998
 )
 
  
loadcmp("Z13a_batch_CS1.Rc")


flst<-dir(pattern='.RData',path=out_dir,full.names=T)

load(flst[1])
pars.Z13a<-rbind(MCMC$bestpar,MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar)

for (i in 2:length(flst)){
load(flst[i])
if((min(sign(diag(cor_ss_nh4_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NH4$bestpar)
if((min(sign(diag(cor_ss_no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_no3_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NO3$bestpar)
if((min(sign(diag(cor_ss_nh4no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4no3_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NH4NO3$bestpar)

if((min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))),na.rm=T)>0) & (!is.na(min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))))))) pars.Z13a<-rbind(pars.Z13a,MCMC$bestpar)
}

pars0<-apply(pars.Z13a,2,median)

niter_1=150000;
niter_2=150000;
burnin.length=50000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=200;
ufactor=200;
lfactorf=500;
ufactorf=200;

out_dir='.//out//SUSF-60-ZHANG2013a'

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

rndgrp<-floor(runif(30,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)

no_cores <- detectCores() - 0
ncl<-max(min(no_cores,lgrpPal),2)

loadcmp("Z13a_batch_CS1.Rc")

###############################################################################################################################################
###############################################################################################################################################
lib_nm<-'zhang2013a.dll'
mod_type='Z13a'
Lb3<-c('NH4','NO3','Urea')

niter_1=10000;
niter_2=10000;
burnin.length=5000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=1000;
ufactor=1000;
lfactorf=1000;
ufactorf=1000;

out_dir='.//out_pre//SUSF-90-ZHANG2013a'

#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://N15CS//run"
setwd(path0)

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

#obs_ao<-read.csv('.//input//nh4no3_fcbs.csv')
#ini_val<-read.csv('.//input//Ini_value_fcbs.csv')

#obs_ao<-read.csv('.//input//nh4no3_nbfb.csv')
#ini_val<-read.csv('.//input//ini_value_nbfb.csv')

obs_ao<-read.csv('.//input//nh4no3_cbssuc_90raw.csv')
ini_val<-read.csv('.//input//Ini_value_cbssuc_90raw.csv')

# obs_ao<-read.csv('.//input//nh4no3_cbssuc_60raw.csv')
# ini_val<-read.csv('.//input//Ini_value_cbssuc_60raw.csv')

#obs_ao<-read.csv('.//input//nh4no3_cbssuc_stgavg.csv')
#ini_val<-read.csv('.//input//Ini_value_cbssuc_stgavg.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15'))]
obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

#rndgrp<-c(-1234,747,1979,470,3000,6000,9000,11310)
rndgrp<-floor(rnorm(90)*10000)
rndgrp<-floor(runif(45,-20000,20000))
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

##pars0 for N15EXP-CBS,zhang2013a
 pars0 = c(
   C_NRec_2_NH4 = 0.031923477,
   K_NH4_2_NRec = 0.0008489,
   C_NH4_2_NLab = 0.028623952,
   K_NLab_2_NH4 = 0.000239385,
   K_NO3_2_NRec = 0.000131491,
   C_NRec_2_NO3 = 0.056204438,
   K_NH4ads_2_NH4 = 0.002753537,
   K_NH4_2_NH4ads = 0.000915349,
   C_NO3_2_NH4 =  0.001863671,
   K_NH4_2_NO3 = 0.003059998
 )
 
 
loadcmp("Z13a_batch_CS1.Rc")


flst<-dir(pattern='.RData',path=out_dir,full.names=T)

load(flst[1])
pars.Z13a<-rbind(MCMC$bestpar,MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar)

for (i in 2:length(flst)){
load(flst[i])
if((min(sign(diag(cor_ss_nh4_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NH4$bestpar)
if((min(sign(diag(cor_ss_no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_no3_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NO3$bestpar)
if((min(sign(diag(cor_ss_nh4no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4no3_obs)))))) pars.Z13a<-rbind(pars.Z13a,MCMC.NH4NO3$bestpar)

if((min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))),na.rm=T)>0) & (!is.na(min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))))))) pars.Z13a<-rbind(pars.Z13a,MCMC$bestpar)
}

pars0<-apply(pars.Z13a,2,median)


niter_1=150000;
niter_2=150000;
burnin.length=50000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=200;
ufactor=200;
lfactorf=500;
ufactorf=200;

out_dir='.//out//SUSF-90-ZHANG2013a'

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

rndgrp<-floor(runif(30,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)

no_cores <- detectCores() - 0
ncl<-max(min(no_cores,lgrpPal),2)

loadcmp("Z13a_batch_CS1.Rc")