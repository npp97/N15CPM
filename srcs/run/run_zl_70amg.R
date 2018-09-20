#run all for Zhaoting soil water 90%

library(compiler)
library(doParallel)
library(doBy)
library(FME)
library(hydromad)
library(missForest)

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
lib_nm<-'zhang2013amg.dll'
mod_type='Z13amg'
niter_1=10000;
niter_2=10000;
burnin.length=5000;
maxit.sce0=5000;
maxit.sce1=5000;
lfactor=1000;
ufactor=1000;
lfactorf=1000;
ufactorf=1000;
Lb3<-c('NH4','NO3')
out_dir='.//out_pre_zl-70-ZHANG2013amg'
#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://N15CSg_zl//run"
setwd(path0)

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

obs_ao<-read.csv('.//input//nh4no3_zl_imp.csv')
ini_val<-read.csv('.//input//Ini_value_zl.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15','N2O','N2O_15'))]
set.seed(1000);obs_al<-missForest(obs_al, maxiter = 100, ntree = 2000)$ximp

obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

rndgrp<-floor(runif(30,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)

no_cores <- detectCores() - 0
ncl<-max(min(no_cores,lgrpPal),2)

#TODO::prepare y0_label and obs_label


loadcmp("subs_2gas.Rc")
#

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
   C_NO3sto_2_NO3 = 0.004020091,
   K_NO3_2_N2O  = 0.0001,
   K_NH4_2_N2O  = 0.00001,
   K_NRec_2_N2O = 0.000001
 )
 
loadcmp("Z13amg15_batch_zl.Rc")
 

flst<-dir(pattern='.RData',path=out_dir,full.names=T)

load(flst[1])
pars.Z13am<-rbind(MCMC$bestpar,MCMC.NH4$bestpar,MCMC.NO3$bestpar)

for (i in 2:length(flst)){
load(flst[i])
if((min(sign(diag(cor_ss_nh4_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4_obs)))))) pars.Z13am<-rbind(pars.Z13am,MCMC.NH4$bestpar)
if((min(sign(diag(cor_ss_no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_no3_obs)))))) pars.Z13am<-rbind(pars.Z13am,MCMC.NO3$bestpar)

if((min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_no3_obs))),na.rm=T)>0) & 
(!is.na(min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_no3_obs))))))) pars.Z13am<-rbind(pars.Z13am,MCMC$bestpar)
}

pars0<-apply(pars.Z13am,2,median)


niter_1=50000;
niter_2=50000;
burnin.length=25000;
maxit.sce0=50000;
maxit.sce1=50000;
lfactor=200;
ufactor=200;
lfactorf=500;
ufactorf=200;

out_dir='.//out_zl_70-ZHANG2013amg'

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

loadcmp("Z13amg15_batch_CS1.Rc")
###############################################################
niter_1=100000;
niter_2=100000;
burnin.length=25000;
lfactor=2;
ufactor=2;

path0<-"D://N15CSg_zl//run"
setwd(path0)

flst<-dir(pattern="*.RData",path=out_dir,full.names=T)

ii_yq_1<-grep('YQ-1',flst)
ii_yq_2<-grep('YQ-2',flst)

ii_jq_1<-grep('JQ-1',flst)
ii_jq_2<-grep('JQ-2',flst)

ii_mh_1<-grep('MH-2',flst)
ii_mh_2<-grep('MH-2',flst)

rst_yq_1<-findBST(flst[ii_yq_1])
rst_yq_2<-findBST(flst[ii_yq_2])

rst_jq_1<-findBST(flst[ii_jq_1])
rst_jq_2<-findBST(flst[ii_jq_2])

rst_mh_1<-findBST(flst[ii_mh_1])
rst_mh_2<-findBST(flst[ii_mh_2])

rst_yq_1f<-getpar_fnal_z13amg(arst=rst_yq_1,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_yq_2f<-getpar_fnal_z13amg(arst=rst_yq_2,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

rst_jq_1f<-getpar_fnal_z13amg(arst=rst_jq_1,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_jq_2f<-getpar_fnal_z13amg(arst=rst_jq_2,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

rst_mh_1f<-getpar_fnal_z13amg(arst=rst_mh_1,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_mh_2f<-getpar_fnal_z13amg(arst=rst_mh_2,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)


out_dir_f='.//out_f_zl-70-ZHANG2013amg'

if(file.exists(out_dir_f)){file.rename(out_dir_f,paste(out_dir_f,substr(file.info(out_dir_f)$ctime,1,10),sep='_'))};
dir.create(out_dir_f)

setwd(out_dir_f)

gen_fig_final(rst_yq_1,'YQ-1-70')
gen_fig_final(rst_yq_2,'YQ-2-70')

gen_fig_final(rst_jq_1,'JQ-1-70')
gen_fig_final(rst_jq_2,'JQ-2-70')

gen_fig_final(rst_mh_1,'MH-1-70')
gen_fig_final(rst_mh_2,'MH-2-70')

 save(
	    list = c("rst_yq_1","rst_yq_2","rst_yq_1f","rst_yq_2f","rst_jq_1","rst_jq_2","rst_jq_1f","rst_jq_2f","rst_mh_1",
      "rst_mh_2","rst_mh_1f","rst_mh_2f"),
	    file = paste("out_f_dxal",'-',mod_type,'.RData', sep = '')
	  )

gen_fig_final(rst_yq_1f,'YQ-1-70f')
gen_fig_final(rst_yq_2f,'YQ-2-70f')

gen_fig_final(rst_jq_1f,'JQ-1-70f')
gen_fig_final(rst_jq_2f,'JQ-2-70f')

gen_fig_final(rst_mh_1f,'MH-1-70f')
gen_fig_final(rst_mh_2f,'MH-2-70f')

setwd(path0)
