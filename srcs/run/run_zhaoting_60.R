#run all for Zhaoting soil water 90%

library(compiler)
library(doParallel)
library(doBy)
library(FME)
#library(hydromad)
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
lib_nm<-'zhang2013am.dll'
mod_type='Z13am'
niter_1=2000;
niter_2=2000;
burnin.length=500;
maxit.sce0=5000;
maxit.sce1=5000;
lfactor=1000;
ufactor=1000;
lfactorf=1000;
ufactorf=1000;
Lb3<-c('NH4','NO3','Urea')
out_dir='.//out_pre_SUSF-CBS-60-ZHANG2013am'
#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://N15CSg//run"
setwd(path0)

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

#obs_ao<-read.csv('.//input//nh4no3_fcbs.csv')
#ini_val<-read.csv('.//input//Ini_value_fcbs.csv')

#obs_ao<-read.csv('.//input//nh4no3_nbfb100.csv')
#ini_val<-read.csv('.//input//ini_value_nbfb100.csv')

# obs_ao<-read.csv('.//input//nh4no3_cbssuc_90raw.csv')
# ini_val<-read.csv('.//input//Ini_value_cbssuc_90raw.csv')

obs_ao<-read.csv('.//input//nh4no3_cbssuc_60raw.csv')
ini_val<-read.csv('.//input//Ini_value_cbssuc_60raw.csv')

#obs_ao<-read.csv('.//input//nh4no3_cbssuc_stgavg.csv')
#ini_val<-read.csv('.//input//Ini_value_cbssuc_stgavg.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15'))]
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


loadcmp("subs_all.Rc")
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
   K_NO3sto_2_NO3 = 0.004020091
 )
 
loadcmp("Z13am_batch_CS1.Rc")
 
flst<-dir(pattern='.RData',path=out_dir,full.names=T)

load(flst[1])
pars.Z13am<-rbind(MCMC$bestpar,MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar)

for (i in 2:length(flst)){
load(flst[i])
if((min(sign(diag(cor_ss_nh4_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4_obs)))))) pars.Z13am<-rbind(pars.Z13am,MCMC.NH4$bestpar)
if((min(sign(diag(cor_ss_no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_no3_obs)))))) pars.Z13am<-rbind(pars.Z13am,MCMC.NO3$bestpar)
if((min(sign(diag(cor_ss_nh4no3_obs)),na.rm=T)>0)&(!is.na(min(sign(diag(cor_ss_nh4no3_obs)))))) pars.Z13am<-rbind(pars.Z13am,MCMC.NH4NO3$bestpar)

if((min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))),na.rm=T)>0) & 
(!is.na(min(sign(c(diag(cor_s4_nh4_obs),diag(cor_s4_nh4no3_obs),diag(cor_s4_no3_obs))))))) pars.Z13am<-rbind(pars.Z13am,MCMC$bestpar)
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

out_dir='.//out_SUSF-CBS-60-ZHANG2013am'

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

loadcmp("Z13am_batch_CS1.Rc")

#################################################################
niter_1=100000;
niter_2=100000;
burnin.length=25000;
lfactor=2;
ufactor=2;

path0<-"D://N15CSg//run"
setwd(path0)

flst<-dir(pattern="*.RData",path=out_dir,full.names=T)

ii_11<-grep('-1_1-',flst)
ii_12<-grep('-1_2-',flst)
ii_13<-grep('-1_3-',flst)

ii_31<-grep('-3_1-',flst)
ii_32<-grep('-3_2-',flst)
ii_33<-grep('-3_3-',flst)

ii_41<-grep('-4_1-',flst)
ii_42<-grep('-4_2-',flst)
ii_43<-grep('-4_3-',flst)

ii_51<-grep('-5_1-',flst)
ii_52<-grep('-5_2-',flst)
ii_53<-grep('-5_3-',flst)

ii_61<-grep('-6_1-',flst)
ii_62<-grep('-6_1-',flst)
ii_63<-grep('-6_1-',flst)

rst_11<-findBST(flst[ii_11])
rst_12<-findBST(flst[ii_12])
rst_13<-findBST(flst[ii_13])

rst_31<-findBST(flst[ii_31])
rst_32<-findBST(flst[ii_32])
rst_33<-findBST(flst[ii_33])

rst_41<-findBST(flst[ii_41])
rst_42<-findBST(flst[ii_42])
rst_43<-findBST(flst[ii_43])

rst_51<-findBST(flst[ii_51])
rst_52<-findBST(flst[ii_52])
rst_53<-findBST(flst[ii_53])

rst_61<-findBST(flst[ii_61])
rst_62<-findBST(flst[ii_62])
rst_63<-findBST(flst[ii_63])


rst_11_60f<-getpar_fnal_z13am(arst=rst_11,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_12_60f<-getpar_fnal_z13am(arst=rst_12,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_13_60f<-getpar_fnal_z13am(arst=rst_13,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

rst_31_60f<-getpar_fnal_z13am(arst=rst_31,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_32_60f<-getpar_fnal_z13am(arst=rst_32,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_33_60f<-getpar_fnal_z13am(arst=rst_33,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

rst_41_60f<-getpar_fnal_z13am(arst=rst_41,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_42_60f<-getpar_fnal_z13am(arst=rst_42,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_43_60f<-getpar_fnal_z13am(arst=rst_43,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

rst_51_60f<-getpar_fnal_z13am(arst=rst_51,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_52_60f<-getpar_fnal_z13am(arst=rst_52,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_53_60f<-getpar_fnal_z13am(arst=rst_53,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

rst_61_60f<-getpar_fnal_z13am(arst=rst_61,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_62_60f<-getpar_fnal_z13am(arst=rst_62,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_63_60f<-getpar_fnal_z13am(arst=rst_63,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)


out_dir_f='.//out_f_SUSF-CBS-60-ZHANG2013amg'

if(file.exists(out_dir_f)){file.rename(out_dir_f,paste(out_dir_f,substr(file.info(out_dir_f)$ctime,1,10),sep='_'))};
dir.create(out_dir_f)

setwd(out_dir_f)

gen_fig(rst_11_60f,'SUSF-1-1-60')
gen_fig(rst_12_60f,'SUSF-1-2-60')
gen_fig(rst_13_60f,'SUSF-1-3-60')

gen_fig(rst_31_60f,'SUSF-3-1-60')
gen_fig(rst_32_60f,'SUSF-3-2-60')
gen_fig(rst_33_60f,'SUSF-3-3-60')

gen_fig(rst_41_60f,'SUSF-4-1-60')
gen_fig(rst_42_60f,'SUSF-4-2-60')
gen_fig(rst_43_60f,'SUSF-4-3-60')

gen_fig(rst_51_60f,'SUSF-5-1-60')
gen_fig(rst_52_60f,'SUSF-5-2-60')
gen_fig(rst_53_60f,'SUSF-5-3-60')

gen_fig(rst_61_60f,'SUSF-6-1-60')
gen_fig(rst_62_60f,'SUSF-6-2-60')
gen_fig(rst_63_60f,'SUSF-6-3-60')

 save(
	    list = c("rst_11_60f","rst_12_60f","rst_13_60f","rst_31_60f","rst_32_60f","rst_33_60f","rst_41_60f","rst_42_60f","rst_43_60f",
      "rst_51_60f","rst_52_60f","rst_53_60f","rst_61_60f","rst_62_60f","rst_63_60f"),
	    file = paste("out_f_SUSF-CBS-60",'-',mod_type,'.RData', sep = '')
	  )

gen_fig_final(rst_11_60f,'SUSF-1-1-60f')
gen_fig_final(rst_12_60f,'SUSF-1-2-60f')
gen_fig_final(rst_13_60f,'SUSF-1-3-60f')

gen_fig_final(rst_31_60f,'SUSF-3-1-60f')
gen_fig_final(rst_32_60f,'SUSF-3-2-60f')
gen_fig_final(rst_33_60f,'SUSF-3-3-60f')

gen_fig_final(rst_41_60f,'SUSF-4-1-60f')
gen_fig_final(rst_42_60f,'SUSF-4-2-60f')
gen_fig_final(rst_43_60f,'SUSF-4-3-60f')

gen_fig_final(rst_51_60f,'SUSF-5-1-60f')
gen_fig_final(rst_52_60f,'SUSF-5-2-60f')
gen_fig_final(rst_53_60f,'SUSF-5-3-60f')

gen_fig_final(rst_61_60f,'SUSF-6-1-60f')
gen_fig_final(rst_62_60f,'SUSF-6-2-60f')
gen_fig_final(rst_63_60f,'SUSF-6-3-60f')

setwd(path0)
############################################################### AVG

obs_ao<-read.csv('.//input//nh4no3_cbssuc_stgavg60.csv')
ini_val<-read.csv('.//input//Ini_value_cbssuc_stgavg60.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15'))]
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

loadcmp("Z13am_batch_CS1.Rc")

###############################################################
niter_1=100000;
niter_2=100000;
burnin.length=25000;
lfactor=2;
ufactor=2;

path0<-"D://N15CSg//run"
setwd(path0)

flst<-dir(pattern="*.RData",path=out_dir,full.names=T)

ii_1<-grep('IDT-1',flst)
ii_3<-grep('IDT-3',flst)
ii_4<-grep('IDT-4',flst)
ii_5<-grep('IDT-5',flst)
ii_6<-grep('IDT-6',flst)

rst_1<-findBST(flst[ii_1])
rst_3<-findBST(flst[ii_3])
rst_4<-findBST(flst[ii_4])
rst_5<-findBST(flst[ii_5])
rst_6<-findBST(flst[ii_6])


rst_1_60f<-getpar_fnal_z13am(arst=rst_1,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_3_60f<-getpar_fnal_z13am(arst=rst_3,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_4_60f<-getpar_fnal_z13am(arst=rst_4,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_5_60f<-getpar_fnal_z13am(arst=rst_5,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_6_60f<-getpar_fnal_z13am(arst=rst_6,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

out_dir_f='.//out_f_SUSF-CBS-60AVG-ZHANG2013amg'

if(file.exists(out_dir_f)){file.rename(out_dir_f,paste(out_dir_f,substr(file.info(out_dir_f)$ctime,1,10),sep='_'))};
dir.create(out_dir_f)

setwd(out_dir_f)

gen_fig(rst_1,'SUSF-1-60AVG')
gen_fig(rst_3,'SUSF-3-60AVG')
gen_fig(rst_4,'SUSF-4-60AVG')
gen_fig(rst_5,'SUSF-5-60AVG')
gen_fig(rst_6,'SUSF-6-60AVG')

 save(
	    list = c("rst_1_60f","rst_3_60f","rst_4_60f","rst_5_60f","rst_6_60f"),
	    file = paste("out_f_SUSF-CBS-60AVG",'-',mod_type,'.RData', sep = '')
	  )

gen_fig_final(rst_1_60f,'SUSF-1-60AVGf')
gen_fig_final(rst_3_60f,'SUSF-3-60AVGf')
gen_fig_final(rst_4_60f,'SUSF-4-60AVGf')
gen_fig_final(rst_5_60f,'SUSF-5-60AVGf')
gen_fig_final(rst_6_60f,'SUSF-6-60AVGf')

setwd(path0)
