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
lib_nm1<-'zhang2013am.dll'
lib_nm2<-'mull2014.dll'
mod_type='M14'
niter_1=15000;
niter_2=15000;
burnin.length=1500;
maxit.sce0=5000;
maxit.sce1=5000;
lfactor=1000;
ufactor=1000;
lfactorf=1000;
ufactorf=1000;
Lb3<-c('NH4','NO3','Urea')
out_dir=".//out_pre_SUSF-CBS-90A-M14"
#path0<-"/home/junhui/research/N15CS/run"
path0<-"D://workdir//N15CPM//run"
setwd(path0)

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

obs_ao<-read.csv('.//input//nh4no3_cbssuc_stgavg90A.csv')
ini_val<-read.csv('.//input//Ini_value_cbssuc_stgavg90A.csv')

lvlid<-as.character(levels(as.factor(obs_ao$ID)))
nlvlid<-nlevels(as.factor(obs_ao$ID))

obs_al<-obs_ao[,which(names(obs_ao) %in% c('ID','lable','time','NH4', 'NO3', 'NH4_15', 'NO3_15','N2O','N2O_15'))]
set.seed(1000);obs_al<-missForest(obs_al, maxiter = 100, ntree = 2000)$ximp

obs_al$NH4_15<-obs_al$NH4*obs_al$NH4_15/100;
obs_al$NO3_15<-obs_al$NO3*obs_al$NO3_15/100;

obs_al$ID<-as.character(obs_al$ID)
ini_val$ID <- as.character(ini_val$ID)

rndgrp<-floor(runif(8,-20000,20000))
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

pars0.M14=c(
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
   K_NH4_2_NO2nit = 0.044,
   K_NO2nit_2_NO3 = 5.26,   
   K_NO3_2_NO3sto = 0.0005222072,
   K_NO3sto_2_NO3 = 0.004020091,
   K_NO3_2_NO2den = 0.000057,
   K_NRec_2_NO2org = 0.0000003,
   K_NO2org_2_N2Oorg = 0.00001,
   K_NO2nit_2_N2Onit = 0.18,
   K_NO2den_2_N2Oden = 0.021,
   K_NO2denNRec_2_N2Ocod = 0.0000006,
   K_N2Oorg_2_N2 = 0.000001,
   K_N2Ocod_2_N2 = 0.00006,
   K_N2Oden_2_N2 = 0.20,
   K_N2Onit_2_N2 = 0.18,
   K_N2Oorg_2_N2OorgE = 0.0036,
   K_N2Ocod_2_N2OcodE = 0.033,
   K_N2Oden_2_N2OdenE = 0.0048,
   K_N2Onit_2_N2OnitE = 0.00021
 )
 
 pars0.z13am=c(
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
 
loadcmp("m14_batch_CS1.Rc")
 

flst<-dir(pattern='.RData',path=out_dir,full.names=T)

load(flst[1])
pars.M14<-rbind(pars1.M14.nh4no3,pars2.M14.nh4no3,pars2.M14.no3,pars1.M14.nh4)
pars.Z13am<-rbind(MCMC0.NH4$bestpar,MCMC0.NO3$bestpar)
pars0.M14<-apply(pars.M14,2,median)
pars0.z13am<-apply(pars.Z13am,2,median)

niter_1=30000;
niter_2=30000;
burnin.length=5000;
lfactor=200;
ufactor=200;
lfactorf=500;
ufactorf=200;

out_dir='.//out_SUSF-CBS-90A-M14'

if(file.exists(out_dir)){file.rename(out_dir,paste(out_dir,substr(file.info(out_dir)$ctime,1,10),sep='_'))};
dir.create(out_dir)

rndgrp<-floor(runif(8,-20000,20000))
#rndgrp<-floor(seq(-10000,10000,length.out = 5))
grpPal<-merge(lvlid,rndgrp)
names(grpPal)<-c("ID","rndseed")
grpPal$ID<-as.character(grpPal$ID)
lgrpPal<-nrow(grpPal)

no_cores <- detectCores() - 0
ncl<-max(min(no_cores,lgrpPal),2)

loadcmp("m14_batch_CS1.Rc")

###############################################################
niter_1=50000;
niter_2=50000;
burnin.length=5000;
lfactor=2;
ufactor=2;

flst<-dir(pattern="*.RData",path=out_dir,full.names=T)
path0<-"D://N15CSg//run"
setwd(path0)

flst<-dir(pattern="*.RData",path=out_dir,full.names=T)

ii_1<-grep('SUCC-CBS-1',flst)
ii_3<-grep('SUCC-CBS-3',flst)
ii_4<-grep('SUCC-CBS-4',flst)
ii_5<-grep('SUCC-CBS-5',flst)
ii_6<-grep('SUCC-CBS-6',flst)

rst_1<-findBST(flst[ii_1])
rst_3<-findBST(flst[ii_3])
rst_4<-findBST(flst[ii_4])
rst_5<-findBST(flst[ii_5])
rst_6<-findBST(flst[ii_6])

rst_1_60f<-getpar_fnal_z13amg(arst=rst_1,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_3_60f<-getpar_fnal_z13amg(arst=rst_3,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_4_60f<-getpar_fnal_z13amg(arst=rst_4,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_5_60f<-getpar_fnal_z13amg(arst=rst_5,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)
rst_6_60f<-getpar_fnal_z13amg(arst=rst_6,fntrt_p=3,lib_nm,lfactor=lfactor,ufactor=ufactor,niter_1=niter_1,niter_2=niter_2,burninlength=burninlength)

out_dir_f='.//out_f_SUSF-CBS-60A-ZHANG2013amg'

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
