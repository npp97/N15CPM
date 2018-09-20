
cl<-makeCluster(ncl)
registerDoParallel(cl)

pars<-foreach(i=1:lgrpPal,.packages=c('FME','compiler','doBy'),.combine=rbind) %dopar%
{
  rngjk=grpPal$rndseed[i]
  set.seed(rngjk)  
  dyn.load(lib_nm1)
  dyn.load(lib_nm2)
# for (i in 1:nlvlid){
   sink(paste(out_dir,'//rst_',grpPal$ID[i],"-CS-",rngjk,'-',mod_type,'.txt',sep=''))

    ii_nh4 <- which((obs_al$ID %in% grpPal$ID[i])&(obs_al$lable %in% Lb3[1]))
    jk.nh4 <- which((ini_val$ID %in% grpPal$ID[i])&(ini_val$lable %in% Lb3[1]))  
    obs_nh4 <- obs_al[ii_nh4, 3:9]
    
    print(paste("-----------", grpPal$ID[i], "-----------"))
    
    #------------------Inti STATE VARIABLES
    y0.nh4 = c(NLab = ini_val$iNLab[jk.nh4],
               NLab_15 = ini_val$iNLab_15[jk.nh4]*0.01*ini_val$iNLab[jk.nh4],
               NH4 = ini_val$iNH4[jk.nh4],
               NH4_15 = ini_val$iNH4_15[jk.nh4]*0.01*ini_val$iNH4[jk.nh4],
               NO3 =ini_val$iNO3[jk.nh4],
               NO3_15 =ini_val$iNO3_15[jk.nh4]*0.01*ini_val$iNO3[jk.nh4],
               NH4ads = ini_val$iNH4ads[jk.nh4],
               NH4ads_15 = ini_val$iNH4ads_15[jk.nh4]*0.01*ini_val$iNH4ads[jk.nh4],
               NRec = ini_val$iNRec[jk.nh4],
               NRec_15 = ini_val$iNRec_15[jk.nh4]*0.01*ini_val$iNRec[jk.nh4],
               NO2 = ini_val$NO2[jk.nh4],
               NO2_15 = ini_val$NO2[jk.nh4]*0.01*ini_val$NO2_15[jk.nh4],
               NO2org =  ini_val$NO2[jk.nh4]*0.07,
               NO2den =  ini_val$NO2[jk.nh4]*0.50,
               NO2nit =  ini_val$NO2[jk.nh4]*0.43,
               NO2org_15 = ini_val$NO2[jk.nh4]*0.01*ini_val$NO2_15[jk.nh4]*0.07,
               NO2den_15 = ini_val$NO2[jk.nh4]*0.01*ini_val$NO2_15[jk.nh4]*0.50,
               NO2nit_15 = ini_val$NO2[jk.nh4]*0.01*ini_val$NO2_15[jk.nh4]*0.43,
               N2O = 0,
               N2O_15 = 0,
               N2Oorg = 0,
               N2Ocod = 0,
               N2Onit = 0,
               N2Oden = 0,
               N2Oorg_15 = 0,
               N2Ocod_15 = 0,
               N2Onit_15 = 0,
               N2Oden_15 = 0,
               EN2Oorg = 0,
               EN2Ocod = 0,
               EN2Onit = 0,
               EN2Oden = 0,
               EN2Oorg_15 = 0,
               EN2Ocod_15 = 0,
               EN2Onit_15 = 0,
               EN2Oden_15 = 0,
               N2 = 0,
               N2_15 = 0
     );
    nh40.nh4<-coef(lm(NH4~time,data=obs_nh4[1:2,]))[1]
    no30.nh4<-coef(lm(NO3~time,data=obs_nh4[1:2,]))[1]
    nh4_15_0.nh4<-coef(lm(NH4_15~time,data=obs_nh4[1:2,]))[1]
    no3_15_0.nh4<-coef(lm(NO3_15~time,data=obs_nh4[1:2,]))[1]
    n2o0.nh4<-coef(lm(N2O~time,data=obs_nh4[1:2,]))[1]
    n2o_15_0.nh4<-coef(lm(N2O_15~time,data=obs_nh4[1:2,]))[1]
    
    #nh4_14_0.nh4<-coef(lm(NH4_14~time,data=obs_nh4[1:2,]))[1]
    #no3_14_0.nh4<-coef(lm(NO3_14~time,data=obs_nh4[1:2,]))[1]
    
    obs_nh40<-c(time=0,NH4=nh40.nh4,NO3=no30.nh4,NH4_15=nh4_15_0.nh4,NO3_15=no3_15_0.nh4,N2O=n2o0.nh4,N2O_15=n2o_15_0.nh4);#,NH4_14=nh4_14_0.nh4,NO3_14=no3_14_0.nh4);
    obs_nh4<-rbind(obs_nh40,obs_nh4);
	nn.nh4<-get_n2oc15(time0=obs_nh4$time,n2o=obs_nh4$N2O,n2o_15=obs_nh4$N2O_15)
	obs_nh4[,'N2O']=nn.nh4[,'N2O'];obs_nh4[,'N2O_15']=nn.nh4[,'N2O_15'];
	

    y0.nh4['NH4']<-obs_nh4[1,'NH4']
    y0.nh4['NO3']<-obs_nh4[1,'NO3']; 
    y0.nh4['NH4_15']<-obs_nh4[1,'NH4_15'];
    y0.nh4['NO3_15']<-obs_nh4[1,'NO3_15'];
    y0.nh4['N2O']<-obs_nh4[1,"N2O"];
    y0.nh4['N2O_15']<-obs_nh4[1,"N2O_15"];
    
    
    ii_no3 <- which((obs_al$ID %in% grpPal$ID[i])&(obs_al$lable %in% Lb3[2]))
    jk.no3 <- which((ini_val$ID %in% grpPal$ID[i])&(ini_val$lable %in% Lb3[2]))  
    obs_no3 <- obs_al[ii_no3, 3:9]
    
    #------------------Inti STATE VARIABLES
    y0.no3 = c(NLab = ini_val$iNLab[jk.no3],
               NLab_15 = ini_val$iNLab_15[jk.no3]*0.01*ini_val$iNLab[jk.no3],
               NH4 = ini_val$iNH4[jk.no3],
               NH4_15 = ini_val$iNH4_15[jk.no3]*0.01*ini_val$iNH4[jk.no3],
               NO3 =ini_val$iNO3[jk.no3],
               NO3_15 =ini_val$iNO3_15[jk.no3]*0.01*ini_val$iNO3[jk.no3],
               NH4ads = ini_val$iNH4ads[jk.no3],
               NH4ads_15 = ini_val$iNH4ads_15[jk.no3]*0.01*ini_val$iNH4ads[jk.no3],
               NRec = ini_val$iNRec[jk.no3],
               NRec_15 = ini_val$iNRec_15[jk.no3]*0.01*ini_val$iNRec[jk.no3],
               NO2 = ini_val$NO2[jk.no3],
               NO2_15 = ini_val$NO2[jk.no3]*0.01*ini_val$NO2_15[jk.no3],
               NO2org =  ini_val$NO2[jk.no3]*0.07,
               NO2den =  ini_val$NO2[jk.no3]*0.50,
               NO2nit =  ini_val$NO2[jk.no3]*0.43,
               NO2org_15 = ini_val$NO2[jk.no3]*0.01*ini_val$NO2_15[jk.no3]*0.07,
               NO2den_15 = ini_val$NO2[jk.no3]*0.01*ini_val$NO2_15[jk.no3]*0.50,
               NO2nit_15 = ini_val$NO2[jk.no3]*0.01*ini_val$NO2_15[jk.no3]*0.43,
               N2O = 0,
               N2O_15 = 0,
               N2Oorg = 0,
               N2Ocod = 0,
               N2Onit = 0,
               N2Oden = 0,
               N2Oorg_15 = 0,
               N2Ocod_15 = 0,
               N2Onit_15 = 0,
               N2Oden_15 = 0,
               EN2Oorg = 0,
               EN2Ocod = 0,
               EN2Onit = 0,
               EN2Oden = 0,
               EN2Oorg_15 = 0,
               EN2Ocod_15 = 0,
               EN2Onit_15 = 0,
               EN2Oden_15 = 0,
               N2 = 0,
               N2_15 = 0
    );
	
    nh40.no3<-coef(lm(NH4~time,data=obs_no3[1:2,]))[1]
    no30.no3<-coef(lm(NO3~time,data=obs_no3[1:2,]))[1]
    nh4_15_0.no3<-coef(lm(NH4_15~time,data=obs_no3[1:2,]))[1]
    no3_15_0.no3<-coef(lm(NO3_15~time,data=obs_no3[1:2,]))[1]
    n2o0.no3<-coef(lm(N2O~time,data=obs_no3[1:2,]))[1]
    n2o_15_0.no3<-coef(lm(N2O_15~time,data=obs_no3[1:2,]))[1]
    
    obs_no30<-c(time=0,NH4=nh40.no3,NO3=no30.no3,NH4_15=nh4_15_0.no3,NO3_15=no3_15_0.no3,N2O=n2o0.no3,N2O_15=n2o_15_0.no3);#,NH4_14=nh4_14_0.no3,NO3_14=no3_14_0.no3)
    obs_no3<-rbind(obs_no30,obs_no3);
	nn.no3<-get_n2oc15(time0=obs_no3$time,n2o=obs_no3$N2O,n2o_15=obs_no3$N2O_15)
	obs_no3[,'N2O']=nn.no3[,'N2O'];obs_no3[,'N2O_15']=nn.no3[,'N2O_15'];
	
    y0.no3['NH4']<-obs_no3[1,'NH4']
    y0.no3['NO3']<-obs_no3[1,'NO3']; 
    y0.no3['NH4_15']<-obs_no3[1,'NH4_15'];
    y0.no3['NO3_15']<-obs_no3[1,'NO3_15'];
    y0.no3['N2O']<-obs_no3[1,"N2O"];
	y0.no3['N2O_15']<-obs_no3[1,'N2O_15'];
    
    ii_nh4no3 <- which((obs_al$ID %in% grpPal$ID[i])&(obs_al$lable %in% Lb3[3]))
    jk.nh4no3 <- which((ini_val$ID %in% grpPal$ID[i])&(ini_val$lable %in% Lb3[3]))  
    obs_nh4no3 <- obs_al[ii_nh4no3, 3:9]
    
    #------------------Inti STATE VARIABLES
    y0.nh4no3 = c(NLab = ini_val$iNLab[jk.nh4no3],
                  NLab_15 = ini_val$iNLab_15[jk.nh4no3]*0.01*ini_val$iNLab[jk.nh4no3],
                  NH4 = ini_val$iNH4[jk.nh4no3],
                  NH4_15 = ini_val$iNH4_15[jk.nh4no3]*0.01*ini_val$iNH4[jk.nh4no3],
                  NO3 =ini_val$iNO3[jk.nh4no3],
                  NO3_15 =ini_val$iNO3_15[jk.nh4no3]*0.01*ini_val$iNO3[jk.nh4no3],
                  NH4ads = ini_val$iNH4ads[jk.nh4no3],
                  NH4ads_15 = ini_val$iNH4ads_15[jk.nh4no3]*0.01*ini_val$iNH4ads[jk.nh4no3],
                  NRec = ini_val$iNRec[jk.nh4no3],
                  NRec_15 = ini_val$iNRec_15[jk.nh4no3]*0.01*ini_val$iNRec[jk.nh4no3],
                  NO2 = ini_val$NO2[jk.nh4no3],
                  NO2_15 = ini_val$NO2[jk.nh4no3]*0.01*ini_val$NO2_15[jk.nh4no3],
                  NO2org =  ini_val$NO2[jk.nh4no3]*0.07,
                  NO2den =  ini_val$NO2[jk.nh4no3]*0.50,
                  NO2nit =  ini_val$NO2[jk.nh4no3]*0.43,
                  NO2org_15 = ini_val$NO2[jk.nh4no3]*0.01*ini_val$NO2_15[jk.nh4no3]*0.07,
                  NO2den_15 = ini_val$NO2[jk.nh4no3]*0.01*ini_val$NO2_15[jk.nh4no3]*0.50,
                  NO2nit_15 = ini_val$NO2[jk.nh4no3]*0.01*ini_val$NO2_15[jk.nh4no3]*0.43,
                  N2O = 0,
                  N2O_15 = 0,
                  N2Oorg = 0,
                  N2Ocod = 0,
                  N2Onit = 0,
                  N2Oden = 0,
                  N2Oorg_15 = 0,
                  N2Ocod_15 = 0,
                  N2Onit_15 = 0,
                  N2Oden_15 = 0,
                  EN2Oorg = 0,
                  EN2Ocod = 0,
                  EN2Onit = 0,
                  EN2Oden = 0,
                  EN2Oorg_15 = 0,
                  EN2Ocod_15 = 0,
                  EN2Onit_15 = 0,
                  EN2Oden_15 = 0,
                  N2 = 0,
                  N2_15 = 0              
     );
    
    nh40.nh4no3<-coef(lm(NH4~time,data=obs_nh4no3[1:2,]))[1]
    no30.nh4no3<-coef(lm(NO3~time,data=obs_nh4no3[1:2,]))[1]
    nh4_15_0.nh4no3<-coef(lm(NH4_15~time,data=obs_nh4no3[1:2,]))[1]
    no3_15_0.nh4no3<-coef(lm(NO3_15~time,data=obs_nh4no3[1:2,]))[1]
    n2o0.nh4no3<-coef(lm(N2O~time,data=obs_nh4no3[1:2,]))[1]
    n2o_15_0.nh4no3<-coef(lm(N2O_15~time,data=obs_nh4no3[1:2,]))[1]
    
    #nh4_14_0.nh4no3<-coef(lm(NH4_14~time,data=obs_nh4no3[1:2,]))[1]
    #no3_14_0.nh4no3<-coef(lm(NO3_14~time,data=obs_nh4no3[1:2,]))[1]
    
    obs_nh4no30<-c(time=0,NH4=nh40.nh4no3,NO3=no30.nh4no3,NH4_15=nh4_15_0.nh4no3,NO3_15=no3_15_0.nh4no3,N2O=n2o0.nh4no3,N2O_15=n2o_15_0.nh4no3);#,NH4_14=nh4_14_0.nh4no3,NO3_14=no3_14_0.nh4no3);
    obs_nh4no3<-rbind(obs_nh4no30,obs_nh4no3);
	nn.nh4no3<-get_n2oc15(time0=obs_nh4no3$time,n2o=obs_nh4no3$N2O,n2o_15=obs_nh4no3$N2O_15)
	obs_nh4no3[,'N2O']=nn.nh4no3[,'N2O'];obs_nh4no3[,'N2O_15']=nn.nh4no3[,'N2O_15'];
	
	y0.nh4no3['NH4']<-obs_nh4no3[1,'NH4']
    y0.nh4no3['NO3']<-obs_nh4no3[1,'NO3']; 
    y0.nh4no3['NH4_15']<-obs_nh4no3[1,'NH4_15'];
    y0.nh4no3['NO3_15']<-obs_nh4no3[1,'NO3_15'];
    y0.nh4no3['N2O']<-obs_nh4no3[1,"N2O"];
    y0.nh4no3['N2O_15']<-obs_nh4no3[1,'N2O_15'];

     times <- obs_nh4[, 'time']
    
    ### MCMC simulation
    #----------------------------------------------------------------------------
    upper<-lower<-pars0.m07s;
    lower1= pars0.m07s/lfactor;
    #lower[]<-0;
    upper1= pars0.m07s*ufactor;
    #upper=apply(data.frame(rep(1,length(pars0)),upper),1,min);
     #print("---------------NH4-N15")
    try({
    MCMC0.NH4<-try(mdMCMCM07sB(ff=M07sAcostL,pars=pars0.m07s,y0=y0.nh4[1:10],obs=obs_nh4[,1:5],niter1=niter_1*2,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
    #print(MCMC.NH4$bestpar)
    MCMC0.NO3<-try(mdMCMCM07sB(ff=M07sAcostL,pars=pars0.m07s,y0=y0.no3[1:10],obs=obs_no3[,1:5],niter1=niter_1*2,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
    })    

####
####
####
if((class(MCMC0.NH4) %in% "modMCMC")|(class(MCMC0.NO3) %in% "modMCMC")){

if((class(MCMC0.NH4) %in% "modMCMC")&(class(MCMC0.NO3) %in% "modMCMC"))
 {
	 #One by One 
##############
############## NH4-->NO3-->(1) M14 and --> (2) NH4NO3 --> M14
##############
	  pars0.m07s<-MCMC0.NH4$bestpar;
      pars.fix.nh4.m07s=c(
      MCMC0.NH4$bestpar["C_NRec_2_NH4"],
      MCMC0.NH4$bestpar["C_NH4_2_NRec"],
      MCMC0.NH4$bestpar["C_NH4_2_NLab"],
      MCMC0.NH4$bestpar["K_NLab_2_NH4"],
      MCMC0.NH4$bestpar["K_NH4ads_2_NH4"],
      MCMC0.NH4$bestpar["K_NH4_2_NH4ads"]
      )
	  
      ii.fix<-which(names(pars0.m07s) %in% names(pars.fix.nh4.m07s))
      ii.v<-which(!(names(pars0.m07s) %in% names(pars.fix.nh4.m07s)))

      lower1=pars0.m07s[ii.v]/lfactor;
      upper1=pars0.m07s[ii.v]*ufactor;
	  
      MCMC01.NO3<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars0.m07s[ii.v],pars0=pars0.m07s,pars_fix=pars0.m07s[ii.fix],obs=obs_no3[,1:5],niter1=niter_1,y0=y0.no3[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
      
	  pars1.m07s.no3<-pars0.m07s;
	  pars1.m07s.no3[ii.v]<-MCMC01.NO3$bestpar;
	   
	  lower1=pars1.m07s.no3/lfactor;
	  upper1=pars1.m07s.no3*ufactor;
	  
	  MCMC01.NH4NO3<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars1.m07s.no3,pars0=pars1.m07s.no3,pars_fix=pars1.m07s.no3,obs=obs_nh4no3[,1:5],niter1=niter_1,y0=y0.no3[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	  
      pars0.M14.no3<-pars0.M14
      pars0.M14.no3["C_NRec_2_NH4"] = as.numeric(pars1.m07s.no3["C_NRec_2_NH4"]);
      pars0.M14.no3["C_NH4_2_NRec"]  =  as.numeric(pars1.m07s.no3["C_NH4_2_NRec"]);
      pars0.M14.no3["C_NH4_2_NLab"] = as.numeric(pars1.m07s.no3["C_NH4_2_NLab"]);
      pars0.M14.no3["K_NLab_2_NH4"] = as.numeric(pars1.m07s.no3["K_NLab_2_NH4"]);
      pars0.M14.no3["K_NH4ads_2_NH4"]  = as.numeric(pars1.m07s.no3["K_NH4ads_2_NH4"]);
      pars0.M14.no3["K_NH4_2_NH4ads"] = as.numeric(pars1.m07s.no3["K_NH4_2_NH4ads"]);
      pars0.M14.no3["C_NO3_2_NRec"] = as.numeric(pars1.m07s.no3["C_NO3_2_NRec"]);
      pars0.M14.no3["C_NO3_2_NH4"] = as.numeric(pars1.m07s.no3["C_NO3_2_NH4"]);
      pars0.M14.no3["C_NRec_2_NO3"] = as.numeric(pars1.m07s.no3["C_NRec_2_NO3"]);
      
      ii.fix<-which(names(pars0.M14.no3) %in% names(pars1.m07s.no3))
      ii.v<-which(!(names(pars0.M14.no3) %in% names(pars1.m07s.no3)))
      
      lower1=pars0.M14.no3[ii.v]/lfactor;
      upper1=pars0.M14.no3[ii.v]*ufactor;
      
#	  lower1=pars0.M14/lfactor ##Tmp
#	  upper1=pars0.M14*ufactor ##
	  
  
	  MCMC4_2.NO3<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.no3[ii.v],pars0=pars0.M14.no3,pars_fix=pars0.M14.no3[ii.fix],obs=obs_nh4,niter1=niter_1,y0=y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
  
	  pars2.M14.no3<-pars0.M14.no3;
	  pars2.M14.no3[ii.v]<-MCMC4_2.NO3$bestpar;
	  
	  
	  #
	  pars0.M14.nh4no3<-pars0.M14
	  pars0.M14.nh4no3["C_NRec_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NH4"]);
	  pars0.M14.nh4no3["C_NH4_2_NRec"]  =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NRec"]);
	  pars0.M14.nh4no3["C_NH4_2_NLab"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NLab"]);
	  pars0.M14.nh4no3["K_NLab_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NLab_2_NH4"]);
	  pars0.M14.nh4no3["K_NH4ads_2_NH4"]  =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4ads_2_NH4"]);
	  pars0.M14.nh4no3["K_NH4_2_NH4ads"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4_2_NH4ads"]);
	  pars0.M14.nh4no3["C_NO3_2_NRec"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NRec"]);
	  pars0.M14.nh4no3["C_NO3_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NH4"]);
	  pars0.M14.nh4no3["C_NRec_2_NO3"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NO3"]);
	  
	  ii.fix<-which(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar))
	  ii.v<-which(!(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar)))
	  
	  lower1=pars0.M14.nh4no3[ii.v]/lfactor;
	  upper1=pars0.M14.nh4no3[ii.v]*ufactor;
	  
	  
	  MCMC4_2.NH4NO3<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.nh4no3[ii.v],pars0=pars0.M14.nh4no3,pars_fix=pars0.M14.nh4no3[ii.fix],obs=obs_nh4no3,niter1=niter_1,y0=y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	  
	  pars2.M14.nh4no3<-pars0.M14.no3;
	  pars2.M14.nh4no3[ii.v]<-MCMC4_2.NH4NO3$bestpar;
	  
###########
########### NO3-->NH4-->(1) M14 and --> (2) NH4NO3 --> M14
###########    
     pars0.m07s<-MCMC0.NO3$bestpar;
    
     pars.fix.no3.m07s<-c(
      MCMC0.NO3$bestpar["C_NO3_2_NRec"],
      MCMC0.NO3$bestpar["C_NRec_2_NO3"],
      MCMC0.NO3$bestpar["C_NO3_2_NH4"],
     );
	 
    ii.fix<-which(names(pars0.m07s) %in% names(pars.fix.no3.m07s))
    ii.v<-which(!(names(pars0.m07s) %in% names(pars.fix.no3.m07s)))

    lower1=pars0.m07s[ii.v]/lfactor;
    upper1=pars0.m07s[ii.v]*ufactor;
    
    MCMC01.NH4<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars0.m07s[ii.v],pars0=pars0.m07s,pars_fix=pars0.m07s[ii.fix],obs=obs_nh4[,1:5],niter1=niter_1,y0=y0.nh4[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));

	pars1.m07s.nh4<-pars0.m07s;
	pars1.m07s.nh4[ii.v]<-MCMC01.NH4$bestpar;
	
    lower1=pars1.m07s.nh4/lfactor;
    upper1=pars1.m07s.nh4*ufactor;
    
    MCMC02.NH4NO3<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars1.m07s.nh4,pars0=pars1.m07s.nh4,pars_fix=pars1.m07s.nh4,obs=obs_nh4no3[,1:5],niter1=niter_1,y0=y0.no3[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
#---------------------   
    pars0.M14.nh4<-pars0.M14
    pars0.M14.nh4["C_NRec_2_NH4"] =  as.numeric(pars1.m07s.nh4["C_NRec_2_NH4"]);
    pars0.M14.nh4["C_NH4_2_NRec"]  =  as.numeric(pars1.m07s.nh4["C_NH4_2_NRec"]);
    pars0.M14.nh4["C_NH4_2_NLab"] =  as.numeric(pars1.m07s.nh4["C_NH4_2_NLab"]);
    pars0.M14.nh4["K_NLab_2_NH4"] =  as.numeric(pars1.m07s.nh4["K_NLab_2_NH4"]);
    pars0.M14.nh4["K_NH4ads_2_NH4"]  =  as.numeric(pars1.m07s.nh4["K_NH4ads_2_NH4"]);
    pars0.M14.nh4["K_NH4_2_NH4ads"] =  as.numeric(pars1.m07s.nh4["K_NH4_2_NH4ads"]);
    pars0.M14.nh4["C_NO3_2_NRec"] =  as.numeric(pars1.m07s.nh4["C_NO3_2_NRec"]);
    pars0.M14.nh4["K_NO3_2_NH4"] =  as.numeric(pars1.m07s.nh4["K_NO3_2_NH4"]);
    pars0.M14.nh4["C_NRec_2_NO3"] =  as.numeric(pars1.m07s.nh4["C_NRec_2_NO3"]);

    ii.fix<-which(names(pars0.M14.nh4) %in% names(pars1.m07s.nh4))
    ii.v<-which(!(names(pars0.M14.nh4) %in% names(pars1.m07s.nh4)))
    
    lower1=pars0.M14.nh4[ii.v]/lfactor;
    upper1=pars0.M14.nh4[ii.v]*ufactor;
    
    MCMC4_1.NH4<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.nh4[ii.v],pars0=pars0.M14.nh4,pars_fix=pars0.M14.nh4[ii.fix],obs=obs_nh4,niter1=niter_1,y0=y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars1.M14.nh4<-pars0.M14.nh4;
	pars1.M14.nh4[ii.v]<-MCMC4_1.NH4$bestpar;
	
	#
	pars0.M14.nh4no3<-pars0.M14
	pars0.M14.nh4no3["C_NRec_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NH4"]);
	pars0.M14.nh4no3["C_NH4_2_NRec"]  =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NRec"]);
	pars0.M14.nh4no3["C_NH4_2_NLab"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NLab"]);
	pars0.M14.nh4no3["K_NLab_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NLab_2_NH4"]);
	pars0.M14.nh4no3["K_NH4ads_2_NH4"]  =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4ads_2_NH4"]);
	pars0.M14.nh4no3["K_NH4_2_NH4ads"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4_2_NH4ads"]);
	pars0.M14.nh4no3["C_NO3_2_NRec"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NRec"]);
	pars0.M14.nh4no3["C_NO3_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NH4"]);
	pars0.M14.nh4no3["C_NRec_2_NO3"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NO3"]);
	
	ii.fix<-which(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar))
	ii.v<-which(!(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar)))
	
	lower1=pars0.M14.nh4no3[ii.v]/lfactor;
	upper1=pars0.M14.nh4no3[ii.v]*ufactor;
	
	MCMC4_1.NH4NO3<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.nh4no3[ii.v],pars0=pars0.M14.nh4no3,pars_fix=pars0.M14.nh4no3[ii.fix],obs=obs_nh4no3,niter1=niter_1,y0=y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));

	pars1.M14.nh4no3<-pars0.M14.no3;
	pars1.M14.nh4no3[ii.v]<-MCMC4_1.NH4NO3$bestpar;
	
}
 #######################################################################################################################################################################################################
if((class(MCMC0.NH4) %in% "modMCMC")&(!(class(MCMC0.NO3) %in% "modMCMC")))
{
	pars0.m07s<-MCMC0.NH4$bestpar;
	pars.fix.nh4.m07s=c(
			MCMC0.NH4$bestpar["C_NRec_2_NH4"],
			MCMC0.NH4$bestpar["C_NH4_2_NRec"],
			MCMC0.NH4$bestpar["C_NH4_2_NLab"],
			MCMC0.NH4$bestpar["K_NLab_2_NH4"],
			MCMC0.NH4$bestpar["K_NH4ads_2_NH4"],
			MCMC0.NH4$bestpar["K_NH4_2_NH4ads"]
	)
	
	ii.fix<-which(names(pars0.m07s) %in% names(pars.fix.nh4.m07s))
	ii.v<-which(!(names(pars0.m07s) %in% names(pars.fix.nh4.m07s)))
	
	lower1=pars0.m07s[ii.v]/lfactor;
	upper1=pars0.m07s[ii.v]*ufactor;
	
	MCMC01.NO3<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars0.m07s[ii.v],pars0=pars0.m07s,pars_fix=pars0.m07s[ii.fix],obs=obs_no3[,1:5],niter1=niter_1,y0=y0.no3[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars1.m07s.no3<-pars0.m07s;
	pars1.m07s.no3[ii.v]<-MCMC01.NO3$bestpar;
	
	lower1=pars1.m07s.no3/lfactor;
	upper1=pars1.m07s.no3*ufactor;
	
	MCMC01.NH4NO3<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars1.m07s.no3,pars0=pars1.m07s.no3,pars_fix=pars1.m07s.no3,obs=obs_nh4no3[,1:5],niter1=niter_1,y0=y0.no3[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars0.M14.no3<-pars0.M14
	pars0.M14.no3["C_NRec_2_NH4"] = as.numeric(pars1.m07s.no3["C_NRec_2_NH4"]);
	pars0.M14.no3["C_NH4_2_NRec"]  =  as.numeric(pars1.m07s.no3["C_NH4_2_NRec"]);
	pars0.M14.no3["C_NH4_2_NLab"] = as.numeric(pars1.m07s.no3["C_NH4_2_NLab"]);
	pars0.M14.no3["K_NLab_2_NH4"] = as.numeric(pars1.m07s.no3["K_NLab_2_NH4"]);
	pars0.M14.no3["K_NH4ads_2_NH4"]  = as.numeric(pars1.m07s.no3["K_NH4ads_2_NH4"]);
	pars0.M14.no3["K_NH4_2_NH4ads"] = as.numeric(pars1.m07s.no3["K_NH4_2_NH4ads"]);
	pars0.M14.no3["C_NO3_2_NRec"] = as.numeric(pars1.m07s.no3["C_NO3_2_NRec"]);
	pars0.M14.no3["C_NO3_2_NH4"] = as.numeric(pars1.m07s.no3["C_NO3_2_NH4"]);
	pars0.M14.no3["C_NRec_2_NO3"] = as.numeric(pars1.m07s.no3["C_NRec_2_NO3"]);
	
	ii.fix<-which(names(pars0.M14.no3) %in% names(pars1.m07s.no3))
	ii.v<-which(!(names(pars0.M14.no3) %in% names(pars1.m07s.no3)))
	
	lower1=pars0.M14.no3[ii.v]/lfactor;
	upper1=pars0.M14.no3[ii.v]*ufactor;
	
#	  lower1=pars0.M14/lfactor ##Tmp
#	  upper1=pars0.M14*ufactor ##
	
	
	MCMC4_2.NO3<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.no3[ii.v],pars0=pars0.M14.no3,pars_fix=pars0.M14.no3[ii.fix],obs=obs_nh4,niter1=niter_1,y0=y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars2.M14.no3<-pars0.M14.no3;
	pars2.M14.no3[ii.v]<-MCMC4_2.NO3$bestpar;
	
	
	#
	pars0.M14.nh4no3<-pars0.M14
	pars0.M14.nh4no3["C_NRec_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NH4"]);
	pars0.M14.nh4no3["C_NH4_2_NRec"]  =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NRec"]);
	pars0.M14.nh4no3["C_NH4_2_NLab"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NLab"]);
	pars0.M14.nh4no3["K_NLab_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NLab_2_NH4"]);
	pars0.M14.nh4no3["K_NH4ads_2_NH4"]  =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4ads_2_NH4"]);
	pars0.M14.nh4no3["K_NH4_2_NH4ads"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4_2_NH4ads"]);
	pars0.M14.nh4no3["C_NO3_2_NRec"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NRec"]);
	pars0.M14.nh4no3["C_NO3_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NH4"]);
	pars0.M14.nh4no3["C_NRec_2_NO3"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NO3"]);
	
	ii.fix<-which(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar))
	ii.v<-which(!(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar)))
	
	lower1=pars0.M14.nh4no3[ii.v]/lfactor;
	upper1=pars0.M14.nh4no3[ii.v]*ufactor;
	
	
	MCMC4_2.NH4NO3<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.nh4no3[ii.v],pars0=pars0.M14.nh4no3,pars_fix=pars0.M14.nh4no3[ii.fix],obs=obs_nh4no3,niter1=niter_1,y0=y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars2.M14.nh4no3<-pars0.M14.no3;
	pars2.M14.nh4no3[ii.v]<-MCMC4_2.NH4NO3$bestpar;
	
}
#################################################################################################################################################################################	
if((!(class(MCMC0.NH4) %in% "modMCMC"))&(class(MCMC0.NO3) %in% "modMCMC"))
{
	pars0.m07s<-MCMC0.NO3$bestpar;
	
	pars.fix.no3.m07s<-c(
			MCMC0.NO3$bestpar["C_NO3_2_NRec"],
			MCMC0.NO3$bestpar["C_NRec_2_NO3"],
			MCMC0.NO3$bestpar["C_NO3_2_NH4"],
	);
	
	ii.fix<-which(names(pars0.m07s) %in% names(pars.fix.no3.m07s))
	ii.v<-which(!(names(pars0.m07s) %in% names(pars.fix.no3.m07s)))
	
	lower1=pars0.m07s[ii.v]/lfactor;
	upper1=pars0.m07s[ii.v]*ufactor;
	
	MCMC01.NH4<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars0.m07s[ii.v],pars0=pars0.m07s,pars_fix=pars0.m07s[ii.fix],obs=obs_nh4[,1:5],niter1=niter_1,y0=y0.nh4[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars1.m07s.nh4<-pars0.m07s;
	pars1.m07s.nh4[ii.v]<-MCMC01.NH4$bestpar;
	
	lower1=pars1.m07s.nh4/lfactor;
	upper1=pars1.m07s.nh4*ufactor;
	
	MCMC02.NH4NO3<-try(mdMCMCM07s(ff=M07sAcostL,pars_v=pars1.m07s.nh4,pars0=pars1.m07s.nh4,pars_fix=pars1.m07s.nh4,obs=obs_nh4no3[,1:5],niter1=niter_1,y0=y0.no3[1:10],times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
#---------------------   
	pars0.M14.nh4<-pars0.M14
	pars0.M14.nh4["C_NRec_2_NH4"] =  as.numeric(pars1.m07s.nh4["C_NRec_2_NH4"]);
	pars0.M14.nh4["C_NH4_2_NRec"]  =  as.numeric(pars1.m07s.nh4["C_NH4_2_NRec"]);
	pars0.M14.nh4["C_NH4_2_NLab"] =  as.numeric(pars1.m07s.nh4["C_NH4_2_NLab"]);
	pars0.M14.nh4["K_NLab_2_NH4"] =  as.numeric(pars1.m07s.nh4["K_NLab_2_NH4"]);
	pars0.M14.nh4["K_NH4ads_2_NH4"]  =  as.numeric(pars1.m07s.nh4["K_NH4ads_2_NH4"]);
	pars0.M14.nh4["K_NH4_2_NH4ads"] =  as.numeric(pars1.m07s.nh4["K_NH4_2_NH4ads"]);
	pars0.M14.nh4["C_NO3_2_NRec"] =  as.numeric(pars1.m07s.nh4["C_NO3_2_NRec"]);
	pars0.M14.nh4["K_NO3_2_NH4"] =  as.numeric(pars1.m07s.nh4["K_NO3_2_NH4"]);
	pars0.M14.nh4["C_NRec_2_NO3"] =  as.numeric(pars1.m07s.nh4["C_NRec_2_NO3"]);
	
	ii.fix<-which(names(pars0.M14.nh4) %in% names(pars1.m07s.nh4))
	ii.v<-which(!(names(pars0.M14.nh4) %in% names(pars1.m07s.nh4)))
	
	lower1=pars0.M14.nh4[ii.v]/lfactor;
	upper1=pars0.M14.nh4[ii.v]*ufactor;
	
	MCMC4_1.NH4<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.nh4[ii.v],pars0=pars0.M14.nh4,pars_fix=pars0.M14.nh4[ii.fix],obs=obs_nh4,niter1=niter_1,y0=y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars1.M14.nh4<-pars0.M14.nh4;
	pars1.M14.nh4[ii.v]<-MCMC4_1.NH4$bestpar;
	
	#
	pars0.M14.nh4no3<-pars0.M14
	pars0.M14.nh4no3["C_NRec_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NH4"]);
	pars0.M14.nh4no3["C_NH4_2_NRec"]  =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NRec"]);
	pars0.M14.nh4no3["C_NH4_2_NLab"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NH4_2_NLab"]);
	pars0.M14.nh4no3["K_NLab_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NLab_2_NH4"]);
	pars0.M14.nh4no3["K_NH4ads_2_NH4"]  =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4ads_2_NH4"]);
	pars0.M14.nh4no3["K_NH4_2_NH4ads"] =  as.numeric(MCMC01.NH4NO3$bestpar["K_NH4_2_NH4ads"]);
	pars0.M14.nh4no3["C_NO3_2_NRec"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NRec"]);
	pars0.M14.nh4no3["C_NO3_2_NH4"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NO3_2_NH4"]);
	pars0.M14.nh4no3["C_NRec_2_NO3"] =  as.numeric(MCMC01.NH4NO3$bestpar["C_NRec_2_NO3"]);
	
	ii.fix<-which(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar))
	ii.v<-which(!(names(pars0.M14.nh4no3) %in% names(MCMC01.NH4NO3$bestpar)))
	
	lower1=pars0.M14.nh4no3[ii.v]/lfactor;
	upper1=pars0.M14.nh4no3[ii.v]*ufactor;
	
	MCMC4_1.NH4NO3<-try(mdMCMC14(ff=M14costL,pars_v=pars0.M14.nh4no3[ii.v],pars0=pars0.M14.nh4no3,pars_fix=pars0.M14.nh4no3[ii.fix],obs=obs_nh4no3,niter1=niter_1,y0=y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
	
	pars1.M14.nh4no3<-pars0.M14.no3;
	pars1.M14.nh4no3[ii.v]<-MCMC4_1.NH4NO3$bestpar;
	
}
	
	#Parameters
	#pars1.M14.nh4no3
	#pars2.M14.nh4no3
	#pars2.M14.no3
	#pars1.M14.nh4
	#MCMC objects
	#MCMC4_1.NH4NO3
	#MCMC4_1.NH4
	#MCMC4_2.NH4NO3
	#MCMC4_1.NO3
	
	M14_1.nh4no3 <- (M14s(pars1.M14.nh4no3, y0.nh4no3,times))
	M14_2.nh4no3 <- M14s(pars2.M14.nh4no3,y0.nh4no3,times)
	M14_1.nh4<-M14s(pars1.M14.nh4,y0.nh4,times)
	M14_2.no3<-M14s(pars2.M14.no3,y0.no3,times)
		  
	cor_m14_1_nh4no3_obs<-cor(M14_1.nh4no3[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')],obs_nh4no3[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')])
	cor_m14_2_nh4no3_obs<-cor(M14_2.nh4no3[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')],obs_nh4no3[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')])
	cor_m14_1_nh4_obs<-cor(M14_1.nh4[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')],obs_nh4[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')])
	cor_m14_2_no3_obs<-cor(M14_2.no3[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')],obs_no3[,c('NH4','NH4_15','NO3','NO3_15','N2O','N2O_15')])
	  
	  print('ss.nh4 vs obs_nh4')
	  print(c('M14_1.nh4no3',format(diag(cor_m14_1_nh4no3_obs),digits=4)))
	  print(c('M14_2.nh4no3',format(diag(cor_m14_2_nh4no3_obs),digits=4)))
	  print(c('M14_1.nh4',format(diag(cor_m14_1_nh4_obs),digits=4)))
	  print(c('M14_2.no3',format(diag(cor_m14_2_no3_obs),digits=4)))
	  
	  #iq<-sample.int(n=nrow(MCMC$pars),size=50000)
	  var2sav<-c("y0.no3","y0.nh4","y0.nh4no3","obs_nh4","obs_no3","obs_nh4no3", "pars1.M14.nh4no3","pars2.M14.nh4no3","pars2.M14.no3","pars1.M14.nh4",
			  "MCMC4_1.NH4NO3","MCMC4_1.NH4","MCMC4_2.NH4NO3","MCMC4_1.NO3","MCMC0.NH4","MCMC0.NO3","M14_2.no3","M14_1.nh4","M14_2.nh4no3","M14_1.nh4no3",
			 "cor_m14_1_nh4no3_obs","cor_m14_2_nh4no3_obs","cor_m14_1_nh4_obs","cor_m14_2_no3_obs")
	  save(
	    list = var2sav[(var2sav %in% ls())],
	    file = paste(out_dir,"//", grpPal$ID[i],"-CS-",rngjk,'-',mod_type,'.RData', sep = '')
	  )
}	  
     dyn.unload(lib_nm1)
	 dyn.unload(lib_nm2)
	 
    sink()
     # c(lvlid[i], MCMC$bestpar)
    #-----------------------------------------------------------------------------
}

stopCluster(cl)

