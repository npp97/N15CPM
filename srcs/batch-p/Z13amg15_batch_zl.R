
cl<-makeCluster(ncl)
registerDoParallel(cl)

pars<-foreach(i=1:lgrpPal,.packages=c('FME','compiler','doBy'),.combine=rbind) %dopar% {
  rngjk=grpPal$rndseed[i]
  set.seed(rngjk)  
  dyn.load(lib_nm)
# for (i in 1:nlvlid){
    sink(paste(out_dir,'//rst_',grpPal$ID[i],"-CS-",rngjk,'-',mod_type,'.txt',sep=''))

    ii_nh4 <- which((obs_al$ID %in% grpPal$ID[i])&(obs_al$lable %in% Lb3[1]))
    jk.nh4 <- which((ini_val$ID %in% grpPal$ID[i])&(ini_val$lable %in% Lb3[1]))  
    obs_nh4 <- obs_al[ii_nh4, 3:ncol(obs_al)]
    nn.nh4<-get_n2oc15(time0=obs_nh4$time,n2o=obs_nh4$N2O,n2o_15=obs_nh4$N2O_15)
    obs_nh4[,'N2O']=nn.nh4[,'N2O'];obs_nh4[,'N2O_15']=nn.nh4[,'N2O_15'];

#    obs_nh4<-get_xx(obs_nh4);
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
               NO3sto = ini_val$iNO3sto[jk.nh4],
               NO3sto_15 = ini_val$iNO3sto_15[jk.nh4]*0.01*ini_val$iNO3sto[jk.nh4],
               N2O = 0,
               N2O_15 = 0
      );
    nh40.nh4<-coef(lm(NH4~time,data=obs_nh4[1:2,]))[1]
    no30.nh4<-coef(lm(NO3~time,data=obs_nh4[1:2,]))[1]
    nh4_15_0.nh4<-coef(lm(NH4_15~time,data=obs_nh4[1:2,]))[1]
    no3_15_0.nh4<-coef(lm(NO3_15~time,data=obs_nh4[1:2,]))[1]
    n2o0.nh4<-coef(lm(N2O~time,data=obs_nh4[1:2,]))[1]
    n2o_15_0.nh4<-coef(lm(N2O_15~time,data=obs_nh4[1:2,]))[1]

    #nh4_14_0.nh4<-coef(lm(NH4_14~time,data=obs_nh4[1:2,]))[1]
    #no3_14_0.nh4<-coef(lm(NO3_14~time,data=obs_nh4[1:2,]))[1]
    
    obs_nh40<-c(time=0,NH4=nh40.nh4,NO3=no30.nh4,NH4_15=nh4_15_0.nh4,NO3_15=no3_15_0.nh4,N2O=max(0,n2o0.nh4),N2O_15=max(0,n2o_15_0.nh4));#,NH4_14=nh4_14_0.nh4,NO3_14=no3_14_0.nh4);
    obs_nh4<-rbind(obs_nh40,obs_nh4);

    y0.nh4['NH4']<-obs_nh4[1,'NH4']
    y0.nh4['NO3']<-obs_nh4[1,'NO3']; 
    y0.nh4['NH4_15']<-obs_nh4[1,'NH4_15'];
    y0.nh4['NO3_15']<-obs_nh4[1,'NO3_15'];
 #   y0.nh4['N2O']<-obs_nh4[1,"N2O"];
 #   y0.nh4['N2O_15']<-obs_nh4[1,"N2O_15"];
    
    
    ii_no3 <- which((obs_al$ID %in% grpPal$ID[i])&(obs_al$lable %in% Lb3[2]))
    jk.no3 <- which((ini_val$ID %in% grpPal$ID[i])&(ini_val$lable %in% Lb3[2]))  
    obs_no3 <- obs_al[ii_no3, 3:ncol(obs_al)]
    nn.no3<-get_n2oc15(time0=obs_no3$time,n2o=obs_no3$N2O,n2o_15=obs_no3$N2O_15)
    obs_no3[,'N2O']=nn.no3[,'N2O'];obs_no3[,'N2O_15']=nn.no3[,'N2O_15'];
 #   obs_no3<-get_xx(obs_no3);

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
               NO3sto = ini_val$iNO3sto[jk.no3],
               NO3sto_15 = ini_val$iNO3sto_15[jk.no3]*0.01*ini_val$iNO3sto[jk.no3],
               N2O = 0,
               N2O_15 = 0
    );
    nh40.no3<-coef(lm(NH4~time,data=obs_no3[1:2,]))[1]
    no30.no3<-coef(lm(NO3~time,data=obs_no3[1:2,]))[1]
    nh4_15_0.no3<-coef(lm(NH4_15~time,data=obs_no3[1:2,]))[1]
    no3_15_0.no3<-coef(lm(NO3_15~time,data=obs_no3[1:2,]))[1]
    n2o0.no3<-coef(lm(N2O~time,data=obs_no3[1:2,]))[1]
    n2o_15_0.no3<-coef(lm(N2O_15~time,data=obs_no3[1:2,]))[1]

    #nh4_14_0.no3<-coef(lm(NH4_14~time,data=obs_no3[1:2,]))[1]
    #no3_14_0.no3<-coef(lm(NO3_14~time,data=obs_no3[1:2,]))[1]

    obs_no30<-c(time=0,NH4=nh40.no3,NO3=no30.no3,NH4_15=nh4_15_0.no3,NO3_15=no3_15_0.no3,N2O=max(0,n2o0.no3),N2O_15=max(0,n2o_15_0.no3));#,NH4_14=nh4_14_0.no3,NO3_14=no3_14_0.no3)
    obs_no3<-rbind(obs_no30,obs_no3);
    y0.no3['NH4']<-obs_no3[1,'NH4']
    y0.no3['NO3']<-obs_no3[1,'NO3']; 
    y0.no3['NH4_15']<-obs_no3[1,'NH4_15'];
    y0.no3['NO3_15']<-obs_no3[1,'NO3_15'];
#    y0.no3['N2O']<-obs_no3[1,"N2O"];
#    y0.no3['N2O_15']<-obs_no3[1,"N2O_15"];
       
#--------------Parameters
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
    # 									
    # 
    times <- obs_nh4[, 'time']
    
    ### MCMC simulation
    #----------------------------------------------------------------------------
    pp<-abs(pars0);
    upper<-lower<-pars0;
    lower1= pars0/lfactor;
    #lower[]<-0;
    upper1= pars0*ufactor;
    #lower1[c('alf_NH4_N2O','alf_NO3_N2O')] = 0;
    #upper1[c('alf_NH4_N2O','alf_NO3_N2O')] = 0.7;
    #upper=apply(data.frame(rep(1,length(pars0)),upper),1,min);
     #print("---------------NH4-N15")
    try({
    MCMC.NH4<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=obs_nh4,niter1=niter_1,niter2=niter_2,y0=y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower1,upper=upper1,burninlength=burnin.length));
    #print(MCMC.NH4$bestpar)
    if(class(MCMC.NH4) %in% "modMCMC"){
           
      pp["C_NRec_2_NH4"] = MCMC.NH4$bestpar["C_NRec_2_NH4"];
      pp["Vmax_NH4_2_NRec"]  = MCMC.NH4$bestpar["Vmax_NH4_2_NRec"];
      pp["Km_NH4_2_NRec"] = MCMC.NH4$bestpar["Km_NH4_2_NRec"];
      pp["Vmax_NH4_2_NLab"] = MCMC.NH4$bestpar["Vmax_NH4_2_NLab"];
      pp["Km_NH4_2_NLab"] = MCMC.NH4$bestpar["Km_NH4_2_NLab"];
      pp["Vmax_NLab_2_NH4"] = MCMC.NH4$bestpar["Vmax_NLab_2_NH4"];
      pp["Km_NLab_2_NH4"] = MCMC.NH4$bestpar["Km_NLab_2_NH4"];
      pp["K_NH4ads_2_NH4"]  = MCMC.NH4$bestpar["K_NH4ads_2_NH4"];
      pp["K_NH4_2_NH4ads"] = MCMC.NH4$bestpar["K_NH4_2_NH4ads"];
      pp["Vmax_NH4_2_NO3"]  = MCMC.NH4$bestpar["Vmax_NH4_2_NO3"];
      pp["Km_NH4_2_NO3"] = MCMC.NH4$bestpar["Km_NH4_2_NO3"];
      pp["K_NH4_2_N2O"] = MCMC.NH4$bestpar["K_NH4_2_N2O"];

   
     #de.nh4<-SCEoptim(Z13amg15costL,pp,obs=obs_nh4,y0=y0.nh4,times=times,lower=lower1,upper=upper1,control=list(maxit=maxit.sce1)) 
     #print("-----------SUMMARY of Fitted Parameters of NH4-N15------------------")
      #ig <- which(diff(MCMC.NH4$SS) <= 0)
      #print(summary(as.mcmc(MCMC.NH4$pars[ig, ])))
      print("-----------Best Fitted Parameters of NH4-N15----------------------")
        print(MCMC.NH4$bestpar);
      #print(de.nh4$par)
       #print(MCMC.NH4$bestfunp);
    }
    })
    
    #print("---------------NO3-N15")

    lower2=pp/lfactor;
    upper2=pp*ufactor;
    #lower2[c('alf_NH4_N2O','alf_NO3_N2O')] = 0;
    #upper2[c('alf_NH4_N2O','alf_NO3_N2O')] = 0.7;
   
    
    try({
    MCMC.NO3<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=obs_no3,niter1=niter_1,niter2=niter_2,y0=y0.no3,times,wvar0=0.1,updatecov=100,lower=lower2,upper=upper2,burninlength=burnin.length));
    #print(MCMC.NO3$bestpar)
    
    if(class(MCMC.NO3) %in% "modMCMC"){
      #print("-----------SUMMARY of Fitted Parameters of NO3-N15------------------")
      # ig <- which(diff(MCMC.NO3$SS) <= 0)
      # print(summary(as.mcmc(MCMC.NO3$pars[ig, ])))
      pp["Km_NO3_2_NRec"] = MCMC.NO3$bestpar["Km_NO3_2_NRec"];
      pp["Vmax_NO3_2_NRec"] = MCMC.NO3$bestpar["Vmax_NO3_2_NRec"];
      pp["C_NRec_2_NO3"] = MCMC.NO3$bestpar["C_NRec_2_NO3"];
      pp["Vmax_NO3_2_NH4"] = MCMC.NO3$bestpar["Vmax_NO3_2_NH4"];
      pp["Km_NO3_2_NH4"] = MCMC.NO3$bestpar["Km_NO3_2_NH4"];
      pp["C_NO3sto_2_NO3"] = MCMC.NO3$bestpar["C_NO3sto_2_NO3"];
      pp["K_NO3_2_NO3sto"] = MCMC.NO3$bestpar["K_NO3_2_NO3sto"];
      pp["K_NO3_2_N2O"] = MCMC.NO3$bestpar["K_NO3_2_N2O"];
      pp["K_NRec_2_N2O"] = MCMC.NO3$bestpar["K_NRec_2_N2O"];
#
      print("-----------Best Fitted Parameters  of NO3-N15----------------------")
      print(MCMC.NO3$bestpar);
      #print(MCMC.NO3$bestfunp);
      #de.no3<-DEoptim(fn=Z13amg15costL,lower=MCMC.NO3$bestpar/100,upper=MCMC.NO3$bestpar*100,obs=obs_no3,y0=y0.no3,times,DEoptim.control(NP=150,itermax=1000,CR=0.9))
      #de.no3<-SCEoptim(Z13amg15costL,pp,obs=obs_no3,y0=y0.no3,times=times,lower=lower2,upper=upper2,control=list(maxit=maxit.sce1)) 
      
      #print(de.no3$par)
    }
    })
    
      
    lower4=pp/lfactor;
    upper4=pp*ufactor;
    try({
    MCMC.NH4<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=obs_nh4,niter1=niter_1,niter2=niter_2,y0=y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower4,upper=upper4,burninlength=burnin.length));
    #print(MCMC.NH4$bestpar)
    if(class(MCMC.NH4) %in% "modMCMC"){
           
      pp["C_NRec_2_NH4"] = MCMC.NH4$bestpar["C_NRec_2_NH4"];
      pp["Vmax_NH4_2_NRec"]  = MCMC.NH4$bestpar["Vmax_NH4_2_NRec"];
      pp["Km_NH4_2_NRec"] = MCMC.NH4$bestpar["Km_NH4_2_NRec"];
      pp["Vmax_NH4_2_NLab"] = MCMC.NH4$bestpar["Vmax_NH4_2_NLab"];
      pp["Km_NH4_2_NLab"] = MCMC.NH4$bestpar["Km_NH4_2_NLab"];
      pp["Vmax_NLab_2_NH4"] = MCMC.NH4$bestpar["Vmax_NLab_2_NH4"];
      pp["Km_NLab_2_NH4"] = MCMC.NH4$bestpar["Km_NLab_2_NH4"];
      pp["K_NH4ads_2_NH4"]  = MCMC.NH4$bestpar["K_NH4ads_2_NH4"];
      pp["K_NH4_2_NH4ads"] = MCMC.NH4$bestpar["K_NH4_2_NH4ads"];
      pp["Vmax_NH4_2_NO3"]  = MCMC.NH4$bestpar["Vmax_NH4_2_NO3"];
      pp["Km_NH4_2_NO3"] = MCMC.NH4$bestpar["Km_NH4_2_NO3"];
      pp["K_NH4_2_N2O"] = MCMC.NH4$bestpar["K_NH4_2_N2O"];

   
     #de.nh4<-SCEoptim(Z13amg15costL,pp,obs=obs_nh4,y0=y0.nh4,times=times,lower=lower1,upper=upper1,control=list(maxit=maxit.sce1)) 
     #print("-----------SUMMARY of Fitted Parameters of NH4-N15------------------")
      #ig <- which(diff(MCMC.NH4$SS) <= 0)
      #print(summary(as.mcmc(MCMC.NH4$pars[ig, ])))
      print("-----------Best Fitted Parameters of NH4-N15----------------------")
      print(MCMC.NH4$bestpar);
      #print(de.nh4$par)
       #print(MCMC.NH4$bestfunp);
    }
    })  
      #print("-----------SUMMARY of Fitted Parameters of Total by SCE-UA------------------")
      
     # print(ss$par)
     # # try({
     #  MCMC<-mdMCMC2a(ff=M07costL,pp1=pp,obs=obs_nh4no3,niter1 = niter_1,niter2=niter_2,y0=y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower4,upper=upper4,burninlength=burnin.length);
     #    #})
     #  if(class(MCMC) %in% "modMCMC"){
     #    #print("-----------SUMMARY of Fitted Parameters of Total------------------")
     #    #ig <- which(diff(MCMC$SS) <= 0)
     #    #print(summary(as.mcmc(MCMC$pars[ig, ])))
     #    print("-----------Best Fitted Parameters of Total using NH4NO3----------------------")
     #    print(MCMC$bestpar);
     #    #print(MCMC$bestfunp);
     #  }
     #  })
     # 
        # #print("-----------SUMMARY of Fitted Parameters of Total by SCE-UA------------------")
      # lower4=apply(rbind(MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar,MCMC$bestpar),2,FUN=min)/lfactorf;
      # upper4=apply(rbind(MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar,MCMC$bestpar),2,FUN=max)*ufactorf;
    
    #PUT　ALL INTO ONE 
      try({
  #    lower4=apply(rbind(MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar,MCMC$bestpar,MCMC.TNN$bestpar),2,FUN=min)/lfactorf;
   #   upper4=apply(rbind(MCMC.NH4$bestpar,MCMC.NO3$bestpar,MCMC.NH4NO3$bestpar,MCMC$bestpar,MCMC.TNN$bestpar),2,FUN=median)*ufactorf;      
    
	  lower4=MCMC.NH4$bestpar/lfactorf;
	  upper4=MCMC.NH4$bestpar*ufactorf;
	  #lower4[c('alf_NH4_N2O','alf_NO3_N2O')] = 0;
    #upper4[c('alf_NH4_N2O','alf_NO3_N2O')] = 0.7;
   
	  #MCMC1<-mdMCMC_3ml(ff=Z13amcost_3l,pp1=pp,lower=lower4,upper=upper4,obs_nh4=obs_nh4,obs_no3=obs_no3,obs_n4n3=obs_nh4no3,y0_nh4=y0.nh4,y0_no3=y0.no3,y0_n4n3=y0.nh4no3,times,var0=0.5,wvar0=0.1,niter1=niter_1,niter2=niter_2,updatecov=100,burninlength=burnin.length)
	  MCMC<-modMCMC(f=Z13amg15cost_2l,p=pp,y0_nh4=y0.nh4,y0_no3=y0.no3,times=times,obs_nh4=obs_nh4,obs_no3=obs_no3,niter=niter_2,wvar0=0.1,updatecov=100,lower=lower4,upper=upper4,burninlength=burnin.length);

	  	 # 	  MCMC<-modMCMC(f=Z13amcost_3l,p=pp,y0_nh4=y0.nh4,y0_no3=y0.no3,y0_n4n3=y0.nh4no3,times=times,obs_nh4=obs_nh4,obs_no3=obs_no3,obs_n4n3=obs_nh4no3,niter=niter_2,var0=0.5,wvar0=0.1,updatecov=100,lower=lower4,upper=upper4,burninlength=burnin.length);
 
        #})
      if(class(MCMC) %in% "modMCMC"){
        print("-----------SUMMARY of Fitted Parameters of Total 6 curves------------------")
        #ig <- which(diff(MCMC$SS) <= 0)
        #print(summary(as.mcmc(MCMC$pars[ig, ])))
        print("-----------Best Fitted Parameters of Total----------------------")
        print(MCMC$bestpar);
        #de.6<-SCEoptim(M07cost_3ll,MCMC$bestpar,lower=lower4,upper=upper4,y0_nh4=y0.nh4,y0_no3=y0.no3,y0_n4n3=y0.nh4no3,times=times,
        #               obs_nh4=obs_nh4,obs_no3=obs_no3,obs_n4n3=obs_nh4no3,control=list(maxit=maxit.sce1))
      # de.6<-SCEoptim(Z13amg15cost_3ll,pp,lower=lower4,upper=upper4,obs_n4n3=obs_nh4no3,obs_nh4=obs_nh4,obs_no3=obs_no3,y0_n4n3=y0.nh4no3,y0_nh4=y0.nh4,y0_no3=y0.no3,times,control=list(maxit=maxit.sce1))
        
        #print(de.6$par)
        
         #print(MCMC$bestfunp);
      }
	  
	  ss.nh4<-Z13amg15(MCMC$bestpar,y0.nh4,times)
	  ss.no3<-Z13amg15(MCMC$bestpar,y0.no3,times)

#	  sss.nh4no3 <- (M07(y0.nh4no3,de.6$par, times))
#	  sss.nh4<-(M07(y0.nh4,de.6$par,times))
#	  sss.no3<-(Z13am(y0.no3,de.6$par,times))
	
#	  s3.nh4no3 <- (M07(y0.nh4no3,de.nh4no3$par, times))
#	  s3.nh4<-(M07(y0.nh4,de.nh4$par,times))
#	  s3.no3<-(Z13am(y0.no3,de.no3$par,times))
	  
#	  sss.nh4no3 <- M07(y0.nh4no3,MCMC.TNN$bestpar, times)
#	  sss.nh4<-M07(y0.nh4,MCMC.TNN$bestpar,times)
#	  sss.no3<-M07(y0.no3,MCMC.TNN$bestpar,times)
	  
	  s4.nh4<-Z13amg15(MCMC.NH4$bestpar,y0.nh4,times)
	  s4.no3<-Z13amg15(MCMC.NO3$bestpar,y0.no3,times)
	  
	  cor_ss_nh4_obs<-cor_z13amg15(obs_nh4,MCMC$bestpar,y0.nh4,times)
	  cor_ss_no3_obs<-cor_z13amg15(obs_no3,MCMC$bestpar,y0.no3,times)

	  cor_s4_nh4_obs<-cor_z13amg15(obs_nh4,MCMC.NH4$bestpar,y0.nh4,times)
	  cor_s4_no3_obs<-cor_z13amg15(obs_no3,MCMC.NO3$bestpar,y0.no3,times)
	  
    # cor_sss_nh4_obs<-cor_z13amg15(obs_nh4,de.6$par,y0.nh4,times)
	  # cor_sss_no3_obs<-cor_z13amg15(obs_no3,de.6$par,y0.no3,times)
	  # cor_sss_nh4no3_obs<-cor_z13amg15(obs_nh4no3,de.6$par,y0.nh4no3,times)

	  # cor_s3_nh4_obs<-cor_z13amg15(obs_nh4,de.nh4$par,y0.nh4,times)
	  # cor_s3_no3_obs<-cor_z13amg15(obs_no3,de.no3$par,y0.no3,times)
	  # cor_s3_nh4no3_obs<-cor_z13amg15(obs_nh4no3,de.nh4no3$par,y0.nh4no3,times)

#	  cor_sss_nh4_obs<-cor(sss.nh4[,c('NH4','NH4_15','NO3','NO3_15')],obs_nh4[,c('NH4','NH4_15','NO3','NO3_15')])
#	  cor_sss_no3_obs<-cor(sss.no3[,c('NH4','NH4_15','NO3','NO3_15')],obs_no3[,c('NH4','NH4_15','NO3','NO3_15')])
#	  cor_sss_nh4no3_obs<-cor(sss.nh4no3[,c('NH4','NH4_15','NO3','NO3_15')],obs_nh4no3[,c('NH4','NH4_15','NO3','NO3_15')])
	  
#	  cor_s3_nh4_obs<-cor(s3.nh4[,c('NH4','NH4_15','NO3','NO3_15')],obs_nh4[,c('NH4','NH4_15','NO3','NO3_15')])
#	  cor_s3_no3_obs<-cor(s3.no3[,c('NH4','NH4_15','NO3','NO3_15')],obs_no3[,c('NH4','NH4_15','NO3','NO3_15')])
#	  cor_s3_nh4no3_obs<-cor(s3.nh4no3[,c('NH4','NH4_15','NO3','NO3_15')],obs_nh4no3[,c('NH4','NH4_15','NO3','NO3_15')])
	  
	  
	  
	  print("--ss.nh4-----")
	  print(ss.nh4)
	  print(obs_nh4)
#	  print(sss.nh4)
	  print(s4.nh4)
	  print('ss.nh4 vs obs_nh4')
	  print(c('mcmc.6cv',format(diag(cor_ss_nh4_obs),digits=4)))
	  print(c('mcmc.ibm',format(diag(cor_s4_nh4_obs),digits=4)))
	  #print(c('sce.6cv',format(diag(cor_sss_nh4_obs),digits=4)))
	  #print(c('sce.ibm',format(diag(cor_s3_nh4_obs),digits=4)))
	  
	  
	  print("ss.no3------------")
	  print(ss.no3)
	  print(obs_no3)
#	  print(sss.no3)
	  print(s4.no3)
	  print('ss.no3 vs obs_no3')
	  print(c('mcmc.6cv',format(diag(cor_ss_no3_obs),digits=4)))
	  print(c('mcmc.ibm',format(diag(cor_s4_no3_obs),digits = 4)))
	  #print(c('sce.6cv',format(diag(cor_sss_no3_obs),digits=4)))
	  #print(c('sce.ibm',format(diag(cor_s3_no3_obs),digits=4)))
	  

	  #iq<-sample.int(n=nrow(MCMC$pars),size=50000)
	  var2sav<-c("ss.nh4","ss.no3","s4.nh4","s4.no3","y0.no3","y0.nh4","MCMC","cor_ss_nh4_obs",
	             #"sss.nh4","sss.no3","sss.nh4no3","s3.nh4","s3.no3","s3.nh4no3","cor_sss_no3_obs","cor_sss_nh4no3_obs","cor_s3_nh4_obs",
	             #"cor_s3_no3_obs","cor_s3_nh4no3_obs","cor_sss_nh4_obs","cor_s3_nh4_obs",
	             "cor_ss_no3_obs","cor_s4_nh4_obs","cor_s4_no3_obs",
	             "obs_nh4","obs_no3","ss.nh4","ss.no3","obs_al", "ini_val","MCMC.NH4","MCMC.NO3")
	  save(
	    list = var2sav[(var2sav %in% ls())],
	    file = paste(out_dir,"//", grpPal$ID[i],"-CS-",rngjk,'-',mod_type,'.RData', sep = '')
	  )
	  
      })
     dyn.unload(lib_nm)
    sink()
     # c(lvlid[i], MCMC$bestpar)
    #-----------------------------------------------------------------------------
}

stopCluster(cl)

