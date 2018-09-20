
  fn_n_rseed<-function(fn_str,fntrt_p=3){
      require(stringr);
      fn_str<-basename(fn_str);
      aposi<-as.data.frame(str_locate_all(fn_str,'-'));
      fn_trt<-substr(fn_str,1,aposi[fntrt_p,1]-1);
      fn_rseed<-as.numeric(substr(fn_str,aposi[nrow(aposi)-1,1]+1,aposi[nrow(aposi),1]-1));
      fn_n_rseed<-list(fn_trt=fn_trt,fn_rseed=fn_rseed);
      return(fn_n_rseed)
  }
  
  getpar_fnal_z13amg<-function(arst,fntrt_p=3,lib_nm,lfactor=10,ufactor=10,niter_1=50000,niter_2=50000,burninlength=10000){
      dyn.load(lib_nm);
     times=arst$obs_nh4[,'time'];
      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NH4.bstr$bestpar
      upper=pp*ufactor;
      lower=pp/lfactor;
      MCMC.NH4<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=arst$obs_nh4,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));


      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NO3.bstr$bestpar
      upper=pp*ufactor;
      lower=pp/lfactor;
      MCMC.NO3<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=arst$obs_no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NH4NO3.bstr$bestpar
      upper=pp*ufactor;
      lower=pp/lfactor;
      MCMC.NH4NO3<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=arst$obs_nh4no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

	  s4.nh4no3 <- Z13amg15(MCMC.NH4NO3$bestpar,y0.nh4no3, times)
	  s4.nh4<-Z13amg15(MCMC.NH4$bestpar,y0.nh4,times)
	  s4.no3<-Z13amg15(MCMC.NO3$bestpar,y0.no3,times)
	  
	  cor_s4_nh4_obs<-cor_z13amg15(obs_nh4,MCMC.NH4$bestpar,y0.nh4,times)
	  cor_s4_no3_obs<-cor_z13amg15(obs_no3,MCMC.NO3$bestpar,y0.no3,times)
	  cor_s4_nh4no3_obs<-cor_z13amg15(obs_nh4no3,MCMC.NH4NO3$bestpar,y0.nh4no3,times)
	  dyn.unload(lib_nm); 
 	  
    var2sav<-list(arst,
                MCMC.NH4.fin=MCMC.NH4,MCMC.NO3.fin=MCMC.NO3,MCMC.NH4NO3.fin=MCMC.NH4NO3,
                s4.nh4=s4.nh4,s4.no3=s4.no3,s4.nh4no3=s4.nh4no3,
	            cor_s4_nh4_obs=cor_s4_nh4_obs,cor_s4_no3_obs=cor_s4_no3_obs,cor_s4_nh4no3_obs=cor_s4_nh4no3_obs)

     return(var2sav)
  }
 
 
  getpar_fnal_z13am<-function(arst,fntrt_p=3,lib_nm,lfactor=10,ufactor=10,niter_1=50000,niter_2=50000,burninlength=10000){
      dyn.load(lib_nm);
      times=arst$obs_nh4[,'time'];
      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NH4.bstr$bestpar
      upper=pp*ufactor;
      lower=pp/lfactor;
      MCMC.NH4<-try(mdMCMC2a(ff=Z13amcost,pp1=pp,obs=arst$obs_nh4,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));


      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NO3.bstr$bestpar
      upper=pp*ufactor;
      lower=pp/lfactor;
      MCMC.NO3<-try(mdMCMC2a(ff=Z13amcost,pp1=pp,obs=arst$obs_no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NH4NO3.bstr$bestpar
      upper=pp*ufactor;
      lower=pp/lfactor;
      MCMC.NH4NO3<-try(mdMCMC2a(ff=Z13amcost,pp1=pp,obs=arst$obs_nh4no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

	  s4.nh4no3 <- Z13amg15(MCMC.NH4NO3$bestpar,y0.nh4no3, times)
	  s4.nh4<-Z13amg15(MCMC.NH4$bestpar,y0.nh4,times)
	  s4.no3<-Z13amg15(MCMC.NO3$bestpar,y0.no3,times)
	  
	  cor_s4_nh4_obs<-cor_z13am(obs_nh4,MCMC.NH4$bestpar,y0.nh4,times)
	  cor_s4_no3_obs<-cor_z13am(obs_no3,MCMC.NO3$bestpar,y0.no3,times)
	  cor_s4_nh4no3_obs<-cor_z13am(obs_nh4no3,MCMC.NH4NO3$bestpar,y0.nh4no3,times)
	 
 	  var2sav<-list(arst,
                MCMC.NH4.fin=MCMC.NH4,MCMC.NO3.fin=MCMC.NO3,MCMC.NH4NO3.fin=MCMC.NH4NO3,
                s4.nh4=s4.nh4,s4.no3=s4.no3,s4.nh4no3=s4.nh4no3,
	            cor_s4_nh4_obs=cor_s4_nh4_obs,cor_s4_no3_obs=cor_s4_no3_obs,cor_s4_nh4no3_obs=cor_s4_nh4no3_obs)
     dyn.unload(lib_nm);
     return(var2sav)
  }
 