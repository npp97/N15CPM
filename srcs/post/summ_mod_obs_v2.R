library(coda);

findBST<-function(flst){
    load(flst[1])

    ibest.nh4no3<-ibest.no3<-ibest.nh4<-1;
    ibestr.nh4no3<-ibestr.no3<-ibestr.nh4<-1;

    if (length(which(is.na(s4.nh4[,'NH4'])))>1){
        load(flst[2]);
        ibest.nh4no3<-ibest.no3<-ibest.nh4<-2;
        ibestr.nh4no3<-ibestr.no3<-ibestr.nh4<-2;
    }
    s4.nh4.fb<-rbind(s4.nh4)#,s4.nh4)
    s4.no3.fb<-rbind(s4.no3)#,s4.no3)
    s4.nh4no3.fb<-rbind(s4.nh4no3)#,s4.nh4no3)

    best.nh4<-s4.nh4
    bestf.nh4<-MCMC.NH4$bestfunp
    bestr.nh4<-s4.nh4
    bestfr.nh4<-sum(diag(cor_s4_nh4_obs[,-1])^2)

    best.no3<-s4.no3
    bestf.no3<-MCMC.NO3$bestfunp
    bestr.no3<-s4.no3
    bestfr.no3<-sum(diag(cor_s4_no3_obs[,-1])^2)

    best.nh4no3<-ss.nh4no3
    bestf.nh4no3<-MCMC.NH4NO3$bestfunp
    bestr.nh4no3<-s4.nh4no3
    bestfr.nh4no3<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)

    #obs_nh4.fb<-obs_nh4
    #obs_no3.fb<-obs_no3
    #obs_nh4no3.fb <- obs_nh4no3

for (i in 2:length(flst)){
    load(flst[i])
    ibest.nh4no3<-ibest.no3<-ibest.nh4<-i;
    ibestr.nh4no3<-ibestr.no3<-ibestr.nh4<-i;

    s4.nh4.fb<-rbind(s4.nh4.fb,s4.nh4)#,s4.nh4)
    s4.no3.fb<-rbind(s4.no3.fb,s4.no3)#,s4.no3)
    s4.nh4no3.fb<-rbind(s4.nh4no3.fb,s4.nh4no3)#,s4.nh4no3)
    if (length(which(is.na(s4.nh4[,'NH4'])))>1) next;
    if (length((is.na(c(diag(cor_s4_nh4_obs),diag(cor_s4_no3_obs),diag(cor_s4_nh4no3_obs)))))>1) next;

    if(bestf.nh4>MCMC.NH4$bestfunp){
        best.nh4<-s4.nh4;
        bestf.nh4<-MCMC.NH4$bestfunp;
        ibest.nh4<-i;
    }
    if(bestf.no3>MCMC.NO3$bestfunp){
        best.no3<-s4.no3;
        bestf.no3<-MCMC.NO3$bestfunp;
        ibest.no3<-i;
    }
    if(bestf.nh4no3>MCMC.NH4NO3$bestfunp){
        best.nh4no3<-s4.nh4no3;
        bestf.nh4no3<-MCMC.NH4NO3$bestfunp;
        ibest.nh4no3<-i;
    }
##
    tmp<-sum(diag(cor_s4_nh4_obs[,-1])^2)
    if(bestfr.nh4<tmp){
        bestr.nh4<-s4.nh4;
        bestfr.nh4<tmp;
        ibestr.nh4<-i;
    }
    tmp<-sum(diag(cor_s4_no3_obs[,-1])^2)
    if(bestfr.no3<tmp){
        bestr.no3<-s4.no3;
        bestf.no3<-tmp;
        ibestr.no3<-i;
    }
    tmp<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)
    if(bestfr.nh4no3<tmp){
        bestr.nh4no3<-s4.nh4no3;
        bestfr.nh4no3<-tmp;
        ibestr.nh4no3<-i;
    }
}
    load(flst[ibest.nh4]);
        MCMC.NH4.bst<-MCMC.NH4;
    load(flst[ibest.no3]);
        MCMC.NO3.bst<-MCMC.NO3;
    load(flst[ibest.nh4no3])
        MCMC.NH4NO3.bst<-MCMC.NH4NO3;

    load(flst[ibestr.nh4]);
        MCMC.NH4.bstr<-MCMC.NH4;
    load(flst[ibestr.no3]);
        MCMC.NO3.bstr<-MCMC.NO3;
    load(flst[ibestr.nh4no3])
        MCMC.NH4NO3.bstr<-MCMC.NH4NO3;
    #group variables to return
    arst<-list(MCMC.NH4.bst=MCMC.NH4.bst,MCMC.NO3.bst=MCMC.NO3.bst,MCMC.NH4NO3.bst=MCMC.NH4NO3.bst,
    MCMC.NH4.bstr=MCMC.NH4.bstr,MCMC.NO3.bstr=MCMC.NO3.bstr,MCMC.NH4NO3.bstr=MCMC.NH4NO3.bstr,
    obs_nh4=obs_nh4,obs_no3=obs_no3,obs_nh4no3=obs_nh4no3,
    y0.nh4=y0.nh4,y0_no3=y0.no3,y0.nh4no3=y0.nh4no3,
    fn.best.nh4=flst[ibest.nh4],fn.best.no3=flst[ibest.no3],fn.best.nh4no3=flst[ibest.nh4no3],
    fn.bestr.nh4=flst[ibestr.nh4],fn.bestr.no3=flst[ibestr.no3],fn.bestr.nh4no3=flst[ibestr.nh4no3], 
    best.nh4=as.data.frame(best.nh4),best.no3=as.data.frame(best.no3),best.nh4no3=as.data.frame(best.nh4no3),
    bestr.nh4=as.data.frame(bestr.nh4),bestr.no3=as.data.frame(bestr.no3),bestr.nh4no3=as.data.frame(bestr.nh4no3),
    s4.nh4=as.data.frame(s4.nh4.fb),s4.no3=as.data.frame(s4.no3.fb),s4.nh4no3=as.data.frame(s4.nh4no3.fb));

return(arst)    
}


valid_plotMCMC<-function(obs,mod_a){
mod_a<-mod_a[,names(obs)];
  tn<-names(mod_a);
  nob<-ncol(obs);
  yl<-apply(rbind(mod_a,obs),FUN=range,2,na.rm=T)
  nm<-ncol(mod_a);
  n<-floor(sqrt(nm-1));
  m<-floor((nm-1)/(n-1));
  times<-obs[,'time'];
  par0<-par();
#  win.graph();
  par(mfcol=c(3,2), mar=c(3,3,2,0.5))
  if(nrow(mod_a)>nrow(obs)){
         for (i in 2:(length(tn))){
            plot(x=mod_a[,'time'],y=mod_a[,i],type='p',pch='+',cex=2,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
            if(i<=(nob)) points(obs[,'time'],obs[,i],type='b',cex=2,pch=19,col=2);
         };
  } else{
        for (i in 2:(length(tn))){
            plot(x=mod_a[,'time'],y=mod_a[,i],type='l',lwd=1.5,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
            title(xlab='Time(Hour)',ylab=paste(tn[i],'(mg/kg)',sep=''),line=2)
        if(i<=(nob)) points(obs[,'time'],obs[,i],type='p',cex=2,pch=19,col=2);
  };

 }	
  return(1);
}

gen_fig<-function(arst,char_trt){

    pdf(paste('valid_best_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,best.nh4);
        valid_plotMCMC(obs_no3,best.no3);
        valid_plotMCMC(obs_nh4no3,best.nh4no3);
    })
    dev.off()

    pdf(paste('valid_bestr_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,bestr.nh4);
        valid_plotMCMC(obs_no3,bestr.no3);
        valid_plotMCMC(obs_nh4no3,bestr.nh4no3);
    })
    dev.off()

    pdf(paste('valid_s4_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,s4.nh4);
        valid_plotMCMC(obs_no3,s4.no3);
        valid_plotMCMC(obs_nh4no3,s4.nh4no3);
    })
    dev.off()

    pdf(paste('denMCMC_pars_bestr_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        densplot(as.mcmc(MCMC.NH4.bstr$pars[diff(MCMC.NH4.bstr$SS)<0,]))
        traceplot(as.mcmc(MCMC.NH4.bstr$pars[diff(MCMC.NH4.bstr$SS)<0,]))
        plot(as.mcmc(MCMC.NH4.bstr$pars[diff(MCMC.NH4.bstr$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NO3.bstr$pars[diff(MCMC.NO3.bstr$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NH4NO3.bstr$pars[diff(MCMC.NH4NO3.bstr$SS)<0,]),trace=T,smooth=T)  
    })
    dev.off()

    pdf(paste('denMCMC_pars_best_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        densplot(as.mcmc(MCMC.NH4.bst$pars[diff(MCMC.NH4.bst$SS)<0,]))
        traceplot(as.mcmc(MCMC.NH4.bst$pars[diff(MCMC.NH4.bst$SS)<0,]))
        plot(as.mcmc(MCMC.NH4.bst$pars[diff(MCMC.NH4.bst$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NO3.bst$pars[diff(MCMC.NO3.bst$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NH4NO3.bst$pars[diff(MCMC.NH4NO3.bst$SS)<0,]),trace=T,smooth=T)  
    })
    dev.off()
}


path0<-'D://N15CSg//run//out_NBFB-CBS-ZHANG2013amg_2018-02-19'
#path0<-'D://N15CSg//run//out_NBFBF-NE-60-ZHANG2013am'#
#path0<-'D://N15CSg//run//out_NBFBF-NE-60-ZHANG2013am'
#path0<-'D://N15CS1//run//out_SUSF-6090AVG-ZHANG2013a'
setwd(path0)

flst<-dir(pattern="*.RData")

ii_FB<-grep('FB',flst)
ii_NB<-grep('NB',flst)

rst_FB<-findBST(flst[ii_FB])
rst_NB<-findBST(flst[ii_NB])

gen_fig(rst_FB,'FB-CBM-100')
gen_fig(rst_NB,'NB-CBM-100')

########################################
path0<-'D://N15CSg//run//out_NBFB-CBS-ZHANG2013amg'

path0<-'D://N15CSg//run//out_NBFB-CBS-ZHANG2013am'
setwd(path0)

flst<-dir(pattern="*.RData")

ii_FB_CBM_60<-grep('FB-CBM-60',flst)
ii_NB_CBM_60<-grep('NB-CBM-60',flst)
ii_FB_CBM_100<-grep('FB-CBM-100',flst)
ii_NB_CBM_100<-grep('NB-CBM-100',flst)
ii_FF_CBM_60<-grep('FF-CBM-60',flst)
ii_NF_CBM_60<-grep('NF-CBM-60',flst)


rst_FB_CBM_60<-findBST(flst[ii_FB_CBM_60])
rst_NB_CBM_60<-findBST(flst[ii_NB_CBM_60])
rst_FB_CBM_100<-findBST(flst[ii_FB_CBM_100])
rst_NB_CBM_100<-findBST(flst[ii_NB_CBM_100])

rst_FF_CBM_60<-findBST(flst[ii_FF_CBM_60])
rst_NF_CBM_60<-findBST(flst[ii_NF_CBM_60])

gen_fig(rst_FB_CBM_60,'FB-CBM-60')
gen_fig(rst_NB_CBM_60,'NB-CBM-60')
gen_fig(rst_FB_CBM_100,'FB-CBM-100')
gen_fig(rst_NB_CBM_100,'NB-CBM-100')
gen_fig(rst_FF_CBM_60,'FF-CBM-60')
gen_fig(rst_NF_CBM_60,'NF-CBM-60')

#############################
path0<-'D://N15CSg//run//out_NBFB-CBS-ZHANG2013amg'
path0<-'D://N15CSg//run//out_NBFBF-NE-60-ZHANG2013am'
path0<-'D://N15CSg//run//out_pre_NBFB-NE-100-ZHANG2013amg';#FB-CBM,FF-XXM,NB-CBM,NF-XXM
#path0<-'D://N15CS//run//Naddition_ne__2018-01-09-ZHANG2013AM'

setwd(path0)

flst<-dir(pattern="*.RData")

ii_FB_CBM<-grep('FB-CBM',flst)
ii_NB_CBM<-grep('NB-CBM',flst)
ii_FF_XXM<-grep('FF-XXM',flst)
ii_NF_XXM<-grep('NF-XXM',flst)


rst_FB_CBM<-findBST(flst[ii_FB_CBM])
rst_NB_CBM<-findBST(flst[ii_NB_CBM])
rst_FF_XXM<-findBST(flst[ii_FF_XXM])
rst_NF_XXM<-findBST(flst[ii_NF_XXM])

gen_fig(rst_FB_CBM,'FB-CBM')
gen_fig(rst_NB_CBM,'NB-CBM')
gen_fig(rst_FF_XXM,'FF-XXM')
gen_fig(rst_NF_XXM,'NF-XXM')


###########################
path0<-'D://N15CSg//run//out_NBFB-CBS-ZHANG2013amg'
path0<-'D://N15CSg//run//out_pre_SUSF-CBS-60-ZHANG2013am'
#path0<-'D://N15CSg//run//out_pre_NBFB-NE-100-ZHANG2013amg';#FB-CBM,FF-XXM,NB-CBM,NF-XXM
#path0<-'D://N15CS//run//Naddition_ne__2018-01-09-ZHANG2013AM'

setwd(path0)

flst<-dir(pattern="*.RData")

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

gen_fig(rst_11,'SUSF-1-1')
gen_fig(rst_12,'SUSF-1-2')
gen_fig(rst_13,'SUSF-1-3')

gen_fig(rst_31,'SUSF-3-1')
gen_fig(rst_32,'SUSF-3-2')
gen_fig(rst_33,'SUSF-3-3')

gen_fig(rst_41,'SUSF-4-1')
gen_fig(rst_42,'SUSF-4-2')
gen_fig(rst_43,'SUSF-4-3')

gen_fig(rst_51,'SUSF-5-1')
gen_fig(rst_52,'SUSF-5-2')
gen_fig(rst_53,'SUSF-5-3')

gen_fig(rst_61,'SUSF-6-1')
gen_fig(rst_62,'SUSF-6-2')
gen_fig(rst_63,'SUSF-6-3')
