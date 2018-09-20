
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
  win.graph();
  par(mfcol=c(3,2), mar=c(2,3,2,0.5))
  if(nrow(mod_a)>nrow(obs)){
         for (i in 2:(length(tn))){
            plot(x=mod_a[,'time'],y=mod_a[,i],type='p',pch='+',cex=2,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
            if(i<=(nob)) points(obs[,'time'],obs[,i],type='b',cex=2,pch=19,col=2);
         };
  } else{
        for (i in 2:(length(tn))){
            plot(x=mod_a[,'time'],y=mod_a[,i],type='l',lwd=1.5,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
        if(i<=(nob)) points(obs[,'time'],obs[,i],type='p',cex=2,pch=19,col=2);
  };

 }	
  return(1);
}



path0<-'D://N15CSg//run//out_NBFB-CBS-ZHANG2013amg'
#path0<-'D://N15CSg//run//out_NBFBF-NE-60-ZHANG2013am'#
#path0<-'D://N15CSg//run//out_NBFBF-NE-60-ZHANG2013am'
#path0<-'D://N15CS1//run//out_SUSF-6090AVG-ZHANG2013a'
setwd(path0)

flst<-dir(pattern="*.RData")

ii_FB<-grep('FB',flst)
ii_NB<-grep('NB',flst)


load(flst[ii_FB[1]])

if (length(which(is.na(s4.nh4[,'NH4'])))>1) load(flst[ii_FB[2]]);
s4.nh4.fb<-rbind(s4.nh4)#,s4.nh4)
s4.no3.fb<-rbind(s4.no3)#,s4.no3)
s4.nh4no3.fb<-rbind(s4.nh4no3)#,s4.nh4no3)

best.nh4.fb<-s4.nh4
bestf.nh4.fb<-MCMC.NH4$bestfunp
bestr.nh4.fb<-s4.nh4
bestfr.nh4.fb<-sum(diag(cor_s4_nh4_obs[,-1])^2)

best.no3.fb<-s4.no3
bestf.no3.fb<-MCMC.NO3$bestfunp
bestr.no3.fb<-s4.no3
bestfr.no3.fb<-sum(diag(cor_s4_no3_obs[,-1])^2)

best.nh4no3.fb<-ss.nh4no3
bestf.nh4no3.fb<-MCMC.NH4NO3$bestfunp
bestr.nh4no3.fb<-s4.nh4no3
bestfr.nh4no3.fb<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)

obs_nh4.fb<-obs_nh4
obs_no3.fb<-obs_no3
obs_nh4no3.fb <- obs_nh4no3

for (i in 2:length(ii_FB)){
    load(flst[ii_FB[i]])

    s4.nh4.fb<-rbind(s4.nh4.fb,s4.nh4)#,s4.nh4)
    s4.no3.fb<-rbind(s4.no3.fb,s4.no3)#,s4.no3)
    s4.nh4no3.fb<-rbind(s4.nh4no3.fb,s4.nh4no3)#,s4.nh4no3)
    if (length(which(is.na(s4.nh4[,'NH4'])))>1) next;
    if (length((is.na(c(diag(cor_s4_nh4_obs),diag(cor_s4_no3_obs),diag(cor_s4_nh4no3_obs)))))>1) next;

    if(bestf.nh4.fb>MCMC.NH4$bestfunp){
        best.nh4.fb<-s4.nh4;
        bestf.nh4.fb<-MCMC.NH4$bestfunp;
    }
    if(bestf.no3.fb>MCMC.NO3$bestfunp){
        best.no3.fb<-s4.no3;
        bestf.no3.fb<-MCMC.NO3$bestfunp;
    }
    if(bestf.nh4no3.fb>MCMC.NH4NO3$bestfunp){
        best.nh4no3.fb<-s4.nh4no3;
        bestf.nh4no3.fb<-MCMC.NH4NO3$bestfunp;
    }
##
    tmp<-sum(diag(cor_s4_nh4_obs[,-1])^2)
    if(bestfr.nh4.fb<tmp){
        bestr.nh4.fb<-s4.nh4;
        bestfr.nh4.fb<tmp;
    }
    tmp<-sum(diag(cor_s4_no3_obs[,-1])^2)
    if(bestfr.no3.fb<tmp){
        bestr.no3.fb<-s4.no3;
        bestf.no3.fb<-tmp;
    }
    tmp<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)
    if(bestfr.nh4no3.fb<tmp){
        bestr.nh4no3.fb<-s4.nh4no3;
        bestfr.nh4no3.fb<-tmp;
    }

}


load(flst[ii_NB[1]])
if (length(which(is.na(s4.nh4[,'NH4'])))>1) load(flst[ii_FB[2]]);
s4.nh4.nb<-rbind(s4.nh4)#,s4.nh4)
s4.no3.nb<-rbind(s4.no3)#,s4.no3)
s4.nh4no3.nb<-rbind(s4.nh4no3)#,s4.nh4no3)

best.nh4.nb<-s4.nh4
bestf.nh4.nb<-MCMC.NH4$bestfunp
bestr.nh4.nb<-s4.nh4
bestfr.nh4.nb<-sum(diag(cor_s4_nh4_obs[,-1])^2)

best.no3.nb<-s4.no3
bestf.no3.nb<-MCMC.NO3$bestfunp
bestr.no3.nb<-s4.no3
bestfr.no3.nb<-sum(diag(cor_s4_no3_obs[,-1])^2)

best.nh4no3.nb<-s4.nh4no3
bestf.nh4no3.nb<-MCMC.NH4NO3$bestfunp
bestr.nh4no3.nb<-s4.nh4no3
bestfr.nh4no3.nb<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)

obs_nh4.nb<-obs_nh4
obs_no3.nb<-obs_no3
obs_nh4no3.nb <- obs_nh4no3

for (i in 2:length(ii_NB)){
    load(flst[ii_NB[i]])
    s4.nh4.nb<-rbind(s4.nh4.nb,s4.nh4)#,s4.nh4)
    s4.no3.nb<-rbind(s4.no3.nb,s4.no3)#,s4.no3)
    s4.nh4no3.nb<-rbind(s4.nh4no3.nb,s4.nh4no3)#,s4.nh4no3)
    if (length(which(is.na(s4.nh4[,'NH4'])))>1) next;
    if (length((is.na(c(diag(cor_s4_nh4_obs),diag(cor_s4_no3_obs),diag(cor_s4_nh4no3_obs)))))>1) next;

	
    if(bestf.nh4.nb>MCMC.NH4$bestfunp){
        best.nh4.nb<-s4.nh4;
        bestf.nh4.nb<-MCMC.NH4$bestfunp;
    }
    if(bestf.no3.nb>MCMC.NO3$bestfunp){
        best.no3.nb<-s4.no3;
        bestf.no3.nb<-MCMC.NO3$bestfunp;
    }
    if(bestf.nh4no3.nb>MCMC.NH4NO3$bestfunp){
        best.nh4no3.nb<-s4.nh4no3;
        bestf.nh4no3.nb<-MCMC.NH4NO3$bestfunp;
    }

    
    tmp<-sum(diag(cor_s4_nh4_obs[,-1])^2)
    if(bestfr.nh4.nb<tmp){
        bestr.nh4.nb<-s4.nh4;
        bestfr.nh4.nb<tmp;
    }
    tmp<-sum(diag(cor_s4_no3_obs[,-1])^2)
    if(bestfr.no3.nb<tmp){
        bestr.no3.nb<-s4.no3;
        bestf.no3.nb<-tmp;
    }
    tmp<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)
    if(bestfr.nh4no3.nb<tmp){
        bestr.nh4no3.nb<-s4.nh4no3;
        bestfr.nh4no3.nb<-tmp;
    }
}


valid_plotMCMC(obs_nh4.fb,as.data.frame(best.nh4.fb))
valid_plotMCMC(obs_nh4.nb,as.data.frame(best.nh4.nb))
valid_plotMCMC(obs_no3.fb,as.data.frame(best.no3.fb))
valid_plotMCMC(obs_no3.nb,as.data.frame(best.no3.nb))
valid_plotMCMC(obs_nh4no3.fb,as.data.frame(best.nh4no3.fb))
valid_plotMCMC(obs_nh4no3.nb,as.data.frame(best.nh4no3.nb))


valid_plotMCMC(obs_nh4.fb,as.data.frame(bestr.nh4.fb))
valid_plotMCMC(obs_nh4.nb,as.data.frame(bestr.nh4.nb))
valid_plotMCMC(obs_no3.fb,as.data.frame(bestr.no3.fb))
valid_plotMCMC(obs_no3.nb,as.data.frame(bestr.no3.nb))
valid_plotMCMC(obs_nh4no3.fb,as.data.frame(bestr.nh4no3.fb))
valid_plotMCMC(obs_nh4no3.nb,as.data.frame(bestr.nh4no3.nb))


valid_plotMCMC(obs_nh4.fb,as.data.frame(s4.nh4.fb))
valid_plotMCMC(obs_nh4.nb,as.data.frame(s4.nh4.nb))
valid_plotMCMC(obs_no3.fb,as.data.frame(s4.no3.fb))
valid_plotMCMC(obs_no3.nb,as.data.frame(s4.no3.nb))
valid_plotMCMC(obs_nh4no3.fb,as.data.frame(s4.nh4no3.fb))
valid_plotMCMC(obs_nh4no3.nb,as.data.frame(s4.nh4no3.nb))






# library(doBy)

# win.graph()
# par(mfrow=c(3,2))

# plot(NH4~time,data=as.data.frame(ss.nh4.fb))
# points(obs_nh4.fb[,'NH4']~obs_nh4.fb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(NH4_15~time,data=as.data.frame(ss.nh4.fb))
# points(obs_nh4.fb[,'NH4_15']~obs_nh4.fb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(NO3~time,data=as.data.frame(ss.nh4.fb))
# points(obs_nh4.fb[,'NO3']~obs_nh4.fb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(NO3_15~time,data=as.data.frame(ss.nh4.fb))
# points(obs_nh4.fb[,'NO3_15']~obs_nh4.fb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(N2O~time,data=as.data.frame(ss.nh4.fb))
# points(obs_nh4.fb[,'N2O']~obs_nh4.fb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(N2O_15~time,data=as.data.frame(ss.nh4.fb))
# points(obs_nh4.fb[,'N2O_15']~obs_nh4.fb[,'time'],type='b',pch=19,cex=2,col=2)


# win.graph()
# par(mfrow=c(3,2))

# plot(NH4~time,data=as.data.frame(ss.nh4.nb))
# points(obs_nh4.nb[,'NH4']~obs_nh4.nb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(NH4_15~time,data=as.data.frame(ss.nh4.nb))
# points(obs_nh4.nb[,'NH4_15']~obs_nh4.nb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(NO3~time,data=as.data.frame(ss.nh4.nb))
# points(obs_nh4.nb[,'NO3']~obs_nh4.nb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(NO3_15~time,data=as.data.frame(ss.nh4.nb))
# points(obs_nh4.nb[,'NO3_15']~obs_nh4.nb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(N2O~time,data=as.data.frame(ss.nh4.nb))
# points(obs_nh4.nb[,'N2O']~obs_nh4.nb[,'time'],type='b',pch=19,cex=2,col=2)
# plot(N2O_15~time,data=as.data.frame(ss.nh4.nb))
# points(obs_nh4.nb[,'N2O_15']~obs_nh4.nb[,'time'],type='b',pch=19,cex=2,col=2)


# win.graph()
# par(mfrow=c(3,2))

# plot(obs_nh4.fb[,'NH4']~obs_nh4.fb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.fb[,'NH4'],best.nh4.fb[,'NH4']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='NH4')
# points(best.nh4.fb[,'NH4']~best.nh4.fb[,'time'],pch='+',cex=2)
# plot(obs_nh4.fb[,'NH4_15']~obs_nh4.fb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.fb[,'NH4_15'],best.nh4.fb[,'NH4_15']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='NH4_15')
# points(best.nh4.fb[,'NH4_15']~best.nh4.fb[,'time'],pch='+',cex=2)
# plot(obs_nh4.fb[,'NO3']~obs_nh4.fb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.fb[,'NO3'],best.nh4.fb[,'NO3']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='NO3')
# points(best.nh4.fb[,'NO3']~best.nh4.fb[,'time'],pch='+',cex=2)
# plot(obs_nh4.fb[,'NO3_15']~obs_nh4.fb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.fb[,'NO3_15'],best.nh4.fb[,'NO3_15'])*c(0.95,1.05),na.rm=T)),xlab='Time',ylab='NO3_15')
# points(best.nh4.fb[,'NO3_15']~best.nh4.fb[,'time'],pch='+',cex=2)
# plot(obs_nh4.fb[,'N2O']~obs_nh4.fb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.fb[,'N2O'],best.nh4.fb[,'N2O']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='N2O')
# points(best.nh4.fb[,'N2O']~best.nh4.fb[,'time'],pch='+',cex=2)
# plot(obs_nh4.fb[,'N2O_15']~obs_nh4.fb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.fb[,'N2O_15'],best.nh4.fb[,'N2O_15']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='N2O_15')
# points(best.nh4.fb[,'N2O_15']~best.nh4.fb[,'time'],pch='+',cex=2)


# win.graph()
# par(mfrow=c(3,2))

# plot(obs_nh4.nb[,'NH4']~obs_nh4.nb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.nb[,'NH4'],best.nh4.nb[,'NH4']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='NH4')
# points(best.nh4.nb[,'NH4']~best.nh4.nb[,'time'],pch='+',cex=2)
# plot(obs_nh4.nb[,'NH4_15']~obs_nh4.nb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.nb[,'NH4_15'],best.nh4.nb[,'NH4_15']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='NH4_15')
# points(best.nh4.nb[,'NH4_15']~best.nh4.nb[,'time'],pch='+',cex=2)
# plot(obs_nh4.nb[,'NO3']~obs_nh4.nb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.nb[,'NO3'],best.nh4.nb[,'NO3']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='NO3')
# points(best.nh4.nb[,'NO3']~best.nh4.nb[,'time'],pch='+',cex=2)
# plot(obs_nh4.nb[,'NO3_15']~obs_nh4.nb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.nb[,'NO3_15'],best.nh4.nb[,'NO3_15'])*c(0.95,1.05),na.rm=T)),xlab='Time',ylab='NO3_15')
# points(best.nh4.nb[,'NO3_15']~best.nh4.nb[,'time'],pch='+',cex=2)
# plot(obs_nh4.nb[,'N2O']~obs_nh4.nb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.nb[,'N2O'],best.nh4.nb[,'N2O']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='N2O')
# points(best.nh4.nb[,'N2O']~best.nh4.nb[,'time'],pch='+',cex=2)
# plot(obs_nh4.nb[,'N2O_15']~obs_nh4.nb[,'time'],type='b',pch=19,cex=1,col=2,ylim=c(range(c(obs_nh4.nb[,'N2O_15'],best.nh4.nb[,'N2O_15']),na.rm=T))*c(0.95,1.05),xlab='Time',ylab='N2O_15')
# points(best.nh4.nb[,'N2O_15']~best.nh4.nb[,'time'],pch='+',cex=2)


