library(FME)


parms_gen_lhs<-function (PARAMETERS, NUMSAMPLES, PMIN, PMAX, ALGORITHM) 
{
	NUMPARAMS <- length(PARAMETERS)
	ALGORITHM <- tolower(ALGORITHM)
	if (ALGORITHM == "optimum") {
		design <- lhs::optimumLHS(NUMSAMPLES, NUMPARAMS, 
				maxSweeps = 2, eps = 0.1)
	}
	else {
		design <- lhs::randomLHS(NUMSAMPLES, NUMPARAMS)
	}
	for (k in 1:NUMSAMPLES) {
		for (l in 1:NUMPARAMS) {
			lhc_max <- PMAX[l]
			lhc_min <- PMIN[l]
			value <- (design[k, l] * (lhc_max - lhc_min)) + lhc_min
			design[k, l] <- value
		}
	}
	colnames(design) <- c(PARAMETERS)
	return(design)
}


#define C_NRec_2_NH4    parms[0]
#define Vmax_NH4_2_NRec parms[1]
#define Km_NH4_2_NRec   parms[2]
#define Vmax_NH4_2_NLab parms[3]
#define Km_NH4_2_NLab   parms[4]
#define Vmax_NLab_2_NH4 parms[5]
#define Km_NLab_2_NH4   parms[6]
#define Vmax_NO3_2_NRec parms[7] 
#define Km_NO3_2_NRec   parms[8]
#define C_NRec_2_NO3    parms[9] 
#define K_NH4ads_2_NH4  parms[10] 
#define K_NH4_2_NH4ads  parms[11]
#define Vmax_NO3_2_NH4  parms[12] 
#define Km_NO3_2_NH4    parms[13]
#define Vmax_NH4_2_NO3  parms[14]
#define Km_NH4_2_NO3    parms[15]
#define K_NO3_2_NO3sto  parms[16]
#define K_NO3sto_2_NO3  parms[17]
#define K_NO3_2_N2O     parms[18]
#define K_NH4_2_N2O     parms[19]
#define K_NRec_2_N2O    parms[20]
#define alf_NH4_N2O  
#define alf_NO3_N2O
library(coda);

findBST<-function(flst){
  
    load(flst[1])

    bestf0.nh4<-MCMC.NH4$bestfunp;
    bestf0.no3<-MCMC.NO3$bestfunp;
#   bestf0.nh4no3<-MCMC.NH4NO3$bestfunp;

  for (i in 2:length(flst)){
    load(flst[i]);
    bestf0.nh4<-rbind(bestf0.nh4,MCMC.NH4$bestfunp);
    bestf0.no3<-rbind(bestf0.no3,MCMC.NO3$bestfunp);
#   bestf0.nh4no3<-rbind(bestf0.nh4no3,MCMC.NH4NO3$bestfunp);
  } 

    bestf0th.nh4<-quantile(bestf0.nh4,c(0.95));
    bestf0th.no3<-quantile(bestf0.no3,c(0.95));
#   bestf0th.nh4no3<-quantile(bestfZ13amg15cost_2l0.nh4no3,c(0.95));

    load(flst[1])

    ibest.nh4no3<-ibest.no3<-ibest.nh4<-1;
    ibestr.nh4no3<-ibestr.no3<-ibestr.nh4<-1;

    if (length(which(is.na(s4.nh4[,'NH4'])))>1){
        load(flst[2]);
        ibest.nh4no3<-ibest.no3<-ibest.nh4<-2;
        ibestr.nh4no3<-ibestr.no3<-ibestr.nh4<-2;
    }

    s4.nh4.fb<-(s4.nh4)#,s4.nh4)
    s4.no3.fb<-(s4.no3)#,s4.no3)
#   s4.nh4no3.fb<-(s4.nh4no3)#,s4.nh4no3)c


    par.nh4.fb<-(as.data.frame(MCMC.NH4$bestpar))
    par.no3.fb<-(as.data.frame(MCMC.NO3$bestpar))
#    par.nh4no3.fb<-(as.data.frame(MCMC.NH4NO3$bestpar))

    best.nh4<-s4.nh4
    bestf.nh4<-MCMC.NH4$bestfunp
    bestr.nh4<-s4.nh4
    bestfr.nh4<-sum(diag(cor_s4_nh4_obs[,-1])^2)

    best.no3<-s4.no3
    bestf.no3<-MCMC.NO3$bestfunp
    bestr.no3<-s4.no3
    bestfr.no3<-sum(diag(cor_s4_no3_obs[,-1])^2)

#    best.nh4no3<-ss.nh4no3
#    bestf.nh4no3<-MCMC.NH4NO3$bestfunp
#    bestr.nh4no3<-s4.nh4no3
#    bestfr.nh4no3<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)

    #obs_nh4.fb<-obs_nh4
    #obs_no3.fb<-obs_no3
    #obs_nh4no3.fb <- obs_nh4no3

for (i in 2:length(flst)){
    load(flst[i])
    
    if (length(which(is.na(s4.nh4[,'NH4'])))>1) next;
    if (length(which(is.na(c(diag(cor_s4_nh4_obs),diag(cor_s4_no3_obs)))))>1) next;
    if ((MCMC.NH4$bestfunp>bestf0th.nh4) | (MCMC.NO3$bestfunp > bestf0th.no3)) next;
    ibest.nh4no3<-ibest.no3<-ibest.nh4<-i;
    ibestr.nh4no3<-ibestr.no3<-ibestr.nh4<-i;

    s4.nh4.fb<-rbind(s4.nh4.fb,s4.nh4)#,s4.nh4)
    s4.no3.fb<-rbind(s4.no3.fb,s4.no3)#,s4.no3)
#   s4.nh4no3.fb<-rbind(s4.nh4no3.fb,s4.nh4no3)#,s4.nh4no3)
   
    par.nh4.fb<-cbind(par.nh4.fb,as.data.frame(MCMC.NH4$bestpar))
    par.no3.fb<-cbind(par.no3.fb,as.data.frame(MCMC.NO3$bestpar))
#   par.nh4no3.fb<-cbind(par.nh4no3.fb,as.data.frame(MCMC.NH4NO3$bestpar))

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
    # if(bestf.nh4no3>MCMC.NH4NO3$bestfunp){
        # best.nh4no3<-s4.nh4no3;
        # bestf.nh4no3<-MCMC.NH4NO3$bestfunp;
        # ibest.nh4no3<-i;
    # }
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
    # tmp<-sum(diag(cor_s4_nh4no3_obs[,-1])^2)
    # if(bestfr.nh4no3<tmp){
        # bestr.nh4no3<-s4.nh4no3;
        # bestfr.nh4no3<-tmp;
        # ibestr.nh4no3<-i;
    # }
}
    load(flst[ibest.nh4]);
        MCMC.NH4.bst<-MCMC.NH4;
    load(flst[ibest.no3]);
        MCMC.NO3.bst<-MCMC.NO3;
#    load(flst[ibest.nh4no3])
#        MCMC.NH4NO3.bst<-MCMC.NH4NO3;

    load(flst[ibestr.nh4]);
        MCMC.NH4.bstr<-MCMC.NH4;
    load(flst[ibestr.no3]);
        MCMC.NO3.bstr<-MCMC.NO3;
#    load(flst[ibestr.nh4no3])
#        MCMC.NH4NO3.bstr<-MCMC.NH4NO3;
    par.nh4=t(par.nh4.fb);par.no3=t(par.no3.fb);#par.nh4no3=t(par.nh4no3.fb);
    row.names(par.nh4)<-1:nrow(par.nh4);
    row.names(par.no3)<-1:nrow(par.no3);
#    row.names(par.nh4no3)<-1:nrow(par.nh4no3);
   
    par.asumm.nh4<-summary(as.mcmc(par.nh4));
    par.asumm.no3<-summary(as.mcmc(par.no3));
#    par.asumm.nh4no3<-summary(as.mcmc(par.nh4no3));    

    par.asumm.nh4<-cbind(as.data.frame(par.asumm.nh4$statistics),as.data.frame(par.asumm.nh4$quantiles))
    par.asumm.no3<-cbind(as.data.frame(par.asumm.no3$statistics),as.data.frame(par.asumm.no3$quantiles))
#    par.asumm.nh4no3<-cbind(as.data.frame(par.asumm.nh4no3$statistics),as.data.frame(par.asumm.nh4no3$quantiles))

    par.dsumm.nh4<-parDist(par.nh4)$qbs;
    par.dsumm.no3<-parDist(par.no3)$qbs;
#    par.dsumm.nh4no3<-parDist(par.nh4no3)$qbs;

    par.mcsumm.nh4<-parBBoot(par.nh4,nboot=4000)$qbs;
    par.mcsumm.no3<-parBBoot(par.no3,nboot=4000)$qbs;
#    par.mcsumm.nh4no3<-parBBoot(par.nh4no3,nboot=4000)$qbs;

    #group variables to return
    arst<-list(MCMC.NH4.bst=MCMC.NH4.bst,MCMC.NO3.bst=MCMC.NO3.bst,#MCMC.NH4NO3.bst=MCMC.NH4NO3.bst,
    MCMC.NH4.bstr=MCMC.NH4.bstr,MCMC.NO3.bstr=MCMC.NO3.bstr,#MCMC.NH4NO3.bstr=MCMC.NH4NO3.bstr,
    par.nh4=par.nh4,par.no3=par.no3,#par.nh4no3=par.nh4no3,
    par.asumm.nh4=par.asumm.nh4,par.asumm.no3=par.asumm.no3,#par.asumm.nh4no3=par.asumm.nh4no3,
    par.dumm.nh4=par.dsumm.nh4,par.dsumm.no3=par.dsumm.no3,#par.dsumm.nh4no3=par.dsumm.nh4no3,
    par.mcsumm.nh4=par.mcsumm.nh4,par.mcsumm.no3=par.mcsumm.no3,#par.mcsumm.nh4no3=par.mcsumm.nh4no3,
    obs_nh4=obs_nh4,obs_no3=obs_no3,#obs_nh4no3=obs_nh4no3,
    y0.nh4=y0.nh4,y0.no3=y0.no3,#y0.nh4no3=y0.nh4no3,
    bestf0.nh4=bestf0.nh4,bestf0.no3=bestf0.no3,#bestf0.nh4no3=bestf0.nh4no3,
    fn.best.nh4=flst[ibest.nh4],fn.best.no3=flst[ibest.no3],#fn.best.nh4no3=flst[ibest.nh4no3],
    fn.bestr.nh4=flst[ibestr.nh4],fn.bestr.no3=flst[ibestr.no3],#fn.bestr.nh4no3=flst[ibestr.nh4no3], 
    best.nh4=as.data.frame(best.nh4),best.no3=as.data.frame(best.no3),#best.nh4no3=as.data.frame(best.nh4no3),
    bestr.nh4=as.data.frame(bestr.nh4),bestr.no3=as.data.frame(bestr.no3),#bestr.nh4no3=as.data.frame(bestr.nh4no3),
    s4.nh4=as.data.frame(s4.nh4.fb),s4.no3=as.data.frame(s4.no3.fb)#,s4.nh4no3=as.data.frame(s4.nh4no3.fb)
	);

return(arst)    
}

rangEst<-function(parMat){
  require(bayesboot);
  parMat<-as.data.frame(parMat);
  q95<-q05<-as.data.frame(matrix(NA,nrow=1,ncol=ncol(parMat)));
  names(q95)<-names(q05)<-names(parMat);
  l.p<-ncol(parMat)
  for (i in 1:l.p){
    tmp<-bayesboot(data = parMat[,i], statistic = weighted.mean, use.weights = TRUE);
    stmp<-summary(tmp)$value[3:4];
    q05[i]<-stmp[1];
    q95[i]<-stmp[2];
  }
  return(rbind(lo=q05,upp=q95))
}


parBBoot<-function(parMat,nboot){
  require(bayesboot);
  parMat<-as.data.frame(parMat);
  bs<-as.data.frame(matrix(NA,nrow=nboot,ncol=ncol(parMat)));
  qbs<- bs<-as.data.frame(matrix(NA,nrow=9,ncol=ncol(parMat)));
  row.names(qbs)<-c("mean","sd","hdi.low","hdi.high","q2.5%","q25%","median","q75%","q97.5%")
  names(bs)<-names(qbs)<-names(parMat);
  l.p<-ncol(parMat)
  for (i in 1:l.p){
    tmp<-bayesboot(data = parMat[,i], R=nboot, weighted.mean, use.weights = TRUE);
    bs[1:nrow(tmp),i]<-as.data.frame(tmp);
    qbs[,i]<-as.numeric(summary(tmp)$value);
  }
  return(list(bs=bs,qbs=qbs))
}

myfindPeaks <-function (x, thresh=0.05, span=0.25, lspan=0.05, noisey=TRUE)
{
  n <- length(x)
  y <- x
  mu.y.loc <- y
  if(noisey)
  {
    mu.y.loc <- (x[1:(n-2)] + x[2:(n-1)] + x[3:n])/3
    mu.y.loc <- c(mu.y.loc[1], mu.y.loc, mu.y.loc[n-2])
  }
  y.loess <- loess(x~I(1:n), span=span)
  y <- y.loess[[2]]
  sig.y <- var(y.loess$resid, na.rm=TRUE)^0.5
  DX.1 <- sign(diff(mu.y.loc, na.pad = FALSE))
  pks <- which(diff(DX.1, na.pad = FALSE) < 0 & DX.1[-(n-1)] > 0) + 1
  out <- pks
  if(noisey)
  {
    n.w <- floor(lspan*n/2)
    out <- NULL
    for(pk in pks)
    {
      inner <- (pk-n.w):(pk+n.w)
      outer <- c((pk-2*n.w):(pk-n.w),(pk+2*n.w):(pk+n.w))
      mu.y.outer <- mean(y[outer])
      if(!is.na(mu.y.outer)) 
        if (mean(y[inner])-mu.y.outer > thresh*sig.y) out <- c(out, pk)
    }
  }
  out
}

Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }

  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

parDist<-function(parMat){
  require(fitdistrplus);
  parMat<-as.data.frame(parMat);
  bs<-as.data.frame(matrix(NA,nrow=10000,ncol=ncol(parMat)));
  qbs<- bs<-as.data.frame(matrix(NA,nrow=8,ncol=ncol(parMat)));
  row.names(qbs)<-c("mean","sd","peaks","q2.5%","q25%","median","q75%","q97.5%")
  names(bs)<-names(qbs)<-names(parMat);
  l.p<-ncol(parMat)
  nmdist<-c('norm','lnorm','gamma','weibull','cauchy')
  
  for (i in 1:l.p){
    lnormf<-try(fitdist(data = parMat[,i], "lnorm",method='mle'));
    gammaf<-try(fitdist(data = parMat[,i], "gamma", lower = c(0, 0),method='mle'));
    weibullf<-try(fitdist(data = parMat[,i], "weibull",lower = c(0, 0),method='mle'));
    normf<-try(fitdist(data = parMat[,i], "norm",method='mle'));
    cauchyf<-fitdist(data=parMat[,i],'cauchy')
  
    l.fdr<-list(normf,lnormf,gammaf,weibullf,cauchyf);
    ca<-0;
    for (k in 1:length(l.fdr)){if(class(l.fdr[[k]]) %in% 'try-error') ca<-c(ca,k);};ca<-ca[ca>0];
    if(length(ca)>0){l.fdr[ca]<-NULL;nmdist<-nmdist[-ca];}
    
    gst<-gofstat(l.fdr);  
    idist<-Mode(apply(as.data.frame(gst[c('ad','cvm','ks')]),2,which.min));
    
     p=c(0.025,0.25,0.5,0.75,0.975);
      if(nmdist[idist] %in% 'lnorm') {tmp=rlnorm(n=10000,meanlog=lnormf$estimate['meanlog'],sd=lnormf$estimate['sdlog']);
                      stmp=qlnorm(p=p,meanlog=lnormf$estimate['meanlog'],sd=lnormf$estimate['sdlog'])}
      if(nmdist[idist] %in% 'gamma') {tmp=rgamma(n=10000,shape=gammaf$estimate['shape'],rate=gammaf$estimate['rate']);
                      stmp<-qgamma(p=p,shape=gammaf$estimate['shape'],rate=gammaf$estimate['rate'])}
      if(nmdist[idist] %in% 'weibull') {tmp=rweibull(n=10000,shape=weibullf$estimate['shape'],scale=weibullf$estimate['scale']);
                      stmp=qweibull(p=p,shape=weibullf$estimate['shape'],scale=weibullf$estimate['scale']);}
      if(nmdist[idist] %in% 'norm') {tmp=rnorm(n=10000,mean=normf$estimate['mean'],sd=normf$estimate['sd']);
                      stmp=qnorm(p=p,mean=normf$estimate['mean'],sd=normf$estimate['sd'])}
      if(nmdist[idist] %in% 'cauchy') {tmp=rcauchy(n=10000,location=cauchyf$estimate['location'],scale=cauchyf$estimate['scale']);
                      stmp=qnorm(p=p,location=cauchyf$estimate['location'],scale=cauchyf$estimate['scale'])}
  

    bs[1:length(tmp),i]<-as.data.frame(tmp);
    dtmp<-density(tmp);mpk<-dtmp$x[which.max(dtmp$y)];
    qbs[1:8,i]<-c(mean(tmp),sd(tmp),mpk,stmp)
  }
  
  return(list(bs=bs,qbs=qbs))
}

rangEstr<-function(parMat){
  parMat<-as.data.frame(parMat);
  q95<-q05<-as.data.frame(matrix(NA,nrow=1,ncol=ncol(parMat)));
  names(q95)<-names(q05)<-names(parMat);
  l.p<-ncol(parMat)
  for (i in 1:l.p){
    stmp<-quantile(parMat[,i],c(0.05,0.95));
    q05[i]<-stmp[1];
    q95[i]<-stmp[2];
  }
  return(list(lo=q05,upp=q95))
}

valid_plotMCMC<-function(obs,mod_a){
  mod_a<-as.data.frame(mod_a);
  obs<-as.data.frame(obs);

  mod_a<-mod_a[,names(obs)];
  tn<-names(mod_a);
  nob<-ncol(obs);
  yl<-apply(rbind(mod_a,obs),FUN=range,2,na.rm=T)
  nm<-ncol(mod_a);
  n<-floor(sqrt(nm-1));
  m<-floor((nm-1)/(n-1));
  times<-obs[,'time'];
#  par0<-par();
#  win.graph();
  par(mfcol=c(3,2), mar=c(3,3,2,0.5))
  if(nrow(mod_a)>nrow(obs)){
         for (i in 2:(length(tn))){
            plot(x=mod_a[,'time'],y=mod_a[,i],type='l',lwd=2,pch='+',cex=2,xlab='',ylab='',ylim=c(yl[1,i]*0.9,yl[2,i]*1.1));title(tn[i]);
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
#        valid_plotMCMC(obs_nh4no3,best.nh4no3);
    })
    dev.off()

    pdf(paste('valid_bestr_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,bestr.nh4);
        valid_plotMCMC(obs_no3,bestr.no3);
#        valid_plotMCMC(obs_nh4no3,bestr.nh4no3);
    })
    dev.off()

    pdf(paste('valid_s4_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,s4.nh4);
        valid_plotMCMC(obs_no3,s4.no3);
#        valid_plotMCMC(obs_nh4no3,s4.nh4no3);
    })
    dev.off()

    pdf(paste('denMCMC_pars_bestr_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        densplot(as.mcmc(MCMC.NH4.bstr$pars[diff(MCMC.NH4.bstr$SS)<0,]))
        traceplot(as.mcmc(MCMC.NH4.bstr$pars[diff(MCMC.NH4.bstr$SS)<0,]))
        plot(as.mcmc(MCMC.NH4.bstr$pars[diff(MCMC.NH4.bstr$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NO3.bstr$pars[diff(MCMC.NO3.bstr$SS)<0,]),trace=T,smooth=T)  
#        plot(as.mcmc(MCMC.NH4NO3.bstr$pars[diff(MCMC.NH4NO3.bstr$SS)<0,]),trace=T,smooth=T)  
    })
    dev.off()

    pdf(paste('denMCMC_pars_best_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        densplot(as.mcmc(MCMC.NH4.bst$pars[diff(MCMC.NH4.bst$SS)<0,]))
        traceplot(as.mcmc(MCMC.NH4.bst$pars[diff(MCMC.NH4.bst$SS)<0,]))
        plot(as.mcmc(MCMC.NH4.bst$pars[diff(MCMC.NH4.bst$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NO3.bst$pars[diff(MCMC.NO3.bst$SS)<0,]),trace=T,smooth=T)  
#        plot(as.mcmc(MCMC.NH4NO3.bst$pars[diff(MCMC.NH4NO3.bst$SS)<0,]),trace=T,smooth=T)  
    })
    dev.off()

    pdf(paste('denMCMC_pars_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        anh4<-parBBoot(par.nh4,nboot=4000);
        ano3<-parBBoot(par.no3,nboot=4000);
#        anh4no3<-parBBoot(par.nh4no3,nboot=4000);
             
        densplot(as.mcmc(anh4$bs))
        plot(as.mcmc(anh4$bs),trace=T,smooth=T)
        densplot(as.mcmc(ano3$bs))
        plot(as.mcmc(ano3$bs),trace=T,smooth=T)
#        densplot(as.mcmc(anh4no3$bs))
#        plot(as.mcmc(anh4no3$bs),trace=T,smooth=T)
    })
    dev.off()
}


gen_fig_final<-function(arst,char_trt){

    pdf(paste('valid_final_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,s4.nh4);
        valid_plotMCMC(obs_no3,s4.no3);
#        valid_plotMCMC(obs_nh4no3,s4.nh4no3);
    })
    dev.off()

    pdf(paste('denMCMC_pars_final_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        densplot(as.mcmc(MCMC.NH4.fin$pars[diff(MCMC.NH4.fin$SS)<0,]))
        traceplot(as.mcmc(MCMC.NH4.fin$pars[diff(MCMC.NH4.fin$SS)<0,]))
        plot(as.mcmc(MCMC.NH4.fin$pars[diff(MCMC.NH4.fin$SS)<0,]),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NO3.fin$pars[diff(MCMC.NO3.fin$SS)<0,]),trace=T,smooth=T)  
#        plot(as.mcmc(MCMC.NH4NO3.fin$pars[diff(MCMC.NH4NO3.fin$SS)<0,]),trace=T,smooth=T)  
    })
    dev.off()

    
    pdf(paste('denMCMC_pars_final_all_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        densplot(as.mcmc(MCMC.NH4.fin$pars))
        traceplot(as.mcmc(MCMC.NH4.fin$pars))
        plot(as.mcmc(MCMC.NH4.fin$pars),trace=T,smooth=T)  
        plot(as.mcmc(MCMC.NO3.fin$pars),trace=T,smooth=T)  
#        plot(as.mcmc(MCMC.NH4NO3.fin$pars),trace=T,smooth=T)  
    })
    dev.off()

    pdf(paste('denMCMC_pars_mcmc_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        anh4<-parBBoot(MCMC.NH4.fin$pars,nboot=4000);
        ano3<-parBBoot(MCMC.NO3.fin$pars,nboot=4000);
#        anh4no3<-parBBoot(MCMC.NH4NO3.fin$pars,nboot=4000);
           
        densplot(as.mcmc(anh4$bs))
        plot(as.mcmc(anh4$bs),trace=T,smooth=T)
        densplot(as.mcmc(ano3$bs))
        plot(as.mcmc(ano3$bs),trace=T,smooth=T)
#        densplot(as.mcmc(anh4no3$bs))
#        plot(as.mcmc(anh4no3$bs),trace=T,smooth=T)
    })
    dev.off()

    pdf(paste('denMCMC_pars_dist_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        anh4<-parDist(MCMC.NH4.fin$pars);
        ano3<-parDist(MCMC.NO3.fin$pars);
#        anh4no3<-parDist(MCMC.NH4NO3.fin$pars);
           
        densplot(as.mcmc(anh4$bs))
        plot(as.mcmc(anh4$bs),trace=T,smooth=T)
        densplot(as.mcmc(ano3$bs))
        plot(as.mcmc(ano3$bs),trace=T,smooth=T)
#        densplot(as.mcmc(anh4no3$bs))
#        plot(as.mcmc(anh4no3$bs),trace=T,smooth=T)
    })
    dev.off()
 
}


gen_fig2<-function(arst,char_trt){

    pdf(paste('valid_final_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,bestr.nh4);
        valid_plotMCMC(obs_no3,bestr.no3);
#        valid_plotMCMC(obs_nh4no3,bestr.nh4no3);
    })
    dev.off()

       
    pdf(paste('denMCMC_pars_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        anh4<-parBBoot(par.nh4,nboot=4000);
        ano3<-parBBoot(par.no3,nboot=4000);
#        anh4no3<-parBBoot(par.nh4no3,nboot=4000);
           
        densplot(as.mcmc(anh4$bs))
        plot(as.mcmc(anh4$bs),trace=T,smooth=T)
        densplot(as.mcmc(ano3$bs))
        plot(as.mcmc(ano3$bs),trace=T,smooth=T)
#        densplot(as.mcmc(anh4no3$bs))
#        plot(as.mcmc(anh4no3$bs),trace=T,smooth=T)
    })
    dev.off()
}


gen_fig_final2<-function(arst,char_trt){

    pdf(paste('valid_final_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        valid_plotMCMC(obs_nh4,s4.nh4);
        valid_plotMCMC(obs_no3,s4.no3);
#        valid_plotMCMC(obs_nh4no3,s4.nh4no3);
    })
    dev.off()

       
    pdf(paste('denMCMC_pars_',char_trt,'.pdf',sep=''),onefile=T)
    with(arst,{
        anh4<-parBBoot(MCMC.NH4.fin$pars,nboot=4000);
        ano3<-parBBoot(MCMC.NO3.fin$pars,nboot=4000);
#        anh4no3<-parBBoot(MCMC.NH4NO3.fin$pars,nboot=4000);
           
        densplot(as.mcmc(anh4$bs))
        plot(as.mcmc(anh4$bs),trace=T,smooth=T)
        densplot(as.mcmc(ano3$bs))
        plot(as.mcmc(ano3$bs),trace=T,smooth=T)
#        densplot(as.mcmc(anh4no3$bs))
#        plot(as.mcmc(anh4no3$bs),trace=T,smooth=T)
    })
    dev.off()
}

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
      rng0<-rangEstr(arst$par.nh4)
      lu<-apply(rbind(rng0$lo,rng0$upp,pp),2,range)
      upper=lu[2,];
      lower=lu[1,];
  MCMC.NH4<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=arst$obs_nh4,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));


      fn_n_seed<-fn_n_rseed(arst$fn.bestr.no3,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NO3.bstr$bestpar;
      rng0<-rangEstr(arst$par.no3)
      lu<-apply(rbind(rng0$lo,rng0$upp,pp),2,range)
      upper=lu[2,];
      lower=lu[1,];

  MCMC.NO3<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=arst$obs_no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

      # fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4no3,fntrt_p);
      # set.seed(fn_n_seed$fn_rseed);
      # pp<-arst$MCMC.NH4NO3.bstr$bestpar
     
      # rng0<-rangEstr(arst$par.nh4no3)
      # lu<-apply(rbind(rng0$lo,rng0$upp,pp),2,range)
      # upper=lu[2,];
      # lower=lu[1,];
  # MCMC.NH4NO3<-try(mdMCMC2a(ff=Z13amg15cost,pp1=pp,obs=arst$obs_nh4no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

	  # s4.nh4no3 <- Z13amg15(MCMC.NH4NO3$bestpar,arst$y0.nh4no3, times)
	  s4.nh4<-Z13amg15(MCMC.NH4$bestpar,arst$y0.nh4,times)
	  s4.no3<-Z13amg15(MCMC.NO3$bestpar,arst$y0.no3,times)
	  
	  cor_s4_nh4_obs<-cor_z13amg15(arst$obs_nh4,MCMC.NH4$bestpar,arst$y0.nh4,times)
	  cor_s4_no3_obs<-cor_z13amg15(arst$obs_no3,MCMC.NO3$bestpar,arst$y0.no3,times)
#	  cor_s4_nh4no3_obs<-cor_z13amg15(arst$obs_nh4no3,MCMC.NH4NO3$bestpar,arst$y0.nh4no3,times)
	 
      par.summ.nh4<-summary(as.mcmc(MCMC.NH4$pars[diff(MCMC.NH4$SS)<0,]));
      par.summ.no3<-summary(as.mcmc(MCMC.NO3$pars[diff(MCMC.NO3$SS)<0,]));
#      par.summ.nh4no3<-summary(as.mcmc(MCMC.NH4NO3$pars[diff(MCMC.NH4NO3$SS)<0,]));    
   
      par.summ.nh4<-cbind(as.data.frame(par.summ.nh4$statistics),as.data.frame(par.summ.nh4$quantiles))
      par.summ.no3<-cbind(as.data.frame(par.summ.no3$statistics),as.data.frame(par.summ.no3$quantiles))
#      par.summ.nh4no3<-cbind(as.data.frame(par.summ.nh4no3$statistics),as.data.frame(par.summ.nh4no3$quantiles))

      par.asumm.nh4<-summary(as.mcmc(MCMC.NH4$pars));
      par.asumm.no3<-summary(as.mcmc(MCMC.NO3$pars));
#      par.asumm.nh4no3<-summary(as.mcmc(MCMC.NH4NO3$pars));    

      par.asumm.nh4<-cbind(as.data.frame(par.asumm.nh4$statistics),as.data.frame(par.asumm.nh4$quantiles))
      par.asumm.no3<-cbind(as.data.frame(par.asumm.no3$statistics),as.data.frame(par.asumm.no3$quantiles))
#      par.asumm.nh4no3<-cbind(as.data.frame(par.asumm.nh4no3$statistics),as.data.frame(par.asumm.nh4no3$quantiles))

      par.dsumm.nh4<-parDist(MCMC.NH4$pars)$qbs;
      par.dsumm.no3<-parDist(MCMC.NO3$pars)$qbs;
#      par.dsumm.nh4no3<-parDist(MCMC.NH4NO3$pars)$qbs;

      par.mcsumm.nh4<-parBBoot(MCMC.NH4$pars,nboot=4000)$qbs;
      par.mcsumm.no3<-parBBoot(MCMC.NO3$pars,nboot=4000)$qbs;
#      par.mcsumm.nh4no3<-parBBoot(MCMC.NH4NO3$pars,nboot=4000)$qbs;


     var2sav<-list(arst=arst,
                MCMC.NH4.fin=MCMC.NH4,MCMC.NO3.fin=MCMC.NO3,#MCMC.NH4NO3.fin=MCMC.NH4NO3,
                par.summ.nh4=par.summ.nh4,par.summ.no3=par.summ.no3,#par.summ.nh4no3=par.summ.nh4no3,
                par.asumm.nh4=par.asumm.nh4,par.asumm.no3=par.asumm.no3,#par.asumm.nh4no3=par.asumm.nh4no3,
                par.dumm.nh4=par.dsumm.nh4,par.dsumm.no3=par.dsumm.no3,#par.dsumm.nh4no3=par.dsumm.nh4no3,
                par.mcsumm.nh4=par.mcsumm.nh4,par.mcsumm.no3=par.mcsumm.no3,#par.mcsumm.nh4no3=par.mcsumm.nh4no3,
                s4.nh4=s4.nh4,s4.no3=s4.no3,#s4.nh4no3=s4.nh4no3,
                obs_nh4=arst$obs_nh4,obs_no3=arst$obs_no3,#obs_nh4no3=arst$obs_nh4no3,
                y0.nh4=arst$y0.nh4,y0.no3=arst$y0.no3,#y0.nh4no3=arst$y0.nh4no3, 
	            cor_s4_nh4_obs=cor_s4_nh4_obs,cor_s4_no3_obs=cor_s4_no3_obs#,cor_s4_nh4no3_obs=cor_s4_nh4no3_obs
				)
     dyn.unload(lib_nm);
     return(var2sav)
  }
 
 
  getpar_fnal_z13am<-function(arst,fntrt_p=3,lib_nm,lfactor=10,ufactor=10,niter_1=50000,niter_2=50000,burninlength=10000){
      dyn.load(lib_nm);
      times=arst$obs_nh4[,'time'];
      fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NH4.bstr$bestpar
      rng0<-rangEstr(arst$par.nh4)
      lu<-apply(rbind(rng0$lo,rng0$upp,pp),2,range)
      upper=lu[2,];
      lower=lu[1,];
      MCMC.NH4<-try(mdMCMC2a(ff=Z13amcostL,pp1=pp,obs=arst$obs_nh4,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));


      fn_n_seed<-fn_n_rseed(arst$fn.bestr.no3,fntrt_p);
      set.seed(fn_n_seed$fn_rseed);
      pp<-arst$MCMC.NO3.bstr$bestpar
      rng0<-rangEstr(arst$par.no3)
      lu<-apply(rbind(rng0$lo,rng0$upp,pp),2,range)
      upper=lu[2,];
      lower=lu[1,];
     MCMC.NO3<-try(mdMCMC2a(ff=Z13amcostL,pp1=pp,obs=arst$obs_no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

      # fn_n_seed<-fn_n_rseed(arst$fn.bestr.nh4no3,fntrt_p);
      # set.seed(fn_n_seed$fn_rseed);
      # pp<-arst$MCMC.NH4NO3.bstr$bestpar
      # rng0<-rangEstr(arst$par.nh4no3)
      # lu<-apply(rbind(rng0$lo,rng0$upp,pp),2,range)
      # upper=lu[2,];
      # lower=lu[1,];
    # MCMC.NH4NO3<-try(mdMCMC2a(ff=Z13amcostL,pp1=pp,obs=arst$obs_nh4no3,niter1=niter_1,niter2=niter_2,y0=arst$y0.nh4no3,times,wvar0=0.1,updatecov=100,lower=lower,upper=upper,burninlength=burnin.length));

#	  s4.nh4no3 <- Z13am(arst$y0.nh4no3,MCMC.NH4NO3$bestpar, times)
	  s4.nh4<-Z13am(arst$y0.nh4,MCMC.NH4$bestpar,times)
	  s4.no3<-Z13am(arst$y0.no3,MCMC.NO3$bestpar,times)
	  
	  cor_s4_nh4_obs<-cor_z13am(arst$obs_nh4,MCMC.NH4$bestpar,arst$y0.nh4,times)
	  cor_s4_no3_obs<-cor_z13am(arst$obs_no3,MCMC.NO3$bestpar,arst$y0.no3,times)
#	  cor_s4_nh4no3_obs<-cor_z13am(arst$obs_nh4no3,MCMC.NH4NO3$bestpar,arst$y0.nh4no3,times)
	  
      par.summ.nh4<-summary(as.mcmc(MCMC.NH4$pars[diff(MCMC.NH4$SS)<0,]));
      par.summ.no3<-summary(as.mcmc(MCMC.NO3$pars[diff(MCMC.NO3$SS)<0,]));
#      par.summ.nh4no3<-summary(as.mcmc(MCMC.NH4NO3$pars[diff(MCMC.NH4NO3$SS)<0,]));    

      par.summ.nh4<-cbind(as.data.frame(par.summ.nh4$statistics),as.data.frame(par.summ.nh4$quantiles))
      par.summ.no3<-cbind(as.data.frame(par.summ.no3$statistics),as.data.frame(par.summ.no3$quantiles))
#      par.summ.nh4no3<-cbind(as.data.frame(par.summ.nh4no3$statistics),as.data.frame(par.summ.nh4no3$quantiles))  

      par.asumm.nh4<-summary(as.mcmc(MCMC.NH4$pars));
      par.asumm.no3<-summary(as.mcmc(MCMC.NO3$pars));
#      par.asumm.nh4no3<-summary(as.mcmc(MCMC.NH4NO3$pars));    

      par.asumm.nh4<-cbind(as.data.frame(par.asumm.nh4$statistics),as.data.frame(par.asumm.nh4$quantiles))
      par.asumm.no3<-cbind(as.data.frame(par.asumm.no3$statistics),as.data.frame(par.asumm.no3$quantiles))
#      par.asumm.nh4no3<-cbind(as.data.frame(par.asumm.nh4no3$statistics),as.data.frame(par.asumm.nh4no3$quantiles))

      par.dsumm.nh4<-parDist(MCMC.NH4$pars)$qbs;
      par.dsumm.no3<-parDist(MCMC.NO3$pars)$qbs;
#      par.dsumm.nh4no3<-parDist(MCMC.NH4NO3$pars)$qbs;

      par.mcsumm.nh4<-parBBoot(MCMC.NH4$pars,nboot=4000)$qbs;
      par.mcsumm.no3<-parBBoot(MCMC.NO3$pars,nboot=4000)$qbs;
#      par.mcsumm.nh4no3<-parBBoot(MCMC.NH4NO3$pars,nboot=4000)$qbs;

 	  var2sav<-list(arst=arst,
                MCMC.NH4.fin=MCMC.NH4,MCMC.NO3.fin=MCMC.NO3,#MCMC.NH4NO3.fin=MCMC.NH4NO3,
                par.summ.nh4=par.summ.nh4,par.summ.no3=par.summ.no3,#par.summ.nh4no3=par.summ.nh4no3,
                par.asumm.nh4=par.asumm.nh4,par.asumm.no3=par.asumm.no3,#par.asumm.nh4no3=par.asumm.nh4no3,
                par.dumm.nh4=par.dsumm.nh4,par.dsumm.no3=par.dsumm.no3,#par.dsumm.nh4no3=par.dsumm.nh4no3,
                par.mcsumm.nh4=par.mcsumm.nh4,par.mcsumm.no3=par.mcsumm.no3,#par.mcsumm.nh4no3=par.mcsumm.nh4no3, 
                s4.nh4=s4.nh4,s4.no3=s4.no3,#s4.nh4no3=s4.nh4no3,
                obs_nh4=arst$obs_nh4,obs_no3=arst$obs_no3,#obs_nh4no3=arst$obs_nh4no3,
                y0.nh4=arst$y0.nh4,y0.no3=arst$y0.no3,#y0.nh4no3=arst$y0.nh4no3, 
	            cor_s4_nh4_obs=cor_s4_nh4_obs,cor_s4_no3_obs=cor_s4_no3_obs#,cor_s4_nh4no3_obs=cor_s4_nh4no3_obs
				)
     dyn.unload(lib_nm);
     return(var2sav)
  }
 
get_n2oc15<-function(time0,n2o,n2o_15){
  a2<-approx(x=time0,y=n2o,xout=seq(min(time0),max(time0)));
  a2$y=cumsum(a2$y);
  cn2o<-approx(x=a2$x,y=a2$y,xout=time0)$y;
  a3<-approx(x=time0,y=n2o*n2o_15/100,xout=seq(min(time0),max(time0)));
  a3$y<-cumsum(a3$y);
  cn2o_15<-approx(x=a3$x,y=a3$y,xout=time0)$y;
  return(cbind(N2O=cn2o,N2O_15=cn2o_15));
}


get_n2oc<-function(time0,n2o){
  a2<-approx(x=time0,y=n2o,xout=seq(min(time0),max(time0)));
  a2$y=cumsum(a2$y);
  cn2o<-approx(x=a2$x,y=a2$y,xout=time0)$y;
  return(cbind(N2O=cn2o));
}


get_n2oca<-function(times,n2o,n2o_15){
  a2<-approx(x=times,y=n2o,xout=seq(min(times),max(times),1));
  a2$y=cumsum(a2$y);
  cn2o<-approx(x=a2$x,y=a2$y,xout=times)$y;
  yy=n2o*n2o_15/100;
  a3<-approx(x=times,y=yy,xout=seq(min(times),max(times),1),yleft=yy[1],yright=yy[length(n2o)]);
  a3$y<-cumsum(a3$y);
  cn2o_15<-approx(x=a3$x,y=a3$y,xout=times)$y;
  return(cbind(N2O=cn2o,N2O_15=cn2o_15));
}


Z13amg15 <- function (pars,y0,times) {
  times0<-times;
  times=seq(min(times),max(times),by=1);
 
  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4);
 
  # N2O.out<-out[,'NH4_2_N2O']+out[,'NO3_2_N2O']+out[,'NRec_2_N2O'];
  # N2O_15.out <- out[,'NH4_2_N2O_15'] + out[,'NO3_2_N2O_15'] + out[,'NRec_2_N2O_15'];
  out <-out[times %in% times0,];
 
  # out<-cbind(out,N2O=N2O.out,N2O_15=N2O_15.out)
  return(out)
}


Z13amg <- function (pars,y0,times) {
  times0<-times;
  times=seq(min(times),max(times),by=1);
  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4);
 
  # N2O.out<-out[,'NH4_2_N2O']+out[,'NO3_2_N2O']+out[,'NRec_2_N2O'];
  # N2O_15.out <- out[,'NH4_2_N2O_15'] + out[,'NO3_2_N2O_15'] + out[,'NRec_2_N2O_15'];

  out <-out[times %in% times0,];
  # out<-cbind(out,N2O=N2O.out,N2O_15=N2O_15.out)
  return(out)
}



Z13am15 <- function (pars,y0,times) {
  times0<-times;
  times=seq(min(times),max(times),by=1);
  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013am',nout=4);
  out <-out[times %in% times0,];

  # N2O.out<-out[,'NH4_2_N2O']+out[,'NO3_2_N2O']+out[,'NRec_2_N2O'];
  # N2O_15.out <- out[,'NH4_2_N2O_15'] + out[,'NO3_2_N2O_15'] + out[,'NRec_2_N2O_15'];

  # out<-cbind(out,N2O=N2O.out,N2O_15=N2O_15.out)
  return(out)
}


Z13am <- function (pars,y0,times) {

  times0<-times;
  times=seq(min(times),max(times),by=1);

  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013am',nout=4);
  out <-out[times %in% times0,];
 
  # N2O.out<-out[,'NH4_2_N2O']+out[,'NO3_2_N2O']+out[,'NRec_2_N2O'];
  # N2O_15.out <- out[,'NH4_2_N2O_15'] + out[,'NO3_2_N2O_15'] + out[,'NRec_2_N2O_15'];

  # out<-cbind(out,N2O=N2O.out,N2O_15=N2O_15.out)
  return(out)
}

cor_z13amg15<-function(obs,parms,y0,times){
 
  out<-Z13amg15(parms,y0,times);

  return(cor(obs,out[,names(obs)]));
}


cor_z13amg<-function(obs,parms,y0,times){
  
  out<-Z13amg(parms,y0,times);

  return(cor(obs,out[,names(obs)]));
}

cor_z13am15<-function(obs,parms,y0,times){
  
  out<-Z13am15(parms,y0,times);

  return(cor(obs,out[,names(obs)]));
}


cor_z13am<-function(obs,parms,y0,times){
  
  out<-Z13am(y0,parms,times);

  return(cor(obs,out[,names(obs)]));
}

M14<-function(pars,y0,times){
  times0<-times;
  times=seq(min(times),max(times),by=1);
  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='ode45',dllname='mull2014',nout=9);
  out <-out[times %in% times0,];
  return(out)
}

M14cost <- function (pars_v,y0,pars0,pars_fix,times,obs) {
  iipar_fix<-which(names(pars0) %in% names(pars_fix));
  pars0[iipar_fix]<-pars_fix;
  iipar_v<-which(names(pars0) %in% names(pars_v));
  pars0[iipar_v]<-pars_v;
  
  out<-M14(pars0,y0,times);
  cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F);
  return(cost_1);
}

M14costL <- function (pars_v,y0,pars0,pars_fix,times,obs) {
  iipar_fix<-which(names(pars0) %in% names(pars_fix));
  pars0[iipar_fix]<-pars_fix;
  iipar_v<-which(names(pars0) %in% names(pars_v));
  pars0[iipar_v]<-pars_v;
  
  out<-M14(pars0,y0,times);  
  cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)$minlogp;
  return(cost_1);
}


Z13amAcost <- function (p=pars_v,y0,pars0,pars_fix,times,obs) {
  iipar_fix<-which(names(pars0) %in% names(pars_fix))
  pars0[iipar_fix]<-pars_fix;
  iipar_v<-which(names(pars0) %in% names(p));
  pars0[iipar_v]<-p;
  
  out <- as.data.frame(ode(y=y0,parms=pars0,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_1 <- modCost(model=out, obs=obs, weight='std',scaleVar =F);
  return(cost_1);
}

Z13amAcostL <- function (p=pars_v,y0,pars0,pars_fix,times,obs) {
  iipar_fix<-which(names(pars0) %in% names(pars_fix));
  pars0[iipar_fix]<-pars_fix;
  iipar_v<-which(names(pars0) %in% names(p));
  pars0[iipar_v]<-p;
  
  out <- as.data.frame(ode(y=y0,parms=pars0,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_1 <- modCost(model=out, obs=obs, weight='std',scaleVar =F)$minlogp;
  return(cost_1);
}


Z13amAcostV <- function (p=pars_v,y0,pars0,pars_fix,times,obs) {
  iipar_fix<-which(names(pars0) %in% names(pars_fix));
  pars0[iipar_fix]<-pars_fix;
  iipar_v<-which(names(pars0) %in% names(p));
  pars0[iipar_v]<-p;
  
  out1 <- as.data.frame(ode(y=y0,parms=pars0,times=times,method='ode45',func='derivs',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_1 <- modCost(model=out1, obs=obs, weight='std',scaleVar =F)$model;
  return(cost_1)
}


mdMCMC14<-function(ff=M14cost,pars_v,pars0,pars_fix,obs,niter1=35000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
  MCMC1 <- modMCMC(f=ff, p=pars_v, y=y0,pars0=pars0,pars_fix=pars_fix,niter=niter1,times=times,obs=obs,ntrydr = 3,
                   wvar0=0.1,updatecov=updatecov,lower=lower, upper=upper,burninlength=burninlength);
  if(class(MCMC1) %in% "modMCMC"){
	Nt<-as.numeric(raftery.diag(as.mcmc(MCMC1$pars))$resmatrix[2]);
	if((MCMC1$naccepted<=100)|(niter1<Nt)){
		niter=max(c(Nt,50000));
		MCMC1 <- modMCMC(f=ff,p=MCMC1$bestpar,y=y0,pars0=pars0,pars_fix=pars_fix,niter=Nt,times=times,obs=obs,ntrydr = 3,
                     wvar0=0.1, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
	}
  }else
  {
	  MCMC1 <- modMCMC(f=ff, p=pars_v, y=y0,pars0=pars0,pars_fix=pars_fix,niter=niter1,times=times,obs=obs,ntrydr = 3,
			  wvar0=0.1,updatecov=updatecov/10,lower=lower, upper=upper,burninlength=burninlength);
	  if(class(MCMC1) %in% "modMCMC"){
		  Nt<-as.numeric(raftery.diag(as.mcmc(MCMC1$pars))$resmatrix[2]);
		  if((MCMC1$naccepted<=100)|(niter1<Nt)){
			  niter=max(c(Nt,50000));
			  MCMC1 <- modMCMC(f=ff,p=MCMC1$bestpar,y=y0,pars0=pars0,pars_fix=pars_fix,niter=Nt,times=times,obs=obs,ntrydr = 3,
					  wvar0=0.1, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
		  }
	  }
	}
  return(MCMC1);
}


M14costB <- function (pars,y0,times,obs) {
	out<-M14(pars,y0,times);
	cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)
	return(cost_1)
}

M14costBL <- function (pars,y0,times,obs) {
	out<-M14(pars,y0,times);  
	cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)$minlogp
	return(cost_1)
}

M14costBV <- function (pars,y0,times,obs) {
	out<-M14(pars,y0,times);
	cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)$model;
	return(cost_1)
}

mdMCMC14B<-function(ff=M14costBL,pars,obs,niter1=35000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
	MCMC1 <- modMCMC(f=ff, p=pars, y=y0,niter=niter1,times=times,obs=obs,ntrydr = 3,
			wvar0=0.1,updatecov=updatecov,lower=lower, upper=upper,burninlength=burninlength);

	if(class(MCMC1) %in% "modMCMC"){
		Nt<-as.numeric(raftery.diag(as.mcmc(MCMC1$pars))$resmatrix[2]);
		if((MCMC1$naccepted<=100)|(niter1<Nt)){
			niter=max(c(Nt,50000));
			MCMC1 <- modMCMC(f=ff,p=MCMC1$bestpar,y=y0,niter=niter,times=times,obs=obs,ntrydr = 3,
					wvar0=0.1, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
		}
	}
	return(MCMC1);
}

mdMCMCZ13am<-function(ff=Z13amAcost,pars_v,pars0,pars_fix,obs,niter1=35000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
  MCMC1 <- modMCMC(f=ff,p=pars_v,y=y0,pars0=pars0,pars_fix=pars_fix,niter=niter1,times=times,obs=obs,ntrydr = 3,
                   wvar0=0.1,updatecov=updatecov,lower=lower, upper=upper,burninlength=burninlength);

 if(class(MCMC1) %in% "modMCMC"){
	Nt<-as.numeric(raftery.diag(as.mcmc(MCMC1$pars))$resmatrix[2]);
	
	if((MCMC1$naccepted<=100)|(niter1<Nt)){
		niter=max(c(Nt,50000));
		MCMC1 <- modMCMC(f=ff,p=MCMC1$bestpar,y=y0,pars0=pars0,pars_fix=pars_fix,niter=niter,times=times,obs=obs,ntrydr = 3,
                     wvar0=0.1, updatecov=updatecov,lower=lower, upper=upper,burninlength=10000);
	}
  }
  
  return(MCMC1);
}


Z13amcostL <- function (pars,y0,times,obs) {
	times0<-times;
	times=seq(min(times),max(times),by=1);
	
	out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='ode45',dllname='zhang2013am',nout=4);
	cost_1 <- modCost(model=out[times %in% times0,], obs=obs, weight='std',scaleVar =F)$minlogp
	return(cost_1)
}

Z13amcostV <- function (pars,y0,times,obs) {
	times0<-times;
	times=seq(min(times),max(times),by=1);
	
	out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='ode45',dllname='zhang2013am',nout=4);
	cost_1 <- modCost(model=out[times %in% times0,], obs=obs, weight='std',scaleVar =F)$model
	return(cost_1)
}
Z13am<-function(y0,pars,times){
	times0<-times;
	times=seq(min(times),max(times),by=1);
	
	out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='ode45',dllname='zhang2013am',nout=4);
	return(out[times %in% times0,]);
}


mdMCMCZ13amB<-function(ff=Z13amcostL,pars,y0,obs,niter1=35000,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
	MCMC1 <- modMCMC(f=ff,p=pars,y=y0,niter=niter1,times=times,obs=obs,ntrydr = 3,
			wvar0=0.1,updatecov=updatecov,lower=lower, upper=upper,burninlength=burninlength);

	if(class(MCMC1) %in% "modMCMC"){
		Nt<-as.numeric(raftery.diag(as.mcmc(MCMC1$pars))$resmatrix[2]);
		if((MCMC1$naccepted<=100)|(niter1<Nt)){
			MCMC1 <- modMCMC(f=ff,p=MCMC1$bestpar,y=y0,niter=Nt,times=times,obs=obs,ntrydr = 3,
					wvar0=0.1, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
		}
	}
	return(MCMC1);
}

Z13amgcost <- function (pars,y0,times,obs) {
  out<-Z13amg15(pars,y0,times);
  # obs$NH4_2_N2O = pars['alf_NH4_N2O']*obs$N2O;
  # obs$NO3_2_N2O = pars['alf_NO3_N2O']*obs$N2O;
  # obs$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O;
  
#  obs$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs$N2O_15;
#  obs$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs$N2O_15;
#  obs$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O_15;
   
  cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)
  return(cost_1)
}

Z13amgcostL <- function (pars,y0,times,obs) {
  #out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',dllname='zhang2013amg',nout=4);
  out<-Z13amg15(pars,y0,times);  
# obs$NH4_2_N2O = pars['alf_NH4_N2O']*obs$N2O; obs$NO3_2_N2O = pars['alf_NO3_N2O']*obs$N2O;obs$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O;
#  obs$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs$N2O_15;obs$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs$N2O_15;obs$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O_15;
 
  cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)$minlogp
  return(cost_1)
}



Z13amgcost_3l <- function (pars,y0_nh4,y0_no3,y0_n4n3,times,obs_nh4,obs_no3,obs_n4n3) {

  out_nh4 <- as.data.frame(Z13amg15(pars,y0_nh4,times));
		#obs_nh4$NH4_2_N2O = pars['alf_NH4_N2O']*obs_nh4$N2O; obs_nh4$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_nh4$N2O; obs_nh4$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O;
#		obs_nh4$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_nh4$N2O_15;obs_nh4$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_nh4$N2O_15;obs_nh4$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O_15;
  cost_nh4 <- modCost(model=out_nh4[,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(Z13amg15(pars,y0_no3,times));
#		obs_no3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_no3$N2O; obs_no3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_no3$N2O; obs_no3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O;
#		obs_no3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_no3$N2O_15;obs_no3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_no3$N2O_15;obs_no3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O_15;
  cost_no3 <- modCost(model=out_no3[,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  out_n4n3 <-  as.data.frame(Z13amg15(pars,y0_n4n3,times));
#		obs_n4n3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_n4n3$N2O; obs_n4n3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_n4n3$N2O; obs_n4n3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O;
#		obs_n4n3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_n4n3$N2O_15;obs_n4n3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_n4n3$N2O_15;obs_n4n3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O_15;
  cost_n4n3 <- modCost(model=out_n4n3[,names(obs_n4n3)], obs=obs_n4n3,weight='std',scaleVar =F,cost=cost_no3);
  cost_1<-cost_n4n3
  return(cost_1)
}


Z13amgcost_3ll <- function (pars,y0_nh4,y0_no3,y0_n4n3,times,obs_nh4,obs_no3,obs_n4n3) { 
 
  out_nh4 <- as.data.frame(Z13amg15(pars,y0_nh4,times));
		#obs_nh4$NH4_2_N2O = pars['alf_NH4_N2O']*obs_nh4$N2O; obs_nh4$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_nh4$N2O; obs_nh4$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O;
#		obs_nh4$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_nh4$N2O_15;obs_nh4$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_nh4$N2O_15;obs_nh4$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O_15;
  cost_nh4 <- modCost(model=out_nh4[,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(Z13amg15(pars,y0_no3,times));
#		obs_no3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_no3$N2O; obs_no3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_no3$N2O; obs_no3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O;
#		obs_no3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_no3$N2O_15;obs_no3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_no3$N2O_15;obs_no3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O_15;
  cost_no3 <- modCost(model=out_no3[,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  out_n4n3 <-  as.data.frame(Z13amg15(pars,y0_n4n3,times));
#		obs_n4n3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_n4n3$N2O; obs_n4n3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_n4n3$N2O; obs_n4n3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O;
#		obs_n4n3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_n4n3$N2O_15;obs_n4n3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_n4n3$N2O_15;obs_n4n3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O_15;
  cost_n4n3 <- modCost(model=out_n4n3[,names(obs_n4n3)], obs=obs_n4n3,weight='std',scaleVar =F,cost=cost_no3);
  cost_1<-cost_n4n3
  return(cost_1$minlogp)
}


Z13amgcost_2l <- function (pars,y0_nh4,y0_no3,times,obs_nh4,obs_no3) {

  out_nh4 <- as.data.frame(Z13amg15(pars,y0_nh4,times));
  cost_nh4 <- modCost(model=out_nh4[,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(Z13amg15(pars,y0_no3,times));
  cost_no3 <- modCost(model=out_no3[,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  cost_1<-cost_no3
  return(cost_1)
}


Z13amgcost_2ll <- function (pars,y0_nh4,y0_no3,times,obs_nh4,obs_no3) { 
 
  out_nh4 <- as.data.frame(Z13amg15(pars,y0_nh4,times));
  cost_nh4 <- modCost(model=out_nh4[,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(Z13amg15(pars,y0_no3,times));
  cost_no3 <- modCost(model=out_no3[,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  cost_1<-cost_no3
  return(cost_1$minlogp)
}


Z13amg15cost <- function (pars,y0,times,obs) {
  times0<-times;
  times=seq(min(times),max(times),by=1);

  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4);
  out<-out[times %in% times0,];
  # obs$NH4_2_N2O = pars['alf_NH4_N2O']*obs$N2O;
  # obs$NO3_2_N2O = pars['alf_NO3_N2O']*obs$N2O;
  # obs$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O;
  
  # obs$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs$N2O_15;
  # obs$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs$N2O_15;
  # obs$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O_15;

  # out[,'N2O'] = out[,'NRec_2_N2O']+out[,'NH4_2_N2O']+out[,'NO3_2_N2O']
  # out[,'N2O_15'] = out[,'NRec_2_N2O_15']+out[,'NH4_2_N2O_15']+out[,'NO3_2_N2O_15']
 
  cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)
  return(cost_1)
}

Z13amg15costL <- function (pars,y0,times,obs) {
  times0<-times;
  times=seq(min(times),max(times),by=1);

  out <- ode(y=y0,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4);
  out<-out[times %in% times0,];
  # obs$NH4_2_N2O = pars['alf_NH4_N2O']*obs$N2O; obs$NO3_2_N2O = pars['alf_NO3_N2O']*obs$N2O;obs$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O;
  # obs$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs$N2O_15;obs$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs$N2O_15;obs$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs$N2O_15;
 
 #  out=cbind(out,'N2O'= out[,'NRec_2_N2O']+out[,'NH4_2_N2O']+out[,'NO3_2_N2O'])
 #  out=cbind(out,'N2O_15' = out[,'NRec_2_N2O_15']+out[,'NH4_2_N2O_15']+out[,'NO3_2_N2O_15'])

  cost_1 <- modCost(model=out[,names(obs)], obs=obs, weight='std',scaleVar =F)$minlogp
  return(cost_1)
}


Z13amg15cost_3l <- function (pars,y0_nh4,y0_no3,y0_n4n3,times,obs_nh4,obs_no3,obs_n4n3) {
  times0<-times;
  times=seq(min(times),max(times),by=1);

  out_nh4 <- as.data.frame(ode(y=y0_nh4,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
#		obs_nh4$NH4_2_N2O = pars['alf_NH4_N2O']*obs_nh4$N2O; obs_nh4$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_nh4$N2O; obs_nh4$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O;
#		obs_nh4$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_nh4$N2O_15;obs_nh4$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_nh4$N2O_15;obs_nh4$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O_15;
    # out_nh4=cbind(out_nh4,'N2O'= out_nh4[,'NRec_2_N2O']+out_nh4[,'NH4_2_N2O']+out_nh4[,'NO3_2_N2O'])
    # out_nh4=cbind(out_nh4,'N2O_15' = out_nh4[,'NRec_2_N2O_15']+out_nh4[,'NH4_2_N2O_15']+out_nh4[,'NO3_2_N2O_15'])

  cost_nh4 <- modCost(model=out_nh4[times %in% times0,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(ode(y=y0_no3,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
#		obs_no3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_no3$N2O; obs_no3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_no3$N2O; obs_no3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O;
#		obs_no3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_no3$N2O_15;obs_no3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_no3$N2O_15;obs_no3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O_15;
    # out_no3=cbind(out_no3,'N2O' = out_no3[,'NRec_2_N2O']+out_no3[,'NH4_2_N2O']+out_no3[,'NO3_2_N2O'])
    # out_no3=cbind(out_no3,'N2O_15' = out_no3[,'NRec_2_N2O_15']+out_no3[,'NH4_2_N2O_15']+out_no3[,'NO3_2_N2O_15'])
  cost_no3 <- modCost(model=out_no3[times %in% times0,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  out_n4n3 <- as.data.frame(ode(y=y0_n4n3,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
#		obs_n4n3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_n4n3$N2O; obs_n4n3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_n4n3$N2O; obs_n4n3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O;
#		obs_n4n3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_n4n3$N2O_15;obs_n4n3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_n4n3$N2O_15;obs_n4n3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O_15;
    # out_n4n3=cbind(out_n4n3,'N2O' = out_n4n3[,'NRec_2_N2O']+out_n4n3[,'NH4_2_N2O']+out_n4n3[,'NO3_2_N2O'])
    # out_n4n3=cbind(out_n4n3,'N2O_15' = out_n4n3[,'NRec_2_N2O_15']+out_n4n3[,'NH4_2_N2O_15']+out_n4n3[,'NO3_2_N2O_15'])

  cost_n4n3 <- modCost(model=out_n4n3[times %in% times0,names(obs_n4n3)], obs=obs_n4n3,weight='std',scaleVar =F,cost=cost_no3);
#  cost_1<-(cost_nh4$model+cost_no3$model+cost_n4n3$model);
  cost_1<-cost_n4n3
  return(cost_1)
}


Z13amg15cost_3ll <- function (pars,y0_nh4,y0_no3,y0_n4n3,times,obs_nh4,obs_no3,obs_n4n3) { 
  times0<-times;
  times=seq(min(times),max(times),by=1);
 
  out_nh4 <- as.data.frame(ode(y=y0_nh4,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
#		obs_nh4$NH4_2_N2O = pars['alf_NH4_N2O']*obs_nh4$N2O; obs_nh4$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_nh4$N2O; obs_nh4$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O;
#		obs_nh4$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_nh4$N2O_15;obs_nh4$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_nh4$N2O_15;obs_nh4$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_nh4$N2O_15;
    # out_nh4=cbind(out_nh4,'N2O'= out_nh4[,'NRec_2_N2O']+out_nh4[,'NH4_2_N2O']+out_nh4[,'NO3_2_N2O'])
    # out_nh4=cbind(out_nh4,'N2O_15' = out_nh4[,'NRec_2_N2O_15']+out_nh4[,'NH4_2_N2O_15']+out_nh4[,'NO3_2_N2O_15'])

  cost_nh4 <- modCost(model=out_nh4[times %in% times0,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(ode(y=y0_no3,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
#		obs_no3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_no3$N2O; obs_no3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_no3$N2O; obs_no3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O;
#		obs_no3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_no3$N2O_15;obs_no3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_no3$N2O_15;obs_no3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_no3$N2O_15;
    # out_no3=cbind(out_no3,'N2O' = out_no3[,'NRec_2_N2O']+out_no3[,'NH4_2_N2O']+out_no3[,'NO3_2_N2O'])
    # out_no3=cbind(out_no3,'N2O_15' = out_no3[,'NRec_2_N2O_15']+out_no3[,'NH4_2_N2O_15']+out_no3[,'NO3_2_N2O_15'])
  cost_no3 <- modCost(model=out_no3[times %in% times0,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  out_n4n3 <- as.data.frame(ode(y=y0_n4n3,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
#		obs_n4n3$NH4_2_N2O = pars['alf_NH4_N2O']*obs_n4n3$N2O; obs_n4n3$NO3_2_N2O = pars['alf_NO3_N2O'] * obs_n4n3$N2O; obs_n4n3$NRec_2_N2O = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O;
#		obs_n4n3$NH4_2_N2O_15 = pars['alf_NH4_N2O']*obs_n4n3$N2O_15;obs_n4n3$NO3_2_N2O_15 = pars['alf_NO3_N2O']*obs_n4n3$N2O_15;obs_n4n3$NRec_2_N2O_15 = max(0,(1-pars['alf_NH4_N2O']-pars['alf_NO3_N2O']))*obs_n4n3$N2O_15;
    # out_n4n3=cbind(out_n4n3,'N2O' = out_n4n3[,'NRec_2_N2O']+out_n4n3[,'NH4_2_N2O']+out_n4n3[,'NO3_2_N2O'])
    # out_n4n3=cbind(out_n4n3,'N2O_15' = out_n4n3[,'NRec_2_N2O_15']+out_n4n3[,'NH4_2_N2O_15']+out_n4n3[,'NO3_2_N2O_15'])

  cost_n4n3 <- modCost(model=out_n4n3[times %in% times0,names(obs_n4n3)], obs=obs_n4n3,weight='std',scaleVar =F,cost=cost_no3);
#
  cost_1<-cost_n4n3
  return(cost_1$minlogp)
}


Z13amg15cost_2l <- function (pars,y0_nh4,y0_no3,times,obs_nh4,obs_no3) {
  times0<-times;
  times=seq(min(times),max(times),by=1);

  out_nh4 <- as.data.frame(ode(y=y0_nh4,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
  cost_nh4 <- modCost(model=out_nh4[times %in% times0,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(ode(y=y0_no3,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
  cost_no3 <- modCost(model=out_no3[times %in% times0,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);
  cost_1<-cost_no3
  return(cost_1)
}


Z13amg15cost_2ll <- function (pars,y0_nh4,y0_no3,times,obs_nh4,obs_no3) { 
  times0<-times;
  times=seq(min(times),max(times),by=1);
 
  out_nh4 <- as.data.frame(ode(y=y0_nh4,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
  cost_nh4 <- modCost(model=out_nh4[times %in% times0,names(obs_nh4)], obs=obs_nh4, weight='std',scaleVar =F)

  out_no3 <- as.data.frame(ode(y=y0_no3,parms=pars,times=times,func='derivs',initfunc='initmod',method='lsoda',dllname='zhang2013amg',nout=4));
  cost_no3 <- modCost(model=out_no3[times %in% times0,names(obs_no3)], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);

  cost_1<-cost_no3
  return(cost_1$minlogp)
}


Z13amcost_3l <- function (pars,y0_nh4,y0_no3,y0_n4n3,times,obs_nh4,obs_no3,obs_n4n3) {
  times0<-times;
  times=seq(min(times),max(times),by=1);

  out_nh4 <- as.data.frame(ode(y=y0_nh4,parms=pars,times=times,func='derivs',method='lsoda',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_nh4 <- modCost(model=out_nh4[times %in% times0,], obs=obs_nh4, weight='std',scaleVar =F)
  out_no3 <- as.data.frame(ode(y=y0_no3,parms=pars,times=times,func='derivs',method='lsoda',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_no3 <- modCost(model=out_no3[times %in% times0,], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);
  out_n4n3 <- as.data.frame(ode(y=y0_n4n3,parms=pars,times=times,func='derivs',method='lsoda',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_n4n3 <- modCost(model=out_n4n3[times %in% times0,], obs=obs_n4n3,weight='std',scaleVar =F,cost=cost_no3);
#  cost_1<-(cost_nh4$model+cost_no3$model+cost_n4n3$model);
  cost_1<-cost_n4n3
  return(cost_1)
}


Z13amcost_3ll <- function (pars,y0_nh4,y0_no3,y0_n4n3,times,obs_nh4,obs_no3,obs_n4n3) {
  out_nh4 <- as.data.frame(ode(y=y0_nh4,parms=pars,times=times,func='derivs',method='lsoda',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_nh4 <- modCost(model=out_nh4[times %in% times0,], obs=obs_nh4, weight='std',scaleVar =F)
  out_no3 <- as.data.frame(ode(y=y0_no3,parms=pars,times=times,func='derivs',method='lsoda',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_no3 <- modCost(model=out_no3[times %in% times0,], obs=obs_no3, weight='std',scaleVar =F,cost=cost_nh4);
  out_n4n3 <- as.data.frame(ode(y=y0_n4n3,parms=pars,times=times,func='derivs',method='lsoda',initfunc='initmod',dllname='zhang2013am',nout=4));
  cost_n4n3 <- modCost(model=out_n4n3[times %in% times0,], obs=obs_n4n3,weight='std',scaleVar =F,cost=cost_no3);
  #  cost_1<-(cost_nh4$model+cost_no3$model+cost_n4n3$model);
  cost_1<-cost_n4n3
  return(cost_1$minlogp)
}

mdMCMC<-function(ff=M07cost,pp,obs,niter=35000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
  Fit <- modFit(f=ff,p=pp,y=y0,times=times,obs=obs,lower=lower,upper=upper);
  pp <- Fit$par;
  
  MCMC1 <- modMCMC(f=ff, p=pp, niter=15000, y=y0,times=times,obs=obs,ntrydr = 3,
                   wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
  
  Cov0<-cov(MCMC1$pars)*2.4^2/(length(pp));
  #lower=(summary(MCMC1)[5,])*0.5;
  upper=summary(MCMC1)[7,];
  pp=unlist(summary(MCMC1)[6,]);
  MCMC <- modMCMC(f=ff, p=MCMC1$bestpar, niter=niter, y=y0,times=times,obs=obs,jump=Cov0,ntrydr = 3,
                  wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper, burninlength=burninlength);
  
  return(MCMC);
}


mdMCMC1<-function(ff=M07cost,pp1,obs,niter1=35000,niter2=50000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
  Fit <- modFit(f=ff,p=pp1,y=y0,times=times,obs=obs,lower=lower,upper=upper,method="Pseudo");
  pp <- Fit$par;
  
  MCMC1 <- modMCMC(f=ff, p=pp, niter=niter1, y=y0,times=times,obs=obs,ntrydr = 3,
                   wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
  
  Cov0<-cov(MCMC1$pars)*2.4^2/(length(pp));
 # lower=(summary(MCMC1)[5,])*0.5;
  upper=summary(MCMC1)[7,];
  pp=unlist(summary(MCMC1)[6,]);
  MCMC <-modMCMC(f=ff, p=MCMC1$bestpar, niter=niter2, y=y0,times=times,obs=obs,jump=Cov0,ntrydr = 3,
                  wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper, burninlength=burninlength);
  if(MCMC$naccepted<=10){
    Fit <- modFit(f=ff,p=MCMC$bestpar,y=y0,times=times,obs=obs,lower=lower,upper=upper,method="Pseudo",control=list(numiter=50000));
    pp <- Fit$par;
    
    MCMC <- modMCMC(f=ff, p=pp, niter=niter2, y=y0,times=times,obs=obs,ntrydr = 3,
                     wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
  }
  return(MCMC);
}

mdMCMC2<-function(ff=M07cost,pp1,obs,niter1=35000,niter2=50000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
  MCMC1 <- modMCMC(f=ff, p=pp1, niter=niter1, y=y0,times=times,obs=obs,ntrydr = 3,
                   wvar0=0,updatecov=updatecov,lower=lower, upper=upper,burninlength=burninlength);
  if(MCMC1$naccepted<=10){
   
    MCMC1 <- modMCMC(f=ff, p=MCMC1$bestpar, niter=niter2, y=y0,times=times,obs=obs,ntrydr = 3,
                    wvar0=0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
  }
  Cov0<-cov(MCMC1$pars)*2.4^2/(length(pp));
  # lower=(summary(MCMC1)[5,])*0.5;
  pp=unlist(summary(MCMC1)[6,]);
  MCMC <-modMCMC(f=ff, p=pp, niter=niter2, y=y0,times=times,obs=obs,jump=Cov0,ntrydr = 3,
                 wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper, burninlength=burninlength);
  if(MCMC$naccepted<=10){
   MCMC <- modMCMC(f=ff, p=MCMC$bestpar, niter=niter2*2, y=y0,times=times,obs=obs,ntrydr = 3,var0=0.1,
                    wvar0=wvar0, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
  }
  return(MCMC);
}

mdMCMC2a<-function(ff=M07cost,pp1,pars0,pars.fix,obs,niter1=35000,niter2=50000,y0=y0,times,wvar0=0.1,updatecov=100,lower,upper,burninlength=10000)
{
  MCMC1 <- modMCMC(f=ff, p=pp1, niter=niter1, y=y0,times=times,obs=obs,ntrydr = 3,
                   wvar0=0.1,updatecov=updatecov,lower=lower, upper=upper,burninlength=burninlength);
  if(class(MCMC1) %in% "modMCMC"){
	 Nt<-as.numeric(raftery.diag(as.mcmc(MCMC1$pars))$resmatrix[2]);
	 if((MCMC1$naccepted<=100)|(niter1<Nt)){
		MCMC1 <- modMCMC(f=ff, p=MCMC1$bestpar, niter=niter2, y=y0,times=times,obs=obs,ntrydr = 3,
                     wvar0=0.1, updatecov=updatecov,lower=lower, upper=upper,burninlength=1000);
	}
  }
  return(MCMC1);
}
