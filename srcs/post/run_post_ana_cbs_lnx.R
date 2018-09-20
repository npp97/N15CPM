library(compiler)
#path0<-"D://N15CS//run"
path0<-"//home//junhui//research//N15CS//run"
#path0<-"D://N15CS//run"
path0a<-"//home//junhui//research//N15CS//run//"
#path0a<-"D://N15CS//run"

setwd(path0)

#exp_code<-"_NE";cd_lth<-3;path1<-"out_cbs_2017-04-26_mm_ne_z2013am"
#exp_code<-"_SUSF_AVG";cd_lth<-2;path1<-"SUSF-6090-AVG-CS-ZHANG2013AM"
exp_code<-"_SUSF_60";cd_lth<-2;path1<-"SUSF-90-CS-ZHANG2013AM"

path_dl<-"/home/junhui/research/N15CS/run";
#path_dl<-"D://N15CS//run"
dll_name="zhang2013am.so"
#dll_name<-"zhang2013am.dll"
#trt_ssp=5;trt_esp=11; #FOR SUSF-6090-AVG-CS-ZHANG2013AM
trt_ssp=1;trt_esp=14;
repnum=10000;
model_Type='zhang2013am'

loadcmp(paste(path0,'//post_ana_cbs.Rc',sep=''))


#exp_code<-"_NE";cd_lth<-3;path1<-"out_cbs_2017-04-26_mm_ne_z2013am"
#exp_code<-"_SUSF_AVG";cd_lth<-2;path1<-"SUSF-6090-AVG-CS-ZHANG2013AM"
exp_code<-"_Naddtion_NECS";cd_lth<-2;path1<-"Naddition_ne__2018-01-09-ZHANG2013AM"

path_dl<-"/home/junhui/research/N15CS/run";
#path_dl<-"D://N15CS//run"
dll_name="zhang2013am.so";
#dll_name<-"zhang2013am.dll"
#trt_ssp=5;trt_esp=11; #FOR SUSF-6090-AVG-CS-ZHANG2013AM
trt_ssp=1;trt_esp=9;
repnum=10000;
model_Type='zhang2013am'

loadcmp(paste(path0,'//post_ana_cbs.Rc',sep=''))
