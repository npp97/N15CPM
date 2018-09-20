require(compiler)
path0<-"D://N15CSg_zl//srcs"
#path0<-"/home/junhui/research/N15CS/srcs"
setwd(path0)

 system('rm ./model/*.o')
  system('rm ./model/*.dll')
  system('rm ./batch-p/*.Rc')
  system('rm ./subs/*.Rc ')
  system('rm ./subs/*.so')
  

  system('R CMD SHLIB ./model/zhang2013amg.c');
   
  cmpfile('./subs/subs_2gas.R');
  
  cmpfile('./batch-p/Z13amg15_batch_zl.R');
  
  cmpfile('./post/post_ana_cbs.R')

  system('rm ./model/*.o')
  system('mv ./model/*.dll ../run/')
  system('mv ./model/*.so ../run/')
  system('mv ./batch-p/*.Rc ../run/')
  system('mv ./subs/*.Rc ../run/')
  system('mv ./post/*.Rc ../run/')
  

