require(compiler)
path0<-"D://workdir//N15CPM//srcs"
#path0<-"/home/junhui/research/N15CS/srcs"
setwd(path0)

 system('rm ./model/*.o')
  system('rm ./model/*.dll')
  system('rm ./batch-p/*.Rc')
  system('rm ./subs/*.Rc ')
  system('rm ./subs/*.so')
  


 # system('R CMD SHLIB ./model/muller2003.c');
  #system('R CMD SHLIB ./model/muller2003f.c'); 
  #system('R CMD SHLIB ./model/zhang2013.c');
  system('R CMD SHLIB ./model/zhang2013a.c');
  #system('R CMD SHLIB ./model/zhang2013e.c');
  #system('R CMD SHLIB ./model/zhang2013am.c');
  
  #system('R CMD SHLIB ./model/muller2003a.c');
  #system('R CMD SHLIB ./model/muller2003b.c');  
  #system('R CMD SHLIB ./model/muller2003c.c');  
  #system('R CMD SHLIB ./model/muller2003_45.c');
  #system('R CMD SHLIB ./model/zhu2016.c');  
  
  cmpfile('./subs/subs_all.R');
  
  #cmpfile('./batch-p/muller2003m_batch_p2nc.R');
  #cmpfile('./batch-p/muller2003_batch_CS.R');
  #cmpfile('./batch-p/muller2003_batch_CS1.R');  
  #cmpfile('./batch-p/muller2003_batch_p.R');
  cmpfile('./batch-p/Z13a_batch_CS1.R');
  
  #cmpfile('./batch-p/zhang2013a_batch_p.r');
  #cmpfile('./batch-p/zhang2013_batch_p.r');
  #cmpfile('./batch-p/muller2003a_batch_p.R');
  #cmpfile('./batch-p/muller2003b_batch_p.R');
  #cmpfile('./batch-p/muller2003c_batch_p.R');
  
  #cmpfile('./batch-p/muller2003_batch_p1.R');
  #cmpfile('./batch-p/muller2003f_batch_p1.R');
  #cmpfile('./batch-p/zhang2013a_batch_p2nc.r');
  #cmpfile('./batch-p/zhang2013am_batch_p2nc.R');

  #cmpfile('./batch-p/zhang2013am_batch_CS1.R');
  cmpfile('./post/post_ana_cbs.R')
  #cmpfile('./batch-p/muller2003_45_batch_p1.R');
  #cmpfile('./batch-p/zhu2016_batch_p2.r');

  system('rm ./model/*.o')
  system('mv ./model/*.dll ../run/')
  system('mv ./model/*.so ../run/')
    system('mv ./batch-p/*.Rc ../run/')
  system('mv ./subs/*.Rc ../run/')
  system('mv ./post/*.Rc ../run/')
  

