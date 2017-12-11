###====================
#=========ERGM========
#====================


###my.ergm is an ERGM object inherited from ergm package--available via statnet
VIF.ERGM<-function(my.ergm){
  
  ###simulate from posterior distribution of ERGM--toggle nsim for more robustness/less computation time
  m2<-simulate(my.ergm,statsonly=TRUE,nsim=1000)
  
  
  cor.mat<-cor(m2) #calculate correlation matrix
  corr5<-cor.mat[-c(1),-c(1)] ##omit edges term
  VIFS<-matrix(0,nr=1,nc=ncol(corr5))
  
  for(i in 1:ncol(corr5)){
    
    gvec<-as.vector(corr5[-c(i),i]) ##create vector of correlations between covariate of interest and other covariates in the model
    tgvec<-t(gvec)            
    xcor<-solve(corr5[-c(i),-c(i)]) ##create square matrix of correlations between covariates in the model other than the one of interest
    Rsq<-tgvec%*%xcor%*%gvec
    VIFS[1,i]<-1/(1-Rsq)
  }
  ##Columns are covariates as they appear in the SAOM object
  
  VIFS
}

