##Measure to detect collinearity in ERGM
#Based on Duxbury (2018) in Sociological Methods and Research


###====================
#=========ERGM========
#====================


###my.ergm is an ERGM object inherited from ergm package--available via statnet
VIF.ERGM<-function(my.ergm){
  require(ergm)

  cor.mat<-cov2cor(my.ergm$covar) #calculate correlation matrix
  corr5<-cor.mat[-c(1),-c(1)] ##omit edges term
  
  corr5<-corr5[!is.na(corr5[1:nrow(corr5)]),]
  corr5<-corr5[,which(!is.na(corr5[1,1:ncol(corr5)]))]
  
  VIFS<-matrix(0,nr=1,nc=ncol(corr5))
  
  for(i in 1:ncol(corr5)){
    
    gvec<-as.vector(corr5[-c(i),i]) ##create vector of correlations between covariate of interest and other covariates in the model
    tgvec<-t(gvec)            
    xcor<-solve(corr5[-c(i),-c(i)]) ##create square matrix of correlations between covariates in the model other than the one of interest
    Rsq<-tgvec%*%xcor%*%gvec
    VIFS[1,i]<-1/(1-Rsq)
  }
  ##Columns are covariates as they appear in the SAOM object
  colnames(VIFS)<-names(my.ergm$coef[-c(1)])
 message("Higher values indicate greater correlation.\nVIF > 20 is concerning, VIF > 100 indicates severe collinearity.")
  VIFS
}


VIF.ERGM(my.ergm)









#================================
###=========Example==============
#=================================

##Not run


##simulate network
library(igraph)
A<-erdos.renyi.game(75,.1,type="gnp")

##convert to sna object
a<-get.adjacency(A)
a<-as.matrix(a)

library(statnet)
a<-as.network(a,directed=FALSE)


##run ERGM
m1<- ergm( a ~ edges + kstar(2)+meandeg+
             gwesp(decay=.7, fixed=TRUE),
           control=control.ergm(seed=40))

##evaluate collinearity
VIF.ERGM(m1)



           
  

