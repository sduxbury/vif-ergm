#############Study one: VIFs per network configuration

###########

library(igraph)
library(statnet)
library(doParallel, quiet=TRUE)


####FUNCTIONS#####

####convert igraph object to sna network--igraph object generated in sim.set() function below
netgen<-function(n){
  
  n<-get.adjacency(n)
  n<-as.matrix(n)
  n<-as.network(n,directed=FALSE)}

###assign node-level attributes
set.att<-function(x){
  attr1<-c(rnorm((numnodes/3),2),rnorm((numnodes/3),0),rnorm((numnodes/3),-2))
  attr2<-c(rnorm((numnodes/3),2),rnorm((numnodes/3),0),rnorm((numnodes/3),-2))
  attr3<-c(rnorm((numnodes/3),2),rnorm((numnodes/3),0),rnorm((numnodes/3),-2))
  
  set.vertex.attribute(x,"attr1",attr1)
  set.vertex.attribute(x,"attr2",attr2)
  set.vertex.attribute(x,"attr3",attr3)}

###simulate random network, convert to sna, assign node-level attributes
sim.set<-function(){
  A<-erdos.renyi.game(75, 0.1, type="gnp")
  a<-netgen(A)
  a<-set.att(a)
}



###run ERGM, simulate distribution of networks calculate correlation matrix
ergm.func<-function(a){
  
  cl<-detectCores()##detects cores for parallel processing--optional
  
  m1<- ergm( a ~ edges + nodecov('attr1') + nodecov('attr2') + nodecov('attr3') +
               degcor,
             control=control.ergm(parallel=cl, parallel.type="PSOCK", 
                                  force.main=TRUE,main.method="MCMLE",MCMLE.density.guard.min=10000))
  m2<-simulate(m1,statsonly=TRUE,nsim=1000)
  corr5<-cor(m2)}


###combines functions to simulate network, assign node level attributes, run an ERGM, and calculate the correlation matrix
big.func<-function(){
  n<-sim.set()
  n<-try(ergm.func(n))}


###calculate VIF for the first variable
vif1<-function(x){
  gvec<-x[3:5,2]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[3:5,3:5]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)}

# VIF for the 2nd attribute
vif2<-function(x){
  gvec<-x[c(2,4:5),3]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,4:5),c(2,4:5)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)}

# VIF for the 3rd attribute
vif3<-function(x){
  gvec<-x[c(2,3,5),4]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,5),c(2,3,5)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)}

# VIF for the 4th term, gwesp
vif4<-function(x){gvec<-x[c(2:4),5]
gvec<-as.vector(gvec)
tgvec<-t(gvec)
xcorr<-x[c(2:4),c(2:4)]
xcor<-solve(xcorr)
Rsq<-tgvec%*%xcor%*%gvec
swvifs<-1/(1-Rsq)}






# random network, density=0.1, size=75, vmr=.5

numnodes<-75
density<-.1

setwd("X:/VIFS")

numnodes<-75
density<-.1
swvifs<-matrix(0,nr=1000,nc=4)


for (i in 1:1000){
  
  repeat{
    repeat{
      repeat{
        
        corr5<-big.func()
        if(!is(corr5, "try-error") & all(!is.na(corr5))) break   
        
      }
      
      test1<-try(solve(corr5[c(2:4),c(2:4)]))
      test2<-try(solve(corr5[c(3:5),c(3:5)]))
      
      if(!is(test1, "try-error") & !is(test2, "try-error")) break}
    
    swvifs[i,1]<-try(vif1(corr5)) 
    swvifs[i,2]<-try(vif2(corr5))
    swvifs[i,3]<-try(vif3(corr5))
    swvifs[i,4]<-try(vif4(corr5))
    
    ####ensures that VIFs will yield meaningful values, i.e. that all values in the correlation matrix are between -1 and 1
      ##Correlations may be >1 or <-1 when the likelihood function is completely unsolvable
    if(all(!is.na(swvifs[i,])) & all(swvifs[i,]>=1) & !is(swvifs[i,], "try-error")) break
  }
  
  save(swvifs, file='ExampleVIF.csv')
}



























########Study two: Variance change in ERGM reruns

##########
library(igraph)
library(statnet)
library(doParallel, quiet=TRUE)

####FUNCTIONS#####

##convert igraph network to sna
netgen<-function(n){
  
  n<-get.adjacency(n)
  n<-as.matrix(n)
  n<-as.network(n,directed=FALSE)}

####assign node-level variables to network
set.att<-function(x){
  attr1<-c(rnorm((numnodes/3),2),rnorm((numnodes/3),0),rnorm((numnodes/3),-2))
  attr2<-c(rnorm((numnodes/3),2),rnorm((numnodes/3),0),rnorm((numnodes/3),-2))
  attr3<-c(rnorm((numnodes/3),2),rnorm((numnodes/3),0),rnorm((numnodes/3),-2))
  
  set.vertex.attribute(x,"attr1",attr1)
  set.vertex.attribute(x,"attr2",attr2)
  set.vertex.attribute(x,"attr3",attr3)}

###simulate network, convert network, assign node-level attributes
sim.set<-function(){
  A<-erdos.renyi.game(75,0.1,type="gnp")
  a<-netgen(A)
  a<-set.att(a)
}

###run ERGM on a network
ergm.func<-function(a){
  
  cl<-detectCores()###necessary for parallel processing
  
  m1<- ergm( a ~ edges + nodecov('attr1') + nodecov('attr2') + nodecov('attr3') +
               gwesp(decay=.7, fixed=TRUE)+degreepopularity+degcor,
             control=control.ergm(parallel=cl, parallel.type="PSOCK", 
                                  force.main=TRUE,main.method="MCMLE",MCMLE.density.guard.min=10000))}

###generate correlation matrix
corr.gen<-function(z){
  m2<-simulate(z,statsonly=TRUE,nsim=1000)
  corr5<-cor(m2)}


###calculate VIFs
vif1<-function(x){
  gvec<-x[3:7,2]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[3:7,3:7]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)}

vif2<-function(x){
  gvec<-x[c(2,4:7),3]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,4:7),c(2,4:7)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)
}


vif3<-function(x){
  gvec<-x[c(2,3,5:7),4]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,5:7),c(2,3,5:7)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)
  
}


vif4<-function(x){
  gvec<-x[c(2:4,6,7),5]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,4,6,7),c(2,3,4,6,7)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)
  
}

vif5<-function(x){
  gvec<-x[c(2:5,7),6]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2:5,7),c(2:5,7)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)
}

vif6<-function(x){
  gvec<-x[c(2:6),7]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2:6),c(2:6)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  swvifs<-1/(1-Rsq)
}



###return VIFs in a matrix
vifmat<-function(x){
  vifx<-matrix(0, nr=1, nc=6)
  vifx[1,1]<-vif1(x)
  vifx[1,2]<-vif2(x)
  vifx[1,3]<-vif3(x)
  vifx[1,4]<-vif4(x)
  vifx[1,5]<-vif5(x)
  vifx[1,6]<-vif6(x)
  vifx
}
##################################















########################density=0.1##############


mega<-1
numnodes<-75


for (i in 1:500){
  repeat{
    repeat{
      
      repeat{
        a<-sim.set()
        m1<-try(ergm.func(a))
        
        if(!is(m1, "try-error")) break
      }        
      repeat{
        m2<-try(ergm.func(a))
        if(!is(m2, "try-error")) break
        
      }
      
      m1corr<-corr.gen(m1)
      m2corr<-corr.gen(m2)
      
      test1<-try(solve(m1corr[c(2:6),c(2:6)]))
      test2<-try(solve(m1corr[c(3:7),c(3:7)]))
      
      test3<-try(solve(m2corr[c(2:6),c(2:6)]))
      test4<-try(solve(m2corr[c(3:7),c(3:7)]))
      
      if(!is(test1, "try-error") & !is(test2, "try-error") & 
         !is(test3, "try-error") & !is(test4, "try-error")) break}
    
    vifs1<-vifmat(m1corr)
    vifs2<-vifmat(m2corr)
    
    if(all(!is.na(vifs1)) & all(vifs1>=1) & !is(vifs1, "try-error") &
       all(!is.na(vifs2)) & all(vifs2>=1) & !is(vifs2, "try-error")) break}
  
  loopid<-as.vector(rep(i+1, each=6))
  dens<-as.vector(rep(0.1, each=6))
  
  
  vifschange<-(vifs2-vifs1)/abs(vifs1)
  vifschange<-as.vector(vifschange)
  
  
  vifs1<-as.vector(t(vifs1))
  vifs2<-as.vector(t(vifs2))
  lengthv<-length(vifs1)
  
  se1<-sqrt(diag(m1$covar)) ###standard errors for model
  se1<-se1[2:7]
  
  se2<-sqrt(diag(m2$covar))
  se2<-se2[2:7]
  
  sechange<-(se2-se1)/abs(se1)
  sechange<-as.vector(sechange)
  
  varchange<-(se2^2-se1^2)/abs(se1^2)##change in variance
  
  se1<-as.vector(se1)
  se2<-as.vector(se2)
  
  
  endoa<-as.vector(c(0,0,0,1,2,3)) ##identifier for each endogenous predictor, 1:gwesp; 2:degpop; 3:degcor
  
  
  badhes1<-try(solve(m1$hessian))
  warntest<-if(is(badhes1, "try-error")){1} else{0}
  badhes2<-try(solve(m2$hessian))
  warntest2<-if(is(badhes2, "try-error")){1} else{0}
  
  warntest<-as.vector(rep(warntest, each=6))
  warntest2<-as.vector(rep(warntest2, each=6))
  
  piece<-cbind(loopid, warntest, warntest2,
               endoa,vifs1,vifs2,vifschange, dens,
               se1,se2,sechange,varchange)
  
  
  
  mega<-if(is.data.frame(mega)){as.data.frame(rbind(mega,piece))} else(as.data.frame(piece))
  save(mega, file='ExampleVIF.SEchange.csv')
}

#################################################