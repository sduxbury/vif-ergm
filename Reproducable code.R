###Reproducable code


###For multilevel regression
  ##Load the total dataset
load("ranvifssmall.RData")
  ##data for only random networks
load("randomnetworkvifs.RData")
  ##data for only scale-free networks
load("sfnetworkvifs.RData")
  ##data for small-world networks
load("swnetworkvifs.RData")


##Multilevel regression
library(lme4)
library(piecewiseSEM)
library(lmtest)

##baseline models
nullmodel<-lmer(log(vifs)^0.1~1+(1|modid), data=ranvifset, na.action=na.exclude)
summary(nullmodel)                

m2lm<-lm(log(vifs)^0.1~dens+VMR+topog+gwesp+degcor+degpop+degpopgwespr+degpopdegcor+
  gwespdegcor+threevar, data=ranvifset, na.action=na.exclude)


###table 2
m2a<-lmer(log(vifs)^0.1~dens+VMR+topog+gwesp+degcor+degpop+degpopgwespr+degpopdegcor+
            gwespdegcor+threevar+(1|modid), data=ranvifset, na.action=na.exclude)
summ2a<-summary(m2a) ##model 1
rm2a<-rsquared(m2a)

##HLM is the correct modelling choice, even though variance component=0
lrtest(m2a,m2lm)

m3<-lmer(log(vifs)^0.1~VMR+dens+topog+gwesp+degcor+degpop+degpopgwespr+degpopdegcor+
           gwespdegcor+threevar+gwesp:dens+degcor:dens+degpop:dens+degpopgwespr:dens+degpopdegcor:dens+
           gwespdegcor:dens+threevar:dens+(1|modid), data=ranvifset, na.action=na.exclude)
summ3<-summary(m3)##model 2
rm3<-rsquared(m3)



####model fit
library(lmerTest)

#p values for random effects
randm2a<-rand(m2a)
randm3<-rand(m3)








######Table 3
##randomnetwork
mrandom<-lmer(log(vifs)^0.1~VMR+dens+gwesp+degcor+degpop+degpopgwespr+degpopdegcor+
                gwespdegcor+threevar+gwesp:dens+degcor:dens+degpop:dens+degpopgwespr:dens+degpopdegcor:dens+
                gwespdegcor:dens+threevar:dens+(1|modid), data=randomnetworks, na.action=na.exclude)
summrandom<-summary(mrandom)
rmrandom<-rsquared(mrandom)
randrandom<-rand(mrandom)

###Scale free
mscalefree<-lmer(log(vifs)^0.1~VMR+dens+gwesp+degcor+degpop+degpopgwespr+degpopdegcor+
                   gwespdegcor+threevar+gwesp:dens+degcor:dens+degpop:dens+degpopgwespr:dens+degpopdegcor:dens+
                   gwespdegcor:dens+threevar:dens+(1|modid), data=sfnetworks, na.action=na.exclude)
summsf<-summary(mscalefree)
rmscalefree<-rsquared(mscalefree)
randscalefree<-rand(mscalefree)

###small world
msmallworld<-lmer(log(vifs)^0.1~VMR+dens+gwesp+degcor+degpop+degpopgwespr+degpopdegcor+
                    gwespdegcor+threevar+gwesp:dens+degcor:dens+degpop:dens+degpopgwespr:dens+degpopdegcor:dens+
                    gwespdegcor:dens+threevar:dens+(1|modid), data=swnetworks, na.action=na.exclude)
summsw<-summary(msmallworld)
rmsmallworld<-rsquared(msmallworld)
randsmallworld<-rand(msmallworld)





#####"Interpreting ERGM VIFs" --variance fluctuations results
load("main files for SE change.RData")

###total mega is VIFs for each variable
###loopsmega is data for each model


###Table 4
##focal variable
cor(log(totalmega$vifsmean)^0.1, log(totalmega$vartot)) ###0.712
cor(log(totalmega$vifsmean)^0.1, log(totalmega$setot)) ##0.652

#model maximum variance fluctuation
cor(log(loopsmega$twomaxvifsperloop)^0.1, log(loopsmega$maxvartot))
cor(log(loopsmega$twomaxvifsperloop)^0.1, log(loopsmega$maxsetot))

#mean model variance fluctuation
cor(log(loopsmega$twomaxvifsperloop)^0.1, log(loopsmega$meanvartot))
cor(log(loopsmega$twomaxvifsperloop)^0.1, log(loopsmega$meansetot))



##Table 5
#############focal var
##variance
mean(totalmega$vartot[which(totalmega$vifsmean>=20 & totalmega$vifsmean<=150)]) ###0.65
mean(totalmega$vartot[which(totalmega$vifsmean>=20)]) ##0.55
mean(totalmega$vartot[which(totalmega$vifsmean>=150)]) ##0.196


##standard errors
mean(totalmega$setot[which(totalmega$vifsmean<=20)]) ##0.000

mean(totalmega$setot[which(totalmega$vifsmean>=20 & totalmega$vifsmean<=150)]) ###0.23
mean(totalmega$setot[which(totalmega$vifsmean>=20)]) ##0.20
mean(totalmega$setot[which(totalmega$vifsmean>=150)]) ##0.10




###models
###max
mean(loopsmega$maxvartot[which(loopsmega$twomaxvifsperloop>=20 &loopsmega$twomaxvifsperloop<=150)])###0.63
mean(loopsmega$maxvartot[which(loopsmega$twomaxvifsperloop<=20)])###0.57
mean(loopsmega$maxvartot[which(loopsmega$twomaxvifsperloop>=150)])###0.48


mean(loopsmega$maxsetot[which(loopsmega$twomaxvifsperloop>=20 &loopsmega$twomaxvifsperloop<=150)])###0.23
mean(loopsmega$maxsetot[which(loopsmega$twomaxvifsperloop>=20)])###0.22
mean(loopsmega$maxsetot[which(loopsmega$twomaxvifsperloop>=150)])###0.20



####mean

mean(loopsmega$meanvartot[which(loopsmega$twomaxvifsperloop>=20 &loopsmega$twomaxvifsperloop<=150)])###0.25
mean(loopsmega$meanvartot[which(loopsmega$twomaxvifsperloop>=20)])###0.22
mean(loopsmega$meanvartot[which(loopsmega$twomaxvifsperloop>=150)])###0.14


mean(loopsmega$meansetot[which(loopsmega$twomaxvifsperloop>=20 &loopsmega$twomaxvifsperloop<=150)])###0.0.09
mean(loopsmega$meansetot[which(loopsmega$twomaxvifsperloop>=20)])###0.08
mean(loopsmega$meansetot[which(loopsmega$twomaxvifsperloop>=150)])###0.06








###empirical example: Faux mesa high
##Table 6

library(statnet)
data(faux.mesa.high)

##model 1
a1<-ergm(faux.mesa.high~edges+nodematch("Race")+nodematch("Sex")+nodematch("Grade"), 
         control=control.ergm(force.main=TRUE, main.method="MCMLE"))

m2<-simulate(a1,statsonly=TRUE,nsim=10000)


###calculate the correlation matrix
corr5<-cor(m2)


viff1<-function(x){
  gvec<-x[3:4,2] 
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[3:4,3:4] 
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)}


##VIF for second variable
viff2<-function(x){
  gvec<-x[c(2,4),3]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,4),c(2,4)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
}


viff3<-function(x){
  gvec<-x[c(2,3),4]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3),c(2,3)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
  
}

###print VIFS
va1<-viff1(corr5)
va2<-viff2(corr5)
va3<-viff3(corr5)




##model 2

a2<-ergm(faux.mesa.high~edges+nodematch("Race")+nodematch("Sex")+nodematch("Grade")+gwdegree(0.7, fixed=TRUE),
         control=control.ergm(force.main=TRUE, main.method="MCMLE"))


m2<-simulate(a2,statsonly=TRUE,nsim=10000)


###calculate the correlation matrix
corr5<-cor(m2)


viff1<-function(x){
  gvec<-x[3:5,2] 
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[3:5,3:5] 
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)}


##VIF for second variable
viff2<-function(x){
  gvec<-x[c(2,4:5),3]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,4:5),c(2,4:5)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
}


viff3<-function(x){
  gvec<-x[c(2,3,5),4]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,5),c(2,3,5)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
  
}

viff4<-function(x){
  gvec<-x[c(2,3,4),5]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,4),c(2,3,4)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
  
}



vb1<-viff1(corr5)
vb2<-viff2(corr5)
vb3<-viff3(corr5)
vb4<-viff4(corr5)








###model 3

a4<-ergm(faux.mesa.high~edges+nodematch("Race")+nodematch("Sex")+nodematch("Grade")+gwesp(decay=0.7, fixed=TRUE)+
           gwdegree(0.7, fixed=TRUE),
         control=control.ergm(force.main=TRUE))




m2<-simulate(a4,statsonly=TRUE,nsim=10000)


###calculate the correlation matrix
corr5<-cor(m2)




viff1<-function(x){
  gvec<-x[3:6,2] ###This vector should be the column for the variable you are examining and all rows excluding the edges term and variable you are examining
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[3:6,3:6] ###square matrix that excludes the edges term and the variable you are interested in
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)}


##VIF for second variable
viff2<-function(x){
  gvec<-x[c(2,4:6),3]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,4:6),c(2,4:6)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
}


viff3<-function(x){
  gvec<-x[c(2,3,5:6),4]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,5:6),c(2,3,5:6)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
  
}


viff4<-function(x){
  gvec<-x[c(2:4,6),5]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2,3,4,6),c(2,3,4,6)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
  
}

viff5<-function(x){
  gvec<-x[c(2:5),6]
  gvec<-as.vector(gvec)
  tgvec<-t(gvec)
  xcorr<-x[c(2:5),c(2:5)]
  xcor<-solve(xcorr)
  Rsq<-tgvec%*%xcor%*%gvec
  vifs<-1/(1-Rsq)
}






vc1<-viff1(corr5)
vc2<-viff2(corr5)
vc3<-viff3(corr5)
vc4<-viff4(corr5)
vc5<-viff5(corr5)











###supplementary analysis: VIFs in large networks
##load large networks
  ##random networkds
load("vifslargenetworksupplement.RData")
  ##scale-free networks
load("sfvifslargenetworksupplement.RData")
  ##small world networks
load("swvifslargenetworksupplement.RData")

###load small networks of same conditions
  #random
load("vifs16.csv")
  ##small-world
load("swvifs16a.csv")
  ##scale-free
load("sfvifs16a.csv")

###note "vifslarge" refers to large network, "vifs" is small network

##scalefree
mean(sfvifslarge)
sd(sfvifslarge)

mean(sfvifs)
sd(sfvifs)

t.test(sfvifslarge,sfvifs)

###random
mean(vifslarge$vifs)
sd(vifslarge$vifs)

mean(vifs)
sd(vifs)

t.test(vifslarge$vifs,vifs)

##small-world
mean(swvifslarge$vifs)
sd(swvifslarge$vifs)

mean(swvifs)
sd(swvifs)

t.test(swvifslarge,swvifs)

###total
totallarge<-cbind(as.vector(sfvifslarge),swvifslarge$vifs,vifslarge$vifs)
totalsmall<-rbind(sfvifs,swvifs,vifs)

mean(totallarge)
sd(totallarge)

mean(totalsmall)
sd(totalsmall)

t.test(totallarge,totalsmall)

