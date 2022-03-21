rm(list=ls())
######################import package######################
# install.packages("deSolve","NLRoot","mvtnorm","truncnorm","expm","fda","MASS","plyr","dplyr")
library(deSolve)
library(NLRoot)
library(mvtnorm)
library(truncnorm)
library(expm)
library(fda)
library(MASS)
library(plyr)
library(dplyr)
##################################set seed####################################################
set.seed(123)

####################eval.basis#########################
source("my_bsplineS.R")
source("my_getbasismatrix.R")
source("my_eval.basis.R")
##################SNTO#####################
source("SNTO_knots.R")
################## generate design ##########################
#input parameter
N1 = 120;N2 = 10
N=N0;min_knots=5
sr = 0.5
source("Design.R")
# load("realdataU2_20_N138.Rdata")
# load("realdataU1_200_N138.Rdata")
###########################Initialization-chapter5############################
source("ode_equation2_obj2.R")

RealData = read.csv("BlowfilesData.csv",header = F)
N0 = nrow(RealData)
RealData_names=c("T","yt")
names(RealData) = RealData_names
Time=RealData[,1]
yt=RealData[,2]
yt[121]=1
lyt=log(yt)

v0=abs(round(rnorm(1, 0, 5),2))
gamma0=abs(round(runif(1, 0, 50),2))
p0=abs(round(rnorm(1, 0, 5),2))

N=N0
sr = 0.5;delta = 10^(-8);lambda=5
N_list=c(20,30,40,50,60,70,80,90,100,110,120,136)

##run
R= SNTO_knots(lambda,N_list,N,Time,U1,U2,N1,N2,delta)

###########################output knots with min{CMINn}############################
CMINn=vector()
for (i in 1:length(N_list)){
  CMINn[i]=R[[i]][N_list[i]+1]
}
n=N_list[which.min(CMINn)]
print(n)
Kbest=R[[which.min(CMINn)]]
print(Kbest[-(n+1)])
