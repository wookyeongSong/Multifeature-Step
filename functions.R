##########################################################################################################
##### Multi-feature Clustering of Step Data using Multivariate Functional Principal Component Analysis
##### Built-in Functions. Automatically loaded by "simulation_code.R" 
##### Date : 2020/11/18
##### Author : Wookyeong Song
###########################################################################################################
###########################################################################################################
library(fdapace)
library(fda)
library(cluster)
library(factoextra)
library(fields)
library(lubridate)
library(doBy)
library(plyr)
library(combinat)
library(NbClust)

## load step data
load("C:\\Users\\user\\Desktop\\WK\\class\\intern\\hsoh(step)\\step_79.RData")

step = full_step2[-1,]

## Standardized Cumulative Sum Function
## input : step - p x N matrix (Ex. p = 1440, N= 21394)
convert_amount <- function(step) {
  
  data_amount = apply(step,2,cumsum)
  
  std_data_amount = data_amount/mean(data_amount[nrow(step),])
  
  std_data_amount
}

## Standardized Ordered Quantile Slope Function
## input : step - p x N matrix (Ex. p = 1440, N= 21394)
##         K - number of quantile Q1 
convert_intensity <-function(step,K) {
  
  intmatrix <- matrix(0,nrow=K+1,ncol=ncol(step))
  for(i in 1:ncol(step)) {
    total_act=sum(step[,i])
    s_t=cumsum(step[,i])
    for(j in 2:K) {
      intmatrix[j,i]<-min(which(s_t>(j-1)* total_act /K  ) )
    }
    intmatrix[(K+1),i]<- nrow(step)
  }
  
  data_intensity = matrix(NA,nrow=nrow(step), ncol=ncol(step))
  num = as.integer(nrow(step)/K)
  for(i in 1:ncol(data_intensity)) {
    for(j in 1:K) {
      if( intmatrix[j+1,i]!=  intmatrix[j,i] ){
        data_intensity[c( (intmatrix[j ,i]+1 ):  intmatrix[j+1 ,i])  ,i]<- sum(step[,i])/(K*(intmatrix[(j+1),i]-intmatrix[j,i])  )
      }
    }
  }
  
  data_intensity = apply(data_intensity,2,sort)
  
  std_data_intensity = data_intensity/mean(data_intensity[nrow(step),])
  
  std_data_intensity
}

## Standardized Mean Score Function
## input : step - p x N matrix (Ex. p = 1440, N= 21394)
##         K - number of quantile Q2 
convert_pattern <- function(step,K) {
  
  pj=c(0:K)/K
  
  patmatrix <- matrix(nrow=(K+1),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    patmatrix[1,i]<- 0
  }
  for(i in 1:ncol(step)) {
    for(j in 2:K) {
      s_t=cumsum(sort(step[,i]))
      patmatrix[j,i]<-which(s_t == s_t[ s_t>= ((j-1)*sum(step[,i])/K)  ][1])  [1]  ##t value that have sc(T)=j
    }
  }
  for(i in 1:ncol(step)) {
    patmatrix[(K+1),i]<- nrow(step)
  }
  
  ### data_pattern_score
  data_pattern_score=matrix(nrow=nrow(step),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    for(j in 1:K) {
      data_pattern_score[which(patmatrix[(j+1),i]>=rank(step[,i])  &  rank(step[,i])>patmatrix[j,i]),i]<-(j-1)
    }
  }
  
  ### data_pattern
  data_pattern=matrix(nrow=nrow(step),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    for(j in 1:K) {
      data_pattern[((nrow(step)/K)*j-(nrow(step)/K-1)):((nrow(step)/K)*j),i]<-mean(data_pattern_score[((nrow(step)/K)*j-(nrow(step)/K-1)):((nrow(step)/K)*j),i])
    }
  }
  
  sort_data_pattern=apply(data_pattern,2,sort)
  
  std_data_pattern= data_pattern/mean(sort_data_pattern[nrow(step),])
  
  std_data_pattern
}

## Construct Multi-Feature Step Data
## input : step - p x N matrix (Ex. p = 1440, N= 21394)
##         Q1 - number of quantiles for the ordered slope function
##         Q2 - number of quantiles for the mean score function
MultiFeat_data <- function(step,Q1,Q2) {
  
  data_amount = convert_amount(step)
  data_intensity = convert_intensity(step,Q1)
  data_pattern = convert_pattern(step,Q2)
  
  MultiFeat_data_array = array(data=NA,dim=c(nrow(step),ncol(step),3))
  
  MultiFeat_data_array[,,1]<-data_amount
  MultiFeat_data_array[,,2]<-data_intensity
  MultiFeat_data_array[,,3]<-data_pattern
  
  MultiFeat_data_array
}

## Obtain MFPCA score
## input : MultiFeat_data_array - Multi-Feature step data from MultiFeat_data function
##         threshold - percentages of total variance from MFPC
MultiFeat_MFPCA = function(MultiFeat_data_array,threshold) {
  
  xdim = dim(MultiFeat_data_array[,,1])[1]
  ## Cubic B-spline basis
  fdabasis = create.bspline.basis(rangeval=c(0,xdim),norder=4,breaks = seq(0,xdim,by=6))
  fdatime = seq(1, xdim, by=1)
  fdafd = smooth.basis(fdatime, MultiFeat_data_array, fdabasis)$fd
  
  stepPca = pca.fd(fdafd,nharm=10)
  
  cum_varprop <- vector()
  for(q in 1:10){
    cum_varprop[q]=sum(stepPca$varprop[1:q])
  }
  
  cum_varprop
  num_PCscore=which(cum_varprop>threshold)[1]
  
  fun_PCscore<-stepPca$scores[,1:num_PCscore,1]+stepPca$scores[,1:num_PCscore ,2]+stepPca$scores[, 1:num_PCscore,3]
  
  fun_PCscore
}

## Finding best match cluster when obtaining CCR
## input : x - true label
##         y - our clustering result
best_match_factor=function(x,y){
  x=factor(x)
  A=NA
  for(k in 1:length(permn(length(table(x))))){
    new_b= as.factor(y)
    newy=mapvalues(new_b, from = levels(new_b), to =  as.character( permn(length(table(x)))[k][[1]] ) )
    A=c(A,  length(which(newy==x) ) )
  }
  A=A[-1]
  newy=mapvalues(new_b, from = levels(new_b), to =  as.character( permn(length(table(x)))[which.max(A)][[1]] ) )
  return(as.numeric(as.character(newy)) )
}


show_result = function(ccr_km,)


## Functions below are used to draw figures
convert_nonorder_intensity <-function(step,K) {

  intmatrix <- matrix(0,nrow=K+1,ncol=ncol(step))
  
  for(i in 1:ncol(step)) {
    total_act=sum(step[,i])
    s_t=cumsum(step[,i])
    for(j in 2:K) {
      intmatrix[j,i]<-min(which(s_t>(j-1)* total_act /K  ) )
    }
    intmatrix[(K+1),i]<- nrow(step)
  }
  
  data_intensity = matrix(NA,nrow=nrow(step), ncol=ncol(step))
  num = as.integer(nrow(step)/K)
  for(i in 1:ncol(data_intensity)) {
    for(j in 1:K) {
      if( intmatrix[j+1,i]!=  intmatrix[j,i] ){
        data_intensity[c( (intmatrix[j ,i]+1 ):  intmatrix[j+1 ,i])  ,i]<- sum(step[,i])/(K*(intmatrix[(j+1),i]-intmatrix[j,i])  )
      }
    }
  }
  
  data_intensity
}

convert_nonstandard_intensity <-function(step,K) {
  K=K
  
  intmatrix <- matrix(0,nrow=K+1,ncol=ncol(step))
  
  for(i in 1:ncol(step)) {
    total_act=sum(step[,i])
    s_t=cumsum(step[,i])
    for(j in 2:K) {
      intmatrix[j,i]<-min(which(s_t>(j-1)* total_act /K  ) )
    }
    intmatrix[(K+1),i]<- nrow(step)
  }
  
  
  data_intensity = matrix(NA,nrow=nrow(step), ncol=ncol(step))
  num = as.integer(nrow(step)/K)
  for(i in 1:ncol(data_intensity)) {
    for(j in 1:K) {
      if( intmatrix[j+1,i]!=  intmatrix[j,i] ){
        data_intensity[c( (intmatrix[j ,i]+1 ):  intmatrix[j+1 ,i])  ,i]<- sum(step[,i])/(K*(intmatrix[(j+1),i]-intmatrix[j,i])  )
      }
    }
  }
  
  sort_data_intensity = apply(data_intensity,2,sort)
  
  sort_data_intensity
}


convert_nonstandard_pattern <- function(step,K) {
  
  patmatrix <- matrix(nrow=(K+1),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    patmatrix[1,i]<- 0
  }
  for(i in 1:ncol(step)) {
    for(j in 2:K) {
      s_t=cumsum(sort(step[,i]))
      patmatrix[j,i]<-which(s_t == s_t[ s_t>= ((j-1)*sum(step[,i])/K)  ][1])  [1]  ##t value that have sc(T)=j
    }
  }
  for(i in 1:ncol(step)) {
    patmatrix[(K+1),i]<- nrow(step)
  }
  
  ### data_pattern_score
  data_pattern_score=matrix(nrow=nrow(step),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    for(j in 1:K) {
      data_pattern_score[which(patmatrix[(j+1),i]>=rank(step[,i])  &  rank(step[,i])>patmatrix[j,i]),i]<-(j-1)
    }
  }
  
  ### data_pattern
  data_pattern=matrix(nrow=nrow(step),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    for(j in 1:K) {
      data_pattern[((nrow(step)/K)*j-(nrow(step)/K-1)):((nrow(step)/K)*j),i]<-mean(data_pattern_score[((nrow(step)/K)*j-(nrow(step)/K-1)):((nrow(step)/K)*j),i])
    }
  }

  data_pattern
}

