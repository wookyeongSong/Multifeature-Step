##########################################################################################################
##### Multi-feature Clustering of Step Data using Multivariate Functional Principal Component Analysis
##### Simulation Codes
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
library(funFEM)
library(funHDDC)
library(mclust)
library(combinat)

source('C:\\Users\\user\\Multifeature-Step\\functions.R')


###########################################################################################################
###########################################################################################################
###########################################################################################################
## Simulation 1
## Step-like Simulation Data - Amount
## input : k - seed number (1~100)
###########################################################################################################
###########################################################################################################
###########################################################################################################
step_like_simul_amount <-function(k) {
  set.seed(k)
  x_low <- as.integer(rnorm(100,150,15))
  x_mid <- as.integer(rnorm(100,250,15))
  x_high <- as.integer(rnorm(100,350,15))
  
  sample_low_afternoon <- matrix(0,nrow=1440,ncol = 100)
  inten_prop_afternoon = c(0,0.02,0.04,0.08,0.15,0.25,0.25,0.1,0.06,0.03,0.02,0)
  tmp1 = inten_prop_afternoon[c(1:2,11:12)]
  tmp2 = inten_prop_afternoon[5:8]
  tmp3 = inten_prop_afternoon[c(3:4,9:10)]
  
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_prop_afternoon <- c(tmp1[1:2],tmp3[1:2],tmp2,tmp3[3:4],tmp1[3:4])
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_low[i]*inten_prop_afternoon[j]),replace = FALSE)
      sample_low_afternoon[walk_location,i] <-as.integer(rexp(as.integer(x_low[i]*inten_prop_afternoon[j]),1/32.5))
    }
  }
  
  sample_mid_afternoon <- matrix(0,nrow=1440,ncol = 100)
  inten_prop_afternoon = c(0,0.02,0.04,0.08,0.15,0.25,0.25,0.1,0.06,0.03,0.02,0)
  tmp1 = inten_prop_afternoon[c(1:2,11:12)]
  tmp2 = inten_prop_afternoon[5:8]
  tmp3 = inten_prop_afternoon[c(3:4,9:10)]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_prop_afternoon <- c(tmp1[1:2],tmp3[1:2],tmp2,tmp3[3:4],tmp1[3:4])
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_mid[i]*inten_prop_afternoon[j]),replace = FALSE)
      sample_mid_afternoon[walk_location,i] <-as.integer(rexp(as.integer(x_mid[i]*inten_prop_afternoon[j]),1/32.5))
    }
  }
  
  sample_high_afternoon <- matrix(0,nrow=1440,ncol = 100)
  inten_prop_afternoon = c(0,0.02,0.04,0.08,0.15,0.25,0.25,0.1,0.06,0.03,0.02,0)
  tmp1 = inten_prop_afternoon[c(1:2,11:12)]
  tmp2 = inten_prop_afternoon[5:8]
  tmp3 = inten_prop_afternoon[c(3:4,9:10)]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_prop_afternoon <- c(tmp1[1:2],tmp3[1:2],tmp2,tmp3[3:4],tmp1[3:4])
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_high[i]*inten_prop_afternoon[j]),replace = FALSE)
      sample_high_afternoon[walk_location,i] <-as.integer(rexp(as.integer(x_high[i]*inten_prop_afternoon[j]),1/32.5))
    }
  }
  
  step_like_amount = cbind(sample_low_afternoon,sample_mid_afternoon,sample_high_afternoon)
  
  step_like_amount
}


true = c(rep(1,100),rep(2,100),rep(3,100))

ccr <- function(result, true) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == result)/length(result))
}

ccr_kmeans_amount=vector()
rand_kmeans_amount=vector()

ccr_pam_amount=vector()
rand_pam_amount=vector()

ccr_funfem_amount=vector()
rand_funfem_amount=vector()

ccr_funhddc_amount=vector()
rand_funhddc_amount=vector()

## Perform simulation
for (k in 1:100) {
  step_like_amount = step_like_simul_amount(k)
  
  MultiFeat_data_array = MultiFeat_data(step_like_amount,8,3)
  
  funPCscore = MultiFeat_MFPCA(MultiFeat_data_array,0.9)
  
  ##Kmeans
  num_cl=3
  MCluster <- kmeans(funPCscore, centers= num_cl, iter.max = 100000)
  
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_kmeans_amount[k] <- ccr(true,cluster_result)
  rand_kmeans_amount[k] <- adjustedRandIndex(true,cluster_result)
  
  ##PAM
  num_cl=3
  MCluster <- pam(funPCscore, num_cl, metric="manhattan")
  
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_pam_amount[k] <- ccr(true,cluster_result)
  rand_pam_amount[k] <- adjustedRandIndex(true,cluster_result)
  
  ## FunFEM
  fdabasis = create.bspline.basis(rangeval=c(0,1440),norder=4,nbasis=21)
  fdatime = seq(1, 1440, by=1)
  sim_fdafd = smooth.basis(fdatime, step_like_amount, fdabasis)$fd
  
  res <- funFEM(sim_fdafd,model='all',K=3)
  cluster_result = best_match_factor(true,res$cls)
  
  ccr_funfem_amount[k] <- ccr(true,cluster_result)
  rand_funfem_amount[k] <- adjustedRandIndex(true,cluster_result)
  
  ##FunHDDC
  fdabasis = create.fourier.basis(rangeval=c(0,1440),period = 1440,nbasis=21)
  fdatime = seq(1, 1440, by=1)
  sim_fdafd = smooth.basis(fdatime, step_like_amount, fdabasis)$fd
  
  res <- funHDDC(sim_fdafd,model=c( 'AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),K=3)
  cluster_result = best_match_factor(true,res$class)
  
  ccr_funhddc_amount[k] <- ccr(true,cluster_result)
  rand_funhddc_amount[k] <- adjustedRandIndex(true,cluster_result)
}

## result
step_sim_amount = cbind(ccr_kmeans_amount,rand_kmeans_amount,ccr_pam_amount,rand_pam_amount,ccr_funfem_amount,rand_funfem_amount,ccr_funhddc_amount,rand_funhddc_amount)
apply(step_sim_amount,2,mean)
apply(step_sim_amount,2,sd)

###########################################################################################################
###########################################################################################################
###########################################################################################################
## Simulation 2
## Step-like Simulation Data - Intensity
## input : k - seed number (1~100)
###########################################################################################################
###########################################################################################################
###########################################################################################################
step_like_simul_intensity <-function(k) {
  set.seed(k)
  x_intensity <- as.integer(rnorm(300,150,5))
  
  sample_low_evening <- matrix(0,nrow=1440,ncol = 100)
  inten_low_evening = sort(c(0.05,0.05,0.05,0.05,0.075,0.075,0.075,0.075,0.125,0.125,0.125,0.125))
  tmp1 = inten_low_evening[1:4]
  tmp2 = inten_low_evening[5:8]
  tmp3 = inten_low_evening[9:12]
  
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_low_evening <- c(tmp1,tmp2,tmp3)
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_intensity[i]*inten_low_evening[j]),replace = FALSE)
      sample_low_evening[walk_location,i] <-as.integer(rexp(as.integer(x_intensity[i]*inten_low_evening[j]),1/32.5))
    }
  }
  
  sample_mid_evening <- matrix(0,nrow=1440,ncol = 100)
  inten_mid_evening = c(0,0.1,0.1,0.0,0.15,0.15,0.0,0.0,0.0,0.0,0.25,0.25)
  tmp1 = inten_mid_evening[1:4]
  tmp2 = inten_mid_evening[5:8]
  tmp3 = inten_mid_evening[9:12]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_mid_evening <- c(tmp1,tmp2,tmp3)
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_intensity[100+i]*inten_mid_evening[j]),replace = FALSE)
      sample_mid_evening[walk_location,i] <-as.integer(rexp(as.integer(x_intensity[100+i]*inten_mid_evening[j]),1/32.5))
    }
  }
  
  sample_high_evening <- matrix(0,nrow=1440,ncol = 100)
  inten_high_evening = c(0,0.0,0.2,0.0,0.0,0.0,0.3,0.0,0.0,0.0,0.0,0.5)
  tmp1 = inten_high_evening[1:4]
  tmp2 = inten_high_evening[5:8]
  tmp3 = inten_high_evening[9:12]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_high_evening <- c(tmp1,tmp2,tmp3)
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_intensity[200+i]*inten_high_evening[j]),replace = FALSE)
      sample_high_evening[walk_location,i] <-as.integer(rexp(as.integer(x_intensity[200+i]*inten_high_evening[j]),1/32.5))
    }
  }
  
  step_like_intensity = cbind(sample_low_evening,sample_mid_evening,sample_high_evening)
  
  step_like_intensity
}

true = c(rep(1,100),rep(2,100),rep(3,100))

ccr <- function(result, true) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == result)/length(result))
}

ccr_kmeans_intensity=vector()
rand_kmeans_intensity=vector()

ccr_pam_intensity=vector()
rand_pam_intensity=vector()

ccr_funfem_intensity=vector()
rand_funfem_intensity=vector()

ccr_funhddc_intensity=vector()
rand_funhddc_intensity=vector()


for (k in 1:100) {
  step_like_intensity = step_like_simul_intensity(k)
  
  MultiFeat_data_array = MultiFeat_data(step_like_intensity,8,3)
  
  funPCscore = MultiFeat_MFPCA(MultiFeat_data_array,0.9)
  
  ##Kmeans
  num_cl=3
  MCluster <- kmeans(funPCscore, centers= num_cl, iter.max = 100000)
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_kmeans_intensity[k] <- ccr(true,cluster_result)
  rand_kmeans_intensity[k] <- adjustedRandIndex(true,cluster_result)
  
  ##PAM
  num_cl=3
  MCluster <- pam(funPCscore, num_cl, metric="manhattan")
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_pam_intensity[k] <- ccr(true,cluster_result)
  rand_pam_intensity[k] <- adjustedRandIndex(true,cluster_result)
  
  ## FunFEM
  fdabasis = create.bspline.basis(rangeval=c(0,1440),norder=4,nbasis=21)
  fdatime = seq(1, 1440, by=1)
  sim_fdafd = smooth.basis(fdatime, step_like_intensity, fdabasis)$fd
  
  res <- funFEM(sim_fdafd,model='all',K=3)
  cluster_result = best_match_factor(true,res$cls)

  ccr_funfem_intensity[k] <- ccr(true,cluster_result)
  rand_funfem_intensity[k] <- adjustedRandIndex(true,cluster_result)
  
  ##FunHDDC
  fdabasis = create.fourier.basis(rangeval=c(0,1440),period = 1440,nbasis=21)
  fdatime = seq(1, 1440, by=1)
  sim_fdafd = smooth.basis(fdatime, step_like_intensity, fdabasis)$fd
  
  res <- funHDDC(sim_fdafd,model=c( 'AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),K=3)
  cluster_result = best_match_factor(true,res$class)
  
  ccr_funhddc_intensity[k] <- ccr(true,cluster_result)
  rand_funhddc_intensity[k] <- adjustedRandIndex(true,cluster_result)
  
}

step_sim_intensity = cbind(ccr_kmeans_intensity,rand_kmeans_intensity,ccr_pam_intensity,rand_pam_intensity,ccr_funfem_intensity,rand_funfem_intensity,ccr_funhddc_intensity,rand_funhddc_intensity)
apply(step_sim_intensity,2,mean)
apply(step_sim_intensity,2,sd)


###########################################################################################################
###########################################################################################################
###########################################################################################################
## Simulation 3
## Step-like Simulation Data - Pattern
## input : k - seed number (1~100)
###########################################################################################################
###########################################################################################################
###########################################################################################################
step_like_simul_pattern <-function(k) {
  set.seed(k)
  x_pattern <- as.integer(rnorm(300,250,15))
  sample_mid_morning <- matrix(0,nrow=1440,ncol = 100)
  inten_prop_morning = c(0.12,0.13,0.1,0.1,0.1,0.1,0.07,0.08,0.05,0.05,0.05,0.05)
  tmp1 = inten_prop_morning[1:4]
  tmp2 = inten_prop_morning[5:8]
  tmp3 = inten_prop_morning[9:12]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_prop_morning <- c(tmp1,tmp2,tmp3)
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_pattern[i]*inten_prop_morning[j]),replace = FALSE)
      sample_mid_morning[walk_location,i] <-as.integer(rexp(as.integer(x_pattern[i]*inten_prop_morning[j]),1/20))
    }
  }
  
  sample_mid_afternoon <- matrix(0,nrow=1440,ncol = 100)
  inten_prop_afternoon = c(0.1,0.07,0.08,0.1,0.1,0.13,0.12,0.1,0.05,0.05,0.05,0.05)
  tmp1 = inten_prop_afternoon[1:4]
  tmp2 = inten_prop_afternoon[5:8]
  tmp3 = inten_prop_afternoon[9:12]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_prop_afternoon <- c(tmp1,tmp2,tmp3)
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_pattern[100+i]*inten_prop_afternoon[j]),replace = FALSE)
      sample_mid_afternoon[walk_location,i] <-as.integer(rexp(as.integer(x_pattern[100+i]*inten_prop_afternoon[j]),1/20))
    }
  }
  
  sample_mid_evening <- matrix(0,nrow=1440,ncol = 100)
  inten_prop_evening = sort(inten_prop_morning)
  tmp1 = inten_prop_evening[1:4]
  tmp2 = inten_prop_evening[5:8]
  tmp3 = inten_prop_evening[9:12]
  for (i in 1:100) {
    tmp1 <- tmp1[sample(4,4,replace=FALSE)]
    tmp2 <- tmp2[sample(4,4,replace=FALSE)]
    tmp3 <- tmp3[sample(4,4,replace=FALSE)]
    inten_prop_evening <- c(tmp1,tmp2,tmp3) 
    for (j in 1:12) {
      walk_location = sample((120*j-119):(120*j),as.integer(x_pattern[200+i]*inten_prop_evening[j]),replace = FALSE)
      sample_mid_evening[walk_location,i] <-as.integer(rexp(as.integer(x_pattern[200+i]*inten_prop_evening[j]),1/20))
    }
  }
  step_like_pattern = cbind(sample_mid_morning,sample_mid_afternoon,sample_mid_evening)
  
  step_like_pattern
}

true = c(rep(1,100),rep(2,100),rep(3,100))

ccr <- function(result, true) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == result)/length(result))
}

ccr_kmeans_pattern=vector()
rand_kmeans_pattern=vector()

ccr_pam_pattern=vector()
rand_pam_pattern=vector()

ccr_funfem_pattern=vector()
rand_funfem_pattern=vector()

ccr_funhddc_pattern=vector()
rand_funhddc_pattern=vector()

for (k in 1:100) {
  step_like_pattern = step_like_simul_pattern(k)
  
  MultiFeat_data_array = MultiFeat_data(step_like_pattern,8,3)
  
  funPCscore = MultiFeat_MFPCA(MultiFeat_data_array,0.9)
  
  ##Kmeans
  num_cl=3
  MCluster <- kmeans(funPCscore, centers= num_cl, iter.max = 100000)
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_kmeans_pattern[k] <- ccr(true,cluster_result)
  rand_kmeans_pattern[k] <- adjustedRandIndex(true,cluster_result)
  
  ##PAM
  num_cl=3
  MCluster <- pam(funPCscore, num_cl, metric="manhattan")
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_pam_pattern[k] <- ccr(true,cluster_result)
  rand_pam_pattern[k] <- adjustedRandIndex(true,cluster_result)
  
  ## FunFEM
  fdabasis = create.bspline.basis(rangeval=c(0,1440),norder=4,nbasis=21)
  fdatime = seq(1, 1440, by=1)
  sim_fdafd = smooth.basis(fdatime, step_like_pattern, fdabasis)$fd
  
  res <- funFEM(sim_fdafd,model='all',K=3)
  cluster_result = best_match_factor(true,res$cls)
  
  ccr_funfem_pattern[k] <- ccr(true,cluster_result)
  rand_funfem_pattern[k] <- adjustedRandIndex(true,cluster_result)
  
  ##FunHDDC
  fdabasis = create.fourier.basis(rangeval=c(0,1440),period = 1440,nbasis=21)
  fdatime = seq(1, 1440, by=1)
  sim_fdafd = smooth.basis(fdatime, step_like_pattern, fdabasis)$fd
  
  res <- funHDDC(sim_fdafd,model = c('AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),K=3)
  cluster_result = best_match_factor(true,res$class)
  
  ccr_funhddc_pattern[k] <- ccr(true,cluster_result)
  rand_funhddc_pattern[k] <- adjustedRandIndex(true,cluster_result)
  
}
step_sim_pattern = cbind(ccr_kmeans_pattern,rand_kmeans_pattern,ccr_pam_pattern,rand_pam_pattern,ccr_funfem_pattern,rand_funfem_pattern,ccr_funhddc_pattern,rand_funhddc_pattern)
apply(step_sim_pattern,2,mean)
apply(step_sim_pattern,2,sd)



###########################################################################################################
###########################################################################################################
###########################################################################################################
## Simulation 4
## Sinusoidal with differnet amount 
## input : k - seed number (1~100)
###########################################################################################################
###########################################################################################################
###########################################################################################################
sinusoidal_amount <- function(k) {
  set.seed(k)
  time = seq(0,1023/1024,length.out=1024)
  
  sinusoidal1 = matrix(0,nrow=1024,ncol=50)
  for (i in 1:1024) {
    sinusoidal1[i,]= 0.5*abs(rnorm(50,sin(5*time[i]),0.3))
  }
  
  sinusoidal2 = matrix(0,nrow=1024,ncol=50)
  for (i in 1:1024) {
    sinusoidal2[i,]= 1*abs(rnorm(50,sin(5*time[i]),0.3))
  }
  
  sinusoidal3 = matrix(0,nrow=1024,ncol=50)
  for (i in 1:1024) {
    sinusoidal3[i,]= 1.5*abs(rnorm(50,sin(5*time[i]),0.3))
  }
  
  sinusoidal4 = matrix(0,nrow=1024,ncol=50)
  for (i in 1:1024) {
    sinusoidal4[i,]= 2*abs(rnorm(50,sin(5*time[i]),0.3))
  }
  sinusoidal = cbind(sinusoidal1,sinusoidal2,sinusoidal3,sinusoidal4)
  
  for (j in 1:200) {
    sinusoidal[,j] = sinusoidal[,j]*rbinom(1024,1,0.2)
  }
  
  sinusoidal
}


true = c(rep(1,50),rep(2,50),rep(3,50),rep(4,50))

ccr <- function(result, true) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == result)/length(result))
}

ccr_kmeans_sinamount=vector()
rand_kmeans_sinamount=vector()

ccr_pam_sinamount=vector()
rand_pam_sinamount=vector()

ccr_funfem_sinamount=vector()
rand_funfem_sinamount=vector()

ccr_funhddc_sinamount=vector()
rand_funhddc_sinamount=vector()


for (k in 1:100) {
  sinusoidal = sinusoidal_amount(k)
  
  MultiFeat_data_array = MultiFeat_data(sinusoidal,8,4)
  
  funPCscore = MultiFeat_MFPCA(MultiFeat_data_array,0.9)
  
  ##Kmeans
  num_cl=4
  MCluster <- kmeans(funPCscore, centers= num_cl, iter.max = 100000)
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_kmeans_sinamount[k] <- ccr(true,cluster_result)
  rand_kmeans_sinamount[k] <- adjustedRandIndex(true,cluster_result)
  
  ##PAM
  num_cl=4
  MCluster <- pam(funPCscore, num_cl, metric="manhattan")
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_pam_sinamount[k] <- ccr(true,cluster_result)
  rand_pam_sinamount[k] <- adjustedRandIndex(true,cluster_result)
  
  ## FunFEM
  fdabasis = create.bspline.basis(rangeval=c(0,1024),norder=4,nbasis=21)
  fdatime = seq(1, 1024, by=1)
  sim_sinusoidal_fdafd = smooth.basis(fdatime, sinusoidal, fdabasis)$fd
  res = funFEM(sim_sinusoidal_fdafd,model='all',K=4)
  
  cluster_result = best_match_factor(true,res$cls)
  
  ccr_funfem_sinamount[k] <- ccr(true,cluster_result)
  rand_funfem_sinamount[k] <- adjustedRandIndex(true,cluster_result)
  
  
  ##FunHDDC
  fdabasis = create.fourier.basis(rangeval=c(0,1024),period = 1024)
  fdatime = seq(1, 1024, by=1)
  sim_sinusoidal_fdafd = smooth.basis(fdatime, sinusoidal, fdabasis)$fd
  
  res <- funHDDC(sim_sinusoidal_fdafd,model=c( 'AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),K=4)
  cluster_result = best_match_factor(true,res$class)
  
  ccr_funhddc_sinamount[k] <- ccr(true,cluster_result)
  rand_funhddc_sinamount[k] <- adjustedRandIndex(true,cluster_result)
  
}

sinusoidal_sim_amount = cbind(ccr_kmeans_sinamount,rand_kmeans_sinamount,ccr_pam_sinamount,rand_pam_sinamount,ccr_funfem_sinamount,rand_funfem_sinamount,ccr_funhddc_sinamount,rand_funhddc_sinamount)
apply(sinusoidal_sim_amount,2,mean)
apply(sinusoidal_sim_amount,2,sd)

###########################################################################################################
###########################################################################################################
###########################################################################################################
## Simulation 5
## Doppler Shift with differnet pattern 
## input : k - seed number (1~100)
###########################################################################################################
###########################################################################################################
###########################################################################################################
doppler_pattern <- function(k) {
  set.seed(k)
  time = seq(0,511/512,length.out=512)
  
  doppler1 = matrix(0,nrow=512,ncol=50)
  doppler1[1,] = rnorm(50,0.6,0.05)
  for (i in 2:512) {
    doppler1[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i])),0.05)
  }
  
  doppler2 = matrix(0,nrow=512,ncol=50)
  for (i in 1:512) {
    doppler2[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-1/3)),0.05)
  }
  
  doppler3 = matrix(0,nrow=512,ncol=50)
  for (i in 1:512) {
    doppler3[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-2/3)),0.05)
  }
  
  doppler4 = matrix(0,nrow=512,ncol=50)
  for (i in 1:512) {
    doppler4[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-1)),0.05)
  }
  
  doppler = cbind(doppler1,doppler2,doppler3,doppler4)
  
  for (j in 1:200) {
    doppler[,j] = doppler[,j]*rbinom(512,1,0.2)
  }
  
  doppler
}

true = c(rep(1,50),rep(2,50),rep(3,50),rep(4,50))

ccr <- function(result, true) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == result)/length(result))
}

ccr_kmeans_dopshif=vector()
rand_kmeans_dopshif=vector()

ccr_pam_dopshif=vector()
rand_pam_dopshif=vector()

ccr_funfem_dopshif=vector()
rand_funfem_dopshif=vector()

ccr_funhddc_dopshif=vector()
rand_funhddc_dopshif=vector()

for (k in 1:100) {
  doppler = doppler_pattern(k)
  
  MultiFeat_data_array = MultiFeat_data(doppler,8,4)
  
  funPCscore = MultiFeat_MFPCA(MultiFeat_data_array,0.9)
  
  ##Kmeans
  num_cl=4
  MCluster <- kmeans(funPCscore, centers= num_cl, iter.max = 100000)
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_kmeans_dopshif[k] <- ccr(true,cluster_result)
  rand_kmeans_dopshif[k] <- adjustedRandIndex(true,cluster_result)
  
  ##PAM
  num_cl=4
  MCluster <- pam(funPCscore, num_cl, metric="manhattan")
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_pam_dopshif[k] <- ccr(true,cluster_result)
  rand_pam_dopshif[k] <- adjustedRandIndex(true,cluster_result)
  
  ## FunFEM
  fdabasis = create.bspline.basis(rangeval=c(0,512),norder=4,nbasis=21)
  fdatime = seq(1, 512, by=1)
  sim_doppler_fdafd = smooth.basis(fdatime, doppler, fdabasis)$fd
  res = funFEM(sim_doppler_fdafd,model='all',K=4)
  
  cluster_result = best_match_factor(true,res$cls)
  
  ccr_funfem_dopshif[k] <- ccr(true,cluster_result)
  rand_funfem_dopshif[k] <- adjustedRandIndex(true,cluster_result)
  
  ##FunHDDC
  fdabasis = create.fourier.basis(rangeval=c(0,512),period = 512)
  fdatime = seq(1, 512, by=1)
  sim_doppler_fdafd = smooth.basis(fdatime, doppler, fdabasis)$fd
  
  res <- funHDDC(sim_doppler_fdafd,model=c('AkjBkQkDk','AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),K=4)
  
  cluster_result = best_match_factor(true,res$class)
  
  ccr_funhddc_dopshif[k] <- ccr(true,cluster_result)
  rand_funhddc_dopshif[k] <- adjustedRandIndex(true,cluster_result)
  
}

doppler_sim_pattern = cbind(ccr_kmeans_dopshif,rand_kmeans_dopshif,ccr_pam_dopshif,rand_pam_dopshif,ccr_funfem_dopshif,rand_funfem_dopshif,ccr_funhddc_dopshif,rand_funhddc_dopshif)
apply(doppler_sim_pattern,2,mean)
apply(doppler_sim_pattern,2,sd)




###########################################################################################################
###########################################################################################################
###########################################################################################################
## Simulation 6
## Doppler shift - Intensity
## input : k - seed number (1~100)
###########################################################################################################
###########################################################################################################
###########################################################################################################
doppler_intensity <-function(k) {
  set.seed(k)
  time = seq(0,511/512,length.out=512)
  
  doppler = matrix(0,nrow=512,ncol=200)
  doppler1 = matrix(0,nrow=512,ncol=50)
  for (i in 1:512) {
    doppler1[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-0.1)),0.05)
  }
  
  for (j in 1:50) {
    doppler[,j] = doppler1[,j]*rbinom(512,1,0.2)
  }
  
  #doppler2 = matrix(0,nrow=512,ncol=50)
  #for (i in 1:512) {
  #  doppler2[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-0.1)),0.5)
  #}
  
  for (j in 1:50) {
    tmp1 = sample(c(1,2,3,4),2)
    tmp2 = sample(c(5,6,7,8),2)
    #tmp3 = sample(c(9,10,11,12),2)
    #tmp4 = sample(c(13,14,15,16),2)
    doppler[(64*tmp1[1]-63):(64*tmp1[1]),(50+j)] = doppler1[(64*tmp1[1]-63):(64*tmp1[1]),j]*rbinom(64,1,0.4)
    doppler[(64*tmp1[2]-63):(64*tmp1[2]),(50+j)] = doppler1[(64*tmp1[2]-63):(64*tmp1[2]),j]*rbinom(64,1,0.4)
    doppler[(64*tmp2[1]-63):(64*tmp2[1]),(50+j)] = doppler1[(64*tmp2[1]-63):(64*tmp2[1]),j]*rbinom(64,1,0.4)
    doppler[(64*tmp2[2]-63):(64*tmp2[2]),(50+j)] = doppler1[(64*tmp2[2]-63):(64*tmp2[2]),j]*rbinom(64,1,0.4)
    #doppler[(32*tmp3[1]-31):(32*tmp3[1]),(50+j)] = doppler1[(32*tmp3[1]-31):(32*tmp3[1]),j]*rbinom(32,1,0.4)
    #doppler[(32*tmp3[2]-31):(32*tmp3[2]),(50+j)] = doppler1[(32*tmp3[2]-31):(32*tmp3[2]),j]*rbinom(32,1,0.4)
    #doppler[(32*tmp4[1]-31):(32*tmp4[1]),(50+j)] = doppler1[(32*tmp4[1]-31):(32*tmp4[1]),j]*rbinom(32,1,0.4)
    #doppler[(32*tmp4[2]-31):(32*tmp4[2]),(50+j)] = doppler1[(32*tmp4[2]-31):(32*tmp4[2]),j]*rbinom(32,1,0.4)
  }
  
  #doppler3 = matrix(0,nrow=512,ncol=50)
  #for (i in 1:512) {
  #  doppler3[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-0.1)),0.5)
  #}
  
  for (j in 1:50) {
    tmp1 = sample(c(1,2,3,4),2)
    tmp2 = sample(c(5,6,7,8),2)
    #tmp3 = sample(c(9,10,11,12),2)
    #tmp4 = sample(c(13,14,15,16),2)
    doppler[(64*tmp1[1]-63):(64*tmp1[1]),(100+j)] = doppler1[(64*tmp1[1]-63):(64*tmp1[1]),j]*rbinom(64,1,0.6)
    doppler[(64*tmp1[2]-63):(64*tmp1[2]),(100+j)] = doppler1[(64*tmp1[2]-63):(64*tmp1[2]),j]*rbinom(64,1,0.2)
    doppler[(64*tmp2[1]-63):(64*tmp2[1]),(100+j)] = doppler1[(64*tmp2[1]-63):(64*tmp2[1]),j]*rbinom(64,1,0.6)
    doppler[(64*tmp2[2]-63):(64*tmp2[2]),(100+j)] = doppler1[(64*tmp2[2]-63):(64*tmp2[2]),j]*rbinom(64,1,0.2)
    #doppler[(32*tmp3[1]-31):(32*tmp3[1]),(100+j)] = doppler1[(32*tmp3[1]-31):(32*tmp3[1]),j]*rbinom(32,1,0.6)
    #doppler[(32*tmp3[2]-31):(32*tmp3[2]),(100+j)] = doppler1[(32*tmp3[2]-31):(32*tmp3[2]),j]*rbinom(32,1,0.2)
    #doppler[(32*tmp4[1]-31):(32*tmp4[1]),(100+j)] = doppler1[(32*tmp4[1]-31):(32*tmp4[1]),j]*rbinom(32,1,0.6)
    #doppler[(32*tmp4[2]-31):(32*tmp4[2]),(100+j)] = doppler1[(32*tmp4[2]-31):(32*tmp4[2]),j]*rbinom(32,1,0.2)
  }
  
  #doppler4 = matrix(0,nrow=512,ncol=50)
  #for (i in 1:512) {
  #  doppler4[i,]=rnorm(50,0.6+0.6*sqrt(time[i]*(1-time[i]))*sin(2.1*pi/(time[i]-0.1)),0.5)
  #}
  
  for (j in 1:50) {
    tmp1 = sample(c(1,2,3,4),1)
    tmp2 = sample(c(5,6,7,8),1)
    #tmp3 = sample(c(9,10,11,12),1)
    #tmp4 = sample(c(13,14,15,16),1)
    doppler[(64*tmp1[1]-63):(64*tmp1[1]),(150+j)] = doppler1[(64*tmp1[1]-63):(64*tmp1[1]),j]*rbinom(64,1,0.8)
    doppler[(64*tmp2[1]-63):(64*tmp2[1]),(150+j)] = doppler1[(64*tmp2[1]-63):(64*tmp2[1]),j]*rbinom(64,1,0.8)
  }
  
  doppler
}

true = c(rep(1,50),rep(2,50),rep(3,50),rep(4,50))

ccr <- function(result, true) {
  if(length(result) != length(true)) {warning('result vector and true vector has different size')}
  return(sum(true == result)/length(result))
}

ccr_kmeans_dopinten=vector()
rand_kmeans_dopinten=vector()

ccr_pam_dopinten=vector()
rand_pam_dopinten=vector()

ccr_funfem_dopinten=vector()
rand_funfem_dopinten=vector()

ccr_funhddc_dopinten=vector()
rand_funhddc_dopinten=vector()

for (k in 1:100) {
  doppler = doppler_intensity(k)
  
  MultiFeat_data_array = MultiFeat_data(doppler,8,2)
  
  funPCscore = MultiFeat_MFPCA(MultiFeat_data_array,0.9)
  
  ##Kmeans
  num_cl=4
  MCluster <- kmeans(funPCscore, centers= num_cl, iter.max = 100000)
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_kmeans_dopinten[k] <- ccr(true,cluster_result)
  rand_kmeans_dopinten[k] <- adjustedRandIndex(true,cluster_result)
  
  ##PAM
  num_cl=4
  MCluster <- pam(funPCscore, num_cl, metric="manhattan")
  cluster_result = best_match_factor(true,MCluster$cluster)
  
  ccr_pam_dopinten[k] <- ccr(true,cluster_result)
  rand_pam_dopinten[k] <- adjustedRandIndex(true,cluster_result)
  
  ## FunFEM
  fdabasis = create.bspline.basis(rangeval=c(0,512),norder=4,nbasis=21)
  fdatime = seq(1, 512, by=1)
  sim_doppler_fdafd = smooth.basis(fdatime, doppler, fdabasis)$fd
  res = funFEM(sim_doppler_fdafd,model='all',K=4)
  
  cluster_result = best_match_factor(true,res$cls)
  
  ccr_funfem_dopinten[k] <- ccr(true,cluster_result)
  rand_funfem_dopinten[k] <- adjustedRandIndex(true,cluster_result)
  
  ##FunHDDC
  fdabasis = create.fourier.basis(rangeval=c(0,512),period = 512)
  fdatime = seq(1, 512, by=1)
  sim_doppler_fdafd = smooth.basis(fdatime, doppler, fdabasis)$fd
  
  res <- funHDDC(sim_doppler_fdafd,model=c('AkjBkQkDk','AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'),K=4)
  
  cluster_result = best_match_factor(true,res$class)
  
  ccr_funhddc_dopinten[k] <- ccr(true,cluster_result)
  rand_funhddc_dopinten[k] <- adjustedRandIndex(true,cluster_result)
  
}

doppler_sim_intensity = cbind(ccr_kmeans_dopinten,rand_kmeans_dopinten,ccr_pam_dopinten,rand_pam_dopinten,ccr_funfem_dopinten,rand_funfem_dopinten,ccr_funhddc_dopinten,rand_funhddc_dopinten)
apply(doppler_sim_intensity,2,mean)
apply(doppler_sim_intensity,2,sd)

