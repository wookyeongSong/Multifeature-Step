convert_amount <- function(step) {
  data_amount<-matrix(nrow=nrow(step), ncol=ncol(step))
  for(i in 1:ncol(step)){
    data_amount[,i]= cumsum(step[,i])
  }
  return(data_amount/mean(data_amount[nrow(step),]))
}

convert_intensity <-function(step,K) {
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
  
  sort_data_intensity = matrix(NA,nrow=nrow(step), ncol=ncol(step))
  for(i in 1:ncol(data_intensity)) {
    sort_data_intensity[,i]<-sort(data_intensity[,i])
  }
  
  return(sort_data_intensity/mean(sort_data_intensity[nrow(step),]))
}


convert_pattern <- function(step,K) {
  new_patmatrix <- matrix(0,nrow= nrow(step),ncol=ncol(step))
  
  K=K
  pj=c(0:K)/K
  
  
  # for(i in 1:ncol(step)) {
  # s_pj= quantile(rank(step[which(step[,i]!=0),i]) , probs=pj) # quantile(c(1:nrow(step)) , probs=pj)
  # index=which(step[,i]!=0)
  # for(t in 1: length(index)){
  
  # new_patmatrix[index[t],i]=min(which(rank(step[which(step[,i]!=0),i] )[t] <= s_pj))-1
  
  # }
  
  # }
  
  
  
  
  
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
  
  sort_data_pattern=matrix(nrow=nrow(step),ncol=ncol(step))
  for(i in 1:ncol(step)) {
    sort_data_pattern[,i] = sort(data_pattern[,i])
  }
  
  return(data_pattern/mean(sort_data_pattern[nrow(step),]))
}

convert_nonorder_intensity <-function(step,K) {
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
  
  return(data_intensity)
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
  
  sort_data_intensity = matrix(NA,nrow=nrow(step), ncol=ncol(step))
  for(i in 1:ncol(data_intensity)) {
    sort_data_intensity[,i]<-sort(data_intensity[,i])
  }
  
  return(sort_data_intensity)
}


convert_nonstandard_pattern <- function(step,K) {
  new_patmatrix <- matrix(0,nrow= nrow(step),ncol=ncol(step))
  
  K=K
  pj=c(0:K)/K
  
  
  # for(i in 1:ncol(step)) {
  # s_pj= quantile(rank(step[which(step[,i]!=0),i]) , probs=pj) # quantile(c(1:nrow(step)) , probs=pj)
  # index=which(step[,i]!=0)
  # for(t in 1: length(index)){
  
  # new_patmatrix[index[t],i]=min(which(rank(step[which(step[,i]!=0),i] )[t] <= s_pj))-1
  
  # }
  
  # }
  
  
  
  
  
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

  
  return(data_pattern)
}


best_match_factor=function(x,y){
  x=factor(x)
  A=NA
  for(k in 1:length(permn(length(table(x)))) ){
    new_b= as.factor(y)
    newy=mapvalues(new_b, from = levels(new_b), to =  as.character( permn(length(table(x)))[k][[1]] ) )
    A=c(A,  length(which(newy==x) ) )
  }
  A=A[-1]
  newy=mapvalues(new_b, from = levels(new_b), to =  as.character( permn(length(table(x)))[which.max(A)][[1]] ) )
  return(as.numeric(as.character(newy)) )
}
