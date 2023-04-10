
####
#### LOADING, FACTOR AND COMMON COMPONENT
TFM_est=function(x,r,method='PE',tol=1e-4,maxiter=100){
  x <- aperm(x,c(2:length(dim(x)),1))
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]

  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  ans.init <- init.tensor(x,r,norm.true=TRUE)

  ddd=dd[-d]
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter+2)

  x.hat <- ttl(x.tnsr,lapply(ans.init$Q,t),d.seq)
  x.hat <- ttl(x.hat,ans.init$Q,d.seq)

  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q
  Ft <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
  Ft.all <- apply(Ft@data,c(1:(d-1)),sum)
  fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm


  if(method=="IE"){
    ans.Q <- ans.init$Q
  }else if(method=="PE"){
    for(i in 1:(d-1)){
      x.new <- aperm(ttl(x.tnsr,lapply(ans.init$Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
      ans.iter <- init.tensor(x.new,c(r[i],r[-i]),oneside.true=TRUE,norm.true=FALSE)
      ans.Q[[i]] <- ans.iter$Q[[1]]
    }
  }else{
    while((dis > tol) & (iiter < maxiter)){
      for(i in 1:(d-1)){
        if(method=="iPE"){
          x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
          ans.iter <- init.tensor(x.new,c(r[i],r[-i]),oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else if(method=="HUBER"){
          weight<-w.T(x,ans.Q)
          x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
          ans.iter <- huber.tensor(x.new,c(r[i],r[-i]),weight,oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else{
          stop('Wrong method !')
        }
      }
        ddd=dd[-d]

        x.hat <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
        x.hat <- ttl(x.hat,ans.Q,d.seq)
        Ft <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
        Ft.all <- apply(Ft@data,c(1:(d-1)),sum)

        fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
        dis <- abs(fnorm.resid[iiter+1] - fnorm.resid[iiter])
        if(iiter==1){
          Qfirst <- ans.Q
          x.hat.first <- x.hat@data
        }
        iiter <- iiter + 1
      }
    }

  fnorm.resid <- fnorm.resid[fnorm.resid != 0]
  fnorm.resid <- fnorm.resid^2
  x0 <- matrix(x,prod(dd[-d]))
  x0 <- t(scale(t(x0),scale=FALSE) )
  x0 <- array(x0,dd)
  return(list("Ft"=aperm(Ft@data,c(d,1:(d-1))),"Ft.all"=Ft.all,"Q"=ans.Q,"x.hat"=aperm(x.hat@data,c(d,1:(d-1))),"niter"=iiter,"fnorm.resid"=fnorm.resid[iiter]))
}



####
#### FACTOR NUMBER
TFM_FN = function(x,r=NULL,method='PE',tol=1e-4,maxiter=100){
  x <- aperm(x,c(2:length(dim(x)),1))
  dd <- dim(x)
  d <- length(dd) # d >= 2
  d.seq <- 1:(d-1)
  n <- dd[d]
  x.tnsr <- as.tensor(x)
  tnsr.norm <- fnorm(x.tnsr)
  factor.num <- array(NA, c(d-1,maxiter+1))
  ddd<-dd[-d]

  if(is.null(r)){
    for (i in 1:(d-1)){
      r[i] = ceiling(dd[i]/3)
    }
  }
  ans.init <- init.tensor(x,r,norm.true=TRUE)
  for(i in 1:(d-1)){
      factor.num[i,1]=tensor.er(ans.init$lambda[[i]],ddd[i],ddd[-i],n)
  }
  r.init <- factor.num[,1]
  iiter <- 1
  dis <- 1
  fnorm.resid <- rep(0,maxiter+1)
  x.hat <- ttl(x.tnsr,lapply(ans.init$Q,t),d.seq)
  x.hat <- ttl(x.hat,ans.init$Q,d.seq)

  fnorm.resid[1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
  ans.Q <- ans.init$Q

  if(method=='IE'){
    iiter=iiter+1
    factor.num[,iiter]=r
  }else {
    for(j in 1:(d-1)){r[j] = min(dd[j],r[j]+1)}
    while((dis > tol) & (iiter < maxiter)){
      for(i in 1:(d-1)){
        if(method=="PE"){
          x.new <- aperm(ttl(x.tnsr,lapply(ans.init$Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
          ans.iter <- init.tensor(x.new,c(r[i],r[-i]),oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else if(method=="iPE"){
          x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
          ans.iter <- init.tensor(x.new,c(r[i],r[-i]),oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else if(method=="HUBER"){
          weight<-w.T(x,ans.Q)
          x.new <- aperm(ttl(x.tnsr,lapply(ans.Q[-i],t),ms=d.seq[-i])@data,c(i,d.seq[-i],d))
          ans.iter <- huber.tensor(x.new,c(r[i],r[-i]),weight,oneside.true=TRUE,norm.true=FALSE)
          ans.Q[[i]] <- ans.iter$Q[[1]]
        }else{
          stop('Wrong estimation method input !')
        }

        factor.num[i,1+iiter]=tensor.er(ans.iter$lambda[[1]],ddd[i],ddd[-i],n)
        r[i]=factor.num[i,1+iiter]

      }

      x.hat <- ttl(x.tnsr,lapply(ans.Q,t),d.seq)
      x.hat <- ttl(x.hat,ans.Q,d.seq)

      fnorm.resid[iiter+1] <- fnorm(x.tnsr-x.hat)/tnsr.norm
      dis <- abs(fnorm.resid[iiter+1] - fnorm.resid[iiter])
      if(iiter==1){
        Qfirst <- ans.Q
        x.hat.first <- x.hat@data
      }
      iiter <- iiter + 1
    }
  }
  # factor.num[,,maxiter]=factor.num[,,iiter]
  factor.num[,maxiter]=factor.num[,iiter-1]
  fnorm.resid <- fnorm.resid[fnorm.resid != 0]

  # label the factor number path
  path = t(factor.num[,1:(iiter)])
  path.rowname = c()
  for(ii in 1:iiter){path.rowname <- c(path.rowname,paste('iteration ',ii-1,sep=''))}
  path.colname = c()
  for(ii in 1:(d-1)){path.colname <- c(path.colname,paste('mode ',ii,sep=''))}
  rownames(path)=path.rowname
  colnames(path)=path.colname

  path
  factor.num[,maxiter]
  # return(list("path"=t(factor.num[,penalty,1:(iiter)]),"factor.num"=factor.num[,penalty,maxiter]))
  return(list("path"=path,"factor.num"=factor.num[,maxiter]))
}




###LOADING, FACTOR AND COMMON COMPONENT HELPER
init.tensor <- function(x,r,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # TIPUP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_d * n
  # if oneside.true==TRUE, then only compute the one sided column space,
  # not the other sides, this option is useful for the iterative method
  dd <- dim(x)
  d <- length(dd) # d >= 2
  n <- dd[d]
  dd.prod <- prod(dd) / n
  x.matrix <- matrix(x,ncol=n)
  if(oneside.true==TRUE){
    k=1
  } else{
    k=d-1
  }
  ans.M <- ans.Q <- ans.lambda <- NULL

  for(i in 1:k){
    x.mat <- array(x.matrix,c(dd[-d],n))
    Omega <- tensor(x.mat,x.mat,c(1:d)[-i],c(1:d)[-i])/(n*dd.prod/dd[i])
    M.temp <- Omega %*% t(Omega)
    M.temp <- M.temp / dd.prod * dd[i]

    ans.eig <- eigen(M.temp)
    ans.M <- c(ans.M,list(M.temp))
    ans.Q <- c(ans.Q,list(ans.eig$vectors[,1:r[i],drop=FALSE]))
    ans.lambda <- c(ans.lambda,list(ans.eig$values))
  }
  norm.percent <- NULL
  x.hat <- NULL
  if(norm.true==TRUE){
    x.tnsr <- as.tensor(x)
    #x.hat <- get.hat(x.tnsr,ans.Q,1:k)
    x.hat <- ttl(x.tnsr,lapply(ans.Q,t),1:k)
    x.hat <- ttl(x.hat,ans.Q,1:k)
    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("M"=ans.M,"Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}



#### FACTOR NUMBER HELPER
tensor.er<-function(reigen,p1,p2,n){
  if(length(p2)>1){
    p2=prod(p2)
  }
  p=p1*p2
  m1=ceiling(p1/3)

  lambda.p1<-reigen[p1:1]
  #cumlambda.p1<-cumsum(lambda.p1)
  #cumlambda.p1<-cumlambda.p1[(p1-1):1]

  #ratio
  #ratio<-lambda.p1[(p1-1):(p1-m1)]/lambda.p1[p1:(p1-m1+1)]
  ratio<-reigen[1:m1]/reigen[2:(1+m1)]
  #factor.p1<-which.min(ratio)
  factor.p1<-which.max(ratio)
  #which.max(eigval1[1:k1]/(eigval1[2:(1+k1)]))

  factor.p1
}



###
### HUBER HELPER

###tau
tau.est<-function(Y,A){
  t=c()
  dd <- dim(Y)
  d <- length(dd) # d >= 2
  n <- dd[d]
  dd.prod <- prod(dd) / n
  d.seq<-1:(d-1)
  Y.tnsr <- as.tensor(Y)
  Ft <- ttl(Y.tnsr,lapply(A,t),d.seq)
  Gt <- matrix((Y.tnsr-ttl(Ft,A,d.seq))@data,dd.prod,n)
  for (i in 1:n) {
    t[i]=fnorm(as.tensor(Gt[,i]))
  }
  t=median(t)
  return(t)
}


##weights
w.T<-function(Y,A){
  w=c()
  t=c()
  dd <- dim(Y)
  d <- length(dd) # d >= 2
  n <- dd[d]
  dd.prod <- prod(dd) / n
  d.seq<-1:(d-1)
  Y.tnsr <- as.tensor(Y)
  Ft <- ttl(Y.tnsr,lapply(A,t),d.seq)
  Gt <- matrix((Y.tnsr-ttl(Ft,A,d.seq))@data,dd.prod,n)
  for (i in 1:n) {
    t[i]=fnorm(as.tensor(Gt[,i]))
  }
  tau=median(t)
  for (i in 1:n) {
    if(t[i]<=tau){
      w[i]=1/2
    }else{
      w[i]=tau/2/t[i]
    }
  }
  return(w)
}


huber.tensor <- function(x,r,weight,oneside.true=FALSE,norm.true=FALSE){
  # x: tensor of any dimension
  # TIPUP initialization (one step)
  # x: d1 * d2 * d3 * ... * d_d * n
  # if oneside.true==TRUE, then only compute the one sided column space,
  # not the other sides, this option is useful for the iterative method
  dd <- dim(x)
  d <- length(dd) # d >= 2
  n <- dd[d]
  dd.prod <- prod(dd) / n
  x.matrix <- matrix(x,ncol=n)
  x.weight <- matrix(0,dd.prod,n)
  for (i in 1:n){
    x.weight[,i]<-sqrt(weight[i])*x.matrix[,i]
  }
  if(oneside.true==TRUE){
    k=1
  } else{
    k=d-1
  }
  ans.M <- ans.Q <- ans.lambda <- NULL


  for(i in 1:k){
    x.mat <- array(x.weight,c(dd[-d],n))
    Omega <- tensor(x.mat,x.mat,c(1:d)[-i],c(1:d)[-i])/(n)
    M.temp <- Omega %*% t(Omega)
    M.temp <- M.temp / dd.prod * dd[i]

    ans.eig <- eigen(M.temp)
    ans.M <- c(ans.M,list(M.temp))
    ans.Q <- c(ans.Q,list(ans.eig$vectors[,1:r[i],drop=FALSE]))
    ans.lambda <- c(ans.lambda,list(ans.eig$values))
  }
  norm.percent <- NULL
  x.hat <- NULL
  if(norm.true==TRUE){
    x.tnsr <- as.tensor(x)
    #x.hat <- get.hat(x.tnsr,ans.Q,1:k)
    x.hat <- ttl(x.tnsr,lapply(ans.Q,t),1:k)
    x.hat <- ttl(x.hat,ans.Q,1:k)
    norm.percent <- fnorm(x.tnsr-x.hat)/fnorm(x.tnsr)
    x.hat <- x.hat@data
  }
  list("M"=ans.M,"Q"=ans.Q,"lambda"=ans.lambda,"norm.percent"=norm.percent,"x.hat"=x.hat)
}

