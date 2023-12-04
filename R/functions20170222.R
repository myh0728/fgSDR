################################################
##                                            ##
##  Gauss elimination to "row"-echelon forms  ##
##                                            ##
################################################

Gauss_rref <- function(A)
{
  M <- A
  m <- nrow(A)
  n <- ncol(A)
  k_max <- min(m,n)
  npivot <- 0

  for (k in 1:k_max)
  {
    i_max <- which(abs(M[k:m,k])==max(abs(M[k:m,k])))[1]+k-1

    if (M[i_max,k]==0)
      next

    M[c(k,i_max),] <- M[c(i_max,k),]
    npivot <- npivot+1

    for (i in 1:m)
    {
      if (i==npivot)
      next

      M[i,] <- M[i,]-M[npivot,]*M[i,k]/M[npivot,k]
    }

    M[npivot,] <- M[npivot,]/M[npivot,k]
  }

  return(M)
}

########################
##                    ##
##  Kernel functions  ##
##                    ##
########################

# kernel function

K2_Epanechnikov <- function(u)
{
  3/4*(1-u^2)*(abs(u)<=1)
}

K2_Biweight <- function(u)
{
  15/16*(1-u^2)^2*(abs(u)<=1)
}

K4_Biweight <- function(u)
{
  (105/64)*(1-3*(u^2))*((1-u^2)^2)*(abs(u)<=1)
}

dK4_Biweight <- function(u)
{
  (105/32)*u*(1-u^2)*(9*u^2-5)*(abs(u)<=1)
}

dK2_Biweight <- function(u)
{
  -15/4*u*(1-u^2)*(abs(u)<=1)
}

###########################################
##                                       ##
##  Local constant smoothing estimators  ##
##                                       ##
###########################################

### kernel matrix

# Xi: covariates, c(n_s,p_s) matrix
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing matrix at different evaluation points, c(n_s,k_s) matrix

K.matrix <- function(Xi,x,K,h)
{
  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]
  k_s <- dim(x)[1]

  Xh <- t(t(Xi)/h)
  xh <- t(t(x)/h)
  KXh <- K(Xh[rep(1:n_s,times=k_s),]-xh[rep(1:k_s,each=n_s),])
  dim(KXh) <- c(n_s*k_s,p_s)
  KXh <- apply(KXh,1,prod)
  dim(KXh) <- c(n_s,k_s)
  return(KXh)
}

### kernel matrix with leave-one-out cross validation

# Xi: covariates, c(n_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing matrix at observed points, c(n_s,n_s) matrix

K.matrix.cv <- function(Xi,K,h)
{
  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]

  Xh <- t(t(Xi)/h)
  KXh <- K(Xh[rep(1:n_s,times=n_s),]-Xh[rep(1:n_s,each=n_s),])
  dim(KXh) <- c(n_s*n_s,p_s)
  KXh <- apply(KXh,1,prod)
  dim(KXh) <- c(n_s,n_s)
  diag(KXh) <- 0
  return(KXh)
}

### kernel smoother

# Xi: covariates, c(n_s,p_s) matrix
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoother at different evaluation points, c(n_s,k_s) matrix

K.smoother <- function(Xi,x,K,h)
{
  K.m <- K.matrix(Xi,x,K,h)
  Nn <- t(K.m)
  Nd <- colSums(K.m)
  smoother <- t(Nn*(Nd!=0)/(Nd+(Nd==0)))
  return(smoother)
}

### kernel smoother with leave-one-out cross validation

# Xi: covariates, c(n_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoother at observed points, c(n_s,n_s) matrix

K.smoother.cv <- function(Xi,K,h)
{
  K.m <- K.matrix.cv(Xi,K,h)
  Nn <- t(K.m)
  Nd <- colSums(K.m)
  smoother <- t(Nn*(Nd!=0)/(Nd+(Nd==0)))
  return(smoother)
}

### kernel smoothing estimator for density

# Xi: covariates, c(n_s,p_s) matrix
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for density at different evaluation points, c(k_s) vector

K.density <- function(Xi,x,K,h)
{
  K.m <- K.matrix(Xi,x,K,h)
  dhat <- colSums(K.m)
  dhat <- dhat*(dhat>=0)
  return(dhat)
}

### kernel smoothing estimator for density with leave-one-out cross validation

# Xi: covariates, c(n_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for density at observed points, c(n_s) vector

K.density.cv <- function(Xi,K,h)
{
  K.m <- K.matrix.cv(Xi,K,h)
  dhat <- colSums(K.m)
  dhat <- dhat*(dhat>=0)
  return(dhat)
}

### kernel smoothing estimator for mean

# Xi: covariates, c(n_s,p_s) matrix
# Yik: responses indexed by another variable, c(n_s,y_k) matrix
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for mean at different evaluation points and index values, c(k_s,y_k) matrix

K.mean <- function(Xi,Yik,x,K,h)
{
  K.s <- K.smoother(Xi,x,K,h)
  mhat <- t(K.s) %*% Yik
  return(mhat)
}

### kernel smoothing estimator for mean with leave-one-out cross validation

# Xi: covariates, c(n_s,p_s) matrix
# Yik: responses indexed by another variable, c(n_s,y_k) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for mean at observed points and index values, c(n_s,y_k) matrix

K.mean.cv <- function(Xi,Yik,K,h)
{
  K.s <- K.smoother.cv(Xi,K,h)
  mhat <- t(K.s) %*% Yik
  return(mhat)
}

### kernel smoothing estimator for distribution

# Xi: covariates, c(n_s,p_s) matrix
# Dik: indicator responses at different evaluation points, c(n_s,y_k) matrix
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for distribution at different evaluation points, c(k_s,y_k) matrix

K.dist <- function(Xi,Dik,x,K,h)
{
  K.s <- K.smoother(Xi,x,K,h)
  Fhat <- t(K.s) %*% Dik
  Fhat <- Fhat*(Fhat>0)*(Fhat<1)+(Fhat>=1)
  return(Fhat)
}


### kernel smoothing estimator for distribution with leave-one-out cross validation

# Xi: covariates, c(n_s,p_s) matrix
# Dik: indicator responses at different evaluation points, c(n_s,y_k) matrix
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for distribution at observed points, c(n_s,y_k) matrix

K.dist.cv <- function(Xi,Dik,K,h)
{
  K.s <- K.smoother.cv(Xi,K,h)
  Fhat <- t(K.s) %*% Dik
  Fhat <- Fhat*(Fhat>0)*(Fhat<1)+(Fhat>=1)
  return(Fhat)
}


### kernel smoothing estimator for distribution with discrete covariate

# Xi: covariates, c(n_s,p_s) matrix
# Dik: indicator responses at different evaluation points, c(n_s,y_k) matrix
# Zi: univariate discrete variable, c(n_s) vector
# x: evaluation points, c(k_s,p_s) matrix
# z: evaluation points, c(k_s) vector
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for distribution at different evaluation points, c(k_s,y_k) matrix

K.dist.z <- function(Xi,Dik,Zi,x,z,K,h)
{
  K.m <- K.matrix(Xi,x,K,h)
  CZ <- outer(Zi,z,FUN="==")
  nr <- K.m*CZ
  Nn <- t(nr) %*% Dik
  Nd <- colSums(nr)
  Fhat <- Nn*(Nd!=0)/(Nd+(Nd==0))
  Fhat <- Fhat*(Fhat>0)*(Fhat<1)+(Fhat>=1)
  return(Fhat)
}

### kernel smoothing estimator for distribution with discrete covariate and leave-one-out cross validation

# Xi: covariates, c(n_s,p_s) matrix
# Dik: indicator responses at different evaluation points, c(n_s,y_k) matrix
# Zi: univariate discrete variable, c(n_s) vector
# K: kernel function
# h: bandwidth

# value: kernel smoothing estimator for distribution at different evaluation points, c(k_s,y_k) matrix

K.dist.z.cv <- function(Xi,Dik,Zi,K,h)
{
  K.m <- K.matrix.cv(Xi,K,h)
  CZ <- outer(Zi,Zi,FUN="==")
  nr <- K.m*CZ
  Nn <- t(nr) %*% Dik
  Nd <- colSums(nr)
  Fhat <- Nn*(Nd!=0)/(Nd+(Nd==0))
  Fhat <- Fhat*(Fhat>0)*(Fhat<1)+(Fhat>=1)
  return(Fhat)
}

### Kaplan-Meier type estimator for conditional survival

# Xi: covariates, c(n_s,p_s) matrix
# Yi: observed time, c(n_s) vector
# di: non-censoring status, c(n_s) vector
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth

KMNWE <- function(Xi,Yi,di,x,K,h)
{
  t_points <- sort(unique(Yi[di==1]))
  l_s <- length(t_points)

  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]
  k_s <- dim(x)[1]

  Xh <- t(t(Xi)/h)
  xh <- t(t(x)/h)
  KXh <- K(Xh[rep(1:n_s,times=k_s),]-xh[rep(1:k_s,each=n_s),])
  dim(KXh) <- c(n_s*k_s,p_s)
  KXh <- apply(KXh,1,prod)
  dim(KXh) <- c(n_s,k_s)

  Rhat <- t(KXh) %*% outer(Yi,t_points,FUN=">=")
  dH1hat <- t(KXh) %*% (outer(Yi,t_points,FUN="==")*di)
  dLambda <- dH1hat*(Rhat!=0)/(Rhat+(Rhat==0))
  dLambda <- dLambda*(dLambda>=0)*(dLambda<=1)+(dLambda>1)
  Shat <- apply(1-dLambda,1,cumprod)

  return(list(jump=t_points,survival=Shat))
}

### Kaplan-Meier type estimator for conditional survival with weights (for bootstrapping)

# Xi: covariates, c(n_s,p_s) matrix
# Yi: observed time, c(n_s) vector
# di: non-censoring status, c(n_s) vector
# x: evaluation points, c(k_s,p_s) matrix
# K: kernel function
# h: bandwidth
# w: weights with sum equal to 1, c(n_s) vector

KMNWE.w <- function(Xi,Yi,di,x,K,h,w)
{
  t_points <- sort(unique(Yi[di==1]))
  l_s <- length(t_points)

  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]
  k_s <- dim(x)[1]

  Xh <- t(t(Xi)/h)
  xh <- t(t(x)/h)
  KXh <- K(Xh[rep(1:n_s,times=k_s),]-xh[rep(1:k_s,each=n_s),])
  dim(KXh) <- c(n_s*k_s,p_s)
  KXh <- apply(KXh,1,prod)
  dim(KXh) <- c(n_s,k_s)

  Rhat <- t(KXh) %*% (outer(Yi,t_points,FUN=">=")*w)
  dH1hat <- t(KXh) %*% (outer(Yi,t_points,FUN="==")*di*w)
  dLambda <- dH1hat*(Rhat!=0)/(Rhat+(Rhat==0))
  dLambda <- dLambda*(dLambda>=0)*(dLambda<=1)+(dLambda>1)
  Shat <- apply(1-dLambda,1,cumprod)

  return(list(jump=t_points,survival=Shat))
}

##############################################
##                                          ##
##  Cross-validated least squares criteria  ##
##                                          ##
##############################################

# leave-one-out cross validation criterion for mean regression

# Xi: covariates, c(n_s,p_s) matrix
# Yi: response, c(n_s) vector
# K: kernel function
# h: bandwidth

# value: estimated prediction risk, real value

CV.mean <- function(Xi,Yi,K,h)
{
  PR <- mean((Yi-K.mean.cv(Xi,Yi,K,h))^2)
  return(PR)
}

# leave-one-out cross validation criterion for distribution regression

# Xi: covariates, c(n_s,p_s) matrix
# Dik: responses at different evaluation points, c(n_s,y_k) matrix
# w: weights, c(y_k) vector
# K: kernel function
# h: bandwidth

# value: estimated prediction risk, real value

CV.dist <- function(Xi,Dik,w,K,h)
{
 APR <- mean(colSums(t((Dik-K.dist.cv(Xi,Dik,K,h))^2)*w))
 return(APR)
}

# integrated sum of square with score and information under cross-validation for distribution SDR regression
# (on Stiefel manifold)

# Xij: difference between covariates, c(n_s,n_s,p_s) matrix
# B: basis matrix, c(p_s,d_s) matrix
# Dik: responses at data points, c(n_s,k_s) matrix
# w: weights with respect to integral, c(k_s) vector
# K: kernel function
# dK: 1st derivative of kernel function
# h: bandwidth

# value:
# $value: integration sum of square, real value
# $score: sample score vectors
# $information: information matrix

CV.dist.nt <- function(Xij,B,Dik,w,K,dK,h)
{
  n_s <- dim(Xij)[1]
  p_s <- dim(Xij)[3]
  d_s <- dim(B)[2]
  k_s <- dim(Dik)[2]

  dim(Xij) <- c(n_s^2,p_s)
  XBhi <- (Xij %*% B)/h
  KXBhij.v <- K(XBhi)
  dim(KXBhij.v) <- c(n_s,n_s,d_s)
  KXBhij <- apply(KXBhij.v,c(1,2),prod)
  diag(KXBhij) <- 0

  Nn <- KXBhij %*% Dik
  Dn <- rowSums(KXBhij)
  Fhat <- Nn*(Dn!=0)/(Dn+(Dn==0))
  Fhat <- Fhat*(Fhat>0)*(Fhat<1)+(Fhat>=1)

  cv <- mean(colSums(t((Dik-Fhat)^2)*w))

  dKXBhij.v <- dK(XBhi)/h
  dim(dKXBhij.v) <- c(n_s,n_s,d_s)
  if (d_s>=2)
  {
    for (k.d in 1:d_s)
    {
      dKXBhij.v[,,k.d] <- dKXBhij.v[,,k.d]*apply(KXBhij.v[,,-k.d],c(1,2),prod)
    }
  }

  F1hat <- array(0,c(n_s,k_s,p_s,d_s))

  for (k1 in 1:p_s)
  {
    for (k2 in 1:d_s)
    {
      dBK <- dKXBhij.v[,,k2]*Xij[,k1]
      F1hat[,,k1,k2] <- (dBK %*% Dik-Fhat*rowSums(dBK))*(Dn!=0)/(Dn+(Dn==0))
    }
  }

  Shati <- F1hat*as.vector(Dik-Fhat)*(-2)
  Shati <- apply(aperm(Shati,c(2,1,3,4))*w,c(2,3,4),sum)

  Vhat <- array(0,c(p_s*d_s,p_s*d_s))

  dim(F1hat) <- c(n_s,k_s,p_s*d_s)

  for (k1 in 1:(p_s*d_s))
  {
    for (k2 in 1:(p_s*d_s))
    {
      Vhat[k1,k2] <- 2*sum(colMeans(F1hat[,,k1]*F1hat[,,k2])*w)
    }
  }

  results <- list(value=cv,score=Shati,information=Vhat)
  return(results)
}

# integrated sum of square with score and information under cross-validation for distribution SDR regression
# (in local coordinate)

# Xij: difference between covariates, c(n_s,n_s,p_s) matrix
# B: basis matrix, c(p_s-d_s,d_s) matrix
# Dik: responses at data points, c(n_s,k_s) matrix
# w: weights with respect to integral, c(k_s) vector
# K: kernel function
# dK: 1st derivative of kernel function
# h: bandwidth

# value:
# $value: integration sum of square, real value
# $score: sample score vectors
# $information: information matrix

CV.dist.nt.local <- function(Xij,B,Dik,w,K,dK,h)
{
  n_s <- dim(Xij)[1]
  p_s <- dim(Xij)[3]
  d_s <- dim(B)[2]
  k_s <- dim(Dik)[2]

  dim(Xij) <- c(n_s^2,p_s)
  XBhi <- (Xij %*% rbind(diag(d_s),B))/h
  KXBhij.v <- K(XBhi)
  dim(KXBhij.v) <- c(n_s,n_s,d_s)
  KXBhij <- apply(KXBhij.v,c(1,2),prod)
  diag(KXBhij) <- 0

  Nn <- KXBhij %*% Dik
  Dn <- rowSums(KXBhij)
  Fhat <- Nn*(Dn!=0)/(Dn+(Dn==0))
  Fhat <- Fhat*(Fhat>0)*(Fhat<1)+(Fhat>=1)

  cv <- mean(colSums(t((Dik-Fhat)^2)*w))

  dKXBhij.v <- dK(XBhi)/h
  dim(dKXBhij.v) <- c(n_s,n_s,d_s)
  if (d_s>=2)
  {
    for (k.d in 1:d_s)
    {
      dKXBhij.v[,,k.d] <- dKXBhij.v[,,k.d]*apply(KXBhij.v[,,-k.d],c(1,2),prod)
    }
  }

  F1hat <- array(0,c(n_s,k_s,p_s-d_s,d_s))

  for (k1 in (d_s+1):p_s)
  {
    for (k2 in 1:d_s)
    {
      dBK <- dKXBhij.v[,,k2]*Xij[,k1]
      F1hat[,,k1-d_s,k2] <- (dBK %*% Dik-Fhat*rowSums(dBK))*(Dn!=0)/(Dn+(Dn==0))
    }
  }

  Shati <- F1hat*as.vector(Dik-Fhat)*(-2)
  Shati <- apply(aperm(Shati,c(2,1,3,4))*w,c(2,3,4),sum)

  Vhat <- array(0,c((p_s-d_s)*d_s,(p_s-d_s)*d_s))

  dim(F1hat) <- c(n_s,k_s,(p_s-d_s)*d_s)

  for (k1 in 1:((p_s-d_s)*d_s))
  {
    for (k2 in 1:((p_s-d_s)*d_s))
    {
      Vhat[k1,k2] <- 2*sum(colMeans(F1hat[,,k1]*F1hat[,,k2])*w)
    }
  }

  results <- list(value=cv,score=Shati,information=Vhat)
  return(results)
}

# leave-one-out cross validation criterion for distribution regression with discrete covariate

# Xi: covariates, c(n_s,p_s) matrix
# Dik: responses at different evaluation points, c(n_s,y_k) matrix
# Zi: univariate discrete variable, c(n_s) vector
# w: weights, c(y_k) vector
# K: kernel function
# h: bandwidth

# value: estimated prediction risk, real value

CV.dist.z <- function(Xi,Dik,Zi,w,K,h)
{
  APR <- mean(colSums(t((Dik-K.dist.z.cv(Xi,Dik,Zi,K,h))^2)*w))
  return(APR)
}

#############################
##                         ##
##  Generalized Cox model  ##
##                         ##
#############################

# leave-one-out cross validation criterion for survival regression under multivariate baseline PH model

# y_k is the length of observed times, the y_k dimension should be in the order of observed times

# Nik: observed failure process, c(n_s,y_k) matrix
# dNik: incident failure indicator, c(n_s,y_k) matrix
# YeIik: risk factor, c(n_s,y_k) matrix
# S0hatik: conditional mean of risk factor, c(n_s,y_k) matrix
# Xi: covariates, c(n_s,p_s) matrix
# w: weights, c(y_k) vector
# K: kernel function
# h: bandwidth

# value: estimated prediction risk, real value

CV.surv.MBPH <- function(Nik,dNik,YeIik,S0hatik,Xi,w,K,h)
{
  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]

  KXh <- K(Xi[rep(1:n_s,times=n_s),]-Xi[rep(1:n_s,each=n_s),]/h)
  dim(KXh) <- c(n_s,n_s,p_s)
  KXh <- apply(KXh,c(1,2),prod)
  diag(KXh) <- 0

  Nn <- KXh %*% dNik
  Dn <- KXh %*% YeIik
  dLambda <- Nn*(Dn!=0)/(Dn+(Dn==0))
  mik <- t(apply(S0hatik*dLambda,1,cumsum))
  mik <- mik*(mik>0)*(mik<=1)+(mik>1)
  APR <- mean(colSums(t((Nik-mik)^2)*w))

  return(APR)
}

# leave-one-out cross validation criterion for survival regression under multivariate baseline PH model
# with scores and information matrix (on Stiefel manifold)

# y_k is the length of observed times, the y_k dimension should be in the order of observed times

# Nik: observed failure process, c(n_s,y_k) matrix
# dNik: incident failure indicator, c(n_s,y_k) matrix
# YeIik: risk factor, c(n_s,y_k) matrix
# S0hatik: conditional mean of risk factor, c(n_s,y_k) matrix
# Xij: difference between covariates, c(n_s,n_s,p_s) matrix
# B: basis matrix, c(p_s,d_s) matrix
# w: weights, c(y_k) vector
# K: kernel function
# dK: 1st derivative of kernel function
# h: bandwidth

# value: estimated prediction risk, real value

CV.surv.MBPH.nt <- function(Nik,dNik,YeIik,S0hatik,Xij,B,w,K,dK,h)
{
  n_s <- dim(Xij)[1]
  p_s <- dim(Xij)[3]
  d_s <- dim(B)[2]
  y_k <- dim(Nik)[2]

  dim(Xij) <- c(n_s^2,p_s)
  XBhi <- (Xij %*% B)/h
  KXBhij.v <- K(XBhi)
  dim(KXBhij.v) <- c(n_s,n_s,d_s)
  KXBhij <- apply(KXBhij.v,c(1,2),prod)
  diag(KXBhij) <- 0

  Nn <- KXBhij %*% dNik
  Dn <- KXBhij %*% YeIik
  dLambda <- Nn*(Dn!=0)/(Dn+(Dn==0))
  mik <- t(apply(S0hatik*dLambda,1,cumsum))
  mik <- mik*(mik>0)*(mik<=1)+(mik>1)
  cv <- mean(colSums(t((Nik-mik)^2)*w))

  dKXBhij.v <- dK(XBhi)/h
  dim(dKXBhij.v) <- c(n_s,n_s,d_s)
  if (d_s>=2)
  {
    for (k.d in 1:d_s)
    {
      dKXBhij.v[,,k.d] <- dKXBhij.v[,,k.d]*apply(KXBhij.v[,,-k.d],c(1,2),prod)
    }
  }

  L1hat <- array(0,c(n_s,y_k,p_s,d_s))

  for (k1 in 1:p_s)
  {
    for (k2 in 1:d_s)
    {
      dBK <- dKXBhij.v[,,k2]*Xij[,k1]
      L1hat[,,k1,k2] <- t(apply(S0hatik*(dBK %*% dNik-dLambda*(dBK %*% YeIik))*(Dn!=0)/(Dn+(Dn==0)),1,cumsum))
    }
  }

  Shati <- L1hat*as.vector(Nik-mik)*(-2)
  Shati <- apply(aperm(Shati,c(2,1,3,4))*w,c(2,3,4),sum)

  Vhat <- array(0,c(p_s*d_s,p_s*d_s))

  dim(L1hat) <- c(n_s,y_k,p_s*d_s)

  for (k1 in 1:(p_s*d_s))
  {
    for (k2 in 1:(p_s*d_s))
    {
      Vhat[k1,k2] <- 2*sum(colMeans(L1hat[,,k1]*L1hat[,,k2])*w)
    }
  }

  results <- list(value=cv,score=Shati,information=Vhat)
  return(results)
}

# leave-one-out cross validation criterion for survival regression under multivariate baseline PH model
# with scores and information matrix (in local coordinate)

# y_k is the length of observed times, the y_k dimension should be in the order of observed times

# Nik: observed failure process, c(n_s,y_k) matrix
# dNik: incident failure indicator, c(n_s,y_k) matrix
# YeIik: risk factor, c(n_s,y_k) matrix
# S0hatik: conditional mean of risk factor, c(n_s,y_k) matrix
# Xij: difference between covariates, c(n_s,n_s,p_s) matrix
# B: basis matrix, c(p_s-d_s,d_s) matrix
# w: weights, c(y_k) vector
# K: kernel function
# dK: 1st derivative of kernel function
# h: bandwidth

# value: estimated prediction risk, real value

CV.surv.MBPH.nt.local <- function(Nik,dNik,YeIik,S0hatik,Xij,B,w,K,dK,h)
{
  n_s <- dim(Xij)[1]
  p_s <- dim(Xij)[3]
  d_s <- dim(B)[2]
  y_k <- dim(Nik)[2]

  dim(Xij) <- c(n_s^2,p_s)
  XBhi <- (Xij %*% rbind(diag(d_s),B))/h
  KXBhij.v <- K(XBhi)
  dim(KXBhij.v) <- c(n_s,n_s,d_s)
  KXBhij <- apply(KXBhij.v,c(1,2),prod)
  diag(KXBhij) <- 0

  Nn <- KXBhij %*% dNik
  Dn <- KXBhij %*% YeIik
  dLambda <- Nn*(Dn!=0)/(Dn+(Dn==0))
  mik <- t(apply(S0hatik*dLambda,1,cumsum))
  mik <- mik*(mik>0)*(mik<=1)+(mik>1)
  cv <- mean(colSums(t((Nik-mik)^2)*w))

  dKXBhij.v <- dK(XBhi)/h
  dim(dKXBhij.v) <- c(n_s,n_s,d_s)
  if (d_s>=2)
  {
    for (k.d in 1:d_s)
    {
      dKXBhij.v[,,k.d] <- dKXBhij.v[,,k.d]*apply(KXBhij.v[,,-k.d],c(1,2),prod)
    }
  }

  L1hat <- array(0,c(n_s,y_k,p_s-d_s,d_s))

  for (k1 in (d_s+1):p_s)
  {
    for (k2 in 1:d_s)
    {
      dBK <- dKXBhij.v[,,k2]*Xij[,k1]
      L1hat[,,k1-d_s,k2] <- t(apply(S0hatik*(dBK %*% dNik-dLambda*(dBK %*% YeIik))*(Dn!=0)/(Dn+(Dn==0)),1,cumsum))
    }
  }

  Shati <- L1hat*as.vector(Nik-mik)*(-2)
  Shati <- apply(aperm(Shati,c(2,1,3,4))*w,c(2,3,4),sum)

  Vhat <- array(0,c((p_s-d_s)*d_s,(p_s-d_s)*d_s))

  dim(L1hat) <- c(n_s,y_k,(p_s-d_s)*d_s)

  for (k1 in 1:((p_s-d_s)*d_s))
  {
    for (k2 in 1:((p_s-d_s)*d_s))
    {
      Vhat[k1,k2] <- 2*sum(colMeans(L1hat[,,k1]*L1hat[,,k2])*w)
    }
  }

  results <- list(value=cv,score=Shati,information=Vhat)
  return(results)
}

# Multivariate baseline proportional hazards model

MBCoxph <- function(Yi,di,Zi,Xi,K,h,beta_n)
{
  n_s <- length(Yi)
  pz_s <- dim(Zi)[2]
  px_s <- dim(Xi)[2]

  jumps <- sort(unique(Yi))
  k_s <- length(jumps)
  w_s <- colMeans(outer(Yi,jumps,FUN="=="))

  dNi <- outer(Yi,jumps,FUN="==")*di
  SYi <- outer(Yi,jumps,FUN=">=")
  eIi <- as.vector(exp(Zi %*% beta_n))

  KXij <- K.matrix(Xi,Xi,K,h)

  dNXi <- KXij %*% dNi

  S0i <- KXij %*% (SYi*eIi)

  S1i <- array(0,c(n_s,k_s,pz_s))
  S2i <- array(0,c(n_s,k_s,pz_s,pz_s))

  for (p1 in 1:pz_s)
  {
    S1i[,,p1] <- KXij %*% (SYi*Zi[,p1]*eIi)

    for (p2 in 1:pz_s)
    {
      S2i[,,p1,p2] <- KXij %*% (SYi*Zi[,p1]*Zi[,p2]*eIi)
    }
  }

  S0i.inv <- (S0i!=0)/(S0i+(S0i==0))

  a <- as.vector(S1i)*as.vector(S0i.inv)
  dim(a) <- c(n_s,k_s,pz_s)

  b <- aperm(a,c(1,3,2))
  b <- -aperm(b-as.vector(Zi),c(1,3,2))

  dLambdai <- dNXi*S0i.inv*eIi
  dLambdai <- dLambdai*(dLambdai>0)

  c <- SYi*dLambdai

  dMi <- dNi-c

  Shat <- as.vector(b)*as.vector(dMi)
  dim(Shat) <- c(n_s,k_s,pz_s)
  Shat <- apply(aperm(Shat,c(2,1,3))*w_s,c(2,3),sum)
  #Shat <- apply(Shat,2,mean)

  V1 <- array(0,c(n_s,k_s,pz_s,pz_s))
  V2 <- array(0,c(n_s,k_s,pz_s,pz_s))

  for (p1 in 1:pz_s)
  {
    for (p2 in 1:pz_s)
    {
      V1[,,p1,p2] <- a[,,p1]*a[,,p2]-S2i[,,p1,p2]*S0i.inv
      V2[,,p1,p2] <- b[,,p1]*b[,,p2]
    }
  }

  Vhat_reduced <- -as.vector(V2)*as.vector(c)
  dim(Vhat_reduced) <- c(n_s,k_s,pz_s,pz_s)
  Vhat_reduced <- apply(aperm(Vhat_reduced,c(2,1,3,4))*w_s,c(2,3,4),sum)
  Vhat_reduced <- apply(Vhat_reduced,c(2,3),mean)

  Vhat_complete <- as.vector(V1)*as.vector(dMi)
  dim(Vhat_complete) <- c(n_s,k_s,pz_s,pz_s)
  Vhat_complete <- apply(aperm(Vhat_complete,c(2,1,3,4))*w_s,c(2,3,4),sum)
  Vhat_complete <- apply(Vhat_complete,c(2,3),mean)
  Vhat_complete <- Vhat_complete+Vhat_reduced

  results <- list(score=Shat,I.complete=Vhat_complete,I.reduced=Vhat_reduced,hazard=dLambdai,time=jumps)
  return(results)
}

######################################
###                                ###
### Optimization over Grassmaniann ###
###                                ###
######################################

ntStep <- function(y, f, fy, fyy, ltol = 1e-3) {

  n <- nrow(y)
  p <- ncol(y)

  ty <- t(y)
  pit <- diag(1, n) - y %*% ty
  tpit <- t(pit)
  tpitpit <- tpit %*% pit

  b <- -pit %*% fy

  dim(b) <- c(n*p, 1)
  A <- matrix(0, n*p, n*p)

  B <- ty %*% fy;
  B <- (B + t(B))/2 # Ensure that the Hessian is symmetric

  for (j in 1:p) {
    for (l in 1:p) {
      flj <- fyy[(j-1)*n + (1:n), (l-1)*n + (1:n)]
      a1 <- tpit %*% flj %*% pit
      a2 <- B[l,j] * tpitpit
      A[(j-1)*n + (1:n), (l-1)*n + (1:n)] <- a1 - a2
    }
  }

  l.control <- abs(eigen(A, only.values = TRUE)$values)
  lam <- min(l.control)
  Lam <- max(l.control)

  if (abs(lam)<=ltol) A <- A+diag(ltol+lam,n*p)

  if (Lam<1e10)
  {
    v <- solve(A, b)
  }else
  {
    v <- matrix(0,n,p)
    warning('Hessian not calculated')
  }

  dim(v) <- c(n, p)
  v <- pit %*% v
  return(v)

}

geod <- function(y, v, r = 1) {

  if (max(abs(t(y) %*% v)) > 1e-6) warning('v not a tangent vector')

  v <- (diag(1,nrow(y)) - y%*%t(y)) %*% v
  v.svd <- svd(v)
  svd.dim <- length(v.svd$d)
  move <- cbind(y %*% v.svd$v, v.svd$u) %*% rbind(diag(cos(v.svd$d*r),svd.dim,svd.dim),
                                                  diag(sin(v.svd$d*r),svd.dim,svd.dim)) %*% t(v.svd$v)
  return(move)
}

######################################
###                                ###
### Sufficient dimenison reduction ###
###                                ###
######################################

# gradient kernel sufficient dimension reduction

# KernelDeriv()
#
# Arguments
#  Xi:  explanatory variables (input data)
#  Yi:  response variables (teaching data)
#  SGX:  bandwidth (deviation) parameter in Gaussian kernel for X
#  SGY:  bandwidth (deviation) parameter in Gaussian kernel for Y
#  EPS:  regularization coefficient
#
# Return value(s)
#  results
#    $basis: basis matrix
#    $values: eigenvalues

KernelDeriv <- function(Xi,Yi,SGX,SGY,EPS=1e-5)
{
  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]

  I <- diag(n_s)

  sx2 <- 2*SGX*SGX
  sy2 <- 2*SGY*SGY

  # Gram matrix of X
  ab <- Xi %*% t(Xi)
  aa <- as.matrix(diag(ab))
  D <- aa[,rep(1,times=n_s)]
  xx <- pmax(D+t(D)-2*ab,matrix(0,n_s,n_s))
  Kx <- exp(-xx/sx2)

  # Gram matrix of Y
  ab <- Yi %*% t(Yi)
  aa <- as.matrix(diag(ab))
  D <- aa[,rep(1,times=n_s)]
  yy <- pmax(D+t(D)-2*ab,matrix(0,n_s,n_s))
  Ky <- exp(-yy/sy2)

  # Derivative of k(X_i, x) w.r.t. x
  Dx <- Xi[rep(1:n_s,n_s),]
  dim(Dx) <- c(n_s,n_s,p_s)
  Xij <- Dx-aperm(Dx,c(2,1,3))
  Xij <- Xij/SGX/SGX
  H <- Xij*as.vector(Kx)

  #compute  sum_i H(X_i)'*Kx^-1*Ky*Kx^-1*H(X_i)

  Kx.inv <- solve(Kx+n_s*EPS*I)
  F_n <- as.vector(Kx.inv %*% Ky %*% Kx.inv)
  Hm <- H
  dim(Hm) <- c(n_s,n_s*p_s)
  HH <- t(Hm) %*% Hm
  dim(HH) <- c(n_s,p_s,n_s,p_s)
  HHm <- aperm(HH,c(1,3,2,4))
  dim(HHm) <- c(n_s*n_s,p_s,p_s)
  R <- apply(HHm*F_n,c(2,3),mean)
  RR <- eigen(R)

  results <- list(basis=RR$vectors,values=RR$values)
  return(results)
}

# MedianDist()
#
#  Computing the median of the distances from a data matrix.
#
# Arguments
#  Xi:  data matrix
#
# Return value(s)
#  s: median of the pairwise distances \|X_i - X_j\|

MedianDist <- function(Xi)
{
  n_s <- dim(Xi)[1]
  p_s <- dim(Xi)[2]
  ab <- Xi %*% t(Xi)
  aa <- diag(ab)
  Dx <- aa[rep(1:n_s,times=n_s)]+aa[rep(1:n_s,each=n_s)]
  dim(Dx) <- c(n_s,n_s)
  Dx <- Dx-2*ab
  Dx <- as.vector(Dx-diag(diag(Dx)))
  s <- sqrt(median(Dx[Dx!=0]))
  return(s)
}

# cumSIR

cumSIR <- function(Xi,Yi)
{
  n_s <- dim(Yi)[1]
  k_s <- dim(Yi)[2]
  p_s <- dim(Xi)[2]

  Di <- matrix(1,n_s,n_s)
  for (k in 1:k_s)
  {
    Di <- Di*outer(Yi[,k],Yi[,k],FUN="<=")
  }

  Xi.cs <- t(t(Xi)-colMeans(Xi))

  my <- t(Xi.cs) %*% Di/n_s
  K.cSIR <- my %*% t(my)/n_s

  RR <- eigen(K.cSIR)
  results <- list(basis=solve(var(Xi),RR$vectors),values=RR$values)
  return(results)
}

########################################
###                                  ###
### Fantope projection and selection ###
###                                  ###
########################################

FPS <- function(S,d,lambda,rho,eps,iter.max)
{
  p <- dim(S)[1]
  Yt <- matrix(0,p,p)
  Ut <- matrix(0,p,p)

  for (iter in 1:iter.max)
  {
    pf <- eigen(Yt-Ut+S/rho)
    c_points <- sort(unique(c(pf$values,pf$values-1)),decreasing=TRUE)

    root_d <- c_points[1]
    value_d <- sum(pmin(pmax(pf$values-root_d,0),1))

    if (value_d<d)
    {
      for (k.search in 2:length(c_points))
      {
        root_t <- c_points[k.search]
        value_t <- sum(pmin(pmax(pf$values-root_t,0),1))

        if (value_t<=d)
        {
          root_d <- root_t
          value_d <- value_t
        }

        if (value_t>d)
        {
          root_d <- root_d-(root_d-root_t)*(d-value_d)/(value_t-value_d)
          value_d <- sum(pmin(pmax(pf$values-root_d,0),1))
        }

        #print(c(root_d,value_d))

        if (value_t>=d)
          break
      }
    }

    Xt <- pf$vectors %*% diag(pmin(pmax(pf$values-root_d,0),1)) %*% t(pf$vectors)
    temp <- Xt+Ut
    Yt.new <- pmax(abs(temp)-lambda/rho,0)*((temp>=0)-(temp<0))
    Ut.new <- Ut+Xt-Yt.new

    stop.criterion <- max(sum((Xt-Yt.new)^2),rho^2*sum((Yt.new-Yt)^2))

    Yt <- Yt.new
    Ut <- Ut.new

    #print(iter)

    if (stop.criterion<=(d*eps^2))
      break
  }

  return(Yt)
}

#############################
###                       ###
### functional regression ###
###                       ###
#############################

# functional cumSIR
# for dense predictor

cumSIR.Xt <- function(Xti,Yi,KN=3,t_points)
{
  n_s <- length(Yi)

  y_points <- sort(unique(Yi))
  y_l <- length(y_points)
  y_weights <- colMeans(outer(Yi,y_points,FUN="=="))
  Dyi <- outer(Yi,y_points,FUN="<=")

  Xti.c <- t(t(Xti)-colMeans(Xti))

  m.ty <- t(Xti.c) %*% Dyi/n_s
  M.st <- m.ty %*% (t(m.ty)*y_weights)
  RR <- eigen(M.st)
  Gx <- t(Xti.c) %*% Xti.c/n_s

  A <- eigen(Gx)
  Gx.inv <- A$vectors[,1:KN] %*% diag(1/A$values[1:KN]) %*% t(A$vectors[,1:KN])
  beta_hat <- Gx.inv %*% RR$vectors
  beta_hat <- t(t(beta_hat)/sqrt(colSums(beta_hat^2)))

  results <- list(matrix=M.st,basis=beta_hat,values=RR$values,points=t_points)
  return(results)
}

# functional gradient-MSDR
# for dense predictor

fgMSDR.Xt <- function(Xti,Yi,KN=3,t_points,v_min,v_max,dv)
{
  n_s <- length(Yi)
  t_l <- length(t_points)
  c_points <- seq(v_min,v_max,dv)
  c_l <- length(c_points)

  vXt <- outer(Xti,c_points,FUN="*")
  Xtv <- exp(vXt*1i)
  dim(Xtv) <- c(n_s,t_l*c_l)

  rho.y <-colMeans(Xtv*Yi)
  Rx <- t(Xtv) %*% Conj(Xtv)/n_s

  A <- eigen(Rx)
  Rx.inv <- A$vectors[,1:KN] %*% diag(1/A$values[1:KN]) %*% t(A$vectors[,1:KN])
  linear.represent <- as.vector(Rx.inv %*% as.vector(rho.y))
  coe <- outer(Conj(linear.represent),linear.represent,FUN="*")*Rx
  dim(coe) <- c(t_l,c_l,t_l,c_l)
  coe <- apply(aperm(coe,c(2,1,3,4))*c_points,c(2,3,4),sum)*dv
  Mst <- Re(apply(aperm(coe,c(3,1,2))*c_points,c(2,3),sum)*dv)
  RR <- eigen(Mst %*% t(Mst))

  results <- list(matrix=Mst,basis=RR$vectors,values=RR$values,points=t_points)
  return(results)
}

# functional gradient-SDR
# for dense predictor

fgSDR.Xt <- function(Xti,Yi,KN=3,t_points,v_min=-1,v_max=1,dv=0.01)
{
  n_s <- length(Yi)
  t_l <- length(t_points)
  c_points <- seq(v_min,v_max,dv)
  c_l <- length(c_points)

  vXt <- outer(Xti,c_points,FUN="*")
  Xtv <- exp(vXt*1i)
  dim(Xtv) <- c(n_s,t_l*c_l)

  y_points <- sort(unique(Yi))
  y_l <- length(y_points)
  y_weights <- colMeans(outer(Yi,y_points,FUN="=="))
  Dy <- outer(Yi,y_points,FUN="<=")

  rho.y <- t(Dy) %*% Xtv/n_s
  Rx <- t(Xtv) %*% Conj(Xtv)/n_s

  A <- eigen(Rx)
  Rx.inv <- A$vectors[,1:KN] %*% diag(1/A$values[1:KN]) %*% t(A$vectors[,1:KN])
  linear.represent <- Rx.inv %*% t(rho.y)

  coe <- linear.represent %*% (Conj(t(linear.represent))*y_weights)*Rx
  dim(coe) <- c(t_l,c_l,t_l,c_l)
  coe <- apply(aperm(coe,c(2,1,3,4))*c_points,c(2,3,4),sum)*dv
  Mst <- Re(apply(aperm(coe,c(3,1,2))*c_points,c(2,3),sum)*dv)
  RR <- eigen(Mst %*% t(Mst))

  results <- list(matrix=Mst,basis=RR$vectors,values=RR$values,points=t_points)
  return(results)
}

# functional cumSIR
# for sparse predictor

# Xij: predictors, c(sum(Ni),3) matrix
#      column 1: subject labels
#      column 2: time points
#      column 3: predictor values

cumSIR.Xij <- function(Xij,Ni,Yi,t_points,KN=3,K,h.mu,h.m,h.Rx)
{
  n_s <- length(Yi)
  ind <- c(0,cumsum(Ni))

  mut <- as.vector(K.mean(Xi=as.matrix(Xij[,2]),Yik=as.matrix(Xij[,3]),x=as.matrix(t_points),K,h.mu))

  y_points <- sort(unique(Yi))
  y_l <- length(y_points)
  y_weights <- colMeans(outer(Yi,y_points,FUN="=="))
  Diy <- outer(Yi,y_points,FUN="<=")

  a <- K.mean(Xi=as.matrix(Xij[,2]),Yik=Diy[rep(1:n_s,times=Ni),]*Xij[,3],x=as.matrix(t_points),K,h.m)
  b <- outer(mut,colMeans(Diy),FUN="*")
  m.ty <- a-b

  M.st <- m.ty %*% (t(m.ty)*y_weights)
  RR <- eigen(M.st)

  c.d <- K.matrix(Xi=as.matrix(Xij[,2]),x=as.matrix(t_points),K,h.Rx)
  c.n <- c.d*Xij[,3]

  Rx.n <- 0
  Rx.d <- 0

  for (i in 1:n_s)
  {
    Rx.d <- Rx.d+t(c.d[rep((ind[i]+1):ind[i+1],each=Ni[i]),]) %*% c.d[rep((ind[i]+1):ind[i+1],times=Ni[i]),]
    Rx.n <- Rx.n+t(c.n[rep((ind[i]+1):ind[i+1],each=Ni[i]),]) %*% c.n[rep((ind[i]+1):ind[i+1],times=Ni[i]),]
  }

  Rx <- Rx.n*(Rx.d!=0)/(Rx.d+(Rx.d==0))
  A <- eigen(Rx)
  Rx.inv <- A$vectors[,1:KN] %*% diag(1/A$values[1:KN]) %*% t(A$vectors[,1:KN])

  beta_hat <- Rx.inv %*% RR$vectors
  beta_hat <- t(t(beta_hat)/sqrt(colSums(beta_hat^2)))

  results <- list(matrix=M.st,basis=beta_hat,values=RR$values,points=t_points)
  return(results)
}

# functional gradient-SDR
# for sparse predictor

fgSDR.Xij <- function(Xij,Ni,Yi,KN=3,t_points,v_min=-1,v_max=1,dv=0.01,K=K2_Biweight,h.rho=NULL,h.Rx=NULL)
{
  n_s <- length(Yi)
  ind <- c(0,cumsum(Ni))
  t_l <- length(t_points)
  c_points <- seq(v_min,v_max,dv)
  c_l <- length(c_points)

  y_points <- sort(unique(Yi))
  y_l <- length(y_points)
  y_weights <- colMeans(outer(Yi,y_points,FUN="=="))
  Diy <- outer(Yi,y_points,FUN="<=")

  eX <- exp(outer(Xij[,3],c_points,FUN="*")*1i)
  rho.response <- Diy[rep(1:n,times=Ni),rep(1:y_l,times=c_l)]*eX[,rep(1:c_l,each=y_l)]
  if (is.null(h.rho))
  {
    cv.h <- function(h)
    {
      cv <- mean((rho.response-K.mean.cv(Xi=as.matrix(Xij[,2]),Yik=rho.response,K,h.rho=h))^2)
      return(cv)
    }
    h.rho <- nlminb(1, objective = cv.h)$par
    rho.y <- K.mean(Xi=as.matrix(Xij[,2]),Yik=rho.response,x=as.matrix(t_points),K,h.rho)
  }else
  {
    rho.y <- K.mean(Xi=as.matrix(Xij[,2]),Yik=rho.response,x=as.matrix(t_points),K,h.rho)
  }

  #dim(rho.y) <- c(t_l,y_l*c_l)

  if (is.null(h.Rx))
  {
    h.Rx <- h.rho
  }
  c.d <- K.matrix(Xi=as.matrix(Xij[,2]),x=as.matrix(t_points),K,h.Rx)
  c.n <- c.d[,rep(1:t_l,times=c_l)]*eX[,rep(1:c_l,each=t_l)]
  #dim(c.n) <- c(sum(Ni),t_l*c_l)

  Rx.n <- 0
  Rx.d <- 0

  for (i in 1:n_s)
  {
    Rx.d <- Rx.d+t(c.d[rep((ind[i]+1):ind[i+1],each=Ni[i]),]) %*% c.d[rep((ind[i]+1):ind[i+1],times=Ni[i]),]
    Rx.n <- Rx.n+t(c.n[rep((ind[i]+1):ind[i+1],each=Ni[i]),]) %*% c.n[rep((ind[i]+1):ind[i+1],times=Ni[i]),]
  }

  dim(Rx.n) <- c(t_l,c_l,t_l,c_l)
  Rx.n <- aperm(Rx.n,c(1,3,2,4))
  dim(Rx.d) <- NULL
  Rx <- Rx.n*(Rx.d!=0)/(Rx.d+(Rx.d==0))
  Rx <- aperm(Rx,c(1,3,2,4))
  dim(Rx) <- c(t_l*c_l,t_l*c_l)
  A <- eigen(Rx)
  Rx.inv <- A$vectors[,1:KN] %*% diag(1/A$values[1:KN]) %*% t(A$vectors[,1:KN])
  dim(rho.y) <- c(t_l,y_l,c_l)
  rho.y <- aperm(rho.y,c(1,3,2))
  dim(rho.y) <- c(t_l*c_l,y_l)
  linear.represent <- Rx.inv %*% rho.y

  coe <- linear.represent %*% (Conj(t(linear.represent))*y_weights)*Rx
  dim(coe) <- c(t_l,c_l,t_l,c_l)
  coe <- apply(aperm(coe,c(2,1,3,4))*c_points,c(2,3,4),sum)*dv
  Mst <- Re(apply(aperm(coe,c(3,1,2))*c_points,c(2,3),sum)*dv)
  RR <- eigen(Mst %*% t(Mst))

  results <- list(matrix=Mst,basis=RR$vectors,values=RR$values,points=t_points)
  return(results)
}

