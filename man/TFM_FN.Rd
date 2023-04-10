\name{TFM_FN}
\alias{TFM_FN}
\title{Estimation Factor Numbers via Eigenvalue-Ratio Criterion
}
\description{
This function is to estimate factor numbers via eigenvalue-ratio criterion corresponding to initial estimation without projection, one-step projection estimation, iterative projection estimation and iterative weighted projection estimation by Huber loss.}
\usage{
TFM_FN(x, r = NULL, method = "PE", tol = 1e-04, maxiter = 100)
}
\arguments{
  \item{x}{
\eqn{T \times p_1 \times \cdots \times p_K} tensor-valued time series.}
  \item{r}{
input rank of the factor tensor.}
  \item{method}{
character string, specifying the type of the factor estimation method to be used. \describe{
 \item{\code{"IE",}}{Initial estimation, without projection.}
 \item{\code{"PE",}}{One-step projection estimation.}
 \item{\code{"iPE",}}{Iterative projection estimation.}
 \item{\code{"HUBER",}}{Iterative weighted projection estimation based on huber loss function.}
}}
  \item{tol}{
tolerance in terms of the Frobenius norm.}
  \item{maxiter}{
maximum number of iterations if error stays above \code{tol}.}
}
\details{See Barigozzi et al. (2022) and Barigozzi et al. (2023) for details.}
\value{
return a list containing the following:\describe{
\item{\code{path}}{a \eqn{K \times (\rm{niter}+1)} matrix of the estimated Tucker rank of the factor process as a path of the maximum number of iteration (\eqn{\rm{niter}}) used. The \eqn{i}-th column is the estimated rank \eqn{\hat r_1, \hat r_2, \cdots, \hat r_K} at \eqn{(i-1)}-th iteration.}
\item{\code{factor.num}}{final solution of the estimated Tucker rank of the factor process \eqn{\hat r_1, \hat r_2, \cdots, \hat r_K}.}
}
}
\references{
Barigozzi M, He Y, Li L, Trapani L. Robust Estimation of Large Factor Models for Tensor-valued Time Series. <arXiv:2206.09800> 

Barigozzi M, He Y, Li L, Trapani L. Statistical Inference for Large-dimensional Tensor Factor Model by Iterative Projection. <arXiv:2303.18163> 
}
\author{
Matteo Barigozzi, Yong He, Lingxiao Li,  Lorenzo Trapani.}

\examples{
library(rTensor)
set.seed(1234)
p <- c(12,16,20) # dimensions of tensor time series
r <- c(3,4,5)  # dimensions of factor series
A<-list()
Q<-list()
for(i in 1:3){
  A[[i]]<-matrix(rnorm(p[i]*r[i],0,1),p[i],r[i])
  Q[[i]]=eigen(A[[i]]\%*\%t(A[[i]]))$vectors
}
T<-100
F<-array(NA,c(T,r))
E<-array(NA,c(T,p))
S<-array(NA,c(T,p))
X<-array(NA,c(T,p))
for(t in 1:T){
  F[t,,,]<-array(rnorm(prod(r),0,1),r)
  E[t,,,]<-array(rnorm(prod(p),0,1),p)
  S[t,,,]<-ttl(as.tensor(F[t,,,]),A,c(1,2,3))@data
  X[t,,,]<-S[t,,,]+E[t,,,]
}
rank<-TFM_FN(X,r=NULL,method='PE')
}
