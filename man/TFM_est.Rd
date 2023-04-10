\name{TFM_est}
\alias{TFM_est}
\title{Estimation of Factor Model for High-Dimensional Tensor Time Series}
\description{This function is to estimate the tensor factor model via four different methods, namely the initial estimation without initial (IE), one-step projection estimation (PE), iterative projection estimation (iPE) and iterative weighted projection estimation by Huber loss (HUBER).
}
\usage{
TFM_est(x, r, method = "PE", tol = 1e-04, maxiter = 100)
}

\arguments{
  \item{x}{\eqn{T \times p_1 \times \cdots \times p_K} tensor-valued time series.}
  \item{r}{input rank of the factor tensor.}
  \item{method}{character string, specifying the type of the estimation method to be used. \describe{
 \item{\code{"IE",}}{Initial estimation, without projection.}
 \item{\code{"PE",}}{One-step projection estimation.}
 \item{\code{"iPE",}}{Iterative projection estimation.}
 \item{\code{"HUBER",}}{Iterative weighted projection estimation based on huber loss function.}
}}
  \item{tol}{tolerance in terms of the Frobenius norm.}
  \item{maxiter}{maximum number of iterations if error stays above \code{tol}.}
}
\details{See Barigozzi et al. (2022) and Barigozzi et al. (2023) for details.}
\value{return a list containing the following:\describe{
\item{\code{Ft}}{estimated factor processes of dimension \eqn{T \times r_1 \times r_2 \times \cdots \times r_K}.}
\item{\code{Ft.all}}{Summation of factor processes over time, of dimension \eqn{r_1,r_2,\cdots,r_K}.}
\item{\code{Q}}{a list of estimated factor loading matrices \eqn{Q_1,Q_2,\cdots,Q_K}. }
\item{\code{x.hat}}{fitted signal tensor, of dimension \eqn{T \times p_1 \times p_2 \times \cdots \times p_K}.}
\item{\code{niter}}{number of iterations.}
\item{\code{fnorm.resid}}{Frobenius norm of residuals, divide the Frobenius norm of the original tensor.}
}
}
\references{
Barigozzi M, He Y, Li L, Trapani L. Robust Estimation of Large Factor Models for Tensor-valued Time Series. <arXiv:2206.09800> 

Barigozzi M, He Y, Li L, Trapani L. Statistical Inference for Large-dimensional Tensor Factor Model by Iterative Projection. <arXiv:2303.18163>
}
\author{Matteo Barigozzi, Yong He, Lingxiao Li,  Lorenzo Trapani.}
\examples{
library(rTensor)
set.seed(1234)
p <- c(12,16,20)  # dimensions of tensor time series
r <- c(3,4,5)  # dimensions of factor series
A<-list()
Q<-list()
for(i in 1:3){
  A[[i]]<-matrix(rnorm(p[i]*r[i],0,1),p[i],r[i])
  Q[[i]]<-eigen(A[[i]]\%*\%t(A[[i]]))$vectors
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
result <- TFM_est(X,r,method='PE')  
Q.hat<-result$Q
Ft.hat <- result$Ft
}

