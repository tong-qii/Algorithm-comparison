## functions 
ase2 <- function(A, dim) {
  require(irlba)
  
  #diag(A) <- rowSums(A) / (nrow(A)-1) # diagaug line
  
  if(nrow(A) >= 400){
    A.svd <- irlba(A, nu = dim, nv = dim)
    A.svd.values <- A.svd$d[1:dim]
    A.svd.vectors <- A.svd$v[,1:dim,drop=F]
    if(dim == 1)
      A.coords <- sqrt(A.svd.values) * A.svd.vectors
    else
      A.coords <- A.svd.vectors %*% diag(sqrt(A.svd.values))
  } else{
    A.svd <- svd(A)
    A.svd.values <- A.svd$d[1:dim]
    if(dim == 1)
      A.coords <- matrix(A.svd$v[,1,drop=F] * sqrt(A.svd$d[1]),ncol=dim)
    else
      A.coords <- A.svd$v[,1:dim] %*% diag(sqrt(A.svd$d[1:dim]))
  }
  
  return(list(X=A.coords,D=A.svd.values))
}

## Procrustes alignment
procrustes<-function(X,Y, type = "I"){
  if(type == "C"){
    X <- X/norm(X, type = "F")*sqrt(nrow(X))
    Y <- Y/norm(Y, type = "F")*sqrt(nrow(Y))
  }
  if(type == "D"){
    tX <- rowSums(X^2)
    tX[tX <= 1e-15] <- 1
    tY <- rowSums(Y^2)
    tY[tY <= 1e-15] <- 1
    X <- X/sqrt(tX)
    Y <- Y/sqrt(tY)
  }
  
  tmp <- t(X) %*% Y
  tmp.svd <- svd(tmp)
  W <- tmp.svd$u %*% t(tmp.svd$v)
  return(list(error = norm(X%*%W - Y, type = "F"), W = W))
}

# distance
euklDist <- function(A,B){
  dd<-sqrt(apply(array(apply(B,1,function(x){(x-t(A))^2}),c(ncol(A),nrow(A),nrow(B))),2:3,sum))
  return(dd)}

## Creates the Classical OMNI matrix
buildOmni <- function(Alist)
{
  require(Matrix)
  
  m <- length(Alist)
  nvec <- sapply(Alist, nrow)
  nsum <- c(0, cumsum(nvec))
  
  omni <- as.matrix(bdiag(Alist))
  for (i in 1:(m-1)) {
    irng <- (nsum[i]+1):nsum[i+1]
    for (j in (i+1):m) {
      jrng <- (nsum[j]+1):nsum[j+1]
      omni[irng,jrng] <- (as.matrix(Alist[[i]]) + as.matrix(Alist[[j]]))
    }
  }
  omni <- (omni + t(omni)) / 2
  return(omni)
}

# elbow
getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE) {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006. 
  
  #  if (is.unsorted(-d))
  
  
  if (is.matrix(dat)) {
    d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
    d <- sort(dat,decreasing=TRUE)
  }
  
  if (!is.logical(threshold))
    d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
    stop(paste("d must have elements that are larger than the threshold ",
               threshold), "!", sep="")
  
  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
    mu1 <- mean(d[1:q])
    mu2 <- mean(d[-(1:q)])              # = NaN when q = p
    sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
      (p - 1 - (q < p))
    lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
      sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
    q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
    if (is.matrix(dat)) {
      sdv <- d # apply(dat,2,sd)
      plot(sdv,type="b",xlab="dim",ylab="stdev")
      points(q,sdv[q],col=2,pch=19)
    } else {
      plot(dat, type="b")
      points(q,dat[q],col=2,pch=19)
    }
  }
  
  return(q)
}








