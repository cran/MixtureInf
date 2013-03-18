## This function computes \tilde B_{22} under Exponential mixture.
## \tilde B_{22} has been standardized to the correlation matrix.
## 
## To calculate \tilde B_{22}, the interval (0,\infty) has been divided into N=10000 subintervals 
##  with the probability in each subinterval being 1/N.
## 
## In practice, you may change N to be a large value to see its effect.
## In Li and Chen (2009), N=400000.
## 
## alpha:  mixing proportions.
## theta:  mixing parameters.
tb2.exp <- function(alpha, theta) {

  m0 <- length(alpha)
  theta <- theta/sum(alpha*theta)

  N <- 10000
  quan <- matrix((0:(N-1)+0.5)/N, ncol=1)
  x <- as.numeric(apply(quan, 1, qexpmix, alpha=alpha, theta=theta))

  dexps <- sapply(1/theta, dexp, x=x)

  y <- (matrix(rep(x, times=m0), ncol=m0) -
        rep(theta, each=length(x))) /
          rep(theta^2, each=length(x)) * dexps
  
  z <- (matrix(rep(x^2, times=m0), ncol=m0) -
        (rep(4 * theta, each=length(x)) *
         matrix(rep(x, times=m0), ncol=m0)) +
        rep(2*theta^2, each=length(x))) /
          rep(theta^4, each=length(x)) * dexps / 2
  
  delta <- dexps - dexps[,m0]
  
  pdf <- as.vector(dexps %*% alpha)
  
  bi <- cbind(delta[,-m0], y, z) / pdf

  B <- crossprod(bi)

  B11 <- B[1:(2*m0-1), 1:(2*m0-1)]
  B12 <- B[1:(2*m0-1), (2*m0):(3*m0-1)]
  B22 <- B[(2*m0):(3*m0-1), (2*m0):(3*m0-1)]

  tB22 <- B22 - crossprod(B12, solve(B11)) %*% B12
  diagB22 <- diag(diag(tB22)^(-1/2), m0, m0)

  corr <- diagB22 %*% tB22 %*% diagB22

  round(corr, 3)
}


## This function computes the MLE of  mixing distribution under null hypothesis for Exponential mixture. 
## 
## x:   data. 
## m0:  order under null hypothesis.
phi0.exp <- function(x, m0) {	
  
  if (m0>1) {

    len <- 16  # number of initial values chosen for the EM-algorithm.
    eps <- 1e-5        # tolerance value for the EM-algorithm to stop.

    output <- matrix(nrow=len, ncol=2*m0+1)

    for (i in 1:len) {
      alpha <- runif(m0, 0, 1)
      alpha <- alpha/sum(alpha)
      theta <- sort(runif( m0, min(x), max(x) ))

      pdf.component <- sapply(1/theta, dexp, x=x) *
        rep(alpha, each=length(x)) + 1e-100/m0
      pdf <- rowSums(pdf.component)
      ln0 <- sum(log(pdf))

      err <- 1
      while (err>eps) {
        w <- pdf.component/pdf
        alpha <- colMeans(w)
        theta <- colSums(w*x) / colSums(w)

        pdf.component <- sapply(1/theta, dexp, x=x) *
          rep(alpha, each=length(x)) + 1e-100/m0
        pdf <- rowSums(pdf.component)
        ln1 <- sum(log(pdf))
        err <- abs(ln1-ln0)
        ln0 <- ln1
      }
      output[i,] <- c(alpha, theta, ln1)
    }

    index <- which.max(output[,2*m0+1])
    out <- output[index,]
    alpha <- out[1:m0]
    theta <- out[(m0+1):(2*m0)]
    ln <- out[2*m0+1]

    index <- order(theta)
    alpha0 <- alpha[index]
    theta0 <- theta[index]

    list(alpha=alpha0,
         theta=theta0,
         loglik=ln)
  }
  else{
    theta0 <- mean(x)
    list(alpha=1,
         theta=theta0,
         loglik=sum(dexp(x, rate=1/theta0, log=TRUE)))
  }

}


## This function computes the EM-test statistic for Exponential mixture. 
## 
## x:        data. 
## outnull:  output from phi0.exp function. 
## C:        tuning parameter for EM-test procedure
## len:      number of initial values for the EM-algorithms.
## eps:      tolence value for the EM-algorithm to stop. 
emstat.exp <- function(x, outnull, C, len=5, eps=1e-5) {

  theta0 <- outnull$theta
  m0 <- length(theta0)	

  eta <- rep(0, m0+1)
  eta[1] <- min(x)
  eta[m0+1] <- max(x)

  if (m0>1) {
    for (i in 2:m0) {
      eta[i] <- (theta0[i-1]+theta0[i])/2
    }
  }

  
  bbeta <- c()

  for (h in 1:m0) {
    bbeta <- rbind(cbind(bbeta, rep(0.1, 3^{h-1} ) ),
                   cbind(bbeta, rep(0.3, 3^{h-1} ) ),
                   cbind(bbeta, rep(0.5, 3^{h-1} ) ) )
  }

  mnk <- matrix(nrow=3^m0, ncol=3)

  for (j in 1:(3^m0)) {

    beta <- bbeta[j,]

    output <- matrix(nrow=len, ncol=4*m0+1)

    for (i in 1:len) {
      alpha <- runif(m0, 0, 1)
      alpha <- alpha / sum(alpha)

      alpha1 <- alpha*beta
      alpha2 <- alpha*(1-beta)

      theta1 <- rep(0, m0)
      theta2 <- rep(0, m0)

      for (l in 1:m0) {
        theta1[l] <- runif(1, eta[l], eta[l+1])
        theta2[l] <- runif(1, eta[l], eta[l+1])
      }

      pdf.part1 <- sapply(1/theta1, dexp, x=x)
      pdf.part2 <- sapply(1/theta2, dexp, x=x)

      pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
      pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

      pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

      ln0 <- sum(log(pdf))

      err <- 1

      while (err > eps) {

        w1 <- pdf.component1/pdf
        w2 <- pdf.component2/pdf

        alpha <- colMeans(w1+w2)
        alpha1 <- alpha*beta
        alpha2 <- alpha*(1-beta)

        theta1 <- colSums(w1*x) / colSums(w1)
        theta2 <- colSums(w2*x) / colSums(w2)

        for (l in 1:m0) {
          theta1[l] <- max(min(theta1[l], eta[l+1]), eta[l])
          theta2[l] <- max(min(theta2[l], eta[l+1]), eta[l])
        }

        pdf.part1 <- sapply(1/theta1, dexp, x=x)
        pdf.part2 <- sapply(1/theta2, dexp, x=x)

        pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
        pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

        pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

        ln1 <- sum(log(pdf))

        err <- abs(ln1-ln0)
        ln0 <- ln1
      }
      output[i,] <- c(alpha, beta, theta1, theta2, ln1)
    }

    index <- which.max(output[,4*m0+1])

    alpha <- output[index, 1:m0]
    beta <- output[index, (m0+1):(2*m0)]
    theta1 <- output[index, (2*m0+1):(3*m0)]
    theta2 <- output[index, (3*m0+1):(4*m0)]

    outstat <- rep(0, 3)

    for (k in 1:3) {
      alpha1 <- alpha*beta
      alpha2 <- alpha*(1-beta)

      pdf.part1 <- sapply(1/theta1, dexp, x=x)
      pdf.part2 <- sapply(1/theta2, dexp, x=x)

      pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
      pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

      pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

      ln <- sum(log(pdf))

      pen <- C*sum(log(1-abs(1-2*beta)))

      outstat[k] <- pen+ln	

      if (k<3) {
        w1 <- pdf.component1/pdf
        w2 <- pdf.component2/pdf

        for (h in 1:m0) {
          if (sum(w1[,h])/(sum(w1[,h])+sum(w2[,h]))<=0.5) {
            beta[h] <- min((sum(w1[,h])+C)/(sum(w1[,h])+sum(w2[,h])+C), 0.5)
          }
          else {
            beta[h] <- max((sum(w1[,h]))/(sum(w1[,h])+sum(w2[,h])+C), 0.5)
          }
        }

        alpha <- colMeans(w1+w2)
        alpha1 <- alpha*beta
        alpha2 <- alpha*(1-beta)

        theta1 <- colSums(w1*x) / colSums(w1)
        theta2 <- colSums(w2*x) / colSums(w2)
      }
    }
    mnk[j,] <- outstat
  }
  emnk <- apply(2*(mnk-outnull$loglik), 2, max)

  ifelse(emnk < 0, 0, emnk)
}


## This function provides (1) the MLE of mixing distribution under H0,
##                        (2) the EM-test statistics,
##                        (3) the corresponding p-value,  
##  for testing H_0:m=m0 versus H_A:m>m_0 under Exponential mixture. 
## 
## x:        data. 
## m0:       order under null hypothesis.
## C:        tuning parameter for EM-test procedure; 
##           if not provided, the recommended value will be used.   
## len:      number of initial values for the EM-algorithms.
## eps:      tolence value for the EM-algorithm to stop. 
emtest.exp <- function(x, m0, C=NULL, len=5, eps=1e-5) {

  meanx <- mean(x)
  x <- x/meanx
  outnull <- phi0.exp(x, m0)
  n <- length(x)

  if (m0==1) {

    ##Choice of C##
    if (is.null(C)) C <- exp(0.74+82/n)/(1+exp(0.74+82/n))

    a1 <- 0.5 
    ah <- c(1-a1, a1)
  }

  if (m0==2) {

    tb2 <- tb2.exp(outnull$alpha, outnull$theta)

    ##Choice of C##

    if (is.null(C)) C <- exp(2.3-8.5*tb2[1, 2])/(1+exp(2.3-8.5*tb2[1, 2]) )

    ah <- c(0.5-acos(tb2[1, 2])/2/pi, 0.5, acos(tb2[1, 2])/2/pi)
  }

  if (m0==3) {
    tb2 <- tb2.exp(outnull$alpha, outnull$theta)
	
    ##Choice of C##

    if (m0==3) C <- exp(2.2-5.9*tb2[1, 2]-5.9*tb2[2, 3])/(1+exp(2.3-5.9*tb2[1, 2]-5.9*tb2[2, 3])) 

    a0 <- 0.5-acos(tb2[1, 2])/4/pi-acos(tb2[1, 3])/4/pi-acos(tb2[2, 3])/4/pi
    a2 <- 0.5-a0

    w123 <- (tb2[1, 2]-tb2[1, 3]*tb2[2, 3])/sqrt(1-tb2[1, 3]^2)/sqrt(1-tb2[2, 3]^2)
    w132 <- (tb2[1, 3]-tb2[1, 2]*tb2[3, 2])/sqrt(1-tb2[1, 2]^2)/sqrt(1-tb2[3, 2]^2)
    w231 <- (tb2[2, 3]-tb2[2, 1]*tb2[3, 1])/sqrt(1-tb2[2, 1]^2)/sqrt(1-tb2[3, 1]^2)

    a1 <- 0.75-acos(w123)/4/pi-acos(w132)/4/pi-acos(w231)/4/pi

    a3 <- 0.5-a1

    ah <- c(a0, a1, a2, a3)
  }

  if (m0 >= 4) {
    tb2 <- tb2.exp(outnull$alpha, outnull$theta)

    if (is.null(C)) {
      if (m0 >= 4) C <- 0.5 
    }

    ah <- emtest.thm3(tb2, N=10000, tol=1e-8)
  }

  emnk <- emstat.exp(x, outnull, C, len=len, eps=eps)
  pvalue <- colSums(ah * sapply(emnk, pchisq, df=0:m0, lower.tail=FALSE))

  output <- list(family='exponential',
                 m0=m0,
                 alpha=outnull$alpha,
                 theta=outnull$theta*meanx,
                 C=C,
                 emstat=emnk,
                 ah=ah,
                 pvalue=pvalue)
  class(output) <- 'emtest'

  output
}


