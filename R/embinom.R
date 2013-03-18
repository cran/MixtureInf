## This function computes \tilde B_{22}  under Binomial mixture.
##
## \tilde B_{22} has been standardized to the correlation matrix.
##
## alpha:  mixing proportions.
## theta:  mixing parameters.
## size :  size parameter for Binomial mixture.
tb2.binom<-function(alpha, theta, size) {

  m0 <- length(alpha)
  x <- 0:size

  dbinoms <- sapply(theta, dbinom, x=x, size=size)

  y <- (matrix(rep(x, times=m0), ncol=m0) /
        rep(theta, each=length(x)) -
        (matrix(rep(size - x, times=m0), ncol=m0)) /
          rep(1-theta, each=length(x))) * dbinoms

  z <- (matrix(rep(x * (x-1), times=m0), ncol=m0) /
        rep(theta^2, each=length(x)) -
        (matrix(rep(2 * x * (size - x), times=m0), ncol=m0)) /
        rep(theta * (1-theta), each=length(x)) +
        matrix(rep((size-x) * (size-x-1), times=m0), ncol=m0) /
        rep((1-theta)^2, each=length(x))) * dbinoms / 2

  delta <- dbinoms - dbinoms[,m0]
  pdf <- as.vector(dbinoms %*% alpha)

  bi <- cbind(delta[,-m0], y, z) / pdf

  B <- crossprod(bi, diag(pdf)) %*% bi

  B11 <- B[1:(2*m0-1), 1:(2*m0-1)]
  B12 <- B[1:(2*m0-1), (2*m0):(3*m0-1)]
  B22 <- B[(2*m0):(3*m0-1), (2*m0):(3*m0-1)]

  tB22 <- B22 - crossprod(B12, solve(B11)) %*% B12
  diagB22 <- diag(diag(tB22)^(-1/2), m0, m0)

  corr <- diagB22 %*% tB22 %*% diagB22

  round(corr, 3)
}


## This function computes a_h values in Theorem 3 of Li and Chen (2009) 
##
## tb:  covariance matrix \tilde B_{22}.
## N:   number of repetitions.
## tol: a value below which is judged as 0.
emtest.thm3<-function(tb, N=10000, tol=1e-8) {

  ## Need to first install the R package quadprog.  
  
  if (!require("quadprog")) {
      stop("Package 'quadprog' is required for testing order > 3.")
  }
  m0 <- nrow(tb)
  
  eig.tb <- eigen(tb)
  tb2 <- eig.tb$vectors %*% diag(sqrt(eig.tb$values)) %*% t(eig.tb$vectors)

  output <- c()
  for (i in 1:N) {
    out <- quadprog::solve.QP(Dmat=tb, dvec=tb2%*%rnorm(m0, 0, 1), Amat=diag(rep(1, m0)), bvec=rep(0, m0))
    output <- c(output, sum(out$solution>tol))
  }
  ah <- output
  table(ah)/N
}


## This function computes the MLE of mixing distribution under null hypothesis for Binomial mixture.
##
## x:      data; can be either a vector or a matrix with the 1st column being the observed values 
##        and the 2nd column being the corresponding frequencies. 
## m0:     order under null hypothesis.
## size :  size parameter for Binomial mixture.
phi0.binom<-function(x, m0, size) {

  x <- as.matrix(x)
  if (ncol(x) == 1) {
    y <- as.matrix(table(x))
    count <- as.numeric(rownames(y))
    freq <- y[, 1]
    x <- cbind(count, freq)
  }

  count <- as.numeric(x[, 1])
  freq <- as.numeric(x[, 2])

  if (m0>1) {

    len <- 16  # number of initial values chosen for the EM-algorithm.
    eps <- 1e-5       # tolerance value for the EM-algorithm to stop. 
    output <- matrix(nrow=len, ncol=2*m0+1)
	
    for (i in 1:len) {
      alpha <- runif(m0, 0, 1)
      alpha <- alpha / sum(alpha)
      theta <- sort(runif(m0, 0, 1))

      pdf.component <- sapply(theta, dbinom, x=count, size=size) *
        rep(alpha, each=length(count)) + 1e-100/m0
      pdf <- rowSums(pdf.component)
      ln0 <- sum(freq*log(pdf))

      err <- 1
      while (err > eps) {
        w <- pdf.component / pdf
        alpha <- colSums(freq*w) / sum(freq)
        theta <- colSums(freq*w*count) / colSums(freq*w) / size

        pdf.component <- sapply(theta, dbinom, x=count, size=size) *
          rep(alpha, each=length(count)) + 1e-100/m0
        pdf <- rowSums(pdf.component)

        ln1 <- sum(freq*log(pdf))
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
  else {
    theta0 <- sum(freq*count)/sum(freq)/size
    list(alpha=1,
         theta=theta0,
         loglik=sum(freq*dbinom(count, size, theta0, log=TRUE)))
  }
}


## This function computes the EM-test statistics for Binomial mixture. 
##
## x:        data; can be either a vector or a matrix with the 1st column being the observed values 
##          and the 2nd column being the corresponding frequency. 
## outnull:  output from phi0.binom function. 
## size :    size parameter for Binomial mixture.
## C:        tuning parameter for EM-test procedure. 
## len:      number of initial values chosen for the EM-algorithm.
## eps:      tolence value for the EM-algorithm to stop. 
emstat.binom<-function(x, outnull, size, C, len=5, eps=1e-5) {

  x <- as.matrix(x)
  if (ncol(x) < 2) {
    y <- as.matrix(table(x))
    count <- as.numeric(rownames(y))
    freq <- y[,1]
    x <- cbind(count, freq)
  }	

  count <- as.numeric(x[,1])
  freq <- as.numeric(x[,2])

  theta0 <- outnull$theta
  m0 <- length(theta0)	

  ## Calculate eta_h's
  eta <- rep(0, m0+1)
  eta[1] <- 0
  eta[m0+1] <- 1
  if (m0>1) {
    for (i in 2:m0) {
      eta[i] <- (theta0[i-1]+theta0[i])/2
    }
  }

  ## Calculate the collection of beta	
  bbeta <- c()

  for (h in 1:m0) {
    bbeta <- rbind(cbind(bbeta, rep(0.1, 3^{h-1} ) ),
                   cbind(bbeta, rep(0.3, 3^{h-1} ) ),
                   cbind(bbeta, rep(0.5, 3^{h-1} ) ) )
  }

  mnk <- matrix(nrow=3^m0, ncol=3)

  ## For each beta, calculate the statistic m_n^{(k)}
  for (j in 1:(3^m0)) {

    output <- matrix(nrow=len, ncol=4*m0+1)
    beta <- bbeta[j,]

    for (i in 1:len) {
      alpha <- runif(m0, 0, 1)
      alpha <- alpha/sum(alpha)

      alpha1 <- alpha*beta
      alpha2 <- alpha*(1-beta)

      theta1 <- rep(0, m0)
      theta2 <- rep(0, m0)

      for (l in 1:m0) {
        theta1[l] <- runif(1, eta[l], eta[l+1])
        theta2[l] <- runif(1, eta[l], eta[l+1])
      }

      pdf.part1 <- sapply(theta1, dbinom, x=count, size=size)
      pdf.part2 <- sapply(theta2, dbinom, x=count, size=size)

      pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
      pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

      pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

      ln0 <- sum(freq*log(pdf))

      err <- 1

      while (err > eps) {

        w1 <- pdf.component1/pdf
        w2 <- pdf.component2/pdf

        alpha <- colSums(freq*w1 + freq*w2) / sum(freq)
        alpha1 <- alpha*beta
        alpha2 <- alpha*(1-beta)

        theta1 <- colSums(freq*w1*count) / colSums(freq*w1) /size
        theta2 <- colSums(freq*w2*count) / colSums(freq*w2) / size

        for (l in 1:m0) {
          theta1[l] <- max(min(theta1[l], eta[l+1]), eta[l])
          theta2[l] <- max(min(theta2[l], eta[l+1]), eta[l])
        }

        pdf.part1 <- sapply(theta1, dbinom, x=count, size=size)
        pdf.part2 <- sapply(theta2, dbinom, x=count, size=size)

        pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
        pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

        pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

        ln1 <- sum(freq*log(pdf))

        err <- abs(ln1-ln0)
        ln0 <- ln1
      }
      output[i,] <- c(alpha, beta, theta1, theta2, ln1)
    }

    index <- which.max(output[, 4*m0+1])

    ## The following are parameters in phi(beta_0)^{(1)}
    alpha <- output[index, 1:m0]
    beta <- output[index, (m0+1):(2*m0)]
    theta1 <- output[index, (2*m0+1):(3*m0)]
    theta2 <- output[index, (3*m0+1):(4*m0)]


    ## Iteration starts from here###
    outstat <- rep(0, 3)

    for (k in 1:3) {
      alpha1 <- alpha*beta
      alpha2 <- alpha*(1-beta)

      pdf.part1 <- sapply(theta1, dbinom, x=count, size=size)
      pdf.part2 <- sapply(theta2, dbinom, x=count, size=size)

      pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
      pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

      pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

      ln <- sum(freq*log(pdf))

      pen <- C*sum(log(1-abs(1-2*beta)))

      outstat[k] <- pen+ln	

      if (k < 3) {
        w1 <- pdf.component1/pdf
        w2 <- pdf.component2/pdf

        for (h in 1:m0) {
          if (sum(freq*w1[,h])/(sum(freq*w1[,h])+sum(freq*w2[,h])<=0.5)) {
            beta[h] <- min((sum(freq*w1[,h])+C)/(sum(freq*w1[,h])+sum(freq*w2[,h])+C), 0.5)
          }
          else {
            beta[h] <- max((sum(freq*w1[,h]))/(sum(freq*w1[,h])+sum(freq*w2[,h])+C), 0.5)
          }
        }

        alpha <- colSums(freq*w1 + freq*w2) / sum(freq)
        alpha1 <- alpha*beta
        alpha2 <- alpha*(1-beta)

        theta1 <- colSums(freq*w1*count) / colSums(freq*w1) /size
        theta2 <- colSums(freq*w2*count) / colSums(freq*w2) / size
      }
    }
    mnk[j,] <- outstat
  }

  emnk <- apply(2*(mnk-outnull$loglik), 2, max)

  ifelse(emnk < 0, 0, emnk)
}


## This function provides (1) the MLE of mixing distribution under H0,
##                       (2) the EM-test statistics,
##                       (3) the corresponding p-value,
## for testing H_0:m=m0 versus H_A:m>m_0 under Binomial mixture. 
##
## x:        data; can be either a vector or a matrix with the 1st column being the observed values 
##          and the 2nd column being the corresponding frequency. 
## m0:       order under null hypothesis.
## size :    size parameter for Binomial mixture.
## C:        tuning parameter for EM-test procedure; 
##          if not provided, the recommended value will be used.   
## len:      number of initial values chosen for the EM-algorithm.
## eps:      tolence value for the EM-algorithm to stop.  
emtest.binom<-function(x, m0, size, C=NULL, len=5, eps=1e-5) {

  x <- as.matrix(x)
  if (ncol(x) == 1) {
    y <- as.matrix(table(x))
    count <- as.numeric(rownames(y))
    freq <- y[,1]
    x <- cbind(count, freq)
  }

  freq <- as.numeric(x[,2])

  n <- sum(freq)	
  
  outnull <- phi0.binom(x, m0=m0, size=size)

  if (m0==1) {

    if (is.null(C)) C <- 0.54

    theta0 <- outnull$theta
    a1 <- 0.5
    ah <- c(1-a1, a1)
  }

  if (m0==2) {

    tb2 <- tb2.binom(outnull$alpha, outnull$theta, size)

    ## Choice of C##

    if (is.null(C)) C <- 0.5*exp(5-10.6*tb2[1, 2]-123/n)/(1+exp(5-10.6*tb2[1, 2]-123/n) )

    ah <- c(0.5-acos(tb2[1, 2])/2/pi, 0.5, acos(tb2[1, 2])/2/pi)

  }

  if (m0==3) {
    tb2 <- tb2.binom(outnull$alpha, outnull$theta, size)
	
    ## Choice of C##

    if (m0==3) C <- 0.5*exp(3.3-5.5*tb2[1, 2]-5.5*tb2[2, 3]-165/n)/(1+exp(3.3-5.5*tb2[1, 2]-5.5*tb2[2, 3]-165/n)) 

    a0 <- 0.5-acos(tb2[1, 2])/4/pi-acos(tb2[1, 3])/4/pi-acos(tb2[2, 3])/4/pi
    a2 <- 0.5-a0

    w123 <- (tb2[1, 2]-tb2[1, 3]*tb2[2, 3])/sqrt(1-tb2[1, 3]^2)/sqrt(1-tb2[2, 3]^2)
    w132 <- (tb2[1, 3]-tb2[1, 2]*tb2[3, 2])/sqrt(1-tb2[1, 2]^2)/sqrt(1-tb2[3, 2]^2)
    w231 <- (tb2[2, 3]-tb2[2, 1]*tb2[3, 1])/sqrt(1-tb2[2, 1]^2)/sqrt(1-tb2[3, 1]^2)

    a1 <- 0.75-acos(w123)/4/pi-acos(w132)/4/pi-acos(w231)/4/pi

    a3 <- 0.5-a1

    ah <- c(a0, a1, a2, a3)
  }

  if (m0>=4) {
    tb2 <- tb2.binom(outnull$alpha, outnull$theta, size)

    if (is.null(C)) {
      if (m0>=4) C <- 0.5 
    }

    ah <- emtest.thm3(tb2, N=10000, tol=1e-8)
  }


  emnk <- emstat.binom(x, outnull, size, C, len=len, eps=eps)
  pvalue <- colSums(ah * sapply(emnk, pchisq, df=0:m0, lower.tail=FALSE))

  output <- list(family='binomial',
                 m0=m0,
                 alpha=outnull$alpha,
                 theta=outnull$theta,
                 C=C,
                 emstat=emnk,
                 ah=ah,
                 pvalue=pvalue)
  class(output) <- 'emtest'

  output
}




