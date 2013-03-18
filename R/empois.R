## This function computes \tilde B_{22} under Poisson mixture.
## \tilde B_{22} has been standardized to the correlation matrix.
## 
## alpha:  mixing proportion vector.
## theta:  mixing parameter vector.
tb2.pois <- function(alpha, theta) {
  
  m0 <- length(alpha)

  size <- min(qpois(1-(1e-100), max(theta)), 200)

  x <- 0:size


  dpoises <- sapply(theta, dpois, x=x)

  y <- (matrix(rep(x, times=m0), ncol=m0) /
        rep(theta, each=length(x)) - 1) * dpoises

  z <- (matrix(rep(x * (x-1), times=m0), ncol=m0) /
        rep(theta^2, each=length(x)) -
        (matrix(rep(2 * x, times=m0), ncol=m0)) /
        rep(theta, each=length(x)) + 1) * dpoises / 2

  delta <- dpoises - dpoises[,m0]
  pdf <- as.vector(dpoises %*% alpha)

  bi <- cbind(delta[,1:(m0-1)], y, z)/pdf

  B <- crossprod(bi, diag(pdf))%*%bi

  B11 <- B[1:(2*m0-1), 1:(2*m0-1)]
  B12 <- B[1:(2*m0-1),(2*m0):(3*m0-1)]
  B22 <- B[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]

  tB22 <- B22 - crossprod(B12, solve(B11)) %*% B12
  diagB22 <- diag(diag(tB22)^(-1/2), m0, m0)

  corr <- diagB22%*%tB22%*%diagB22

  round(corr, 3)
}


## This function computes the MLE of mixing distribution under null hypothesis for Poisson mixture. 
## 
## x:   data; can be either a vector or a matrix with the 1st column being the observed values 
##      and the 2nd column being the corresponding frequency. 
## m0:  order under null hypothesis.
phi0.pois <- function(x, m0) {

  x <- as.matrix(x)
  if (dim(x)[2] == 1) {
    y <- as.matrix(table(x))
    count <- as.numeric(rownames(y))
    freq <- y[,1]
    x <- cbind(count, freq)
  }


  count <- as.numeric(x[,1])
  freq <- as.numeric(x[,2])

  if (m0>1) {
    len <- 25  # number of initial values chosen for the EM-algorithm.
    eps <- 1e-6         # tolerance value for the EM-algorithm to stop. 

    output <- matrix(nrow=len, ncol=2*m0+1)

    for (i in 1:len) {
      alpha <- runif(m0, 0, 1)
      alpha <- alpha/sum(alpha)
      theta <- sort(runif( m0, 0, max(count)))

      pdf.component <- sapply(theta, dpois, x=count) *
        rep(alpha, each=length(count)) + 1e-100/m0
      pdf <- rowSums(pdf.component)
      ln0 <- sum(freq*log(pdf))

      err <- 1
      while(err>eps) {
        w <- pdf.component/pdf
        alpha <- colSums(freq*w) / sum(freq)
        theta <- colSums(freq*w*count) / colSums(freq*w)

        pdf.component <- sapply(theta, dpois, x=count) *
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
    theta0 <- sum(count*freq)/sum(freq)
    list(alpha=1,
         theta=theta0,
         loglik=sum(freq*dpois(count, theta0, log=TRUE)) )
  }

}


## This function computes the EM-test statistics for Poisson mixture. 
## 
## x:        data; can be either a vector or a matrix with the 1st column being the observed values 
##           and the 2nd column being the corresponding frequency. 
## outnull:  output from phi0.pois function. 
## C:        tuning parameter for EM-test procedure
## len:      number of initial values for the EM-algorithms.
## eps:      tolence value for the EM-algorithm to stop.  
emstat.pois <- function(x, outnull, C, len=5, eps=1e-5) {

  x <- as.matrix(x)
  if (dim(x)[2] == 1) {
    y <- as.matrix(table(x))
    count <- as.numeric(rownames(y))
    freq <- y[,1]
    x <- cbind(count, freq)
  }

  count <- as.numeric(x[,1])
  freq <- as.numeric(x[,2])

  theta0 <- outnull$theta
  m0 <- length(theta0)	


  eta <- rep(0, m0+1)
  eta[1]=0
  eta[m0+1]=max(count)

  if (m0>1) {
    for (i in 2:m0) {
      eta[i]=(theta0[i-1]+theta0[i])/2
    }
  }

  
  bbeta <- c()

  for (h in 1:m0) {
    bbeta <- rbind(cbind(bbeta, rep(0.1, 3^{h-1} ) ),
                   cbind(bbeta, rep(0.3, 3^{h-1} ) ),
                   cbind(bbeta, rep(0.5, 3^{h-1} ) ) )
  }

  mnk <- matrix(nrow=3^m0, ncol=3)

  ## For each beta, calculate the statistic m_n^{(k)}
  for (j in 1:(3^m0)) {

    beta <- bbeta[j,]
    output <- matrix(nrow=len, ncol=4*m0+1)

    for (i in 1:len) {
      alpha <- runif(m0, 0, 1)
      alpha <- alpha/sum(alpha)

      alpha1 <- alpha*beta
      alpha2 <- alpha*(1-beta)

      theta1 <- rep(0, m0)
      theta2 <- rep(0, m0)

      for (l in 1:m0) {
        theta1[l]=runif(1, eta[l], eta[l+1])
        theta2[l]=runif(1, eta[l], eta[l+1])
      }

      pdf.part1 <- sapply(theta1, dpois, x=count)
      pdf.part2 <- sapply(theta2, dpois, x=count)

      pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
      pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

      pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

      ln0 <- sum(freq*log(pdf))

      err <- 1

      while(err>eps) {

        w1 <- pdf.component1/pdf
        w2 <- pdf.component2/pdf

        alpha <- colSums(freq*w1+freq*w2) / sum(freq)
        alpha1 <- alpha*beta
        alpha2 <- alpha*(1-beta)

        theta1 <- colSums(freq*w1*count) / colSums(freq*w1)
        theta2 <- colSums(freq*w2*count) / colSums(freq*w2)

        for (l in 1:m0) {
          theta1[l]=max(min(theta1[l], eta[l+1]), eta[l])
          theta2[l]=max(min(theta2[l], eta[l+1]), eta[l])
        }

        pdf.part1 <- sapply(theta1, dpois, x=count)
        pdf.part2 <- sapply(theta2, dpois, x=count)

        pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
        pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

        pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

        ln1 <- sum(freq*log(pdf))

        err <- abs(ln1-ln0)
        ln0 <- ln1
      }
      output[i,] <- c(alpha, beta, theta1, theta2, ln1)
    }

    index <- which.max(output[,4*m0+1])

    alpha <- output[index, 1:m0]
    beta <- output[index,(m0+1):(2*m0)]
    theta1 <- output[index,(2*m0+1):(3*m0)]
    theta2 <- output[index,(3*m0+1):(4*m0)]

    outstat <- rep(0, 3)

    for (k in 1:3) {
      alpha1 <- alpha*beta
      alpha2 <- alpha*(1-beta)

      pdf.part1 <- sapply(theta1, dpois, x=count)
      pdf.part2 <- sapply(theta2, dpois, x=count)

      pdf.component1 <- pdf.part1 * rep(alpha1, each=nrow(pdf.part1)) + 1e-100/m0
      pdf.component2 <- pdf.part2 * rep(alpha2, each=nrow(pdf.part2)) + 1e-100/m0

      pdf <- rowSums(pdf.component1) + rowSums(pdf.component2)

      ln <- sum(freq*log(pdf))

      pen <- C*sum(log(1-abs(1-2*beta)))

      outstat[k]=pen+ln	

      if (k<3) {
        w1 <- pdf.component1/pdf
        w2 <- pdf.component2/pdf

        for (h in 1:m0) {
          if (sum(freq*w1[,h])/(sum(freq*w1[,h])+sum(freq*w2[,h])<=0.5)) {
            beta[h]=min((sum(freq*w1[,h])+C)/(sum(freq*w1[,h])+sum(freq*w2[,h])+C), 0.5)
          }
          else {
            beta[h]=max((sum(freq*w1[,h]))/(sum(freq*w1[,h])+sum(freq*w2[,h])+C), 0.5)
          }
        }

        alpha <- colSums(freq*w1+freq*w2) / sum(freq)
        alpha1 <- alpha*beta
        alpha2 <- alpha*(1-beta)

        theta1 <- colSums(freq*w1*count) / colSums(freq*w1)
        theta2 <- colSums(freq*w2*count) / colSums(freq*w2)
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
##  for testing H_0:m=m0 versus H_A:m>m_0 under Poisson mixture. 
## 
## x:        data; can be either a vector or a matrix with the 1st column being the observed values 
##           and the 2nd column being the corresponding frequency. 
## m0:       order under null hypothesis.
## C:        tuning parameter for EM-test procedure; 
##           if not provided, the recommended value will be used.   
## len:      number of initial values for the EM-algorithms.
## eps:      tolence value for the EM-algorithm to stop. 
emtest.pois <- function(x, m0, C=NULL, len=5, eps=1e-5) {

  x <- as.matrix(x)
  if (dim(x)[2] == 1) {
    y <- as.matrix(table(x))
    count <- as.numeric(rownames(y))
    freq <- y[,1]
    x <- cbind(count, freq)
  }

  freq <- as.numeric(x[,2])

  n <- sum(freq)	

  outnull <- phi0.pois(x, m0)

  if (m0 == 1) {
    if (is.null(C)) C <- 0.54
    a1 <- 0.5
    ah <- c(1-a1, a1)
  }

  if (m0 == 2) {

    tb2 <- tb2.pois(outnull$alpha, outnull$theta)

    if (is.null(C)) C <- 0.5*exp(5-10.6*tb2[1, 2]-123/n)/(1+exp(5-10.6*tb2[1, 2]-123/n) )

    ah <- c(0.5-acos(tb2[1, 2])/2/pi, 0.5, acos(tb2[1, 2])/2/pi)
  }

  if (m0 == 3) {
    tb2 <- tb2.pois(outnull$alpha, outnull$theta)
	
    if (m0 == 3) C <- 0.5*exp(3.3-5.5*tb2[1, 2]-5.5*tb2[2, 3]-165/n)/(1+exp(3.3-5.5*tb2[1, 2]-5.5*tb2[2, 3]-165/n)) 

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
    tb2 <- tb2.pois(outnull$alpha, outnull$theta)

    if (is.null(C)) {
      if (m0 >= 4) C <- 0.5 
    }

    ah <- emtest.thm3(tb2, N=10000, tol=1e-8)
  }

  emnk <- emstat.pois(x, outnull, C, len=len, eps=eps)
  p1 <- sum(ah*pchisq(emnk[1], 0:m0, lower.tail = F))
  p2 <- sum(ah*pchisq(emnk[2], 0:m0, lower.tail = F))
  p3 <- sum(ah*pchisq(emnk[3], 0:m0, lower.tail = F))

  pvalue <- c(p1, p2, p3)
  
  output <- list(family='poisson',
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
