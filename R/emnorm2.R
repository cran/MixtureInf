### The following five functions: over1, over2, over3, formula2 and formula3 are ###
### used to calculate the an values for m0=2 and 3 ###
over1 <- function(alpi, mui, sigi, alpj, muj, sigj) {

  if (sigi == sigj) {
    delta <- abs(mui-muj)/sigi

    out <- pnorm(-delta/2 + log( alpj / alpi)/delta, 0, 1)
  }

  if (sigi>sigj) {
    ncp <- (mui-muj)*sigi/(sigi^2-sigj^2)

    value <- sigj^2*(mui-muj)^2/(sigj^2-sigi^2)^2-sigj^2/(sigi^2-sigj^2)*log(alpi^2*sigj^2/alpj^2/sigi^2 )
    sqrt.value <- sqrt(max(value, 0))

    out <- pnorm(sqrt.value-ncp, 0, 1)-pnorm(-sqrt.value-ncp, 0, 1)
  }

  if (sigi<sigj) {

    ncp <- (mui-muj)*sigi/(sigi^2-sigj^2)

    value <- sigj^2*(mui-muj)^2/(sigj^2-sigi^2)^2-sigj^2/(sigi^2-sigj^2)*log(alpi^2*sigj^2/alpj^2/sigi^2 )

    sqrt.value <- sqrt(max(value, 0))


    out <- 1-pnorm(sqrt.value-ncp, 0, 1)+pnorm(-sqrt.value-ncp, 0, 1)

  }

  out
}


over2 <- function(par) {
  
  alp1 <- par$alpha[1]
  alp2 <- par$alpha[2]
  mu1 <- par$mu[1]
  mu2 <- par$mu[2]
  sig1 <- par$sigma[1]
  sig2 <- par$sigma[2]

  part1 <- over1(alp1, mu1, sig1, alp2, mu2, sig2)
  part2 <- over1(alp2, mu2, sig2, alp1, mu1, sig1)
  (part1+part2)/2
}


over3 <- function(par) {
  
  alp1 <- par$alpha[1]
  alp2 <- par$alpha[2]
  alp3 <- par$alpha[3]


  mu1 <- par$mu[1]
  mu2 <- par$mu[2]
  mu3 <- par$mu[3]


  sig1 <- par$sigma[1]
  sig2 <- par$sigma[2]
  sig3 <- par$sigma[3]


  part1 <- over1(alp1, mu1, sig1, alp2, mu2, sig2)
  part2 <- over1(alp2, mu2, sig2, alp1, mu1, sig1)
  w12 <- (part1+part2)/2

  part3 <- over1(alp2, mu2, sig2, alp3, mu3, sig3)
  part4 <- over1(alp3, mu3, sig3, alp2, mu2, sig2)
  w23 <- (part3+part4)/2
  c(w12, w23)
}


formula2 <- function(par, n) {
  overlap <- over2(par)
  tover <- log(overlap/(1-overlap))
  an <- 0.35*exp(-1.859    -0.577*tover -60.453/n  )/(1+exp(-1.859    -0.577*tover -60.453/n))
  an
}


formula3 <- function(par, n) {
  overlap <- over3(par)
  tover1 <- log(overlap[1]/(1-overlap[1]))
  tover2 <- log(overlap[2]/(1-overlap[2]))

  an2 <- -1.602-0.240*tover1-0.240*tover2-130.394/n
  an <- 0.35*exp(an2)/(1+exp(an2))
  an
}


### The penaly function on \sigma^2 ###
pn <- function(sigma, sigma0, an) {
  -an*( sigma0^2/sigma^2 + log(sigma^2/sigma0^2) -1)
}


### The following three functions are used to  ###
### calculate the EM-test statistics and their p-values for testing order=1 ###
max2.norm <- function(x, beta, eps=1e-6) {
  mu <- mean(x)
  sn <- mean((x-mu)^2)

  output <- c()

  for (i in 1:25) {
    mu1 <- runif(1, min(x), max(x))
    mu2 <- runif(1, min(x), max(x))
    sig1 <- runif(1, 0.01, 1)*sn
    sig2 <- runif(1, 0.01, 1)*sn
    alp1 <- beta
    alp2 <- 1-beta
    ln0 <- sum(log(alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2)+1e-50 ))
    pen0 <- pn(sig1, sn, 0.25)+pn(sig2, sn, 0.25)
    pln0 <- ln0+pen0
    err <- 1

    while (err>eps) {
      pdf <- alp1*dnorm(x, mu1, sig1)+ alp2*dnorm(x, mu2, sig2)+1e-50
      
      w1 <- ( alp1*dnorm(x, mu1, sig1)+1e-50/2 )/pdf
      w2 <- ( alp2*dnorm(x, mu2, sig2)+1e-50/2 )/pdf

      mu1 <- sum(w1*x)/sum(w1)
      mu2 <- sum(w2*x)/sum(w2)

      sig1 <- sqrt((sum(w1*(x-mu1)^2)+2*0.25*sn) / (sum(w1)+2*0.25))
      sig2 <- sqrt((sum(w2*(x-mu2)^2)+2*0.25*sn) / (sum(w2)+2*0.25))
      ln1 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2)+1e-50 ))
      pen1 <- pn(sig1, sn, 0.25)+pn(sig2, sn, 0.25)
      pln1 <- ln1+pen1
      err <- abs(pln1-pln0)
      pln0 <- pln1
    }
    output <- rbind(output, c(beta, mu1, mu2, sig1, sig2, pln1))
  }
  len <- nrow(output)
  index <- (1:len)[output[,5] >= max(output[,5])][1]
  outpar <- output[index,]
  outpar
}


emite2.norm <- function(x, pars) {
  sn <- mean((x-mean(x))^2)

  beta   <- pars[1]
  mu1    <- pars[2]
  mu2    <- pars[3]
  sig1   <- pars[4]
  sig2   <- pars[5]
  alp1   <- beta
  alp2   <- 1-beta

  output <- rep(0, 3)
  output[1] <- sum(log(alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2)+1e-50 ))+pn(sig1, sn, 0.25)+pn(sig2, sn, 0.25)+log(1-abs(2*beta-1))

  for ( i in 1:2) {
    pdf <- alp1*dnorm(x, mu1, sig1)+ alp2*dnorm(x, mu2, sig2)+1e-50
    
    w1 <- ( alp1*dnorm(x, mu1, sig1)+1e-50/2 )/pdf
    w2 <- ( alp2*dnorm(x, mu2, sig2)+1e-50/2 )/pdf

    mu1 <- sum(w1*x)/sum(w1)
    mu2 <- sum(w2*x)/sum(w2)

    sig1 <- sqrt((sum(w1*(x-mu1)^2)+2*0.25*sn) / (sum(w1)+2*0.25))
    sig2 <- sqrt((sum(w2*(x-mu2)^2)+2*0.25*sn) / (sum(w2)+2*0.25))

    if (sum(w1)/(sum(w1)+sum(w2))<=0.5) {
      beta <- min((sum(w1)+1)/(sum(w1)+sum(w2)+1), 0.5)
    }
    else {
      beta <- max((sum(w1))/(sum(w1)+sum(w2)+1), 0.5)
    }

    alp1   <- beta
    alp2   <- 1-beta

    output[i+1] <- sum(log(alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2)+1e-50 ))+pn(sig1, sn, 0.25)+pn(sig2, sn, 0.25)+log(1-abs(2*beta-1))
  }
  output
}


emtest1.norm <- function(x, eps=1e-6) {
  output <- c()

  for (beta in c(0.1, 0.3, 0.5)) {
    out <- max2.norm(x, beta, eps=eps)[1:5]
    output <- rbind(output, emite2.norm(x, out))
  }

  out.alt <- apply(output, 2, max)

  mu0 <- mean(x)
  sig0 <- sqrt( mean( (x-mu0)^2 ) )
  out.null <- sum(log(dnorm(x, mu0, sig0)))

  emstat <- 2*(out.alt-out.null)

  list(alpha=1,
       theta=list(mu=mu0, sigma=sig0),
       emstat=emstat,
       pvalue=pchisq(emstat, 2, lower.tail=F))
}


### The following four functions are used to  ###
### calculate the EM-test statistics and pvalues for testing order=2 ###
phi2.norm <- function(x, eps=1e-6) {

  sn <- sqrt( mean((x-mean(x))^2) )
  n <- length(x)

  output <- c()

  len <- 25  # number of initial values chosen for the EM-algorithm.
  eps <- 1e-5       # tolerance value for the EM-algorithm to stop.
  output <- matrix(nrow=len, ncol=8)
  
  for (i in 1:len) {
    alp1 <- runif(1, 0, 1)
    alp2 <- 1-alp1
    mu1 <- runif(1, min(x), max(x))
    mu2 <- runif(1, min(x), max(x))
    sig1 <- runif(1, 0.01, 1)*sn
    sig2 <- runif(1, 0.01, 1)*sn
    ln0 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2)+1e-50 ))
    pen0 <- pn(sig1, sn, 1/n)+pn(sig2, sn, 1/n)
    pln0 <- ln0+pen0

    err <- 1
    while (err>eps) {
      pdf <- alp1*dnorm(x, mu1, sig1)+ alp2*dnorm(x, mu2, sig2)+1e-50
      
      w1 <- ( alp1*dnorm(x, mu1, sig1)+1e-50/2 )/pdf
      w2 <- ( alp2*dnorm(x, mu2, sig2)+1e-50/2 )/pdf

      alp1 <- mean(w1)
      alp2 <- mean(w2)

      mu1 <- sum(w1*x)/sum(w1+1e-50/2)
      mu2 <- sum(w2*x)/sum(w2+1e-50/2)

      sig1 <- sqrt( ( sum(w1*(x-mu1)^2)+2*sn/n )/(sum(w1)+2/n ) )
      sig2 <- sqrt( ( sum(w2*(x-mu2)^2)+2*sn/n )/(sum(w2)+2/n ) )

      ln1 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2)+1e-50 ))
      pen1 <- pn(sig1, sn, 1/n)+pn(sig2, sn, 1/n)
      pln1 <- ln1+pen1

      err <- abs(pln1-pln0)
      pln0 <- pln1
    }

    output[i,] <- c(alp1, alp2, mu1, mu2, sig1, sig2, ln1, pln1)
  }

  index <- which.max(output[,8])
  out <- output[index,]
  if (out[3] > out[4]) {
    out[1:2] <- out[2:1]
    out[3:4] <- out[4:3]
    out[5:6] <- out[6:5]
  }

  list(alpha=out[1:2], mu=out[3:4], sigma=out[5:6], loglik=out[7])
}


max4.norm <- function(x, an, beta, par0, eps=1e-6) {
  mid <- mean(par0$mu)

  sig01 <- par0$sigma[1]
  sig02 <- par0$sigma[2]


  len <- 25
  output <- matrix(nrow=len, ncol=13)
  for (i in 1:len) {
    alpha <- par0$alpha
    mu1 <- runif(2, min(x), mid)
    mu2 <- runif(2, mid, max(x))
    mu <- c(mu1, mu2)
    sig1 <- runif(2, 0.25*sig01, 2*sig01)
    sig2 <- runif(2, 0.25*sig02, 2*sig02)
    sigma <- c(sig1, sig2)

    beta1 <- beta[1]
    beta2 <- beta[2]

    alp1 <- alpha[1]*beta1
    alp2 <- alpha[1]*(1-beta1)
    alp3 <- alpha[2]*beta2
    alp4 <- alpha[2]*(1-beta2)

    mu1 <- mu[1]
    mu2 <- mu[2]
    mu3 <- mu[3]
    mu4 <- mu[4]

    sig1 <- sigma[1]
    sig2 <- sigma[2]
    sig3 <- sigma[3]
    sig4 <- sigma[4]

    ln0 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4)+1e-50))

    pen0 <- pn(sig1, sig01, an)+pn(sig2, sig01, an)+pn(sig3, sig02, an)+pn(sig4, sig02, an)

    pln0 <- ln0+pen0

    err <- 1
    while (err>eps) {
      pdf <- alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) +1e-50

      w1 <- (alp1*dnorm(x, mu1, sig1)+1e-50/4)/pdf
      w2 <- (alp2*dnorm(x, mu2, sig2)+1e-50/4)/pdf
      w3 <- (alp3*dnorm(x, mu3, sig3)+1e-50/4)/pdf
      w4 <- (alp4*dnorm(x, mu4, sig4)+1e-50/4)/pdf

      alpha[1] <- mean(w1+w2)
      alpha[2] <- mean(w3+w4)

      alp1 <- alpha[1]*beta1
      alp2 <- alpha[1]*(1-beta1)
      alp3 <- alpha[2]*beta2
      alp4 <- alpha[2]*(1-beta2)

      mu1 <- sum(w1*x)/sum(w1)
      mu2 <- sum(w2*x)/sum(w2)
      mu1 <- min(mu1, mid)
      mu2 <- min(mu2, mid)

      mu3 <- sum(w3*x)/sum(w3)
      mu4 <- sum(w4*x)/sum(w4)

      mu3 <- max(mu3, mid)
      mu4 <- max(mu4, mid)

      sig1 <- sqrt( ( sum(w1*(x-mu1)^2)+2*an*sig01^2 )/(sum(w1)+2*an ) )
      sig2 <- sqrt( ( sum(w2*(x-mu2)^2)+2*an*sig01^2 )/(sum(w2)+2*an ) )
      sig3 <- sqrt( ( sum(w3*(x-mu3)^2)+2*an*sig02^2 )/(sum(w3)+2*an ) )
      sig4 <- sqrt( ( sum(w4*(x-mu4)^2)+2*an*sig02^2 )/(sum(w4)+2*an ) )

      ln0 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) +1e-50))
      pen0 <- pn(sig1, sig01, an)+pn(sig2, sig01, an)+pn(sig3, sig02, an)+pn(sig4, sig02, an)

      pln1 <- ln0+pen0
      err <- abs(pln1-pln0)
      pln0 <- pln1
    }

    output[i,] <- c(beta, alpha, mu1, mu2, mu3, mu4, sig1, sig2, sig3, sig4, pln1)
  }

  index <- which.max(output[,13])
  output[index,]
}


emite4.norm <- function(x, par, par0, an) {
  beta1 <- par[1]
  beta2 <- par[2]

  alpha <- par[3]

  alp1 <- alpha*beta1
  alp2 <- alpha*(1-beta1)
  alp3 <- (1-alpha)*beta2
  alp4 <- (1-alpha)*(1-beta2)

  mu1 <- par[5]
  mu2 <- par[6]
  mu3 <- par[7]
  mu4 <- par[8]

  sigma1 <- par[9]
  sigma2 <- par[10]
  sigma3 <- par[11]
  sigma4 <- par[12]

  sigma10 <- par0$sigma[1]
  sigma20 <- par0$sigma[2]

  ln <- sum(log(alp1*dnorm(x, mu1, sigma1)+ alp2*dnorm(x, mu2, sigma2)+alp3*dnorm(x, mu3, sigma3)+alp4*dnorm(x, mu4, sigma4)+1e-100))

  pen1 <- pn(sigma1, sigma10, an)+pn(sigma2, sigma10, an)+pn(sigma3, sigma20, an)+pn(sigma4, sigma20, an)

  pen2 <-  log(1-abs(1-2*beta1))+ log(1-abs(1-2*beta2))

  pln <- ln+pen1+pen2

  k <- 2

  output <- rep(0, k+1)
  output[1] <- pln

  for (i in 1:k) {

    pdf <- alp1*dnorm(x, mu1, sigma1)+ alp2*dnorm(x, mu2, sigma2)+alp3*dnorm(x, mu3, sigma3)+alp4*dnorm(x, mu4, sigma4)+1e-100
    
    w1 <- ( alp1*dnorm(x, mu1, sigma1)+1e-100/4 )/pdf
    w2 <- ( alp2*dnorm(x, mu2, sigma2)+1e-100/4 )/pdf
    w3 <- ( alp3*dnorm(x, mu3, sigma3)+1e-100/4 )/pdf
    w4 <- ( alp4*dnorm(x, mu4, sigma4)+1e-100/4 )/pdf

    alpha <- mean(w1+w2)

    mu1 <- sum(w1*x)/sum(w1)
    mu2 <- sum(w2*x)/sum(w2)
    mu3 <- sum(w3*x)/sum(w3)
    mu4 <- sum(w4*x)/sum(w4)

    sigma1 <- sqrt( ( sum(w1*(x-mu1)^2)+2*an*sigma10^2 )/(sum(w1)+2*an ) )
    sigma2 <- sqrt( ( sum(w2*(x-mu2)^2)+2*an*sigma10^2 )/(sum(w2)+2*an ) )
    sigma3 <- sqrt( ( sum(w3*(x-mu3)^2)+2*an*sigma20^2 )/(sum(w3)+2*an ) )
    sigma4 <- sqrt( ( sum(w4*(x-mu4)^2)+2*an*sigma20^2 )/(sum(w4)+2*an ) )

    if (sum(w1)/(sum(w1)+sum(w2))<=0.5) {
      beta1 <- min((sum(w1)+1)/(sum(w1)+sum(w2)+1), 0.5)
    }
    else {
      beta1 <- max((sum(w1))/(sum(w1)+sum(w2)+1), 0.5)
    }


    if (sum(w3)/(sum(w3)+sum(w4))<=0.5) {
      beta2 <- min( (sum(w3)+1)/(sum(w3)+sum(w4)+1), 0.5 )
    }
    else {
      beta2 <- max( (sum(w3))/(sum(w3)+sum(w4)+1), 0.5 )
    }

    alp1 <- alpha*beta1
    alp2 <- alpha*(1-beta1)
    alp3 <- (1-alpha)*beta2
    alp4 <- (1-alpha)*(1-beta2)

    ln <- sum(log(alp1*dnorm(x, mu1, sigma1)+ alp2*dnorm(x, mu2, sigma2)+alp3*dnorm(x, mu3, sigma3)+alp4*dnorm(x, mu4, sigma4)+1e-100))

    pen1 <- pn(sigma1, sigma10, an)+pn(sigma2, sigma10, an)+pn(sigma3, sigma20, an)+pn(sigma4, sigma20, an)

    pen2 <- log(1-abs(1-2*beta1))+log(1-abs(1-2*beta2))

    pln <- ln+pen1+pen2

    output[i+1] <- pln
  }
  output
}


emtest2.norm <- function(x, eps=1e-6) {
  par0 <- phi2.norm(x, eps=eps)
  ln0 <- par0$loglik
  an <- formula2(par0, length(x))

  output <- c()

  for (beta1 in c(0.1, 0.3, 0.5)) {
    for (beta2 in c(0.1, 0.3, 0.5)) {
      
      beta <- c(beta1, beta2)

      par1 <- max4.norm(x, an, beta, par0, eps=eps)

      out <- emite4.norm(x, par1, par0, an=an)

      output <- rbind(output, out)
    }
  }

  output <- 2*(output-ln0)

  emstat <- apply(output, 2, max)

  list(alpha=par0$alpha,
       theta=list(mu=par0$mu, sigma=par0$sigma),
       emstat=emstat,
       pvalue=pchisq(emstat, 4, lower.tail=F))
}


### The following four functions are used to ###
### Calculate the EM-test statistics and P-values for testing order=3 ###
phi3.norm <- function(x, eps=1e-6) {
  sn <- sqrt( mean((x-mean(x))^2) )
  n <- length(x)

  len <- 25
  output <- matrix(nrow=len, ncol=11)

  for (i in 1:len) {
    alpha <- runif(3, 0, 1)
    alpha <- alpha/sum(alpha)
    mu <- runif(3, min(x), max(x))
    sigma <- runif(3, 0.1, 1)*sn

    alp1 <- alpha[1]
    alp2 <- alpha[2]
    alp3 <- alpha[3]

    mu1 <- mu[1]
    mu2 <- mu[2]
    mu3 <- mu[3]

    sig1 <- sigma[1]
    sig2 <- sigma[2]
    sig3 <- sigma[3]

    ln0 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+1e-50))
    pen0 <- pn(sig1, sn, 1/n)+pn(sig2, sn, 1/n)+pn(sig3, sn, 1/n)
    pln0 <- ln0+pen0
    err <- 1

    while (err>eps) {
      pdf <- alp1*dnorm(x, mu1, sig1)+ alp2*dnorm(x, mu2, sig2)+alp3*dnorm(x, mu3, sig3)+1e-50
      
      w1 <- ( alp1*dnorm(x, mu1, sig1)+1e-50/3 )/pdf
      w2 <- ( alp2*dnorm(x, mu2, sig2)+1e-50/3 )/pdf
      w3 <- ( alp3*dnorm(x, mu3, sig3)+1e-50/3 )/pdf

      alp1 <- mean(w1)
      alp2 <- mean(w2)
      alp3 <- mean(w3)

      mu1 <- sum(w1*x)/sum(w1+1e-50/3)
      mu2 <- sum(w2*x)/sum(w2+1e-50/3)
      mu3 <- sum(w3*x)/sum(w3+1e-50/3)

      sig1 <- sqrt( ( sum(w1*(x-mu1)^2)+2*sn/n )/(sum(w1)+2/n ) )
      sig2 <- sqrt( ( sum(w2*(x-mu2)^2)+2*sn/n )/(sum(w2)+2/n ) )
      sig3 <- sqrt( ( sum(w3*(x-mu3)^2)+2*sn/n )/(sum(w3)+2/n ) )

      ln1 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+1e-50))
      pen1 <- pn(sig1, sn, 1/n)+pn(sig2, sn, 1/n)+pn(sig3, sn, 1/n)
      pln1 <- ln1+pen1
      err <- abs(pln1-pln0)
      pln0 <- pln1
    }
    output[i,] <- c(alp1, alp2, alp3, mu1, mu2, mu3, sig1, sig2, sig3, ln1, pln1)
  }

  index <- which.max(output[,11])
  outpar <- output[index,]

  alpha <- outpar[1:3]
  mu <- outpar[4:6]
  sigma <- outpar[7:9]

  index <- sort(mu, index.return = TRUE)$ix
  alpha <- alpha[index]
  mu <- mu[index]
  sigma <- sigma[index]

  ln <- outpar[10]

  list(alpha=alpha, mu=mu, sigma=sigma, loglik=ln)
}


max6.norm <- function(x, an, beta, par0, eps=1e-6) {
  
  mid1 <- mean(par0$mu[1:2])
  mid2 <- mean(par0$mu[2:3])

  sig01 <- par0$sigma[1]
  sig02 <- par0$sigma[2]
  sig03 <- par0$sigma[3]


  len <- 25
  output <- matrix(nrow=len, ncol=19)
  for (i in 1:len) {
    alpha <- par0$alpha
    mu1 <- runif(2, min(x), mid1)
    mu2 <- runif(2, mid1, mid2)
    mu3 <- runif(2, mid2, max(x))
    mu <- c(mu1, mu2, mu3)
    sig1 <- runif(2, 0.25*sig01, 2*sig01)
    sig2 <- runif(2, 0.25*sig02, 2*sig02)
    sig3 <- runif(2, 0.25*sig03, 2*sig03)
    sigma <- c(sig1, sig2, sig3)

    beta1 <- beta[1]
    beta2 <- beta[2]
    beta3 <- beta[3]

    alp1 <- alpha[1]*beta1
    alp2 <- alpha[1]*(1-beta1)
    alp3 <- alpha[2]*beta2
    alp4 <- alpha[2]*(1-beta2)
    alp5 <- alpha[3]*beta3
    alp6 <- alpha[3]*(1-beta3)

    mu1 <- mu[1]
    mu2 <- mu[2]
    mu3 <- mu[3]
    mu4 <- mu[4]
    mu5 <- mu[5]
    mu6 <- mu[6]

    sig1 <- sigma[1]
    sig2 <- sigma[2]
    sig3 <- sigma[3]
    sig4 <- sigma[4]
    sig5 <- sigma[5]
    sig6 <- sigma[6]

    ln0 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) + alp5*dnorm(x, mu5, sig5) + alp6*dnorm(x, mu6, sig6) +1e-50))

    pen0 <- pn(sig1, sig01, an)+pn(sig2, sig01, an)+pn(sig3, sig02, an)+pn(sig4, sig02, an)+pn(sig5, sig03, an)+pn(sig6, sig03, an)

    pln0 <- ln0+pen0
    err <- 1

    while (err>eps) {
      
      pdf <- alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) + alp5*dnorm(x, mu5, sig5) + alp6*dnorm(x, mu6, sig6) +1e-50

      w1 <- (alp1*dnorm(x, mu1, sig1)+1e-50/6)/pdf
      w2 <- (alp2*dnorm(x, mu2, sig2)+1e-50/6)/pdf
      w3 <- (alp3*dnorm(x, mu3, sig3)+1e-50/6)/pdf
      w4 <- (alp4*dnorm(x, mu4, sig4)+1e-50/6)/pdf
      w5 <- (alp5*dnorm(x, mu5, sig5)+1e-50/6)/pdf
      w6 <- (alp6*dnorm(x, mu6, sig6)+1e-50/6)/pdf

      alpha[1] <- mean(w1+w2)
      alpha[2] <- mean(w3+w4)
      alpha[3] <- mean(w5+w6)

      alp1 <- alpha[1]*beta1
      alp2 <- alpha[1]*(1-beta1)
      alp3 <- alpha[2]*beta2
      alp4 <- alpha[2]*(1-beta2)
      alp5 <- alpha[3]*beta3
      alp6 <- alpha[3]*(1-beta3)

      mu1 <- sum(w1*x)/sum(w1)
      mu2 <- sum(w2*x)/sum(w2)
      mu1 <- min(mu1, mid1)
      mu2 <- min(mu2, mid1)

      mu3 <- sum(w3*x)/sum(w3)
      mu4 <- sum(w4*x)/sum(w4)

      mu3 <- min(max(mu3, mid1), mid2)
      mu4 <- min(max(mu4, mid1), mid2)

      mu5 <- sum(w5*x)/sum(w5)
      mu6 <- sum(w6*x)/sum(w6)

      mu5 <- max(mid2, mu5)
      mu6 <- max(mid2, mu6)

      sig1 <- sqrt( ( sum(w1*(x-mu1)^2)+2*an*sig01^2 )/(sum(w1)+2*an ) )
      sig2 <- sqrt( ( sum(w2*(x-mu2)^2)+2*an*sig01^2 )/(sum(w2)+2*an ) )
      sig3 <- sqrt( ( sum(w3*(x-mu3)^2)+2*an*sig02^2 )/(sum(w3)+2*an ) )
      sig4 <- sqrt( ( sum(w4*(x-mu4)^2)+2*an*sig02^2 )/(sum(w4)+2*an ) )
      sig5 <- sqrt( ( sum(w5*(x-mu5)^2)+2*an*sig03^2 )/(sum(w5)+2*an ) )
      sig6 <- sqrt( ( sum(w6*(x-mu6)^2)+2*an*sig03^2 )/(sum(w6)+2*an ) )

      ln0 <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) + alp5*dnorm(x, mu5, sig5) + alp6*dnorm(x, mu6, sig6) +1e-50))

      pen0 <- pn(sig1, sig01, an)+pn(sig2, sig01, an)+pn(sig3, sig02, an)+pn(sig4, sig02, an)+pn(sig5, sig03, an)+pn(sig6, sig03, an)
      pln1 <- ln0+pen0
      err <- abs(pln1-pln0)
      pln0 <- pln1
    }



    output[i,] <- c(beta, alpha, mu1, mu2, mu3, mu4, mu5, mu6, sig1, sig2, sig3, sig4, sig5, sig6, pln1)
  }

  index <- which.max(output[,19])
  output[index,]
}


emite6.norm <- function(x, par, par0, an) {
  
  beta <- par[1:3]
  alpha <- par[4:6]
  mu <- par[7:12]
  sigma <- par[13:18]
  
  beta1 <- beta[1]
  beta2 <- beta[2]
  beta3 <- beta[3]

  alp1 <- alpha[1]*beta1
  alp2 <- alpha[1]*(1-beta1)
  alp3 <- alpha[2]*beta2
  alp4 <- alpha[2]*(1-beta2)
  alp5 <- alpha[3]*beta3
  alp6 <- alpha[3]*(1-beta3)

  mu1 <- mu[1]
  mu2 <- mu[2]
  mu3 <- mu[3]
  mu4 <- mu[4]
  mu5 <- mu[5]
  mu6 <- mu[6]

  sig1 <- sigma[1]
  sig2 <- sigma[2]
  sig3 <- sigma[3]
  sig4 <- sigma[4]
  sig5 <- sigma[5]
  sig6 <- sigma[6]

  sig01 <- par0$sigma[1]
  sig02 <- par0$sigma[2]
  sig03 <- par0$sigma[3]

  ln <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) + alp5*dnorm(x, mu5, sig5) + alp6*dnorm(x, mu6, sig6) +1e-50))

  pen1 <- pn(sig1, sig01, an)+pn(sig2, sig01, an)+pn(sig3, sig02, an)+pn(sig4, sig02, an)+pn(sig5, sig03, an)+pn(sig6, sig03, an)

  pen2 <- log(1-abs(1-2*beta1))+log(1-abs(1-2*beta2))+log(1-abs(1-2*beta3))

  pln <- ln+pen1+pen2

  k <- 2
  output <- rep(0, k+1)
  output[1] <- pln

  for (i in 1:k) {

    pdf <- alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) + alp5*dnorm(x, mu5, sig5) + alp6*dnorm(x, mu6, sig6) +1e-50

    w1 <- (alp1*dnorm(x, mu1, sig1)+1e-50/6)/pdf
    w2 <- (alp2*dnorm(x, mu2, sig2)+1e-50/6)/pdf
    w3 <- (alp3*dnorm(x, mu3, sig3)+1e-50/6)/pdf
    w4 <- (alp4*dnorm(x, mu4, sig4)+1e-50/6)/pdf
    w5 <- (alp5*dnorm(x, mu5, sig5)+1e-50/6)/pdf
    w6 <- (alp6*dnorm(x, mu6, sig6)+1e-50/6)/pdf

    alpha[1] <- mean(w1+w2)
    alpha[2] <- mean(w3+w4)
    alpha[3] <- mean(w5+w6)


    mu1 <- sum(w1*x)/sum(w1)
    mu2 <- sum(w2*x)/sum(w2)
    mu3 <- sum(w3*x)/sum(w3)
    mu4 <- sum(w4*x)/sum(w4)
    mu5 <- sum(w5*x)/sum(w5)
    mu6 <- sum(w6*x)/sum(w6)


    sig1 <- sqrt( ( sum(w1*(x-mu1)^2)+2*an*sig01^2 )/(sum(w1)+2*an ) )
    sig2 <- sqrt( ( sum(w2*(x-mu2)^2)+2*an*sig01^2 )/(sum(w2)+2*an ) )
    sig3 <- sqrt( ( sum(w3*(x-mu3)^2)+2*an*sig02^2 )/(sum(w3)+2*an ) )
    sig4 <- sqrt( ( sum(w4*(x-mu4)^2)+2*an*sig02^2 )/(sum(w4)+2*an ) )
    sig5 <- sqrt( ( sum(w5*(x-mu5)^2)+2*an*sig03^2 )/(sum(w5)+2*an ) )
    sig6 <- sqrt( ( sum(w6*(x-mu6)^2)+2*an*sig03^2 )/(sum(w6)+2*an ) )

    if (sum(w1)/(sum(w1)+sum(w2))<=0.5) {
      beta1 <- min((sum(w1)+1)/(sum(w1)+sum(w2)+1), 0.5)
    }
    else {
      beta1 <- max((sum(w1))/(sum(w1)+sum(w2)+1), 0.5)
    }


    if (sum(w3)/(sum(w3)+sum(w4))<=0.5) {
      beta2 <- min( (sum(w3)+1)/(sum(w3)+sum(w4)+1), 0.5 )
    }
    else {
      beta2 <- max( (sum(w3))/(sum(w3)+sum(w4)+1), 0.5 )
    }

    if (sum(w5)/(sum(w5)+sum(w6))<=0.5) {
      beta3 <- min( (sum(w5)+1)/(sum(w5)+sum(w6)+1), 0.5 )
    }
    else {
      beta3 <- max( (sum(w5))/(sum(w5)+sum(w6)+1), 0.5 )
    }


    alp1 <- alpha[1]*beta1
    alp2 <- alpha[1]*(1-beta1)
    alp3 <- alpha[2]*beta2
    alp4 <- alpha[2]*(1-beta2)
    alp5 <- alpha[3]*beta3
    alp6 <- alpha[3]*(1-beta3)

    ln <- sum(log( alp1*dnorm(x, mu1, sig1) + alp2*dnorm(x, mu2, sig2) + alp3*dnorm(x, mu3, sig3)+alp4*dnorm(x, mu4, sig4) + alp5*dnorm(x, mu5, sig5) + alp6*dnorm(x, mu6, sig6) +1e-50))

    pen1 <- pn(sig1, sig01, an)+pn(sig2, sig01, an)+pn(sig3, sig02, an)+pn(sig4, sig02, an)+pn(sig5, sig03, an)+pn(sig6, sig03, an)

    pen2 <-  log(1-abs(1-2*beta1))+ log(1-abs(1-2*beta2))+ log(1-abs(1-2*beta3))

    pln <- ln+pen1+pen2

    output[i+1] <- pln
  }
  output
}


emtest3.norm <- function(x, eps=1e-6) {
  
  par0 <- phi3.norm(x, eps=eps)
  ln0 <- par0$loglik
  an <- formula3(par0, length(x))

  output <- c()

  for (beta1 in c(0.1, 0.3, 0.5)) {
    for (beta2 in c(0.1, 0.3, 0.5)) {
      for (beta3 in c(0.1, 0.3, 0.5)) {	
        beta <- c(beta1, beta2, beta3)
        par1 <- max6.norm(x, an=an, beta=beta, par0=par0, eps=eps)

        out <- emite6.norm(x, par1, par0, an=an)
        output <- rbind(output, out)
      }
    }
  }
  output <- 2*(output-ln0)

  emstat <- apply(output, 2, max)

  list(alpha=par0$alpha,
       theta=list(mu=par0$mu, sigma=par0$sigma),
       emstat=emstat,
       pvalue=pchisq(emstat, 6, lower.tail=F))
}


#### emtest.norm(x, m, eps=1e-6) can be used to test order=1, 2, 3. ###
emtest.norm2 <- function(x, m, eps=1e-6) {
  if (m == 1) 
    out <- emtest1.norm(x, eps=eps)
  
  if (m == 2) 
    out <- emtest2.norm(x, eps=eps)
  
  if (m == 3) 
    out <- emtest3.norm(x, eps=eps)

  if (m >= 4) {
    stop("Currently, we can only test for order=1, 2, or 3.")
  }

  output <- list(family='normal_unequalvar',
                 m0=m,
                 alpha=out$alpha,
                 theta=out$theta,
                 emstat=out$emstat,
                 pvalue=out$pvalue)
  class(output) <- 'emtest'

  output
} 

