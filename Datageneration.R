library(mvtnorm)
Datageneration = function(n, beta, vari, corr, Dist = 'mzNormal')
{
   beta0 <- beta.0  <- beta
   d <- p <- length(beta.0)
   sigmax <- matrix(corr, d, d) + diag(1-corr, d)
   if( Dist == 'mzNormal' ){X  <- rmvnorm(n, rep(0, d), sigmax)}
   if( Dist == 'nzNormal' ){X  <- rmvnorm(n, rep(1, d), sigmax)}
   if( Dist == 'ueNormal' )
   {
      sigmax <- matrix(corr, d, d) + diag(1/((1:d)^(-0.1))-corr, d)
      X  <- rmvnorm(n, rep(0, d), sigmax)
   }
   if( Dist == 'T3' ){X  <- rmvt(n, sigma = sigmax, df = 3)}
   if( Dist == 'T2' ){X  <- rmvt(n, sigma = sigmax, df = 2)}
   if( Dist == 'T1' ){X  <- rmvt(n, sigma = sigmax, df = 1)}
   if( Dist == 'mixNormal' )
   {
      X = matrix(0, n, d)
      U=runif(n)
      X[which(U>0.5),]  = rmvnorm( sum(U > 0.5), rep(1, d), sigmax)
      X[which(U<=0.5),]  = rmvnorm( sum(U <= 0.5), rep(-1, d), sigmax) 
   }
   if( Dist == 'EXP' )
   {
      X = c()
      for(i in 1:d)
         X = cbind(rexp(n, 1/2), X)
   }
   X = cbind(1, X)
   beta0 <- c(1, beta.0)
   p <- d + 1
   Y  <- c(X %*% beta0) + rnorm(n, 0, vari)
   list(Y = Y, X = X)
}   
