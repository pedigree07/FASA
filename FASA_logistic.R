#########################################################################
################### MLE
#########################################################################
MLE <- function(x, y, w, maxiter = 100, tol = 1e-5) 
{
   beta <- rep(0, ncol(x))
   msg <- c()
   for(iter in 1:maxiter) 
   {
      p <- c(1 - 1 / (1 + exp(x %*% beta)))
      H <- t(x) %*% (p * (1 - p) * w * x)
      S <- colSums((y - p) * w * x)
      tryCatch({
         shs <- NA
         shs <- solve(H, S) },
         error=function(e){cat("\n ERROR :", iter, conditionMessage(e), "\n")})
      if (is.na(shs[1])) 
      {
         msg <- "Not converge"
         beta <- NA
         break
      }
      beta.new <- beta + shs
      diff = sum((beta.new - beta)^2)
      beta  <- beta.new
      if(diff < tol)
      {
         msg <- "converge"
         break
      }
      iter  <- iter
   }
   list(beta=beta, message=msg, iter=iter, diff=diff)
}
#########################################################################
################### Algorithm 1
#########################################################################
FASA_logistic1 = function(Y, X, subsize, pilot.ind, r1, r2)
{
   n = dim(X)[1]; p = dim(X)[2]; k = log2(n)

   idx = pilot.ind; r0 = length(pilot.ind)
   x.samp <- X[idx,]
   y.samp <- Y[idx]
   w = n
   fit.uni = MLE(x = x.samp, y = y.samp, w = w)
   beta.uni <- fit.uni$beta
   if (anyNA(beta.uni)) 
   {
      result <- list(time = time.FAop, beta.FAop=beta.uni, sigma.FAop=sigma.FAop, msg="first stage not converge")
   }else{
      
   P.uni  = 1 - 1 / (1 + exp(c(X %*% beta.uni)))
   w.uni = P.uni * (1 - P.uni)
   sqW.X = X * sqrt(w.uni)
   Dd = (2*rbinom(n,1,0.5)-1)/sqrt(r1);
   D.sqW.X = sqW.X * Dd; 
   S = rep(0, n)
   S[sample(1:n, r1, F)] = 1
   SX = hadamard_cpp(S, D.sqW.X)
   H_candidate = runif(p * r2)
   First_H = which(H_candidate <= 1/6); Second_H = which(5/6 < H_candidate) 
   H_candidate[First_H] = sqrt(3/r2); H_candidate[Second_H] = -sqrt(3/r2)
   H_candidate[ -c(First_H, Second_H)] = 0
   Pi_2 = matrix(H_candidate,p,r2)        
   
   SAOP = solve(t(SX) %*% SX, Pi_2)
         
   P.FAop <- sqrt((Y - P.uni)^2 * rowSums((X %*% SAOP)^2))
   P.FAop <- P.FAop / sum(P.FAop)
         
   idx.FAop <- sample(1:n, subsize, T, P.FAop)
   x.FAop <- X[c(idx.FAop, idx),]
   y.FAop <- Y[c(idx.FAop, idx)]
   pinv.FAop <- c(1 / P.FAop[idx.FAop], rep(n,r0))
   fit.FAop <- MLE(x=x.FAop, y=y.FAop, w=pinv.FAop)
         
   Tsubsize <- length(pinv.FAop)
   beta.FAop <- fit.FAop$beta
   
     
   p.FAop  <- 1 - 1 / (1 + exp(c(x.FAop %*% fit.FAop$beta)))
   w.FAop <- p.FAop * (1 - p.FAop)
   M.FAop <- solve(t(x.FAop) %*% (x.FAop * (w.FAop * pinv.FAop))) * Tsubsize * n
   Vc.FAop <- t(x.FAop) %*% (x.FAop * (y.FAop-p.FAop)^2 * pinv.FAop^2) / Tsubsize^2 / n^2
   V.FAop <- M.FAop %*% Vc.FAop %*% M.FAop
   sigma.FAop <- diag(V.FAop)
   msg <- c(fit.uni$message, fit.FAop$message)
   result <- list(beta.FAop=beta.FAop, sigma.FAop=sigma.FAop, msg=msg)
   }
   return(result)
}
#########################################################################
################### Algorithm 2
#########################################################################
FASA_logistic2 = function(Y, X, subsize, pilot.ind, r2)
{
   n = dim(X)[1]; p = dim(X)[2]; k = log2(n)

   idx = pilot.ind; r0 = length(pilot.ind)
   x.samp <- X[idx,]
   y.samp <- Y[idx]
   w = n
   fit.uni = MLE(x = x.samp, y = y.samp, w = w)
   beta.uni <- fit.uni$beta
   if (anyNA(beta.uni)) 
   {
      result <- list(time = time.FAop, beta.FAop=beta.uni, sigma.FAop=sigma.FAop, msg="first stage not converge")
   }else{

   H_candidate = runif(p * r2)
   First_H = which(H_candidate <= 1/6); Second_H = which(5/6 < H_candidate) 
   H_candidate[First_H] = sqrt(3/r2); H_candidate[Second_H] = -sqrt(3/r2)
   H_candidate[ -c(First_H, Second_H)] = 0
   Pi_2 = matrix(H_candidate,p,r2)        
         
   SX = X[idx,]
         
   P.uni  = 1 - 1 / (1 + exp(c(X %*% beta.uni)))
   w.uni = P.uni * (1 - P.uni)
   SAOP = solve(t(SX) %*% (SX * c(w.uni[idx])), Pi_2)
         
   P.FAop <- sqrt((Y - P.uni)^2 * rowSums((X%*%SAOP)^2))
   P.FAop <- P.FAop / sum(P.FAop)
         
   idx.FAop <- sample(1:n, subsize, T, P.FAop)
   x.FAop <- X[c(idx.FAop, idx),]
   y.FAop <- Y[c(idx.FAop, idx)]
   pinv.FAop <- c(1 / P.FAop[idx.FAop], rep(n,r0))
   fit.FSAop <- MLE(x=x.FAop, y=y.FAop, w=pinv.FAop)
         
   Tsubsize <- length(pinv.FAop)
   beta.FSAop <- fit.FSAop$beta
   
   p.FAop  <- 1 - 1 / (1 + exp(c(x.FAop %*% fit.FSAop$beta)))
   w.FAop <- p.FAop * (1 - p.FAop)
   M.FAop <- solve(t(x.FAop) %*% (x.FAop * (w.FAop * pinv.FAop))) * Tsubsize * n
   Vc.FAop <- t(x.FAop) %*% (x.FAop * (y.FAop-p.FAop)^2 * pinv.FAop^2) / Tsubsize^2 / n^2
   V.FAop <- M.FAop %*% Vc.FAop %*% M.FAop
   sigma.FSAop <- diag(V.FAop)
   msg <- c(fit.uni$message, fit.FSAop$message)
   result <- list(beta.FSAop=beta.FSAop, sigma.FSAop=sigma.FSAop, msg=msg)
   }
   return(result)
}