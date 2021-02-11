#############################################################
# Algorithm 1
############################################################# 
FASA_linear1 = function(Y, X, subsize, pilot.ind, r1, r2)
{   
   n = dim(X)[1]; p = dim(X)[2]; k = log2(n)
   H_candidate = runif(p * r2)
   First_H = which(H_candidate <= 1/6); Second_H = which(5/6 < H_candidate) 
   H_candidate[First_H] = sqrt(3/r2); H_candidate[Second_H] = -sqrt(3/r2)
   H_candidate[ -c(First_H, Second_H)] = 0
   Pi_2 = matrix(H_candidate,p,r2)        
      
   Dd = (2*rbinom(n,1,0.5)-1)/sqrt(r1);
   DX = X * Dd; 
   S = rep(0, n)
   S[sample(1:n, r1, F)] = 1
   SX = hadamard_cpp(S, DX)
      
   SAOP = solve(t(SX) %*% SX, Pi_2)
      
   ##U_idx <- sample(1:n, r0, T)
   U_idx = pilot.ind; r0 = length(pilot.ind)
   x.U <- X[U_idx,]
   y.U <- Y[U_idx] 
   pilot.est = solve(t(x.U)%*%x.U,t(x.U)%*%y.U)
   Epsil_r0 <- Y - X%*% pilot.est
   
   Vari <- sum(Epsil_r0^2)/(n)
   Second <- (Epsil_r0^2 - Vari)^2/(4 * n^2 * Vari)
   p.slev <- sqrt(Epsil_r0^2 * rowSums((X %*% SAOP)^2) + Second)      
   p.slev <- p.slev/sum(p.slev)
   idx.lv <- sample(1:n, subsize, T, p.slev)
   x.lv <- X[c(idx.lv,U_idx),]
   y.lv <- Y[c(idx.lv,U_idx)]
   pi.i <- c(p.slev[idx.lv],rep(1/n,r0))
   
   eig.lv <- eigen( t(x.lv) %*% (1/pi.i * x.lv))
   iI.lv <- eig.lv$vectors %*% (t(eig.lv$vectors)/eig.lv$values)
      
   iI.lv %*% t(x.lv) %*% (1/ pi.i * y.lv)
   beta.FastAop <- iI.lv %*% t(x.lv) %*% (1/ pi.i * y.lv)
   
   Epsil_r = y.lv - c(x.lv%*%beta.FastAop)
   Inv.sigma_x = iI.lv * ( n * (subsize + r0) ) 
   sigma_pi = ( t(x.lv) %*% (x.lv*(Epsil_r^2)/(pi.i^2)) ) / (( n * (subsize + r0) )^2) 
   sigma.FastAop = diag(Inv.sigma_x %*% sigma_pi %*% Inv.sigma_x)
   vari = sum(Epsil_r^2/pi.i)/( sum(1/pi.i) )
   list(beta.FastAop = beta.FastAop, sigma.FastAop = sigma.FastAop, vari.FastAop = vari)
}
#############################################################
# Algorithm 2. 
############################################################# 
FASA_linear2 = function(Y, X, subsize, pilot.ind, r2)
{   
   n = dim(X)[1]; p = dim(X)[2]; k = log2(n)
   H_candidate = runif(p * r2)
   First_H = which(H_candidate <= 1/6); Second_H = which(5/6 < H_candidate) 
   H_candidate[First_H] = sqrt(3/r2); H_candidate[Second_H] = -sqrt(3/r2)
   H_candidate[ -c(First_H, Second_H)] = 0
   Pi_2 = matrix(H_candidate,p,r2)        
      
   U_idx = pilot.ind; r0 = length(pilot.ind)
   x.U = X[U_idx,]
   SAOP = solve(t(x.U) %*% (x.U/r0),Pi_2)
   
   y.U <- Y[U_idx] 
   pilot.est = solve(t(x.U)%*%x.U,t(x.U)%*%y.U)
   Epsil_r0 <- Y - X%*% pilot.est
      
   Vari <- sum(Epsil_r0^2)/(n)
   Second <- (Epsil_r0^2 - Vari)^2/(4 * n^2 * Vari)
      
   p.slev <- sqrt(Epsil_r0^2 * rowSums((X %*% SAOP)^2) + Second)
   p.slev <- p.slev/sum(p.slev)
   idx.lv <- sample(1:n, subsize, T, p.slev)
   x.lv <- X[c(idx.lv,U_idx),]
   y.lv <- Y[c(idx.lv,U_idx)]
   pi.i <- c(p.slev[idx.lv],rep(1/n,r0))

   eig.lv <- eigen( t(x.lv) %*% (1/pi.i * x.lv))
   iI.lv <- eig.lv$vectors %*% (t(eig.lv$vectors)/eig.lv$values)
      
   iI.lv %*% t(x.lv) %*% (1/ pi.i * y.lv)
   beta.FastSAop <- iI.lv %*% t(x.lv) %*% (1/ pi.i * y.lv)
   
   Epsil_r = y.lv - c(x.lv%*%beta.FastSAop)

   Inv.sigma_x = iI.lv * ( n * (subsize + r0) ) 
   sigma_pi = ( t(x.lv) %*% (x.lv*(Epsil_r^2)/(pi.i^2)) ) / (( n * (subsize + r0) )^2) 

   sigma.FastSAop = diag(Inv.sigma_x %*% sigma_pi %*% Inv.sigma_x)
   vari = sum(Epsil_r^2/pi.i)/( sum(1/pi.i) )
   list(beta.FastSAop = beta.FastSAop, sigma.FastSAop = sigma.FastSAop, vari.FastAop = vari)
}
