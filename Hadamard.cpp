#include <RcppArmadillo.h>
#include <time.h>
#include <stdlib.h>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat  hadamard_cpp(LogicalVector S, mat v){
   int n = v.n_rows;
   if (n==1){
      return(v);
   } else { 
   LogicalVector S1 = S[seq(0,n/2-1)];
   LogicalVector S2 = S[seq(n/2,n-1)];
   bool FLogic = is_true(any(S1)); 
   bool SLogic = is_true(any(S2)); 
   mat v1 = v.rows(0,n/2-1);
   mat v2 = v.rows(n/2,n-1);   
   if( FLogic == TRUE && SLogic == TRUE)
   {  
      mat mu = hadamard_cpp(S1,v1+v2);
      mat ml = hadamard_cpp(S2,v1-v2);
      return(join_cols(mu,ml));
   }else if(FLogic == TRUE && SLogic == FALSE){   
      mat mu = hadamard_cpp(S1, v1+v2);
      return(mu);
   }else if(FLogic == FALSE && SLogic == TRUE){
      mat ml = hadamard_cpp(S2,v1-v2);
      return(ml);   
   }else{ 
      mat M;
      return M;
      }
   }
}   
   