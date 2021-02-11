source('Datageneration.R')
Rcpp::sourceCpp('Hadamard.cpp')
source("FASA_linear.R")
###################################################################
######################### Data Generation #########################
###################################################################
n <- 2^17
beta  <- c(rep(3,10))
Data = Datageneration(n, beta, vari = 3, corr = 0.5, Dist = 'mzNormal'); 

######################### Pilot sample for using algorithm 1 and 2 ###########
pilot.ind <- sample(1:n, 400, T)
######################### Algorithm1 for linear model ########################
FASA1 = FASA_linear1(Data$Y, Data$X, subsize = 500, pilot.ind = pilot.ind, r1 = 1000, r2 = 10)
######################### Algorithm2 for linear model ########################
FASA2 = FASA_linear2(Data$Y, Data$X, subsize = 500, pilot.ind = pilot.ind, r2 = 10)
