#!/usr/bin/env Rscript
path = getwd();
source(paste(path, "/install.r", sep=""));
source(paste(path, "/heston2.R", sep=""));

#  Heston Model Calibration
# 
#  Options2011: set of options with on index
#  
#  Original author: Jonathan Frei, 2015
#  Author: 0xAlbert, 2023
# 

#  K is a vector of strikes 
#  PC is a vector containing either a 1 or 2 (1 for call, 2 for put)
#  Price is a vector of market prices of the options 
#  S is a vector of spot prices
#  T is a vector of maturities (must be the same unit as r and q)
#  dateO is the starting date of the option
#  r is a vector of interest rates
#  q is a vector of dividend rates


# PC==1 for calls, PC==2 is for puts
# dateO == min(dateO) selects the first dates only (Wednesdays here)
data <- readMat("file.mat");
indx = which(data$PC == 2 & data$dateO == min(data$dateO));

K <- data$K;
S <- data$S;
T <- data$T;
t <- data$t;
r <- data$r;
q <- data$q;
PC <- data$PC;
Price <- data$Price;

# % reducing vectors to selected options
K = K[indx];

S = S[indx];
T = T[indx];
t = replicate(length(T), 0);
r = r[indx];
q = q[indx];
PC = PC[indx];
MarketPrice = Price[indx];

# tic;
# Tip: Smaller upper boundaries of [100 100 1-eps 100 100 ] instead of [Inf Inf 1-eps Inf Inf ] and trying 
# different starting parameters might help if the calibration gets stuck in a local minimum.
startparameters = c(
    0.2^2, 
    0.1, 
    -0.1, 
    0.25^2, 
    1
);

getPrice <- 
function(            
    startparameters
) {
    Heston1993KahlJaeckelLordRev3(
        PC, 
        S, 
        K, 
        T, 
        t, 
        r, 
        q,  
        startparameters[1], # v0 initial variance
        startparameters[2], # theta long run mean variance
        startparameters[3], # kappa mean reversion speed of  volatility
        startparameters[4], # sigma volatility of volatility
        startparameters[5], # rho correlation between returns volatility
        # alphas = c(0) # alpha can be a vector supplied by the user, otherwise the function attempts to find a payoff-dependent optimal alpha
    );
}

tic;

tryCatch({

    controlpars <- nls.lm.control(
        ftol = sqrt(.Machine$double.eps),
        ptol = sqrt(.Machine$double.eps), 
        gtol = 0, 
        diag = list(), 
        epsfcn = 0,
        factor = 100, 
        maxfev = integer(), 
        maxiter = 50, 
        nprint = 0
    )

    xopt = nls.lm(startparameters, 
        lower=c(eps(), eps(), -1 + eps(), eps(),  eps()), 
        upper=c(100, 100, 1 - eps(), 100, 100), 
        # upper=c(Inf, Inf, 1 - eps(), Inf, Inf), 
        getPrice,
        control = controlpars
    );
    
    print(xopt);

}, error = function(e) {
    if (length(e) > 1) {
        print(e);
    }
})   

toc;
