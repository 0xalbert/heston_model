#!/usr/bin/env Rscript
path = getwd()
source(paste(path, "/install.r", sep=""))

Heston1993KahlJaeckelLordRev3 <- 
function(PC, S, K, T, t, r, q, v0, theta, rho, kappa, sigma, alphas=NA) {
    #  Heston pricing function based on the implementation suggested by
    #  Roger Lord and Chrisitan Kahl in "Optimal Fourier inversion in 
    #  semi-analytical option pricing"
    # 
    #  Input: (PC till q can be vectorized)
    #        PC: 1 for Calls, 2 for Puts
    #        S: Spot
    #        K: Strike
    #        T: Maturity
    #        t: start date
    #        r: interest rate
    #        q: dividend
    #        v0: initial variance
    #        theta: long run mean variance
    #        kappa: mean reversion speed of  volatility
    #        sigma: volatility of volatility
    #        rho: correlation between returns volatility
    #        alpha: alpha can be a vector supplied by the user, otherwise the
    #        function attempts to find a payoff-dependent optimal alpha
    # 
    #  Output: Price for each option, optionally generated alphas
    # 
    #  Usage: Heston1993KahlJaeckelLordRev3(1, 100, 100, 20,0, 0.05, 0.0, 
    #         0.00003, 0.00003,-0.3, 0.5, 0.0008)
    # 
    #  Original author: Jonathan Frei, 2015
    #  Author: 0xAlbert, 2023
    #  

    print(glue::glue("got v0 {v0}  theta {theta}  rho {rho}  kappa {kappa}  sigma {sigma}  alphas {alphas}"))

    # force column vector
    PC <- as.numeric(PC)
    S <- as.numeric(S)
    K <- as.numeric(K)
    T <- as.numeric(T)
    t <- as.numeric(t)
    r <- as.numeric(r)
    q <- as.numeric(q)

    nos <- length(S)
    prices <- rep(NA, nos)
    tau <- T - t
    mu <- (r - q)
    F <- S * exp(mu * tau)

    if(is.na(alphas)) {
        alphas <- rep(NA, nos)
    } else if(length(alphas) == 1) {
        alphas <- rep(alphas, nos)
    }
    

    alpha0 <- 0.75

    for(ind in 1:nos) {

      getPsi <- function(alpha0) {
        psi(
          alpha0, 
          K = K[ind], 
          F = F[ind], 
          kappa = kappa, 
          theta = theta, 
          rho = rho, 
          sigma = sigma, 
          tau = tau[ind], 
          v0 = v0
        )
      }
      # print(alphas)

      if(is.na(alphas[ind])) {
        tryCatch({

          alphas[ind] <- fzero(
            getPsi,
            alpha0
          )$fval
          
        }, error = function(e) {
          alphas[ind] <- alpha0
        })
      }
    
      print(glue::glue("alpha: {alphas[ind]} K: {K[ind]} F: {F[ind]} kappa: {kappa} theta: {theta} rho: {rho} sigma: {sigma} tau: {tau[ind]} v0: {v0}"))
          
        tryCatch({

          integral <- integrate(
              phi,
              0,
              Inf,
              subdivisions = 5000,
              abs.tol = 1e-5, #rel.tol = .Machine$double.eps^.05,
              rel.tol = 1e-7,  #abs.tol = .Machine$double.eps^.01,
              stop.on.error = TRUE,
              K = K[ind], 
              alpha = alphas[ind], 
              F = F[ind], 
              kappa = kappa, 
              theta = theta, 
              rho = rho, 
              sigma = sigma, 
              tau = tau[ind], 
              v0 = v0
          )$value;

        }, error = function(e) {
            # if (length(e) > 1) {
            #   print(e)
            # }
        })

        if (!is.numeric(integral)) {
          integral <- 0
        }
        
        prices[ind] <- Ralpha(F[ind], K[ind], alphas[ind]) + 1 / pi * integral
        # print(integral)

        if (PC[ind] == 2) {
          prices[ind] <- prices[ind] + K[ind] * exp(-r[ind] * tau[ind])-S[ind] * exp(-q[ind]*tau[ind])
        }
    }

    # # Optionally return alphas 
    # # list(alphas = alphas, prices = prices)
    return(prices)
}

# f <- function(x) {
#    if (all(Im(z <- zapsmall(x)) == 0 )) as.numeric(z) else x
# }

psi <- function(alpha, K, F, kappa, theta, rho, sigma, tau, v0) {
    k <- log(K)
    p <- -alpha * k + 0.5 * log(phi(-(alpha + 1) * 1i, K, alpha, F, kappa, theta, rho, sigma, tau, v0) ^2)
}

Ralpha <- function(F, K, alpha) { 
    r = F * (alpha <= 0) - K * (alpha <= -1) - 0.5 * (F * (alpha == 0) -K * (alpha == -1));
}

cf <- function(u, F, kappa, theta, rho, sigma, tau, v0) { 
    f = log(F);
    c = exp(1i * u * f + A(u, kappa, theta, rho, sigma, tau) + Bv(u, rho, sigma, kappa, tau) * v0);
}

phi <- function(v, K, alpha, F, kappa, theta, rho, sigma, tau, v0) {
    k = log(K);
    y <- Re(exp(-1i * (v - 1i * alpha) * k) * cf(v - 1i * (alpha + 1), F, kappa, theta, rho, sigma, tau, v0) / ((v - 1i * (alpha + 1)) * (v - 1i * alpha)));
}

Bv <- function(u, rho, sigma, kappa, tau) { 
  b = ((beta(u,rho,sigma,kappa) - D(u, rho, sigma, kappa)) * (1 - exp(-D(u, rho, sigma, kappa)  * tau)))/(sigma^2 * (1 - G(u, rho, sigma, kappa) * exp(-D(u, rho, sigma, kappa)*tau)));
  return(b)
}

A <- function(u, kappa, theta, rho, sigma, tau) { 
    f <- beta(u, rho, sigma, kappa) - D(u, rho, sigma, kappa)
    c <- f * tau
    d <- c - (2 * log(phi2(u, rho, sigma, kappa, tau)))
    e <- kappa * theta * d
    a <- e / (sigma^2)
}

phi2 <- function(u, rho, sigma, kappa, tau) { 
  p <- (G(u, rho, sigma, kappa) * exp(-D(u, rho, sigma, kappa) * tau) -1) / (G(u, rho, sigma, kappa) -1)
  return(p)
}

G <- function(u, rho, sigma, kappa) { 
    g = (beta(u, rho, sigma, kappa) -D(u, rho, sigma, kappa)) / (beta(u, rho, sigma, kappa) + D(u, rho, sigma, kappa));
}

D <- function(u, rho, sigma, kappa) { 
    d = beta(u, rho, sigma, kappa) ^2
    x <- d * -4 * alphahat(u) * gamma(sigma);
    x <- sqrt(x);
}

alphahat <- function(u) { 
    x <- -0.5 * u 
    a = x * (1i + u);
}

beta <- function(u, rho, sigma, kappa) { 
    b = kappa - rho * sigma * u * 1i;
}

gamma <- function(sigma) {
    x <- 0.5 * sigma
    y = x^2;
}

Heston1993KahlJaeckelLordRev3(
  2, #1 Put, 2 Call
  100, #S
  150, #K
  180, #F
  0, 
  0, 
  0, 
  0.05, 
  0.0, 
  0.00003, 
  0.00003, 
  -0.3, 
  0.5
)

# data <- readMat('file.mat')

# str(data)
# print(typeof(data))
# dfs <- lapply(data, data.frame, stringsAsFactors = FALSE)
# df <- bind_rows(df) 
# bind_cols(dfs) 


# write_xlsx(df, paste(path, "/file.xls", sep=""))
