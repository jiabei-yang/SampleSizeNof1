# Find J (number of pts per sequence) for ppobs graph
calc_j_popavg <- function(K, L, max.J, psbl.seq, intcpt, slp, 
                          corr.mtrx, sigma.resid.err, sigma.mu, sigma.tau, sigma.mutau,
                          type.1.err, power, min.diff) {
  
  Sigma.eps <- sigma.resid.err^2 * corr.mtrx
  
  if (is.na(max.J)) {
    
    actual.power <- 0
    J            <- 0
    
    while (actual.power < power) {
      J            <- J + 1
      actual.power <- calc_power(K, L, J, psbl.seq, intcpt, slp, 
                                 Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                                 min.diff, type.1.err)
    }
    
  } else { # !is.na(max.I)
    
    actual.power <- 1
    J            <- max.J + 1
    
    while (actual.power > power) {
      
      J <- J - 1
      # prev.power <- actual.power
      
      # when I = 0, it means actual.power for I = 1 is greater than power
      # so break the loop
      if (J == 0) {
        break
      }
      
      actual.power <- calc_power(K, L, J, psbl.seq, intcpt, slp, 
                                 Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                                 min.diff, type.1.err)
      
    }
    
    # Plus 1 back to reach the power needed
    J <- J + 1
    
    # assign prev.power greater than prespecified power to actual power for output
    # actual.power <- prev.power
    
  } # !is.na(max.I)
  
  return(c(K, L, J))
}



# Find L (number of measurements per period) for pts graph
calc_l_popavg <- function(J, K, max.L, gen.max.KL, psbl.seq, intcpt, slp,
                          sigma.resid.err, sigma.mu, sigma.tau, sigma.mutau,
                          corstr, rho, type.1.err, power, min.diff) {
  
  if (is.na(max.L)) {
    
    if (slp == "common") {
      
      actual.power <- 0
      L            <- 0
      
      while ((actual.power < power) & ((L + 1) <= floor(gen.max.KL / K))) {
        
        L <- L + 1
        
        corr.mtrx <- gen_corr_str(K, L, corstr, rho)
        Sigma.eps <- sigma.resid.err^2 * corr.mtrx
        
        actual.power <- calc_power(K, L, J, psbl.seq, intcpt, slp, 
                                   Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                                   min.diff, type.1.err)
        
        if (L %% 50 == 0) { print(c(K, L, J))}
        
      }
      
      # If L goes out of range, let L be NA
      if (actual.power < power) {
        L <- NA
      }
      
    } else { # (slp == "random") start from the largest L
      
      L <- floor(gen.max.KL / K) 
      corr.mtrx <- gen_corr_str(K, L, corstr, rho)
      Sigma.eps <- sigma.resid.err^2 * corr.mtrx
      
      actual.power <- calc_power(K, L, J, psbl.seq, intcpt, slp, 
                                 Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                                 min.diff, type.1.err)
      
      if (actual.power < power) {
        L <- NA
      } else {
        
        while(actual.power > power) {
          
          L <- L - 1
          
          # when I = 0, it means actual.power for I = 1 is greater than power
          # so break the loop
          if (L == 0) {
            break
          }
          
          corr.mtrx <- gen_corr_str(K, L, corstr, rho)
          Sigma.eps <- sigma.resid.err^2 * corr.mtrx
          
          actual.power <- calc_power(K, L, J, psbl.seq, intcpt, slp, 
                                     Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                                     min.diff, type.1.err)
          
          if (L %% 50 == 0) { print(c(K, L, J))}
          
        } # while
        
        # Plus 1 back to reach the power needed
        L <- L + 1
        
      }
      
    }
    
  } else { # !is.na(max.L)
    
    actual.power <- 1
    L            <- max.L + 1
    
    while (actual.power > power) {
      
      L <- L - 1
      # prev.power <- actual.power
      
      # when I = 0, it means actual.power for I = 1 is greater than power
      # so break the loop
      if (L == 0) {
        break
      }
      
      corr.mtrx <- gen_corr_str(K, L, corstr, rho)
      Sigma.eps <- sigma.resid.err^2 * corr.mtrx
      
      actual.power <- calc_power(K, L, J, psbl.seq, intcpt, slp, 
                                 Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                                 min.diff, type.1.err)
      
      if (L %% 50 == 0) { print(c(K, L, J))}
      
    }
    
    # Plus 1 back to reach the power needed
    L <- L + 1
    
    # assign prev.power greater than prespecified power to actual power for output
    # actual.power <- prev.power
    
  }
  
  return(c(K, L, J))
}



# calculate power for 1) different analysis and 2) randomization schemes when pop.avg treatment effect is of interest
calc_power <- function(K, L, J, psbl.seq, intcpt, slp, 
                       Sigma.eps, sigma.mu, sigma.tau, sigma.mutau, 
                       min.diff, type.1.err) {
  
  if (psbl.seq == "Alternating Sequences") {
    
    if ((intcpt == "fixed") & (slp == "common")) {
      actual.power <- calc_power_alter_fixintcpt_fixslp(K, L, J, Sigma.eps, min.diff, type.1.err)
      
    } else if ((intcpt == "fixed") & (slp == "random")) {
      D <- sigma.tau^2
      actual.power <- calc_power_alter_fixintcpt_randslp(K, L, J, Sigma.eps, D, min.diff, type.1.err)
      
    } else if ((intcpt == "random") & (slp == "common")) {
      D <- sigma.mu^2
      actual.power <- calc_power_alter_randintcpt_fixslp(K, L, J, Sigma.eps, D, min.diff, type.1.err)
      
    } else if ((intcpt == "random") & (slp == "random")) {
      D <- matrix(c(sigma.mu^2, sigma.mutau, sigma.mutau, sigma.tau^2), ncol = 2)
      actual.power <- calc_power_alter_randintcpt_randslp(K, L, J, Sigma.eps, D, min.diff, type.1.err)
    }
    
  } else if (psbl.seq == "Pairwise Randomization") {
    
    if ((intcpt == "fixed") & (slp == "common")) {
      actual.power <- calc_power_pair_fixintcpt_fixslp(K, L, J, Sigma.eps, min.diff, type.1.err)
      
    } else if ((intcpt == "fixed") & (slp == "random")) {
      D <- sigma.tau^2
      actual.power <- calc_power_pair_fixintcpt_randslp(K, L, J, Sigma.eps, D, min.diff, type.1.err)
      
    } else if ((intcpt == "random") & (slp == "common")) {
      D <- sigma.mu^2
      actual.power <- calc_power_pair_randintcpt_fixslp(K, L, J, Sigma.eps, D, min.diff, type.1.err)
      
    } else if ((intcpt == "random") & (slp == "random")) {
      D <- matrix(c(sigma.mu^2, sigma.mutau, sigma.mutau, sigma.tau^2), ncol = 2)
      actual.power <- calc_power_pair_randintcpt_randslp(K, L, J, Sigma.eps, D, min.diff, type.1.err)
    }
    
  }
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) alternating sequences 
# 2) fixed intercept; fixed slope
calc_power_alter_fixintcpt_fixslp <- function(K, L, J, Sigma.eps, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  Sigma.ij       <- Sigma.eps
  inv.Sigma.ij   <- solve(Sigma.ij)
  
  if ((K %% 2) == 0) {
    
    # For sum in the variance estimator
    for (j in 1:J){
      
      # 2 is for the alternating sequences
      X.ij.0      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(0, 1), K / 2), each = L)), ncol = J * 2 + 1)
      X.ij.0[, j] <- 1
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(1, 0), K / 2), each = L)), ncol = J * 2 + 1)
      X.ij.1[, j + J] <- 1
      
      sum.sigma.beta <- sum.sigma.beta + 
        t(X.ij.0) %*% inv.Sigma.ij %*% X.ij.0 +
        t(X.ij.1) %*% inv.Sigma.ij %*% X.ij.1
      
    }
    
  } else { # ((K %% 2) != 0)
    
    for (j in 1:J) {
      
      # 2 is for the alternating sequences
      X.ij.0      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = J * 2 + 1)
      X.ij.0[, j] <- 1
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = J * 2 + 1)
      X.ij.1[, j + J] <- 1
      
      sum.sigma.beta <- sum.sigma.beta + 
        t(X.ij.0) %*% inv.Sigma.ij %*% X.ij.0 +
        t(X.ij.1) %*% inv.Sigma.ij %*% X.ij.1
      
    }
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(rep(0, 2 * J), 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) alternating sequences 
# 2) random intercept; fixed slope
calc_power_alter_randintcpt_fixslp <- function(K, L, J, Sigma.eps, D, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  
  Z.ij         <- matrix(rep(1, K * L), ncol = 1)
  Sigma.ij     <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
  inv.Sigma.ij <- solve(Sigma.ij)
  
  if ((K %% 2) == 0) {
    
    # For sum in the variance estimator
    # 2 is for the alternating sequences
    X.ij.0      <- matrix(c(rep(1, K * L), rep(rep(c(0, 1), K / 2), each = L)), ncol = 2)
    X.ij.1      <- matrix(c(rep(1, K * L), rep(rep(c(1, 0), K / 2), each = L)), ncol = 2)
    
    # do not need J loop 
    # since will be the same for all participants assigned to the same sequence
    sum.sigma.beta <- sum.sigma.beta + 
      J * (t(X.ij.0) %*% inv.Sigma.ij %*% X.ij.0) +
      J * (t(X.ij.1) %*% inv.Sigma.ij %*% X.ij.1)
    
  } else { # ((K %% 2) != 0)
    
    # 2 is for the alternating sequences
    X.ij.0      <- matrix(c(rep(1, K * L), rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = 2)
    X.ij.1      <- matrix(c(rep(1, K * L), rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = 2)
    
    sum.sigma.beta <- sum.sigma.beta + 
      J * (t(X.ij.0) %*% inv.Sigma.ij %*% X.ij.0) +
      J * (t(X.ij.1) %*% inv.Sigma.ij %*% X.ij.1)
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(0, 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) alternating sequences 
# 2) fixed intercept; random slope
calc_power_alter_fixintcpt_randslp <- function(K, L, J, Sigma.eps, D, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  
  if ((K %% 2) == 0) {
    
    # For sum in the variance estimator
    for (j in 1:J){
      
      # 2 is for the alternating sequences
      X.ij.0      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(0, 1), K / 2), each = L)), ncol = J * 2 + 1)
      X.ij.0[, j] <- 1
      
      Z.ij.0     <- matrix(rep(rep(c(0, 1), K / 2), each = L), ncol = 1)
      Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(1, 0), K / 2), each = L)), ncol = J * 2 + 1)
      X.ij.1[, j + J] <- 1
      
      Z.ij.1     <- matrix(rep(rep(c(1, 0), K / 2), each = L), ncol = 1)
      Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
      
      sum.sigma.beta <- sum.sigma.beta + 
        t(X.ij.0) %*% solve(Sigma.ij.0) %*% X.ij.0 +
        t(X.ij.1) %*% solve(Sigma.ij.1) %*% X.ij.1
      
    }
    
  } else { # ((K %% 2) != 0)
    
    for (j in 1:J) {
      
      # 2 is for the alternating sequences
      X.ij.0      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = J * 2 + 1)
      X.ij.0[, j] <- 1
      
      Z.ij.0     <- matrix(c(rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = 1)
      Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = J * 2 + 1)
      X.ij.1[, j + J] <- 1
      
      Z.ij.1     <- matrix(c(rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = 1)
      Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
      
      sum.sigma.beta <- sum.sigma.beta + 
        t(X.ij.0) %*% solve(Sigma.ij.0) %*% X.ij.0 +
        t(X.ij.1) %*% solve(Sigma.ij.1) %*% X.ij.1
      
    }
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(rep(0, 2 * J), 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) alternating sequences 
# 2) random intercept; random slope
calc_power_alter_randintcpt_randslp <- function(K, L, J, Sigma.eps, D, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  
  if ((K %% 2) == 0) {
    
    # For sum in the variance estimator
    # 2 is for the alternating sequences
    X.ij.0     <- matrix(c(rep(1, K * L), rep(rep(c(0, 1), K / 2), each = L)), ncol = 2)
    Z.ij.0     <- X.ij.0
    Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
    
    X.ij.1     <- matrix(c(rep(1, K * L), rep(rep(c(1, 0), K / 2), each = L)), ncol = 2)
    Z.ij.1     <- X.ij.1
    Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
    
    sum.sigma.beta <- sum.sigma.beta + 
      J * (t(X.ij.0) %*% solve(Sigma.ij.0) %*% X.ij.0) +
      J * (t(X.ij.1) %*% solve(Sigma.ij.1) %*% X.ij.1)
    
  } else { # ((K %% 2) != 0)
    
    # 2 is for the alternating sequences
    X.ij.0     <- matrix(c(rep(1, K * L), rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = 2)
    Z.ij.0     <- X.ij.0
    Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
    
    X.ij.1     <- matrix(c(rep(1, K * L), rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = 2)
    Z.ij.1     <- X.ij.1
    Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
    
    sum.sigma.beta <- sum.sigma.beta + 
      J * (t(X.ij.0) %*% solve(Sigma.ij.0) %*% X.ij.0) +
      J * (t(X.ij.1) %*% solve(Sigma.ij.1) %*% X.ij.1)
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(0, 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) pairwise randomization
# 2) fixed intercept; fixed slope
calc_power_pair_fixintcpt_fixslp <- function(K, L, J, Sigma.eps, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  Sigma.ij       <- Sigma.eps
  inv.Sigma.ij   <- solve(Sigma.ij)
  
  if ((K %% 2) == 0) { 
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), K/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    # For sum in the variance estimator
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq[seq(2, K, 2)] <- abs(trt.seq[seq(2, K, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      for (j in 1:J){
        
        X.ij <- matrix(c(rep(0, K * L * J * I), trt.seq), ncol = J * I + 1)
        X.ij[, J * (i-1) + j] <- 1
        
        sum.sigma.beta <- sum.sigma.beta + 
          t(X.ij) %*% inv.Sigma.ij %*% X.ij 
      } # for j loop
      
    } # for i loop
    
  } else { # ((K %% 2) != 0)
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), (K + 1)/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq <- trt.seq[-length(trt.seq)]
      trt.seq[seq(2, K-1, 2)] <- abs(trt.seq[seq(2, K-1, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      for (j in 1:J) {
        
        X.ij <- matrix(c(rep(0, K * L * J * I), trt.seq), ncol = J * I + 1)
        X.ij[, J * (i-1) + j] <- 1
        
        sum.sigma.beta <- sum.sigma.beta + 
          t(X.ij) %*% inv.Sigma.ij %*% X.ij 
        
      } # for j loop
      
    } # for i loop
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(rep(0, I * J), 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) pairwise randomization
# 2) random intercept; fixed slope
calc_power_pair_randintcpt_fixslp <- function(K, L, J, Sigma.eps, D, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  
  Z.ij         <- matrix(rep(1, K * L), ncol = 1)
  Sigma.ij     <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
  inv.Sigma.ij <- solve(Sigma.ij)
  
  if ((K %% 2) == 0) {
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), K/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq[seq(2, K, 2)] <- abs(trt.seq[seq(2, K, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
      
      # do not need J loop 
      # since will be the same for all participants assigned to the same sequence
      sum.sigma.beta <- sum.sigma.beta + 
        J * (t(X.ij) %*% inv.Sigma.ij %*% X.ij) 
      
    } # for i loop
    
  } else { # ((K %% 2) != 0)
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), (K + 1)/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq <- trt.seq[-length(trt.seq)]
      trt.seq[seq(2, K-1, 2)] <- abs(trt.seq[seq(2, K-1, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
      
      # do not need J loop 
      # since will be the same for all participants assigned to the same sequence
      sum.sigma.beta <- sum.sigma.beta + 
        J * (t(X.ij) %*% inv.Sigma.ij %*% X.ij) 
      
    } # for i loop
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(0, 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) pairwise randomization
# 2) fixed intercept; random slope
calc_power_pair_fixintcpt_randslp <- function(K, L, J, Sigma.eps, D, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  
  if ((K %% 2) == 0) {
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), K/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq[seq(2, K, 2)] <- abs(trt.seq[seq(2, K, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      for (j in 1:J){
        
        X.ij <- matrix(c(rep(0, K * L * J * I), trt.seq), ncol = J * I + 1)
        X.ij[, J * (i-1) + j] <- 1
        
        Z.ij <- matrix(trt.seq, ncol = 1)
        Sigma.ij <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
        
        sum.sigma.beta <- sum.sigma.beta + 
          t(X.ij) %*% solve(Sigma.ij) %*% X.ij 
      } # for j loop
      
    } # for i loop
    
  } else { # ((K %% 2) != 0)
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), (K + 1)/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq <- trt.seq[-length(trt.seq)]
      trt.seq[seq(2, K-1, 2)] <- abs(trt.seq[seq(2, K-1, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      for (j in 1:J) {
        
        X.ij <- matrix(c(rep(0, K * L * J * I), trt.seq), ncol = J * I + 1)
        X.ij[, J * (i-1) + j] <- 1
        
        Z.ij <- matrix(trt.seq, ncol = 1)
        Sigma.ij <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
        
        sum.sigma.beta <- sum.sigma.beta + 
          t(X.ij) %*% solve(Sigma.ij) %*% X.ij 
        
      } # for j loop
      
    } # for i loop
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(rep(0, I * J), 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}



# Calculate the power for 
# 1) pairwise randomization
# 2) random intercept; random slope
calc_power_pair_randintcpt_randslp <- function(K, L, J, Sigma.eps, D, min.diff, type.1.err) {
  
  sum.sigma.beta <- 0
  
  if ((K %% 2) == 0) {
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), K/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq[seq(2, K, 2)] <- abs(trt.seq[seq(2, K, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
      Z.ij <- X.ij
      Sigma.ij <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
      
      sum.sigma.beta <- sum.sigma.beta + 
        J * (t(X.ij) %*% solve(Sigma.ij) %*% X.ij)
      
    } # for i loop
    
  } else { # ((K %% 2) != 0)
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), (K + 1)/2)
    seq.pair <- expand.grid(seq.pair)
    
    I <- nrow(seq.pair)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq <- trt.seq[-length(trt.seq)]
      trt.seq[seq(2, K-1, 2)] <- abs(trt.seq[seq(2, K-1, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
      Z.ij <- X.ij
      Sigma.ij <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
      
      sum.sigma.beta <- sum.sigma.beta + 
        J * (t(X.ij) %*% solve(Sigma.ij) %*% X.ij)

    } # for i loop
    
  }
  
  Var.theta <- solve(sum.sigma.beta)
  C         <- c(0, 1)
  se.tau    <- sqrt(t(C) %*% Var.theta %*% C)
  
  actual.power <- pnorm(- min.diff / se.tau - qnorm(1 - type.1.err / 2)) + pnorm(min.diff / se.tau - qnorm(1 - type.1.err / 2))
  
  return(actual.power)
  
}

