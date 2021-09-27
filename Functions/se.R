calc_se_niv <- function(K, L, psbl.seq, sigma.resid.err, rho, corr.str) {
  
  corr.mtrx <- gen_corr_str(K, L, corr.str, rho)
  Sigma.eps <- sigma.resid.err^2 * corr.mtrx
  inv.Sigma.eps <- solve(Sigma.eps)
  
  C <- c(0, 1)
  
  # design matrix X can be different based on different sequences
  # which may result in different standard error
  se <- NULL
  if (psbl.seq == "Alternating Sequences") {
    
    if ((K %% 2) == 0) {
      
      # 1st scenario: start with 0
      X.0 <- matrix(c(rep(1, K * L), rep(rep(c(0, 1), K / 2), each = L)), ncol = 2)
      tmp.se <- sqrt(t(C) %*% solve(t(X.0) %*% inv.Sigma.eps %*% X.0) %*% C)
      
      # only "01" in the 1st block needed 
      # because we know the following periods in alternating sequences
      se <- rbind(se, 
                  c("01", tmp.se))  
      
      # 2nd scenario: start with 1
      X.1 <- matrix(c(rep(1, K * L), rep(rep(c(1, 0), K / 2), each = L)), ncol = 2)
      tmp.se <- sqrt(t(C) %*% solve(t(X.1) %*% inv.Sigma.eps %*% X.1) %*% C)
      
      se <- rbind(se, 
                  c("10", tmp.se))  
      
    } else { # ((K %% 2) != 0)
      
      # 1st scenario: start with 0
      X.0 <- matrix(c(rep(1, K * L), rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = 2)
      tmp.se <- sqrt(t(C) %*% solve(t(X.0) %*% inv.Sigma.eps %*% X.0) %*% C)
      
      se <- rbind(se, 
                  c("01", tmp.se))  
      
      # 2nd scenario: start with 1
      X.1 <- matrix(c(rep(1, K * L), rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = 2)
      tmp.se <- sqrt(t(C) %*% solve(t(X.1) %*% inv.Sigma.eps %*% X.1) %*% C)
      
      se <- rbind(se, 
                  c("10", tmp.se))  
      
    }
    
  } else if (psbl.seq == "Pairwise Randomization") {
    
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
        
        # same design matrix for all j's
        # only need to calculate for different sequences
        X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
        
        tmp.se <- sqrt(t(C) %*% solve(t(X.ij) %*% inv.Sigma.eps %*% X.ij) %*% C)
        
        se <- rbind(se, 
                    c(paste(seq.pair[i, ], collapse = ""), tmp.se)) 
        
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
        
        # same design matrix for all j's
        # only need to calculate for different sequences
        X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
        
        tmp.se <- sqrt(t(C) %*% solve(t(X.ij) %*% inv.Sigma.eps %*% X.ij) %*% C)
        
        se <- rbind(se, 
                    c(paste(seq.pair[i, ], collapse = ""), tmp.se)) 
        
      } # for i loop
      
    } # else ((K %% 2) != 0)
    
  } # else if (psbl.seq == "Pairwise Randomization")
  
  se <- matrix(se[!duplicated(se[, 2]), ], ncol = 2)
  
  return(se)
  
}



calc_se_shrk <- function(K, L, J, psbl.seq, intcpt, 
                         sigma.resid.err, rho, corr.str, sigma.mu, sigma.tau, sigma.mutau, var.rand.eff) {
  
  # var.rand.eff: take into account the variation in random effect
  
  corr.mtrx <- gen_corr_str(K, L, corr.str, rho)
  Sigma.eps <- sigma.resid.err^2 * corr.mtrx
  
  if (psbl.seq == "Alternating Sequences") {
    
    if (intcpt == "fixed") {
      C.b <- 1
      D   <- sigma.tau^2
      se  <- calc_se_shrk_alter_fixintcpt(K, L, J, Sigma.eps, D, C.b, var.rand.eff)
      
    } else if (intcpt == "random") {
      C.b <- c(0, 1)
      D   <- matrix(c(sigma.mu^2, sigma.mutau, sigma.mutau, sigma.tau^2), ncol = 2)
      se  <- calc_se_shrk_alter_randintcpt(K, L, J, Sigma.eps, D, C.b, var.rand.eff)
      
    }
    
  } else if (psbl.seq == "Pairwise Randomization") {
    
    if (intcpt == "fixed") {
      C.b <- 1
      D   <- sigma.tau^2
      se  <- calc_se_shrk_pair_fixintcpt(K, L, J, Sigma.eps, D, C.b, var.rand.eff)
      
    } else if (intcpt == "random") {
      C.b <- c(0, 1)
      D   <- matrix(c(sigma.mu^2, sigma.mutau, sigma.mutau, sigma.tau^2), ncol = 2)
      se  <- calc_se_shrk_pair_randintcpt(K, L, J, Sigma.eps, D, C.b, var.rand.eff)
    }
    
  }
  
  se <- matrix(se[!duplicated(se[, 3]), ], ncol = 3)
  
  return(se)
  
}


calc_se_shrk_genFormula <- function(C.theta, C.b, inv.sum.sigma.beta, D, Z.ij, X.ij, inv.Sigma.ij, var.rand.eff) {
  
  if (var.rand.eff) {
    tmp.se <- sqrt(t(C.theta) %*% inv.sum.sigma.beta %*% C.theta + 
                     t(C.b) %*% (D - D %*% t(Z.ij) %*% inv.Sigma.ij %*% Z.ij %*% D + 
                                   D %*% t(Z.ij) %*% inv.Sigma.ij %*% X.ij %*% inv.sum.sigma.beta %*% t(X.ij) %*% inv.Sigma.ij %*% Z.ij %*% D) %*% C.b - 
                     2 * (t(C.theta) %*% inv.sum.sigma.beta %*% t(X.ij) %*% inv.Sigma.ij %*% Z.ij %*% D  %*% C.b))
  } else {
    tmp.se <- sqrt(t(C.theta) %*% inv.sum.sigma.beta %*% C.theta + 
                     t(C.b) %*% D %*% t(Z.ij) %*% 
                     (inv.Sigma.ij - inv.Sigma.ij %*% X.ij %*% inv.sum.sigma.beta %*% t(X.ij) %*% inv.Sigma.ij) %*%
                     Z.ij %*% D %*% C.b)
  }
  
  return(tmp.se)
}


calc_se_shrk_alter_fixintcpt <- function(K, L, J, Sigma.eps, D, C.b, var.rand.eff) {
  
  I <- 2
  C.theta <- c(rep(0, I * J), 1)
  
  sum.sigma.beta <- 0
  
  # se can be different based on different design matrices
  se <- NULL
  if ((K %% 2) == 0) {
    
    # For sum in the variance estimator
    for (j in 1:J){
      
      # 2 is for the alternating sequences
      X.ij.0      <- matrix(c(rep(0, K * L * J * I), rep(rep(c(0, 1), K / 2), each = L)), ncol = J * I + 1)
      X.ij.0[, j] <- 1
      
      Z.ij.0     <- matrix(rep(rep(c(0, 1), K / 2), each = L), ncol = 1)
      Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * I), rep(rep(c(1, 0), K / 2), each = L)), ncol = J * I + 1)
      X.ij.1[, j + J] <- 1
      
      Z.ij.1     <- matrix(rep(rep(c(1, 0), K / 2), each = L), ncol = 1)
      Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
      
      sum.sigma.beta <- sum.sigma.beta + 
        t(X.ij.0) %*% solve(Sigma.ij.0) %*% X.ij.0 +
        t(X.ij.1) %*% solve(Sigma.ij.1) %*% X.ij.1
      
    }
    
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    # Different pts in the trials have different design matrix
    # calculate se for each pt
    for (j in 1:J) {
      
      # sequence "01"
      X.ij.0      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(0, 1), K / 2), each = L)), ncol = J * 2 + 1)
      X.ij.0[, j] <- 1
      
      Z.ij.0     <- matrix(rep(rep(c(0, 1), K / 2), each = L), ncol = 1)
      Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
      
      inv.Sigma.ij.0 <- solve(Sigma.ij.0)
      
      tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                        Z.ij.0, X.ij.0, inv.Sigma.ij.0, 
                                        var.rand.eff)
      se <- rbind(se, c(j, "01", tmp.se))
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(1, 0), K / 2), each = L)), ncol = J * 2 + 1)
      X.ij.1[, j + J] <- 1
      
      Z.ij.1     <- matrix(rep(rep(c(1, 0), K / 2), each = L), ncol = 1)
      Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
      
      inv.Sigma.ij.1 <- solve(Sigma.ij.1)
      
      tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                        Z.ij.1, X.ij.1, inv.Sigma.ij.1, 
                                        var.rand.eff)
      se <- rbind(se, c(j, "10", tmp.se))
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
    
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    # Different pts in the trials have different design matrix
    # calculate se for each pt
    for (j in 1:J) {
      
      # sequence "01"
      X.ij.0      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = J * 2 + 1)
      X.ij.0[, j] <- 1
      
      Z.ij.0     <- matrix(c(rep(rep(c(0, 1), floor(K / 2)), each = L), rep(0, L)), ncol = 1)
      Sigma.ij.0 <- Z.ij.0 %*% D %*% t(Z.ij.0) + Sigma.eps
      
      inv.Sigma.ij.0 <- solve(Sigma.ij.0)
      
      tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                        Z.ij.0, X.ij.0, inv.Sigma.ij.0, 
                                        var.rand.eff)
      se <- rbind(se, c(j, "01", tmp.se))
      
      X.ij.1      <- matrix(c(rep(0, K * L * J * 2), rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = J * 2 + 1)
      X.ij.1[, j + J] <- 1
      
      Z.ij.1     <- matrix(c(rep(rep(c(1, 0), floor(K / 2)), each = L), rep(1, L)), ncol = 1)
      Sigma.ij.1 <- Z.ij.1 %*% D %*% t(Z.ij.1) + Sigma.eps
      
      inv.Sigma.ij.1 <- solve(Sigma.ij.1)
      
      tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                        Z.ij.1, X.ij.1, inv.Sigma.ij.1, 
                                        var.rand.eff)
      se <- rbind(se, c(j, "10", tmp.se))
    } # for j loop
    
  } # else ((K %% 2) != 0)
  
  return(se)
  
}



calc_se_shrk_alter_randintcpt <- function(K, L, J, Sigma.eps, D, C.b, var.rand.eff) {
  
  I <- 2
  C.theta <- c(0, 1)
  
  sum.sigma.beta <- 0
  # se can be different based on different design matrices
  se <- NULL
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
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    # individual design matrices
    inv.Sigma.ij.0 <- solve(Sigma.ij.0)
    
    tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                      Z.ij.0, X.ij.0, inv.Sigma.ij.0, 
                                      var.rand.eff)
    se <- rbind(se, c(1, "01", tmp.se))
    
    inv.Sigma.ij.1 <- solve(Sigma.ij.1)
    
    tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                      Z.ij.1, X.ij.1, inv.Sigma.ij.1, 
                                      var.rand.eff)
    se <- rbind(se, c(1, "10", tmp.se))
    
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
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    # individual design matrices
    inv.Sigma.ij.0 <- solve(Sigma.ij.0)
    
    tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                      Z.ij.0, X.ij.0, inv.Sigma.ij.0, 
                                      var.rand.eff)
    se <- rbind(se, c(1, "01", tmp.se))
    
    inv.Sigma.ij.1 <- solve(Sigma.ij.1)
    
    tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                      Z.ij.1, X.ij.1, inv.Sigma.ij.1, 
                                      var.rand.eff)
    se <- rbind(se, c(1, "10", tmp.se))
    
  } # else ((K %% 2) != 0)
  
  return(se)
}



calc_se_shrk_pair_fixintcpt <- function(K, L, J, Sigma.eps, D, C.b, var.rand.eff) {
  
  I       <- 2^(ceiling(K / 2))
  C.theta <- c(rep(0, I * J), 1)
  
  sum.sigma.beta <- 0
  
  # se can be different based on different design matrices
  se <- NULL
  if ((K %% 2) == 0) { 
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), K/2)
    seq.pair <- expand.grid(seq.pair)
    
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
    
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    # Different pts in the trials have different design matrix
    # calculate se for each pt
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
        
        inv.Sigma.ij <- solve(Sigma.ij)
        
        tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                          Z.ij, X.ij, inv.Sigma.ij, 
                                          var.rand.eff)
        se <- rbind(se, 
                    c(j, paste(seq.pair[i, ], collapse = ""), tmp.se))
        
      } # for j loop
      
    } # for i loop
    
  } else { # if ((K %% 2) == 0) 
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), (K + 1)/2)
    seq.pair <- expand.grid(seq.pair)
    
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
    
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    # Different pts in the trials have different design matrix
    # calculate se for each pt
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
        
        inv.Sigma.ij <- solve(Sigma.ij)
        
        tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                          Z.ij, X.ij, inv.Sigma.ij, 
                                          var.rand.eff)
        # tmp.term.1 <- t(C.theta) %*% inv.sum.sigma.beta %*% C.theta
        # tmp.term.2 <- t(C.b) %*% D %*% t(Z.ij) %*% inv.Sigma.ij %*% Z.ij %*% D %*% C.b
        # tmp.term.3 <- t(C.b) %*% D %*% t(Z.ij) %*% inv.Sigma.ij %*% X.ij %*% inv.sum.sigma.beta %*% t(X.ij) %*% inv.Sigma.ij %*% Z.ij %*% D %*% C.b
        # tmp.sigma.ij <- Sigma.ij
        # tmp.inv.sigma.ij <- inv.Sigma.ij
        
        se <- rbind(se, 
                    c(j, paste(seq.pair[i, ], collapse = ""), tmp.se))
        
      } # for j loop
      
    } # for i loop
    
  } # else ((K %% 2) == 0) 
  
  return(se)
  
}



calc_se_shrk_pair_randintcpt <- function(K, L, J, Sigma.eps, D, C.b, var.rand.eff) {
  
  I       <- 2^(ceiling(K / 2))
  C.theta <- c(0, 1)
  
  sum.sigma.beta <- 0
  # se can be different based on different design matrices
  se <- NULL
  if ((K %% 2) == 0) {
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), K/2)
    seq.pair <- expand.grid(seq.pair)
    
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
    
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq[seq(2, K, 2)] <- abs(trt.seq[seq(2, K, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
      Z.ij <- X.ij
      Sigma.ij <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
      
      # individual design matrices
      inv.Sigma.ij <- solve(Sigma.ij)
      tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                        Z.ij, X.ij, inv.Sigma.ij, 
                                        var.rand.eff)
      se <- rbind(se, 
                  c(1, paste(seq.pair[i, ], collapse = ""), tmp.se))
      
    } # for i loop
    
  } else { # ((K %% 2) != 0)
    
    # Generate the possible sequences through generating the first treatment period in each block
    seq.pair <- rep(list(c(0, 1)), (K + 1)/2)
    seq.pair <- expand.grid(seq.pair)
    
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
    
    inv.sum.sigma.beta <- solve(sum.sigma.beta)
    
    for (i in 1:I) {
      
      # Generate the treatment sequence in the design matrix
      trt.seq <- unlist(rep(seq.pair[i, ], each = 2))
      trt.seq <- trt.seq[-length(trt.seq)]
      trt.seq[seq(2, K-1, 2)] <- abs(trt.seq[seq(2, K-1, 2)] - 1)
      trt.seq <- rep(trt.seq, each = L)
      
      X.ij <- matrix(c(rep(1, K * L), trt.seq), ncol = 2)
      Z.ij <- X.ij
      Sigma.ij <- Z.ij %*% D %*% t(Z.ij) + Sigma.eps
      
      # individual design matrices
      inv.Sigma.ij <- solve(Sigma.ij)
      tmp.se <- calc_se_shrk_genFormula(C.theta, C.b, inv.sum.sigma.beta, D, 
                                        Z.ij, X.ij, inv.Sigma.ij, 
                                        var.rand.eff)
      se <- rbind(se, 
                  c(1, paste(seq.pair[i, ], collapse = ""), tmp.se))
    } # for i loop
    
  } # else ((K %% 2) != 0)
  
  return(se)
  
}
