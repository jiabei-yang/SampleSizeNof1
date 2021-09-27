gen_corr_str <- function(K, L, corstr, corr.resid.err) {
  
  if (corstr == "Independent"){
    
    V <- diag(K * L)
    
  } else if (corstr == "Exchangeable"){
    
    V <- matrix(corr.resid.err, nrow = K * L, ncol = K * L)
    diag(V) <- 1
    
  } else if (corstr == "AR-1") { # AR-1
    
    times <- 1:(K * L)
    H <- abs(outer(times, times, "-"))
    V <- corr.resid.err^H
    
  }
  return(V)
  
}

# calculate I
# May not need this function
calc_I_popavg <- function(K, psbl.seq) {

  if (psbl.seq == "Alternating Sequences") {

    I = 2

  } else if (psbl.seq == "Pairwise Randomization") {

    if ((K %% 2) == 0) {
      I <- 2^(K/2)
    } else {
      I <- 2^((K + 1)/2)
    }
  }

  return(I)
}

optimize_dsgn <- function(data.fig) {
  
  # Input data.frame with columns K, L, I and CorrStr
  # This function finds the optimized K, L, I combinations
  # Since the data is generated from different combinations of K and L
  # need to take care of K, I and L, I combination
  
  # for the same correlation structure, same L and I, retain the minimum K
  data.fig.opt <- data.fig %>%
    arrange(intcpt, slp, J, K, L)
  data.fig.opt <- data.fig.opt[!duplicated(data.fig.opt[, c("intcpt", "slp", "J", "K")]), ]
  
  # for the same correlation structure, same K and I, retain the minimum L
  data.fig.opt <- data.fig.opt %>%
    arrange(intcpt, slp, J, L, K)
  data.fig.opt <- data.fig.opt[!duplicated(data.fig.opt[, c("intcpt", "slp", "J", "L")]), ]
  
  data.fig.opt <- data.fig.opt %>%
    arrange(intcpt, slp, K, L, J)
  data.fig.opt <- data.fig.opt[!duplicated(data.fig.opt[, c("intcpt", "slp", "K", "L")]), ]
  
  # data.fig.opt <- data.fig.opt %>%
  #   arrange(intcpt, slp, K, L, J)
  
  return(data.fig.opt)
}



# Calculate the values for producing the figure
gen_data_fig_ppobs <- function(data.fig.opt) {
  
  data.fig.opt <- data.fig.opt %>%
    mutate(KL = K * L)  %>%
    mutate(IJ   = I * J) %>%
    mutate(IJKL = IJ * KL)
  
  # data.fig.opt <- data.fig.opt %>%
  #   mutate(CorrStr = fct_explicit_na(CorrStr))
  
  summ.data.fig.opt <- data.fig.opt %>% 
    group_by(intcpt, slp, KL) %>%
    summarise(min.IJKL  = min(IJKL),
              mean.IJKL = mean(IJKL),
              max.IJKL  = max(IJKL),
              min.IJ    = min(IJ),
              mean.IJ   = mean(IJ),
              max.IJ    = max(IJ))
  
  summ.data.fig.opt <- summ.data.fig.opt %>%
    left_join(data.fig.opt, by = c("intcpt" = "intcpt", "slp" = "slp", "KL" = "KL"))
  
  summ.data.fig.opt <- data.frame(summ.data.fig.opt)
  
  return(summ.data.fig.opt)
  
}



gen_data_fig_pts <- function(data.fig.opt) {
  
    data.fig.opt <- data.fig.opt %>%
      mutate(KL = K * L)  %>%
      mutate(IJ   = I * J) %>%
      mutate(IJKL = IJ * KL)
    # data.fig.opt <- data.fig.opt %>%
    #   mutate(CorrStr = fct_explicit_na(CorrStr))
    
    summ.data.fig.opt <- data.fig.opt %>% 
      group_by(intcpt, slp, IJ) %>%
      summarise(min.IJKL  = min(IJKL),
                mean.IJKL = mean(IJKL),
                max.IJKL  = max(IJKL),
                min.KL    = min(KL),
                mean.KL   = mean(KL),
                max.KL    = max(KL))
    
    summ.data.fig.opt <- summ.data.fig.opt %>%
      left_join(data.fig.opt, by = c("intcpt" = "intcpt", "slp" = "slp", "IJ" = "IJ"))
    
    summ.data.fig.opt <- data.frame(summ.data.fig.opt)
    
  return(summ.data.fig.opt)
  
}

# gen.chr.seq <- function(seq, psbl.seq, K) {
#   
#   seq <- str_split(seq, "", simplify = T)
#   seq <- as.numeric(seq)
#   
#   if (psbl.seq == "Alternating Sequences") {
#     
#     seq <- rep(seq, ceiling(K/2))
#   } else {
#     
#     seq <- rep(seq, each = 2)
#     
#     if ((K %% 2) == 0) {
#       seq[seq(2, K, 2)] <- abs(seq[seq(2, K, 2)] - 1)
#     } else {
#       seq[seq(2, K - 1, 2)] <- abs(seq[seq(2, K - 1, 2)] - 1)
#     }
#     
#   }
#  
#   seq <- seq[1:K]
#   seq <- paste(seq, collapse = "")
#   
#   return(seq)
# }
