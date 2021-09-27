#!/usr/bin/env Rscript
library(optparse)
library(doParallel)

registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# Read in the arguments from command line
option_list <- list(make_option(c("-q", "--psblseq"), type="character"),
                    make_option(c("-s", "--corrstr"), type="character"),
                    make_option(c("-a", "--alpha"), type="double"),
                    make_option(c("-b", "--beta"), type="double"),
                    make_option(c("-e", "--sigmaepsilonsq"), type="double"),
                    make_option(c("-r", "--rho"), type="double"),
                    make_option(c("-m", "--sigmamusq"), type="double"),
                    make_option(c("-t", "--sigmatausq"), type = "double"),
                    make_option(c("-c", "--sigmamutau"), type = "double"),
                    make_option(c("-d", "--mindiff"), type="double"),
                    make_option(c("-f", "--filename"), type="character"))
# a, b, c, d, e, f, m, q, r, s, t
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)

# functions
source("../Functions/gen.R")
source("../Functions/Power.R")

intcpt  <- c("fixed", "random")
slp     <- c("common", "random")

scnrs <- expand.grid(slp, intcpt)
scnrs           <- scnrs[, 2:1]
colnames(scnrs) <- c("intcpt", "slp")

# parameters
range.ppobs     <- c(1, 300)

CorrStr         <- opt$corrstr
corr            <- opt$rho
sigma.resid.err <- sqrt(opt$sigmaepsilonsq)
type.1.err      <- opt$alpha
power           <- 1 - opt$beta
min.diff        <- opt$mindiff
sigma.tau       <- sqrt(opt$sigmatausq)
sigma.mu        <- sqrt(opt$sigmamusq)
sigma.mutau     <- opt$sigmamutau
psbl.seq        <- opt$psblseq

# CorrStr         <- "AR-1"
# corr            <- 0.4
# sigma.resid.err <- 2
# type.1.err      <- 0.05
# power           <- 0.8
# min.diff        <- 1
# sigma.tau       <- 1
# sigma.mu        <- 2
# sigma.mutau     <- 1
# psbl.seq        <- "Pairwise Randomization"

data.fig <- foreach(scnrs.i = 1:nrow(scnrs)) %dopar% {
  
  tmp.data.fig <- NULL
  
  # Record the previous min.I for last k to save computation time
  prev.max.J <- NULL
  
  # At least 2 periods per participant
  for (k in 2:range.ppobs[2]) {
    # for (k in 4:range.ppobs[2]) {
    
    if ((k %% 10) == 0) {
      print(k) 
    }
    
    I <- ifelse((k %% 2) == 0, 2^(k/2), 2^((k+1)/2))
    
    # Find earlier J values, new J should be equal or smaller since we increase k and l
    psbl.l       <- ceiling(range.ppobs[1] / k):floor(range.ppobs[2] / k)
    max.J        <- rep(NA, length(psbl.l))
    names(max.J) <- as.character(psbl.l)
    
    if (!is.null(prev.max.J)) {
      max.J[names(max.J) %in% names(prev.max.J)] <- prev.max.J[names(prev.max.J) %in% names(max.J)]
    }
    
    # IF psbl.l contains only 1 and max.J is 1, then can skip the following loops
    lth.psbl.l <- length(psbl.l)
    
    # break if the only one element in psbl.l is 1 and max.J for this element is also 1
    # this will not be optimal so nothing saved to data.fig
    # no use if check at the end of the loop
    # if ((lth.psbl.l == 1) & (psbl.l[1] == 1)) {
    #   if (max.J == 1) {
    #     break
    #   }
    # }
    
    for (l in psbl.l) {
      # for (l in 2:100) {
      
      TmpCorrMtrx <- gen_corr_str(K      = k, 
                                  L      = l, 
                                  corstr = CorrStr, 
                                  rho    = corr)
      max.J.ind <- which(names(max.J) == as.character(l))
      
      # When equal to na and not equal to 1, use function to calculate;
      # otherwise assign to 1 directly
      if (!is.na(max.J[max.J.ind])) {
        if (max.J[max.J.ind] == 1) {
          
          TmpDsgn <- c(k, l, 1)
          
        } else {
          
          TmpDsgn   <- calc_j_popavg(K               = k, 
                                     L               = l, 
                                     max.J           = max.J[max.J.ind], 
                                     psbl.seq        = psbl.seq,
                                     intcpt          = scnrs$intcpt[scnrs.i],
                                     slp             = scnrs$slp[scnrs.i],
                                     corr.mtrx       = TmpCorrMtrx, 
                                     sigma.resid.err = sigma.resid.err,
                                     sigma.mu        = sigma.mu,
                                     sigma.tau       = sigma.tau,
                                     sigma.mutau     = sigma.mutau,
                                     type.1.err      = type.1.err, 
                                     power           = power, 
                                     min.diff        = min.diff)
        }
      } else {
        
        # tmp <- max.J
        # K               = k
        # L               = l
        # max.J           = max.J[max.J.ind]
        # psbl.seq        = psbl.seq
        # intcpt          = scnrs$intcpt[scnrs.i]
        # slp             = scnrs$slp[scnrs.i]
        # corr.mtrx       = TmpCorrMtrx
        # sigma.resid.err = sigma.resid.err
        # sigma.mu        = sigma.mu
        # sigma.tau       = sigma.tau
        # sigma.mutau     = sigma.mutau
        # type.1.err      = type.1.err
        # power           = power
        # min.diff        = min.diff
        # max.J <- tmp
        
        TmpDsgn   <- calc_j_popavg(K               = k, 
                                   L               = l, 
                                   max.J           = max.J[max.J.ind], 
                                   psbl.seq        = psbl.seq,
                                   intcpt          = scnrs$intcpt[scnrs.i],
                                   slp             = scnrs$slp[scnrs.i],
                                   corr.mtrx       = TmpCorrMtrx, 
                                   sigma.resid.err = sigma.resid.err,
                                   sigma.mu        = sigma.mu,
                                   sigma.tau       = sigma.tau,
                                   sigma.mutau     = sigma.mutau,
                                   type.1.err      = type.1.err, 
                                   power           = power, 
                                   min.diff        = min.diff)
      }
      
      max.J[max.J.ind]     <- TmpDsgn[3]
      
      # If l the last one, then cannot assign to the next value
      if (l != psbl.l[lth.psbl.l]) {
        
        # need to see if the next value is NA
        if (is.na(max.J[max.J.ind + 1])) {
          max.J[max.J.ind + 1] <- TmpDsgn[3]
        } else {
          
          # if the next value is not NA, only assign if the value is greater than the current value
          if (max.J[max.J.ind + 1] > TmpDsgn[3]) {
            max.J[max.J.ind + 1] <- TmpDsgn[3]
          }
        }
      }
      
      tmp.data.fig <- rbind(tmp.data.fig, c(TmpDsgn, I,
                                            as.character(scnrs$intcpt[scnrs.i]), 
                                            as.character(scnrs$slp[scnrs.i]),
                                            CorrStr))
      
      # if I == 1, no need to run more since for larger K and L, I will be 1 as well
      if (TmpDsgn[3] == 1) {
        
        # can break directly since we will delete the later observations anyway
        max.J[max.J.ind:lth.psbl.l] <- 1
        break
        
      }
      
    } # for l loop
    
    prev.max.J <- max.J
    
    # if both L and J are 1, greater K will lead to L and J being 1 as well
    # therefore break k loop
    nrow.data.fig <- nrow(tmp.data.fig)
    if ((tmp.data.fig[nrow.data.fig, 2] == 1) & (tmp.data.fig[nrow.data.fig, 3] == 1)) {
      break
    }
    
  } # for k loop
  
  print(scnrs[scnrs.i, ])
  
  tmp.data.fig                      
}

save(data.fig, file = paste0("../Data/", opt$filename))