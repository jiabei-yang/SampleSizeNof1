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
gen.max.KL <- 365 * 3
range.pts     <- c(2, 150)

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
  
  j <- 1
  k <- 2
  
  I <- ifelse((k %% 2) == 0, 2^(k/2), 2^((k+1)/2))
  
  # first generate prev.max.L for j = 1
  # will be used for j = 2, 3, ...
  TmpDsgn   <- calc_l_popavg(J               = j,
                             K               = k, 
                             max.L           = NA, 
                             gen.max.KL      = gen.max.KL, 
                             psbl.seq        = psbl.seq,
                             intcpt          = scnrs$intcpt[scnrs.i],
                             slp             = scnrs$slp[scnrs.i],
                             sigma.resid.err = sigma.resid.err,
                             sigma.mu        = sigma.mu,
                             sigma.tau       = sigma.tau,
                             sigma.mutau     = sigma.mutau,
                             corstr          = CorrStr, 
                             rho             = corr,
                             type.1.err      = type.1.err, 
                             power           = power, 
                             min.diff        = min.diff)
  
  prev.max.L <- TmpDsgn[2]
  names(prev.max.L) <- k
  
  tmp.data.fig <- rbind(tmp.data.fig, c(TmpDsgn, I, 
                                        as.character(scnrs$intcpt[scnrs.i]), 
                                        as.character(scnrs$slp[scnrs.i]),
                                        CorrStr))
  
  last.L <- TmpDsgn[2]
  k <- k + 1
  I <- ifelse((k %% 2) == 0, 2^(k/2), 2^((k+1)/2))
  
  while (((last.L > 1) | is.na(last.L)) & 
         (k <= gen.max.KL) & 
         (I <= (range.pts[2] / j))) {
    
    # if the last.L * k goes out of gen.max.KL, last.L go back to NA
    if (!is.na(last.L)) {
      if (last.L > floor(gen.max.KL / k)) {
        last.L <- NA
      }
    }
    
    TmpDsgn   <- calc_l_popavg(J               = j,
                               K               = k, 
                               max.L           = last.L, 
                               gen.max.KL      = gen.max.KL, 
                               psbl.seq        = psbl.seq,
                               intcpt          = scnrs$intcpt[scnrs.i],
                               slp             = scnrs$slp[scnrs.i],
                               sigma.resid.err = sigma.resid.err,
                               sigma.mu        = sigma.mu,
                               sigma.tau       = sigma.tau,
                               sigma.mutau     = sigma.mutau,
                               corstr          = CorrStr, 
                               rho             = corr,
                               type.1.err      = type.1.err, 
                               power           = power, 
                               min.diff        = min.diff)
    
    prev.max.L <- c(prev.max.L, TmpDsgn[2])
    names(prev.max.L)[length(prev.max.L)] <- k
    
    last.L     <- TmpDsgn[2]
    
    tmp.data.fig <- rbind(tmp.data.fig, c(TmpDsgn, I, 
                                          as.character(scnrs$intcpt[scnrs.i]), 
                                          as.character(scnrs$slp[scnrs.i]),
                                          CorrStr))
    
    k <- k + 1
    I <- ifelse((k %% 2) == 0, 2^(k/2), 2^((k+1)/2))
  }
  
  # The maximum number of J should be range.pts/2
  # since the minimum number of different sequences is 2
  for (j in 2:(range.pts[2] / 2)) {
    # for (j in 8:(range.pts[2] / 2)) {
    
    if ((j %% 5) == 0) {
      print(j) 
    }
    
    max.L     <- prev.max.L
    lth.max.L <- length(max.L)
    
    # prepare while loop
    k         <- 2
    max.L.ind <- which(names(max.L) == k)
    I         <- 2
    
    # stop iterating if 1) last.L drops to 1, 2) k goes out of range and 3) # of participants goes out of range
    while (((max.L[max.L.ind] > 1) | is.na(max.L[max.L.ind])) & 
           (k <= min(gen.max.KL, as.numeric(names(max.L)[lth.max.L]))) &
           (I <= (range.pts[2] / j))) {
      
      max.L.ind <- which(names(max.L) == k)
      last.L <- max.L[max.L.ind]
      
      TmpDsgn   <- calc_l_popavg(J               = j,
                                 K               = k, 
                                 max.L           = last.L, 
                                 gen.max.KL      = gen.max.KL, 
                                 psbl.seq        = psbl.seq,
                                 intcpt          = scnrs$intcpt[scnrs.i],
                                 slp             = scnrs$slp[scnrs.i],
                                 sigma.resid.err = sigma.resid.err,
                                 sigma.mu        = sigma.mu,
                                 sigma.tau       = sigma.tau,
                                 sigma.mutau     = sigma.mutau,
                                 corstr          = CorrStr, 
                                 rho             = corr,
                                 type.1.err      = type.1.err, 
                                 power           = power, 
                                 min.diff        = min.diff)
      
      max.L[max.L.ind] <- TmpDsgn[2]
      
      # If k the last one, then cannot assign to the next value
      # if the current value is not NA
      if ((max.L.ind != lth.max.L) & (!is.na(max.L[max.L.ind]))) {
        
        # assign if the next value is NA
        if (is.na(max.L[max.L.ind + 1])) {
          
          # if the L * (k+1) goes out of gen.max.KL, next L go back to NA
          # in other words, only assign when L * (k+1) <= gen.max.KL
          if (TmpDsgn[2] <= floor(gen.max.KL / as.numeric(names(max.L)[max.L.ind + 1]))) {
            max.L[max.L.ind + 1] <- TmpDsgn[2]
          }
          
        } else {
          # if the next value is not NA, only assign if the value is greater than the current value
          if (max.L[max.L.ind + 1] > TmpDsgn[2]) {
            max.L[max.L.ind + 1] <- TmpDsgn[2]
          }
        }
        
      }
      
      tmp.data.fig <- rbind(tmp.data.fig, c(TmpDsgn, I, 
                                            as.character(scnrs$intcpt[scnrs.i]), 
                                            as.character(scnrs$slp[scnrs.i]),
                                            CorrStr))
      
      k <- k + 1
      I <- ifelse((k %% 2) == 0, 2^(k/2), 2^((k+1)/2))
      
    } # while k loop
    
    # do not need if else because the result is the same
    # if (length(which(max.L == 1)) > 1) {
    prev.max.L <- max.L[1:(which(max.L == 1)[1])]
    # } else {
    #   prev.max.L <- max.L
    # }
    
    # DO WE NEED TO DECIDE WHETHER THE REMAINING ONE ELEMENT IS 1?
    # No, since always shrinking length
    if (length(prev.max.L) == 1) {
      break
    }
    
  } # for j loop
  
  print(scnrs[scnrs.i, ])
  
  tmp.data.fig
}

save(data.fig, file = paste0("../Data/", opt$filename))


