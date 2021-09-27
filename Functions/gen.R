cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



gen_corr_str <- function(K, L, corstr, rho) {
  
  if (corstr == "Independent"){
    
    V <- diag(K * L)
    
  } else if (corstr == "Exchangeable"){
    
    V <- matrix(rho, nrow = K * L, ncol = K * L)
    diag(V) <- 1
    
  } else if (corstr == "AR-1") { # AR-1
    
    times <- 1:(K * L)
    H <- abs(outer(times, times, "-"))
    V <- rho^H
    
  }
  return(V)
  
}



prep_data_fig <- function(data.fig, x.axis) {
  
  data.fig <- data.frame(rbind(data.fig[[1]], 
                               data.fig[[2]],
                               data.fig[[3]],
                               data.fig[[4]]))
  colnames(data.fig) <- c("K", "L", "J", "I", "intcpt", "slp", "CorrStr")
  data.fig <- data.fig %>%
    filter(!is.na(L))
  data.fig <- data.fig %>%
    mutate(K = as.integer(levels(K))[K],
           L = as.integer(levels(L))[L],
           J = as.integer(levels(J))[J],
           I = as.integer(levels(I))[I])
  data.fig <- optimize_dsgn(data.fig)
  
  if(x.axis == "KL") {
    data.fig <- gen_data_fig_ppobs(data.fig)
  } else { # x.axis == "IJ"
    data.fig <- gen_data_fig_pts(data.fig)
  }
  
  data.fig <- data.fig %>%
    mutate(intcpt = factor(intcpt, 
                           levels = c("fixed", "random"),
                           labels = c("Fixed Intercepts", "Random Intercepts")),
           slp    = factor(slp,
                           levels = c("common", "random"),
                           labels = c("Common Slope", "Random Slopes")))
  
  data.fig <- data.fig %>%
    mutate(Model = paste0(intcpt, "-", slp)) %>%
    mutate(Model = factor(Model,
                          levels = c("Fixed Intercepts-Common Slope", "Random Intercepts-Common Slope", 
                                     "Fixed Intercepts-Random Slopes", "Random Intercepts-Random Slopes")))
  
  return(data.fig)
}





optimize_dsgn <- function(data.fig) {
  
  # Input data.frame with columns K, L, J, I, intcpt, slp and CorrStr
  # This function finds the optimized K, L, J combinations
  # Since for a given K, I will be the same
  # need to take care of K, J and L, J combination
  
  # for the same correlation structure, same L and I, retain the minimum K
  data.fig.opt <- data.fig %>%
    arrange(intcpt, slp, CorrStr, L, J, K)
  data.fig.opt <- data.fig.opt[!duplicated(data.fig.opt[, c("L", "J", "intcpt", "slp", "CorrStr")]), ]
  
  # for the same correlation structure, same K and I, retain the minimum L
  data.fig.opt <- data.fig.opt %>%
    arrange(intcpt, slp, CorrStr, K, J, L)
  data.fig.opt <- data.fig.opt[!duplicated(data.fig.opt[, c("K", "J", "intcpt", "slp", "CorrStr")]), ]
  
  # for the same correlation structure, same K and L, retain the minimum I
  data.fig.opt <- data.fig.opt %>%
    arrange(intcpt, slp, CorrStr, K, L, J)
  data.fig.opt <- data.fig.opt[!duplicated(data.fig.opt[, c("K", "L", "intcpt", "slp", "CorrStr")]), ]
  
  data.fig.opt <- data.fig.opt %>%
    arrange(intcpt, slp, CorrStr, K, L, J)
  
  return(data.fig.opt)
}



# for figure with x-axis being KL
gen_data_fig_ppobs <- function(data.fig.opt) {
  
  data.fig.opt <- data.fig.opt %>%
    mutate(KL   = K * L)  %>%
    mutate(IJ   = I * J) %>%
    mutate(IJKL = IJ * KL)
  
  summ.data.fig.opt <- data.fig.opt %>% 
    group_by(intcpt, slp, CorrStr, KL) %>%
    summarise(min.IJKL  = min(IJKL),
              mean.IJKL = mean(IJKL),
              max.IJKL  = max(IJKL),
              min.IJ    = min(IJ),
              mean.IJ   = mean(IJ),
              max.IJ    = max(IJ))
  
  summ.data.fig.opt <- summ.data.fig.opt %>%
    left_join(data.fig.opt, by = c("intcpt" = "intcpt", "slp" = "slp", "CorrStr" = "CorrStr", "KL" = "KL"))
  
  summ.data.fig.opt <- data.frame(summ.data.fig.opt)
  
  return(summ.data.fig.opt)
  
}



# for figure with x-axis being IJ
gen_data_fig_pts<- function(data.fig.opt) {
  
  data.fig.opt <- data.fig.opt %>%
    mutate(KL   = K * L)  %>%
    mutate(IJ   = I * J) %>%
    mutate(IJKL = IJ * KL)
  
  summ.data.fig.opt <- data.fig.opt %>% 
    group_by(intcpt, slp, CorrStr, IJ) %>%
    summarise(min.IJKL  = min(IJKL),
              mean.IJKL = mean(IJKL),
              max.IJKL  = max(IJKL),
              min.KL    = min(KL),
              mean.KL   = mean(KL),
              max.KL    = max(KL))
  
  summ.data.fig.opt <- summ.data.fig.opt %>%
    left_join(data.fig.opt, by = c("intcpt" = "intcpt", "slp" = "slp", "CorrStr" = "CorrStr", "IJ" = "IJ"))
  
  summ.data.fig.opt <- data.frame(summ.data.fig.opt)
  
  return(summ.data.fig.opt)
  
}



# default to 2*2 panels: facet + grid.arrange(nrow = 2)
plot.app.popavg.params <- function(data.fig, x.label, label.x = c(F, T), x.breaks, y.label, log.y,
                                   facet.x = T, legend.position = "right", nrow.panel = 2) {

  p1 <- ggplot(data.fig$KL, aes(x = params, y = mean.IJKL, group = KL)) + 
    geom_point(aes(color = KL, shape = KL), size = 2) +
    geom_line(aes(color = KL, linetype = KL)) +
    geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = KL), alpha = 0.3) 
  p1 <- plot.app.popavg.params.sep(p = p1,
                                   legend.title = "KL",
                                   x.breaks     = x.breaks,
                                   if.x.label   = label.x[1], 
                                   x.label      = x.label[1],
                                   log.y        = log.y[1],
                                   facet.x      = facet.x,
                                   legend.position = legend.position)
 
  p2 <- ggplot(data.fig$IJ, aes(x = params, y = mean.IJKL, group = IJ)) + 
    geom_point(aes(color = IJ, shape = IJ), size = 2) +
    geom_line(aes(color = IJ, linetype = IJ)) +
    geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = IJ), alpha = 0.3)
  p2 <- plot.app.popavg.params.sep(p = p2,
                                   legend.title = "IJ",
                                   x.breaks     = x.breaks,
                                   if.x.label   = label.x[2], 
                                   x.label      = x.label[2],
                                   log.y        = log.y[2],
                                   facet.x      = facet.x,
                                   legend.position = legend.position)
  
  # legend.shared <- get_legend(p1)
  # p1 <- p1 + theme(legend.position="none")
  # p2 <- p2 + theme(legend.position="none")

  grid.arrange(p1, p2, nrow = nrow.panel,
               left = y.label)

}



plot.app.popavg.params.sep <- function(p, legend.title, x.breaks, if.x.label, x.label, log.y, 
                                       facet.x, legend.position) {
  
  p <- p + 
    scale_x_continuous(breaks = x.breaks) +
    # scale_y_continuous(limits = y.limits,
    #                    breaks = y.breaks) +
    scale_linetype_manual(name = legend.title,
                          values = c("dashed", "dotdash", "dotted", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")) +
    scale_color_manual(name = legend.title,
                       values = cbPalette[c(7, 5, 3, 4, 1, 6, 2, 8)]) +
    scale_fill_manual(name = legend.title,
                      values = cbPalette[c(7, 5, 3, 4, 1, 6, 2, 8)]) +
    scale_shape_manual(name = legend.title,
                       values = 1:8) +
    theme_bw() 
  
  if (if.x.label) {
    p <- p +
      xlab(x.label) + 
      theme(legend.position  = legend.position,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = NA),
            axis.title.y     = element_blank()) 
  } else {
    p <- p +
      theme(legend.position  = legend.position,
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = NA),
            axis.title       = element_blank()) 
  }
  
  if (facet.x) {
    p <- p + facet_grid(. ~ slp)
  }
  
  if (log.y) {
    p <- p + scale_y_log10()
  }
  
  p
}
