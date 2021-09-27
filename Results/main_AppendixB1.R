library(dplyr)
library(ggplot2)
library(gridExtra)

source("../Functions/gen.R")

# ppobs
load("../Data/main/mainPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
# identical(data.fig[[1]], data.fig[[3]])
# identical(data.fig[[2]], data.fig[[4]])
# 
# which(data.fig[[1]][314] != data.fig[[3]][314])

summ.data.fig.ppobs.opt <- prep_data_fig(data.fig,
                                         x.axis = "KL")

# Pts
load("../Data/main/mainPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.data.fig.pts.opt <- prep_data_fig(data.fig,
                                       x.axis = "IJ")

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

p.main.ppobs <- ggplot(summ.data.fig.ppobs.opt %>%
               filter((KL < 75) & (intcpt == "Fixed Intercepts")), aes(x = KL, y = mean.IJKL, group = slp,
                                      text = paste("Number of measurements per participant: ", KL,
                                                   "<br>Average number of measurements across trials: ", round(mean.IJKL, 2)))) +
  geom_point(aes(color = slp, shape = slp), size = 1.5) + 
  geom_line(aes(color = slp, linetype = slp)) + 
  geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = slp), alpha = 0.3) +
  # scale_y_log10() +
  # facet_grid(slp ~ intcpt) +
  scale_y_continuous(limits = c(70, 1400),
                     breaks = seq(250, 1250, 500)) +
  scale_linetype_manual(values = c("dashed", "dotdash", "dotted", "longdash")) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  xlab("Number of measurements per participant") +
  ylab("Total number of measurements across trials") +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill     = guide_legend(ncol = 2, byrow = T, title = "Slope"),
         shape    = guide_legend(ncol = 2, byrow = T, title = "Slope"),
         linetype = guide_legend(ncol = 2, byrow = T, title = "Slope"),
         color    = guide_legend(ncol = 2, byrow = T, title = "Slope"))

p.main.pts <- ggplot(summ.data.fig.pts.opt %>%
               filter((intcpt == "Fixed Intercepts")), aes(x = IJ, y = mean.IJKL, group = slp)) +
  geom_point(aes(color = slp, shape = slp), size = 1.5) + 
  geom_line(aes(color = slp, linetype = slp)) + 
  geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = slp), alpha = 0.3) +
  scale_y_continuous(limits = c(70, 1400),
                     breaks = seq(250, 1250, 500)) +
  # facet_grid(slp ~ intcpt) +
  scale_linetype_manual(values = c("dashed", "dotdash", "dotted", "longdash")) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  xlab("Number of participants") +
  ylab("Total number of measurements across trials") +
  # ggtitle(paste0("Optimized designs when pairwise randomization is used")) +
  # labs(color    = "Correlation structure",
  #      shape    = "Correlation structure",
  #      fill     = "Correlation structure",
  #      linetype = "Correlation structure") +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill     = guide_legend(ncol = 2, byrow = T, title = "Slope"),
         shape    = guide_legend(ncol = 2, byrow = T, title = "Slope"),
         linetype = guide_legend(ncol = 2, byrow = T, title = "Slope"),
         color    = guide_legend(ncol = 2, byrow = T, title = "Slope"))

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend.shared.main <- get_legend(p.main.ppobs)
p.main.ppobs <- p.main.ppobs + theme(legend.position="none")
p.main.pts   <- p.main.pts + theme(legend.position="none")

# Figure 2
grid.arrange(p.main.ppobs, p.main.pts, legend.shared.main, ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1, 2), rep(3, 2)),
             widths = c(1/2, 1/2), heights = c(9/10, 1/10))

# Figure S1
p.app.ppobs <- ggplot(summ.data.fig.ppobs.opt %>%
               filter((KL < 75)), aes(x = KL, y = mean.IJKL, group = Model,
                                      text = paste("Number of measurements per participant: ", KL,
                                                   "<br>Average number of measurements across trials: ", round(mean.IJKL, 2)))) +
  geom_point(aes(color = Model, shape = Model), size = 1.5) + 
  geom_line(aes(color = Model, linetype = Model)) + 
  geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = Model), alpha = 0.3) +
  # scale_y_log10() +
  # facet_grid(slp ~ intcpt) +
  scale_y_continuous(limits = c(70, 1400),
                     breaks = seq(250, 1250, 500)) +
  scale_linetype_manual(values = c("dashed", "dotdash", "dotted", "longdash")) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  xlab("Number of measurements per participant") +
  ylab("Total number of measurements across trials") +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill     = guide_legend(ncol = 2, byrow = T),
         shape    = guide_legend(ncol = 2, byrow = T),
         linetype = guide_legend(ncol = 2, byrow = T),
         color    = guide_legend(ncol = 2, byrow = T))

p.app.pts <- ggplot(summ.data.fig.pts.opt, aes(x = IJ, y = mean.IJKL, group = Model)) +
  geom_point(aes(color = Model, shape = Model), size = 1.5) + 
  geom_line(aes(color = Model, linetype = Model)) + 
  geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = Model), alpha = 0.3) +
  scale_y_continuous(limits = c(70, 1400),
                     breaks = seq(250, 1250, 500)) +
  # facet_grid(slp ~ intcpt) +
  scale_linetype_manual(values = c("dashed", "dotdash", "dotted", "longdash")) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  xlab("Number of participants") +
  ylab("Total number of measurements across trials") +
  # ggtitle(paste0("Optimized designs when pairwise randomization is used")) +
  # labs(color    = "Correlation structure",
  #      shape    = "Correlation structure",
  #      fill     = "Correlation structure",
  #      linetype = "Correlation structure") +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  guides(fill     = guide_legend(ncol = 2, byrow = T),
         shape    = guide_legend(ncol = 2, byrow = T),
         linetype = guide_legend(ncol = 2, byrow = T),
         color    = guide_legend(ncol = 2, byrow = T))

legend.shared.app <- get_legend(p.app.ppobs)
p.app.ppobs <- p.app.ppobs + theme(legend.position="none")
p.app.pts   <- p.app.pts + theme(legend.position="none")

grid.arrange(p.app.ppobs, p.app.pts, legend.shared.app, ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1, 2), rep(3, 2)),
             widths = c(1/2, 1/2), heights = c(8.5/10, 1.5/10))

