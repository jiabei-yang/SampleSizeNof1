library(dplyr)
library(ggplot2)
library(gridExtra)

source("../Functions/gen.R")

# ppobs
load("../Data/main/mainPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.4.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e1_r04_m4_t1_c1_d1.RData")
summ.ppobs.1.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e2_r04_m4_t1_c1_d1.RData")
summ.ppobs.2.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e3_r04_m4_t1_c1_d1.RData")
summ.ppobs.3.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e5_r04_m4_t1_c1_d1.RData")
summ.ppobs.5.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e6_r04_m4_t1_c1_d1.RData")
summ.ppobs.6.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e7_r04_m4_t1_c1_d1.RData")
summ.ppobs.7.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e8_r04_m4_t1_c1_d1.RData")
summ.ppobs.8.04.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")

load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r01_m4_t1_c1_d1.RData")
summ.ppobs.4.01.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r02_m4_t1_c1_d1.RData")
summ.ppobs.4.02.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r03_m4_t1_c1_d1.RData")
summ.ppobs.4.03.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r05_m4_t1_c1_d1.RData")
summ.ppobs.4.05.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r06_m4_t1_c1_d1.RData")
summ.ppobs.4.06.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r07_m4_t1_c1_d1.RData")
summ.ppobs.4.07.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r08_m4_t1_c1_d1.RData")
summ.ppobs.4.08.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ar1_a005_b02_e4_r09_m4_t1_c1_d1.RData")
summ.ppobs.4.09.ar1 <- prep_data_fig(data.fig,
                                     x.axis = "KL")

load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Ind_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.4.04.ind <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB3ResidVar/ResidVarPpobs_Pair_Exch_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.4.04.exch <- prep_data_fig(data.fig,
                                      x.axis = "KL")

summ.ppobs.sgmEpsl <- rbind(summ.ppobs.4.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 4),
                            summ.ppobs.1.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 1),
                            summ.ppobs.2.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 2),
                            summ.ppobs.3.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 3),
                            summ.ppobs.5.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 5),
                            summ.ppobs.6.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 6),
                            summ.ppobs.7.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 7),
                            summ.ppobs.8.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = 8))
# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.ppobs.sgmEpsl %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.sgmEpsl <- summ.ppobs.sgmEpsl %>% 
  filter(KL %in% c(4, 6, 8, 12, 14))
summ.ppobs.sgmEpsl <- summ.ppobs.sgmEpsl %>% 
  mutate(KL = factor(KL, levels = c(4, 6, 8, 12, 14)))

summ.ppobs.rho <- rbind(summ.ppobs.4.04.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.4),
                        summ.ppobs.4.01.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.1),
                        summ.ppobs.4.02.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.2),
                        summ.ppobs.4.03.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.3),
                        summ.ppobs.4.05.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.5),
                        summ.ppobs.4.06.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.6),
                        summ.ppobs.4.07.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.7),
                        summ.ppobs.4.08.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.8),
                        summ.ppobs.4.09.ar1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.9))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.ppobs.rho %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.rho <- summ.ppobs.rho %>% 
  filter(KL %in% c(4, 6, 8, 10, 12))
summ.ppobs.rho <- summ.ppobs.rho %>% 
  mutate(KL = factor(KL, levels = c(4, 6, 8, 10, 12)))

summ.ppobs.corrStr <- rbind(summ.ppobs.4.04.ar1 %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = "AR-1"),
                            summ.ppobs.4.04.ind %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = "Independent"),
                            summ.ppobs.4.04.exch %>% 
                              filter(intcpt == "Fixed Intercepts") %>%
                              mutate(params = "Exchangeable"))
summ.ppobs.corrStr <- summ.ppobs.corrStr %>%
  mutate(params = factor(params, levels = c("AR-1", "Independent", "Exchangeable")))

# pts
load("../Data/main/mainPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.4.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e1_r04_m4_t1_c1_d1.RData")
summ.pts.1.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e2_r04_m4_t1_c1_d1.RData")
summ.pts.2.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e3_r04_m4_t1_c1_d1.RData")
summ.pts.3.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e5_r04_m4_t1_c1_d1.RData")
summ.pts.5.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e6_r04_m4_t1_c1_d1.RData")
summ.pts.6.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e7_r04_m4_t1_c1_d1.RData")
summ.pts.7.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e8_r04_m4_t1_c1_d1.RData")
summ.pts.8.04.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")

load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r01_m4_t1_c1_d1.RData")
summ.pts.4.01.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r02_m4_t1_c1_d1.RData")
summ.pts.4.02.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r03_m4_t1_c1_d1.RData")
summ.pts.4.03.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r05_m4_t1_c1_d1.RData")
summ.pts.4.05.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r06_m4_t1_c1_d1.RData")
summ.pts.4.06.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r07_m4_t1_c1_d1.RData")
summ.pts.4.07.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r08_m4_t1_c1_d1.RData")
summ.pts.4.08.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ar1_a005_b02_e4_r09_m4_t1_c1_d1.RData")
summ.pts.4.09.ar1 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")

load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Ind_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.4.04.ind <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB3ResidVar/ResidVarPts_Pair_Exch_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.4.04.exch <- prep_data_fig(data.fig,
                                    x.axis = "IJ")

summ.pts.sgmEpsl <- rbind(summ.pts.4.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 4),
                          summ.pts.1.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 1),
                          summ.pts.2.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 2),
                          summ.pts.3.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 3),
                          summ.pts.5.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 5),
                          summ.pts.6.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 6),
                          summ.pts.7.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 7),
                          summ.pts.8.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 8))
# still summarize by ij if plot epsilon on x-axis; 
tmp.1 <- summ.pts.sgmEpsl %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.sgmEpsl <- summ.pts.sgmEpsl %>% 
  filter(IJ %in% c(12, 14, 16, 20, 24))
summ.pts.sgmEpsl <- summ.pts.sgmEpsl %>% 
  mutate(IJ = factor(IJ, levels = c(12, 14, 16, 20, 24)))

summ.pts.rho <- rbind(summ.pts.4.04.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.4),
                      summ.pts.4.01.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.1),
                      summ.pts.4.02.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.2),
                      summ.pts.4.03.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.3),
                      summ.pts.4.05.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.5),
                      summ.pts.4.06.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.6),
                      summ.pts.4.07.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.7),
                      summ.pts.4.08.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.8),
                      summ.pts.4.09.ar1 %>% 
                        filter(intcpt == "Fixed Intercepts") %>%
                        mutate(params = 0.9))
# still summarize by ij if plot epsilon on x-axis; 
tmp.1 <- summ.pts.rho %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.rho <- summ.pts.rho %>% 
  filter(IJ %in% c(14, 16, 18, 20, 24))
summ.pts.rho <- summ.pts.rho %>% 
  mutate(IJ = factor(IJ, levels = c(14, 16, 18, 20, 24)))

summ.pts.corrStr <- rbind(summ.pts.4.04.ar1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = "AR-1"),
                          summ.pts.4.04.ind %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = "Independent"),
                          summ.pts.4.04.exch %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = "Exchangeable"))
summ.pts.corrStr <- summ.pts.corrStr %>%
  mutate(params = factor(params, levels = c("AR-1", "Independent", "Exchangeable")))

# Figure S5
p.sgmEpsl <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.sgmEpsl,
                                                    IJ = summ.pts.sgmEpsl), 
                                    x.label  = rep("Residual error variance", 2), 
                                    label.x  = c(F, T),
                                    x.breaks = 1:8,
                                    y.label  = "                 Total number of measurements across trials",
                                    log.y    = c(F, F),
                                    facet.x  = T,
                                    legend.position = "right",
                                    nrow.panel      = 2)

# Figure S6
p.rho <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.rho,
                                                IJ = summ.pts.rho), 
                                x.label  = rep("Correlation", 2), 
                                label.x  = c(F, T),
                                x.breaks = seq(0.1, 0.9, 0.1),
                                y.label  = "                 Total number of measurements across trials",
                                log.y    = c(F, F),
                                facet.x  = T,
                                legend.position = "right",
                                nrow.panel      = 2)

# Correlation structure
# Figure S4
legend.title <- "Correlation Structure"

p.corrStr.ppobs <- ggplot(summ.ppobs.corrStr %>%
                            filter(KL <= 75), aes(x = KL, y = mean.IJKL, group = params)) + 
  xlab("Number of measurements per participant") +
  geom_point(aes(color = params, shape = params), size = 2) +
  geom_line(aes(color = params, linetype = params)) +
  geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = params), alpha = 0.3) +
  facet_grid(. ~ slp) +
  # scale_y_continuous(limits = y.limits,
  #                    breaks = y.breaks) +
  scale_linetype_manual(name = legend.title,
                        values = c("dashed", "dotdash", "dotted", "longdash", "twodash")) +
  scale_color_manual(name = legend.title,
                     values = cbPalette[c(7, 5, 3, 4, 1, 6, 2, 8)]) +
  scale_fill_manual(name = legend.title,
                    values = cbPalette[c(7, 5, 3, 4, 1, 6, 2, 8)]) +
  scale_shape_manual(name = legend.title,
                     values = 1:8) +
  theme_bw() +
  theme(legend.position  = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.title.y     = element_blank()) 

p.corrStr.pts <- ggplot(summ.pts.corrStr %>%
                          filter(max.IJKL <= 1000), aes(x = IJ, y = mean.IJKL, group = params)) + 
  xlab("Number of participants") +
  geom_point(aes(color = params, shape = params), size = 2) +
  geom_line(aes(color = params, linetype = params)) +
  geom_ribbon(aes(ymin = min.IJKL, ymax = max.IJKL, fill = params), alpha = 0.3) +
  facet_grid(. ~ slp) +
  # scale_y_continuous(limits = y.limits,
  #                    breaks = y.breaks) +
  scale_linetype_manual(name = legend.title,
                        values = c("dashed", "dotdash", "dotted", "longdash", "twodash")) +
  scale_color_manual(name = legend.title,
                     values = cbPalette[c(7, 5, 3, 4, 6, 1, 2, 8)]) +
  scale_fill_manual(name = legend.title,
                    values = cbPalette[c(7, 5, 3, 4, 6, 1, 2, 8)]) +
  scale_shape_manual(name = legend.title,
                     values = 1:8) +
  theme_bw() +
  theme(legend.position  = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA),
        axis.title.y     = element_blank()) 

legend.corrStr  <- get_legend(p.corrStr.ppobs)
p.corrStr.ppobs <- p.corrStr.ppobs + theme(legend.position="none")
p.corrStr.pts   <- p.corrStr.pts + theme(legend.position="none")

p.corrStr <- grid.arrange(p.corrStr.ppobs, p.corrStr.pts, legend.corrStr, nrow = 3,
                          left = "                 Total number of measurements across trials",
                          heights = c(5, 5, 1))
