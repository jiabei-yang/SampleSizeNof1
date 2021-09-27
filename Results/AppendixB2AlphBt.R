library(dplyr)
library(ggplot2)
library(gridExtra)

source("../Functions/gen.R")

# ppobs
load("../Data/main/mainPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.005.02 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a0025_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.0025.02 <- prep_data_fig(data.fig,
                                    x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a006_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.006.02 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a007_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.007.02 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a008_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.008.02 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a009_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.009.02 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a01_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.01.02 <- prep_data_fig(data.fig,
                                  x.axis = "KL")

load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a005_b005_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.005.005 <- prep_data_fig(data.fig,
                                    x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a005_b01_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.005.01 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a005_b0125_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.005.0125 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a005_b015_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.005.015 <- prep_data_fig(data.fig,
                                    x.axis = "KL")
load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a005_b0175_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.005.0175 <- prep_data_fig(data.fig,
                                     x.axis = "KL")
# load("../Data/AppendixB2AlphBt/AlphBtPpobs_Pair_Ar1_a005_b025_e4_r04_m4_t1_c1_d1.RData")
# summ.ppobs.005.025 <- prep_data_fig(data.fig,
#                                     x.axis = "KL")

summ.ppobs.alpha <- rbind(summ.ppobs.005.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.05),
                          summ.ppobs.0025.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.025),
                          summ.ppobs.006.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.06),
                          summ.ppobs.007.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.07),
                          summ.ppobs.008.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.08),
                          summ.ppobs.009.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.09),
                          summ.ppobs.01.02 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.1))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.ppobs.alpha %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.alpha <- summ.ppobs.alpha %>% 
  filter(KL %in% c(4, 6, 8, 10, 12))
summ.ppobs.alpha <- summ.ppobs.alpha %>% 
  mutate(KL = factor(KL, levels = c(4, 6, 8, 10, 12)))

summ.ppobs.beta <- rbind(summ.ppobs.005.02 %>% 
                           filter(intcpt == "Fixed Intercepts") %>%
                           mutate(params = 0.2),
                         summ.ppobs.005.005 %>% 
                           filter(intcpt == "Fixed Intercepts") %>%
                           mutate(params = 0.05),
                         summ.ppobs.005.01 %>% 
                           filter(intcpt == "Fixed Intercepts") %>%
                           mutate(params = 0.1),
                         summ.ppobs.005.0125 %>% 
                           filter(intcpt == "Fixed Intercepts") %>%
                           mutate(params = 0.125),
                         summ.ppobs.005.015 %>% 
                           filter(intcpt == "Fixed Intercepts") %>%
                           mutate(params = 0.15),
                         summ.ppobs.005.0175 %>% 
                           filter(intcpt == "Fixed Intercepts") %>%
                           mutate(params = 0.175))

# tmp <- rbind(summ.ppobs.005.02 %>% 
#                filter(intcpt == "Random Intercepts") %>%
#                mutate(params = 0.2),
#              summ.ppobs.005.005 %>% 
#                filter(intcpt == "Random Intercepts") %>%
#                mutate(params = 0.05),
#              summ.ppobs.005.01 %>% 
#                filter(intcpt == "Random Intercepts") %>%
#                mutate(params = 0.1),
#              summ.ppobs.005.0125 %>% 
#                filter(intcpt == "Random Intercepts") %>%
#                mutate(params = 0.125),
#              summ.ppobs.005.015 %>% 
#                filter(intcpt == "Random Intercepts") %>%
#                mutate(params = 0.15),
#              summ.ppobs.005.0175 %>% 
#                filter(intcpt == "Random Intercepts") %>%
#                mutate(params = 0.175))
# 
# plot(summ.ppobs.beta$mean.IJKL, type = "l")
# lines(tmp$mean.IJKL, col = "red")

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.ppobs.beta %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.beta <- summ.ppobs.beta %>% 
  filter(KL %in% c(4, 6, 8, 10, 12))
summ.ppobs.beta <- summ.ppobs.beta %>% 
  mutate(KL = factor(KL, levels = c(4, 6, 8, 10, 12)))

# pts
load("../Data/main/mainPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.005.02 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a0025_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.0025.02 <- prep_data_fig(data.fig,
                                  x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a006_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.006.02 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a007_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.007.02 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a008_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.008.02 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a009_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.009.02 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a01_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.01.02 <- prep_data_fig(data.fig,
                                x.axis = "IJ")

load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a005_b005_e4_r04_m4_t1_c1_d1.RData")
summ.pts.005.005 <- prep_data_fig(data.fig,
                                  x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a005_b01_e4_r04_m4_t1_c1_d1.RData")
summ.pts.005.01 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a005_b0125_e4_r04_m4_t1_c1_d1.RData")
summ.pts.005.0125 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a005_b015_e4_r04_m4_t1_c1_d1.RData")
summ.pts.005.015 <- prep_data_fig(data.fig,
                                  x.axis = "IJ")
load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a005_b0175_e4_r04_m4_t1_c1_d1.RData")
summ.pts.005.0175 <- prep_data_fig(data.fig,
                                   x.axis = "IJ")
# load("../Data/AppendixB2AlphBt/AlphBtPts_Pair_Ar1_a005_b025_e4_r04_m4_t1_c1_d1.RData")
# summ.pts.005.025 <- prep_data_fig(data.fig,
#                                   x.axis = "IJ")

summ.pts.alpha <- rbind(summ.pts.005.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.05),
                        summ.pts.0025.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.025),
                        summ.pts.006.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.06),
                        summ.pts.007.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.07),
                        summ.pts.008.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.08),
                        summ.pts.009.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.09),
                        summ.pts.01.02 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.1))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.pts.alpha %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.alpha <- summ.pts.alpha %>% 
  filter(IJ %in% c(12, 14, 16, 20, 24))
summ.pts.alpha <- summ.pts.alpha %>% 
  mutate(IJ = factor(IJ, levels = c(12, 14, 16, 20, 24)))

summ.pts.beta <- rbind(summ.pts.005.02 %>% 
                         filter(intcpt == "Fixed Intercepts") %>%
                         mutate(params = 0.2),
                       summ.pts.005.005 %>% 
                         filter(intcpt == "Fixed Intercepts") %>%
                         mutate(params = 0.05),
                       summ.pts.005.01 %>% 
                         filter(intcpt == "Fixed Intercepts") %>%
                         mutate(params = 0.1),
                       summ.pts.005.0125 %>% 
                         filter(intcpt == "Fixed Intercepts") %>%
                         mutate(params = 0.125),
                       summ.pts.005.015 %>% 
                         filter(intcpt == "Fixed Intercepts") %>%
                         mutate(params = 0.15),
                       summ.pts.005.0175 %>% 
                         filter(intcpt == "Fixed Intercepts") %>%
                         mutate(params = 0.175))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.pts.beta %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.beta <- summ.pts.beta %>% 
  filter(IJ %in% c(12, 14, 20, 24, 28))
summ.pts.beta <- summ.pts.beta %>% 
  mutate(IJ = factor(IJ, levels = c(12, 14, 20, 24, 28)))


# Figure S2
p.alpha <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.alpha,
                                                  IJ = summ.pts.alpha), 
                                  x.label  = rep("Type I error", 2), 
                                  label.x  = c(F, T),
                                  x.breaks = c(0.025, seq(0.05, 0.1, 0.01)),
                                  y.label  = "                 Total number of measurements across trials",
                                  log.y    = c(F, T),
                                  facet.x  = T,
                                  legend.position = "right",
                                  nrow.panel      = 2)

# Figure S3
p.beta <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.beta,
                                                 # %>% filter(slp == "Common Slope"),
                                                  IJ = summ.pts.beta), 
                                  x.label  = rep("Type II error", 2), 
                                  label.x  = c(F, T),
                                  x.breaks = c(0.05, seq(0.1, 0.2, 0.025), 0.25),
                                  y.label  = "                 Total number of measurements across trials",
                                  log.y    = c(F, T),
                                  facet.x  = T,
                                  legend.position = "right",
                                  nrow.panel      = 2)

