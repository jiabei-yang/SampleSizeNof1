library(dplyr)
library(ggplot2)
library(gridExtra)

source("../Functions/gen.R")

# ppobs
load("../Data/main/mainPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.4.1.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m1_t1_c1_d1.RData")
summ.ppobs.1.1.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m2_t1_c1_d1.RData")
summ.ppobs.2.1.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m3_t1_c1_d1.RData")
summ.ppobs.3.1.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m5_t1_c1_d1.RData")
summ.ppobs.5.1.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m6_t1_c1_d1.RData")
summ.ppobs.6.1.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")

load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t2_c1_d1.RData")
summ.ppobs.4.2.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t3_c1_d1.RData")
summ.ppobs.4.3.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t4_c1_d1.RData")
summ.ppobs.4.4.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t5_c1_d1.RData")
summ.ppobs.4.5.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t6_c1_d1.RData")
summ.ppobs.4.6.1 <- prep_data_fig(data.fig,
                                  x.axis = "KL")

load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c0_d1.RData")
summ.ppobs.4.1.0 <- prep_data_fig(data.fig,
                                  x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c02_d1.RData")
summ.ppobs.4.1.02 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c06_d1.RData")
summ.ppobs.4.1.06 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c14_d1.RData")
summ.ppobs.4.1.14 <- prep_data_fig(data.fig,
                                   x.axis = "KL")
load("../Data/AppendixB4RandEff/RandEffPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c18_d1.RData")
summ.ppobs.4.1.18 <- prep_data_fig(data.fig,
                                   x.axis = "KL")

# random intercepts
# only random intercepts is possible
summ.ppobs.sgmMu <- rbind(summ.ppobs.4.1.1 %>% 
                            filter(intcpt == "Random Intercepts") %>%
                            mutate(params = 4),
                          summ.ppobs.1.1.1 %>% 
                            filter(intcpt == "Random Intercepts") %>%
                            mutate(params = 1),
                          summ.ppobs.2.1.1 %>% 
                            filter(intcpt == "Random Intercepts") %>%
                            mutate(params = 2),
                          summ.ppobs.3.1.1 %>% 
                            filter(intcpt == "Random Intercepts") %>%
                            mutate(params = 3),
                          summ.ppobs.5.1.1 %>% 
                            filter(intcpt == "Random Intercepts") %>%
                            mutate(params = 5),
                          summ.ppobs.6.1.1 %>% 
                            filter(intcpt == "Random Intercepts") %>%
                            mutate(params = 6))

# still summarize by KL if plot sgmMu on x-axis; 
tmp.1 <- summ.ppobs.sgmMu %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.sgmMu <- summ.ppobs.sgmMu %>% 
  filter(KL %in% c(4, 8, 16, 33, 60))
summ.ppobs.sgmMu <- summ.ppobs.sgmMu %>% 
  mutate(KL = factor(KL, levels = c(4, 8, 16, 33, 60)))
# 4, 8, 16, 33, 60

# random slopes
# only for random slopes; checked: the same for random intercepts and fixed intercepts
summ.ppobs.sgmTau <- rbind(summ.ppobs.4.1.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 1),
                           summ.ppobs.4.2.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 2),
                           summ.ppobs.4.3.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 3),
                           summ.ppobs.4.4.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 4),
                           summ.ppobs.4.5.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 5),
                           summ.ppobs.4.6.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 6))

# tmp <- rbind(summ.ppobs.4.1.1 %>% 
#                              filter((intcpt == "Fixed Intercepts") & (slp == "Random Slopes")) %>%
#                              mutate(params = 1),
#                            summ.ppobs.4.2.1 %>% 
#                              filter((intcpt == "Fixed Intercepts") & (slp == "Random Slopes")) %>%
#                              mutate(params = 2),
#                            summ.ppobs.4.4.1 %>% 
#                              filter((intcpt == "Fixed Intercepts") & (slp == "Random Slopes")) %>%
#                              mutate(params = 4),
#                            summ.ppobs.4.6.1 %>% 
#                              filter((intcpt == "Fixed Intercepts") & (slp == "Random Slopes")) %>%
#                              mutate(params = 6))

# still summarize by KL if plot tau on x-axis; 
tmp.1 <- summ.ppobs.sgmTau %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.sgmTau <- summ.ppobs.sgmTau %>% 
  filter(KL %in% c(4, 8, 12, 24))
summ.ppobs.sgmTau <- summ.ppobs.sgmTau %>% 
  mutate(KL = factor(KL, levels = c(4, 8, 12, 24)))

# only random intercepts and random slopes model have correlation between random effects
summ.ppobs.sgmMuTau <- rbind(summ.ppobs.4.1.1 %>% 
                               filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                               mutate(params = 1 / (2 * 1)),
                             summ.ppobs.4.1.0 %>% 
                               filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                               mutate(params = 0),
                             summ.ppobs.4.1.02 %>% 
                               filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                               mutate(params = 0.2 / (2 * 1)),
                             summ.ppobs.4.1.06 %>% 
                               filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                               mutate(params = 0.6 / (2 * 1)),
                             summ.ppobs.4.1.14 %>% 
                               filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                               mutate(params = 1.4 / (2 * 1)),
                             summ.ppobs.4.1.18 %>% 
                               filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                               mutate(params = 1.8 / (2 * 1)))

# still summarize by KL if plot sgmMu on x-axis; 
tmp.1 <- summ.ppobs.sgmMuTau %>%
  group_by(KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(KL) %>%
  summarise(n = n())

summ.ppobs.sgmMuTau <- summ.ppobs.sgmMuTau %>% 
  filter(KL %in% c(4, 9, 16, 32))
summ.ppobs.sgmMuTau <- summ.ppobs.sgmMuTau %>% 
  mutate(KL = factor(KL, levels = c(4, 9, 16, 32)))

# pts
load("../Data/main/mainPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.4.1.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m1_t1_c1_d1.RData")
summ.pts.1.1.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m2_t1_c1_d1.RData")
summ.pts.2.1.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m3_t1_c1_d1.RData")
summ.pts.3.1.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m5_t1_c1_d1.RData")
summ.pts.5.1.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m6_t1_c1_d1.RData")
summ.pts.6.1.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")

load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t2_c1_d1.RData")
summ.pts.4.2.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t3_c1_d1.RData")
summ.pts.4.3.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t4_c1_d1.RData")
summ.pts.4.4.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t5_c1_d1.RData")
summ.pts.4.5.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t6_c1_d1.RData")
summ.pts.4.6.1 <- prep_data_fig(data.fig,
                                x.axis = "IJ")

load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c0_d1.RData")
summ.pts.4.1.0 <- prep_data_fig(data.fig,
                                x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c02_d1.RData")
summ.pts.4.1.02 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c06_d1.RData")
summ.pts.4.1.06 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c14_d1.RData")
summ.pts.4.1.14 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")
load("../Data/AppendixB4RandEff/RandEffPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c18_d1.RData")
summ.pts.4.1.18 <- prep_data_fig(data.fig,
                                 x.axis = "IJ")

# random intercepts
summ.pts.sgmMu <- rbind(summ.pts.4.1.1 %>% 
                          filter(intcpt == "Random Intercepts") %>%
                          mutate(params = 4),
                        summ.pts.1.1.1 %>% 
                          filter(intcpt == "Random Intercepts") %>%
                          mutate(params = 1),
                        summ.pts.2.1.1 %>% 
                          filter(intcpt == "Random Intercepts") %>%
                          mutate(params = 2),
                        summ.pts.3.1.1 %>% 
                          filter(intcpt == "Random Intercepts") %>%
                          mutate(params = 3),
                        summ.pts.5.1.1 %>% 
                          filter(intcpt == "Random Intercepts") %>%
                          mutate(params = 5),
                        summ.pts.6.1.1 %>% 
                          filter(intcpt == "Random Intercepts") %>%
                          mutate(params = 6))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.pts.sgmMu %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.sgmMu <- summ.pts.sgmMu %>% 
  filter(IJ %in% c(12, 16, 20, 24, 30))
summ.pts.sgmMu <- summ.pts.sgmMu %>% 
  mutate(IJ = factor(IJ, levels = c(12, 16, 20, 24, 30)))

# random slopes
summ.pts.sgmTau <- rbind(summ.pts.4.1.1 %>% 
                           filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                           mutate(params = 1),
                         summ.pts.4.2.1 %>% 
                           filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                           mutate(params = 2),
                         summ.pts.4.3.1 %>% 
                           filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                           mutate(params = 3),
                         summ.pts.4.4.1 %>% 
                           filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                           mutate(params = 4),
                         summ.pts.4.5.1 %>% 
                           filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                           mutate(params = 5),
                         summ.pts.4.6.1 %>% 
                           filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                           mutate(params = 6))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.pts.sgmTau %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.sgmTau <- summ.pts.sgmTau %>% 
  filter(IJ %in% c(42, 46, 50, 54))
summ.pts.sgmTau <- summ.pts.sgmTau %>% 
  mutate(IJ = factor(IJ, levels = c(42, 46, 50, 54)))

# only random intercepts and random slopes model have correlation between random effects
summ.pts.sgmMuTau <- rbind(summ.pts.4.1.1 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 1 / (2 * 1)),
                           summ.pts.4.1.0 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 0),
                           summ.pts.4.1.02 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 0.2 / (2 * 1)),
                           summ.pts.4.1.06 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 0.6 / (2 * 1)),
                           summ.pts.4.1.14 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 1.4 / (2 * 1)),
                           summ.pts.4.1.18 %>% 
                             filter((intcpt == "Random Intercepts") & (slp == "Random Slopes")) %>%
                             mutate(params = 1.8 / (2 * 1)))

# still summarize by KL if plot sgmMu on x-axis; 
tmp.1 <- summ.pts.sgmMuTau %>%
  group_by(IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(IJ) %>%
  summarise(n = n())

summ.pts.sgmMuTau <- summ.pts.sgmMuTau %>% 
  filter(IJ %in% c(16, 22, 28, 36, 46))
summ.pts.sgmMuTau <- summ.pts.sgmMuTau %>% 
  mutate(IJ = factor(IJ, levels = c(16, 22, 28, 36, 46)))

# Figure S7
p.sgmMu <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.sgmMu,
                                                  IJ = summ.pts.sgmMu), 
                                  x.label  = rep("Variance of random intercepts", 2), 
                                  label.x  = c(F, T),
                                  x.breaks = 1:6,
                                  y.label  = "                 Total number of measurements across trials",
                                  log.y    = c(F, F),
                                  facet.x  = T,
                                  legend.position = "right",
                                  nrow.panel      = 2)

# Figure S9
p.sgmTau <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.sgmTau,
                                                   IJ = summ.pts.sgmTau), 
                                   x.label  = rep("Variance of random slopes", 2), 
                                   label.x  = c(T, T),
                                   x.breaks = 1:6,
                                   y.label  = "Total number of measurements across trials",
                                   log.y    = c(F, F),
                                   facet.x  = F,
                                   legend.position = "bottom",
                                   nrow.panel = 1)

# Figure S8
p.sgmMuTau <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.sgmMuTau,
                                                     IJ = summ.pts.sgmMuTau), 
                                     x.label  = rep("Correlation of random intercepts and slopes", 2), 
                                     label.x  = c(T, T),
                                     x.breaks = c(0, seq(0.1, 0.9, 0.2)),
                                     y.label  = "Total number of measurements across trials",
                                     log.y    = c(F, F),
                                     facet.x  = F,
                                     legend.position = "bottom",
                                     nrow.panel = 1)

