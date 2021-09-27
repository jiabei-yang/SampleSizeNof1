library(dplyr)
library(ggplot2)
library(gridExtra)

source("../Functions/gen.R")

# ppobs
load("../Data/main/mainPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.ppobs.1 <- prep_data_fig(data.fig,
                              x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d05.RData")
summ.ppobs.05 <- prep_data_fig(data.fig,
                               x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d15.RData")
summ.ppobs.15 <- prep_data_fig(data.fig,
                               x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d2.RData")
summ.ppobs.2 <- prep_data_fig(data.fig,
                              x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d25.RData")
summ.ppobs.25 <- prep_data_fig(data.fig,
                               x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d3.RData")
summ.ppobs.3 <- prep_data_fig(data.fig,
                              x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d35.RData")
summ.ppobs.35 <- prep_data_fig(data.fig,
                               x.axis = "KL")
load("../Data/AppendixB5Delta/DeltaPpobs_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d4.RData")
summ.ppobs.4 <- prep_data_fig(data.fig,
                              x.axis = "KL")

summ.ppobs.delta <- rbind(summ.ppobs.1 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 1),
                          summ.ppobs.05 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 0.5),
                          summ.ppobs.15 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 1.5),
                          summ.ppobs.2 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 2),
                          summ.ppobs.25 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 2.5),
                          summ.ppobs.3 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 3),
                          summ.ppobs.35 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 3.5),
                          summ.ppobs.4 %>% 
                            filter(intcpt == "Fixed Intercepts") %>%
                            mutate(params = 4))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.ppobs.delta %>%
  group_by(slp, KL, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, KL) %>%
  summarise(n = n())

summ.ppobs.delta <- summ.ppobs.delta %>% 
  filter(KL %in% c(4, 6, 8, 12))
summ.ppobs.delta <- summ.ppobs.delta %>% 
  mutate(KL = factor(KL, levels = c(4, 6, 8, 12)))

# pts
load("../Data/main/mainPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d1.RData")
summ.pts.1 <- prep_data_fig(data.fig,
                            x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d05.RData")
summ.pts.05 <- prep_data_fig(data.fig,
                             x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d15.RData")
summ.pts.15 <- prep_data_fig(data.fig,
                             x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d2.RData")
summ.pts.2 <- prep_data_fig(data.fig,
                            x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d25.RData")
summ.pts.25 <- prep_data_fig(data.fig,
                             x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d3.RData")
summ.pts.3 <- prep_data_fig(data.fig,
                            x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d35.RData")
summ.pts.35 <- prep_data_fig(data.fig,
                             x.axis = "IJ")
load("../Data/AppendixB5Delta/DeltaPts_Pair_Ar1_a005_b02_e4_r04_m4_t1_c1_d4.RData")
summ.pts.4 <- prep_data_fig(data.fig,
                            x.axis = "IJ")

summ.pts.delta <- rbind(summ.pts.1 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 1),
                        summ.pts.05 %>%  
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 0.5),
                        summ.pts.15 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 1.5),
                        summ.pts.2 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 2),
                        summ.pts.25 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 2.5),
                        summ.pts.3 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 3),
                        summ.pts.35 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 3.5),
                        summ.pts.4 %>% 
                          filter(intcpt == "Fixed Intercepts") %>%
                          mutate(params = 4))

# still summarize by KL if plot epsilon on x-axis; 
tmp.1 <- summ.pts.delta %>%
  group_by(slp, IJ, params) %>%
  summarise(n = n())
tmp.2 <- tmp.1 %>%
  group_by(slp, IJ) %>%
  summarise(n = n())

summ.pts.delta <- summ.pts.delta %>% 
  filter(IJ %in% c(4, 6, 8, 10, 12))
summ.pts.delta <- summ.pts.delta %>% 
  mutate(IJ = factor(IJ, levels = c(4, 6, 8, 10, 12)))

# Figure S10
p.delta <- plot.app.popavg.params(data.fig = list(KL = summ.ppobs.delta,
                                                  IJ = summ.pts.delta), 
                                  x.label  = rep("Minimal clinically important treatment effect", 2), 
                                  label.x  = c(F, T),
                                  x.breaks = seq(0.5, 4, 0.5),
                                  y.label  = "                 Total number of measurements across trials",
                                  log.y    = c(F, T),
                                  facet.x  = T,
                                  legend.position = "right",
                                  nrow.panel      = 2)
