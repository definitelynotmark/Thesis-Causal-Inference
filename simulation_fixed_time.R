## @knitr simulation_fixed_time_setup
library(tidyverse)
library(magrittr)
theme_set(theme_bw())
options(ggplot2.discrete.colour = function(...) scale_colour_brewer(..., palette = "Set2"))
options(ggplot2.discrete.fill = function(...) scale_fill_brewer(..., palette = "Set2"))

## @knitr simulation_fixed_time_body
# model parameters
n <- 1000
p_L <- 0.3 # prob of L = 1
p_A <- 0.5 # prob of A = 1 in RCT
p_A0 <- 1 / 4 # prob of A = 1 if L = 0
p_A1 <- 3 / 4 # prob of A = 1 if L = 1

shape <- 2 # Weibull shape
rate_L <- 3 # rate mult factor if L = 1
rate_A <- 0.5 # rate mult factor if A = 1

# simulation of data
set.seed(1)
L <- rbinom(n, 1, p_L)
A_RCT <- rbinom(n, 1, p_A)
A_conf <- rbinom(n, 1, ifelse(L, p_A1, p_A0))
T0 <- rweibull(n, shape, 1 / ifelse(L, rate_L, 1))
T1 <- rweibull(n, shape, 1 / ifelse(L, rate_L, 1) / rate_A)

T_RCT <- ifelse(A_RCT, T1, T0)
T_conf <- ifelse(A_conf, T1, T0)

# estimators of ATE
# naive difference between treatment groups
S_assoc <- function(A, T){
  T0 <- T[A == 0]
  T1 <- T[A == 1]
  
  S0 <- Vectorize(function(s) mean(T0 > s))
  S1 <- Vectorize(function(s) mean(T1 > s))
  
  data.frame(t = T, S0 = S0(T), S1 = S1(T)) %>%
    mutate(ATE = S1 - S0) %>%
    arrange(t)
}

# g-formula estimator
S_g <- function(L, A, T){
  T0 <- subset(T, A == 0)
  T1 <- subset(T, A == 1)
  T_L0A0 <- subset(T, A == 0 & L == 0)
  T_L0A1 <- subset(T, A == 1 & L == 0)
  T_L1A0 <- subset(T, A == 0 & L == 1)
  T_L1A1 <- subset(T, A == 1 & L == 1)
  
  p_l <- mean(L)
  
  S0 <- Vectorize(function(t) (1 - p_l) * mean(T_L0A0 > t) + 
                    p_l * mean(T_L1A0 > t))
  S1 <- Vectorize(function(t) (1 - p_l) * mean(T_L0A1 > t) + 
                    p_l * mean(T_L1A1 > t))
  
  data.frame(t = T, S0 = S0(T), S1 = S1(T)) %>%
    mutate(ATE = S1 - S0) %>%
    arrange(t)
}

# IPW estimator
S_IPW <- function(L, A, T){
  T0 <- subset(T, A == 0)
  T1 <- subset(T, A == 1)
  p_0 <- mean(subset(A, L == 0))
  p_1 <- mean(subset(A, L == 1))
  
  p <- function(a, l){
    ifelse(l, ifelse(a, p_1, 1 - p_1), ifelse(a, p_0, 1 - p_0))
  }
  
  S0 <- Vectorize(function(t) mean((A == 0) * (T > t) / p(0, L)))
  S1 <- Vectorize(function(t) mean((A == 1) * (T > t) / p(1, L)))
  
  data.frame(t = T, S0 = S0(T), S1 = S1(T)) %>%
    mutate(ATE = S1 - S0) %>%
    arrange(t)
}

## @knitr simulation_fixed_time_plot
# plot
tstart <- 0
tstop <- 4
trueATE <- data.frame(t = seq(tstart, tstop, length.out = 1000)) %>%
  mutate(ATE = (1 - p_L) * pweibull(t, shape, 1 / rate_A,  lower.tail = FALSE) + 
           p_L * pweibull(t, shape, 1 / rate_L / rate_A,  lower.tail = FALSE) - 
           ((1 - p_L) * pweibull(t, shape,  lower.tail = FALSE) + 
              p_L * pweibull(t, shape, 1 / rate_L,  lower.tail = FALSE) ))

ggplot(mapping = aes(x = t, y = ATE)) +
  scale_x_continuous(limits = c(tstart, tstop)) +
  scale_y_continuous(name = "Average Treatment Effect", labels = scales::percent) +
  geom_line(aes(color = "True"), data = trueATE) +
  geom_step(aes(color = "RCT"), data = S_assoc(A_RCT, T_RCT)) +
  geom_step(aes(color = "Confounded"), data = S_assoc(A_conf, T_conf)) +
  geom_step(aes(color = "g/IPW"), data = S_IPW(L, A_conf, T_conf)) +
  guides(color = guide_legend(NULL), fill = guide_legend(NULL))
