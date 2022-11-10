# Code generating Tables A1 to A4
library(rstan)
options(mc.cores = 4)
rstan_options(auto_write = T)

load('input.RData')

fit <- stan('appendix.stan', data = data, iter = 10000,
            warmup = 5000, thin = 10, include = F,
            pars = c('z_m', 'z_n', 'z_theta', 'tau_m', 'tau_n',
                     'L_omega_m', 'L_omega_n', 'L_omega_theta'))
est <- summary(fit)$summary

library(tidyverse)
library(lme4)
library(optimx)
control <- lmerControl('optimx', optCtrl = list(method = 'nlminb'))

dat <- data.frame(t = as.vector(t(data$t)),
                  ans = as.vector(t(data$y)),
                  c = est[grep('^v_c\\[', rownames(est)), 1],
                  theta = est[grep('^v_theta\\[', rownames(est)), 1],
                  eta = est[grep('^v_eta\\[', rownames(est)), 1]) %>%
        mutate(n = row_number() - 1, id = n %/% data$M + 1,
               item = n %% data$M + 1, n = NULL) %>%
        mutate(log.t = log(t), thetac = theta + c, z = factor(ans))

# full model
summary(lmer(log.t ~ z + thetac + eta + (1 + thetac + eta | id) + (1 | item),
             dat, control = control))

# content only
summary(lmer(log.t ~ z + thetac + (1 + thetac | id) + (1 | item), dat,
        control = control))

# response style only
summary(lmer(log.t ~ z + eta + (1 + eta | id) + (1 | item), dat, control = control))
