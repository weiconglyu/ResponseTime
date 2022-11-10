# Code generating Tables 1 to 4
library(tidyverse)
library(lme4)
library(ltm)
library(optimx)
control.lmer <- lmerControl('optimx', optCtrl = list(method = 'nlminb'))

inv.logit <- function(x) {
    log(x / (1 - x))
}

load('input.RData')

eta <- t(apply(data$y, 1, function(y) {
    proportions(table(factor(y, 1:data$K)))
}))
rs <- sapply(1:data$N, function(n) {
    eta[n, as.integer(data$y[n, ])]
})

y <- t(apply(data$y, 1, function(y) {
    y[data$s == -1] <- 1 + data$K - y[data$s == -1]
    y
}))
prob <- do.call(c, lapply(1:data$D, function(d) {
    yy <- y[, data$d == d]
    fit <- grm(yy, control = list(iter.qN = 1000))
    fitted(fit, yy, 'conditional-probabilities')
}))
content <- sapply(1:data$N, function(n) {
    sapply(1:data$M, function(m) {
        prob[[m]][n, y[n, m]]
    })
})

dat <- data.frame(t = unlist(data$t), y = unlist(data$y),
                  content = as.numeric(t(content)), rs = as.numeric(t(rs))) %>%
    mutate(id = (row_number() - 1) %% data$N + 1,
           item = (row_number() - 1) %/% data$N + 1,
           z = factor(y),
           inv.logit.content = inv.logit(content),
           inv.logit.rs = inv.logit(rs))

# full model
summary(lmer(log(t) ~ inv.logit.content + inv.logit.rs + z +
                      (1 + inv.logit.content + inv.logit.rs | id) + (1 | item),
             dat, control = control.lmer))

# content only
summary(lmer(log(t) ~ inv.logit.content + z +
                      (1 + inv.logit.content | id) + (1 | item),
             dat, control = control.lmer))

# response style only
summary(lmer(log(t) ~ inv.logit.rs + z +
                      (1 + inv.logit.rs | id) + (1 | item),
             dat, control = control.lmer))
