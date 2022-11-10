# Data cleaning
# Raw data can be downloaded from https://openpsychometrics.org/_rawdata/
library(tidyverse)

# MIES
set.seed(428571)

all <- na.omit(read.csv('data.csv', sep = '\t'))
y <- all[, grep('Q[0-9]*A', names(all))]
all <- all[apply(y, 1, sd) > 0 & rowSums(y <= 0) == 0, ]
t <- all[, grep('Q[0-9]*E', names(all))]
all <- all[rowSums(t <= 0 | t >= 60000) == 0, ]
sample <- all[sample.int(nrow(all), 500), ]

y <- sample[, grep('Q[0-9]*A', names(sample))]
z <- t(apply(as.matrix(sample[, grep('Q[0-9]*I', names(sample))]), 1, order))
t <- sample[, grep('Q[0-9]*E', names(sample))]
s <- as.vector(sign(cor(y[, 1], y)))

data <- list(N = nrow(y), M = ncol(y), K = max(y), D = 1, y = y, z = z, d = rep(1, length(y)), t = t, s = s)
save(data, sample, file = 'input_mies.RData')


# IPIP
set.seed(857142)

all <- na.omit(read.csv('data.csv', sep = '\t', na.strings = 'NULL'))
y <- all[, grep('^...[0-9]+$', names(all))]
all <- all[apply(y, 1, sd) > 0 & rowSums(y <= 0) == 0, ]
t <- all[, grep('^...[0-9]+_E$', names(all))]
all <- all[rowSums(t <= 0 | t >= 60000) == 0, ]
sample <- all[sample.int(nrow(all), 500), ]

y <- sample[, grep('^...[0-9]+$', names(sample))]
t <- sample[, grep('^...[0-9]+_E$', names(sample))]
d <- rep(1:5, each = 10)
D <- 5
s <- as.vector(sapply(1:D, function(k) {
    sign(cor(y[, d == k])[1, ])
}))

data <- list(N = nrow(y), M = ncol(y), K = max(y), D = D, y = y, d = d, t = t, s = s)
save(data, sample, file = 'input_ipip.RData')


# DASS
set.seed(142857)

all <- na.omit(read.csv('data.csv', sep = '\t'))
y <- all[, grep('Q[0-9]*A', names(all))]
all <- all[apply(y, 1, sd) > 0 & rowSums(y <= 0) == 0, ]
t <- all[, grep('Q[0-9]*E', names(all))]
all <- all[rowSums(t <= 0 | t >= 60000) == 0, ]
sample <- all[sample.int(nrow(all), 500), ]

y <- sample[, grep('Q[0-9]*A', names(sample))]
z <- t(apply(as.matrix(sample[, grep('Q[0-9]*I', names(sample))]), 1, order))
t <- sample[, grep('Q[0-9]*E', names(sample))]
s <- as.vector(sign(cor(y[, 1], y)))

data <- list(N = nrow(y), M = ncol(y), K = max(y), D = 1, y = y, z = z, d = rep(1, length(y)), t = t, s = s)
save(data, sample, file = 'input_dass.RData')


# FTI
set.seed(285714)

all <- na.omit(read.csv('data.csv', sep = '\t'))
y <- all[, grep('Q[0-9]*A', names(all))]
all <- all[apply(y, 1, sd) > 0 & rowSums(y <= 0) == 0, ]
t <- all[, grep('Q[0-9]*E', names(all))]
all <- all[rowSums(t <= 0 | t >= 60000) == 0, ]
sample <- all[sample.int(nrow(all), 500), ]

y <- sample[, grep('Q[0-9]*A', names(sample))]
z <- t(apply(as.matrix(sample[, grep('Q[0-9]*I', names(sample))]), 1, order))
t <- sample[, grep('Q[0-9]*E', names(sample))]
d <- rep(1:4, each = 14)
D <- 4
s <- as.vector(sapply(1:D, function(k) {
    sign(cor(y[, d == k])[1, ])
}))

data <- list(N = nrow(y), M = ncol(y), K = max(y), D = D, y = y, z = z, d = d, t = t, s = s)
save(data, sample, file = 'input_fti.RData')
