########################
# Assignment 1 -  Team 13
########################

# Set working directory
setwd("C:/Users/felix/Documents/felix/_ku-leuven/_kurse/1_multivariate-statistics/assignment-1")


########################
# Task 1
########################



########################
# Task 2
########################

# Load libraries
library(psych)

# Load data
load("ndata.RData")

# Reverse items for consistent orientation
ndata$hmsfmlsh_r <- 6 - ndata$hmsfmlsh
ndata$rlgatnd_r  <- 8 - ndata$rlgatnd
ndata$pray_r     <- 8 - ndata$pray

########################
# Task 2 a
########################

# Select the 16 items/observed variables into one object
items <- ndata[, c(
  "ppltrst", "pplfair", "pplhlp",
  "freehms", "hmsfmlsh_r", "hmsacld",
  "imbgeco", "imueclt", "imwbcnt",
  "rlgdgr", "rlgatnd_r", "pray_r",
  "eqwrkbg", "eqpolbg", "eqmgmbg", "eqpaybg"
)]

# Check number of missing values for each item
colSums(is.na(items))

# Compute the observed correlation matrix
cormat <- cor(items)

# EFA on correlation matrix with 5 factors and oblique rotation
# So: 5 factors, maximum likelihood estimation, oblique rotation
efa <- fa(cormat, covar = FALSE, nfactors = 5,
          rotate = "oblimin", fm = "ml", n.obs = 2518)

# Show pattern loadings, communalities, common and explained variance, factor correlations, model fit
efa

# Display all structure loadings, including small values (cutoff = 0), rounded to 3 decimal places (digits = 3)
print(efa$Structure, cutoff = 0, digits = 3)

# Display matrix of residual correlations
round(efa$residual, 3)

# Compute proportion of non-redundant residual correlations with absolute values greater than 0.05
sum(abs(efa$residual[upper.tri(efa$residual)]) > 0.05) / (16 * 15 / 2)



















