########################
# Assignment 1 -  Team 13
########################

########################
# Task 1
########################

# Set working directory
setwd("C:/Users/felix/Documents/felix/_ku-leuven/_kurse/1_multivariate-statistics/assignment-1")

# Clear workspace
rm(list = ls())

# Load libraries
library(paran)

# Load data
load("agdata.Rdata")

# Reverse items for consistent orientation
agdata$hmsfmlsh_r <- 6 - agdata$hmsfmlsh
agdata$rlgatnd_r  <- 8 - agdata$rlgatnd
agdata$pray_r     <- 8 - agdata$pray

########################
# Task 1 a
########################

# Select the 16 items/observed variables into one object
items1 <- agdata[, c(
  "ppltrst", "pplfair", "pplhlp",
  "freehms", "hmsfmlsh_r", "hmsacld",
  "imbgeco", "imueclt", "imwbcnt",
  "rlgdgr", "rlgatnd_r", "pray_r",
  "eqwrkbg", "eqpolbg", "eqmgmbg", "eqpaybg"
)]

# Check number of missing values for each item
colSums(is.na(items1))

# Standardize variables
zitems1 <- scale(items1, center = TRUE, scale = TRUE)

# Conduct PCA on standardized variables
pca <- prcomp(zitems1)

# Compute eigenvalues
# PCA produces 16 eigenvalues here: one for each principal component. 
# Each eigenvalue indicates how much of the total variance a component explains.
round(pca$sdev^2, 3)

# Kaiser's criterion: retain components with eigenvalues greater than 1
which(pca$sdev^2 > 1)

# Number of components retained according to Kaiser's criterion
sum(pca$sdev^2 > 1)
# Result: Retain 3 components

# Scree plot of eigenvalues: retain components before elbow
screeplot(pca, type = "lines", main = "Scree Plot")
abline(h = 1, col = "red", lty = 2)
# Result: Retain probably 2 components

# Proportion of explained variance of each component
round(pca$sdev^2 / 16, 3)

# Cumulative proportion of explained variance by components
round(cumsum(pca$sdev^2 / 16), 3)

# Procedure of Horn
# Helps to decide the number of components to be retained
set.seed(258)
paran(zitems1, iterations = 5000,
      graph = TRUE, cfa = FALSE, centile = 0)

########################
# Task 1 b
########################

# Compute component loadings 
# They are correlations between variables and components because variables were standardized.
component_loadings <- pca$rotation %*% diag(pca$sdev)

# Display component loadings and Variance Accounted For (VAF)
# For the two components retained based on procedure of Horn.
round(cbind(PC1 = component_loadings[, 1],
            PC2 = component_loadings[, 2],
            VAF = apply(component_loadings[, 1:2]^2, 1, sum)), 2)

# Compute unstandardized component scores for the 30 countries
component_scores <- predict(pca)

# Standardize component scores: Divide each component score by its standard deviation.
# >>> The standardized components then have mean 0 and standard deviation 1.
standardized_scores <- component_scores %*% diag(1 / pca$sdev)

# Display standardized scores on the first two components
# >> Gives us coordinates of all 30 countries in following biplot.
round(standardized_scores[, 1:2], 2)

# Biplot of the first two principal components
# Countries are represented by standardized component scores.
# Item arrows represent component loadings, i.e. correlations with PC1 and PC2.
# pc.biplot = TRUE uses standardized component scores for countries and component loadings for items!
par(cex = 1)
biplot(pca, pc.biplot = TRUE, cex = 0.8, xlab = "PC1", ylab = "PC2")














########################
# Task 2
########################

# Set working directory
setwd("C:/Users/felix/Documents/felix/_ku-leuven/_kurse/1_multivariate-statistics/assignment-1")

# Clear workspace
rm(list = ls())

# Load libraries
library(GPArotation)
library(lavaan)
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

# Show pattern loadings, total and common variance, factor correlations, model fit
efa

# Display all structure loadings, including small values (cutoff = 0), rounded (digits = 3)
print(efa$Structure, cutoff = 0, digits = 3)

# Compute residual correlations
resid <- cormat - (efa$loadings %*% efa$Phi %*% t(efa$loadings) + diag(efa$uniquenesses))

# Display matrix of residual correlations
round(resid, 3)

# Compute proportion of non-redundant residual correlations with absolute value greater 0.05
n <- sum(ifelse(abs(resid) > 0.05, 1, 0)) / 2
print(n / (16 * 15 / 2))

########################
# Task 2 b
########################

# Center items without standardizing their variances
citems <- scale(items, center = TRUE, scale = FALSE)

# Compute covariance matrix
covmat <- cov(citems)

# Specify CFA model with 5 correlated latent factors for the 5 concepts.
# Each item loads only on its intended factor. 
# The first loading of each factor is fixed to 1 for model identification. 
cfa5 <- '
  # Define measurement model
  social_trust        =~ 1*ppltrst + pplfair + pplhlp
  negative_gays       =~ 1*freehms + hmsfmlsh_r + hmsacld
  positive_immigrants =~ 1*imbgeco + imueclt + imwbcnt
  religiousness       =~ 1*rlgdgr + rlgatnd_r + pray_r
  gender_equality     =~ 1*eqwrkbg + eqpolbg + eqmgmbg + eqpaybg

  # Factor variances are freely estimated
  social_trust        ~~ social_trust
  negative_gays       ~~ negative_gays
  positive_immigrants ~~ positive_immigrants
  religiousness       ~~ religiousness
  gender_equality     ~~ gender_equality

  # Factor covariances
  social_trust        ~~ negative_gays + positive_immigrants + religiousness + gender_equality
  negative_gays       ~~ positive_immigrants + religiousness + gender_equality
  positive_immigrants ~~ religiousness + gender_equality
  religiousness       ~~ gender_equality
'

# Fit CFA model on  covariance matrix
fitcfa5 <- cfa(cfa5, sample.cov = covmat, sample.nobs = 2518)

# Summary of results
summary(fitcfa5, fit.measures = TRUE)

# Fit measures
fitmeasures(fitcfa5, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Show standardized loadings and factor correlations between the 5 factors
standSol <- standardizedSolution(fitcfa5)
standSol


# Convergent validity
# Extract standardized factor loadings
loadings <- subset(standSol, op == "=~")
loadings[, c("lhs", "rhs", "est.std", "pvalue")]

# Average variance extracted for each factor
# Assessing convergent validity:
# First of all: Standardized loadings should be positive and significant.
# AVE should exceed 0.50.
ave <- tapply(loadings$est.std^2, loadings$lhs, mean)
round(ave, 3)



# Discriminant validity
# Extract standardized correlations between different factors
factors <- unique(loadings$lhs)

# This line of the code to create "correlations" was generated using GenAI. 
correlations <- subset(standSol, op == "~~" &
                         lhs %in% factors &
                         rhs %in% factors &
                         lhs != rhs)

correlations[, c("lhs", "rhs", "est.std")]

# Maximum shared variance for each factor
msv <- c(
  social_trust = max(correlations$est.std[
    correlations$lhs == "social_trust" |
      correlations$rhs == "social_trust"]^2),
  
  negative_gays = max(correlations$est.std[
    correlations$lhs == "negative_gays" |
      correlations$rhs == "negative_gays"]^2),
  
  positive_immigrants = max(correlations$est.std[
    correlations$lhs == "positive_immigrants" |
      correlations$rhs == "positive_immigrants"]^2),
  
  religiousness = max(correlations$est.std[
    correlations$lhs == "religiousness" |
      correlations$rhs == "religiousness"]^2),
  
  gender_equality = max(correlations$est.std[
    correlations$lhs == "gender_equality" |
      correlations$rhs == "gender_equality"]^2)
)

# Compare AVE with MSV
report <- cbind(ave[factors], msv[factors])
colnames(report) <- c("AVE", "MSV")
rownames(report) <- factors

round(report, 3)

# Assess discriminant validity
report[, "AVE"] > report[, "MSV"]


########################
# Composite reliability
########################

# Define function to calculate composite reliability from standardized factor loadings
compositereliability <- function(x) {
  A <- (sum(x))^2
  B <- sum(1 - x^2)
  return(A / (A + B))
}

# Calculate composite reliability for each factor
reliab <- c(
  social_trust = compositereliability(
    loadings$est.std[loadings$lhs == "social_trust"]),
  negative_gays = compositereliability(
    loadings$est.std[loadings$lhs == "negative_gays"]),
  positive_immigrants = compositereliability(
    loadings$est.std[loadings$lhs == "positive_immigrants"]),
  religiousness = compositereliability(
    loadings$est.std[loadings$lhs == "religiousness"]),
  gender_equality = compositereliability(
    loadings$est.std[loadings$lhs == "gender_equality"])
)
round(reliab, 3)


########################
# Summary table
########################

summary_table <- data.frame(
  Factor = factors,
  AVE = ave[factors],
  MSV = msv[factors],
  CompositeReliability = reliab[factors]
)

print(summary_table, digits = 3, row.names = FALSE)








