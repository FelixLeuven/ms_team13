########################
# Assignment 1 -  Team 13
########################

########################
# Task 1
########################

load("agdata.Rdata")

#a 

#reversing hmsfmlsh, rlgatnd and pray such that, within each concept, a higher score indicates more of the concept it measures (negative attitude towards gays/lesbians, religiousness)
agdata$hmsfmlsh <- 6 - agdata$hmsfmlsh
agdata$rlgatnd <- 8 - agdata$rlgatnd
agdata$pray    <- 8 - agdata$pray

#select the 16 items, so that no other columns enter the analysis
items = agdata[, c("ppltrst", "pplfair", "pplhlp",
                   "freehms", "hmsfmlsh", "hmsacld",
                   "imbgeco", "imueclt", "imwbcnt",
                   "rlgdgr", "rlgatnd", "pray",
                   "eqwrkbg", "eqpolbg", "eqmgmbg", "eqpaybg")]

#check number of missing values per item 
colSums(is.na(items))

#standardise variables because of different scales of items 
zagdata = scale(items, center=TRUE, scale=TRUE)

#analysis of standardised agdata with PCA 
pcaagdata = prcomp(zagdata)

#eigenvalues
round(pcaagdata$sdev^2,2)

#proportion of explained variance per component
#total variance = 16 (because of standardised variables) 
round(pcaagdata$sdev^2/16,2) 

#cumulative proportion 
round(cumsum(pcaagdata$sdev^2)/16,2) #the first two components explain together 81% of the variance in the data 

#scree plot 
screeplot(pcaagdata, type="lines")
abline(h=1, col="red") 
sum(pcaagdata$sdev^2 > 1)
#elbow at component 2, and little crack at component 3 (scree plot: retain 1, or maybe even 2 components)
#Kaiser's rule: retain 3 components, kaiser's rule

#Horn's procedure
library(paran)

#mean of bootstrapped eigenvalues 
set.seed(589)
paran(zagdata,iterations=5000, graph=TRUE, cfa=FALSE, centile=0) #2 components retained 

#95th percentile (more conservative)
set.seed(897)
paran(zagdata, iterations=5000, graph=TRUE, cfa=FALSE, centile=95) #1 component retained (borderline: could include a second component -> see plot)

#compute component loadings
A = pcaagdata$rotation%*%diag(pcaagdata$sdev)
round(A[,1:2],2)
#PC1 is a general liberal/trusting/secular vs conservative/low-trust/religious dimension 
#the items on negative attitude towards gays/lesbians load in the same direction as the religious items (conservative side)
#PC2 is a gender-equality vs social trust dimension 

#proportion of variance accounted for (VAF) in each item by the first and first two principal components
round(cbind(VAF1 = diag(A[, 1] %*% t(A[, 1])), VAF2 = diag(A[, 1:2] %*% t(A[, 1:2]))), 2)
#VAF1: low VAF value for eqwrkbg (0.17) and the religious items: rlgdgr (0.44) and pray (0.52)
#VAF 2: not the case if you retain the first 2 components (this makes a case for retaining 2 principal components) 

#b

#biplots are not perfect representations because they only consider the first 2 principal components 
#this biplot, however, explains 81% of the variance, so it is a good approximation of the data
#pc.biplot=TRUE: countries are represented by standardised component scores, items by component loadings
biplot(pcaagdata, pc.biplot=TRUE, xlim=c(-3,3), ylim=c(-3,3))

#squared length of vectors = proportion of variance that is accounted for by the first two principal components 
round(apply(A[,1:2]^2,1,sum),2) 
#eqmgmbg has the longest vector -> the two components account for the largest share of eqmgmbg’s variance.
#rlgdgr has the shortest vector -> the two components account for the smallest share of rlgdgr's variance

#angle between 2 vectors provides info about the predicted correlation between the corresponding pair of variables 
#predicted correlation between variables 
round(A[,1:2]%*%t(A[,1:2]),2)
#eqpolbg and eqmgmbg have the highest correlation (0.95) -> countries that see gender equality in politics as positive also see it as positive in business management.
#eqwrkbg and pplhlp have the weakest correlation (0) -> a country’s view on women in paid work says nothing about how helpful its people are considered to be.
#all gender equality items are highly correlated with each other
#eqwrkbg is almost perpendicular to the trust and religious items (predicted correlation close to 0)

#standardised component scores 
Zs = zagdata%*%pcaagdata$rotation%*%diag(1/pcaagdata$sdev)
#print list of countries' scores on PC1 and PC2
round(Zs[,1:2],2)


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


########################
# Task 2 c
########################

# Preparation of centered items and country variable for multi-group SEM
semdata <- data.frame(country = ndata$cntry, citems)

# Specification of multi-group SEM consisting of 5 measurement models and 1 structural model
sem5 <- '
  # 5 measurement models: five latent variables: fixing first loading of each concept to 1
  social_trust        =~ 1*ppltrst + pplfair + pplhlp
  negative_gays       =~ 1*freehms + hmsfmlsh_r + hmsacld
  positive_immigrants =~ 1*imbgeco + imueclt + imwbcnt
  religiousness       =~ 1*rlgdgr + rlgatnd_r + pray_r
  gender_equality     =~ 1*eqwrkbg + eqpolbg + eqmgmbg + eqpaybg

  # 1 structural model: explain negative attitudes towards gays/lesbians
  negative_gays ~ social_trust + positive_immigrants + religiousness + gender_equality
'

# 1) 
# Configural measurement invariance model with country-specific regression coefficients
# Fit configural model: loadings and regression coefficients are freely estimated in BE and PL
# group = "country" >>> multi-group SEM
fitsem_config_1 <- sem(
  sem5,
  data = semdata,
  group = "country"
)

# Display model results and fit measures
summary(fitsem_config_1, fit.measures = TRUE)
# Display standardized estimates for both countries
#standardizedSolution(fitsem_config_1)
# Display main fit measures
#fitmeasures(fitsem_config_1, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# 2)
# Configural measurement invariance model with equal regression coefficients across BE and PL
# Factor loadings are freely estimated, but regression coefficients are constrained to be equal across BE and PL
fitsem_config_2 <- sem(
  sem5,
  data = semdata,
  group = "country",
  group.equal = "regressions"
)

# Display model results and fit measures
summary(fitsem_config_2, fit.measures = TRUE)

# 3)
# Metric measurement invariance model with country-specific regression coefficients
# Factor loadings are constrained to be equal across BE and PL.
# Regression coefficients are freely estimated in each country.
fitsem_metric_3 <- sem(
  sem5,
  data = semdata,
  group = "country",
  group.equal = "loadings"
)

# Display model results and fit measures
summary(fitsem_metric_3, fit.measures = TRUE)


# 4)
# Metric measurement invariance model with equal regression coefficients across BE and PL
# Factor loadings AND regression coefficients are both constrained to be equal across BE and PL.
fitsem_metric_4 <- sem(
  sem5,
  data = semdata,
  group = "country",
  group.equal = c("loadings", "regressions")
)

# Display model results and fit measures
summary(fitsem_metric_4, fit.measures = TRUE)


# Comparing the fit measures of the four estimated models 
# The individual LR tests compare each model against the saturated model.
# Problem: They are sensitive to sample size
# To compare all four models simultaneously, we have to use AIC, BIC, CFI, TLI, RMSEA and SRMR.
# Because not all four models are nested.
fitmeasures_list <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
rbind(
  Model1 = fitMeasures(fitsem_config_1, fitmeasures_list),
  Model2 = fitMeasures(fitsem_config_2, fitmeasures_list),
  Model3 = fitMeasures(fitsem_metric_3, fitmeasures_list),
  Model4 = fitMeasures(fitsem_metric_4, fitmeasures_list)
)

# Comparing the fit measures of the four estimated models 

# In addition, comparison tests where possible to select the best model
# Models 2) and 3) are NOT nested, so they cannot be compared using a standard LR test.
# Hence, to compare all 4 models simultaneaously, we have to use the fit measures reported above.

# Models 1) and 2) are nested: Model 2) constrains the four regression coefficients
# to be equal across BE and PL, whereas Model 1) estimates them freely.
# Both models estimate the factor loadings freely across countries.
# LR test examines whether these equality constraints significantly worsen model fit.
anova(fitsem_config_1, fitsem_config_2)

# Models 1) and 3) are nested: Model 3) constrains the factor loadings
# to be equal across BE and PL, whereas Model 1) estimates them freely.
# LR test examines whether these equality constraints significantly worsen model fit.
anova(fitsem_config_1, fitsem_metric_3)

# Models 1) and 4) are nested: Model 4) constrains both factor loadings
# and regression coefficients to be equal across BE and PL,
# whereas Model 1) estimates both freely.
# LR test examines whether these equality constraints significantly worsen model fit.
anova(fitsem_config_1, fitsem_metric_4)

# Models 2) and 4) are nested: Model 4) constrains the factor loadings
# to be equal across BE and PL, whereas Model 2) estimates them freely.
# Both models constrain the regression coefficients to be equal across countries.
# LR test examines whether these additional equality constraints significantly worsen model fit.
anova(fitsem_config_2, fitsem_metric_4)

# Models 3) and 4) are nested: Model 4) constrains the four regression coefficients
# to be equal across BE and PL, whereas Model 3) estimates them freely.
# Both models constrain the factor loadings to be equal across countries.
# LR test examines whether these additional equality constraints significantly worsen model fit.
anova(fitsem_metric_3, fitsem_metric_4)


################################
# Results of best model: model 1
################################

# Display model results and fit measures
summary(fitsem_config_1, fit.measures = TRUE)

# Main fit measures
fitMeasures(fitsem_config_1, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr"))

# Standardized solution for both countries
standS_model1 <- standardizedSolution(fitsem_config_1)
options(max.print = 10000)
print(standS_model1)

# Extract R sqaured for negative attitudes towards gays and lesbians by country 
# The following line of code was generated using GenAI. 
sapply(lavInspect(fitsem_config_1, "rsquare"), function(x) x["negative_gays"])


# Composite reliability
# Function compositereliability from task 2b used here

# Extract standardized factor loadings from Model 1
loadings_sem <- subset(standS_model1, op == "=~")

# Calculate composite reliability for each factor in Belgium (group 1)
reliab_BE <- c(social_trust = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "social_trust" & loadings_sem$group == 1]),
  negative_gays = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "negative_gays" & loadings_sem$group == 1]),
  positive_immigrants = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "positive_immigrants" & loadings_sem$group == 1]),
  religiousness = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "religiousness" & loadings_sem$group == 1]),
  gender_equality = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "gender_equality" & loadings_sem$group == 1]))

# Calculate composite reliability for each factor in Poland (group 2)
reliab_PL <- c(social_trust = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "social_trust" & loadings_sem$group == 2]),
  negative_gays = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "negative_gays" & loadings_sem$group == 2]),
  positive_immigrants = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "positive_immigrants" & loadings_sem$group == 2]),
  religiousness = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "religiousness" & loadings_sem$group == 2]),
  gender_equality = compositereliability(loadings_sem$est.std[loadings_sem$lhs == "gender_equality" & loadings_sem$group == 2]))

# Display composite reliabilities for Belgium and Poland
round(cbind(Belgium = reliab_BE, Poland = reliab_PL), 3)






