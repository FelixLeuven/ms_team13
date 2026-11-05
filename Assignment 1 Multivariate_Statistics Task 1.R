load("agdata.Rdata")

#a 

#reversing hmsfmlsh, rlgatnd and pray such that, within each concept, a higher score indicates more of the concept it measures (negative attitude towards gays/lesbians, religiousness)
agdata$hmsfmlsh <- 6 - agdata$hmsfmlsh
agdata$rlgatnd <- 8 - agdata$rlgatnd
agdata$pray    <- 8 - agdata$pray

#select the 16 items, so that no other columns enter the analysis
items = agdata[, c("ppltrst", "pplfair", "pplhlp",
                   "freehms", "hmsfmlsh_r", "hmsacld",
                   "imbgeco", "imueclt", "imwbcnt",
                   "rlgdgr", "rlgatnd_r", "pray_r",
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