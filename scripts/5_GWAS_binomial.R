# Title: "Binomial GWAS of binary response"
# Author: Aglaia Szukala
# Short description: Script to run different GWAS models using adegenet and glmnet when ADB damage is coded as a 0/1 trait

# Load libraries
library(nnet)
library(glmnet)
library(adegenet)
library(ggplot2)

# Set working directory
setwd("//bfwsbfile1/Institut2/Aglaia/EscheNot/manuscript/SubmissionDataScripts/data/")

##########               ##########
########## GENOTYPE DATA ##########
##########               ##########

gen <- read.table("4_NewRefAnn_UniqueHits_NoHybridInd_GWAS+LD/Dat_6742SNPs_1037Ind.txt",
                  sep="\t", 
                  header=T, 
                  check.names = F)

gen <- as.matrix(gen)
gen[1:5,1:5]
dim(gen)

# Codify missing values as NA
gen[gen == -1] <- NA

# Estimate amount of NAs in dataset
(sum(is.na(gen)) * 100) / (dim(gen)[1] * dim(gen)[2]) # 0,29% missing data

##########                     ##########
########## IMPUTE MISSING DATA ##########
##########                     ##########

# Impute missing data using the most common genotype at each SNP across all individuals.
gen.imp <-
  apply(gen, 2, function(x)
    replace(x, is.na(x), as.numeric(names(which.max(
      table(x)
    )))))
sum(is.na(gen.imp)) # No NAs
dim(gen.imp)
dim(gen)

##########                ##########
########## PHENOTYPE DATA ##########
##########                ##########

phen <- read.table("5_PhenotypeClimDat/Dat_phenotype_climate.txt")

# Keep only samples as in genotype matrix
phen <- phen[match(colnames(gen.imp),phen$Taxa),]
dim(phen)
head(phen)

# Check that samples Ids are in the same order in phenotype and genotype matrices
phen$Taxa <- as.character(phen$Taxa)
identical(colnames(gen.imp),phen$Taxa) # all good

# Convert tolerance and flushing phenotypes to factors
str(phen)
phen$Tolerance <- as.factor(phen$Tolerance)
phen$Flush_1st <- as.factor(phen$Flush_1st)
phen$Flush_2nd <- as.factor(phen$Flush_2nd)
phen$Field_test <- as.factor(phen$Field_test)

###########                        ##########
########### GWAS with Binary Model ##########
###########                        ##########

# Interpret ADB damage as a binary phenotype. Two cases:
# 1) Take only class 1 (healthy) and 6 (extremely sick)
# 2) Take classes 1 and 2 as healthy (no ADB = 0) and 5 and 6 as sick (with ADB = 1)
# Test both cases

# Build case 1 phenotype data
binaryDat <- phen[phen$Tolerance == "1" | phen$Tolerance == "6",] # keep only class 1 and 6
dim(binaryDat)
head(binaryDat)

# Assign 0 to healthy and 1 to sick by making a new ADB variable
binaryDat$ADB <- "0" # healthy
binaryDat[binaryDat$Tolerance == "6",]$ADB <- "1" # sick
table(binaryDat$Tolerance) # check frequency of each class in tolerance
table(binaryDat$ADB) # check frequency of classes 1 and 0 

# Make dataframe for phenotype
Dat_ADB01 <- as.data.frame(binaryDat[,"ADB"])
head(Dat_ADB01)
rownames(Dat_ADB01) <- binaryDat$Taxa
colnames(Dat_ADB01) <- "ADB"


# write.table(Dat_ADB01,file="5_PhenotypeClimDat/Dat_ADB_binary.txt",
#             col.names = T,
#             row.names =T,
#             quote = F)
# write.table(rownames(Dat_ADB01), file="5_PhenotypeClimDat/Dat_ADB_binary_taxaList.txt",
#             quote = F,
#             row.names = F,
#             col.names = F)

# Build case 2 phenotype data
binaryDat <- phen[phen$Tolerance == "1" | phen$Tolerance == "2" | phen$Tolerance == "5" | phen$Tolerance == "6",]
dim(binaryDat)
head(binaryDat)
binaryDat$ADB <- "0" # healthy
binaryDat[binaryDat$Tolerance == "6" | binaryDat$Tolerance == "5",]$ADB <- "1" # sick
table(binaryDat$Tolerance)
table(binaryDat$ADB)
Dat_ADB01 <- as.data.frame(binaryDat[,"ADB"])
head(Dat_ADB01)
rownames(Dat_ADB01) <- binaryDat$Taxa
colnames(Dat_ADB01) <- "ADB"

# Save these datasets for analysis in Tassel
# write.table(Dat_ADB01,file="5_PhenotypeClimDat/Dat_ADB_binary_case2.txt",
#             col.names = T,
#             row.names =T,
#             quote = F)
# write.table(rownames(Dat_ADB01), file="5_PhenotypeClimDat/Dat_ADB_binary_case_2_taxaList.txt",
#             quote = F,
#             row.names = F,
#             col.names = F)


###########                                            ##########
########### GWAS without and with population structure ##########
###########                                            ##########

# Keep only individual genotypes that are present in binary phenotype table
genADB01 <- gen.imp[, match(rownames(Dat_ADB01),colnames(gen.imp))]
genADB01[1:5,1:5]
dim(genADB01)
genADB01 <- t(genADB01)

# Prepare Data for adegenet and glmnet
X <- as.matrix(genADB01) # genotype matrix
y <- Dat_ADB01$ADB # phenotype matrix

####                                                      ####
#### First run GWAS without correction for stratification ####
####                                                      ####

# 1) Naive association test using a univariate method (Fisher’s exact test)
pval <- apply(X, 2, function(e)
  fisher.test(table(factor(e, levels=c(0,1)), y))$p.value
  )

# Check p-val of most significant SNP and see how many SNPs are below 0.05 before FDR correction
min(pval)
length(which(pval < 0.05))

# Multiple testing correction
pval.corrected.bonf <- p.adjust(pval, 
                                method="bonferroni"
                                ) # bonferroni correction
pval.corrected.fdr <- p.adjust(pval, 
                               method="fdr"
                               ) # pval correction

# How many SNPs are retained?
length(which(pval.corrected.fdr<0.05)) # no significant snp
length(which(pval.corrected.bonf<0.05))

# Produce Manhattan plot
log.pval <- -log10(pval.corrected.fdr)
set.seed(1)
log.pval <- jitter(log.pval, amount=0.5)
plot(log.pval,
     col = transp(azur(5)),
     pch = 19,
     cex = 1.5,
     ylim=c(-0.5, 5),
     main="Fisher's exact test \n(FDR correction)",
     xlab="SNP loci", 
     ylab="FDR-corrected -log10(p-value)"
      )
thresh <- -log10(0.05)
abline(h=thresh, col = "red", lwd=2) # nothing interesting at all...

# 2) Multivariate method: Fit Lasso multinomial logistic regression model
model_glmnet_bin <- cv.glmnet(X, 
                              y, 
                              family = "binomial", # the family argument allows to set a probability distribution
                              alpha = 1,
                              lambda.min.ratio=0.01
) 
plot(model_glmnet_bin)

# The LASSO method generates coefficients for each variable, though the majority of these
# will have been shrunk to zero in the penalization step. We extract these coefficients from the
# model generated by cv.glmnet with the function coef.
beta <- as.vector(t(coef(model_glmnet_bin, 
                         s="lambda.min")
                    )
                  )

# We retrieve the results of the implicit feature selection performed by LASSO by identifying
# the variables that have non-zero coefficients.
res <- which(beta[-1] !=0)
length(res) # no significant SNPs

# Code to extract significant SNPs if I would have some
coefs.model <- beta[-1][res]
# names(coefs.model) <- colnames(snps)[res]

# A standard visualisation tool for the LASSO method is a plot of the fraction of deviance
# explained versus the values of the coefficients retained in the LASSO regression model.

fit <- model_glmnet_bin$glmnet.fit
# Improved plot compared to default
y.pos <- coefs.model-c(0.2, 0.25, 0.5, 0.75)
plot(fit, 
     xvar = "dev", 
     label = FALSE
     )
text(x=1.01, 
     y=y.pos,
     labels=names(coefs.model), 
     col="black", 
     pos=2, 
     cex=0.6
     )
grid()
title("Fraction of deviance explained by LASSO coefficients", 
      line=3
      )

####                                        ####
#### Run GWAS correcting for stratification ####
####                                        ####

# Run PCA and compute euclidean distances between individuals of the trimmed genotypes 
# to see how many PCs I should include for population structure

# PCA of genotypes - plots barplot of eigenvalues automatically
pca <- dudi.pca(X, 
                scale=F
                ) # 7 axes showing larger increments were selected
head(pca$eig)
head(pca$l1)

# Plot PCA
s.label(pca$li, 
        sub="PCA - PC 1 and 2"
        )
add.scatter.eig(pca$eig,5,1,2, 
                ratio=.26, posi="topleft"
                )

# Add a quantitative assessment of individuals clustering via squared Euclidean distances
# between individuals (function dist) and hierarchical clustering with complete linkage
# (hclust) to define tight clusters:
D <- dist(pca$li)^2
clust <- hclust(D, 
                method="complete"
                )

# Plot the distances stored in the dist object D in a heatmap
temp <- as.data.frame(as.matrix(D))
temp <- t(as.matrix(D))
temp <- temp[,ncol(temp):1]
par(mar=c(5.1,4.1,4.1,2.1))
image(x=1:length(rownames(X)), 
      y=1:length(rownames(X)), 
      temp, 
      col=rev(heat.colors(nlevels(as.factor(D)))),
      xaxt="n", 
      yaxt="n",
      xlab="",
      ylab="")
axis(side=2, 
     at=1:length(rownames(X)), 
     lab=rev(rownames(X)), 
     las=2, 
     cex.axis=.46)
axis(side=3, 
     at=1:length(rownames(X)), 
     lab=rev(rownames(X)),  
     las=2, 
     cex.axis=.46)
title("Genetic distances between isolates", 
      outer=TRUE, 
      line=-1)

# Define clusters 
pop <- factor(cutree(clust, 
                     k=7
                     )
              ) #attribute individuals to each one of the chosen clusters
table(pop) # most individuals are in one cluster with a few "outliers", but remeber that the overall amount of explined variance by the first PCs after trimming the data was rather low

# Plot clusters on PCA
s.class(pca$li, 
        fac=pop, 
        col=transp(funky(7)), 
        cpoint=2,
        sub="PCA - axes 1 and 2"
        )
add.scatter.eig(pca$eig,5,1,2, ratio=.18, posi="topleft")

# Perform DAPC-based feature selection to include structure in the models
# Begin the DAPC approach by running cross-validation to help
# select the number of PCs of PCA to retain that will maximize the ability to discriminate
# between the two phenotypic groups.
set.seed(1)
xval1 <- xvalDapc(X, y,
                  n.pca=c(c(1,3), 
                          seq(10, 80, 10)),
                  n.rep=20) # may take a moment...

# Let’s take a look at the object xval1 containing the results of cross-validation
xval1[2:6] # 1 PC seem to be enough to correct phenotypic groups, so basically no structure, but let´s still try with 7 PCs

# The last element of the output of xvalDapc is a dapc object generated with the optimal
# number of PCs, as indicated by RMSE. Store this in an object called dapc1
dapc1 <- xval1[[7]]

# Use the function snpzip to perform feature selection and visualize results; 
# Visualize a discriminant function plot showing the separation of phenotypes based on SNPs
result <- snpzip(X, 
                 dapc1,
                 method="ward", 
                 xval.plot = FALSE,
                 plot = TRUE, 
                 loading.plot = TRUE
                 )
# Clearly, SNPs cannot distinguish the phenotype classe, which overlap
summary(dapc1) # no healthy individual was assigned correctly as expected

# Include population stratification in the GWAS model
# Create a snps dataset corrected for population structure using k components
snps.corrected <- apply(X, 
                        2, 
                        function(e) 
                          residuals(lm(e~pca$li[,1]+pca$li[,2]+pca$li[,3]+pca$li[,4]+pca$li[,5]+pca$li[,6]+pca$li[,7] # add less or more components for testing
               )
               )
               )
# Let’s inspect The corrected SNPs matrix:
dim(snps.corrected)
range(snps.corrected)

# To visually assess whether our correction for population stratification has been successful,
# we can run a second PCA analysis, this time with the corrected SNPs matrix:
pca2 <- dudi.pca(snps.corrected, scale=FALSE, scannf=FALSE, nf=7)
barplot(pca2$eig, main="PCA eigenvalues")
s.class(pca2$li, 
        fac=pop, 
        col=transp(funky(7)), 
        cpoint=2,
        sub="PCA - axes 1 and 2"
)
# The data looks well adjusted

# GWAS
# Instead of Fisher’s exact test use an alternative univariate approach that consists of
# two stages. First, we generate a simple linear model between each column of our corrected
# SNPs matrix and the phenotypic trait. Second, run an analysis of variance (ANOVA) on
# each model generated, specifying a Chi-squared test of association. From this, retrieve
# a p-value for the significance of association between each corrected SNP and the resistance
# phenotype.

pval2 <- numeric(0)
y <- as.numeric(y) # needs numeric phenotype
head(y)
for(i in 1:ncol(snps.corrected)){
  foo <- suppressWarnings(glm(y ~ snps.corrected[,i], 
                              family="binomial"))
  ANOVA <- anova(foo, test="Chisq")
  pval2[i] <- ANOVA$"Pr(>Chi)"[2]
}

min(pval2, na.rm=TRUE)
# 0.0002735301
length(which(pval2 < 0.05)) # previous to FDR correction
pval.corrected.fdr <- p.adjust(pval2, method="fdr")
res <- which(pval.corrected.fdr < 0.05)
length(res) # no significant SNPs

# Multivariate method
# Fit Lasso multinomial logistic regression model for many snps
model_glmnet_bin <- cv.glmnet(snps.corrected, y, 
                              family = "binomial", 
                              alpha = 1,
                              lambda.min.ratio=0.01
) # the family argument allows us to set a probability distribution!
plot(model_glmnet_bin)

# The LASSO method generates coefficients for each variable, though the majority of these
# will have been shrunk to zero in the penalization step. We extract these coefficients from the
# model generated by cv.glmnet with the function coef.
beta <- as.vector(t(coef(model_glmnet_bin, s="lambda.min")))

# We retrieve the results of the implicit feature selection performed by LASSO by identifying
# the variables that have non-zero coefficients.
res <- which(beta[-1] !=0)
length(res) # no significant SNPs

# DAPC
result <- snpzip(snps.corrected, as.factor(y),
                 xval.plot=TRUE, plot=TRUE, loading.plot=TRUE,
                 method="ward")

# Best lambda identified through cross-validation
best_lambda <- model_glmnet_bin$lambda.min
print(best_lambda)

############ Check this ##################
# Fit the final model with the best lambda
final_model <- glmnet(X, y, family = "binomial", alpha = 1, lambda = best_lambda)
print(final_model)

# Extract the coefficients (SNPs with non-zero coefficients are considered significant)
coefficients <- coef(final_model, s = best_lambda)

# Convert the coefficients to a matrix for easy handling
coefficients_matrix <- as.matrix(coefficients) #class 1
head(coefficients_matrix)

# Identify SNPs with non-zero coefficients
significant_snps <- rownames(coefficients_matrix)[coefficients_matrix != 0]
significant_snps <- significant_snps[significant_snps != "(Intercept)"] # Exclude the intercept
significant_snps # Keine significant SNPs using binomial model, The SNPs significant in the multinomial jsut differentiate class 1 and 6 from the rest




