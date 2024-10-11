# Title: "Multinomial GWAS of categorical response"
# Author: Aglaia Szukala
# Short description: Script to run GWAS models using adegenet and glmnet when ADB damage is coded as a 1 to 6 categories
# First and second flushing were also analyzed the same way

# Load libraries
library(nnet)
library(glmnet)
library(adegenet)
library(ggplot2)

# Set working directory
setwd("./data/")

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

# Check values in flush variable
table(phen$Flush_2nd) 
table(phen$Flush_1st)
table(phen$Tolerance)

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
head(phen)

#############                                             #############
############# Fit a multinomial logistic regression model #############
#############                                             #############

# Build data
genotype <- t(gen.imp)
genotype <- as.data.frame(genotype)
genotype[1:5,1:5]
dim(genotype)

# Preparing the data for glmnet: choose the one to be analyze

#2) 1st Flushing
phenfl <- phen[!is.na(phen$Flush_1st),] #remove NAs
dim(phenfl)
head(phenfl)
head(phen)
phenfl2 <- phenfl[phenfl$flushClass == "green",] # keep only correct assignments that are in the "green" class
dim(phenfl2)
head(phenfl2)
genotype<- genotype[rownames(genotype) %in% phenfl2$Taxa,]
X <- as.matrix(genotype)
dim(genotype)
dim(X)
y <- as.data.frame(phenfl2[,"Flush_1st"])
rownames(y) <- phenfl2$Taxa
colnames(y) <- "Flush_1st"
dim(y)
#or
y <- as.data.frame(phenfl2[,"Flush_2nd"])
rownames(y) <- phenfl2$Taxa
colnames(y) <- "Flush_2nd"

###########                                          ##########
########### GWAS accounting for population structure ##########
###########                                          ##########

# Run PCA and compute euclidean distances between individuals of the trimmed genotypes 
# to see how many PCs I should include for population structure

# PCA of genotypes - plots barplot of eigenvalues automatically
pca <- dudi.pca(X, 
                scale=F
) # from a first check I saw that here it makes more sense to keep only 2 axes showing larger increments were selected
head(pca$eig)
head(pca$l1)

# Plot PCA
s.label(pca$li, 
        xax = 1, 
        yax = 7, # check from 2 to 7; at PC 7 there are no clear separation anymore
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
                     k=6 # PCs retained in pca
)
) #attribute individuals to each one of the chosen clusters
table(pop) # cluster 2 contains the western individuals

# Plot clusters on PCA
s.class(pca$li, 
        fac=pop,
        xax = 1, 
        yax = 4,
        col=transp(funky(6)), 
        cpoint=2,
        sub="PCA - axes 1 and 2"
)
add.scatter.eig(pca$eig,5,1,2, 
                ratio=.18, 
                posi="topleft"
                ) # looks ok to me

# Perform DAPC-based feature selection to include structure in the models
# Begin the DAPC approach by running cross-validation to help
# select the number of PCs of PCA to retain that will maximize the ability to discriminate
# between the two phenotypic groups.
set.seed(1)
xval1 <- xvalDapc(X, 
                  as.matrix(y),
                  n.pca=c(c(1,3),seq(10, 80, 10)),
                  n.rep=20) 
        
# Take a look at the object xval1 containing the results of cross-validation
xval1[2:6] # 3 PCs seem optimal to correct phenotypic groups, but let´s still try with values up to 7 PCs

# The last element of the output of xvalDapc is a dapc object generated with the optimal
# number of PCs, as indicated by RMSE. Store this in an object called dapc1
dapc1 <- xval1[[7]]

# Use the function snpzip to perform feature selection and visualize results; 
# Visualize a discriminant function plot showing the separation of phenotypes based on SNPs
result <- snpzip(X, 
                 dapc1,
                 method="ward.D", # or "ward" 
                 xval.plot = FALSE,
                 plot = TRUE, 
                 loading.plot = TRUE
)

summary(dapc1)

###########                                               ##########
########### Include population stratification in the GWAS ##########
###########                                               ##########

# Include pop stratification
# Include population stratification in the model by correcting the snps matrix for structure
# Apply the residuals correction to each SNP column in the SNP matrix X
# Linear model regresses the SNP (e) on the first k principal components from the DAPC or PCA results
# Residuals() extracts the residuals from the linear model. 
# These residuals are essentially the SNP values after the population structure (represented by the PCA components) has been accounted for. 

snps.corrected <- apply(X, 2, function(e)
  residuals(lm(e~pca$li[,1]+pca$li[,2]+pca$li[,3] #+pca$li[,4]+pca$li[,5]+pca$li[,6]
               ))) # try with 2 to 7

# Let’s inspect our corrected SNPs matrix:
dim(snps.corrected)
snps.corrected[1:5,1:5]
range(snps.corrected)

# To visually assess whether our correction for population stratification has been successful,
# we can run a second PCA analysis, this time with the corrected SNPs matrix:
pca2 <- dudi.pca(snps.corrected, 
                 scale=FALSE, 
                 scannf=T
                 )
s.class(pca2$li, 
        xax = 1, 
        yax = 3,
        fac=pop, 
        col=transp(funky(7)), 
        cpoint=2,
        sub="PCA - axes 1 and 2"
) # correction works well

# Run the GWAS model on the corrected snps matrix
# Run cross-validation and plot the returned object
cvfit <- cv.glmnet(snps.corrected,
                   y$Flush_2nd, 
                   family = "multinomial",
                   alpha = 1
                   )

# Extract Best lambda from the cross-validation
best_lambda <- cvfit$lambda.min
best_lambda

# Fit the final model with the best lambda
final_model <- glmnet(snps.corrected, 
                      y$Flush_2nd, 
                      family = "multinomial", 
                      alpha = 1, 
                      lambda = best_lambda
                      )
#plot(final_model)

# Extract the coefficients for the model at the best lambda
coefficients <- coef(final_model, s = best_lambda) # coefficients is a list where each element corresponds to a category of the outcome

# Initialize an empty vector to store significant SNPs
significant_snps <- c()

# Loop over each category of the response variable
for (i in 1:length(coefficients)) {
  # Get the coefficient matrix for this category
  coef_matrix <- as.matrix(coefficients[[i]])
  
  # Find SNPs (or PCs) with non-zero coefficients
  non_zero_indices <- which(coef_matrix != 0)
  
  # Get the SNP names (or PCs)
  snp_names <- rownames(coef_matrix)[non_zero_indices]
  
  # Exclude intercept term from significant SNPs
  snp_names <- snp_names[snp_names != "(Intercept)"]
  
  # Store the SNPs
  significant_snps <- unique(c(significant_snps, snp_names))
}

# significant_snps now contains the list of SNPs with non-zero coefficients
print(significant_snps)
  
