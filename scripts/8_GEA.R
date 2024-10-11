# Title: "Lfmm and RDA"
# Author: Aglaia Szukala
# Short description: Script for latent factor mixed model and redundancy analysis of climatic variables 

# Load libraries
library(lfmm)
library(readxl)
library(factoextra)
library(qvalue)
library(SNPRelate)
library(vegan)
library(GWASTools)
library(lattice)
library(psych)
library(gplots)
library(vcd)
library(qqman)
library(dplyr)

############                                                                ###########
########################## LOAD and PRPARE GENOTYPE DATA ##############################
############                                                                ###########

setwd("./data")

# Genotype data after filtering by maf, missingness, and monomorphic loci
gen <- read.table("4_NewRefAnn_UniqueHits_NoHybridInd_GWAS+LD/Dat_6742SNPs_1037Ind.txt",
                  sep="\t", 
                  header=T, 
                  check.names = F)
gen <- as.matrix(gen)

# Codify missing values as NA
gen[gen == -1] <- NA

# Check amount of missing data
sum(is.na(gen))
(sum(is.na(gen)) * 100) / (dim(gen)[1] * dim(gen)[2]) # 0,29% missing data

# I impute using the most common genotype at each SNP across all individuals.
gen.imp <-
  apply(gen, 2, function(x)
    replace(x, is.na(x), as.numeric(names(which.max(
      table(x)
    )))))
sum(is.na(gen.imp)) # No NAs
dim(gen.imp)

############                                                                      ###########
########################## LOAD and PREPARE ENVIRONMENTAL DATA ##############################
############                                                                      ###########

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

# Predictors must be encoded as numeric
for (i in 9:14) {
  phen[,i] <- as.numeric(phen[,i])
}
head(phen)
str(phen)

# Remove 4 samples for which I am missing coordinates and environmental data
env <- phen[-c(491,493,872,874),]
dim(phen)
dim(env)
gen.imp <- gen.imp[,env$Taxa]
dim(gen.imp)

# Check distribution of environmental predictors
par(mfrow=c(3,2))
for (i in 9:14) {
  hist(na.omit(env[,i]), main = paste0("Distribution of ",colnames(env[i])),
       xlab = paste0(colnames(env[i])))
  #lines(density(na.omit(env[,i])), col="black", lwd=2)
}

# Perform PCA on the predictor variables
pred <- env[, c(9:14)]
head(pred)
dim(pred)
pred.pca <- rda(pred, scale = T) # rda function in vegan
summary(pred.pca)$cont
par(mfrow=c(1,1))
screeplot(pred.pca, main = "Eigenvalues of predictors")

# Correlation between the PC axis and predictors:
round(scores(
  pred.pca,
  choices = 1:3,
  display = "species",
  scaling = 0
),
digits = 3) # first PC correlates with Precipitation, T and Altitude to a similar extent

# Store the synthetic PC axis predictor as pred.PC for use in LFMM
pred.PC <- scores(pred.pca,
                   choices = 1:6,
                   display = "sites",
                   scaling = 0)
head(pred.PC)

############                                       ###########
################### LFMM ridge estimates #####################
############                                       ###########


# Run pca with broken sticks to confirm that K=3 best explains the data
gen.pca <- rda(gen.imp, scale = T)
screeplot(
  gen.pca,
  main = "Screeplot of Genetic Data with Broken Stick",
  bstick = TRUE,
  type = "barplot",
  npcs = 10
)

# retrieve data with chromosome and position info for manhattan plot
snps_info <- read.table(file="4_NewRefAnn_UniqueHits_NoHybridInd_GWAS+LD/Dat_SNPsList_ChrPos.txt")
colnames(snps_info) <- c("chrom", "pos","snpID")
dim(snps_info)
snps_info$chromN <- as.numeric(factor(snps_info$chrom, levels = unique(snps_info$chrom),
                                      labels = c(1:24)))
head(snps_info)

# Define latent factor mixed models
run.lfmm = function(y, x, k) { #add g if you want to add a new lambda parameter
  frax.lfmm <- lfmm_ridge(Y = y, X = x, K = k, it.max = 1000000) # define model using ridge estimates
  #frax.lfmm <- lfmm_lasso(Y = y, X = x, K = k,nozero.prop = 0.01) # define model using lasso estimates
  frax.pv <-
    lfmm_test(
      Y = y,
      X = x,
      lfmm = frax.lfmm,
      calibrate = "gif"
    ) # get test statistics of predictor(s)
  
  print(paste0("this are the calculated GIF: ", frax.pv$gif)) # print out to screen the GIF
  
  # plot histograms of p-values before and after corrections: needs the GIF to be manually readjusted?
  par(mfrow = c(1, 2))
  hist(frax.pv$pvalue, main = "Unadjusted p-values", xlab = "p-value")
  hist(
    frax.pv$calibrated.pvalue[, 1],
    main = paste0("GIF-adjusted p-values (GIF=", round(frax.pv$gif, 2), ")"),
    xlab = "p-value"
  )
  
  # obtain z-scores for health
  zscore <-
    frax.pv$score[, 1]   # zscores for first predictor: health (therefore always chekc that health is the first predictor in the X matrix, otherwise you get outliers for the other predictors)
  
  # Use your own gif values for recalibration
  #new.gif1 <- g
  
  # calculate new calibrated p-values using new GIF and plot histogram
  # adj.pv1 <- pchisq(zscore ^ 2 / new.gif1, df = 1, lower = FALSE)
  # hist(
  #   adj.pv1,
  #   main = paste0("REadjusted p-values (GIF=", new.gif1, ")"),
  #   xlab = "adjusted p-values"
  # )
  
  # apply FDR correction
  frax.qv <- qvalue(frax.pv$calibrated.pvalue)$qvalues
  print(paste0("Number of SNPs with FDR<0.05: ", length(which(frax.qv < 0.05))))
  print(paste0("Number of SNPs with calib p-val<0.00001: ", length(
    which(frax.pv$calibrated.pvalue < 0.00001)
  )))
  
  # qqplots of p-values
  qqplot(
    rexp(length(frax.pv$pvalue), rate = log(10)),-log10(frax.pv$pvalue),
    xlab = "Expected quantile",
    pch = 19,
    cex = .4,
    main = "QQplot of raw p-values"
  )
  abline(0, 1)

  qqplot(
    rexp(length(frax.pv$calibrated.pvalue), rate = log(10)),-log10(frax.pv$calibrated.pvalue),
    xlab = "Expected quantile",
    pch = 19,
    cex = .4,
    main = "QQplot of calibrated p-values"
  )
  abline(0, 1)
  
  # identify which SNPs these are
  frax.FDR.1 <- colnames(y)[which(frax.qv < 0.05)]
  frax.PV.01 <- colnames(y)[which(frax.pv$calibrated.pvalue < 0.0001)]
  
  # plot manhattan plot
  snps_sub <- snps_info[match(rownames(frax.pv$calibrated.pvalue), snps_info$snpID), ]
  snps_sub$p <- frax.pv$calibrated.pvalue[,1] 
  
  #par(mfrow=c(1,2))
  manhattan(
    na.omit(snps_sub[!snps_sub$chromN == 24,]), # this is an not assembled fragment that we can leave out for plotting purpose
    chr = "chromN",
    bp = "pos",
    snp = "snpID",
    p = "p",
    col = c("#ffbc42", "#218380"),
    annotateTop = T,
    genomewideline = F,
    suggestiveline = F,
    #ylim=c(0,6),
    xaxt='n' #,
    #xlab=""
  )
  abline(h = -log10(0.001), col="#780000",lty = 1)
  (GWAS_Bonn_corr_threshold <- -log10(0.05 / nrow(frax.pv$calibrated.pvalue)))
  abline(h = GWAS_Bonn_corr_threshold, col="#780000",lty = 3)
  
  qq(snps_sub$p)
  
  return(list(frax.FDR.1, frax.PV.01))#frax.GIFadjQV.1, 
}

# Run lfmm for annual precipitation, with different K values
run.lfmm(t(gen.imp), pred$AnnPrec,1) # k=1 does not correct well for structure
run.lfmm(t(gen.imp), pred$AnnPrec,3)
run.lfmm(t(gen.imp), pred$AnnPrec,7) # K= 7 captures better the fine structure 
run.lfmm(t(gen.imp), log(pred$AnnPrec),7)

# Use synthetic PC and other predictors 
run.lfmm(t(gen.imp), pred.PC[,1:3],3)
run.lfmm(t(gen.imp), pred$Altitude,3)
run.lfmm(t(gen.imp), log(pred$Altitude),3)
run.lfmm(t(gen.imp), log(pred$Altitude),7)
run.lfmm(t(gen.imp), pred$AnnMeanT,7)


#################                  #################
################# Run RDA analysis #################
#################                  #################

# Perform PCA of genotype matrix
pc_filt <- prcomp(t(gen.imp))

# Add PC 1 to 3 of genotype data as structure proxy
env$PC1 <- pc_filt$x[,1]
env$PC2 <- pc_filt$x[,2]
env$PC3 <- pc_filt$x[,3]
env$PC4 <- pc_filt$x[,4]
env$PC5 <- pc_filt$x[,5]
env$PC6 <- pc_filt$x[,6]
head(env)

# define predictors you wish to use
pred.red <- env[, c(9:11,16:21)]
pred.red$syntheticPC <- pred.PC[,1]
head(pred.red)

# visualize correlations among predictors
pairs.panels(pred.red, scale=F) 

# RDA model
# 1) Simplest model without structure
frax.rda <- 
  rda(t(gen.imp) ~ ., 
      data=pred.red[,c(1:3)], 
      scale=T
      )

# 2) Conditioned RDA, the effects of (structure and) is partialled out
frax.rda <- 
  rda(t(gen.imp) ~ 
        `AnnPrec`+
        `AnnMeanT`+
        `Altitude`+
        Condition(PC1)+
        Condition(PC2),#test different numbers of PCs
      data=pred.red, 
      scale=T
      )

# 3) Final model: Conditioned RDA, the effects of (structure and) is partialled out, + use the synthetic PC as predictor due to collinearity among predictors
frax.rda <- 
  rda(t(gen.imp) ~ 
        `syntheticPC`+
        Condition(PC1)+
        Condition(PC2)+
        Condition(PC3)+
        Condition(PC4)+
        Condition(PC5)+
        Condition(PC6),
      data=pred.red, 
      scale=T
      )

# Extract the equivalent of regression coefficients for each 
# explanatory variable on each canonical axis
coef(frax.rda)

# Calculate the adjusted R2 using:
RsquareAdj(frax.rda) # The constrained ordination explains about 0.1% of variation

# The eigenvalues for the constrained axes reflect the variance explained by each canonical axis:
summary(frax.rda)$concont
summary(eigenvals(frax.rda, model = "constrained"))
screeplot(frax.rda) # visualize this information using a screeplot of the canonical eigenvalues

# Run a formal test of statistical significance of each constrained axis using anova 
# Note: takes quite long! (up to a few hours on large data sets)
anova.cca(frax.rda, by="axis") # or by term
signif.full <- anova.cca(frax.rda, parallel=getOption("mc.cores")) # default is permutation=999
signif.full
# vegan has a simple function for checking Variance Inflation Factors for the predictor variables used in the model. 

# Check collinearity among predictors
vif.cca(frax.rda) # very low values (<5) indicates that multicollinearity among these predictors shouldn’t be a problem

# Extract significant SNPs

# SNP loadings are stored as species in the RDA object
load.rda <- summary(frax.rda)$species # all axes
#load.rda <- summary(frax.rda)$species[,1] # onlyfot precipitation
head(load.rda)

# Look at histograms of the loadings on each RDA axis, see their (relatively normal) distribution.
par(mfrow=c(1,1))
hist(load.rda[,1], main="Loadings on RDA1")
abline(v= (mean(load.rda[,1]) + c(-1, 1) * 3 * sd(load.rda[,1])), 
       col = "#d00000", lty=3) # more stringent
abline(v= (mean(load.rda[,1]) + c(-1, 1) * 2.5 * sd(load.rda[,1])), 
       col = "#9d0208", lty=3) # less stringent treshhold

# Use if you chose model 2
hist(load.rda[,2], main="Loadings on RDA2")
abline(v= (mean(load.rda[,2]) + c(-1, 1) * 3 * sd(load.rda[,2])), 
       col = "#d00000", lty=3) # more stringent
abline(v= (mean(load.rda[,2]) + c(-1, 1) * 2.5 * sd(load.rda[,2])), 
       col = "#9d0208", lty=3) # less stringent treshhold

hist(load.rda[,3], main="Loadings on RDA3") 
abline(v= (mean(load.rda[,3]) + c(-1, 1) * 3 * sd(load.rda[,3])), 
       col = "#d00000", lty=3) # more stringent
abline(v= (mean(load.rda[,3]) + c(-1, 1) * 2.5 * sd(load.rda[,3])), 
       col = "#9d0208", lty=3) # less stringent treshhold

# 2.5 standard deviation cutoff (two-tailed p-value = 0.012).
# A bit more stringent: a 3 standard deviation cutoff (two-tailed p-value = 0.0027)

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
  x[x < lims[1] | x > lims[2]]           # locus names in these tails
}

cand_rd1 <- outliers(load.rda[,1], 3) 
# cand_rd2 <- outliers(load.rda[,2], 3) # use depending on the model chosen
# cand_rd3 <- outliers(load.rda[,3], 3)

# number of candidates
length(cand_rd1) 
names(cand_rd1)
# length(cand_rd2)  
# names(cand_rd2)
# length(cand_rd3) 
# names(cand_rd3)

# If several predictors were used
# how many candidates are shared by more predictors
frax.rda.cand <- c(names(outliers(load.rda[,1],3)),
                   names(outliers(load.rda[,2],3)),
                   names(outliers(load.rda[,3],3))
                  )
frax.rda.cand[duplicated(frax.rda.cand)] # no candidates are shared by more than one predictors
# the two SNPs also found in lfmm are associated with both T and Precipitation gradient

# Plot RDA results 

# Define colors for plotting
# Precipitation
pred.red$colPrec <- "#D86D46" 
pred.red$colPrec[pred.red$AnnPrec >= mean(pred.red$AnnPrec)] <- "#7BB2D9"

# Altitude
pred.red$colAlt <- "#E29578" #lightblue - low altitude
pred.red$colAlt[pred.red$Altitude >= mean(pred.red$Altitude)] <- "#95B2B0" #blue - high altitude

# T
pred.red$colT <- "#E29578" 
pred.red$colT[pred.red$AnnMeanT >= mean(pred.red$AnnMeanT)] <- "#95B2B0"
pred.red$shapeT <- 21 
pred.red$shapeT[pred.red$AnnMeanT >= mean(pred.red$AnnMeanT)] <- 22

# set colors for candidate SNPs
bgcol  <- ifelse(rownames(gen.imp) %in% unique(names(outliers(load.rda[,1], 3))), 'black', '#00000000')
snpcol <- ifelse(rownames(gen.imp) %in% unique(names(outliers(load.rda[,1], 3))), '#B9C211', '#00000000')
# bgcol_T  <- ifelse(rownames(gen.imp) %in% unique(names(outliers(load.rda[,2], 3))), 'black', '#00000000')
# snpcol_T <- ifelse(rownames(gen.imp) %in% unique(names(outliers(load.rda[,2], 3))), '#c1121f', '#00000000')
# bgcol_prec  <- ifelse(rownames(gen.imp) %in% unique(names(outliers(load.rda[,3], 3))), 'black', '#00000000')
# snpcol_prec <- ifelse(rownames(gen.imp) %in% unique(names(outliers(load.rda[,3], 3))), '#c1121f', '#00000000')

par(mfrow=c(1,1))
plot(frax.rda, 
     scaling=3, 
     choices=c(1,2), 
     type="n",
     ylim=c(-1.4,2),
     xlim=c(-1.5,1.5)
)
# the SNPs
points(frax.rda, #plot all the SNPs
       display="species", 
       pch=20, 
       cex=0.7, 
       col="gray", 
       scaling=3,
       choices=c(1,2))
points(frax.rda, #plot individual samples and separate them using colors for high or low precipitation (i.e. higher or lower than mean precipitation)
       display="sites", 
       pch=21, #shape by T 
       cex=1.5, 
       col="black",
       scaling=3,
       lwd = 1.8,
       choices=c(1,2),
       bg=pred.red$colPrec
           ) 
# candidate SNPs prec
points(frax.rda, 
       display="species", 
       pch=24, 
       cex=1.6, 
       col= bgcol, 
       bg=snpcol,
       lwd = 1.8,
       scaling=3,
       choices=c(1,2))
# the predictors
text(frax.rda,
     scaling=3,
     display="bp",
     col="black",
     #cex=1.8,
     choices=c(1,2),
     substitute(paste(bold("")))
)
legend("topright", 
       legend=c("Low precipitation","High precipitation","All SNPs","Candidate SNPs"), #"Altitude < 431m","Altitude > 431m","SNPs","candidate SNPs" #"early flushing","late flushing",
       bty="n", 
       col="black", 
       pch=c(21,21,20,24), 
       cex=1.5, pt.bg=c("#D86D46","#7BB2D9","gray","#B9C211")) #"#508991","#172A3A", #"#fca311","#1d3557","#a7c957","#386641",

# pdf("//bfwsbfile1/Institut2/Aglaia/EscheNot/Results/cRDA/flushing_included/RDA_plot.pdf")
# 
# dev.off()
