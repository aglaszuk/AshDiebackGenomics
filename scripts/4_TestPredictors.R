# Title: "Association amongs predictors"
# Author: Aglaia Szukala
# Short description: Script to check for associations among predictors used in the GEA models

# Load libraries
library(psych) # pairs.panels 
library(gplots) # heatmap

# Set working directory
setwd("./data/")

# Load phenotype data
phen <- read.table("5_PhenotypeClimDat/Dat_phenotype_climate.txt", 
                   sep="\t", 
                   header=T, 
                   check.names = F
)

head(phen)
dim(phen)
str(phen)

# Change Taxa IDs to character
phen$Taxa <- as.character(phen$Taxa) 
phen$Tolerance <- as.character(phen$Tolerance)
phen$Flush_1st <- as.character(phen$Flush_1st)
phen$Flush_2nd <- as.character(phen$Flush_2nd)

# Keep only data for which flushing is in the green and yellow category, the orange individuals died in the year of flushing assessment and should not be included in analyses on flushing behaviour
phen_sub1 <- phen[!phen$flushClass == "orange",]
dim(phen)
dim(phen_sub1)

# Keep only data for which flushing is in the green category, i.e. individuals bei which flushing cannot depend on ADB infection, as they were free of symptoms when flushing was assessed
phen_sub2 <- phen_sub1[!phen_sub1$flushClass == "yellow",]
dim(phen_sub2)

# Visualize correlation between categorical variables with association plot
par(mfrow=c(1,3))
assocplot(table(phen_sub1[,c(6,5)]),
          col = c("grey", "black"),
          xlab = "Damage Class (DC)",
          ylab = "R (germination year)")
assocplot(table(phen_sub1[,c(6,7)]),
          col = c("grey", "black"),
          xlab = "Damage Class (DC)", 
          ylab = "First Flush") 
assocplot(table(phen_sub1[,c(6,8)]),
          col = c("grey", "black"),
          xlab = "Damage Class (DC)", 
          ylab = "Second Flush")

# Only green individuals
assocplot(table(phen_sub2[,c(6,5)]),
          col = c("grey", "black"),
          xlab = "Damage Class (DC)", 
          ylab = "R (germination year)")
assocplot(table(phen_sub2[,c(6,7)]),
          col = c("grey", "black"),
          xlab = "Damage Class (DC)",
          ylab = "First Flush") 
assocplot(table(phen_sub2[,c(6,8)]),
          col = c("grey", "black"),
          xlab = "Damage Class (DC)", 
          ylab = "Second Flush")

# Chisquare test for correlation between categorical variables
chisq.test(table(phen[,c(5,6)])) 
chisq.test(table(phen_sub1[,c(6,7)])) 
chisq.test(table(phen_sub1[,c(6,8)]))
chisq.test(table(phen_sub2[,c(6,7)])) 
chisq.test(table(phen_sub2[,c(6,8)]))

# visualize correlations among continuous predictors
pairs.panels(phen[,c(9:14)], scale=F) 
