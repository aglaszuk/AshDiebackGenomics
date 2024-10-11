# Title: "Principal Component Analysis of common ash genotype data"
# Author: Aglaia Szukala
# Short description: Script to perform PCA and visualize geographic distribution of mother trees

# Load libraries
library(SNPRelate) 
library(gdsfmt) 
library(maps) 
library(ggplot2)  
library(ggrepel)

# Set working directory
setwd("./data/")

#######                   #######
####### 1 - GENOTYPE DATA #######
#######                   #######

# Load genotype data (use 1 or 2)
# 1- Genotype data before filtering by maf, missingness, and monomorphic loci
gen <- read.table("2_Axiom_qualityFilterDat/Dat_11582SNPs_1107Ind.txt", # or test the data after filtering by maf, missingness and monomorphic loci
                  sep="\t", 
                  header=T, 
                  check.names = F)
gen <- as.matrix(gen)

# Remove duplicated individuals that were used as controls in plates
keep <- unique(colnames(gen))
gen <- gen[,select=keep]
dim(gen) # 11582  SNPs and 1096 individuals after removing duplicates

#2- Genotype data after filtering by maf, missingness, and monomorphic loci
gen <- read.table("3_ControlsRemoval_Maf_miss_mono_Filters/Dat_7313SNPs_1096Ind_mafMissMonoFiltered.txt",
                  sep="\t", 
                  header=T, 
                  check.names = F)
gen <- as.matrix(gen)

# Explore and edit the genotype data
gen[1:5,1:5] # SNPs as rows and samples as columns
dim(gen) # 11582 SNPs and 1107 Individuals (controls still included)

# Distribution of genotype classes
table(gen) # -1 are missing values
hist(gen)
gen[gen == -1] <- NA
hist(gen)


#######                               #######
####### 2 - PHENOTYPE / CLIMATIC DATA #######
#######                               #######

# Load and check phenotype data
phen <- read.table("5_PhenotypeClimDat/Dat_phenotype_climate.txt", 
                   sep="\t", 
                   header=T, 
                   check.names = F
)

str(phen)
head(phen)
dim(phen)

# Make sure individual IDs are in the same order in gen and phen
phen$Taxa <- as.character(phen$Taxa) # Change Taxa IDs to character
phen = phen[match(colnames(gen),phen$Taxa),]
identical(colnames(gen), phen$Taxa)

# Plot samples on map
plot(phen[,c(4,3)], pch = 19, cex = .5, 
     xlab = "Longitude (E)", ylab = "Latitude (N)")
map(add = T, interior = T) 
#text(phen[,c(4,3)], labels=phen$Taxa, cex=0.6, font=1)

# Four individuals miss geographic coordinates
phen[is.na(phen$Y_Northing),"Taxa"]

#######                      #######
####### 3 - PCA OF GENOTYPES #######
#######                      #######

# Use package SNPRelate to create a gds file 
snpgdsCreateGeno("gen.gds", genmat = gen,
                 sample.id = colnames(gen), snp.id = rownames(gen),
                 snpfirstdim=TRUE)

genofile <- snpgdsOpen("gen.gds")
#snpgdsClose(genofile)

# Perform PCA
pca <- snpgdsPCA(genofile,autosome.only=FALSE)

# Percent variation explained by PCs
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
barplot(round(pc.percent[1:32], 2), #Fig S3B
        ylab = "Proportion of variance explained by PCs")

# Create a data-frame for plotting including the first 4 PCs
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the third eigenvector
                  EV4 = pca$eigenvect[,4],    # the fourth eigenvector
                  stringsAsFactors = FALSE)

print(tab)
tab <- tab[!tab$sample.id %in% phen[is.na(phen$Y_Northing),"Taxa"],]

# Create a color variable for plotting
tab$y <- phen$Y_Northing[match(tab$sample.id, phen$Taxa)]
tab$x <- phen$X_Easting[match(tab$sample.id, phen$Taxa)]

# Plot PCA
ggplot(tab) +
  geom_point(aes(EV2,EV1, fill=EV2), pch=21, size=4) +
  scale_fill_gradient(low="#ffbc42", high="#d81159") + #"#ffbc42", "#218380", "#d81159"
  theme_classic(base_size = 15) +
  labs(x = paste("PC2 = ",round(pc.percent[1],1),"%",sep=""),
       y = paste("PC1 = ",round(pc.percent[2],1),"%",sep=""))

ggplot(tab) +
  geom_point(aes(EV2,EV1, fill=EV1), pch=21, size=4) +
  scale_fill_gradient(low="#218380", high="#d81159") + #"#ffbc42", "#218380", "#d81159"
  theme_classic(base_size = 15) +
  labs(x = paste("PC2 = ",round(pc.percent[1],1),"%",sep=""),
       y = paste("PC1 = ",round(pc.percent[2],1),"%",sep=""))

#Plot first 4 PC
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], labels=lbls)
