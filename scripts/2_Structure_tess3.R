# Title: "Population structure using TESS3"
# Author: Aglaia Szukala
# Short description: Script to estimate of spatial population structure using the TESS3 method

# Load libraries
library(tess3r)
library(maps)
library(raster)
library(rworldmap)
library(ggplot2)
library(ggrepel)
library(readxl)
library(directlabels)

# Set working directory
setwd("/path/to/working/directory")

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

# Change -1 to NA
gen[gen == -1] <- NA

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

# Extract geographic coordinates
coord <- as.matrix(cbind(phen$Y_Northing,phen$X_Easting))
rownames(coord) <- phen$Taxa
head(coord)

# Plot samples on map
plot(phen[,c(4,3)], pch = 19, cex = .5, 
     xlab = "Longitude (E)", ylab = "Latitude (N)")
map(add = T, interior = T) 
#text(phen[,c(4,3)], labels=phen$Taxa, cex=0.6, font=1)

# Four samples are missing geographic coordinates, these need to be removed to run Tess3
ids_to_remove <- phen[is.na(phen$Y_Northing),"Taxa"]
gen <- gen[,!(colnames(gen) %in% ids_to_remove)]
coord = coord[match(colnames(gen),rownames(coord)),]
dim(gen) # after removing NAs 1092 individuals left
dim(coord)
identical(colnames(gen), rownames(coord))


#######                                #######
####### 3 - GENETIC STRUCTURE in TeSS3 #######
#######                                #######

# Create a tess3 object
tess3.obj <- tess3(X = t(gen), 
                   coord = coord, #needs a coordinate matrix
                   K = 1:10, # K population clusters to test
                   method = "projected.ls", 
                   ploidy = 2, 
                   openMP.core.num = 4,
                   max.iteration = 500, 
                   rep = 50 
) 

# Visual cross-validation of K values, Fig S3A
plot(tess3.obj, 
     pch = 19, 
     col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score", main = "All samples"
)

# Retrieve tess3 Q matrix for several K values
for (i in 2:10) {
  assign(paste0("q.matrix", i),qmatrix(tess3.obj, K = i))
} 

# Produce a barplot of the Q-matrix 
my.colors <- c("#ffbc42", "#218380", "#d81159") # as many colors as K, add colors to plot higher K values
my.palette <- CreatePalette(my.colors)

barplot(q.matrix3, # choose q matrix
        border = NA, 
        space = 0, 
        col.palette = my.palette, 
        xlab = "Individuals", 
        ylab = "Ancestry proportions", 
        main = "Ancestry matrix",
        sort.by.Q = T)

# Order samples by increasing ancestry percentage
barplot.tess3Q = function(height, sort.by.Q = TRUE, col.palette = NULL, palette.length = 9, lab = FALSE, ...){
  Q = height
  rm(height)
  
  if (class(Q)[1] != "tess3Q") {warning("Object Q not of class tess3Q.")}
  ## defines a default color palette (8 colors)
  if (is.null(col.palette)) {
    cat("Use CreatePalette() to define color palettes.\n")
    if (ncol(Q) > 8)
      stop("The default color palette expects less than 9 clusters.")
    TestRequiredPkg("RColorBrewer")
    col.palette = list(
      c(RColorBrewer::brewer.pal(palette.length,"Reds")),
      c(RColorBrewer::brewer.pal(palette.length,"Greens")),
      c(RColorBrewer::brewer.pal(palette.length,"Blues")),
      c(RColorBrewer::brewer.pal(palette.length,"YlOrBr")),
      c(RColorBrewer::brewer.pal(palette.length,"RdPu")),
      c(RColorBrewer::brewer.pal(palette.length,"Greys")),
      c(RColorBrewer::brewer.pal(palette.length,"Purples")),
      c(RColorBrewer::brewer.pal(palette.length,"Oranges"))
    )
  }
  ## colors
  if (palette.length > 3) {
    colpal = sapply(col.palette, FUN = function(x) x[palette.length/2])}
  else {
    colpal = sapply(col.palette, FUN = function(x) x[palette.length])
  }
  
  if (sort.by.Q) {
    gr = apply(Q, MARGIN = 1, which.max)
    gm = max(gr)
    gr.o = order(sapply(1:gm, FUN = function(g) mean(Q[,g])))
    gr = sapply(gr, FUN = function(i) gr.o[i])
    or = order(gr)
    Qm = t(Q[or,])
    class(Qm) = "matrix"
    graphics::barplot(Qm, col = colpal,...)
    return(list(order = or))
  }
  else {
    Qm = t(Q)
    class(Qm) = "matrix"
    graphics::barplot(Qm, col = colpal, ...)
    return(list(order = 1:nrow(Q)))
  }
}

barplot.tess3Q(q.matrix2[order(q.matrix2[,1], decreasing = F),],
               space = 0,border = NA, 
               ylab = "Ancestry proportions",
               sort.by.Q = F,
               col.palette = my.palette)
barplot.tess3Q(q.matrix3[order(q.matrix3[,3], decreasing = F),],
               space = 0,border = NA, 
               ylab = "Ancestry proportions",
               sort.by.Q = F,
               col.palette = my.palette)
