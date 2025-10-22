# Title: "Plot linked genome around SNPs"
# Author: Aglaia Szukala
# Short description: Script plots the area downstream and upstream of SNPs that is in LD based on LD decay estimates 

library(chromPlot)


dat <- read.table(file="6_Data4LDscripts/intervals.txt")
colnames(dat) <- c("Chrom","Start","End","Name")

# Define areas in linkage with SNPs or not
dat$gieStain <- "gneg"
dat$gieStain[!dat$Name == "uncovered"] <- "gpos"

dat$Colors <- "gray80"
dat$Colors[!dat$Name == "uncovered"] <- "black"

# Plot data
par(mai = c(0.0, 0.2, 0.2, 0.2),
    omi = c(0, 0, 0, 0))  
chromPlot(gaps=dat, 
          bands= dat, 
          figCols=5, 
          segLwd=4,
          cex = 1.2
          )