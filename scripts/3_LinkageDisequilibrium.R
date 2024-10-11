# Title: "Linkage decay by chromosome"
# Author: Aglaia Szukala following tutorial at https://eacooper400.github.io/gen8900/exercises/ld-1.html
# Short description: Obtain average LD, recombination rate under drift equilibrium [Hill and Weir (1988)], half decay distance and plot linkage decay

# Load libraries
library(ggplot2)

# Import TASSEL LD output file
ld <- read.delim(file.choose(),
                 stringsAsFactors = FALSE,
                 header=TRUE, 
                 sep = "\t"
                 ) 
# Inspect file
head(ld)
dim(ld)
str(ld)

# Remove NaN sites r2 calculations
ld_sub <- ld[ld$R.2 != "NaN",]
dim(ld_sub)

# Change distances between SNPs to numeric and remove NaN
ld_sub$dist <- as.numeric(ld_sub$Dist_bp)
ld_sub <- ld_sub[ld_sub$dist != "NaN",]
dim(ld_sub)

# Rename R.2
ld_sub$rsq <- ld_sub$R.2
head(ld_sub)

# Keep only informative columns
dat <- ld_sub[,c(1,2,7,8,15:19)]
head(dat)

# 1- Define 500bp intervals using the seq() function
bins=seq(1,max(dat$dist), by=500) #seq() function generates a series of numbers between a given minimum (1) and maximum value, in intervals of a specified size
head(bins)

# 2- Create a vector of 0s that will be replaced by the average r2 value in each bin.
my.means=rep(0, length(bins))

# Create a table with 2 columns: the first column will be the bin values, and the second will be the average r2 value for each bin.
LD.averages=data.frame(bins, my.means)
head(LD.averages)

# Loop through the list of intervals, find the subset of data that corresponds to each interval, and get the mean for that subset of data:

for (i in 1:length(bins)) {
  # Use the subset function to get all of the distance values that fall between the start of the ith bin and that value plus 500 (the interval size).
  data.interval=subset(dat, (dat$dist >= bins[i] & dat$dist < (bins[i]+500))) 
  LD.averages$my.means[i]=mean(data.interval$rsq) 
}

head(LD.averages)

# To estimate the population recombination rate (rho), fit a nonlinear model based on Hill and Weir’s (1988) equation
# N is the number of the genotypes that have the SNP site
hill.weir.eq=(rsq~(((10+(rho*dist))/((2+(rho*dist))*(11+(rho*dist))))*(1+(((3+(rho*dist))*(12+(12*(rho*dist))+((rho*dist)**2)))/(N*(2+(rho*dist))*(11+(rho*dist)))))))

# Rho values range from about 0.5 to 2, start with 0.1
rho.start=0.1

# R’s built in non-linear least squares function (nls). 
# Formula describes the relationship between the above equation, a starting value for the unknown parameter rho, and iteratively tries to find the value for the unknown parameter 
m=nls(formula=hill.weir.eq, data=dat, start=list(rho=rho.start)) 
results.m=summary(m)

results.m

# Test the fit of the model by looking at the correlation between the actual values and the predicted values:
cor(dat$rsq, predict(m))

# Extract the estimate of rho from the model results, and calculate the expected rho values using the same Hill and Weir equation, but now with your new estimate of rho plugged in:
rho.estimate=results.m$parameters[1]
Distance=sort(dat$dist)

# Obtain the expected distribution of the recombination parameter
exp.rsquared=(((10+(rho.estimate*Distance))/
                 ((2+(rho.estimate*Distance))*
                    (11+(rho.estimate*Distance))))*
                (1+(((3+(rho.estimate*Distance))*
                       (12+(12*(rho.estimate*Distance))+
                          ((rho.estimate*Distance)**2)))/
                      (n*(2+(rho.estimate*Distance))*
                         (11+(rho.estimate*Distance))))))


# Obtain max LD distance 
maxld <- max(exp.rsquared,na.rm=TRUE) #using max LD value from the expected distribution
halfdecay = maxld*0.5

# Obtain half LD decay distance
halfdecaydist <- Distance[which.min(abs(exp.rsquared-halfdecay))]
halfdecaydist

# Reorder data by distance
#newdat <- newdat[order(newdat$dat.dist),]

# Calculate recombination rate (frequency)
dat$rec <- dat$dist*rho.estimate
dat$rec
head(dat)

# Plot the fitted curve by adding a line of the expected results to the graph:
plot(x=dat$dist, 
     y=dat$rsq, 
     xlab="Distance between SNPs (bp)", 
     ylab=expression(R^2), 
     pch=20, 
     col=rgb(0,0,0,alpha=0.2), 
     main="Decay of Linkage Disequilibrium",
     xlim=c(0,4e+06)
)
# Add points to the plot to show the interval average values (I will use 2 commands: one to add the points, and then a second to connect the points with a line):
points(x=LD.averages$bins, y=LD.averages$my.means, col="red", pch=20) 
lines(x=LD.averages$bins, y=LD.averages$my.means, col="red", lwd=2)
# Add a Smoothed Line to Show the Trend of Decay
lines(Distance, exp.rsquared, col="purple", lwd=2)
legend(100000,0.99, c("Means", "Expected R2"), lty=c(1,1), col=c("red", "purple"))

# Make plot
ggplot(dat,
       aes(x=dist, y=rsq, color=rec))+ # distance in centimorgan (cM)
  geom_point(pch=20) + 
  scale_colour_gradient(low="#218380", high="#ffbc42") +
  theme_classic() +
  xlab("Distance between SNPs (bp)") +
  ylab(expression(LD ~ (R^2))) +
  labs(color = "Recombination \n frequency") +
  geom_vline(xintercept=halfdecaydist, linetype="dashed", color = "darkred") +
  geom_text(aes(x=6, label=paste0(round(halfdecaydist/1000, digits=0), " Kb"), y=0.5),colour="black")+
  geom_line(color='red',aes(x=Distance, y=exp.rsquared)) # Add a Smoothed Line to Show the Trend of Decay
#+ xlim(0,1.5e+06) #visualize only the lower part of the x axis to zoom into the "interesting" part of the plot


