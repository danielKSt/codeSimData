# We perform some initial testing on pixels, not doing the point process stuff here

# Data can be imported here ----
library("rhdf5")
library("tidyverse")
setwd("/Users/danielks/Library/CloudStorage/OneDrive-NTNU/PhD/Simulated data")
completeSurfaceVortexData <- h5read(file = 'surface_vortices.h5', name = 'surf')
completeSurfaceDivergenceData <- h5read(file = 'surface_vortices.h5', name = 'sdiv')
endTimestep <- dim(completeSurfaceVortexData)[1]

#### Check for "spatial uniformness" of vortices with various length of timesteps ----

find_rate_outliers <- function(sizeOfTimestep, vortexData){
  print(sizeOfTimestep)
  endTimestep <- dim(vortexData)[1]
  numberOfTimesteps <- length(seq(from = 1, to = endTimestep, by = sizeOfTimestep))
  numberOfTimesWithVortexOnPixel <- matrix(integer(256*256), nrow = 256)
  for (i in seq(from = 1, to = endTimestep, by = sizeOfTimestep)) {
    numberOfTimesWithVortexOnPixel <- numberOfTimesWithVortexOnPixel + vortexData[i, 1:256, 1:256]
  }
  rateOfTimesWithVortexOnPixel <- numberOfTimesWithVortexOnPixel/numberOfTimesteps
  
  meanNoVortexOnPixel = mean(numberOfTimesWithVortexOnPixel)
  meanRateOfVortex <- meanNoVortexOnPixel/numberOfTimesteps
  approxStDev <- sqrt(numberOfTimesteps*meanRateOfVortex*(1-meanRateOfVortex))
  
  data <- expand.grid(X = 1:256, Y = 1:256)
  data$vorticesOnPixel <- array(numberOfTimesWithVortexOnPixel)
  pltres <- ggplot(data, mapping = aes(x = X, y = Y, fill = vorticesOnPixel))+
    geom_tile() + 
    ggtitle(sprintf("Distribution of vortices with step size %i", sizeOfTimestep)) +
    theme(legend.position = "bottom")
  
  bigOutlierPixels <- (numberOfTimesWithVortexOnPixel > meanNoVortexOnPixel + 2*approxStDev)
  smallOutlierPixels <- (numberOfTimesWithVortexOnPixel < meanNoVortexOnPixel - 2*approxStDev)
  
  return(list("rateOfOutliers" = (sum(bigOutlierPixels)+sum(smallOutlierPixels))/256^2, "PlotResult" = pltres, "NumberPerPixel" = numberOfTimesWithVortexOnPixel))
}

rateOutliersData <- vector(mode = "list", length = 5)
stepSizes <- c(1, 5, 10, 20, 40)
rateOutliersData <- stepSizes |> lapply(find_rate_outliers, vortexData = completeSurfaceVortexData)

png("Presentation_december23/vortexDistribution.png", res = 100, width = 2500, height = 1545)
rateOutliersData[[1]]$PlotResult
dev.off()

library("gridExtra")
png("Presentation_december23/vortexDistributionmanyStepSizes.png", res = 400, width = 10000, height = 6179)
grid.arrange(rateOutliersData[[1]]$PlotResult, rateOutliersData[[2]]$PlotResult, rateOutliersData[[3]]$PlotResult, 
             rateOutliersData[[4]]$PlotResult, ncol = 2, nrow = 2)
dev.off()

for(i in c(1:5)){
  a <- sprintf("Step size %i: Outlier rate = %.3f", stepSizes[i], rateOutliersData[[i]]$rateOfOutliers)
  print(a)
}

#### Check various correlations for pixels ----

nullToNull <- matrix(integer(256*256), nrow = 256)
nullToOne <- matrix(integer(256*256), nrow = 256)
oneToNull <- matrix(integer(256*256), nrow = 256)
oneToOne <- matrix(integer(256*256), nrow = 256)
for (t in c(2:endTimestep)) {
  nullToNull <- nullToNull + ((1-completeSurfaceVortexData[t, 1:256, 1:256])*(1-completeSurfaceVortexData[t-1, 1:256, 1:256]))
  nullToOne <- nullToOne + (completeSurfaceVortexData[t, 1:256, 1:256]*(1-completeSurfaceVortexData[t-1, 1:256, 1:256]))
  oneToNull <- oneToNull + ((1-completeSurfaceVortexData[t, 1:256, 1:256])*completeSurfaceVortexData[t-1, 1:256, 1:256])
  oneToOne <- oneToOne + (completeSurfaceVortexData[t, 1:256, 1:256]*completeSurfaceVortexData[t-1, 1:256, 1:256])
}
p1 <- oneToOne/(oneToOne+oneToNull)
#p1NansimpleReplace <- p1
#p1NansimpleReplace[which(is.nan(p1))] <- mean(na.omit(p1))
#p1NansimpleReplace[which(p1NansimpleReplace==0)] <- mean(na.omit(p1))
p0 <- nullToOne/(nullToOne+nullToNull)

p1mean <- mean(na.omit(p1))
p0mean <- mean(na.omit(p0))

p1variances <- sqrt((oneToNull + oneToOne)*p1mean*(1-p1mean))
p0variances <- sqrt((nullToNull + nullToOne)*p0mean*(1-p0mean))

p1bigOutliers <- (oneToOne > (oneToOne+oneToNull)*p1mean + 2*p1variances)
p1smallOutliers <- (oneToOne < (oneToOne+oneToNull)*p1mean - 2*p1variances)
p1outliers <- (oneToOne > (oneToOne+oneToNull)*p1mean + 2*p1variances) + (oneToOne < (oneToOne+oneToNull)*p1mean - 2*p1variances)
sum(p1outliers)/256^2

p0bigOutliers <- (nullToOne > (nullToOne+nullToNull)*p0mean + 2*p0variances)
p0smallOutliers <- (nullToOne < (nullToOne+nullToNull)*p0mean - 2*p0variances)
p0outliers <- (nullToOne > (nullToOne+nullToNull)*p0mean + 2*p0variances) + (nullToOne < (nullToOne+nullToNull)*p0mean - 2*p0variances)
sum(p0outliers)/256^2


#We make some plots of the calculated probabilities and the outliers we have found

data <- expand.grid(X = 1:256, Y = 1:256)
data$probFromOne <- array(p1)
data$bigOutliersFromOne <- array(p1bigOutliers)
data$smallOutliersFromOne <- array(p1smallOutliers)
data$outliersFromOne <- array(p1outliers)
data$probFromZero <- array(p0)
data$bigOutliersFromZero <- array(p0bigOutliers)
data$smallOutliersFromZero <- array(p0smallOutliers)
data$outliersFromZero <- array(p0outliers)


pf1 <- ggplot(data, mapping = aes(x = X, y = Y, fill = probFromOne))+
  geom_tile() + 
  theme(legend.position = "bottom")

pf0 <- ggplot(data, mapping = aes(x = X, y = Y, fill = probFromZero))+
  geom_tile() + 
  theme(legend.position = "bottom")

of1 <- ggplot(data, mapping = aes(x = X, y = Y, fill = outliersFromOne))+
  geom_tile() + 
  theme(legend.position = "bottom")

bof1 <- ggplot(data, mapping = aes(x = X, y = Y, fill = bigOutliersFromOne))+
  geom_tile() + 
  theme(legend.position = "bottom")

sof1 <- ggplot(data, mapping = aes(x = X, y = Y, fill = smallOutliersFromOne))+
  geom_tile() + 
  theme(legend.position = "bottom")

of0 <- ggplot(data, mapping = aes(x = X, y = Y, fill = outliersFromZero))+
  geom_tile() + 
  theme(legend.position = "bottom")

bof0 <- ggplot(data, mapping = aes(x = X, y = Y, fill = bigOutliersFromZero))+
  geom_tile() + 
  theme(legend.position = "bottom")

sof0 <- ggplot(data, mapping = aes(x = X, y = Y, fill = smallOutliersFromZero))+
  geom_tile() + 
  theme(legend.position = "bottom")

png(filename = "Presentation_december23/probabilitiesMediumSimpleModel.png", res = 400, width = 2500, height = 1545)
grid.arrange(pf0, pf1, ncol = 2)
dev.off()

png(filename = "Presentation_december23/probabilitiesMediumSimpleModel_2.png", res = 400, width = 2500, height = 3090)
grid.arrange(bof0, bof1, sof0, sof1, ncol = 2)
dev.off()

#### We compare the raw pixels with the beta squared values ----

numberOfVortices <- c(1:endTimestep)
betaSquared <- c(1:endTimestep)

for (i in c(1:endTimestep)) {
  numberOfVortices[i] <- sum(completeSurfaceVortexData[i, 1:256, 1:256])
  betaSquared[i] <- sum(completeSurfaceDivergenceData[i, 1:256, 1:256]^2)
}
numberOfVortices <- (100/max(numberOfVortices))*numberOfVortices
betaSquared <- (100/max(betaSquared))*betaSquared

numberOfVorticesAndBetaSquared <- data.frame(t = c(1:endTimestep), pixelsWithVortex = numberOfVortices, betaSq = betaSquared)
numberOfVorticesAndBetaSquared <- melt(numberOfVorticesAndBetaSquared, id.vars = "t")
ggplot(data = numberOfVorticesAndBetaSquared, aes(x = t, y = value, color = variable))+
  geom_line()



