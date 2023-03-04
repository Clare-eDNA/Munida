# Playing around with Temperature and such

setwd("~/Documents/MunidaTransect/data/TempSal")

install.packages("~/Documents/MunidaTransect/data/TempSal/EcotoneFinder_0.2.0.tar.gz", dependencies=TRUE,repos = NULL, type = "source")
library(EcotoneFinder)
require(ggplot2)
require(colorspace)
require(Rmisc)

metadat<-read.csv("Metadata.csv")

metadat$X<-NULL
metadat$NO3.N<-NULL
metadat$Zone<-NULL

head(metadat)

str(metadat)

data<-na.omit(metadat)
data$lat<-abs(data$lat)

head(data)

## Looking at temperature gradients
?EcotoneFinder
EcoFinder <- EcotoneFinder(data[ , c("Temp","Sal")], data$MetersFrmShore, method = c("cmeans"), groups = 3, m.exp = 2, standardize = "standardize")
EcoFinder

?ggEcotone
Plot1 <- ggEcotone(EcoFinder, slope = NULL, plot.data = TRUE,
                      method = c("cmeans"),
                      col = c("#D33F6A", "#E99A2C", "#E2E6BD"),
                      facet = list(c("data"), c("cmeans")),
                      title = "Species distribution and fuzzy clusters",
                      xlab = "Gradient", ylab = "Membership grades") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme_bw()
Plot1


Plot <- ggEcotone(EcoFinder, slope = NULL, plot.data = FALSE,
                  method = c("cmeans"),
                  col = c("#D33F6A", "#E99A2C", "#E2E6BD"),
                  facet = c("cmeans"),
                  title = "fuzzy clusters and derivatives",
                  xlab = "Gradient", ylab = "Membership grades") +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme_bw()
Plot

## Clustering on Species data

EcoFinder <- EcotoneFinder(data[ , 10:ncol(data)], data$MetersFrmShore, method = "all", groups = 2, m.exp = 2, diversity= "all", na.rm =T,standardize = "hellinger")






#playing around with other stuff

plot<-plot(data$Temp~data$Location)
plot2<-plot(data$Sal~data$Location)

matplot(data$Location, cbind (data$Temp, data$Sal),pch = 19)

