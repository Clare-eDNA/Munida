library(lubridate)
library(readxl)
ALLData<- read_xlsx("Munida_continuous values_only.xlsx", sheet = "ODV format")

Munida<-ALLData[,c(2,4,5,6,7,8,10,14,16,28)]

Munida$Station<-as.factor(Munida$Station)
levels(Munida$Station)
Munida$`mon/day/yr` <- as.Date(Munida$`mon/day/yr`, format="%Y-%m-%d")

Munida <- as.data.frame(Munida)

unique(Munida$`mon/day/yr`)

Munida <- rbind(Munida[Munida$`mon/day/yr` == "2017-02-02",],
      Munida[Munida$`mon/day/yr` == "2017-02-23",])

Munida <- na.omit(Munida)
Munida <- Munida[Munida$Station != 0,]

names(Munida)

#### Ploting temperature & salinity on same graph:
# Salinity:
library(lattice)
xycolSal=list(superpose.symbol = list(col=colorRampPalette(c("lightblue", "blue"))(10),
                                      fill=Munida$`mon/day/yr`, pch=16),
              superpose.line=list(col=colorRampPalette(c("lightblue", "blue"))(10),
                                  fill=Munida$`mon/day/yr`),
  strip.background=list(col="grey90"),
              strip.border=list(col="black"),
              simpleTheme(col=2, lty=1, col.line=c("darkgoldenrod1","blue")))
XYplotSal=xyplot(`Salinity [psu]` ~ `distance [km from TH]` | `mon/day/yr`,
                 data=Munida,
                 group=Munida$`mon/day/yr`,
                 type=c("p","l","g"),
                 xlab=list("Distance (km)", cex=3),
                 ylab=list("Salinity (PSU)", cex=3),
                 ylab.right=list("Temperature (ºC)", cex=3),
                 cex.lab=3,
                 layout=c(1,2),
                 as.table=T,
                 scales=list(alternating=T, y=list(relation="free"), cex=2),
                 par.settings=xycolSal, left.padding=2)
XYplotSal
# Temperature:
xycolTemp=list(superpose.symbol = list(col=colorRampPalette(c("darkgoldenrod1", "brown"))(6),
                                       fill=Munida$`mon/day/yr`, pch=16),
               superpose.line=list(col=colorRampPalette(c("darkgoldenrod1", "brown"))(6),
                                   fill=Munida$`mon/day/yr`),
  strip.background=list(col="grey90"),
               strip.border=list(col="black"))
XYplotTemp=xyplot(`Temperature [oC]` ~ `distance [km from TH]` | `mon/day/yr`,
                  data=Munida,
                  group=Munida$`mon/day/yr`,
                  type=c("p","l","g"),
                  ylab=list("Temperature (ºC)", cex=3),
                  layout=c(1,2),
                  as.table=T,
                  scales=list(alternating=T, y=list(relation="free"), cex=2),
                  par.settings=xycolTemp)
XYplotTemp
# Add the legend afterward...
library(latticeExtra)
doubleYScale(XYplotSal,XYplotTemp,
             use.style=F,
             add.ylab2=F)


### Ecotone:

library(EcotoneFinder)
library(ggplot2)

EcoTempSalCruise1 <- EcotoneFinder(data = Munida[Munida$`mon/day/yr` == "2017-02-02", c(8,9)],
                                  dist = Munida[Munida$`mon/day/yr` == "2017-02-02",]$`distance [km from TH]`,
                                  method = c("cmeans", "fanny"), groups = 4, m.exp = 2,
                                  na.rm =T, standardize = "standardize")

ggEcotone(EcoTempSalCruise1, method = "fanny", col = c("red", "gold", "green", "blue")) +
  theme_bw()+labs(title="2 Feb 2017",x ="KM from Taiaroa head")+ theme(plot.title = element_text(size=30))

EcoTempSalCruise2 <- EcotoneFinder(data = Munida[Munida$`mon/day/yr` == "2017-02-23", c(8,9)],
                                   dist = Munida[Munida$`mon/day/yr` == "2017-02-23",]$`distance [km from TH]`,
                                   method = c("all"), groups = 4, m.exp = 2,
                                   na.rm =T, standardize = "standardize")

ggEcotone(EcoTempSalCruise2, method = "fanny", col = c("red", "gold", "green","blue")) +
  theme_bw()+labs(title="23 Feb 2017", x ="KM from Taiaroa head")+ theme(plot.title = element_text(size=30))


