####################################################
####   Script for analysing derived variables   ####
####################################################

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","dplyr","ggplot2","RColorBrewer","lme4","gridExtra","effects","caper","car", "blighty", "plotrix","reshape2","effects","MuMIn")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

source("CheckResidsFunction.R")

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### read in butterfly and moth datasets ####

## we have three recording levels of data to work with, depending on recording intensity in the distribution datasets
# let's start with the most basic, 'recorded'

# read in butterfly data
butterfly.raw <- read.csv("Data/Derived/BMS_BNM.csv", header = T)
summary(butterfly.raw)

# and moth data
moth.raw <- read.csv("Data/Derived/RIS_NMRS.csv", header = T)
summary(moth.raw)

## merge together, labelling butterflies and moths
butterfly.raw$TAXON <- "Butterfly"
moth.raw$TAXON <- "Moth"

raw.data <- rbind(butterfly.raw,moth.raw)

### now we need to mark out which species are 'northerly-distributed' and can't therefore be used in analyses of Northern Range Margins
## this data is held in the 'species to use' tables
## read in 'species to use' tables
butterfly.to.use <- read.csv("Data/butterflyspeciestouse.csv", header = T)
moth.to.use <- read.csv("Data/mothspeciestouse.csv", header = T)

to.use <- rbind(butterfly.to.use[,c(1,8)], moth.to.use[,c(1,8)])


## merge in details of northerly species
data.dist <- merge(to.use,raw.data)

### we also want details of species' voltinism - 
# at this point just fairly crudely as a categorical variable (i.e. "univoltine"/"multivoltine")

voltinism <- read.csv("Data/voltinism.csv", header = T)

data.volt <- merge(voltinism[,-1], data.dist)

data.volt$VOLTINISM <- as.factor(data.volt$VOLTINISM)
data.volt$STRICT_VOLTINISM <- as.factor(data.volt$STRICT_VOLTINISM)
data.volt$TAXON <- as.factor(data.volt$TAXON)
summary(data.volt)

## finally, for each of the four derived variables, let's generate a 'state' categorical variable for what the trend is
data.volt$PHENO.STATE <- as.factor(ifelse(data.volt$PHENO.P > 0.05, "Non-significant",
                                          ifelse(data.volt$PHENO.SLOPE < 0, "Significantly later emergence",
                                                 "Significantly earlier emergence")))

data.volt$ABUND.STATE <- as.factor(ifelse(data.volt$ABUND.P > 0.05, "Non-significant",
                                          ifelse(data.volt$ABUND.SLOPE < 0, "Significantly declining abundance",
                                                 "Significantly increasing abundance")))

data.volt$DISTRIB.STATE <- as.factor(ifelse(data.volt$DISTRIB.P > 0.05, "Non-significant",
                                          ifelse(data.volt$DISTRIB.SLOPE < 0, "Significantly retracting distribution",
                                                 "Significantly expanding distribution")))

data.volt$MARGIN.STATE <- as.factor(ifelse(data.volt$MARGIN.P > 0.05, "Non-significant",
                                            ifelse(data.volt$MARGIN.SLOPE < 0, "Significantly retracting southwards",
                                                   "Significantly expanding northwards")))


## divide it back into butterflies and moths
data.b <- data.volt[which(data.volt$TAXON=="Butterfly"), ]
data.m <- data.volt[which(data.volt$TAXON=="Moth"), ]


## and also divide it into univoltine and multivoltine
data.uv <- data.volt[which(data.volt$VOLTINISM==1), ]
data.mv <- data.volt[which(data.volt$VOLTINISM==2), ]

## and finally into habitat specialist and generalist
data.hs <- data.volt[which(data.volt$CLASS=="HS"), ]
data.wc <- data.volt[which(data.volt$CLASS=="WC"), ]

### analysis ####
## first let's plot up some histograms of each variable, coloured by significance

hist1 <- ggplot(data.volt)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("All Lepidoptera")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1

# butterflies
hist1.b <- ggplot(data.b)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Butterflies")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.b


# moths
hist1.m <- ggplot(data.m)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Moths")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.m


## abundance
hist2 <- ggplot(data.volt)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2


# butterflies
hist2.b <- ggplot(data.b)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.b


# moths
hist2.m <- ggplot(data.m)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.m


## distribution
hist3 <- ggplot(data.volt)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3



# butterflies
hist3.b <- ggplot(data.b)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.b


# moths
hist3.m <- ggplot(data.m)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.m



## northern range margin
hist4 <- ggplot(data.volt[which(data.volt$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4

# butterflies
hist4.b <- ggplot(data.b[which(data.b$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.b


# moths
hist4.m <- ggplot(data.m[which(data.m$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.m


m.hist <- grid.arrange(hist1.b,hist1,hist1.m,
                       hist2.b,hist2,hist2.m,
                       hist3.b,hist3,hist3.m,
                       hist4.b,hist4,hist4.m,
                       ncol=3)


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

leg.hist1 <- g_legend(hist1.b)
leg.hist2 <- g_legend(hist2.b)
leg.hist3 <- g_legend(hist3.b)
leg.hist4 <- g_legend(hist4.b)

m.hist.a <- grid.arrange(arrangeGrob(hist1.b + theme(legend.position="none"),hist1 + theme(legend.position="none"),hist1.m + theme(legend.position="none"),
                                     hist2.b + theme(legend.position="none"),hist2 + theme(legend.position="none"),hist2.m + theme(legend.position="none"),
                                     hist3.b + theme(legend.position="none"),hist3 + theme(legend.position="none"),hist3.m + theme(legend.position="none"),
                                     hist4.b + theme(legend.position="none"),hist4 + theme(legend.position="none"),hist4.m + theme(legend.position="none"),
                                     ncol=3),
                         arrangeGrob(leg.hist1,leg.hist2,
                                     leg.hist3,leg.hist4,ncol =1),ncol=2, widths=c(10,2))

ggsave("Results/Plots/Histograms.svg", plot = m.hist.a, device = svg, width = 350, height = 350, units = "mm", limitsize = T)

## now also plot histograms of SSI scores, to see whether this might drive the differences between butterflies and moths in the distribution data

## SSI
hist5 <- ggplot(data.volt)+
  geom_histogram(aes(x = SSI, fill = CLASS), bins=15)+
  theme_bw()+
  xlim(-0,2.5)+
  ylim(0,100)+
  xlab("Species Specialisation Index")+
  ylab("Number of species")+
  ggtitle("All Lepidoptera")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Specialisation class")

hist5

# butterflies
hist5.b <- ggplot(data.volt[which(data.volt$TAXON=="Butterfly"), ])+
  geom_histogram(aes(x = SSI, fill = CLASS), bins=15)+
  theme_bw()+
  xlim(0,2.5)+
  ylim(0,100)+
  xlab("Species Specialisation Index")+
  ylab("Number of species")+
  ggtitle("Butterflies")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Specialisation class")

hist5.b


# moths
hist5.m <- ggplot(data.volt[which(data.volt$TAXON=="Moth"), ])+
  geom_histogram(aes(x = SSI, fill = CLASS), bins=15)+
  theme_bw()+
  xlim(0,2.5)+
  ylim(0,100)+
  xlab("Species Specialisation Index")+
  ylab("Number of species")+
  ggtitle("Moths")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Specialisation class")

hist5.m

# arrange and export

leg.hist5 <- g_legend(hist5.b)

m.hist5.a <- grid.arrange(arrangeGrob(hist5.b + theme(legend.position="none"),hist5 + theme(legend.position="none"),hist5.m + theme(legend.position="none"),
                                     ncol=3),
                         leg.hist5,ncol=2, widths=c(10,2))

ggsave("Results/Plots/Final/FigA.svg", plot = m.hist5.a, device = svg, width = 350, height = 100, units = "mm", limitsize = T)


## and now some histograms making the same comparison but between univoltine and multivoltine species

## phenology
# uni
hist1.uv <- ggplot(data.uv)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Univoltine species")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.uv


# multi
hist1.mv <- ggplot(data.mv)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Multivoltine species")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.mv


## abundance
# univoltine
hist2.uv <- ggplot(data.uv)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.uv


# multivoltine
hist2.mv <- ggplot(data.mv)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.mv


## distribution
# univoltine
hist3.uv <- ggplot(data.uv)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.uv


# multivoltine
hist3.mv <- ggplot(data.mv)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.mv



## northern range margin
# univoltine
hist4.uv <- ggplot(data.uv[which(data.uv$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.uv


# moths
hist4.mv <- ggplot(data.mv[which(data.mv$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.mv


m.hist.v <- grid.arrange(hist1.uv,hist1,hist1.mv,
                       hist2.uv,hist2,hist2.mv,
                       hist3.uv,hist3,hist3.mv,
                       hist4.uv,hist4,hist4.mv,
                       ncol=3)


leg.hist1v <- g_legend(hist1.uv)
leg.hist2v <- g_legend(hist2.uv)
leg.hist3v <- g_legend(hist3.uv)
leg.hist4v <- g_legend(hist4.uv)

m.hist.va <- grid.arrange(arrangeGrob(hist1.uv + theme(legend.position="none"),hist1 + theme(legend.position="none"),hist1.mv + theme(legend.position="none"),
                                     hist2.uv + theme(legend.position="none"),hist2 + theme(legend.position="none"),hist2.mv + theme(legend.position="none"),
                                     hist3.uv + theme(legend.position="none"),hist3 + theme(legend.position="none"),hist3.mv + theme(legend.position="none"),
                                     hist4.uv + theme(legend.position="none"),hist4 + theme(legend.position="none"),hist4.mv + theme(legend.position="none"),
                                     ncol=3),
                         arrangeGrob(leg.hist1v,leg.hist2v,
                                     leg.hist3v,leg.hist4v,ncol =1),ncol=2, widths=c(10,2))

ggsave("Results/Plots/VoltinismHistograms.svg", plot = m.hist.va, device = svg, width = 350, height = 350, units = "mm", limitsize = T)


## and finally some histograms making the same comparison but between specialist and generalist species

## phenology
# hs
hist1.hs <- ggplot(data.hs)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Habitat specialists")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.hs


# wc
hist1.wc <- ggplot(data.wc)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Wider countryside generalists")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.wc


## abundance
# hs
hist2.hs <- ggplot(data.hs)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.hs


# wc
hist2.wc <- ggplot(data.wc)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.wc


## distribution
# hs
hist3.hs <- ggplot(data.hs)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.hs


# wc
hist3.wc <- ggplot(data.wc)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.wc



## northern range margin
# hs
hist4.hs <- ggplot(data.hs[which(data.hs$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.hs


# wc
hist4.wc <- ggplot(data.wc[which(data.wc$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.wc


m.hist.s <- grid.arrange(hist1.hs,hist1,hist1.wc,
                         hist2.hs,hist2,hist2.wc,
                         hist3.hs,hist3,hist3.wc,
                         hist4.hs,hist4,hist4.wc,
                         ncol=3)


leg.hist1s <- g_legend(hist1.hs)
leg.hist2s <- g_legend(hist2.hs)
leg.hist3s <- g_legend(hist3.hs)
leg.hist4s <- g_legend(hist4.hs)

m.hist.sa <- grid.arrange(arrangeGrob(hist1.hs + theme(legend.position="none"),hist1 + theme(legend.position="none"),hist1.wc + theme(legend.position="none"),
                                      hist2.hs + theme(legend.position="none"),hist2 + theme(legend.position="none"),hist2.wc + theme(legend.position="none"),
                                      hist3.hs + theme(legend.position="none"),hist3 + theme(legend.position="none"),hist3.wc + theme(legend.position="none"),
                                      hist4.hs + theme(legend.position="none"),hist4 + theme(legend.position="none"),hist4.wc + theme(legend.position="none"),
                                      ncol=3),
                          arrangeGrob(leg.hist1s,leg.hist2s,
                                      leg.hist3s,leg.hist4s,ncol =1),ncol=2, widths=c(10,2))

ggsave("Results/Plots/SpecialisationHistograms.svg", plot = m.hist.sa, device = svg, width = 350, height = 350, units = "mm", limitsize = T)



### now we are ready to start analysing relationships between phenology and other variables

## first, abundance
# start with a basic plot
g1 <- ggplot(data.volt)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g1


# now test the relationship using a linear model

hist(data.volt$ABUND.SLOPE)

model1 <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
             data = data.volt)

summary(model1)
drop1(model1, test = "Chi")
Anova(model1, type = "III")


## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine
hist(data.uv$ABUND.SLOPE)

model1.u <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
             data = data.uv)

summary(model1.u)
drop1(model1.u, test = "Chi")


# predict data to plot lines
newdata1u <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1u$ABUND.SLOPE <- predict(model1.u, newdata = newdata1u, type="response")
newdata1u <- ddply(newdata1u, .(PHENO.SLOPE), numcolwise(mean))

g1.u <- ggplot(data.uv)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata1u,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-0.15,0.1)+
  xlim(-3,2)+
  ggtitle("Univoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g1.u


# multivoltine
hist(data.mv$ABUND.SLOPE)

model1.m <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv)

summary(model1.m)
drop1(model1.m, test = "Chi")


# predict data to plot lines
newdata1m <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1m$ABUND.SLOPE <- predict(model1.m, newdata = newdata1m, type="response")
newdata1m <- ddply(newdata1m, .(PHENO.SLOPE), numcolwise(mean))

g1.m <- ggplot(data.mv)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata1m,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-0.15,0.1)+
  xlim(-3,2)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g1.m



## combine and export the plots
m1 <- grid.arrange(g1.u,g1.m,
                   ncol=2)


leg1 <- g_legend(g1.u)

m1a <- grid.arrange(arrangeGrob(g1.u + theme(legend.position="none"),g1.m + theme(legend.position="none"),
                                     ncol=2),
                         leg1,nrow=2, heights=c(10,2))

ggsave("Results/Plots/Abundance.png", plot = m1a, width = 350, height = 200, units = "mm", limitsize = F)






## second, distribution
# start with a basic plot
g2 <- ggplot(data.volt)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g2

# now test the relationship using a linear model

model2 <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
               data = data.volt)

summary(model2)
drop1(model2, test = "Chi")



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine
hist(data.uv$DISTRIB.SLOPE)

model2.u <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.uv)

summary(model2.u)
drop1(model2.u, test = "Chi")


# predict data to plot lines
newdata2u <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2u$DISTRIB.SLOPE <- predict(model2.u, newdata = newdata2u, type="response")
newdata2u <- ddply(newdata2u, .(PHENO.SLOPE), numcolwise(mean))


g2.u <- ggplot(data.uv)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata2u,
            aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-3,2)+
  ggtitle("Univoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g2.u


# multivoltine
hist(data.mv$DISTRIB.SLOPE)

model2.m <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv)

summary(model2.m)
drop1(model2.m, test = "Chi")


# predict data to plot lines
newdata2m <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2m$DISTRIB.SLOPE <- predict(model2.m, newdata = newdata2m, type="response")
newdata2m <- ddply(newdata2m, .(PHENO.SLOPE), numcolwise(mean))

g2.m <- ggplot(data.mv)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata2m,
            aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-3,2)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g2.m


## combine and export the plots
m2 <- grid.arrange(g2.u,g2.m,
                   ncol=2)


leg2 <- g_legend(g2.u)

m2a <- grid.arrange(arrangeGrob(g2.u + theme(legend.position="none"),g2.m + theme(legend.position="none"),
                                ncol=2),
                    leg2,nrow=2, heights=c(10,2))

ggsave("Results/Plots/Distribution.png", plot = m2a, width = 350, height = 200, units = "mm", limitsize = F)




## third, NRM
# start with a basic plot
g3 <- ggplot(data.volt[which(data.volt$NORTHERLY==F), ])+
  geom_point(aes(y = PHENO.SLOPE, x = MARGIN.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g3

# now test the relationship using a linear model

model3 <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
               data = data.volt[which(data.volt$NORTHERLY==F), ])

summary(model3)
drop1(model3, test = "Chi")



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine
hist(data.uv[which(data.uv$NORTHERLY==F), ]$MARGIN.SLOPE)

model3.u <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.uv[which(data.uv$NORTHERLY==F), ])

summary(model3.u)
drop1(model3.u, test = "Chi")


# predict data to plot lines
newdata3u <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3u$MARGIN.SLOPE <- predict(model3.u, newdata = newdata3u, type="response")
newdata3u <- ddply(newdata3u, .(PHENO.SLOPE), numcolwise(mean))


g3.u <- ggplot(data.uv[which(data.uv$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata3u,
            aes(y = MARGIN.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-3,2)+
  ggtitle("Univoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g3.u


# multivoltine
hist(data.mv[which(data.mv$NORTHERLY==F), ]$MARGIN.SLOPE)

model3.m <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv[which(data.mv$NORTHERLY==F), ])

summary(model3.m)
drop1(model3.m, test = "Chi")


# predict data to plot lines
newdata3m <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3m$MARGIN.SLOPE <- predict(model3.m, newdata = newdata3m, type="response")
newdata3m <- ddply(newdata3m, .(PHENO.SLOPE), numcolwise(mean))


g3.m <- ggplot(data.mv[which(data.mv$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata3m,
            aes(y = MARGIN.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-3,2)+
  ggtitle("Multivoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g3.m




## combine and export the plots
m3 <- grid.arrange(g3.u,g3.m,
                   ncol=2)


leg3 <- g_legend(g3.u)

m3a <- grid.arrange(arrangeGrob(g3.u + theme(legend.position="none"),g3.m + theme(legend.position="none"),
                                ncol=2),
                    leg3,nrow=2, heights=c(10,2))

ggsave("Results/Plots/NorthernRangeMargin.png", plot = m3a, width = 350, height = 200, units = "mm", limitsize = F)



m1all <- grid.arrange(arrangeGrob(g1.u + theme(legend.position="none"),g1.m + theme(legend.position="none"),
                                  g2.u + theme(legend.position="none") + ggtitle(""),g2.m + theme(legend.position="none") + ggtitle(""),
                                  g3.u + theme(legend.position="none") + ggtitle(""),g3.m + theme(legend.position="none") + ggtitle(""),
                                  ncol=2),
                      leg1,ncol=2,widths=c(10,2))


ggsave("Results/Plots/Figure1.svg", plot = m1all, device = svg, width = 350, height = 400, units = "mm", limitsize = F)





## repeat for well- and heavily-recorded restricted data

## well-recorded

# read in butterfly data
butterfly.raw.well <- read.csv("Data/Derived/BMS_BNM_well.csv", header = T)
summary(butterfly.raw.well)

# and moth data
moth.raw.well <- read.csv("Data/Derived/RIS_NMRS_well.csv", header = T)
summary(moth.raw.well)

## merge together, labelling butterflies and moths
butterfly.raw.well$TAXON <- "Butterfly"
moth.raw.well$TAXON <- "Moth"

raw.data.well <- rbind(butterfly.raw.well,moth.raw.well)

### now we need to mark out which species are 'northerly-distributed' and can't therefore be used in analyses of Northern Range Margins
## this data is held in the 'species to use' tables and isn't changed

## merge in details of northerly species
data.dist.well <- merge(to.use,raw.data.well)

### we also want details of species' voltinism - again unchanged
data.volt.well <- merge(voltinism[,-1], data.dist.well)

data.volt.well$VOLTINISM <- as.factor(data.volt.well$VOLTINISM)
summary(data.volt.well)

## finally, for each of the four derived variables, let's generate a 'state' categorical variable for what the trend is
data.volt.well$PHENO.STATE <- as.factor(ifelse(data.volt.well$PHENO.P > 0.05, "Non-significant",
                                          ifelse(data.volt.well$PHENO.SLOPE < 0, "Significantly later emergence",
                                                 "Significantly earlier emergence")))

data.volt.well$ABUND.STATE <- as.factor(ifelse(data.volt.well$ABUND.P > 0.05, "Non-significant",
                                          ifelse(data.volt.well$ABUND.SLOPE < 0, "Significantly declining abundance",
                                                 "Significantly increasing abundance")))

data.volt.well$DISTRIB.STATE <- as.factor(ifelse(data.volt.well$DISTRIB.P > 0.05, "Non-significant",
                                            ifelse(data.volt.well$DISTRIB.SLOPE < 0, "Significantly retracting distribution",
                                                   "Significantly expanding distribution")))

data.volt.well$MARGIN.STATE <- as.factor(ifelse(data.volt.well$MARGIN.P > 0.05, "Non-significant",
                                           ifelse(data.volt.well$MARGIN.SLOPE < 0, "Significantly retracting southwards",
                                                  "Significantly expanding northwards")))


## divide it back into butterflies and moths
data.b.well <- data.volt.well[which(data.volt.well$TAXON=="Butterfly"), ]
data.m.well <- data.volt.well[which(data.volt.well$TAXON=="Moth"), ]


## and also divide it into univoltine and multivoltine
data.uv.well <- data.volt.well[which(data.volt.well$VOLTINISM==1), ]
data.mv.well <- data.volt.well[which(data.volt.well$VOLTINISM==2), ]

### analysis ####
## first let's plot up some histograms of each variable, coloured by significance

## phenology
hist1.well <- ggplot(data.volt.well)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("All Lepidoptera")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.well

# butterflies
hist1.b.well <- ggplot(data.b.well)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Butterflies")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.b.well


# moths
hist1.m.well <- ggplot(data.m.well)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Moths")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.m.well


## abundance
hist2.well <- ggplot(data.volt.well)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.well


# butterflies
hist2.b.well <- ggplot(data.b.well)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.b.well


# moths
hist2.m.well <- ggplot(data.m.well)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.m.well


## distribution
hist3.well <- ggplot(data.volt.well)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.well



# butterflies
hist3.b.well <- ggplot(data.b.well)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.b.well


# moths
hist3.m.well <- ggplot(data.m.well)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.m.well



## northern range margin
hist4.well <- ggplot(data.volt.well[which(data.volt.well$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,15)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.well

# butterflies
hist4.b.well <- ggplot(data.b.well[which(data.b.well$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,15)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.b.well


# moths
hist4.m.well <- ggplot(data.m.well[which(data.m.well$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,15)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.m.well


m.hist.well <- grid.arrange(hist1.b.well,hist1.well,hist1.m.well,
                       hist2.b.well,hist2.well,hist2.m.well,
                       hist3.b.well,hist3.well,hist3.m.well,
                       hist4.b.well,hist4.well,hist4.m.well,
                       ncol=3)


leg.hist1.well <- g_legend(hist1.b.well)
leg.hist2.well <- g_legend(hist2.b.well)
leg.hist3.well <- g_legend(hist3.b.well)
leg.hist4.well <- g_legend(hist4.b.well)

m.hist.a.well <- grid.arrange(arrangeGrob(hist1.b.well + theme(legend.position="none"),hist1.well + theme(legend.position="none"),hist1.m.well + theme(legend.position="none"),
                                     hist2.b.well + theme(legend.position="none"),hist2.well + theme(legend.position="none"),hist2.m.well + theme(legend.position="none"),
                                     hist3.b.well + theme(legend.position="none"),hist3.well + theme(legend.position="none"),hist3.m.well + theme(legend.position="none"),
                                     hist4.b.well + theme(legend.position="none"),hist4.well + theme(legend.position="none"),hist4.m.well + theme(legend.position="none"),
                                     ncol=3),
                         arrangeGrob(leg.hist1.well,leg.hist2.well,
                                     leg.hist3.well,leg.hist4.well,ncol =1),ncol=2, widths=c(10,2))

ggsave("Results/Plots/WellRecorded/Histograms.png", plot = m.hist.a.well, width = 450, height = 350, units = "mm", limitsize = T)

## and now some histograms making the same comparison but between univoltine and multivoltine species

## phenology
# uni
hist1.uv.well <- ggplot(data.uv.well)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Univoltine species")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.uv.well


# multi
hist1.mv.well <- ggplot(data.mv.well)+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  theme_bw()+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Multivoltine species")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.mv.well


## abundance
# univoltine
hist2.uv.well <- ggplot(data.uv.well)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.uv.well


# multivoltine
hist2.mv.well <- ggplot(data.mv.well)+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  theme_bw()+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.mv.well


## distribution
# univoltine
hist3.uv.well <- ggplot(data.uv.well)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.uv.well


# multivoltine
hist3.mv.well <- ggplot(data.mv.well)+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.mv.well



## northern range margin
# univoltine
hist4.uv.well <- ggplot(data.uv.well[which(data.uv.well$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.uv.well


# moths
hist4.mv.well <- ggplot(data.mv.well[which(data.mv.well$NORTHERLY==F), ])+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  theme_bw()+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.mv.well


m.hist.v.well <- grid.arrange(hist1.uv.well,hist1.well,hist1.mv.well,
                         hist2.uv.well,hist2.well,hist2.mv.well,
                         hist3.uv.well,hist3.well,hist3.mv.well,
                         hist4.uv.well,hist4.well,hist4.mv.well,
                         ncol=3)


leg.hist1v.well <- g_legend(hist1.uv.well)
leg.hist2v.well <- g_legend(hist2.uv.well)
leg.hist3v.well <- g_legend(hist3.uv.well)
leg.hist4v.well <- g_legend(hist4.uv.well)

m.hist.va.well <- grid.arrange(arrangeGrob(hist1.uv.well + theme(legend.position="none"),hist1.well + theme(legend.position="none"),hist1.mv.well + theme(legend.position="none"),
                                      hist2.uv.well + theme(legend.position="none"),hist2.well + theme(legend.position="none"),hist2.mv.well + theme(legend.position="none"),
                                      hist3.uv.well + theme(legend.position="none"),hist3.well + theme(legend.position="none"),hist3.mv.well + theme(legend.position="none"),
                                      hist4.uv.well + theme(legend.position="none"),hist4.well + theme(legend.position="none"),hist4.mv.well + theme(legend.position="none"),
                                      ncol=3),
                          arrangeGrob(leg.hist1v.well,leg.hist2v.well,
                                      leg.hist3v.well,leg.hist4v.well,ncol =1),ncol=2, widths=c(10,2))

ggsave("Results/Plots/WellRecorded/VoltinismHistograms.png", plot = m.hist.va.well, width = 350, height = 350, units = "mm", limitsize = T)


### now we are ready to start analysing relationships between phenology and other variables

## first, abundance
# start with a basic plot
g1.well <- ggplot(data.volt.well)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g1.well


# now test the relationship using a linear model

hist(data.volt.well$ABUND.EXP.SLOPE)

model1.well <- lmer(ABUND.EXP.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
               data = data.volt.well)

summary(model1.well)
drop1(model1.well, test = "Chi")



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine
hist(data.uv.well$ABUND.SLOPE)

model1.u.well <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.uv.well)

summary(model1.u.well)
drop1(model1.u.well, test = "Chi")


# predict data to plot lines
newdata1u.well <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1u.well$ABUND.SLOPE <- predict(model1.u.well, newdata = newdata1u.well, type="response")
newdata1u.well <- ddply(newdata1u.well, .(PHENO.SLOPE), numcolwise(mean))

g1.u.well <- ggplot(data.uv.well)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata1u.well,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-0.15,0.1)+
  xlim(-3,2)+
  ggtitle("Univoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g1.u.well


# multivoltine
hist(data.mv.well$ABUND.SLOPE)

model1.m.well <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv.well)

summary(model1.m.well)
drop1(model1.m.well, test = "Chi")


# predict data to plot lines
newdata1m.well <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1m.well$ABUND.SLOPE <- predict(model1.m.well, newdata = newdata1m.well, type="response")
newdata1m.well <- ddply(newdata1m.well, .(PHENO.SLOPE), numcolwise(mean))

g1.m.well <- ggplot(data.mv.well)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata1m.well,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-0.15,0.1)+
  xlim(-3,2)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g1.m.well



## combine and export the plots
m1.well <- grid.arrange(g1.u.well,g1.m.well,
                   ncol=2)


leg1.well <- g_legend(g1.u.well)

m1a.well <- grid.arrange(arrangeGrob(g1.u.well + theme(legend.position="none"),g1.m.well + theme(legend.position="none"),
                                ncol=2),
                    leg1.well,nrow=2, heights=c(10,2))

ggsave("Results/Plots/WellRecorded/Abundance.png", plot = m1a.well, width = 350, height = 200, units = "mm", limitsize = F)






## second, distribution
# start with a basic plot
g2.well <- ggplot(data.volt.well)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g2.well

# now test the relationship using a linear model

model2.well <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
               data = data.volt.well)

summary(model2.well)
drop1(model2.well, test = "Chi")



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine
hist(data.uv.well$DISTRIB.SLOPE)

model2.u.well <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.uv.well)

summary(model2.u.well)
drop1(model2.u.well, test = "Chi")


# predict data to plot lines
newdata2u.well <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2u.well$DISTRIB.SLOPE <- predict(model2.u.well, newdata = newdata2u.well, type="response")
newdata2u.well <- ddply(newdata2u.well, .(PHENO.SLOPE), numcolwise(mean))


g2.u.well <- ggplot(data.uv.well)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata2u.well,
            aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-3,2)+
  ggtitle("Univoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g2.u.well


# multivoltine
hist(data.mv.well$DISTRIB.SLOPE)

model2.m.well <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv.well)

summary(model2.m.well)
drop1(model2.m.well, test = "Chi")


# predict data to plot lines
newdata2m.well <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2m.well$DISTRIB.SLOPE <- predict(model2.m.well, newdata = newdata2m.well, type="response")
newdata2m.well <- ddply(newdata2m.well, .(PHENO.SLOPE), numcolwise(mean))

g2.m.well <- ggplot(data.mv.well)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata2m.well,
            aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-3,2)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g2.m.well


## combine and export the plots
m2.well <- grid.arrange(g2.u.well,g2.m.well,
                   ncol=2)


leg2.well <- g_legend(g2.u.well)

m2a.well <- grid.arrange(arrangeGrob(g2.u.well + theme(legend.position="none"),g2.m.well + theme(legend.position="none"),
                                ncol=2),
                    leg2.well,nrow=2, heights=c(10,2))

ggsave("Results/Plots/WellRecorded/Distribution.png", plot = m2a.well, width = 350, height = 200, units = "mm", limitsize = F)




## third, NRM
# start with a basic plot
g3.well <- ggplot(data.volt.well[which(data.volt.well$NORTHERLY==F), ])+
  geom_point(aes(y = PHENO.SLOPE, x = MARGIN.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g3.well

# now test the relationship using a linear model

model3.well <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
               data = data.volt.well[which(data.volt.well$NORTHERLY==F), ])

summary(model3.well)
drop1(model3.well, test = "Chi")



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine
hist(data.uv.well[which(data.uv.well$NORTHERLY==F), ]$MARGIN.SLOPE)

model3.u.well <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.uv.well[which(data.uv.well$NORTHERLY==F), ])

summary(model3.u.well)
drop1(model3.u.well, test = "Chi")


# predict data to plot lines
newdata3u.well <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3u.well$MARGIN.SLOPE <- predict(model3.u.well, newdata = newdata3u.well, type="response")
newdata3u.well <- ddply(newdata3u.well, .(PHENO.SLOPE), numcolwise(mean))


g3.u.well <- ggplot(data.uv.well[which(data.uv.well$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata3u.well,
            aes(y = MARGIN.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype="dashed")+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-3,2)+
  ggtitle("Univoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g3.u.well


# multivoltine
hist(data.mv.well[which(data.mv.well$NORTHERLY==F), ]$MARGIN.SLOPE)

model3.m.well <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv.well[which(data.mv.well$NORTHERLY==F), ])

summary(model3.m.well)
drop1(model3.m.well, test = "Chi")


# predict data to plot lines
newdata3m.well <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3m.well$MARGIN.SLOPE <- predict(model3.m.well, newdata = newdata3m.well, type="response")
newdata3m.well <- ddply(newdata3m.well, .(PHENO.SLOPE), numcolwise(mean))


g3.m.well <- ggplot(data.mv.well[which(data.mv.well$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata3m.well,
            aes(y = MARGIN.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-3,2)+
  ggtitle("Multivoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g3.m.well




## combine and export the plots
m3.well <- grid.arrange(g3.u.well,g3.m.well,
                   ncol=2)


leg3.well <- g_legend(g3.u.well)

m3a.well <- grid.arrange(arrangeGrob(g3.u.well + theme(legend.position="none"),g3.m.well + theme(legend.position="none"),
                                ncol=2),
                    leg3.well,nrow=2, heights=c(10,2))

ggsave("Results/Plots/WellRecorded/NorthernRangeMargin.png", plot = m3a.well, width = 350, height = 200, units = "mm", limitsize = F)



## heavily-recorded ####

#### important: this is the level of recording we will report in the final manuscript,
# so I may make refinements to the code from here onwards that are not fully replicated above

# read in butterfly data
butterfly.raw.heavy <- read.csv("Data/Derived/BMS_BNM_heavy.csv", header = T)
summary(butterfly.raw.heavy)

# and moth data
moth.raw.heavy <- read.csv("Data/Derived/RIS_NMRS_heavy.csv", header = T)
summary(moth.raw.heavy)

## merge together, labelling butterflies and moths
butterfly.raw.heavy$TAXON <- "Butterfly"
moth.raw.heavy$TAXON <- "Moth"

raw.data.heavy <- rbind(butterfly.raw.heavy,moth.raw.heavy)

### now we need to mark out which species are 'northerly-distributed' and can't therefore be used in analyses of Northern Range Margins
## this data is held in the 'species to use' tables and isn't changed

## merge in details of northerly species
data.dist.heavy <- merge(to.use,raw.data.heavy)

### we also want details of species' voltinism - again unchanged
data.volt.heavy <- merge(voltinism[,-1], data.dist.heavy)

data.volt.heavy$VOLTINISM <- as.factor(data.volt.heavy$VOLTINISM)
summary(data.volt.heavy)

## finally, for each of the four derived variables, let's generate a 'state' categorical variable for what the trend is
data.volt.heavy$PHENO.STATE <- as.factor(ifelse(data.volt.heavy$PHENO.P > 0.05, "Non-significant",
                                               ifelse(data.volt.heavy$PHENO.SLOPE < 0, "Significantly later emergence",
                                                      "Significantly earlier emergence")))

data.volt.heavy$ABUND.STATE <- as.factor(ifelse(data.volt.heavy$ABUND.P > 0.05, "Non-significant",
                                               ifelse(data.volt.heavy$ABUND.SLOPE < 0, "Significantly declining abundance",
                                                      "Significantly increasing abundance")))

data.volt.heavy$DISTRIB.STATE <- as.factor(ifelse(data.volt.heavy$DISTRIB.P > 0.05, "Non-significant",
                                                 ifelse(data.volt.heavy$DISTRIB.SLOPE < 0, "Significantly retracting distribution",
                                                        "Significantly expanding distribution")))

data.volt.heavy$MARGIN.STATE <- as.factor(ifelse(data.volt.heavy$MARGIN.P > 0.05, "Non-significant",
                                                ifelse(data.volt.heavy$MARGIN.SLOPE < 0, "Significantly retracting southwards",
                                                       "Significantly expanding northwards")))


data.volt.heavy$TAXON <- as.factor(data.volt.heavy$TAXON)

## divide it back into butterflies and moths
data.b.heavy <- data.volt.heavy[which(data.volt.heavy$TAXON=="Butterfly"), ]
data.m.heavy <- data.volt.heavy[which(data.volt.heavy$TAXON=="Moth"), ]


## and also divide it into univoltine and multivoltine
data.uv.heavy <- data.volt.heavy[which(data.volt.heavy$VOLTINISM==1), ]
data.mv.heavy <- data.volt.heavy[which(data.volt.heavy$VOLTINISM==2), ]

## and also divide it into specialists and generalists
data.hs.heavy <- data.volt.heavy[which(data.volt.heavy$CLASS=="HS"), ]
data.wc.heavy <- data.volt.heavy[which(data.volt.heavy$CLASS=="WC"), ]

# also export it, as we'll want this in as a supplementary table
# but we don't want the VOLTI columns in there, so let's ditch those out first
data.volt.out <- data.volt.heavy[,c(30,1:3,41,4:6,8:24,31:40)]

# now write it
write.csv(data.volt.out, file = "Results/AllSpecies.csv", row.names = F)


## before starting to analyse, let's summarise what data we have from each scheme

# number of records (i.e. one species at one site on one day)
sum(data.b.heavy$RECORDS)

# number of populations (i.e. one species at one site)
sum(data.b.heavy$SITES)
sum(data.m.heavy$SITES)

# number of individuals (sum of all abundance counts)
sum(data.b.heavy$INDIVIDUALS)
sum(data.m.heavy$INDIVIDUALS)

# number of hectad-level presence records
sum(data.b.heavy$HECTAD.RECORDS)
sum(data.m.heavy$HECTAD.RECORDS)


### analysis ####
## first let's plot up some histograms of each variable, coloured by significance

# before plotting each histogram, let's extract a mean and s.d. of change
# and do a t-test for a difference between butterflies and moths

## phenology
mean(data.volt.heavy$PHENO.SLOPE)
sd(data.volt.heavy$PHENO.SLOPE)/sqrt(nrow(data.volt.heavy))

t.test(data.uv.heavy$PHENO.SLOPE)
t.test(data.mv.heavy$PHENO.SLOPE)

t.test(data.volt.heavy$PHENO.SLOPE)

t.test(data.b.heavy$PHENO.SLOPE, data.m.heavy$PHENO.SLOPE)
t.test(data.uv.heavy$PHENO.SLOPE, data.mv.heavy$PHENO.SLOPE)

shapiro.test(data.b.heavy$PHENO.SLOPE)
shapiro.test(data.m.heavy$PHENO.SLOPE)
shapiro.test(data.uv.heavy$PHENO.SLOPE)
shapiro.test(data.mv.heavy$PHENO.SLOPE)


wilcox.test(data.volt.heavy$PHENO.SLOPE, mu = 0, alternative = "two.sided")

wilcox.test(data.b.heavy$PHENO.SLOPE, data.m.heavy$PHENO.SLOPE)
wilcox.test(data.uv.heavy$PHENO.SLOPE, data.mv.heavy$PHENO.SLOPE)
wilcox.test(data.hs.heavy$PHENO.SLOPE, data.wc.heavy$PHENO.SLOPE)



## phenology
hist1.heavy <- ggplot(data.volt.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("All Lepidoptera")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.heavy

# butterflies
hist1.b.heavy <- ggplot(data.b.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+ 
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Butterflies")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text=element_text(size=11))

hist1.b.heavy


# moths
hist1.m.heavy <- ggplot(data.m.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Moths")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.m.heavy


## abundance
mean(data.volt.heavy$ABUND.SLOPE)
sd(data.volt.heavy$ABUND.SLOPE)/sqrt(nrow(data.volt.heavy))

t.test(data.volt.heavy$ABUND.SLOPE)

t.test(data.b.heavy$ABUND.SLOPE, data.m.heavy$ABUND.SLOPE)
t.test(data.uv.heavy$ABUND.SLOPE, data.mv.heavy$ABUND.SLOPE)

shapiro.test(data.b.heavy$ABUND.SLOPE)
shapiro.test(data.m.heavy$ABUND.SLOPE)
shapiro.test(data.uv.heavy$ABUND.SLOPE)
shapiro.test(data.mv.heavy$ABUND.SLOPE)


wilcox.test(data.volt.heavy$ABUND.SLOPE, mu = 0, alternative = "two.sided")

wilcox.test(data.b.heavy$ABUND.SLOPE, data.m.heavy$ABUND.SLOPE)
wilcox.test(data.uv.heavy$ABUND.SLOPE, data.mv.heavy$ABUND.SLOPE)
wilcox.test(data.hs.heavy$ABUND.SLOPE, data.wc.heavy$ABUND.SLOPE)


hist2.heavy <- ggplot(data.volt.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.heavy


# butterflies
hist2.b.heavy <- ggplot(data.b.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")+
  theme(legend.text=element_text(size=11))

hist2.b.heavy


# moths
hist2.m.heavy <- ggplot(data.m.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.m.heavy


## distribution
mean(data.volt.heavy$DISTRIB.SLOPE)
sd(data.volt.heavy$DISTRIB.SLOPE)/sqrt(nrow(data.volt.heavy))

t.test(data.volt.heavy$DISTRIB.SLOPE)

t.test(data.b.heavy$DISTRIB.SLOPE, data.m.heavy$DISTRIB.SLOPE)
t.test(data.uv.heavy$DISTRIB.SLOPE, data.mv.heavy$DISTRIB.SLOPE)


shapiro.test(data.b.heavy$DISTRIB.SLOPE)
shapiro.test(data.m.heavy$DISTRIB.SLOPE)
shapiro.test(data.uv.heavy$DISTRIB.SLOPE)
shapiro.test(data.mv.heavy$DISTRIB.SLOPE)


wilcox.test(data.volt.heavy$DISTRIB.SLOPE, mu = 0, alternative = "two.sided")

wilcox.test(data.b.heavy$DISTRIB.SLOPE, data.m.heavy$DISTRIB.SLOPE)
wilcox.test(data.uv.heavy$DISTRIB.SLOPE, data.mv.heavy$DISTRIB.SLOPE)
wilcox.test(data.hs.heavy$DISTRIB.SLOPE, data.wc.heavy$DISTRIB.SLOPE)



hist3.heavy <- ggplot(data.volt.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(% recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.heavy



# butterflies
hist3.b.heavy <- ggplot(data.b.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")+
  theme(legend.text=element_text(size=11))

hist3.b.heavy


# moths
hist3.m.heavy <- ggplot(data.m.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.m.heavy



## northern range margin
mean(data.volt.heavy$MARGIN.SLOPE)
sd(data.volt.heavy$MARGIN.SLOPE)/sqrt(nrow(data.volt.heavy))

t.test(data.volt.heavy$MARGIN.SLOPE)

t.test(data.b.heavy$MARGIN.SLOPE, data.m.heavy$MARGIN.SLOPE)
t.test(data.uv.heavy$MARGIN.SLOPE, data.mv.heavy$MARGIN.SLOPE)


shapiro.test(data.b.heavy$MARGIN.SLOPE)
shapiro.test(data.m.heavy$MARGIN.SLOPE)
shapiro.test(data.uv.heavy$MARGIN.SLOPE)
shapiro.test(data.mv.heavy$MARGIN.SLOPE)


wilcox.test(data.volt.heavy$MARGIN.SLOPE, mu = 0, alternative = "two.sided")

wilcox.test(data.b.heavy$MARGIN.SLOPE, data.m.heavy$MARGIN.SLOPE)
wilcox.test(data.uv.heavy$MARGIN.SLOPE, data.mv.heavy$MARGIN.SLOPE)
wilcox.test(data.hs.heavy$MARGIN.SLOPE, data.wc.heavy$MARGIN.SLOPE)



hist4.heavy <- ggplot(data.volt.heavy[which(data.volt.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,15)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.heavy

# butterflies
hist4.b.heavy <- ggplot(data.b.heavy[which(data.b.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,15)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")+
  theme(legend.text=element_text(size=11))

hist4.b.heavy


# moths
hist4.m.heavy <- ggplot(data.m.heavy[which(data.m.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,15)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.m.heavy


m.hist.heavy <- grid.arrange(hist1.b.heavy,hist1.heavy,hist1.m.heavy,
                            hist2.b.heavy,hist2.heavy,hist2.m.heavy,
                            hist3.b.heavy,hist3.heavy,hist3.m.heavy,
                            hist4.b.heavy,hist4.heavy,hist4.m.heavy,
                            ncol=3)


leg.hist1.heavy <- g_legend(hist1.b.heavy)
leg.hist2.heavy <- g_legend(hist2.b.heavy)
leg.hist3.heavy <- g_legend(hist3.b.heavy)
leg.hist4.heavy <- g_legend(hist4.b.heavy)

m.hist.a.heavy <- grid.arrange(arrangeGrob(hist1.b.heavy + theme(legend.position="none") + xlab("\n"),hist1.heavy + theme(legend.position="none") + ylab(" "),hist1.m.heavy + theme(legend.position="none") + xlab("\n") + ylab(" "),
                                          hist2.b.heavy + theme(legend.position="none") + xlab("\n"),hist2.heavy + theme(legend.position="none") + ylab(" "),hist2.m.heavy + theme(legend.position="none") + xlab("\n") + ylab(" "),
                                          hist3.b.heavy + theme(legend.position="none") + xlab("\n"),hist3.heavy + theme(legend.position="none") + ylab(" "),hist3.m.heavy + theme(legend.position="none") + xlab("\n") + ylab(" "),
                                          hist4.b.heavy + theme(legend.position="none") + xlab("\n"),hist4.heavy + theme(legend.position="none") + ylab(" "),hist4.m.heavy + theme(legend.position="none") + xlab("\n") + ylab(" "),
                                          ncol=3),
                              arrangeGrob(leg.hist1.heavy,leg.hist2.heavy,
                                          leg.hist3.heavy,leg.hist4.heavy,ncol =1),ncol=2, widths=c(10,3))

ggsave("Results/Plots/Final/FigB.svg", device=svg, plot = m.hist.a.heavy, width = 380, height = 350, units = "mm", limitsize = T)

## and now some histograms making the same comparison but between univoltine and multivoltine species

## phenology
# uni
hist1.uv.heavy <- ggplot(data.uv.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Univoltine species")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.text=element_text(size=11))

hist1.uv.heavy


# multi
hist1.mv.heavy <- ggplot(data.mv.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Multivoltine species")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.mv.heavy


## abundance
# univoltine
hist2.uv.heavy <- ggplot(data.uv.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")+
  theme(legend.text=element_text(size=11))

hist2.uv.heavy


# multivoltine
hist2.mv.heavy <- ggplot(data.mv.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.mv.heavy


## distribution
# univoltine
hist3.uv.heavy <- ggplot(data.uv.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")+
  theme(legend.text=element_text(size=11))

hist3.uv.heavy


# multivoltine
hist3.mv.heavy <- ggplot(data.mv.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.mv.heavy



## northern range margin
# univoltine
hist4.uv.heavy <- ggplot(data.uv.heavy[which(data.uv.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")+
  theme(legend.text=element_text(size=11))

hist4.uv.heavy


# moths
hist4.mv.heavy <- ggplot(data.mv.heavy[which(data.mv.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.mv.heavy


m.hist.v.heavy <- grid.arrange(hist1.uv.heavy,hist1.heavy,hist1.mv.heavy,
                              hist2.uv.heavy,hist2.heavy,hist2.mv.heavy,
                              hist3.uv.heavy,hist3.heavy,hist3.mv.heavy,
                              hist4.uv.heavy,hist4.heavy,hist4.mv.heavy,
                              ncol=3)


leg.hist1v.heavy <- g_legend(hist1.uv.heavy)
leg.hist2v.heavy <- g_legend(hist2.uv.heavy)
leg.hist3v.heavy <- g_legend(hist3.uv.heavy)
leg.hist4v.heavy <- g_legend(hist4.uv.heavy)

m.hist.va.heavy <- grid.arrange(arrangeGrob(hist1.uv.heavy + theme(legend.position="none") + xlab("\n"),hist1.heavy + theme(legend.position="none") + ylab(""),hist1.mv.heavy + theme(legend.position="none") + xlab("\n") + ylab(""),
                                           hist2.uv.heavy + theme(legend.position="none") + xlab("\n"),hist2.heavy + theme(legend.position="none") + ylab(""),hist2.mv.heavy + theme(legend.position="none") + xlab("\n") + ylab(""),
                                           hist3.uv.heavy + theme(legend.position="none") + xlab("\n"),hist3.heavy + theme(legend.position="none") + ylab(""),hist3.mv.heavy + theme(legend.position="none") + xlab("\n") + ylab(""),
                                           hist4.uv.heavy + theme(legend.position="none") + xlab("\n"),hist4.heavy + theme(legend.position="none") + ylab(""),hist4.mv.heavy + theme(legend.position="none") + xlab("\n") + ylab(""),
                                           ncol=3),
                               arrangeGrob(leg.hist1v.heavy,leg.hist2v.heavy,
                                           leg.hist3v.heavy,leg.hist4v.heavy,ncol =1),ncol=2, widths=c(10,3))

ggsave("Results/Plots/Final/FigC.svg", device = svg, plot = m.hist.va.heavy, width = 380, height = 350, units = "mm", limitsize = T)

## and finally some histograms making the same comparison but between specialist and generalist species

## phenology
# hs
hist1.hs.heavy <- ggplot(data.hs.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+  
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Habitat specialists")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.hs.heavy


# wc
hist1.wc.heavy <- ggplot(data.wc.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+  
  geom_histogram(aes(x = PHENO.SLOPE, fill = PHENO.STATE), bins=15)+
  xlim(-4,4)+
  ylim(0,75)+
  ggtitle("Wider countryside generalists")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Emergence date change")+
  theme(plot.title = element_text(hjust = 0.5))

hist1.wc.heavy


## abundance
# hs
hist2.hs.heavy <- ggplot(data.hs.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+  
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.hs.heavy


# wc
hist2.wc.heavy <- ggplot(data.wc.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = ABUND.SLOPE, fill = ABUND.STATE), bins=15)+
  xlim(-0.25,0.25)+
  ylim(0,75)+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Abundance change")

hist2.wc.heavy


## distribution
# hs
hist3.hs.heavy <- ggplot(data.hs.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.hs.heavy


# wc
hist3.wc.heavy <- ggplot(data.wc.heavy)+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = DISTRIB.SLOPE, fill = DISTRIB.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,75)+
  xlab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Distribution change")

hist3.wc.heavy



## northern range margin
# hs
hist4.hs.heavy <- ggplot(data.hs.heavy[which(data.hs.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.hs.heavy


# wc
hist4.wc.heavy <- ggplot(data.wc.heavy[which(data.wc.heavy$NORTHERLY==F), ])+
  geom_hline(aes(yintercept = 0), colour = "grey80")+
  geom_vline(aes(xintercept = 0), colour = "grey80")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_histogram(aes(x = MARGIN.SLOPE, fill = MARGIN.STATE), bins=15)+
  xlim(-20,20)+
  ylim(0,12)+
  xlab("Change in northern range margin \n(km northwards per year)")+
  ylab("Number of species")+
  scale_fill_brewer(palette = "Paired", drop = F, name = "Range margin change")

hist4.wc.heavy


m.hist.s.heavy <- grid.arrange(hist1.hs.heavy,hist1.heavy,hist1.wc.heavy,
                         hist2.hs.heavy,hist2.heavy,hist2.wc.heavy,
                         hist3.hs.heavy,hist3.heavy,hist3.wc.heavy,
                         hist4.hs.heavy,hist4.heavy,hist4.wc.heavy,
                         ncol=3)


leg.hist1s.heavy <- g_legend(hist1.hs.heavy)
leg.hist2s.heavy <- g_legend(hist2.hs.heavy)
leg.hist3s.heavy <- g_legend(hist3.hs.heavy)
leg.hist4s.heavy <- g_legend(hist4.hs.heavy)

m.hist.sa.heavy <- grid.arrange(arrangeGrob(hist1.hs.heavy + theme(legend.position="none"),hist1.heavy + theme(legend.position="none"),hist1.wc.heavy + theme(legend.position="none"),
                                      hist2.hs.heavy + theme(legend.position="none"),hist2.heavy + theme(legend.position="none"),hist2.wc.heavy + theme(legend.position="none"),
                                      hist3.hs.heavy + theme(legend.position="none"),hist3.heavy + theme(legend.position="none"),hist3.wc.heavy + theme(legend.position="none"),
                                      hist4.hs.heavy + theme(legend.position="none"),hist4.heavy + theme(legend.position="none"),hist4.wc.heavy + theme(legend.position="none"),
                                      ncol=3),
                          arrangeGrob(leg.hist1s.heavy,leg.hist2s.heavy,
                                      leg.hist3s.heavy,leg.hist4s.heavy,ncol =1),ncol=2, widths=c(10,2))

ggsave("Results/Plots/Final/FigC1.svg", plot = m.hist.sa.heavy, device = svg, width = 350, height = 350, units = "mm", limitsize = T)



### now we are ready to start analysing relationships between phenology and other variables

## first, abundance
# start with a basic plot
g1.heavy <- ggplot(data.volt.heavy)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g1.heavy

# do a version using facet wrap to show the four treatments separately
g1.heavy.wrap <- ggplot(data.volt.heavy)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()+
  facet_grid(CLASS ~ VOLTINISM)

g1.heavy.wrap



# now test the relationship using a linear model

hist(data.volt.heavy$ABUND.SLOPE)

model1.heavy.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                    data = data.volt.heavy)

summary(model1.heavy.c)
drop1(model1.heavy.c, test = "Chi")

r.squaredGLMM(model1.heavy.c)

# without class

hist(data.volt.heavy$ABUND.SLOPE)

model1.heavy <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                     data = data.volt.heavy)

summary(model1.heavy)
drop1(model1.heavy, test = "Chi")

r.squaredGLMM(model1.heavy)

# add the lines to that figure
newdata1.heavy <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"), VOLTINISM = factor(c(1,2)), CLASS = factor(c("HS","WC")), ABUND.SLOPE = 0)
newdata1.heavy$ABUND.SLOPE <- predict(model1.heavy, newdata = newdata1.heavy, type="response")
newdata1.heavy <- ddply(newdata1.heavy, .(PHENO.SLOPE,VOLTINISM,CLASS), numcolwise(mean))



g1.heavy <- ggplot(data.volt.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  geom_line(data = newdata1.heavy,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, linetype = VOLTINISM, colour = CLASS),
            stat="identity",size=1)+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  scale_linetype_manual(values = c("dashed","solid"))+
  theme_bw()

g1.heavy


## split the analysis up into the various four classes, and add the model lines
data.uv.hs <- data.uv.heavy[which(data.uv.heavy$CLASS == "HS"), ]
data.mv.hs <- data.mv.heavy[which(data.mv.heavy$CLASS == "HS"), ]
data.uv.wc <- data.uv.heavy[which(data.uv.heavy$CLASS == "WC"), ]
data.mv.wc <- data.mv.heavy[which(data.mv.heavy$CLASS == "WC"), ]


# univoltine, hs
hist(data.uv.hs$ABUND.SLOPE)

model1.u.hs <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                      data = data.uv.hs)

summary(model1.u.hs)
drop1(model1.u.hs, test = "Chi")

chkres(model1.u.hs, data.uv.hs$PHENO.SLOPE, data.uv.hs$TAXON)


# predict data to plot lines
newdata1u.hs <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1u.hs$ABUND.SLOPE <- predict(model1.u.hs, newdata = newdata1u.hs, type="response")
newdata1u.hs <- ddply(newdata1u.hs, .(PHENO.SLOPE), numcolwise(mean))
newdata1u.hs$CLASS <- "HS"

# univoltine, wc
hist(data.uv.wc$ABUND.SLOPE)

model1.u.wc <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.uv.wc)

summary(model1.u.wc)
drop1(model1.u.wc, test = "Chi")

chkres(model1.u.wc, data.uv.wc$PHENO.SLOPE, data.uv.wc$TAXON)

## no class
# univoltine
hist(data.uv.heavy$ABUND.SLOPE)

model1.u.heavy <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.uv.heavy)

summary(model1.u.heavy)
drop1(model1.u.heavy, test = "Chi")

chkres(model1.u.heavy, data.uv.heavy$PHENO.SLOPE, data.uv.heavy$TAXON)




# predict data to plot lines
newdata1u.heavy <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1u.heavy$ABUND.SLOPE <- predict(model1.u.heavy, newdata = newdata1u.heavy, type="response")
newdata1u.heavy <- ddply(newdata1u.heavy, .(PHENO.SLOPE), numcolwise(mean))

g1.u.heavy <- ggplot(data.uv.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-0.15,0.15)+
  xlim(-1,2)+
  ggtitle("Univoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g1.u.heavy


newdata1u.wc <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1u.wc$ABUND.SLOPE <- predict(model1.u.wc, newdata = newdata1u.wc, type="response")
newdata1u.wc <- ddply(newdata1u.wc, .(PHENO.SLOPE), numcolwise(mean))
newdata1u.wc$CLASS <- "WC"

newdata1u.all <- rbind(newdata1u.hs,newdata1u.wc)

g1.u.heavy.c <- ggplot(data.uv.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata1u.all,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, colour = CLASS),
            linetype = "dotted", stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  scale_colour_manual(values = c("Black",NA), name = "Class")+
  theme_bw()+
  ylim(-0.15,0.15)+
  xlim(-1,2)+
  ggtitle("Univoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g1.u.heavy.c


# multivoltine
# multivoltine, hs
hist(data.mv.hs$ABUND.SLOPE)

model1.m.hs <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.mv.hs)

summary(model1.m.hs)
drop1(model1.m.hs, test = "Chi")

chkres(model1.m.hs, data.mv.hs$PHENO.SLOPE, data.mv.hs$TAXON)


# predict data to plot lines
newdata1m.hs <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1m.hs$ABUND.SLOPE <- predict(model1.m.hs, newdata = newdata1m.hs, type="response")
newdata1m.hs <- ddply(newdata1m.hs, .(PHENO.SLOPE), numcolwise(mean))
newdata1m.hs$CLASS <- "HS"

# multivoltine, wc
hist(data.mv.wc$ABUND.SLOPE)

model1.m.wc <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.mv.wc)

summary(model1.m.wc)
drop1(model1.m.wc, test = "Chi")

chkres(model1.m.wc, data.mv.wc$PHENO.SLOPE, data.mv.wc$TAXON)

## no class
# multivoltine
hist(data.mv.heavy$ABUND.SLOPE)

model1.m.heavy <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                 data = data.mv.heavy)

summary(model1.m.heavy)
drop1(model1.m.heavy, test = "Chi")

chkres(model1.m.heavy, data.mv.heavy$PHENO.SLOPE, data.mv.heavy$TAXON)

# predict data to plot lines
newdata1m.heavy <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1m.heavy$ABUND.SLOPE <- predict(model1.m.heavy, newdata = newdata1m.heavy, type="response")
newdata1m.heavy <- ddply(newdata1m.heavy, .(PHENO.SLOPE), numcolwise(mean))

g1.m.heavy <- ggplot(data.mv.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata1m.heavy,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-0.15,0.15)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g1.m.heavy




# predict data to plot lines
newdata1m.wc <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),ABUND.SLOPE = 0)
newdata1m.wc$ABUND.SLOPE <- predict(model1.m.wc, newdata = newdata1m.wc, type="response")
newdata1m.wc <- ddply(newdata1m.wc, .(PHENO.SLOPE), numcolwise(mean))
newdata1m.wc$CLASS <- "WC"

newdata1m.all <- rbind(newdata1m.hs,newdata1m.wc)

g1.m.heavy.c <- ggplot(data.mv.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata1m.wc,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            linetype = "dashed", stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  scale_linetype_manual(values = c("dotted","dashed"), name = "Class")+
  theme_bw()+
  ylim(-0.15,0.15)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g1.m.heavy.c

g1.m.dummy <- ggplot(data.mv.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata1m.all,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, linetype = CLASS),
            stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  scale_linetype_manual(values = c("dotted","dashed"), name = "Class")+
  theme_bw()+
  ylim(-0.15,0.15)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
g1.m.dummy

## combine and export the plots
m1.heavy <- grid.arrange(g1.u.heavy,g1.m.heavy,
                        ncol=2)


leg1.heavy <- g_legend(g1.u.heavy)

m1a.heavy <- grid.arrange(arrangeGrob(g1.u.heavy + theme(legend.position="none"),g1.m.heavy + theme(legend.position="none"),
                                     ncol=2),
                         leg1.heavy,nrow=2, heights=c(10,2))

ggsave("Results/Plots/HeavilyRecorded/Abundance.png", plot = m1a.heavy, width = 350, height = 200, units = "mm", limitsize = F)


## before moving on, I want to test the influence of a possible outlier - the univoltine habitat-specialist Miltochrista miniata (Rosy Footman)

# trim it from the dataframe
data.volt.heavy.out <- data.volt.heavy[which(data.volt.heavy$COMMON_NAME != "Rosy Footman"), ]

g1.heavy.wrap.out <- ggplot(data.volt.heavy.out)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()+
  facet_grid(CLASS ~ VOLTINISM)

g1.heavy.wrap.out


model1.heavy.out.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                     data = data.volt.heavy.out)

summary(model1.heavy.out.c)
drop1(model1.heavy.out.c, test = "Chi")

model1.heavy.out <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                         data = data.volt.heavy.out)

summary(model1.heavy.out)
drop1(model1.heavy.out, test = "Chi")


# and check within the UVHS

data.uv.hs.out <- data.uv.hs[which(data.uv.hs$COMMON_NAME != "Rosy Footman"), ]


model1.u.hs.out <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.uv.hs.out)

summary(model1.u.hs.out)
drop1(model1.u.hs.out, test = "Chi")


## second, distribution
# start with a basic plot
g2.heavy <- ggplot(data.volt.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g2.heavy


# now test the relationship using a linear model

hist(data.volt.heavy$DISTRIB.SLOPE)

model2.heavy.c <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                     data = data.volt.heavy)

summary(model2.heavy.c)
drop1(model2.heavy.c, test = "Chi")

r.squaredGLMM(model2.heavy.c)


model2.heavy <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + VOLTINISM + (1|TAXON),
                     data = data.volt.heavy)

summary(model2.heavy)
drop1(model2.heavy, test = "Chi")

r.squaredGLMM(model2.heavy)


# add the lines to that figure
newdata2.heavy <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"), VOLTINISM = factor(c(1,2)), CLASS = factor(c("HS","WC")), DISTRIB.SLOPE = 0)
newdata2.heavy$DISTRIB.SLOPE <- predict(model2.heavy, newdata = newdata2.heavy, type="response")
newdata2.heavy <- ddply(newdata2.heavy, .(PHENO.SLOPE,VOLTINISM,CLASS), numcolwise(mean))



g2.heavy <- ggplot(data.volt.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  geom_line(data = newdata2.heavy,
            aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, linetype = VOLTINISM, colour = CLASS),
            stat="identity",size=1)+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  scale_linetype_manual(values = c("dashed","solid"))+
  theme_bw()

g2.heavy



# univoltine, hs
hist(data.uv.hs$DISTRIB.SLOPE)

model2.u.hs <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.uv.hs)

summary(model2.u.hs)
drop1(model2.u.hs, test = "Chi")

chkres(model2.u.hs, data.uv.hs$PHENO.SLOPE, data.uv.hs$TAXON)


# predict data to plot lines
newdata2u.hs <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2u.hs$DISTRIB.SLOPE <- predict(model2.u.hs, newdata = newdata2u.hs, type="response")
newdata2u.hs <- ddply(newdata2u.hs, .(PHENO.SLOPE), numcolwise(mean))
newdata2u.hs$CLASS <- "HS"

# univoltine, wc
hist(data.uv.wc$DISTRIB.SLOPE)

model2.u.wc <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.uv.wc)

summary(model2.u.wc)
drop1(model2.u.wc, test = "Chi")

chkres(model2.u.wc, data.uv.wc$PHENO.SLOPE, data.uv.wc$TAXON)


# predict data to plot lines
newdata2u.wc <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2u.wc$DISTRIB.SLOPE <- predict(model2.u.wc, newdata = newdata2u.wc, type="response")
newdata2u.wc <- ddply(newdata2u.wc, .(PHENO.SLOPE), numcolwise(mean))
newdata2u.wc$CLASS <- "WC"

newdata2u.all <- rbind(newdata2u.hs,newdata2u.wc)

g2.u.heavy.c <- ggplot(data.uv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-1,2)+
  ggtitle("Univoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g2.u.heavy.c


g2.u.heavy <- ggplot(data.uv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-1,2)+
  ggtitle("Univoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g2.u.heavy


# multivoltine
# multivoltine, hs
hist(data.mv.hs$DISTRIB.SLOPE)

model2.m.hs <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.mv.hs)

summary(model2.m.hs)
drop1(model2.m.hs, test = "Chi")

chkres(model2.m.hs, data.mv.hs$PHENO.SLOPE, data.mv.hs$TAXON)


# predict data to plot lines
newdata2m.hs <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2m.hs$DISTRIB.SLOPE <- predict(model2.m.hs, newdata = newdata2m.hs, type="response")
newdata2m.hs <- ddply(newdata2m.hs, .(PHENO.SLOPE), numcolwise(mean))
newdata2m.hs$CLASS <- "HS"

# multivoltine, wc
hist(data.mv.wc$DISTRIB.SLOPE)

model2.m.wc <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.mv.wc)

summary(model2.m.wc)
drop1(model2.m.wc, test = "Chi")

chkres(model2.m.wc, data.mv.wc$PHENO.SLOPE, data.mv.wc$TAXON)


# predict data to plot lines
newdata2m.wc <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata2m.wc$DISTRIB.SLOPE <- predict(model2.m.wc, newdata = newdata2m.wc, type="response")
newdata2m.wc <- ddply(newdata2m.wc, .(PHENO.SLOPE), numcolwise(mean))
newdata2m.wc$CLASS <- "WC"

newdata2m.all <- rbind(newdata2m.hs,newdata2m.wc)

g2.m.heavy.c <- ggplot(data.mv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g2.m.heavy.c


g2.m.heavy <- ggplot(data.mv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g2.m.heavy



## combine and export the plots
m2.heavy <- grid.arrange(g2.u.heavy,g2.m.heavy,
                         ncol=2)


leg2.heavy <- g_legend(g2.m.heavy)

m2a.heavy <- grid.arrange(arrangeGrob(g2.u.heavy + theme(legend.position="none"),g2.m.heavy + theme(legend.position="none"),
                                      ncol=2),
                          leg2.heavy,nrow=2, heights=c(10,2))

ggsave("Results/Plots/HeavilyRecorded/Distribution.png", plot = m2a.heavy, width = 350, height = 200, units = "mm", limitsize = F)



## third, NRM
# start with a basic plot
g3.heavy <- ggplot(data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g3.heavy


# now test the relationship using a linear model

hist(data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$MARGIN.SLOPE)

model3.heavy.c <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                     data = data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ])

summary(model3.heavy.c)
drop1(model3.heavy.c, test = "Chi")

r.squaredGLMM(model3.heavy.c)

## no class
model3.heavy <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + VOLTINISM + (1|TAXON),
                       data = data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ])

summary(model3.heavy)
drop1(model3.heavy, test = "Chi")

r.squaredGLMM(model3.heavy)




# add the lines to that figure
newdata3.heavy <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"), VOLTINISM = factor(c(1,2)), CLASS = factor(c("HS","WC")), MARGIN.SLOPE = 0)
newdata3.heavy$MARGIN.SLOPE <- predict(model3.heavy, newdata = newdata3.heavy, type="response")
newdata3.heavy <- ddply(newdata3.heavy, .(PHENO.SLOPE,VOLTINISM,CLASS), numcolwise(mean))



g3.heavy <- ggplot(data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  geom_line(data = newdata3.heavy,
            aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, linetype = VOLTINISM, colour = CLASS),
            stat="identity",size=1)+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  scale_linetype_manual(values = c("dashed","solid"))+
  theme_bw()

g3.heavy



# univoltine, hs
hist(data.uv.hs[which(data.uv.hs$NORTHERLY == F), ]$MARGIN.SLOPE)

model3.u.hs <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.uv.hs[which(data.uv.hs$NORTHERLY == F), ])

summary(model3.u.hs)
drop1(model3.u.hs, test = "Chi")

chkres(model3.u.hs, data.uv.hs[which(data.uv.hs$NORTHERLY == F), ]$PHENO.SLOPE, data.uv.hs[which(data.uv.hs$NORTHERLY == F), ]$TAXON)


# predict data to plot lines
newdata3u.hs <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3u.hs$MARGIN.SLOPE <- predict(model3.u.hs, newdata = newdata3u.hs, type="response")
newdata3u.hs <- ddply(newdata3u.hs, .(PHENO.SLOPE), numcolwise(mean))
newdata3u.hs$CLASS <- "HS"

# univoltine, wc
hist(data.uv.wc[which(data.uv.wc$NORTHERLY == F), ]$MARGIN.SLOPE)

model3.u.wc <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.uv.wc[which(data.uv.wc$NORTHERLY == F), ])

summary(model3.u.wc)
drop1(model3.u.wc, test = "Chi")

chkres(model3.u.wc, data.uv.wc[which(data.uv.wc$NORTHERLY == F), ]$PHENO.SLOPE, data.uv.wc[which(data.uv.wc$NORTHERLY == F), ]$TAXON)


# predict data to plot lines
newdata3u.wc <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3u.wc$MARGIN.SLOPE <- predict(model3.u.wc, newdata = newdata3u.wc, type="response")
newdata3u.wc <- ddply(newdata3u.wc, .(PHENO.SLOPE), numcolwise(mean))
newdata3u.wc$CLASS <- "WC"

newdata3u.all <- rbind(newdata3u.hs,newdata3u.wc)

g3.u.heavy.c <- ggplot(data.uv.heavy[which(data.uv.heavy$NORTHERLY == F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-1,2)+
  ggtitle("Univoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g3.u.heavy.c

g3.u.heavy <- ggplot(data.uv.heavy[which(data.uv.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-1,2)+
  ggtitle("Univoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g3.u.heavy


# multivoltine
# multivoltine, hs
hist(data.mv.hs[which(data.mv.hs$NORTHERLY == F), ]$MARGIN.SLOPE)

# no further analysis since n = 2!
newdata3m.hs <- expand.grid(PHENO.SLOPE = 1, MARGIN.SLOPE = 0, CLASS = "HS")


# multivoltine, wc
hist(data.mv.wc[which(data.mv.wc$NORTHERLY == F), ]$MARGIN.SLOPE)

model3.m.wc <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE + (1|TAXON),
                    data = data.mv.wc[which(data.mv.wc$NORTHERLY == F), ])

summary(model3.m.wc)
drop1(model3.m.wc, test = "Chi")

chkres(model3.m.wc, data.mv.wc[which(data.mv.wc$NORTHERLY == F), ]$PHENO.SLOPE, data.mv.wc[which(data.mv.wc$NORTHERLY == F), ]$TAXON)


# predict data to plot lines
newdata3m.wc <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"),MARGIN.SLOPE = 0)
newdata3m.wc$MARGIN.SLOPE <- predict(model3.m.wc, newdata = newdata3m.wc, type="response")
newdata3m.wc <- ddply(newdata3m.wc, .(PHENO.SLOPE), numcolwise(mean))
newdata3m.wc$CLASS <- "WC"

newdata3m.all <- rbind(newdata3m.hs,newdata3m.wc)

g3.m.heavy.c <- ggplot(data.mv.heavy[which(data.mv.heavy$NORTHERLY == F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata3m.all,
            aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, linetype = CLASS),
            stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  scale_linetype_manual(values = c("blank","dashed"), name = "Class")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g3.m.heavy.c

g3.m.heavy <- ggplot(data.mv.heavy[which(data.mv.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = PHENO.SLOPE, fill = TAXON),shape=24)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-1,2)+
  ggtitle("Multivoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g3.m.heavy


## combine and export the plots
m3.heavy <- grid.arrange(g3.u.heavy,g3.m.heavy,
                         ncol=2)


leg3.heavy <- g_legend(g3.m.heavy)

m3a.heavy <- grid.arrange(arrangeGrob(g3.u.heavy + theme(legend.position="none"),g3.m.heavy + theme(legend.position="none"),
                                      ncol=2),
                          leg3.heavy,nrow=2, heights=c(10,2))

ggsave("Results/Plots/HeavilyRecorded/NRM.png", plot = m3a.heavy, width = 350, height = 200, units = "mm", limitsize = F)




m1all.heavy <- grid.arrange(arrangeGrob(g1.u.heavy + theme(legend.position="none") + xlab("\n"),g1.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                  g2.u.heavy + theme(legend.position="none") + ggtitle("") + xlab("\n"),g2.m.heavy + theme(legend.position="none") + ggtitle("") + xlab("\n") + ylab("\n"),
                                  g3.u.heavy + theme(legend.position="none") + ggtitle(""),g3.m.heavy + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                  ncol=2),
                      leg1.heavy,ncol=2,widths=c(10,2))

m1.heavy.nonabund <- grid.arrange(arrangeGrob(g2.u.heavy + theme(legend.position="none") + xlab("\n"),g2.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                              g3.u.heavy + theme(legend.position="none") + ggtitle(""),g3.m.heavy + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                              ncol=2),
                                  leg1.heavy,ncol=2,widths=c(10,2))

ggsave("Results/Plots/Final/FigD.svg", plot = m1all.heavy, device = svg, width = 250, height = 300, units = "mm", limitsize = F)






### before we move on, let's also test whether Louise Mair's previous observation (that positive abundance trends correlate with positive distribution trends)
# also holds true in our dataset here

#### direct relationship between abundance trend and distribution/NRM trends

## distribution
# start with a basic plot
g4.dist <- ggplot(data.volt.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g4.dist

# now test the relationship using a linear model

model4.dist.c <- lmer(DISTRIB.SLOPE ~ ABUND.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                    data = data.volt.heavy)

summary(model4.dist.c)
drop1(model4.dist.c, test = "Chi")

r.squaredGLMM(model4.dist.c)

## no class

model4.dist <- lmer(DISTRIB.SLOPE ~ ABUND.SLOPE + VOLTINISM + (1|TAXON),
                      data = data.volt.heavy)

summary(model4.dist)
drop1(model4.dist, test = "Chi")

r.squaredGLMM(model4.dist)



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine, hs
hist(data.uv.hs$DISTRIB.SLOPE)

model4.u.hs <- lmer(DISTRIB.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.uv.hs)

summary(model4.u.hs)
drop1(model4.u.hs, test = "Chi")

chkres(model4.u.hs, data.uv.hs$ABUND.SLOPE, data.uv.hs$TAXON)


# predict data to plot lines
newdata4u.hs <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata4u.hs$DISTRIB.SLOPE <- predict(model4.u.hs, newdata = newdata4u.hs, type="response")
newdata4u.hs <- ddply(newdata4u.hs, .(ABUND.SLOPE), numcolwise(mean))
newdata4u.hs$CLASS <- "HS"

# univoltine, wc
hist(data.uv.wc$DISTRIB.SLOPE)

model4.u.wc <- lmer(DISTRIB.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.uv.wc)

summary(model4.u.wc)
drop1(model4.u.wc, test = "Chi")

chkres(model4.u.wc, data.uv.wc$PHENO.SLOPE, data.uv.wc$TAXON)


# predict data to plot lines
newdata4u.wc <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata4u.wc$DISTRIB.SLOPE <- predict(model4.u.wc, newdata = newdata4u.wc, type="response")
newdata4u.wc <- ddply(newdata4u.wc, .(ABUND.SLOPE), numcolwise(mean))
newdata4u.wc$CLASS <- "WC"

newdata4u.all <- rbind(newdata4u.hs,newdata4u.wc)

g4.u.heavy.c <- ggplot(data.uv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata4u.all,
            aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, linetype = CLASS),
            stat="identity",size=1)+
  scale_linetype_manual(values = c("dotted","dashed"), name = "Class")+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-0.1,0.15)+
  ggtitle("Univoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g4.u.heavy.c

# predict data to plot lines
newdata4.dist <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), VOLTINISM = as.factor(c(1,2)), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata4.dist$DISTRIB.SLOPE <- predict(model4.dist, newdata = newdata4.dist, type="response")
newdata4.dist <- ddply(newdata4.dist, .(ABUND.SLOPE), numcolwise(mean))


g4.u.heavy <- ggplot(data.uv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata4.dist,
            aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,15)+
  xlim(-0.1,0.15)+
  ggtitle("Univoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g4.u.heavy

# multivoltine
# multivoltine, hs
hist(data.mv.hs$DISTRIB.SLOPE)

model4.m.hs <- lmer(DISTRIB.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.mv.hs)

summary(model4.m.hs)
drop1(model4.m.hs, test = "Chi")

chkres(model4.m.hs, data.mv.hs$ABUND.SLOPE, data.mv.hs$TAXON)


# predict data to plot lines
newdata4m.hs <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata4m.hs$DISTRIB.SLOPE <- predict(model4.m.hs, newdata = newdata4m.hs, type="response")
newdata4m.hs <- ddply(newdata4m.hs, .(ABUND.SLOPE), numcolwise(mean))
newdata4m.hs$CLASS <- "HS"

# multivoltine, wc
hist(data.mv.wc$DISTRIB.SLOPE)

model4.m.wc <- lmer(DISTRIB.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.mv.wc)

summary(model4.m.wc)
drop1(model4.m.wc, test = "Chi")

chkres(model4.m.wc, data.mv.wc$PHENO.SLOPE, data.mv.wc$TAXON)


# predict data to plot lines
newdata4m.wc <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), TAXON = c("Butterfly","Moth"),DISTRIB.SLOPE = 0)
newdata4m.wc$DISTRIB.SLOPE <- predict(model4.m.wc, newdata = newdata4m.wc, type="response")
newdata4m.wc <- ddply(newdata4m.wc, .(ABUND.SLOPE), numcolwise(mean))
newdata4m.wc$CLASS <- "WC"

newdata4m.all <- rbind(newdata4m.hs,newdata4m.wc)

g4.m.heavy.c <- ggplot(data.mv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata4m.wc,
            aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE),
            linetype = "dashed", stat="identity",size=1)+
  scale_linetype_manual(values = c("dotted","dashed"), name = "Class")+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-0.1,0.15)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g4.m.heavy.c


g4.m.heavy <- ggplot(data.mv.heavy)+
  geom_point(aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata4.dist,
            aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-5,15)+
  xlim(-0.1,0.15)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g4.m.heavy



# combine and export the plots
m4.dist <- grid.arrange(g4.u.heavy, g4.m.heavy,
                        ncol=2)

leg4.dist <- g_legend(g4.u.heavy)

m4a.dist <- grid.arrange(arrangeGrob(g4.u.heavy + theme(legend.position="none"),g4.m.heavy + theme(legend.position="none"),
                                     ncol=2),
                         leg4.dist,nrow=2, heights=c(10,2))


ggsave("Results/Plots/HeavilyRecorded/AbundVSDistribution.png", plot = m4a.dist, width = 350, height = 200, units = "mm", limitsize = F)






## NRM
# start with a basic plot
g4.NRM <- ggplot(data.volt.heavy[which(data.volt.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = ABUND.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g4.NRM

# now test the relationship using a linear model

model4.NRM.c <- lmer(MARGIN.SLOPE ~ VOLTINISM * ABUND.SLOPE * CLASS + (1|TAXON),
                   data = data.volt.heavy[which(data.volt.heavy$NORTHERLY==F), ])

summary(model4.NRM.c)
drop1(model4.NRM.c, test = "Chi")

r.squaredGLMM(model4.NRM.c)


## no class
model4.NRM <- lmer(MARGIN.SLOPE ~ VOLTINISM + ABUND.SLOPE + (1|TAXON),
                     data = data.volt.heavy[which(data.volt.heavy$NORTHERLY==F), ])

summary(model4.NRM)
drop1(model4.NRM, test = "Chi")

r.squaredGLMM(model4.NRM)



## split the analysis up into univoltine and multivoltine, and add the model lines
# univoltine, hs
hist(data.uv.hs[which(data.uv.hs$NORTHERLY==F), ]$MARGIN.SLOPE)

model4.NRM.u.hs <- lmer(MARGIN.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.uv.hs[which(data.uv.hs$NORTHERLY==F), ])

summary(model4.NRM.u.hs)
drop1(model4.NRM.u.hs, test = "Chi")

chkres(model4.NRM.u.hs, data.uv.hs[which(data.uv.hs$NORTHERLY==F), ]$ABUND.SLOPE, data.uv.hs[which(data.uv.hs$NORTHERLY==F), ]$TAXON)



# univoltine, wc
hist(data.uv.wc[which(data.uv.wc$NORTHERLY==F), ]$MARGIN.SLOPE)

model4.NRM.u.wc <- lmer(MARGIN.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.uv.wc[which(data.uv.wc$NORTHERLY==F), ])

summary(model4.NRM.u.wc)
drop1(model4.NRM.u.wc, test = "Chi")

chkres(model4.NRM.u.wc, data.uv.wc[which(data.uv.wc$NORTHERLY==F), ]$ABUND.SLOPE, data.uv.wc[which(data.uv.wc$NORTHERLY==F), ]$TAXON)



# multivoltine
# multivoltine, hs
hist(data.mv.hs[which(data.mv.hs$NORTHERLY==F), ]$MARGIN.SLOPE)


# multivoltine, wc
hist(data.mv.wc[which(data.mv.wc$NORTHERLY==F), ]$MARGIN.SLOPE)

model4.NRM.m.wc <- lmer(MARGIN.SLOPE ~ ABUND.SLOPE + (1|TAXON),
                    data = data.mv.wc[which(data.mv.wc$NORTHERLY==F), ])

summary(model4.NRM.m.wc)
drop1(model4.NRM.m.wc, test = "Chi")

chkres(model4.NRM.m.wc, data.mv.wc[which(data.mv.wc$NORTHERLY==F), ]$ABUND.SLOPE, data.mv.wc[which(data.mv.wc$NORTHERLY==F), ]$TAXON)




# predict data to plot lines
newdata4.nrm.c <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), TAXON = c("Butterfly","Moth"), VOLTINISM = c("1","2"), CLASS = c("HS","WC"), MARGIN.SLOPE = 0)
newdata4.nrm.c$MARGIN.SLOPE <- predict(model4.NRM.c, newdata = newdata4.nrm.c, type="response")
newdata4.nrm.c <- ddply(newdata4.nrm.c, .(ABUND.SLOPE), numcolwise(mean))

# predict data to plot lines
newdata4.nrm <- expand.grid(ABUND.SLOPE = seq(-0.1,0.15,0.01), TAXON = c("Butterfly","Moth"), VOLTINISM = c("1","2"), MARGIN.SLOPE = 0)
newdata4.nrm$MARGIN.SLOPE <- predict(model4.NRM, newdata = newdata4.nrm, type="response")
newdata4.nrm <- ddply(newdata4.nrm, .(ABUND.SLOPE), numcolwise(mean))


# plot figs for univoltine and multivoltine

# univoltine

g4.u.heavy.nrm.c <- ggplot(data.uv.heavy[which(data.uv.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = ABUND.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata4.nrm,
            aes(y = MARGIN.SLOPE, x = ABUND.SLOPE),
            linetype = "solid", stat="identity",size=1)+
  scale_linetype_manual(values = c("dashed","blank"), name = "Class")+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-0.1,0.15)+
  ggtitle("Univoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g4.u.heavy.nrm.c


# multivoltine

g4.m.heavy.nrm.c <- ggplot(data.mv.heavy[which(data.mv.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = ABUND.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = newdata4.nrm,
            aes(y = MARGIN.SLOPE, x = ABUND.SLOPE),
            linetype = "solid", stat="identity",size=1)+
  scale_fill_manual(values = c("Black","White"), name = "Class")+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-0.1,0.15)+
  ggtitle("Multivoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))

g4.m.heavy.nrm.c


g4.u.heavy.nrm <- ggplot(data.uv.heavy[which(data.uv.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = ABUND.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata4.nrm,
            aes(y = MARGIN.SLOPE, x = ABUND.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-0.1,0.15)+
  ggtitle("Univoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g4.u.heavy.nrm

g4.m.heavy.nrm <- ggplot(data.mv.heavy[which(data.mv.heavy$NORTHERLY==F), ])+
  geom_point(aes(y = MARGIN.SLOPE, x = ABUND.SLOPE, fill = TAXON),shape=24)+
  geom_line(data = newdata4.nrm,
            aes(y = MARGIN.SLOPE, x = ABUND.SLOPE),
            colour="black",stat="identity", size = 1)+
  scale_fill_manual(values = c("Black","White"), name = "Taxon")+
  theme_bw()+
  ylim(-10,20)+
  xlim(-0.1,0.15)+
  ggtitle("Multivoltine")+
  ylab("Change in northern range margin \n(km northwards per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))

g4.m.heavy.nrm


# combine and export the plots
m4.NRM <- grid.arrange(g4.u.heavy.nrm, g4.m.heavy.nrm,
                       ncol=2)

leg4.NRM <- g_legend(g4.u.heavy)

m4a.NRM <- grid.arrange(arrangeGrob(g4.u.heavy.nrm + theme(legend.position="none"),g4.m.heavy.nrm + theme(legend.position="none"),
                                    ncol=2),
                        leg4.NRM,nrow=2, heights=c(10,2))


ggsave("Results/Plots/HeavilyRecorded/AbundVSNorthernRangeMargin.png", plot = m4a.NRM, width = 350, height = 200, units = "mm", limitsize = F)


### now we want to create a plot to give a dummy legend, showing linetypes for overall, HS and WC together
data.dummy <- newdata4m.all
data.dummy$CLASS <- ifelse(data.dummy$ABUND.SLOPE < 0, "Overall",data.dummy$CLASS)

data.dummy1 <- data.mv.heavy
data.dummy1$CLASS <- ifelse(data.dummy1$STRICT_VOLTINISM == 0, "Overall",data.dummy1$CLASS)


g4.dummy <- ggplot(data.dummy1)+
  geom_point(aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, fill = CLASS,shape = TAXON))+
  geom_line(data = data.dummy,
            aes(y = DISTRIB.SLOPE, x = ABUND.SLOPE, linetype = CLASS),
            stat="identity",size=1)+
  scale_linetype_manual(values = c("dotted","dashed","solid"), 
                        name = "Class",
                        labels = c("HS","WC","Overall"))+
  scale_fill_manual(values = c("Black","White","Grey"),
                    name = "Class",
                    labels = c("HS","WC","Overall"))+
  scale_shape_manual(values = c(24,25), name = "Taxon")+
  theme_bw()+
  ylim(-5,20)+
  xlim(-0.1,0.15)+
  ggtitle("Multivoltine")+
  ylab("Change in distribution size \n(percentage of recorded hectads occupied per year)")+
  xlab("Change in abundance \n(log(odds ratio) per year)")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(fill = guide_legend(override.aes=list(shape=21)))
  

g4.dummy

leg.dummy <- g_legend(g4.dummy)

## and combine again


m4all <- grid.arrange(arrangeGrob(g4.u.heavy + theme(legend.position="none") + xlab("\n"),g4.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                  g4.u.heavy.nrm + theme(legend.position="none") + ggtitle(""),g4.m.heavy.nrm + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                  ncol=2),
                      leg4.dist,ncol=2,widths=c(10,2))


ggsave("Results/Plots/Final/FigE.svg", plot = m4all, device = svg, width = 250, height = 200, units = "mm", limitsize = F)


# and now combine with the distrib and NRM plots from FigD, 
# in several different trial arrangements (so we can pick the most logical)

m4comb1 <- grid.arrange(arrangeGrob(g2.u.heavy + theme(legend.position="none") + xlab("\n"),g2.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                   g3.u.heavy + theme(legend.position="none") + ggtitle(""),g3.m.heavy + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                   g4.u.heavy + theme(legend.position="none") + ggtitle("") + xlab("\n"),g4.m.heavy + theme(legend.position="none") + ggtitle("") + xlab("\n") + ylab("\n"),
                                   g4.u.heavy.nrm + theme(legend.position="none") + ggtitle(""),g4.m.heavy.nrm + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                   ncol=2),
                        leg4.dist,ncol=2,widths=c(10,2))

ggsave("Results/Plots/Final/FigL1.svg", plot = m4comb1, device = svg, width = 250, height = 400, units = "mm", limitsize = F)



m4comb2 <- grid.arrange(arrangeGrob(g2.u.heavy + theme(legend.position="none") + xlab("\n"),g2.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                    g4.u.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),g4.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                    g3.u.heavy + theme(legend.position="none") + ggtitle(""),g3.m.heavy + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    g4.u.heavy.nrm + theme(legend.position="none") + ggtitle("") + ylab("\n"),g4.m.heavy.nrm + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    ncol=4),
                        leg4.dist,ncol=2,widths=c(10,1))

ggsave("Results/Plots/Final/FigL2.svg", plot = m4comb2, device = svg, width = 500, height = 200, units = "mm", limitsize = F)




m4comb3 <- grid.arrange(arrangeGrob(g2.u.heavy + theme(legend.position="none") + xlab("\n"), g4.u.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                    g2.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"), g4.m.heavy + theme(legend.position="none") + xlab("\n") + ylab("\n"),
                                    g3.u.heavy + theme(legend.position="none") + ggtitle(""), g4.u.heavy.nrm + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    g3.m.heavy + theme(legend.position="none") + ggtitle("") + ylab("\n"), g4.m.heavy.nrm + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    ncol=4),
                        leg4.dist,ncol=2,widths=c(10,1))

ggsave("Results/Plots/Final/FigL3.svg", plot = m4comb3, device = svg, width = 500, height = 200, units = "mm", limitsize = F)


m4comb4 <- grid.arrange(arrangeGrob(g2.u.heavy + theme(panel.grid = element_blank()) + theme(legend.position="none") ,g2.m.heavy + theme(panel.grid = element_blank()) + theme(legend.position="none") + ylab("\n"),
                                    g4.u.heavy + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle(""),g4.m.heavy + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    g3.u.heavy + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle(""),g3.m.heavy + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    g4.u.heavy.nrm + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle(""),g4.m.heavy.nrm + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle("") + ylab("\n"),
                                    ncol=2),
                        leg4.dist,ncol=2,widths=c(10,2))

ggsave("Results/Plots/Final/FigL4.svg", plot = m4comb4, device = svg, width = 250, height = 400, units = "mm", limitsize = F)




### follow-up analyses ####

### we've shown that the same patterns emerge regardless of whether you use recorded, well-recorded or heavily-recorded hectads
# for simplicity, let's do the rest with heavily-recorded only (we have most faith in these)


## given that we find pheno ~ abundance more strongly than pheno ~ distribution or pheno ~ NRM
# and given that previously and here, we see abundance ~ distribution (and possibly abundance ~ NRM)
# we want to test a hierarchical model of distribution ~ (abundance ~ phenology*voltinism)
# we can't specify this directly; instead we need to rephrase that as:
# distribution ~ fitted abundance (from the model abundance ~ phenology*voltinism)

# first let's see how this fitted abundance relates to the actual measured abundance
data.volt.heavy$ABUND.PV.c <- fitted(lmer(ABUND.SLOPE ~ PHENO.SLOPE*VOLTINISM*CLASS + (1|TAXON),
                                        data = data.volt.heavy))

plot(ABUND.SLOPE ~ ABUND.PV.c,
               data = data.volt.heavy)


# there's a correlation there (as you'd hope!) but it doesn't look all that strong

# now try testing the relationship between fitted abundance and distribution
qplot(data.volt.heavy$ABUND.PV.c,data.volt.heavy$DISTRIB.SLOPE,
     colour = data.volt.heavy$TAXON)+
  theme_bw()
  


modelh1.c <- lmer(DISTRIB.SLOPE ~ ABUND.PV.c + (1|TAXON),
               data = data.volt.heavy)


summary(modelh1.c)
drop1(modelh1.c, test = "Chi")

chkres(modelh1.c, data.volt.heavy$ABUND.PV.c, data.volt.heavy$TAXON)

r.squaredGLMM(modelh1.c)



## no class
data.volt.heavy$ABUND.PV <- fitted(lmer(ABUND.SLOPE ~ PHENO.SLOPE*VOLTINISM + (1|TAXON),
                                          data = data.volt.heavy))

plot(ABUND.SLOPE ~ ABUND.PV,
     data = data.volt.heavy)


# there's a correlation there (as you'd hope!) but it doesn't look all that strong

# now try testing the relationship between fitted abundance and distribution
qplot(data.volt.heavy$ABUND.PV,data.volt.heavy$DISTRIB.SLOPE,
      colour = data.volt.heavy$TAXON)+
  theme_bw()



modelh1 <- lmer(DISTRIB.SLOPE ~ ABUND.PV + (1|TAXON),
                  data = data.volt.heavy)


summary(modelh1)
drop1(modelh1, test = "Chi")

chkres(modelh1, data.volt.heavy$ABUND.PV, data.volt.heavy$TAXON)

r.squaredGLMM(modelh1)




# and the relationship between fitted abundance and NRM

qplot(data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$ABUND.PV.c,data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$MARGIN.SLOPE,
      colour = data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$TAXON)+
  theme_bw()

modelh2.c <- lmer(MARGIN.SLOPE ~ ABUND.PV.c + (1|TAXON),
                data = data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ])


summary(modelh2.c)
drop1(modelh2.c, test = "Chi")

chkres(modelh2.c, data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$ABUND.PV.c, data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$TAXON)

r.squaredGLMM(modelh2.c)



## no class
qplot(data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$ABUND.PV,data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$MARGIN.SLOPE,
      colour = data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$TAXON)+
  theme_bw()

modelh2 <- lmer(MARGIN.SLOPE ~ ABUND.PV + (1|TAXON),
               data = data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ])


summary(modelh2)
drop1(modelh2, test = "Chi")

chkres(modelh2, data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$ABUND.PV, data.volt.heavy[which(data.volt.heavy$NORTHERLY == F), ]$TAXON)

r.squaredGLMM(modelh2)



#### intra-specific ####

## we want to see whether these patterns hold when we look, not at the overall trends for each species,
# but rather for the *site-level* trends *within* each species
# i.e. within univoltine species, is abundance trend at each site unrelated to phenology trend,
# and within multivoltine species, do they fare better at sites where their phenology advances more?
# we can further subdivide this to see whether it holds both for multivoltine species that are doing well at species-level
# and also for multivoltine species that are declining at species-level


## we now need to go back to the raw data and, for every combination of site*species, 
# extract a *site-level* trend for each of abundance and phenology
# I've done this in the original data extraction scripts, so I can read this data in here:

intraspec.b <- read.csv("Data/Derived/BMS_intraspec.csv", header = T)
intraspec.m <- read.csv("Data/Derived/RIS_intraspec.csv", header = T)

intraspec <- rbind(intraspec.b,intraspec.m)

## first let's categorize all species into three categories: univoltine, multivoltine increasing, and multivoltine declining

data.volt.heavy$CAT <- ifelse(data.volt.heavy$VOLTINISM==1, "Univoltine",
                              ifelse(data.volt.heavy$ABUND.SLOPE > 0, "Multivoltine increasing", "Multivoltine declining"))

# and bind the species specialisation classes and voltinism into the intraspecific dataset
intraspec <- merge(intraspec, data.volt.heavy[,c(1,4,8)])



## now split the dataset by the same categories
# species-level
spec.uv <- data.volt.heavy[which(data.volt.heavy$VOLTINISM == 1), ]
spec.mv <- data.volt.heavy[which(data.volt.heavy$VOLTINISM == 2), ]

spec.mv.pos <- data.volt.heavy[which(data.volt.heavy$VOLTINISM == 2 & data.volt.heavy$ABUND.SLOPE > 0), ]
spec.mv.neg <- data.volt.heavy[which(data.volt.heavy$VOLTINISM == 2 & data.volt.heavy$ABUND.SLOPE < 0), ]

spec.uv.hs <- spec.uv[which(spec.uv$CLASS == "HS"), ]
spec.uv.wc <- spec.uv[which(spec.uv$CLASS == "WC"), ]
spec.mv.hs <- spec.mv[which(spec.mv$CLASS == "HS"), ]
spec.mv.wc <- spec.mv[which(spec.mv$CLASS == "WC"), ]
spec.mv.pos.hs <- spec.mv.pos[which(spec.mv.pos$CLASS == "HS"), ]
spec.mv.pos.wc <- spec.mv.pos[which(spec.mv.pos$CLASS == "WC"), ]
spec.mv.neg.hs <- spec.mv.neg[which(spec.mv.neg$CLASS == "HS"), ]
spec.mv.neg.wc <- spec.mv.neg[which(spec.mv.neg$CLASS == "WC"), ]


# transect-level
intra.uv <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.uv$COMMON_NAME))), ]
intra.mv <- intraspec[which(intraspec$COMMON_NAME %in% c(levels(droplevels(spec.mv.pos$COMMON_NAME)), levels(droplevels(spec.mv.neg$COMMON_NAME)))), ]
intra.mv.pos <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.pos$COMMON_NAME))), ]
intra.mv.neg <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.neg$COMMON_NAME))), ]

intra.uv.hs <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.uv.hs$COMMON_NAME))), ]
intra.uv.wc <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.uv.wc$COMMON_NAME))), ]
intra.mv.hs <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.hs$COMMON_NAME))), ]
intra.mv.wc <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.wc$COMMON_NAME))), ]
intra.mv.pos.hs <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.pos.hs$COMMON_NAME))), ]
intra.mv.pos.wc <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.pos.wc$COMMON_NAME))), ]
intra.mv.neg.hs <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.neg.hs$COMMON_NAME))), ]
intra.mv.neg.wc <- intraspec[which(intraspec$COMMON_NAME %in% levels(droplevels(spec.mv.neg.wc$COMMON_NAME))), ]


## overarching test, using species as a random effect
hist(intraspec$ABUND.SLOPE)

modeli.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|COMMON_NAME),
               data = intraspec)

summary(modeli.c)
drop1(modeli.c, test = "Chi")

r.squaredGLMM(modeli.c)

## no class
modeli <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|COMMON_NAME),
                 data = intraspec)

summary(modeli)
drop1(modeli, test = "Chi")

r.squaredGLMM(modeli)


# now look within each category at the relationship between phenology and abundance, using species as a random effect

## first for univoltine, hs

modeli.u.hs <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                   data = intra.uv.hs)

summary(modeli.u.hs)
drop1(modeli.u.hs, test = "Chi")

r.squaredGLMM(modeli.u.hs)

newdataiu.hs <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.uv.hs$COMMON_NAME)), ABUND.SLOPE = 0)
newdataiu.hs$ABUND.SLOPE <- predict(modeli.u.hs, newdata = newdataiu.hs, type="response")
newdataiu.hs <- ddply(newdataiu.hs, .(PHENO.SLOPE), numcolwise(median))


## for univoltine, wc

modeli.u.wc <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                    data = intra.uv.wc)

summary(modeli.u.wc)
drop1(modeli.u.wc, test = "Chi")

r.squaredGLMM(modeli.u.wc)

newdataiu.wc <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.uv.wc$COMMON_NAME)), ABUND.SLOPE = 0)
newdataiu.wc$ABUND.SLOPE <- predict(modeli.u.wc, newdata = newdataiu.wc, type="response")
newdataiu.wc <- ddply(newdataiu.wc, .(PHENO.SLOPE), numcolwise(median))

## for multivoltine, hs

modeli.m.hs <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                    data = intra.mv.hs)

summary(modeli.m.hs)
drop1(modeli.m.hs, test = "Chi")

r.squaredGLMM(modeli.m.hs)

newdataim.hs <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.mv.hs$COMMON_NAME)), ABUND.SLOPE = 0)
newdataim.hs$ABUND.SLOPE <- predict(modeli.m.hs, newdata = newdataim.hs, type="response")
newdataim.hs <- ddply(newdataim.hs, .(PHENO.SLOPE), numcolwise(median))


## for multivoltine, wc

modeli.m.wc <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                    data = intra.mv.wc)

summary(modeli.m.wc)
drop1(modeli.m.wc, test = "Chi")

r.squaredGLMM(modeli.m.wc)

newdataiu.wc <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.uv.wc$COMMON_NAME)), ABUND.SLOPE = 0)
newdataiu.wc$ABUND.SLOPE <- predict(modeli.u.wc, newdata = newdataiu.wc, type="response")
newdataiu.wc <- ddply(newdataiu.wc, .(PHENO.SLOPE), numcolwise(median))





# as well as this overall line, we also want an individual line for each species
# because species is a random effect in the model this can't come from the model directly

newdata.specs <- data.frame(PHENO.SLOPE = numeric(),
                            COMMON_NAME = factor(),
                            ABUND.SLOPE = numeric(),
                            N = numeric(),
                            Significance = factor(),
                            Relationship = numeric(),
                            State = factor(),
                            Class = factor(),
                            Voltinism = factor())

for (x in levels(droplevels(intraspec$COMMON_NAME))){
  spec <- intraspec[which(intraspec$COMMON_NAME == x), ]
  n <- nrow(spec)
 
  model.spec <- lm(ABUND.SLOPE ~ PHENO.SLOPE,
                   data = spec)
  
  model.sig <- ifelse(drop1(model.spec, test = "F")[2,6] < 0.05, "Significant","Non-significant")
  model.slope <- summary(model.spec)$coefficients[2,1]
  model.state <- ifelse(model.sig == "Non-significant","Non-significant",
                        ifelse(model.slope > 0, "Positive","Negative"))
  
  
  newdataspec <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = x, ABUND.SLOPE = 0)
  newdataspec$ABUND.SLOPE <- predict(model.spec, newdata = newdataspec, type="response")
  newdataspec$N <- n
  newdataspec$Significance <- model.sig
  newdataspec$Relationship <- model.slope
  newdataspec$State <- model.state
  newdataspec$Class <- as.character(spec$CLASS[[1]])
  newdataspec$Voltinism <- as.character(spec$VOLTINISM[[1]])
  
  
  
  newdata.specs <- rbind(newdata.specs,newdataspec)
}

newdata.specs.uv <- newdata.specs[which(newdata.specs$Voltinism == 1), ]
newdata.specs.mv <- newdata.specs[which(newdata.specs$Voltinism == 2), ]
newdata.specs.mvpos <- newdata.specs.mv[which(newdata.specs.mv$COMMON_NAME %in% levels(droplevels(spec.mv.pos$COMMON_NAME))), ]
newdata.specs.mvneg <- newdata.specs.mv[which(newdata.specs.mv$COMMON_NAME %in% levels(droplevels(spec.mv.neg$COMMON_NAME))), ]


## plot univoltines
# first get an overall line
modeli.uv.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * CLASS + (1|COMMON_NAME),
                  data = intra.uv)

summary(modeli.uv.c)
drop1(modeli.uv.c, test = "Chi")

r.squaredGLMM(modeli.uv.c)



modeli.uv <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                  data = intra.uv)

summary(modeli.uv)
drop1(modeli.uv, test = "Chi")

r.squaredGLMM(modeli.uv)

# no relationship = no line = no need to go further

# then plot

gi.u.c <- ggplot(intra.uv)+
  geom_line(data = newdata.specs.uv,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, group = COMMON_NAME, linetype = Class, colour = Significance),
            stat="identity",size = 0.4)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE),shape=24)+
  theme_bw()+
  scale_linetype_manual(values = c("dotted","dashed"))+
  scale_colour_manual(values = c("grey90","black"))+
  ylim(-0.4,0.4)+
  xlim(-5,5)+
  ggtitle("Univoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

gi.u.c

gi.u <- ggplot(intra.uv)+
  geom_line(data = newdata.specs.uv,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, group = COMMON_NAME, linetype = Significance, colour = Significance),
            stat="identity",size = 0.7)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE),shape=24)+
  theme_bw()+
  scale_linetype_manual(values = c("dotted","solid"))+
  scale_colour_manual(values = c("grey80","grey70"))+
  ylim(-0.4,0.4)+
  xlim(-5,5)+
  ggtitle("Univoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

gi.u


species.intra.uv <- ddply(newdata.specs.uv, .(COMMON_NAME, Significance, Relationship, State, N, Class, Voltinism), numcolwise(mean))
species.intra.uv$Significance <- as.factor(species.intra.uv$Significance)
species.intra.uv$State <- as.factor(species.intra.uv$State)
species.intra.uv$Class <- as.factor(species.intra.uv$Class)
summary(species.intra.uv)



## plot multivoltines
# first get an overall line

modeli.mv.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * CLASS + (1|COMMON_NAME),
                  data = intra.mv)

summary(modeli.mv.c)
drop1(modeli.mv.c, test = "Chi")

r.squaredGLMM(modeli.mv.c)


## no class
modeli.mv <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                  data = intra.mv)

summary(modeli.mv)
drop1(modeli.mv, test = "Chi")

r.squaredGLMM(modeli.mv)

newdataimv.c <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.mv$COMMON_NAME)), CLASS = c("HS","WC"), ABUND.SLOPE = 0)
newdataimv.c$ABUND.SLOPE <- predict(modeli.mv.c, newdata = newdataimv.c, type="response")
newdataimv.c <- ddply(newdataimv.c, .(PHENO.SLOPE), numcolwise(median))


newdataimv <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.mv$COMMON_NAME)), ABUND.SLOPE = 0)
newdataimv$ABUND.SLOPE <- predict(modeli.mv, newdata = newdataimv, type="response")
newdataimv <- ddply(newdataimv, .(PHENO.SLOPE), numcolwise(median))


# plot
gi.mv.c <- ggplot(intra.mv)+
  geom_line(data = newdata.specs.mv,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, group = COMMON_NAME, linetype = Class, colour = Significance),
            stat="identity",size = 0.4)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE),shape=24)+
  geom_line(data = newdataimv.c,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            stat="identity", size = 1.5, linetype = "solid", colour = "royalblue")+
  theme_bw()+
  scale_linetype_manual(values = c("dotted","dashed"))+
  scale_colour_manual(values = c("grey90","black"))+
  ylim(-0.4,0.4)+
  xlim(-5,5)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

gi.mv.c


gi.mv <- ggplot(intra.mv)+
  geom_line(data = newdata.specs.mv,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, group = COMMON_NAME, linetype = Significance, colour = Significance),
            stat="identity",size = 0.7)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE),shape=24)+
  geom_line(data = newdataimv,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype = "solid")+
  theme_bw()+
  scale_linetype_manual(values = c("dotted","solid"))+
  scale_colour_manual(values = c("grey80","grey70"))+
  ylim(-0.4,0.4)+
  xlim(-5,5)+
  ggtitle("Multivoltine")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

gi.mv


species.intra.mv <- ddply(newdata.specs.mv, .(COMMON_NAME, Significance, Relationship, State, N), numcolwise(mean))
species.intra.mv$Significance <- as.factor(species.intra.mv$Significance)

species.intra.mv$State <- as.factor(species.intra.mv$State)
summary(species.intra.mv)


## then multivoltine increasing - we'll no longer bother with class since it was n.s. above
# first get an overall line

modeli.mvpos <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                    data = intra.mv.pos)

summary(modeli.mvpos)
drop1(modeli.mvpos, test = "Chi")

r.squaredGLMM(modeli.mvpos)

newdataimvpos <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.mv.pos$COMMON_NAME)), ABUND.SLOPE = 0)
newdataimvpos$ABUND.SLOPE <- predict(modeli.mvpos, newdata = newdataimvpos, type="response")
newdataimvpos <- ddply(newdataimvpos, .(PHENO.SLOPE), numcolwise(median))


# plot
gi.mvpos <- ggplot(intra.mv.pos)+
  geom_line(data = newdata.specs.mvpos,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, group = COMMON_NAME, linetype = Significance, colour = Significance),
            stat="identity",size = 0.7)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE),shape=24)+
  geom_line(data = newdataimvpos,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype = "solid")+
  theme_bw()+
  scale_linetype_manual(values = c("dotted","solid"))+
  scale_colour_manual(values = c("grey80","grey70"))+
  ylim(-0.4,0.4)+
  xlim(-5,5)+
  ggtitle("Multivoltine, increasing abundance")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

gi.mvpos




## and multivoltine declining

modeli.mvneg <- lmer(ABUND.SLOPE ~ PHENO.SLOPE + (1|COMMON_NAME),
                   data = intra.mv.neg)

summary(modeli.mvneg)
drop1(modeli.mvneg, test = "Chi")

r.squaredGLMM(modeli.mvneg)

newdataimvneg <- expand.grid(PHENO.SLOPE = seq(-5,5,0.1), COMMON_NAME = levels(droplevels(intra.mv.neg$COMMON_NAME)), ABUND.SLOPE = 0)
newdataimvneg$ABUND.SLOPE <- predict(modeli.mvneg, newdata = newdataimvneg, type="response")
newdataimvneg <- ddply(newdataimvneg, .(PHENO.SLOPE), numcolwise(median))


# plot
gi.mvneg <- ggplot(intra.mv.neg)+
  geom_line(data = newdata.specs.mvneg,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, group = COMMON_NAME, linetype = Significance, colour = Significance),
            stat="identity",size = 0.7)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE),shape=24)+
  geom_line(data = newdataimvneg,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE),
            colour="black",stat="identity", size = 1, linetype = "solid")+
  theme_bw()+
  scale_linetype_manual(values = c("dotted","solid"))+
  scale_colour_manual(values = c("grey80","grey70"))+
  ylim(-0.4,0.4)+
  xlim(-5,5)+
  ggtitle("Multivoltine, declining abundance")+
  ylab("Change in abundance \n(log(odds ratio) per year)")+
  xlab("Advance in emergence date \n(days earlier per year)")+
  theme(plot.title = element_text(hjust = 0.5))

gi.mvneg




## combine and export the plots
mi <- grid.arrange(gi.mvpos, gi.mvneg,
                         ncol=2)


gi.ua <- gi.u + theme(legend.position = "bottom")

leg.i <- g_legend(gi.u)
leg.ia <- g_legend(gi.ua)
  
mi.a <- grid.arrange(arrangeGrob(gi.mvpos + theme(panel.grid = element_blank()) + theme(legend.position="none") + ggtitle(" "),gi.mvneg  + theme(panel.grid = element_blank()) + theme(legend.position="none") + ylab(" ") + ggtitle(" "),
                                       ncol=2),
                           leg.ia,nrow=2, heights=c(10,2))



ggsave("Results/Plots/Final/FigF.svg", plot = mi.a, device = svg, width = 250, height = 125, units = "mm", limitsize = F)





## then another version of this figure without multivoltine separated out, but with univoltine included

## combine and export the plots
mi.res <- grid.arrange(gi.u, gi.mv,
                   ncol=2)


mi.res.a <- grid.arrange(arrangeGrob(gi.u + theme(legend.position="none"),gi.mv + theme(legend.position="none") + ylab(" "),
                                 ncol=2),
                     leg.ia,nrow=2, heights=c(10,2))



ggsave("Results/Plots/Final/FigG.svg", plot = mi.res.a, device = svg, width = 250, height = 150, units = "mm", limitsize = F)




#### phylo constraint ####
## now we want to repeat the main phenology vs abundance analysis, accounting for phylogenetic autocorrelation

## read in a tree (built at COI locus in Geneious)

phy <- read.tree(file = "Phylogeny/Lepi_C alignment 2 tree.newick")
plot(phy)

# attach tree to the full dataset
cdat <- comparative.data(data=data.volt.heavy, phy=phy, names.col="TREE_NAME", na.omit=F)
print(cdat)

# split into uni and multivoltine
cdat.u <- subset(cdat, subset= VOLTINISM==1)
print(cdat.u)

cdat.m <- subset(cdat, subset= VOLTINISM==2)
print(cdat.m)


## fit and test the PGLS models - first for abundance

PGLS.ab.c <- pgls(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS, cdat)
summary(PGLS.ab.c)
anova(PGLS.ab.c)
AIC(PGLS.ab.c)

PGLS.ab <- pgls(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM, cdat)
summary(PGLS.ab)
anova(PGLS.ab)
AIC(PGLS.ab)



## next, distribution

PGLS.d.c <- pgls(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS, cdat)
summary(PGLS.d.c)
anova(PGLS.d.c)
AIC(PGLS.d.c)


PGLS.d <- pgls(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM, cdat)
summary(PGLS.d)
anova(PGLS.d)
AIC(PGLS.d)

## now, northern range margin - which needs a further subsetting
cdat.n <- subset(cdat, subset= NORTHERLY==F)
print(cdat.n)

PGLS.nrm.c <- pgls(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS, cdat.n)
summary(PGLS.nrm.c)
anova(PGLS.nrm.c)
AIC(PGLS.nrm.c)


PGLS.nrm <- pgls(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM, cdat.n)
summary(PGLS.nrm)
anova(PGLS.nrm)
AIC(PGLS.nrm)


## next, distribution ~ abundance
PGLS.da.c <- pgls(DISTRIB.SLOPE ~ ABUND.SLOPE * VOLTINISM * CLASS, cdat)
summary(PGLS.da.c)
anova(PGLS.da.c)
AIC(PGLS.da.c)

PGLS.da <- pgls(DISTRIB.SLOPE ~ ABUND.SLOPE, cdat)
summary(PGLS.da)
anova(PGLS.da)
AIC(PGLS.da)


PGLS.nrm.c <- pgls(MARGIN.SLOPE ~ ABUND.SLOPE * VOLTINISM * CLASS, cdat.n)
summary(PGLS.nrm.c)
anova(PGLS.nrm.c)
AIC(PGLS.nrm.c)

PGLS.nrm <- pgls(MARGIN.SLOPE ~ ABUND.SLOPE, cdat.n)
summary(PGLS.nrm)
anova(PGLS.nrm)
AIC(PGLS.nrm)


## next, distribution ~ modelled abundance
cdat$data$FITTED.c <- fitted(PGLS.ab.c)

PGLS.h.c <- lm(DISTRIB.SLOPE ~ FITTED.c, cdat$data)
summary(PGLS.h.c)
drop1(PGLS.h.c, test = "F")


cdat$data$FITTED <- fitted(PGLS.ab)

PGLS.h <- lm(DISTRIB.SLOPE ~ FITTED, cdat$data)
summary(PGLS.h)
drop1(PGLS.h, test = "F")


cdat.n <- subset(cdat, subset= NORTHERLY==F)

PGLS.hnrm.c <- lm(MARGIN.SLOPE ~ FITTED.c, cdat.n$data)
summary(PGLS.hnrm.c)
drop1(PGLS.hnrm.c, test = "F")


PGLS.hnrm <- lm(MARGIN.SLOPE ~ FITTED, cdat.n$data)
summary(PGLS.hnrm)
drop1(PGLS.hnrm, test = "F")



#### strict voltinism categories
## now we want to repeat the main analyses, 
# but using only species where we can be certain that their voltinism category is the same across all sites
# i.e. they are never univoltine, even in the far north; or they are never multivoltine, even in the far south

# create a sub-dataset using only the species for which we're certain of their voltinism
data.volt.heavy.r <- data.volt.heavy[which(data.volt.heavy$STRICT_VOLTINISM == 1), ]
summary(data.volt.heavy.r) # we lose about half of multivoltines, and a few univoltines

data.uv.heavy.r <- data.uv.heavy[which(data.uv.heavy$STRICT_VOLTINISM == 1), ]
summary(data.uv.heavy.r)

data.mv.heavy.r <- data.mv.heavy[which(data.mv.heavy$STRICT_VOLTINISM == 1), ]
summary(data.mv.heavy.r)

## first, abundance
# start with a basic plot
g1.heavy.r <- ggplot(data.volt.heavy.r)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g1.heavy.r


# now test the relationship using a linear model

hist(data.volt.heavy.r$ABUND.SLOPE)

model1.heavy.r.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                     data = data.volt.heavy.r)

summary(model1.heavy.r.c)
drop1(model1.heavy.r.c, test = "Chi")

r.squaredGLMM(model1.heavy.r.c)


model1.heavy.r <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                         data = data.volt.heavy.r)

summary(model1.heavy.r)
drop1(model1.heavy.r, test = "Chi")

r.squaredGLMM(model1.heavy.r)



## out of interest, are trends different for the flexible species in each category?
# create a sub-dataset using only the species for which we're certain of their voltinism
data.volt.heavy.f <- data.volt.heavy[which(data.volt.heavy$STRICT_VOLTINISM == 0), ]
summary(data.volt.heavy.f) # we lose about half of multivoltines, and a lots of univoltines

data.uv.heavy.f <- data.uv.heavy[which(data.uv.heavy$STRICT_VOLTINISM == 0), ]
summary(data.uv.heavy.f)

data.mv.heavy.f <- data.mv.heavy[which(data.mv.heavy$STRICT_VOLTINISM == 0), ]
summary(data.mv.heavy.f)

## first, abundance
# start with a basic plot
g1.heavy.f <- ggplot(data.volt.heavy.f)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g1.heavy.f


# now test the relationship using a linear model

hist(data.volt.heavy.f$ABUND.SLOPE)

model1.heavy.f.c <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                       data = data.volt.heavy.f)

summary(model1.heavy.f.c)
drop1(model1.heavy.f.c, test = "Chi")

r.squaredGLMM(model1.heavy.f.c)


model1.heavy.f <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                       data = data.volt.heavy.f)

summary(model1.heavy.f)
drop1(model1.heavy.f, test = "Chi")

r.squaredGLMM(model1.heavy.f)

## second, distribution
# start with a basic plot
g2.heavy.r <- ggplot(data.volt.heavy.r)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g2.heavy.r

# now test the relationship using a linear model

model2.heavy.r.c <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                     data = data.volt.heavy.r)

summary(model2.heavy.r.c)
drop1(model2.heavy.r.c, test = "Chi")

r.squaredGLMM(model2.heavy.r.c)

model2.heavy.r <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                         data = data.volt.heavy.r)

summary(model2.heavy.r)
drop1(model2.heavy.r, test = "Chi")

r.squaredGLMM(model2.heavy.r)


## flexibles
# start with a basic plot
g2.heavy.f <- ggplot(data.volt.heavy.f)+
  geom_point(aes(y = DISTRIB.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g2.heavy.f

# now test the relationship using a linear model

model2.heavy.f.c <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                       data = data.volt.heavy.f)

summary(model2.heavy.f.c)
drop1(model2.heavy.f.c, test = "Chi")

r.squaredGLMM(model2.heavy.f.c)


model2.heavy.f <- lmer(DISTRIB.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                       data = data.volt.heavy.f)

summary(model2.heavy.f)
drop1(model2.heavy.f, test = "Chi")

r.squaredGLMM(model2.heavy.f)


## third, NRM
# start with a basic plot
g3.heavy.r <- ggplot(data.volt.heavy.r[which(data.volt.heavy.r$NORTHERLY==F), ])+
  geom_point(aes(y = PHENO.SLOPE, x = MARGIN.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g3.heavy.r

# now test the relationship using a linear model

model3.heavy.r.c <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM * CLASS + (1|TAXON),
                     data = data.volt.heavy.r[which(data.volt.heavy.r$NORTHERLY==F), ])

summary(model3.heavy.r.c)
drop1(model3.heavy.r.c, test = "Chi")

r.squaredGLMM(model3.heavy.r.c)


model3.heavy.r <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                         data = data.volt.heavy.r[which(data.volt.heavy.r$NORTHERLY==F), ])

summary(model3.heavy.r)
drop1(model3.heavy.r, test = "Chi")

r.squaredGLMM(model3.heavy.r)


## flexibles
# start with a basic plot
g3.heavy.f <- ggplot(data.volt.heavy.f[which(data.volt.heavy.f$NORTHERLY==F), ])+
  geom_point(aes(y = PHENO.SLOPE, x = MARGIN.SLOPE, shape = TAXON, fill = VOLTINISM))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g3.heavy.f

# now test the relationship using a linear model

model3.heavy.f.c <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM *CLASS + (1|TAXON),
                       data = data.volt.heavy.f[which(data.volt.heavy.f$NORTHERLY==F), ])

summary(model3.heavy.f.c)
drop1(model3.heavy.f.c, test = "Chi")

r.squaredGLMM(model3.heavy.f.c)


model3.heavy.f <- lmer(MARGIN.SLOPE ~ PHENO.SLOPE * VOLTINISM + (1|TAXON),
                       data = data.volt.heavy.f[which(data.volt.heavy.f$NORTHERLY==F), ])

summary(model3.heavy.f)
drop1(model3.heavy.f, test = "Chi")

r.squaredGLMM(model3.heavy.f)



#### early- vs late-emergers ###

# as a check, we want to confirm that these patterns cannot be explained simply by 
# whether species' first generation is early or late in the summer

# I might refine this later, but for now I've simply defined species categorically
# by whether the peak of their first generation is before or after 21st June in an average year ('TIMING')

# so, take the basic analyses above and re-run with TIMING replacing VOLTINISM
# first inspect the correlation between the two...

plot(VOLTINISM ~ TIMING, data.volt.heavy)

# about half of early-emerging species are univoltine, and almost all of late-emerging
# so there is a correlation, but also enough differentiation that we should be able to distinguish the effects


## first, abundance
# start with a basic plot
g1.heavy.el <- ggplot(data.volt.heavy)+
  geom_point(aes(y = ABUND.EXP.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = TIMING))+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  theme_bw()

g1.heavy.el


# now test the relationship using a linear model

hist(data.volt.heavy$ABUND.SLOPE)

model1.heavy.el <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * TIMING * CLASS + (1|TAXON),
                     data = data.volt.heavy)

summary(model1.heavy.el)
drop1(model1.heavy.el, test = "Chi")


# add the lines to that figure
newdata1.heavy.el <- expand.grid(PHENO.SLOPE = seq(-3,2,0.1), TAXON = c("Butterfly","Moth"), TIMING = factor(c("Early","Late")),ABUND.SLOPE = 0)
newdata1.heavy.el$ABUND.SLOPE <- predict(model1.heavy.el, newdata = newdata1.heavy.el, type="response")
newdata1.heavy.el <- ddply(newdata1.heavy.el, .(PHENO.SLOPE,TIMING), numcolwise(mean))



g1.heavy.el <- ggplot(data.volt.heavy)+
  geom_point(aes(y = ABUND.SLOPE, x = PHENO.SLOPE, shape = TAXON, fill = TIMING))+
  geom_line(data = newdata1.heavy.el,
            aes(y = ABUND.SLOPE, x = PHENO.SLOPE, linetype = TIMING),
            colour="black",stat="identity",size=1)+
  scale_shape_manual(values = c(24,25))+
  scale_fill_manual(values = c("white","black"))+
  scale_linetype_manual(values = c("dashed","solid"))+
  theme_bw()

g1.heavy.el


# so it does appear that the early/late timing split can explain the pattern...!
# but does it explain it better or worse than the univoltine/multivoltine split?


model1.heavy.umel <- lmer(ABUND.SLOPE ~ PHENO.SLOPE * TIMING * CLASS + (1|TAXON),
                      data = data.volt.heavy) 

summary(model1.heavy.umel)
drop1(model1.heavy.umel, test = "Chi")

# this suggests that the PHENO.SLOPE:VOLTINISM:CLASS interaction is the better predictor so we can move on

# does TIMING help explain any variation within univoltines or multivoltines?
summary(data.uv.heavy)

## split the analysis up into early and late, and add the model lines
# early
hist(data.uv.heavy$ABUND.SLOPE)

model1.u.heavy <- lmer(ABUND.SLOPE ~ PHENO.SLOPE*TIMING + (1|TAXON),
                       data = data.uv.heavy)

summary(model1.u.heavy)
drop1(model1.u.heavy, test = "Chi")

chkres(model1.u.heavy, data.uv.heavy$PHENO.SLOPE, data.uv.heavy$TAXON)


# multivoltine
hist(data.mv.heavy$ABUND.SLOPE)

model1.m.heavy <- lmer(ABUND.SLOPE ~ PHENO.SLOPE*TIMING + (1|TAXON),
                       data = data.mv.heavy)

summary(model1.m.heavy)
drop1(model1.m.heavy, test = "Chi")

chkres(model1.m.heavy, data.mv.heavy$PHENO.SLOPE, data.mv.heavy$TAXON)

# ...no!



#### sites map figure ####
# let's take this opportunity to map out where all the sites fall and how many species each has of univoltine and multivoltine

# read in details of where the sites are and how many univoltine and multivoltine species they had
sites_BMS <- read.csv("Data/UKBMS_sitecounts.csv", header = T)
sites_RIS <- read.csv("Data/RIS_sitecounts.csv", header = T)


# first just try plotting sites as points on a baseR map
blighty()
points(sites_BMS$HEast, sites_BMS$HNorth)

blighty()
points(sites_RIS$HEast, sites_RIS$HNorth)


# now try building a jazzier map

sites_BMS.vec <- levels(droplevels(sites_BMS$SITENAME))
sites_RIS.vec <- levels(droplevels(sites_RIS$SITE))


svg("Results/Plots/sites.svg", width = 12, height = 12, pointsize = 15, bg = "white")
par( mfrow = c(1,2), oma = c(0,0,2,0))

blighty()

for (x in sites_BMS.vec){
  site <- sites_BMS[which(sites_BMS$SITENAME == x), ]
  floating.pie(xpos = site$HEast[[1]], ypos = site$HNorth[[1]], 
               x = c(site$UNIVOLTINES[[1]], site$MULTIVOLTINES[[1]]),
               radius = 3*sqrt(site$SPECIES[[1]]), startpos = pi/2,
               col = c("goldenrod","royalblue"))
  title(main = "UKBMS transects")
  legend("topleft", legend = c("Univoltine species","Multivoltine species"), fill = c("goldenrod","royalblue"))

}

blighty()

for (x in sites_RIS.vec){
  site <- sites_RIS[which(sites_RIS$SITE == x), ]
  floating.pie(xpos = site$HEast[[1]], ypos = site$HNorth[[1]], 
               x = c(site$UNIVOLTINES[[1]], site$MULTIVOLTINES[[1]]),
               radius = 3*sqrt(site$SPECIES[[1]]), startpos = pi/2,
               col = c("goldenrod","royalblue"))
  title(main = "RIS traps")
  
}

dev.off()







#### step-by-step explanation ####
## now we want to break down our explanation of these patterns and see whether we can demonstrate each step:
# 1. early emergence correlates with hot spring (analysis of GDD up to a fixed date/s)
# 2. in early years, there is a higher proportion of individuals after the first generation peaks (higher abundance 2nd gen and/or additional gens)
# 3. this leads to a stronger performance in the following year
# 4. and subsequently to a northwards advance in the following year

### spring temperatures
# the first step in showing all this is to get hold of GDD5 for all relevant sites up to a set date
# we need a dataframe containing site, grid reference (plus X and Y coordinates), year and Julian day of interest

siteGDD_BMS <- read.csv("Data/UKBMS_sites_byseason_gddJ_5km.csv", header = T)
summary(siteGDD_BMS)


siteGDD_RIS <- read.csv("Data/RIS_sites_byseason_gddJ_5km.csv", header = T)
summary(siteGDD_RIS)

# there are a few different columns in each:
colnames(siteGDD_BMS)
colnames(siteGDD_RIS)

# so rearrange them to match first column order and content, and then colnames
siteGDD_BMS <- siteGDD_BMS[,-6]
siteGDD_RIS <- siteGDD_RIS[,c(1:5,12:15,6:7,16:19)]

colnames(siteGDD_BMS)
colnames(siteGDD_RIS)


colnames(siteGDD_RIS) <- colnames(siteGDD_BMS)

# now merge the two together
site_GDD <- rbind(siteGDD_BMS, siteGDD_RIS)

summary(site_GDD)

# finally we want to cast the GDD estimates for the three julian days out into separate columns (this will help with merging this data into abundance/pheno data)

# we'll need to simplify the julian day column because at the moment it contains adjustments for leap years which will prevent us from casting properly

site_GDD$JULIAN_DAY <- ifelse(site_GDD$JULIAN_DAY %in% c(60,61), 60,
                              ifelse(site_GDD$JULIAN_DAY %in% c(152,153), 152, 244))

site_GDD_wide <- dcast(site_GDD, SITENAME + YEAR + SITENO + EASTING + NORTHING + HEast + HNorth + LONGITUDE + LATITUDE + x + y
                       ~ JULIAN_DAY, value.var = 'gdd5J')


colnames(site_GDD_wide)[c(1,12:14)] <- c("SITE","SPRING","SUMMER","AUTUMN")

# generate a variable for spring temp only (i.e. different to SUMMER, which is currently "1st January up to start of summer")

site_GDD_wide$SPRING.ONLY <- site_GDD_wide$SUMMER - site_GDD_wide$SPRING

# now we need the annual data on the individual populations
# read this in as well

popns_BMS <- read.csv("Data/UKBMS_annual.csv", header = T)
summary(popns_BMS)

popns_RIS <- read.csv("Data/RIS_annual.csv", header = T)
summary(popns_RIS)

# make the columns the same in both
popns_RIS <- popns_RIS[,c(1:6,8,11:21)]


# merge each sheet together
popns_all <- rbind(popns_BMS, popns_RIS)
popns_all_volti <- merge(popns_all, voltinism)
popns_GDD <- merge(popns_all_volti,site_GDD_wide)

summary(popns_GDD)

popns_GDD <- popns_GDD[!is.na(popns_GDD$SPRING.ONLY), ]
popns_GDD$VOLTINISM <- as.factor(popns_GDD$VOLTINISM)

summary(popns_GDD)

### now we can do the analysis of spring temp vs phenology

# we want an overall relationship and ones broken down by categories

hist(popns_GDD$SPRING.ONLY)

modelGDD.spr.c <- lmer(PEAKDAY ~ SPRING.ONLY * VOLTINISM * CLASS + (1|COMMON_NAME) + (1|YEAR),
                     data = popns_GDD)

summary(modelGDD.spr.c)
drop1(modelGDD.spr.c, test = "Chi")

chkres(modelGDD.spr.c, popns_GDD$PEAKDAY, popns_GDD$CLASS)

r.squaredGLMM(modelGDD.spr.c)


modelGDD.spr <- lmer(PEAKDAY ~ SPRING.ONLY * VOLTINISM + (1|COMMON_NAME) + (1|YEAR),
                     data = popns_GDD)

summary(modelGDD.spr)
drop1(modelGDD.spr, test = "Chi")

chkres(modelGDD.spr, popns_GDD$PEAKDAY, popns_GDD$CLASS)

r.squaredGLMM(modelGDD.spr)

# separate tests for univoltines and multivoltines
# split the dataset into univoltine and multivoltine
GDD.uv <- popns_GDD[which(popns_GDD$VOLTINISM == 1), ]
GDD.mv <- popns_GDD[which(popns_GDD$VOLTINISM == 2), ]

# and hs and wc
GDD.uv.hs <- GDD.uv[which(GDD.uv$CLASS == "HS"), ]
GDD.mv.hs <- GDD.mv[which(GDD.mv$CLASS == "HS"), ]
GDD.uv.wc <- GDD.uv[which(GDD.uv$CLASS == "WC"), ]
GDD.mv.wc <- GDD.mv[which(GDD.mv$CLASS == "WC"), ]


# univoltine
modelGDD.u <- lmer(PEAKDAY ~ SPRING.ONLY + (1|COMMON_NAME) + (1|YEAR),
                   data = GDD.uv)

summary(modelGDD.u)
drop1(modelGDD.u, test = "Chi")

# multivoltine
modelGDD.m <- lmer(PEAKDAY ~ SPRING.ONLY + (1|COMMON_NAME) + (1|YEAR),
                   data = GDD.mv)

summary(modelGDD.m)
drop1(modelGDD.m, test = "Chi")

# univoltine hs
modelGDD.u.hs <- lmer(PEAKDAY ~ SPRING.ONLY + (1|COMMON_NAME) + (1|YEAR),
                      data = GDD.uv.hs)

summary(modelGDD.u.hs)
drop1(modelGDD.u.hs, test = "Chi")

# univoltine wc
modelGDD.u.wc <- lmer(PEAKDAY ~ SPRING.ONLY + (1|COMMON_NAME) + (1|YEAR),
                      data = GDD.uv.wc)

summary(modelGDD.u.wc)
drop1(modelGDD.u.wc, test = "Chi")

# multivoltine hs
modelGDD.m.hs <- lmer(PEAKDAY ~ SPRING.ONLY + (1|COMMON_NAME) + (1|YEAR),
                      data = GDD.mv.hs)

summary(modelGDD.m.hs)
drop1(modelGDD.m.hs, test = "Chi")

# multivoltine wc
modelGDD.m.wc <- lmer(PEAKDAY ~ SPRING.ONLY + (1|COMMON_NAME) + (1|YEAR),
                      data = GDD.mv.wc)

summary(modelGDD.m.wc)
drop1(modelGDD.m.wc, test = "Chi")




### we want a single plot showing both univoltine and multivoltine responses to spring-only temps
# now construct the models

errorsGDD <- data.frame(effect(c("SPRING.ONLY*VOLTINISM"),modelGDD.spr))
errorsGDD$VOLTINISM <- as.factor(errorsGDD$VOLTINISM)

# build the plots

fig2d <- ggplot(errorsGDD)+
  geom_line(aes(y = fit, x = SPRING.ONLY, colour = VOLTINISM),
            linetype = "solid", stat="identity",size = 1)+
  geom_ribbon(aes(ymin = lower, ymax = upper, x = SPRING.ONLY, fill = VOLTINISM),
              stat="identity",size = 0.6, alpha = 0.2)+
  theme_bw()+
  scale_colour_manual(name = "Voltinism",
                      values = c("royalblue","goldenrod"),
                      labels = c("Univoltine","Multivoltine"))+
  scale_fill_manual(name = "Voltinism",
                    values = c("royalblue","goldenrod"),
                    labels = c("Univoltine","Multivoltine"))+
  ggtitle(" ")+
  ylab("Average peak day of emergence\n")+
  xlab("Accumulated degrees (GDD5), \nMarch 1st to May 31st")+
  theme(plot.title = element_text(hjust = 0.5))


fig2d



### so, hot spring/summer temps do lead to early emergence
## what is the consequence of that for various components of abundance?

popns_all_volti$VOLTINISM <- as.factor(popns_all_volti$VOLTINISM)
summary(popns_all_volti)

# specifically, we want to test at an intraspecific level whether early emergence (i.e. low PEAKDAY) 
# correlates with any of the following:

# total abundance
# first-generation abundance
# second-generation abundance
# timing of second-generation peak
# ratio of first- to second-generation abundance
# (also testing whether first-generation abundance correlates with second-generation abundance)
# abundance in year t+1

# that's a lot of analyses so let's just get cracking!
# split the dataset into uni and multivoltines

popns.uv <- popns_all_volti[which(popns_all_volti$VOLTINISM == 1), ]
popns.mv <- popns_all_volti[which(popns_all_volti$VOLTINISM == 2), ]

popns.uv.hs <- popns.uv[which(popns.uv$CLASS == "HS"), ]
popns.mv.hs <- popns.mv[which(popns.mv$CLASS == "HS"), ]
popns.uv.wc <- popns.uv[which(popns.uv$CLASS == "WC"), ]
popns.mv.wc <- popns.mv[which(popns.mv$CLASS == "WC"), ]


## total abundance
model.f1.c <- lmer(log(MEAN.ABUND) ~ PEAKDAY * VOLTINISM + PEAKDAY * CLASS + (1|COMMON_NAME) + (1|YEAR),
                 data = popns_all_volti)

summary(model.f1.c)
drop1(model.f1.c, test = "Chi")

r.squaredGLMM(model.f1.c)

# univoltine, hs
model.f1.u.hs <- lmer(log(MEAN.ABUND) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.uv.hs)

summary(model.f1.u.hs)
drop1(model.f1.u.hs, test = "Chi")

# multivoltine, hs
model.f1.m.hs <- lmer(log(MEAN.ABUND) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                      data = popns.mv.hs)

summary(model.f1.m.hs)
drop1(model.f1.m.hs, test = "Chi")


# univoltine, wc
model.f1.u.wc <- lmer(log(MEAN.ABUND) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                      data = popns.uv.wc)

summary(model.f1.u.wc)
drop1(model.f1.u.wc, test = "Chi")

# multivoltine, wc
model.f1.m.wc <- lmer(log(MEAN.ABUND) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                      data = popns.mv.wc)

summary(model.f1.m.wc)
drop1(model.f1.m.wc, test = "Chi")

## no class
model.f1 <- lmer(log(MEAN.ABUND) ~ PEAKDAY * VOLTINISM + (1|COMMON_NAME) + (1|YEAR),
                   data = popns_all_volti)

summary(model.f1)
drop1(model.f1, test = "Chi")

r.squaredGLMM(model.f1)


# univoltine
model.f1.u <- lmer(log(MEAN.ABUND) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.uv)

summary(model.f1.u)
drop1(model.f1.u, test = "Chi")

# multivoltine
model.f1.m <- lmer(log(MEAN.ABUND) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.mv)

summary(model.f1.m)
drop1(model.f1.m, test = "Chi")



## also create a single plot with univoltine and multivoltine overlaid

model.f1 <- lmer(log(MEAN.ABUND) ~ PEAKDAY * VOLTINISM + (1|COMMON_NAME) + (1|YEAR),
                 data = popns_all_volti)

summary(model.f1)
drop1(model.f1, test = "Chi")


errors.f1 <- data.frame(effect(c("PEAKDAY*VOLTINISM"),model.f1, xlevels=365))
errors.f1$fite <- exp(errors.f1$fit)
errors.f1$lowere <- exp(errors.f1$lower)
errors.f1$uppere <- exp(errors.f1$upper)



fig2a <- ggplot(errors.f1)+
  geom_ribbon(aes(ymin = lowere, ymax = uppere, x = PEAKDAY, fill = VOLTINISM),
              stat="identity",size = 0.6, alpha = 0.2)+
  geom_line(aes(y = fite, x = PEAKDAY, colour = VOLTINISM),
            stat="identity",size = 1)+
  theme_bw()+
  scale_colour_manual(name = "Voltinism",
                      values = c("royalblue","goldenrod"),
                      labels = c("Univoltine","Multivoltine"))+
  scale_fill_manual(name = "Voltinism",
                      values = c("royalblue","goldenrod"),
                      labels = c("Univoltine","Multivoltine"))+
  ggtitle(" ")+
  scale_y_log10(limits = c(0.05,2))+
  ylab("Mean abundance per transect/trap \nin given year (log scale)")+
  xlab("Peak day of emergence \n")+
  theme(plot.title = element_text(hjust = 0.5))


fig2a



## first-generation abundance
popns.f2 <- popns_all_volti[which(popns_all_volti$FIRST.GEN != "Fail"), ]
popns.uv.f2 <- popns.uv[which(popns.uv$FIRST.GEN != "Fail"), ]
popns.mv.f2 <- popns.mv[which(popns.mv$FIRST.GEN != "Fail"), ]

popns.f2$FIRST.GEN <- as.numeric(as.character(popns.f2$FIRST.GEN))
popns.uv.f2$FIRST.GEN <- as.numeric(as.character(popns.uv.f2$FIRST.GEN))
popns.mv.f2$FIRST.GEN <- as.numeric(as.character(popns.mv.f2$FIRST.GEN))

hist(log(popns.f2$FIRST.GEN))
hist(log(popns.uv.f2$FIRST.GEN))
hist(log(popns.mv.f2$FIRST.GEN))


model.f2 <- lmer(log(FIRST.GEN) ~ PEAKDAY * VOLTINISM * CLASS + (1|COMMON_NAME) + (1|YEAR),
                 data = popns.f2)

summary(model.f2)
drop1(model.f2, test = "Chi")


## second-generation abundance
popns.f3 <- popns_all_volti[which(popns_all_volti$SECOND.GEN != "Fail" & popns_all_volti$SECOND.GEN != 0), ]
popns.uv.f3 <- popns.uv[which(popns.uv$SECOND.GEN != "Fail" & popns.uv$SECOND.GEN != 0), ]
popns.mv.f3 <- popns.mv[which(popns.mv$SECOND.GEN != "Fail" & popns.uv$SECOND.GEN != 0), ]

popns.f3$SECOND.GEN <- as.numeric(as.character(popns.f3$SECOND.GEN))
popns.uv.f3$SECOND.GEN <- as.numeric(as.character(popns.uv.f3$SECOND.GEN))
popns.mv.f3$SECOND.GEN <- as.numeric(as.character(popns.mv.f3$SECOND.GEN))


hist(log(popns.f3$SECOND.GEN))
hist(log(popns.uv.f3$SECOND.GEN))
hist(log(popns.mv.f3$SECOND.GEN))

# overall
model.f3 <- lmer(log(SECOND.GEN) ~ PEAKDAY * VOLTINISM * CLASS  + (1|COMMON_NAME) + (1|YEAR),
                 data = popns.f3)

summary(model.f3)
drop1(model.f3, test = "Chi")


# timing of second-generation peak
popns.f4 <- popns_all_volti[which(popns_all_volti$SECOND.PEAK != "Fail"), ]
popns.uv.f4 <- popns.uv[which(popns.uv$SECOND.PEAK != "Fail"), ]
popns.mv.f4 <- popns.mv[which(popns.mv$SECOND.PEAK != "Fail"), ]

popns.f4$SECOND.PEAK <- as.numeric(as.character(popns.f4$SECOND.PEAK))
popns.uv.f4$SECOND.PEAK <- as.numeric(as.character(popns.uv.f4$SECOND.PEAK))
popns.mv.f4$SECOND.PEAK <- as.numeric(as.character(popns.mv.f4$SECOND.PEAK))

hist(popns.f4$SECOND.PEAK)
hist(popns.uv.f4$SECOND.PEAK)
hist(popns.mv.f4$SECOND.PEAK)

# overall
model.f4 <- lmer(SECOND.PEAK ~ PEAKDAY * VOLTINISM * CLASS + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.f4)

summary(model.f4)
drop1(model.f4, test = "Chi")



## does first-generation and second-generation abundance correlate?
popns.f5 <- popns.f2[which(popns.f2$SECOND.GEN != "Fail" & popns.f2$SECOND.GEN != 0), ]
popns.uv.f5 <- popns.uv.f2[which(popns.uv.f2$SECOND.GEN != "Fail" & popns.uv.f2$SECOND.GEN != 0), ]
popns.mv.f5 <- popns.mv.f2[which(popns.mv.f2$SECOND.GEN != "Fail" & popns.uv.f2$SECOND.GEN != 0), ]

popns.f5$SECOND.GEN <- as.numeric(as.character(popns.f5$SECOND.GEN))
popns.uv.f5$SECOND.GEN <- as.numeric(as.character(popns.uv.f5$SECOND.GEN))
popns.mv.f5$SECOND.GEN <- as.numeric(as.character(popns.mv.f5$SECOND.GEN))

hist(log(popns.f5$SECOND.GEN))
hist(log(popns.uv.f5$SECOND.GEN))
hist(log(popns.mv.f5$SECOND.GEN))

# overall
model.f5 <- lmer(SECOND.GEN ~ FIRST.GEN + VOLTINISM + CLASS + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.f5)

summary(model.f5)
drop1(model.f5, test = "Chi")



## ratio of first- to second-generation abundance
popns.f6 <- popns_all_volti[which(popns_all_volti$VOLTINISM.RATIO != "Fail" 
                                  & popns_all_volti$VOLTINISM.RATIO != 0
                                  & popns_all_volti$VOLTINISM == 2), ]
popns.uv.f6 <- popns.uv[which(popns.uv$VOLTINISM.RATIO != "Fail" & popns.uv$VOLTINISM.RATIO != 0), ]
popns.mv.f6 <- popns.mv[which(popns.mv$VOLTINISM.RATIO != "Fail" & popns.mv$VOLTINISM.RATIO != 0), ]

popns.f6$VOLTINISM.RATIO <- as.numeric(as.character(popns.f6$VOLTINISM.RATIO))
popns.uv.f6$VOLTINISM.RATIO <- as.numeric(as.character(popns.uv.f6$VOLTINISM.RATIO))
popns.mv.f6$VOLTINISM.RATIO <- as.numeric(as.character(popns.mv.f6$VOLTINISM.RATIO))

hist(log(popns.f6$VOLTINISM.RATIO))
hist(log(popns.uv.f6$VOLTINISM.RATIO))
hist(log(popns.mv.f6$VOLTINISM.RATIO))

# overall
model.f6.c <- lmer(log(VOLTINISM.RATIO) ~ PEAKDAY * CLASS + (1|COMMON_NAME) + (1|YEAR),
                 data = popns.f6)

summary(model.f6.c)
drop1(model.f6.c, test = "Chi")

r.squaredGLMM(model.f6.c)

## no class

model.f6 <- lmer(log(VOLTINISM.RATIO) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                 data = popns.f6)

summary(model.f6)
drop1(model.f6, test = "Chi")

r.squaredGLMM(model.f6)


# hs
popns.f6.hs <- popns.f6[which(popns.f6$CLASS == "HS"), ]
popns.f6.wc <- popns.f6[which(popns.f6$CLASS == "WC"), ]


model.f6.hs <- lmer(log(VOLTINISM.RATIO) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                    data = popns.f6.hs)

summary(model.f6.hs)
drop1(model.f6.hs, test = "Chi")


# wc
model.f6.wc <- lmer(log(VOLTINISM.RATIO) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                    data = popns.f6.wc)

summary(model.f6.wc)
drop1(model.f6.wc, test = "Chi")



## also create a version of g.f6.m with formatting matching the other final plots
summary(model.f6)

errors.f6 <- data.frame(effect(c("PEAKDAY"),model.f6, xlevels=365))
errors.f6$fite <- exp(errors.f6$fit)
errors.f6$lowere <- exp(errors.f6$lower)
errors.f6$uppere <- exp(errors.f6$upper)



fig2c <- ggplot(errors.f6)+
  geom_ribbon(aes(ymin = lowere, ymax = uppere, x = PEAKDAY),
              fill = "goldenrod", stat="identity",size = 0.6, alpha = 0.2)+
  geom_line(aes(y = fite, x = PEAKDAY),
            stat="identity",size = 1, colour = "goldenrod")+
  theme_bw()+
  ggtitle(" ")+
  scale_y_log10(limits = c(0.02,50))+
  ylab("Ratio of 1st to 2nd \ngeneration abundance")+
  xlab("Peak day of emergence \n")+
  theme(plot.title = element_text(hjust = 0.5))

fig2c

  
# abundance in year t+1
# first we need to generate this variable of t+1 abundance

head(popns_all_volti)

all.popns <- levels(droplevels(popns_all_volti$POPULATION))

popns.f7.t1 <- data.frame(YEAR = numeric(),
                          NEXT.YEAR = numeric(),
                          POPULATION = factor())

  
  
for (x in all.popns){
  pop <- popns_all_volti[which(popns_all_volti$POPULATION == x), ]
  
  YEAR <- NULL
  NEXT.YEAR <- NULL
  
  for (y in 1995:2014){
    YEAR <- append(YEAR, y)
    
    next.y <- pop[which(pop$YEAR == (1+y)), 'MEAN.ABUND']
    next.y <- ifelse(length(next.y) == 0, NA, next.y)
    
    NEXT.YEAR <- append(NEXT.YEAR, next.y)
    
  }
  
  out <- data.frame(cbind(YEAR, NEXT.YEAR))
  out$POPULATION <- x
  
  popns.f7.t1 <- rbind(popns.f7.t1,out)
  
}


# and merge it back into the main data
popns.f7 <- merge(popns.f7.t1, popns_all_volti, all.y=T)
summary(popns.f7)

popns.f7.uv <- popns.f7[which(popns.f7$VOLTINISM == 1), ]
popns.f7.mv <- popns.f7[which(popns.f7$VOLTINISM == 2), ]

popns.f7.uv.hs <- popns.f7.uv[which(popns.f7.uv$CLASS == "HS"), ]
popns.f7.mv.hs <- popns.f7.mv[which(popns.f7.mv$CLASS == "HS"), ]
popns.f7.uv.wc <- popns.f7.uv[which(popns.f7.uv$CLASS == "WC"), ]
popns.f7.mv.wc <- popns.f7.mv[which(popns.f7.mv$CLASS == "WC"), ]


hist(log(popns.f7$NEXT.YEAR))


# overall
model.f7.c <- lmer(log(NEXT.YEAR) ~ PEAKDAY * VOLTINISM * CLASS + (1|COMMON_NAME) + (1|YEAR),
                 data = popns.f7)

summary(model.f7.c)
drop1(model.f7.c, test = "Chi")

r.squaredGLMM(model.f7.c)

## no class
model.f7 <- lmer(log(NEXT.YEAR) ~ PEAKDAY * VOLTINISM + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.f7)

summary(model.f7)
drop1(model.f7, test = "Chi")

r.squaredGLMM(model.f7)


# univoltine
model.f7.u <- lmer(log(NEXT.YEAR) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.f7.uv)

summary(model.f7.u)
drop1(model.f7.u, test = "Chi")


# multivoltine
model.f7.m <- lmer(log(NEXT.YEAR) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.f7.mv)

summary(model.f7.m)
drop1(model.f7.m, test = "Chi")



# univoltine, hs
model.f7.u.hs <- lmer(log(NEXT.YEAR) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                   data = popns.f7.uv.hs)

summary(model.f7.u.hs)
drop1(model.f7.u.hs, test = "Chi")

# univoltine, wc
model.f7.u.wc <- lmer(log(NEXT.YEAR) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                      data = popns.f7.uv.wc)

summary(model.f7.u.wc)
drop1(model.f7.u.wc, test = "Chi")


# multivoltine, hs
model.f7.m.hs <- lmer(log(NEXT.YEAR) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                      data = popns.f7.mv.hs)

summary(model.f7.m.hs)
drop1(model.f7.m.hs, test = "Chi")

# multivoltine, wc
model.f7.m.wc <- lmer(log(NEXT.YEAR) ~ PEAKDAY + (1|COMMON_NAME) + (1|YEAR),
                      data = popns.f7.mv.wc)

summary(model.f7.m.wc)
drop1(model.f7.m.wc, test = "Chi")





## also create a single plot with univoltine and multivoltine overlaid
summary(model.f7)

errors.f7 <- data.frame(effect(c("PEAKDAY*VOLTINISM"),model.f7, xlevels=365))
errors.f7$fite <- exp(errors.f7$fit)
errors.f7$lowere <- exp(errors.f7$lower)
errors.f7$uppere <- exp(errors.f7$upper)



fig2b <- ggplot(errors.f7)+
  geom_ribbon(aes(ymin = lowere, ymax = uppere, x = PEAKDAY, fill = VOLTINISM),
              stat="identity",size = 0.6, alpha = 0.2)+
  geom_line(aes(y = fite, x = PEAKDAY, colour = VOLTINISM),
            stat="identity",size = 1)+
  theme_bw()+
  scale_fill_manual(name = "Voltinism",
                    values = c("royalblue","goldenrod"),
                    labels = c("Univoltine","Multivoltine"))+
  scale_colour_manual(name = "Voltinism",
                      values = c("royalblue","goldenrod"),
                      labels = c("Univoltine","Multivoltine"))+
  ggtitle(" ")+
  scale_y_log10(limits = c(0.05,2))+
  ylab("Mean abundance per transect/trap \nin following year (log scale)")+
  xlab("Peak day of emergence \n")+
  theme(plot.title = element_text(hjust = 0.5))

fig2b


### final figures ####
## we want to rearrange some components of figures into new sets for the final manuscript

# Fig 1

fig1a <- grid.arrange(arrangeGrob(g1.u.heavy + theme(legend.position="none") + theme(panel.grid = element_blank()) + xlab("\n"),g1.m.heavy + theme(legend.position="none") + theme(panel.grid = element_blank()) + ylab("\n") + xlab("\n"),ncol=2),
                      leg1.heavy,ncol=2,widths=c(10,2))

fig1b <- grid.arrange(arrangeGrob(gi.u + theme(legend.position="none") + theme(panel.grid = element_blank()) + ggtitle(" "),gi.mv + theme(legend.position="none") + theme(panel.grid = element_blank()) + ggtitle(" ") + ylab("\n"),ncol=2),
                      leg.i,ncol=2, widths=c(10,2))


fig1 <- grid.arrange(arrangeGrob(fig1a,fig1b, nrow = 2))

ggsave("Results/Plots/Final/Fig2update.svg", plot = fig1, device = svg, width = 250, height = 200, units = "mm", limitsize = F)



# Fig 2


fig2a <- fig2a + theme(panel.grid = element_blank())
fig2b <- fig2b + theme(panel.grid = element_blank())
fig2c <- fig2c + theme(panel.grid = element_blank())
fig2d <- fig2d + theme(panel.grid = element_blank())

leg2.final <- g_legend(fig2a)



# export

fig2 <- grid.arrange(arrangeGrob(fig2d + theme(legend.position="none"), fig2c + theme(legend.position="none"),
                                 fig2a + theme(legend.position="none"), fig2b + theme(legend.position="none"),
                                 ncol=2),
                     leg2.final, ncol=2, widths=c(10,2))

ggsave("Results/Plots/Final/Fig3update.svg", plot = fig2, device = svg, width = 250, height = 200, units = "mm", limitsize = F)





### check of national relevance of abundance trends

# there is a question over the relevance of the abundance trends calculated for our species 
# (from, in some cases, a small number of sites)
# to the national picture

# we want to use national-scale data from the UKBMS and the RIS to check this out -
# we'll have to do it separately for butterflies and moths as the format of the two datasets can't be married easily

# first, for the UKBMS: we have a single national index of abundance per species per year

BMS.natabund <- read.csv("Data/Collated_Indices_2016.csv", header = T)

summary(BMS.natabund)

# write a loop over this to extract:
# (i) national trend over full time period
# (ii) national trend over 1995-2014

BMS.species <- levels(droplevels(BMS.natabund$Species))

BMS.nattrends <- data.frame(SCI_NAME = factor(),
                            COMMON_NAME = factor(),
                            FULL.TREND = numeric(),
                            STUDY.TREND = numeric())

for (x in BMS.species){
  spec <- BMS.natabund[which(BMS.natabund$Species == x), ]
  
  SCI_NAME <- x
  COMMON_NAME <- as.character(spec[[1,2]])
  
  mod1 <- lm(Collated.Index ~ Year, data = spec)

  FULL.TREND <- summary(mod1)$coefficients[2,1]
  
  mod2 <- lm(Collated.Index ~ Year, data = spec[which(spec$Year %in% 1995:2014), ])
  
  STUDY.TREND <- summary(mod2)$coefficients[2,1]
  
  out <- cbind(SCI_NAME, COMMON_NAME, FULL.TREND, STUDY.TREND)
  
  BMS.nattrends <- rbind(BMS.nattrends, out)
}

BMS.nattrends$FULL.TREND <- as.numeric(as.character(BMS.nattrends$FULL.TREND))
BMS.nattrends$STUDY.TREND <- as.numeric(as.character(BMS.nattrends$STUDY.TREND))

# now merge our own trends into these

BMS.alltrends <- merge(BMS.nattrends, data.volt.heavy[,c('COMMON_NAME','SITES','ABUND.SLOPE')])

summary(BMS.alltrends)

## now we want to test the correlation between our trends and the national ones

# first for the study period only
plot(STUDY.TREND ~ ABUND.SLOPE, data = BMS.alltrends)

model.st.b <- lm(STUDY.TREND ~ ABUND.SLOPE, data = BMS.alltrends)

summary(model.st.b)
drop1(model.st.b, test = "F")

# and for the full period
plot(FULL.TREND ~ ABUND.SLOPE, data = BMS.alltrends)

model.ft.b <- lm(FULL.TREND ~ ABUND.SLOPE, data = BMS.alltrends)

summary(model.ft.b)
drop1(model.ft.b, test = "F")


## repeat the same test for RIS data - here I've pulled in full trends (period 1967-2016) and can't obtain shorter trends

RIS.nattrends <- read.csv("Data/RIS_nattrends.csv", header = T)

summary(RIS.nattrends)

plot(FULL.TREND ~ ABUND.SLOPE, data = RIS.nattrends)

model.ft.m <- lm(FULL.TREND ~ ABUND.SLOPE, data = RIS.nattrends)

summary(model.ft.m)
drop1(model.ft.m, test = "F")
