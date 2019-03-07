############################################################
####   Script for analysing and plotting climate data   ####
############################################################

### Clear the current workspace (DOESN'T CLEAR LOADED LIBRARIES)
rm(list=ls())

### install if necessary and then load the libraries you need

j <- c("rstudioapi","plyr","dplyr","ggplot2","RColorBrewer","reshape2","gridExtra")

new.packages <- j[!(j %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(j, require, character.only = TRUE)  # loads up any libraries that aren't already loaded

### set working directory to script's saved location

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read in data
climate <- read.csv("Data/UK_climate.csv", header = T)

summary(climate)

# we have monthly, quarterly (seasonal) and annual mean temperature data for the years 1910-2018

# we want to use this to assess the change in climate immediately before and during our study period
prelude <- 1980:1994
study <- 1995:2014
interval1 <- 1995:2006
interval2 <- 2003:2014
postlude <- 2015:2018

# first pick out the data for the years we want
all.years <- climate[which(climate$Year %in% c(prelude,study,postlude)), ]

# now pick out the study years
study.years <- climate[which(climate$Year %in% study), ]

# and plot a line through them
mod1 <- lm(SPR ~ Year, data = study.years)
summary(mod1)
drop1(mod1, test = "F")

preds <- data.frame(Year = c(1995:2014))
preds$SPR <- predict(mod1, newdata = preds)



# and start building a plot
g1 <- ggplot(all.years, aes(x = Year, y = SPR))+
  geom_smooth(linetype = "dashed")+
  theme_bw()+
  geom_point() +
  theme(panel.grid = element_blank())+
  geom_line()+
  geom_line(data = preds, aes(x = Year, y = SPR),
            colour = "darkgreen", size = 1)+
  ylab("Spring temperature \n(degrees Celsius)\n")+
  xlab(" ")+
  ylim(5.5,9.5)+
  xlim(1980,2020)


g1



### for comparison to the mean spring temperatures, let's do the same for the mean GDD up to spring at the set of study sites, 
# as this is the *actual* climate variable we use in the analysis

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

site_GDD_wide$SPRING.ONLY <- site_GDD_wide$SUMMER - site_GDD_wide$SPRING

## this variable is now held in the column called "SPRING", so let's take an annual average across all sites of it

site_GDD_wide_comp <- site_GDD_wide[complete.cases(site_GDD_wide),]

GDD_mean <- ddply(site_GDD_wide_comp, .(YEAR), numcolwise(mean))


## and now prepare the panel

# prep a trend-line
mod2 <- lm(SPRING.ONLY ~ YEAR, data = GDD_mean)
summary(mod2)
drop1(mod2, test = "F")

preds2 <- data.frame(YEAR = c(1995:2014))
preds2$SPRING.ONLY <- predict(mod2, newdata = preds2)


# and start building a plot
g2 <- ggplot(GDD_mean, aes(x = YEAR, y = SPRING.ONLY))+
  geom_smooth(linetype = "dashed")+
  theme_bw() +
  theme(panel.grid = element_blank())+
  geom_point()+
  geom_line()+
  geom_line(data = preds2, aes(x = YEAR, y = SPRING.ONLY),
            colour = "darkgreen", size = 1)+
  ylab("Spring temperature \n(GDD5)")+
  xlab("Year")+
  ylim(200,550)+
  xlim(1980,2020)


g2


# merge the panels

m1 <- grid.arrange(g1,g2)


ggsave("Results/Plots/Climate.svg", plot = m1, device = svg, width = 120, height = 140, units = "mm", limitsize = F)

