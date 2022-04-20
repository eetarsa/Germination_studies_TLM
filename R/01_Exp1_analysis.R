# Experiment 1 ----
rm(list=ls())

# Note: before ET loaded in data, I changed replicate numbers in the data files prior to loading. 
# For all species, there were actually 6 replicates of germination tested at 0 MPA per temp (across 2 trials). 
# So all 0 MPA treatments, rep = 6; everything else, rep = 3;
# could also drop 3 reps to balance design if that's a better approach

# Load dependencies ----
library(GerminaR)
library(tidyverse)
library(emmeans)
library(grid)
library(reshape2)
library(dplyr)
library(agricolae)
library(gtools)
library(lme4)
library(glmmTMB)
#devtools::install_github("onofriAndreaPG/drcSeedGerm")
#devtools::install_github("DoseResponse/drcData")
library(drcData)
library(drcSeedGerm)
library(drc)
library(magic)
library(car)
library(metafor)
library(multcomp)
library(ggplot2); theme_set(theme_bw())
library(broom)
library(knitr)
library(markdown)
library(gtable)

# Change contrast options
options(contrast = c("contr.sum", "contr.poly"))

# Load cleaned (non-cumulative) germination data ----
phau.exp1 <- read.csv("./data/exp1_PHAU_noncum.csv", as.is=T)
scac.exp1 <- read.csv("./data/exp1_SCAC_noncum.csv", as.is=T)
scam.exp1 <- read.csv("./data/exp1_SCAM_noncum.csv", as.is=T)

# Merge & clean data sets ----
all.spp <- rbind(phau.exp1, scac.exp1, scam.exp1)
all.spp$temp <- as.factor(all.spp$temp)
all.spp$species <- as.factor(all.spp$species)
all.spp$MPA <- as.factor(all.spp$MPA)
all.spp$rep <- as.factor(all.spp$rep)

# Add day 0 column full of zero's to initialize germination & reorder dataframe
all.spp$Day.0 <- 0
all.spp <- all.spp[, c(1:5,27,6:26)]

# Remove "Day." for data collection columns
colnames(all.spp) <- c("temp", "species", "MPA", "rep", "N",
                       "0","3","4","5","6","7","8","9","10","11",
                       "12","13","14","15","16","18","20","22","24",
                       "26", "28", "30")

# Investigate values that are negative
all.spp$trt.id <- paste(all.spp$temp, all.spp$species, all.spp$MPA, all.spp$rep)
length(unique(all.spp$trt.id))
dim(all.spp) # check that dimension and id length match

neg.val <- all.spp[rowSums(all.spp[6:27] < 0) > 0, ]
neg.val # 6 instances of negative values - check with Tatiana, input error?

# For now, assume input error and make negative values positive (but check with Tatiana)
all.spp[,6:27] <- abs(all.spp[,6:27])

# Remove trt.id column
all.spp <- all.spp[,1:27]

# Add replicate column that include all treatment factors
all.spp1 <- all.spp
all.spp1[["Rep"]] <- with(all.spp1, interaction(temp, species, MPA, rep))

# Flip data set
spp.long <- pivot_longer(all.spp1, cols="0":"30", names_to = "Day", values_to = "Germ")

# Add time interval
spp.long$Start <- spp.long$Day
spp.long$End <- NA

for (i in 1:nrow(spp.long)){
  if(spp.long$Day[i]=="30"){
    spp.long$End[i] <- "Inf"
  }else{
    spp.long$End[i] <- spp.long$Day[i+1]
  }#ifelse
}#for

spp.long$Start <- as.numeric(spp.long$Start)
spp.long$End <- as.numeric(spp.long$End)

# Calculate cumulative germination
spp.long$CumGerm <- NA
for (i in 1:nrow(spp.long)){
  if(spp.long$Day[i]=="0"){
    spp.long$CumGerm[i] <- spp.long$Germ[i]
  }else{
    spp.long$CumGerm[i] <- spp.long$Germ[i] + spp.long$CumGerm[i-1]
  }#ifelse
}#for

# Calculate cumulative germination proportion
spp.long$CumGRP <- spp.long$CumGerm / spp.long$N

# Turn day into factor and level correctly
spp.long$Day <- as.factor(spp.long$Day)
day.vec <- c(0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30)
spp.long$Day <- factor(spp.long$Day, levels = day.vec)

# Assign remaining ungerminated seeds to "Inf" time interval
for (i in 1:nrow(spp.long)){
  if (spp.long$End[i] == "Inf"){
    spp.long$Germ[i] <- spp.long$N[i] - spp.long$CumGerm[i]
  } 
}

# Plot germination curves - clean this up! ----
ggplot(spp.long, aes(x = Day, y = CumGRP, color=rep, shape=MPA)) + 
  geom_point(size=2.5) + ylim(0,1) +
  geom_line() + facet_wrap(~species*temp) + labs(x ="Days", y = "Germination proportion") +
  scale_colour_grey(start = 0, end = .6) + scale_shape_manual(values = c(0, 16, 4, 6, 15)) +
  theme_bw()

germ.stats <- spp.long %>%
  group_by(temp, species, MPA, Day) %>% 
  dplyr::summarize(Mean  = mean(CumGRP, na.rm = TRUE),
                   SD    = sd(CumGRP, na.rm = TRUE),
                   nsize = sum(!is.na(CumGRP)),
                   SE    = SD/sqrt(nsize))

ggplot(germ.stats, aes(x = Day, y = Mean, group=MPA, color=MPA, shape=MPA)) + 
  geom_line(aes(group=MPA, color=MPA)) + 
  geom_point(size=2.5) + ylim(0,1) +
  facet_wrap(~species*temp) + labs(x ="Days", y = "Germination proportion") +
  scale_colour_grey(start = 0, end = 0.5) + scale_shape_manual(values = c(0, 16, 4, 6, 15)) +
  theme_bw()

ggplot(germ.stats, aes(x = Day, y = Mean, group=MPA, color=MPA, shape=MPA)) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.2) +
  geom_line(aes(group=MPA, color=MPA)) + 
  geom_point(size=2.5) + ylim(0,1) +
  facet_wrap(~species*temp) + labs(x ="Days", y = "Germination proportion") +
  scale_colour_grey(start = 0, end = 0.5) + scale_shape_manual(values = c(0, 16, 4, 6, 15)) +
  theme_bw()

# Model germination curves ----
## Split long data set
phau.long <- filter(spp.long, species == "PHAU")
scac.long <- filter(spp.long, species == "SCAC")
scam.long <- filter(spp.long, species == "SCAM")

## PHAU models ----
### Temp 23 ----
phau.23 <- filter(phau.long, temp == "23")
phau23.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(fixed = c(NA, NA, NA)), 
                       data = phau.23, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(phau23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 23 @ 0 MPA")
summary(phau23.mod.mpa0)
ph23.mpa0 <- as.data.frame(tidy(phau23.mod.mpa0))
ph23.mpa0$species <- "PHAU"
ph23.mpa0$temp <- "23"
ph23.mpa0$MPA <- "0"

phau23.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(phau23.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 23 @ -0.15 MPA")
summary(phau23.mod.mpa15)
ph23.mpa15 <- as.data.frame(tidy(phau23.mod.mpa15))
ph23.mpa15$species <- "PHAU"
ph23.mpa15$temp <- "23"
ph23.mpa15$MPA <- "-0.15"

phau23.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(phau23.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 23 @ -0.3 MPA")
summary(phau23.mod.mpa3)
ph23.mpa3 <- as.data.frame(tidy(phau23.mod.mpa3))
ph23.mpa3$species <- "PHAU"
ph23.mpa3$temp <- "23"
ph23.mpa3$MPA <- "-0.3"

phau23.mod.mpa6 <- drm(Germ ~ Start + End, fct = W1.3,
                       data = phau.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(phau23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 23 @ -0.6 MPA")
summary(phau23.mod.mpa6) 
ph23.mpa6 <- as.data.frame(tidy(phau23.mod.mpa6))
ph23.mpa6$species <- "PHAU"
ph23.mpa6$temp <- "23"
ph23.mpa6$MPA <- "-0.6"
# something is not working with replicate 3 - it continues to estimate max germ at above 100%
# the model does not converge or properly estimate when upper bound is set; tried LL.4, LL.3, W1.4, W2.4 with no luck
# Possible solution - drop 3rd replicate?
phau23.6.nor3 <- filter(phau.23, MPA == -0.6 & rep != 3)
phau23.mod.mpa6.V2 <- drm(Germ ~ Start + End, fct = W2.4(fixed = c(NA, 0, NA, NA)),
                       data = phau23.6.nor3, type = "event", 
                       curveid = rep)
plot(phau23.mod.mpa6.V2, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 23 @ -0.6 MPA")
summary(phau23.mod.mpa6.V2) 
ph23.mpa6 <- as.data.frame(tidy(phau23.mod.mpa6.V2))
ph23.mpa6$species <- "PHAU"
ph23.mpa6$temp <- "23"
ph23.mpa6$MPA <- "-0.6"

phau23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(phau23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 23 @ -1.2 MPA")
summary(phau23.mod.mpa12)
ph23.mpa12 <- as.data.frame(tidy(phau23.mod.mpa12))
ph23.mpa12$species <- "PHAU"
ph23.mpa12$temp <- "23"
ph23.mpa12$MPA <- "-1.2"

### Temp 28 ----
phau.28 <- filter(phau.long, temp == "28")
phau28.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.28, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(phau28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 28 @ 0 MPA")
summary(phau28.mod.mpa0)
ph28.mpa0 <- as.data.frame(tidy(phau28.mod.mpa0))
ph28.mpa0$species <- "PHAU"
ph28.mpa0$temp <- "28"
ph28.mpa0$MPA <- "0"

phau28.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.28, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(phau28.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 28 @ -0.15 MPA")
summary(phau28.mod.mpa15)
ph28.mpa15 <- as.data.frame(tidy(phau28.mod.mpa15))
ph28.mpa15$species <- "PHAU"
ph28.mpa15$temp <- "28"
ph28.mpa15$MPA <- "-0.15"

phau28.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(phau28.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 28 @ -0.3 MPA")
summary(phau28.mod.mpa3)
ph28.mpa3 <- as.data.frame(tidy(phau28.mod.mpa3))
ph28.mpa3$species <- "PHAU"
ph28.mpa3$temp <- "28"
ph28.mpa3$MPA <- "-0.3"

phau28.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(phau28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 28 @ -0.6 MPA")
summary(phau28.mod.mpa6)
ph28.mpa6 <- as.data.frame(tidy(phau28.mod.mpa6))
ph28.mpa6$species <- "PHAU"
ph28.mpa6$temp <- "28"
ph28.mpa6$MPA <- "-0.6"

phau28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.28, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(phau28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 28 @ -1.2 MPA")
summary(phau28.mod.mpa12)
ph28.mpa12 <- as.data.frame(tidy(phau28.mod.mpa12))
ph28.mpa12$species <- "PHAU"
ph28.mpa12$temp <- "28"
ph28.mpa12$MPA <- "-1.2"

### Temp 33 ----
phau.33 <- filter(phau.long, temp == "33")
phau33.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.33, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(phau33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 33 @ 0 MPA")
summary(phau33.mod.mpa0)
ph33.mpa0 <- as.data.frame(tidy(phau33.mod.mpa0))
ph33.mpa0$species <- "PHAU"
ph33.mpa0$temp <- "33"
ph33.mpa0$MPA <- "0"

phau33.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.33, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(phau33.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 33 @ -0.15 MPA")
summary(phau33.mod.mpa15)
ph33.mpa15 <- as.data.frame(tidy(phau33.mod.mpa15))
ph33.mpa15$species <- "PHAU"
ph33.mpa15$temp <- "33"
ph33.mpa15$MPA <- "-0.15"

phau33.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(phau33.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 33 @ -0.3 MPA")
summary(phau33.mod.mpa3)
ph33.mpa3 <- as.data.frame(tidy(phau33.mod.mpa3))
ph33.mpa3$species <- "PHAU"
ph33.mpa3$temp <- "33"
ph33.mpa3$MPA <- "-0.3"

phau33.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(phau33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 33 @ -0.6 MPA")
summary(phau33.mod.mpa6)
ph33.mpa6 <- as.data.frame(tidy(phau33.mod.mpa6))
ph33.mpa6$species <- "PHAU"
ph33.mpa6$temp <- "33"
ph33.mpa6$MPA <- "-0.6"

phau33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.33, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(phau33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 33 @ -1.2 MPA")
summary(phau33.mod.mpa12)
ph33.mpa12 <- as.data.frame(tidy(phau33.mod.mpa12))
ph33.mpa12$species <- "PHAU"
ph33.mpa12$temp <- "33"
ph33.mpa12$MPA <- "-1.2"

### Temp 36 ----
phau.36 <- filter(phau.long, temp == "36")
phau36.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.36, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(phau36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 36 @ 0 MPA")
summary(phau36.mod.mpa0)
ph36.mpa0 <- as.data.frame(tidy(phau36.mod.mpa0))
ph36.mpa0$species <- "PHAU"
ph36.mpa0$temp <- "36"
ph36.mpa0$MPA <- "0"

phau36.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.36, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(phau36.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 36 @ -0.15 MPA")
summary(phau36.mod.mpa15)
ph36.mpa15 <- as.data.frame(tidy(phau36.mod.mpa15))
ph36.mpa15$species <- "PHAU"
ph36.mpa15$temp <- "36"
ph36.mpa15$MPA <- "-0.15"

phau36.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(phau36.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 36 @ -0.3 MPA")
summary(phau36.mod.mpa3)
ph36.mpa3 <- as.data.frame(tidy(phau36.mod.mpa3))
ph36.mpa3$species <- "PHAU"
ph36.mpa3$temp <- "36"
ph36.mpa3$MPA <- "-0.3"

phau36.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = phau.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(phau36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 36 @ -0.6 MPA")
summary(phau36.mod.mpa6)
ph36.mpa6 <- as.data.frame(tidy(phau36.mod.mpa6))
ph36.mpa6$species <- "PHAU"
ph36.mpa6$temp <- "36"
ph36.mpa6$MPA <- "-0.6"

phau36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = phau.36, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(phau36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="PHAU 36 @ -1.2 MPA")
summary(phau36.mod.mpa12)
ph36.mpa12 <- as.data.frame(tidy(phau36.mod.mpa12))
ph36.mpa12$species <- "PHAU"
ph36.mpa12$temp <- "36"
ph36.mpa12$MPA <- "-1.2"

## SCAC models ----
### Temp 23 ----
scac.23 <- filter(scac.long, temp == "23")
scac23.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.23, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scac23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 23 @ 0 MPA")
summary(scac23.mod.mpa0)
s23.mpa0 <- as.data.frame(tidy(scac23.mod.mpa0))
s23.mpa0$species <- "SCAC"
s23.mpa0$temp <- "23"
s23.mpa0$MPA <- "0"

scac23.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scac23.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 23 @ -0.15 MPA")
summary(scac23.mod.mpa15)
s23.mpa15 <- as.data.frame(tidy(scac23.mod.mpa15))
s23.mpa15$species <- "SCAC"
s23.mpa15$temp <- "23"
s23.mpa15$MPA <- "-0.15"

scac23.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scac23.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 23 @ -0.3 MPA")
summary(scac23.mod.mpa3)
s23.mpa3 <- as.data.frame(tidy(scac23.mod.mpa3))
s23.mpa3$species <- "SCAC"
s23.mpa3$temp <- "23"
s23.mpa3$MPA <- "-0.3"

scac23.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scac23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 23 @ -0.6 MPA")
summary(scac23.mod.mpa6)
s23.mpa6 <- as.data.frame(tidy(scac23.mod.mpa6))
s23.mpa6$species <- "SCAC"
s23.mpa6$temp <- "23"
s23.mpa6$MPA <- "-0.6"

scac23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scac23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 23 @ -1.2 MPA")
summary(scac23.mod.mpa12)
s23.mpa12 <- as.data.frame(tidy(scac23.mod.mpa12))
s23.mpa12$species <- "SCAC"
s23.mpa12$temp <- "23"
s23.mpa12$MPA <- "-1.2"

### Temp 28 ----
scac.28 <- filter(scac.long, temp == "28")
scac28.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.28, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scac28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 28 @ 0 MPA")
summary(scac28.mod.mpa0)
s28.mpa0 <- as.data.frame(tidy(scac28.mod.mpa0))
s28.mpa0$species <- "SCAC"
s28.mpa0$temp <- "28"
s28.mpa0$MPA <- "0"

scac28.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.28, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scac28.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 28 @ -0.15 MPA")
summary(scac28.mod.mpa15)
s28.mpa15 <- as.data.frame(tidy(scac28.mod.mpa15))
s28.mpa15$species <- "SCAC"
s28.mpa15$temp <- "28"
s28.mpa15$MPA <- "-0.15"

scac28.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scac28.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 28 @ -0.3 MPA")
summary(scac28.mod.mpa3)
s28.mpa3 <- as.data.frame(tidy(scac28.mod.mpa3))
s28.mpa3$species <- "SCAC"
s28.mpa3$temp <- "28"
s28.mpa3$MPA <- "-0.3"

scac28.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scac28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 28 @ -0.6 MPA")
summary(scac28.mod.mpa6)
s28.mpa6 <- as.data.frame(tidy(scac28.mod.mpa6))
s28.mpa6$species <- "SCAC"
s28.mpa6$temp <- "28"
s28.mpa6$MPA <- "-0.6"

scac28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.28, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scac28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 28 @ -1.2 MPA")
summary(scac28.mod.mpa12)
s28.mpa12 <- as.data.frame(tidy(scac28.mod.mpa12))
s28.mpa12$species <- "SCAC"
s28.mpa12$temp <- "28"
s28.mpa12$MPA <- "-1.2"

### Temp 33 ----
scac.33 <- filter(scac.long, temp == "33")
scac33.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.33, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scac33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 33 @ 0 MPA")
summary(scac33.mod.mpa0)
s33.mpa0 <- as.data.frame(tidy(scac33.mod.mpa0))
s33.mpa0$species <- "SCAC"
s33.mpa0$temp <- "33"
s33.mpa0$MPA <- "0"

scac33.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.33, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scac33.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 33 @ -0.15 MPA")
summary(scac33.mod.mpa15)
s33.mpa15 <- as.data.frame(tidy(scac33.mod.mpa15))
s33.mpa15$species <- "SCAC"
s33.mpa15$temp <- "33"
s33.mpa15$MPA <- "-0.15"

scac33.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scac33.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 33 @ -0.3 MPA")
summary(scac33.mod.mpa3)
s33.mpa3 <- as.data.frame(tidy(scac33.mod.mpa3))
s33.mpa3$species <- "SCAC"
s33.mpa3$temp <- "33"
s33.mpa3$MPA <- "-0.3"

scac33.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scac33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 33 @ -0.6 MPA")
summary(scac33.mod.mpa6)
s33.mpa6 <- as.data.frame(tidy(scac33.mod.mpa6))
s33.mpa6$species <- "SCAC"
s33.mpa6$temp <- "33"
s33.mpa6$MPA <- "-0.6"

scac33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.33, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scac33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 33 @ -1.2 MPA")
summary(scac33.mod.mpa12)
s33.mpa12 <- as.data.frame(tidy(scac33.mod.mpa12))
s33.mpa12$species <- "SCAC"
s33.mpa12$temp <- "33"
s33.mpa12$MPA <- "-1.2"

### Temp 36 ----
scac.36 <- filter(scac.long, temp == "36")
scac36.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.36, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scac36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 36 @ 0 MPA")
summary(scac36.mod.mpa0)
s36.mpa0 <- as.data.frame(tidy(scac36.mod.mpa0))
s36.mpa0$species <- "SCAC"
s36.mpa0$temp <- "36"
s36.mpa0$MPA <- "0"

scac36.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.36, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scac36.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 36 @ -0.15 MPA")
summary(scac36.mod.mpa15)
s36.mpa15 <- as.data.frame(tidy(scac36.mod.mpa15))
s36.mpa15$species <- "SCAC"
s36.mpa15$temp <- "36"
s36.mpa15$MPA <- "-0.15"

scac36.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scac36.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 36 @ -0.3 MPA")
summary(scac36.mod.mpa3)
s36.mpa3 <- as.data.frame(tidy(scac36.mod.mpa3))
s36.mpa3$species <- "SCAC"
s36.mpa3$temp <- "36"
s36.mpa3$MPA <- "-0.3"

scac36.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scac.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scac36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 36 @ -0.6 MPA")
summary(scac36.mod.mpa6)
s36.mpa6 <- as.data.frame(tidy(scac36.mod.mpa6))
s36.mpa6$species <- "SCAC"
s36.mpa6$temp <- "36"
s36.mpa6$MPA <- "-0.6"

scac36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scac.36, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scac36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAC 36 @ -1.2 MPA")
summary(scac36.mod.mpa12)
s36.mpa12 <- as.data.frame(tidy(scac36.mod.mpa12))
s36.mpa12$species <- "SCAC"
s36.mpa12$temp <- "36"
s36.mpa12$MPA <- "-1.2"

## SCAM models ----
### Temp 23 ----
scam.23 <- filter(scam.long, temp == "23")
scam23.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.23, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scam23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 23 @ 0 MPA")
summary(scam23.mod.mpa0)
sm23.mpa0 <- as.data.frame(tidy(scam23.mod.mpa0))
sm23.mpa0$species <- "SCAM"
sm23.mpa0$temp <- "23"
sm23.mpa0$MPA <- "0"

scam23.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scam23.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 23 @ -0.15 MPA")
summary(scam23.mod.mpa15)
sm23.mpa15 <- as.data.frame(tidy(scam23.mod.mpa15))
sm23.mpa15$species <- "SCAM"
sm23.mpa15$temp <- "23"
sm23.mpa15$MPA <- "-0.15"

scam23.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scam23.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 23 @ -0.3 MPA")
summary(scam23.mod.mpa3)
sm23.mpa3 <- as.data.frame(tidy(scam23.mod.mpa3))
sm23.mpa3$species <- "SCAM"
sm23.mpa3$temp <- "23"
sm23.mpa3$MPA <- "-0.3"

scam23.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scam23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 23 @ -0.6 MPA")
summary(scam23.mod.mpa6)
sm23.mpa6 <- as.data.frame(tidy(scam23.mod.mpa6))
sm23.mpa6$species <- "SCAM"
sm23.mpa6$temp <- "23"
sm23.mpa6$MPA <- "-0.6"

scam23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scam23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 23 @ -1.2 MPA")
summary(scam23.mod.mpa12) # NOTHING GERMINATED HERE! 
# Create data frame with 0 for max germination and slope, set t50 (e) = max data collection interval (30 days)
term <- c("b", "b", "b", "d", "d", "d", "e", "e", "e")
curve <- c("1", "2", "3", "1", "2", "3", "1", "2", "3")
estimate <- c(0,0,0,0,0,0,30,30,30)
std.error <- rep(0, 9)
statistic <- rep(0, 9)
p.value <- rep(NA, 9)
species <- rep("SCAM", 9)
temp <- rep("23", 9)
MPA <- rep("-1.2", 9)
sm23.mpa12 <- data.frame(term, curve, estimate, std.error, statistic, p.value,
                         species, temp, MPA)

### Temp 28 ----
scam.28 <- filter(scam.long, temp == "28")
scam28.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.28, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scam28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 28 @ 0 MPA")
summary(scam28.mod.mpa0)
sm28.mpa0 <- as.data.frame(tidy(scam28.mod.mpa0))
sm28.mpa0$species <- "SCAM"
sm28.mpa0$temp <- "28"
sm28.mpa0$MPA <- "0"

scam28.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.28, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scam28.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 28 @ -0.15 MPA")
summary(scam28.mod.mpa15)
sm28.mpa15 <- as.data.frame(tidy(scam28.mod.mpa15))
sm28.mpa15$species <- "SCAM"
sm28.mpa15$temp <- "28"
sm28.mpa15$MPA <- "-0.15"

scam28.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scam28.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 28 @ -0.3 MPA")
summary(scam28.mod.mpa3)
sm28.mpa3 <- as.data.frame(tidy(scam28.mod.mpa3))
sm28.mpa3$species <- "SCAM"
sm28.mpa3$temp <- "28"
sm28.mpa3$MPA <- "-0.3"

scam28.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scam28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 28 @ -0.6 MPA")
summary(scam28.mod.mpa6)
sm28.mpa6 <- as.data.frame(tidy(scam28.mod.mpa6))
sm28.mpa6$species <- "SCAM"
sm28.mpa6$temp <- "28"
sm28.mpa6$MPA <- "-0.6"

scam28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.28, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scam28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 28 @ -1.2 MPA")
summary(scam28.mod.mpa12)
sm28.mpa12 <- as.data.frame(tidy(scam28.mod.mpa12))
sm28.mpa12$species <- "SCAM"
sm28.mpa12$temp <- "28"
sm28.mpa12$MPA <- "-1.2"

### Temp 33 ----
scam.33 <- filter(scam.long, temp == "33")
scam33.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.33, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scam33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 33 @ 0 MPA")
summary(scam33.mod.mpa0)
sm33.mpa0 <- as.data.frame(tidy(scam33.mod.mpa0))
sm33.mpa0$species <- "SCAM"
sm33.mpa0$temp <- "33"
sm33.mpa0$MPA <- "0"

scam33.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.33, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scam33.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 33 @ -0.15 MPA")
summary(scam33.mod.mpa15)
sm33.mpa15 <- as.data.frame(tidy(scam33.mod.mpa15))
sm33.mpa15$species <- "SCAM"
sm33.mpa15$temp <- "33"
sm33.mpa15$MPA <- "-0.15"

scam33.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scam33.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 33 @ -0.3 MPA")
summary(scam33.mod.mpa3)
sm33.mpa3 <- as.data.frame(tidy(scam33.mod.mpa3))
sm33.mpa3$species <- "SCAM"
sm33.mpa3$temp <- "33"
sm33.mpa3$MPA <- "-0.3"

scam33.mod.mpa6 <- drm(Germ ~ Start + End, fct = W1.3(),
                       data = scam.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scam33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 33 @ -0.6 MPA")
summary(scam33.mod.mpa6)
sm33.mpa6 <- as.data.frame(tidy(scam33.mod.mpa6))
sm33.mpa6$species <- "SCAM"
sm33.mpa6$temp <- "33"
sm33.mpa6$MPA <- "-0.6"

scam33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.33, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scam33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 33 @ -1.2 MPA")
summary(scam33.mod.mpa12)
sm33.mpa12 <- as.data.frame(tidy(scam33.mod.mpa12))
sm33.mpa12$species <- "SCAM"
sm33.mpa12$temp <- "33"
sm33.mpa12$MPA <- "-1.2"

### Temp 36 ----
scam.36 <- filter(scam.long, temp == "36")
scam36.mod.mpa0 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.36, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(scam36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 36 @ 0 MPA")
summary(scam36.mod.mpa0)
sm36.mpa0 <- as.data.frame(tidy(scam36.mod.mpa0))
sm36.mpa0$species <- "SCAM"
sm36.mpa0$temp <- "36"
sm36.mpa0$MPA <- "0"

scam36.mod.mpa15 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.36, type = "event", 
                        curveid = rep, subset = c(MPA == "-0.15"))
plot(scam36.mod.mpa15, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 36 @ -0.15 MPA")
summary(scam36.mod.mpa15)
sm36.mpa15 <- as.data.frame(tidy(scam36.mod.mpa15))
sm36.mpa15$species <- "SCAM"
sm36.mpa15$temp <- "36"
sm36.mpa15$MPA <- "-0.15"

scam36.mod.mpa3 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.3"))
plot(scam36.mod.mpa3, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 36 @ -0.3 MPA")
summary(scam36.mod.mpa3)
sm36.mpa3 <- as.data.frame(tidy(scam36.mod.mpa3))
sm36.mpa3$species <- "SCAM"
sm36.mpa3$temp <- "36"
sm36.mpa3$MPA <- "-0.3"

scam36.mod.mpa6 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = scam.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(scam36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 36 @ -0.6 MPA")
summary(scam36.mod.mpa6)
sm36.mpa6 <- as.data.frame(tidy(scam36.mod.mpa6))
sm36.mpa6$species <- "SCAM"
sm36.mpa6$temp <- "36"
sm36.mpa6$MPA <- "-0.6"

scam36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = scam.36, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(scam36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="SCAM 36 @ -1.2 MPA")
summary(scam36.mod.mpa12)
sm36.mpa12 <- as.data.frame(tidy(scam36.mod.mpa12))
sm36.mpa12$species <- "SCAM"
sm36.mpa12$temp <- "36"
sm36.mpa12$MPA <- "-1.2"

# Combine all data sets ----
germ.met.all <- rbind(ph23.mpa0, ph23.mpa15, ph23.mpa3, ph23.mpa6, ph23.mpa12,
                      ph28.mpa0, ph28.mpa15, ph28.mpa3, ph28.mpa6, ph28.mpa12,
                      ph33.mpa0, ph33.mpa15, ph33.mpa3, ph33.mpa6, ph33.mpa12,
                      ph36.mpa0, ph36.mpa15, ph36.mpa3, ph36.mpa6, ph36.mpa12,
                      s23.mpa0, s23.mpa15, s23.mpa3, s23.mpa6, s23.mpa12,
                      s28.mpa0, s28.mpa15, s28.mpa3, s28.mpa6, s28.mpa12,
                      s33.mpa0, s33.mpa15, s33.mpa3, s33.mpa6, s33.mpa12,
                      s36.mpa0, s36.mpa15, s36.mpa3, s36.mpa6, s36.mpa12,
                      sm23.mpa0, sm23.mpa15, sm23.mpa3, sm23.mpa6, sm23.mpa12,
                      sm28.mpa0, sm28.mpa15, sm28.mpa3, sm28.mpa6, sm28.mpa12,
                      sm33.mpa0, sm33.mpa15, sm33.mpa3, sm33.mpa6, sm33.mpa12,
                      sm36.mpa0, sm36.mpa15, sm36.mpa3, sm36.mpa6, sm36.mpa12)

# Remove negative for 'b' term (see reference below): disregard negative, comes from the drc package being rooted in bioassays
# https://www.statforbiology.com/2021/stat_drcte_2-methods/
for (i in 1:nrow(germ.met.all)){
  if (germ.met.all$term[i] == "b") (germ.met.all$estimate[i] <- abs(germ.met.all$estimate[i]))
}

# Define & relevel factors
str(germ.met.all)
germ.met.all$term <- as.factor(germ.met.all$term)
germ.met.all$curve <- as.factor(germ.met.all$curve)
germ.met.all$species <- as.factor(germ.met.all$species)
germ.met.all$temp <- as.factor(germ.met.all$temp)
germ.met.all$MPA <- as.factor(germ.met.all$MPA)
germ.met.all$MPA <- factor(germ.met.all$MPA, levels = c("0", "-0.15", "-0.3", "-0.6", "-1.2"))

ggplot(filter(germ.met.all, term == "d"), aes(x=estimate, y=species)) +
  geom_point() + xlim(0,1) # make sure all d fit between 0 & 1

# Calculate average values
size <- function(X) sum(!is.na(X))
gmet.stats <- germ.met.all %>%
  group_by(term, species, temp, MPA) %>%
  dplyr::summarize_at(vars(estimate), 
                      funs(mean, sd, size))
gmet.stats$SE <- gmet.stats$sd/sqrt(gmet.stats$size)

# Plot germination metrics ----
# b: slope around the inflection point
# d: max germination
# e: time to 50% germination (time to 50% of max germination)
## b plots; slope ----
b.g1 <- ggplot(filter(germ.met.all, term == "b"), aes(x=MPA, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha=0.2, position = position_jitter(width = 0.2, height = 0.1)) +
  facet_wrap(~species*temp); b.g1

# color pallete: ("cadetblue3", "goldenrod", "seagreen4", "orangered3")
b.g2 <- ggplot(filter(gmet.stats, term == "b"), aes(x=MPA, y=mean, fill = temp)) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width = 0.8, size=0.8, position=position_dodge(width=0.8)) +
  geom_point(aes(fill=temp), color="black", pch=21, size=3.5, position=position_dodge(width=0.8)) +
  scale_fill_grey(start=1, end=0) +
  facet_wrap(~species, ncol=1); b.g2

## d plots; slope ----
d.g1 <- ggplot(filter(germ.met.all, term == "d"), aes(x=MPA, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha=0.2, position = position_jitter(width = 0.2, height = 0.1)) +
  facet_wrap(~species*temp) + ylim(0,1); d.g1

# color pallete: ("cadetblue3", "goldenrod", "seagreen4", "orangered3")
d.g2 <- ggplot(filter(gmet.stats, term == "d"), aes(x=MPA, y=mean, fill = temp)) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width = 0.8, size=0.8, position=position_dodge(width=0.8)) +
  geom_point(aes(fill=temp), color="black", pch=21, size=3.5, position=position_dodge(width=0.8)) +
  scale_fill_grey(start=1, end=0) +
  ylim(0,1) + facet_wrap(~species, ncol=1); d.g2

## e plots; slope ----
e.g1 <- ggplot(filter(germ.met.all, term == "e"), aes(x=MPA, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha=0.2, position = position_jitter(width = 0.2, height = 0.1)) +
  facet_wrap(~species*temp); e.g1

# color pallete: ("cadetblue3", "goldenrod", "seagreen4", "orangered3")
e.g2 <- ggplot(filter(gmet.stats, term == "e"), aes(x=MPA, y=mean, fill = temp)) +
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width = 0.8, size=0.8, position=position_dodge(width=0.8)) +
  geom_point(aes(fill=temp), color="black", pch=21, size=3.5, position=position_dodge(width=0.8)) +
  scale_fill_grey(start=1, end=0) +
  facet_wrap(~species, ncol=1); e.g2

# Statistical models ----
## d; max germ ----
d.estimate <- filter(germ.met.all, term == "d")
hist(qlogis(d.estimate$estimate)) # this is of limited use, distribution of response conditional on predictors

for (i in 1:nrow(d.estimate)){
  if (d.estimate$estimate[i] == 0.00000000) d.estimate$estimate[i] <-  0.0005
}

# Add time component to dataset to incorporate time differences in model
d.estimate$time <- NA
for (i in 1:nrow(d.estimate)){
  if (d.estimate$MPA[i] == -1.2 | d.estimate$MPA[i] == -0.6 | d.estimate$curve[i] %in% c(4,5,6)){
    d.estimate$time[i] <-  2
  }else{
    d.estimate$time[i] <-  1
  }
}
d.estimate$time <- as.factor(d.estimate$time)
with(d.estimate,table(MPA,time))

# Look at MPA 0
d.estimate.MPA0 <- filter(d.estimate, MPA == 0)
time.mod <- lmer(qlogis(estimate) ~ (temp * species) + (time * species) +
                   (1|temp:time) + (1|temp:time:species), 
                 data=d.estimate.MPA0, control=lmerControl(optimizer="Nelder_Mead")) 
summary(time.mod)
Anova(time.mod)
emmeans(time.mod, ~ time)

ggplot(d.estimate.MPA0, aes(x=temp, y=estimate, color=time)) +
  geom_jitter(height=0, width=0.1, alpha=0.4) +  
  facet_wrap(~time + species)

ggplot(d.estimate.MPA0, aes(x=as.numeric(temp), y=estimate, color=time)) +
  geom_jitter(height=0, width=0.1, alpha=0.4) +  
  facet_wrap(~species) +
  geom_smooth()

time.mod2 <- lmer(qlogis(estimate) ~ (temp * species) +
                   (1|time) + (1|temp:time) + (1|temp:time:species), 
                 data=d.estimate.MPA0, control=lmerControl(optimizer="Nelder_Mead")) 
summary(time.mod2) # there is a time effect, estimates for MPA0 is higher at T1 relative to T2, we don't know if this applies to other MPAs
Anova(time.mod2)

# Run model!
d.mod1 <- lmer(qlogis(estimate) ~ temp * MPA * species +
               (1|time:temp) + (1|time:temp:MPA:species), 
               data=d.estimate, control=lmerControl(optimizer="Nelder_Mead")) 
# (1|temp:time) = chamber level variation within each time period;
# (1|time:temp:MPA) = variation among sets of boxes within chambers and time periods; 
summary(d.mod1)
Anova(d.mod1)
plot(d.mod1) # outlier.. 
# Try to find outlier
ggplot(d.estimate, aes(x=scale(as.numeric(MPA)), y=scale(qlogis(estimate)), color=species)) +
  geom_point() + facet_wrap(~as.numeric(temp)) # outlier - SCAM temp 23 where no germination occured - don't want to drop this outlier bcs it is an important data point

simulationOutput.d.mod1 <- DHARMa::simulateResiduals(d.mod1, plot=T) 
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$species) 
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$Temp) 

# Try beta distributed model
d.mod2 = glmmTMB::glmmTMB(estimate ~ temp * MPA * species + 
                         (1|time:temp) + (1|time:temp:MPA:species), 
                         data = d.estimate, beta_family(link="logit"))
summary(d.mod2)
Anova(d.mod2)
simulationOutput.d.mod2 <- DHARMa::simulateResiduals(d.mod2, plot=T) 
DHARMa::plotResiduals(simulationOutput.d.mod2, form = d.estimate$species) 
DHARMa::plotResiduals(simulationOutput.d.mod2, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod2, form = d.estimate$Temp)

##Pairwise comparisons for significant interactions (temp*MPA*species)
emmip(d.mod2, species ~ MPA|temp) 
emmip(d.mod2, species ~ temp|MPA)
emmip(d.mod2, temp ~ MPA|species)

d.emm1.1 <- emmeans(d.mod2, pairwise ~ MPA|temp|species, adjust='tukey', type = "response") 
d.emm1.2 <- emmeans(d.mod2, pairwise ~ temp|species|MPA, adjust='tukey', type = "response") 

d.mod1.rg <- update(ref_grid(d.mod1), tran = "logit") # back-transform from logit scale
emmeans(d.mod1.rg, ~ MPA + temp + species, type="response") # By default, 95% CIs
emm.d.mod1.rg <- confint(emmeans(d.mod1.rg, ~ MPA + temp + species, type="response"), adjust = "none", level = 0.68) # specify 68% CIs (roughly equivalent to +/- 1 SE)
emm.d.mod1.rg$adj <- (emm.d.mod1.rg$response)*100
emm.d.mod1.rg$adjLCL <- (emm.d.mod1.rg$lower.CL)*100
emm.d.mod1.rg$adjUCL <- (emm.d.mod1.rg$upper.CL)*100
emm.d.mod1.rg$adjSE <- (emm.d.mod1.rg$SE)*100
plot(emm.d.mod1.rg)

col.pal <- c("tomato3", "cadetblue", "goldenrod")
bw.pal <- c("gray80", "gray10", "gray40")
max.germ.plot <- ggplot(data=emm.d.mod1.rg, aes(x = MPA, y = adj, fill = species)) + 
  geom_errorbar(aes(ymin=adj-adjSE, ymax=adj+adjSE), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = species, color = species), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=species), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Water potential (MPA)") + facet_wrap(~temp) + 
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw()
max.germ.plot <- max.germ.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
max.germ.plot

max.germ.plot2 <- ggplot(data=emm.d.mod1.rg, aes(x = temp, y = adj, fill = species)) + 
  geom_errorbar(aes(ymin=adj-adjSE, ymax=adj+adjSE), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = species, color = species), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=species), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Temperature") + facet_wrap(~MPA) + 
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw()
max.germ.plot2 <- max.germ.plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
max.germ.plot2

col.pal2 <- c("cadetblue","seagreen4","goldenrod","tomato3")
max.germ.plot3 <- ggplot(data=emm.d.mod1.rg, aes(x = MPA, y = adj, fill = temp)) + 
  geom_errorbar(aes(ymin=adj-adjSE, ymax=adj+adjSE), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = temp, color = temp), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=temp), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Water potential") + facet_wrap(~species) + 
  scale_fill_manual(values=col.pal2, name="fill") +
  scale_color_manual(values=col.pal2) + 
  theme_bw()
max.germ.plot3 <- max.germ.plot3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
max.germ.plot3

# pull in observed data to overlay
d.estimate$adj <- d.estimate$estimate * 100

obs.data.pc <- ggplot(data=d.estimate, aes(x = MPA, y = adj, fill=temp)) + 
  geom_point(aes(fill=temp), alpha=0.1, pch = 21, size = 3.5, position = position_dodge(width = 0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Water potential (MPA)") + facet_wrap(~species) +
  scale_fill_manual(values=col.pal2) + theme_bw()
obs.data.pc <- obs.data.pc + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
obs.data.pc <- obs.data.pc + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none")
obs.data.pc <- obs.data.pc + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(), strip.text.x = element_text(size = 8))
obs.data.pc

# Merge graphs
g1 <- ggplotGrob(max.germ.plot3)
g2 <- ggplotGrob(obs.data.pc)

pp <- c(subset(g1$layout, grepl("panel", g1$layout$name), se = t:r))

g <- gtable_add_grob(g1, g2$grobs[grepl("panel", g1$layout$name)], 
                     pp$t, pp$l, pp$b, pp$l)

hinvert_title_grob <- function(grob){
  
  # Swap the widths
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  
  # Fix the justification
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}

# Get the y axis title from g2
index <- which(g2$layout$name == "ylab-l") # Which grob contains the y axis title?   EDIT HERE
ylab <- g2$grobs[[index]]                  # Extract that grob
ylab <- hinvert_title_grob(ylab)           # Swap margins and fix justifications

g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, ylab, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, clip = "off", name = "ylab-r")

index <- which(g2$layout$name == "axis-l-1-1")  # Which grob.    EDIT HERE
yaxis <- g2$grobs[[index]]                    # Extract the grob

ticks <- yaxis$children[[2]]
ticks$widths <- rev(ticks$widths)
ticks$grobs <- rev(ticks$grobs)

plot_theme <- function(p) {
  plyr::defaults(p$theme, theme_get())
}

tml <- plot_theme(g1)$axis.ticks.length   # Tick mark length
ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + tml

ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])

yaxis$children[[2]] <- ticks

g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, yaxis, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, 
                     clip = "off", name = "axis-r")

grid.draw(g)
#ggsave("./output/Exp1_byspp_maxgerm.png", g, width=8, height=5, unit="in", dpi=800)

## e; t50 ----
e.estimate <- filter(germ.met.all, term == "e")
hist(log(e.estimate$estimate))

ggplot(e.estimate, aes(x=species, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter() +
  facet_wrap (~temp*MPA)

for (i in 1:nrow(e.estimate)){
  if (e.estimate$estimate[i] == 0.00000000) e.estimate$estimate[i] <-  0.0005
}

# Add time component to dataset to incorporate time differences in model
e.estimate$time <- NA
for (i in 1:nrow(e.estimate)){
  if (e.estimate$MPA[i] == -1.2 | e.estimate$MPA[i] == -0.6){
    e.estimate$time[i] <-  2
  }else{
    e.estimate$time[i] <-  1
  }
}
e.estimate$time <- as.factor(e.estimate$time)

# Run model!
e.mod1 <- lmer(log(estimate) ~ temp * MPA * species +
                 (1|time:temp) + (1|time:temp:MPA:species), 
               data=e.estimate)  
summary(e.mod1) 
Anova(e.mod1)
plot(e.mod1)

simulationOutput.e.mod1 <- DHARMa::simulateResiduals(e.mod1, plot=T)
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate$species) 
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate$Temp) 

##Pairwise comparisons for significant interactions (temp*MPA*species)
emmip(e.mod1, species ~ MPA|temp) 
emmip(e.mod1, species ~ temp|MPA) 
emmip(e.mod1, temp ~ MPA|species)

emmeans(e.mod1, pairwise ~ MPA|temp|species, adjust='tukey') 

e.mod1.rg <- update(ref_grid(e.mod1), tran = "log") # back-transform from log scale
emmeans(e.mod1.rg, ~ MPA + temp + species, type="response") # By default, 95% CIs; what's up with these df??
emm.e.mod1.rg <- confint(emmeans(e.mod1.rg, ~ MPA + temp + species, type="response"), adjust = "none", level = 0.68) # specify 68% CIs (roughly equivalent to +/- 1 SE)

col.pal <- c("tomato3", "cadetblue", "goldenrod")
bw.pal <- c("gray80", "gray10", "gray40")
t50.plot <- ggplot(data=emm.e.mod1.rg, aes(x = MPA, y = response, fill = species)) + 
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = species, color = species), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=species), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Time to 50% germination (days)") + scale_y_continuous(breaks = seq(0,35,5), limits = c(0,35)) + 
  xlab("Water potential (MPA)") + facet_wrap(~temp) +
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw()
t50.plot <- t50.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
t50.plot

t50.plot2 <- ggplot(data=emm.e.mod1.rg, aes(x = temp, y = response, fill = species)) + 
  geom_errorbar(aes(ymin=response-SE, ymax=response+SE), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = species, color = species), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=species), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Time to 50% germination (days)") + scale_y_continuous(breaks = seq(0,35,5), limits = c(0,35)) + 
  xlab("Temperature") + facet_wrap(~MPA) +
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw()
t50.plot2 <- t50.plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
t50.plot2

# pull in observed data to overlay
obs.data.pc <- ggplot(data=e.estimate, aes(x = temp, y = estimate, fill=species)) + 
  geom_point(aes(fill=species), alpha=0.1, pch = 21, size = 3.5, position = position_dodge(width = 0.5)) +
  ylab("Time to 50% germination (days)") + scale_y_continuous(breaks = seq(0,35,5), limits = c(0,35)) + 
  xlab("Temperature") + facet_wrap(~MPA) +
  scale_fill_manual(values=col.pal) + theme_bw()
obs.data.pc <- obs.data.pc + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
obs.data.pc <- obs.data.pc + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none")
obs.data.pc <- obs.data.pc + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(), strip.text.x = element_text(size = 8))
obs.data.pc

# Merge graphs
g1 <- ggplotGrob(t50.plot2)
g2 <- ggplotGrob(obs.data.pc)

pp <- c(subset(g1$layout, grepl("panel", g1$layout$name), se = t:r))

g <- gtable_add_grob(g1, g2$grobs[grepl("panel", g1$layout$name)], 
                     pp$t, pp$l, pp$b, pp$l)

hinvert_title_grob <- function(grob){
  
  # Swap the widths
  widths <- grob$widths
  grob$widths[1] <- widths[3]
  grob$widths[3] <- widths[1]
  grob$vp[[1]]$layout$widths[1] <- widths[3]
  grob$vp[[1]]$layout$widths[3] <- widths[1]
  
  # Fix the justification
  grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
  grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
  grob$children[[1]]$x <- unit(1, "npc") - grob$children[[1]]$x
  grob
}

# Get the y axis title from g2
index <- which(g2$layout$name == "ylab-l") # Which grob contains the y axis title?   EDIT HERE
ylab <- g2$grobs[[index]]                  # Extract that grob
ylab <- hinvert_title_grob(ylab)           # Swap margins and fix justifications

g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, ylab, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, clip = "off", name = "ylab-r")

index <- which(g2$layout$name == "axis-l-1-1")  # Which grob.    EDIT HERE
yaxis <- g2$grobs[[index]]                    # Extract the grob

ticks <- yaxis$children[[2]]
ticks$widths <- rev(ticks$widths)
ticks$grobs <- rev(ticks$grobs)

plot_theme <- function(p) {
  plyr::defaults(p$theme, theme_get())
}

tml <- plot_theme(g1)$axis.ticks.length   # Tick mark length
ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, "npc") + tml

ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])

yaxis$children[[2]] <- ticks

g <- gtable_add_cols(g, g2$widths[g2$layout[index, ]$l], max(pp$r))
g <- gtable_add_grob(g, yaxis, max(pp$t), max(pp$r) + 1, max(pp$b), max(pp$r) + 1, 
                     clip = "off", name = "axis-r")

grid.draw(g)
#ggsave("./output/Exp1_tempx_t50.png", g, width=8, height=5, unit="in", dpi=800)
