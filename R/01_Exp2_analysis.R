# Experiment 2 ----
rm(list=ls())

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

# Load cleaned (non-cumulative) germination data ----
BL.exp2 <- read.csv("./data/exp2_BearLake_noncum.csv", as.is=T)
GSL.exp2 <- read.csv("./data/exp2_GSL_noncum.csv", as.is=T)
KCH.exp2 <- read.csv("./data/exp2_Kirch_noncum.csv", as.is=T)
NP.exp2 <- read.csv("./data/exp2_Ninepipe_noncum.csv", as.is=T)
PH.exp2 <- read.csv("./data/exp2_Pahran_noncum.csv", as.is=T)

# Merge & clean data sets ----
all.pop <- rbind(BL.exp2, GSL.exp2, KCH.exp2, NP.exp2, PH.exp2)
all.pop$temp <- as.factor(all.pop$temp)
all.pop$pop <- as.factor(all.pop$pop)
all.pop$MPA <- as.factor(all.pop$MPA)
all.pop$rep <- as.factor(all.pop$rep)

# Add day 0 column full of zero's to initialize germination & reorder dataframe
all.pop$Day.0 <- 0
all.pop <- all.pop[, c(1:5,26,6:25)]

# Remove "Day." for data collection columns
colnames(all.pop) <- c("temp", "pop", "MPA", "rep", "N",
                       "0","3","4","5","6","7","8","9","10","11",
                       "12","13","14","15","16","18","20","22","24",
                       "26", "28")

# Investigate values that are negative
all.pop$trt.id <- paste(all.pop$temp, all.pop$pop, all.pop$MPA, all.pop$rep)
length(unique(all.pop$trt.id))
dim(all.pop) # check that dimension and id length match

neg.val <- all.pop[rowSums(all.pop[6:26] < 0) > 0, ]
neg.val # 6 instances of negative values - check with Tatiana, input error?

# For now, assume input error and make negative values positive (but check with Tatiana)
all.pop[,6:26] <- abs(all.pop[,6:26])

# Remove trt.id column
all.pop <- all.pop[,1:26]

# Add replicate column that include all treatment factors
all.pop1 <- all.pop
all.pop1[["Rep"]] <- with(all.pop1, interaction(temp, pop, MPA, rep))

# Flip data set
pop.long <- pivot_longer(all.pop1, cols="0":"28", names_to = "Day", values_to = "Germ")

# Add time interval
pop.long$Start <- pop.long$Day
pop.long$End <- NA

for (i in 1:nrow(pop.long)){
  if(pop.long$Day[i]=="28"){
    pop.long$End[i] <- "Inf"
  }else{
    pop.long$End[i] <- pop.long$Day[i+1]
  }#ifelse
}#for

pop.long$Start <- as.numeric(pop.long$Start)
pop.long$End <- as.numeric(pop.long$End)

# Calculate cumulative germination
pop.long$CumGerm <- NA
for (i in 1:nrow(pop.long)){
  if(pop.long$Day[i]=="0"){
    pop.long$CumGerm[i] <- pop.long$Germ[i]
  }else{
    pop.long$CumGerm[i] <- pop.long$Germ[i] + pop.long$CumGerm[i-1]
  }#ifelse
}#for

# Calculate cumulative germination proportion
pop.long$CumGRP <- pop.long$CumGerm / pop.long$N

# Turn day into factor and level correctly
pop.long$Day <- as.factor(pop.long$Day)
day.vec <- c(0,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28)
pop.long$Day <- factor(pop.long$Day, levels = day.vec)

# Assign remaining ungerminated seeds to "Inf" time interval
for (i in 1:nrow(pop.long)){
  if (pop.long$End[i] == "Inf"){
    pop.long$Germ[i] <- pop.long$N[i] - pop.long$CumGerm[i]
  } 
}

# Plot germination curves - clean this up! ----
ggplot(pop.long, aes(x = Day, y = CumGRP, color=rep, shape=MPA)) + 
  geom_point(size=2.5) + ylim(0,1) +
  geom_line() + facet_wrap(~pop*temp) + labs(x ="Days", y = "Germination proportion") +
  scale_colour_grey(start = 0, end = .6) + scale_shape_manual(values = c(0, 16, 4, 6, 15)) +
  theme_bw()

germ.stats <- pop.long %>%
  group_by(temp, pop, MPA, Day) %>% 
  dplyr::summarize(Mean  = mean(CumGRP, na.rm = TRUE),
                   SD    = sd(CumGRP, na.rm = TRUE),
                   nsize = sum(!is.na(CumGRP)),
                   SE    = SD/sqrt(nsize))

ggplot(germ.stats, aes(x = Day, y = Mean, group=MPA, color=MPA, shape=MPA)) + 
  geom_line(aes(group=MPA, color=MPA)) + 
  geom_point(size=2.5) + ylim(0,1) +
  facet_wrap(~pop*temp) + labs(x ="Days", y = "Germination proportion") +
  scale_colour_grey(start = 0, end = 0.5) + scale_shape_manual(values = c(0, 16, 4, 6, 15)) +
  theme_bw()

ggplot(germ.stats, aes(x = Day, y = Mean, group=MPA, color=MPA, shape=MPA)) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width = 0.2) +
  geom_line(aes(group=MPA, color=MPA)) + 
  geom_point(size=2.5) + ylim(0,1) +
  facet_wrap(~pop*temp) + labs(x ="Days", y = "Germination proportion") +
  scale_colour_grey(start = 0, end = 0.5) + scale_shape_manual(values = c(0, 16, 4, 6, 15)) +
  theme_bw()

# Model germination curves ----
## Split long data set
BL.long <- filter(pop.long, pop == "Bear Lake")
GSL.long <- filter(pop.long, pop == "GSL")
K.long <- filter(pop.long, pop == "Kirch")
NP.long <- filter(pop.long, pop == "Ninepipe")
P.long <- filter(pop.long, pop == "Pahranagat")

## Bear Lake models ----
### Temp 23 ----
BL.23 <- filter(BL.long, temp == "23")
BL.23.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                       data = BL.23, type = "event", 
                       curveid = rep, subset = c(MPA == "0"))
plot(BL.23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 23 @ 0 MPA")
summary(BL.23.mod.mpa0)
BL23.mpa0 <- as.data.frame(tidy(BL.23.mod.mpa0))
BL23.mpa0$pop <- "Bear Lake"
BL23.mpa0$temp <- "23"
BL23.mpa0$MPA <- "0"

BL.23.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                       data = BL.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-0.6"))
plot(BL.23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 23 @ -0.6 MPA")
summary(BL.23.mod.mpa6) 
BL23.mpa6 <- as.data.frame(tidy(BL.23.mod.mpa6))
BL23.mpa6$pop <- "Bear Lake"
BL23.mpa6$temp <- "23"
BL23.mpa6$MPA <- "-0.6"

BL.23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = BL.23, type = "event", 
                        curveid = rep, subset = c(MPA == "-1.2"))
plot(BL.23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 23 @ -1.2 MPA")
summary(BL.23.mod.mpa12)
BL23.mpa12 <- as.data.frame(tidy(BL.23.mod.mpa12))
BL23.mpa12$pop <- "Bear Lake"
BL23.mpa12$temp <- "23"
BL23.mpa12$MPA <- "-1.2"

### Temp 28 ----
BL.28 <- filter(BL.long, temp == "28")
BL.28.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = BL.28, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(BL.28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 28 @ 0 MPA")
summary(BL.28.mod.mpa0)
BL28.mpa0 <- as.data.frame(tidy(BL.28.mod.mpa0))
BL28.mpa0$pop <- "Bear Lake"
BL28.mpa0$temp <- "28"
BL28.mpa0$MPA <- "0"

BL.28.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = BL.28, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(BL.28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 28 @ -0.6 MPA")
summary(BL.28.mod.mpa6) 
BL28.mpa6 <- as.data.frame(tidy(BL.28.mod.mpa6))
BL28.mpa6$pop <- "Bear Lake"
BL28.mpa6$temp <- "28"
BL28.mpa6$MPA <- "-0.6"

BL.28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = BL.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(BL.28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 28 @ -1.2 MPA")
summary(BL.28.mod.mpa12)
BL28.mpa12 <- as.data.frame(tidy(BL.28.mod.mpa12))
BL28.mpa12$pop <- "Bear Lake"
BL28.mpa12$temp <- "28"
BL28.mpa12$MPA <- "-1.2"

### Temp 33 ----
BL.33 <- filter(BL.long, temp == "33")
BL.33.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = BL.33, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(BL.33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 33 @ 0 MPA")
summary(BL.33.mod.mpa0)
BL33.mpa0 <- as.data.frame(tidy(BL.33.mod.mpa0))
BL33.mpa0$pop <- "Bear Lake"
BL33.mpa0$temp <- "33"
BL33.mpa0$MPA <- "0"

BL.33.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = BL.33, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(BL.33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 33 @ -0.6 MPA")
summary(BL.33.mod.mpa6) 
BL33.mpa6 <- as.data.frame(tidy(BL.33.mod.mpa6))
BL33.mpa6$pop <- "Bear Lake"
BL33.mpa6$temp <- "33"
BL33.mpa6$MPA <- "-0.6"

BL.33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = BL.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(BL.33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 33 @ -1.2 MPA")
summary(BL.33.mod.mpa12)
BL33.mpa12 <- as.data.frame(tidy(BL.33.mod.mpa12))
BL33.mpa12$pop <- "Bear Lake"
BL33.mpa12$temp <- "33"
BL33.mpa12$MPA <- "-1.2"

### Temp 36 ----
BL.36 <- filter(BL.long, temp == "36")
BL.36.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = BL.36, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(BL.36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 36 @ 0 MPA")
summary(BL.36.mod.mpa0)
BL36.mpa0 <- as.data.frame(tidy(BL.36.mod.mpa0))
BL36.mpa0$pop <- "Bear Lake"
BL36.mpa0$temp <- "36"
BL36.mpa0$MPA <- "0"

BL.36.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = BL.36, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(BL.36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 36 @ -0.6 MPA")
summary(BL.36.mod.mpa6) 
BL36.mpa6 <- as.data.frame(tidy(BL.36.mod.mpa6))
BL36.mpa6$pop <- "Bear Lake"
BL36.mpa6$temp <- "36"
BL36.mpa6$MPA <- "-0.6"

BL.36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = BL.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(BL.36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="BL 36 @ -1.2 MPA")
summary(BL.36.mod.mpa12)
BL36.mpa12 <- as.data.frame(tidy(BL.36.mod.mpa12))
BL36.mpa12$pop <- "Bear Lake"
BL36.mpa12$temp <- "36"
BL36.mpa12$MPA <- "-1.2"

## GSL models ----
### Temp 23 ----
GSL.23 <- filter(GSL.long, temp == "23")
GSL.23.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = GSL.23, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(GSL.23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 23 @ 0 MPA")
summary(GSL.23.mod.mpa0)
GSL23.mpa0 <- as.data.frame(tidy(GSL.23.mod.mpa0))
GSL23.mpa0$pop <- "GSL"
GSL23.mpa0$temp <- "23"
GSL23.mpa0$MPA <- "0"

GSL.23.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = GSL.23, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(GSL.23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 23 @ -0.6 MPA")
summary(GSL.23.mod.mpa6) 
GSL23.mpa6 <- as.data.frame(tidy(GSL.23.mod.mpa6))
GSL23.mpa6$pop <- "GSL"
GSL23.mpa6$temp <- "23"
GSL23.mpa6$MPA <- "-0.6"

GSL.23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = GSL.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(GSL.23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 23 @ -1.2 MPA")
summary(GSL.23.mod.mpa12) # rep 1 max germ estimation is WAY high; can't resolve with weibull; drop rep

GSL.23.mpa12.norep1 <- filter(GSL.23, MPA == "-1.2" & rep != "1")
GSL.23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                        data = GSL.23.mpa12.norep1, type = "event", 
                        curveid = rep)
plot(GSL.23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 23 @ -1.2 MPA")
summary(GSL.23.mod.mpa12)
GSL23.mpa12 <- as.data.frame(tidy(GSL.23.mod.mpa12))
GSL23.mpa12$pop <- "GSL"
GSL23.mpa12$temp <- "23"
GSL23.mpa12$MPA <- "-1.2"

### Temp 28 ----
GSL.28 <- filter(GSL.long, temp == "28")
GSL.28.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = GSL.28, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(GSL.28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 28 @ 0 MPA")
summary(GSL.28.mod.mpa0)
GSL28.mpa0 <- as.data.frame(tidy(GSL.28.mod.mpa0))
GSL28.mpa0$pop <- "GSL"
GSL28.mpa0$temp <- "28"
GSL28.mpa0$MPA <- "0"

GSL.28.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = GSL.28, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(GSL.28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 28 @ -0.6 MPA")
summary(GSL.28.mod.mpa6) 
GSL28.mpa6 <- as.data.frame(tidy(GSL.28.mod.mpa6))
GSL28.mpa6$pop <- "GSL"
GSL28.mpa6$temp <- "28"
GSL28.mpa6$MPA <- "-0.6"

GSL.28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = GSL.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(GSL.28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 28 @ -1.2 MPA")
summary(GSL.28.mod.mpa12)
GSL28.mpa12 <- as.data.frame(tidy(GSL.28.mod.mpa12))
GSL28.mpa12$pop <- "GSL"
GSL28.mpa12$temp <- "28"
GSL28.mpa12$MPA <- "-1.2"

### Temp 33 ----
GSL.33 <- filter(GSL.long, temp == "33")
GSL.33.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = GSL.33, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(GSL.33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 33 @ 0 MPA")
summary(GSL.33.mod.mpa0)
GSL33.mpa0 <- as.data.frame(tidy(GSL.33.mod.mpa0))
GSL33.mpa0$pop <- "GSL"
GSL33.mpa0$temp <- "33"
GSL33.mpa0$MPA <- "0"

GSL.33.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = GSL.33, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(GSL.33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 33 @ -0.6 MPA")
summary(GSL.33.mod.mpa6) 
GSL33.mpa6 <- as.data.frame(tidy(GSL.33.mod.mpa6))
GSL33.mpa6$pop <- "GSL"
GSL33.mpa6$temp <- "33"
GSL33.mpa6$MPA <- "-0.6"

GSL.33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = GSL.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(GSL.33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 33 @ -1.2 MPA")
summary(GSL.33.mod.mpa12)
GSL33.mpa12 <- as.data.frame(tidy(GSL.33.mod.mpa12))
GSL33.mpa12$pop <- "GSL"
GSL33.mpa12$temp <- "33"
GSL33.mpa12$MPA <- "-1.2"

### Temp 36 ----
GSL.36 <- filter(GSL.long, temp == "36")
GSL.36.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = GSL.36, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(GSL.36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 36 @ 0 MPA")
summary(GSL.36.mod.mpa0)
GSL36.mpa0 <- as.data.frame(tidy(GSL.36.mod.mpa0))
GSL36.mpa0$pop <- "GSL"
GSL36.mpa0$temp <- "36"
GSL36.mpa0$MPA <- "0"

GSL.36.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = GSL.36, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(GSL.36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 36 @ -0.6 MPA")
summary(GSL.36.mod.mpa6) 
GSL36.mpa6 <- as.data.frame(tidy(GSL.36.mod.mpa6))
GSL36.mpa6$pop <- "GSL"
GSL36.mpa6$temp <- "36"
GSL36.mpa6$MPA <- "-0.6"

GSL.36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = GSL.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(GSL.36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="GSL 36 @ -1.2 MPA")
summary(GSL.36.mod.mpa12)
GSL36.mpa12 <- as.data.frame(tidy(GSL.36.mod.mpa12))
GSL36.mpa12$pop <- "GSL"
GSL36.mpa12$temp <- "36"
GSL36.mpa12$MPA <- "-1.2"

## Kirch models ----
### Temp 23 ----
K.23 <- filter(K.long, temp == "23")
K.23.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = K.23, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(K.23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 23 @ 0 MPA")
summary(K.23.mod.mpa0)
K23.mpa0 <- as.data.frame(tidy(K.23.mod.mpa0))
K23.mpa0$pop <- "Kirch"
K23.mpa0$temp <- "23"
K23.mpa0$MPA <- "0"

K.23.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = K.23, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(K.23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 23 @ -0.6 MPA")
summary(K.23.mod.mpa6) 
K23.mpa6 <- as.data.frame(tidy(K.23.mod.mpa6))
K23.mpa6$pop <- "Kirch"
K23.mpa6$temp <- "23"
K23.mpa6$MPA <- "-0.6"

K.23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = K.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(K.23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 23 @ -1.2 MPA")
summary(K.23.mod.mpa12)
K23.mpa12 <- as.data.frame(tidy(K.23.mod.mpa12))
K23.mpa12$pop <- "Kirch"
K23.mpa12$temp <- "23"
K23.mpa12$MPA <- "-1.2"

### Temp 28 ----
K.28 <- filter(K.long, temp == "28")
K.28.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = K.28, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(K.28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 28 @ 0 MPA")
summary(K.28.mod.mpa0)
K28.mpa0 <- as.data.frame(tidy(K.28.mod.mpa0))
K28.mpa0$pop <- "Kirch"
K28.mpa0$temp <- "28"
K28.mpa0$MPA <- "0"

K.28.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = K.28, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(K.28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 28 @ -0.6 MPA")
summary(K.28.mod.mpa6) 
K28.mpa6 <- as.data.frame(tidy(K.28.mod.mpa6))
K28.mpa6$pop <- "Kirch"
K28.mpa6$temp <- "28"
K28.mpa6$MPA <- "-0.6"

K.28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = K.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(K.28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 28 @ -1.2 MPA")
summary(K.28.mod.mpa12)
K28.mpa12 <- as.data.frame(tidy(K.28.mod.mpa12))
K28.mpa12$pop <- "Kirch"
K28.mpa12$temp <- "28"
K28.mpa12$MPA <- "-1.2"

### Temp 33 ----
K.33 <- filter(K.long, temp == "33")
K.33.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = K.33, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(K.33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 33 @ 0 MPA")
summary(K.33.mod.mpa0)
K33.mpa0 <- as.data.frame(tidy(K.33.mod.mpa0))
K33.mpa0$pop <- "Kirch"
K33.mpa0$temp <- "33"
K33.mpa0$MPA <- "0"

K.33.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = K.33, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(K.33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 33 @ -0.6 MPA")
summary(K.33.mod.mpa6) 
K33.mpa6 <- as.data.frame(tidy(K.33.mod.mpa6))
K33.mpa6$pop <- "Kirch"
K33.mpa6$temp <- "33"
K33.mpa6$MPA <- "-0.6"

K.33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = K.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(K.33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 33 @ -1.2 MPA")
summary(K.33.mod.mpa12)
K33.mpa12 <- as.data.frame(tidy(K.33.mod.mpa12))
K33.mpa12$pop <- "Kirch"
K33.mpa12$temp <- "33"
K33.mpa12$MPA <- "-1.2"

### Temp 36 ----
K.36 <- filter(K.long, temp == "36")
K.36.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = K.36, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(K.36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 36 @ 0 MPA")
summary(K.36.mod.mpa0)
K36.mpa0 <- as.data.frame(tidy(K.36.mod.mpa0))
K36.mpa0$pop <- "Kirch"
K36.mpa0$temp <- "36"
K36.mpa0$MPA <- "0"

K.36.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = K.36, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(K.36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 36 @ -0.6 MPA")
summary(K.36.mod.mpa6) 
K36.mpa6 <- as.data.frame(tidy(K.36.mod.mpa6))
K36.mpa6$pop <- "Kirch"
K36.mpa6$temp <- "36"
K36.mpa6$MPA <- "-0.6"

K.36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = K.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(K.36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="K 36 @ -1.2 MPA")
summary(K.36.mod.mpa12)
K36.mpa12 <- as.data.frame(tidy(K.36.mod.mpa12))
K36.mpa12$pop <- "Kirch"
K36.mpa12$temp <- "36"
K36.mpa12$MPA <- "-1.2"

## Ninepipe models ----
### Temp 23 ----
NP.23 <- filter(NP.long, temp == "23")
NP.23.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = NP.23, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(NP.23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 23 @ 0 MPA")
summary(NP.23.mod.mpa0)
NP23.mpa0 <- as.data.frame(tidy(NP.23.mod.mpa0))
NP23.mpa0$pop <- "Ninepipe"
NP23.mpa0$temp <- "23"
NP23.mpa0$MPA <- "0"

NP.23.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = NP.23, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(NP.23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 23 @ -0.6 MPA")
summary(NP.23.mod.mpa6) 
NP23.mpa6 <- as.data.frame(tidy(NP.23.mod.mpa6))
NP23.mpa6$pop <- "Ninepipe"
NP23.mpa6$temp <- "23"
NP23.mpa6$MPA <- "-0.6"

NP.23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = NP.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(NP.23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 23 @ -1.2 MPA")
summary(NP.23.mod.mpa12)
NP23.mpa12 <- as.data.frame(tidy(NP.23.mod.mpa12))
NP23.mpa12$pop <- "Ninepipe"
NP23.mpa12$temp <- "23"
NP23.mpa12$MPA <- "-1.2"

### Temp 28 ----
NP.28 <- filter(NP.long, temp == "28")
NP.28.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = NP.28, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(NP.28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 28 @ 0 MPA")
summary(NP.28.mod.mpa0)
NP28.mpa0 <- as.data.frame(tidy(NP.28.mod.mpa0))
NP28.mpa0$pop <- "Ninepipe"
NP28.mpa0$temp <- "28"
NP28.mpa0$MPA <- "0"

NP.28.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = NP.28, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(NP.28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 28 @ -0.6 MPA")
summary(NP.28.mod.mpa6) 
NP28.mpa6 <- as.data.frame(tidy(NP.28.mod.mpa6))
NP28.mpa6$pop <- "Ninepipe"
NP28.mpa6$temp <- "28"
NP28.mpa6$MPA <- "-0.6"

NP.28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = NP.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(NP.28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 28 @ -1.2 MPA")
summary(NP.28.mod.mpa12)
NP28.mpa12 <- as.data.frame(tidy(NP.28.mod.mpa12))
NP28.mpa12$pop <- "Ninepipe"
NP28.mpa12$temp <- "28"
NP28.mpa12$MPA <- "-1.2"

### Temp 33 ----
NP.33 <- filter(NP.long, temp == "33")
NP.33.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = NP.33, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(NP.33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 33 @ 0 MPA")
summary(NP.33.mod.mpa0)
NP33.mpa0 <- as.data.frame(tidy(NP.33.mod.mpa0))
NP33.mpa0$pop <- "Ninepipe"
NP33.mpa0$temp <- "33"
NP33.mpa0$MPA <- "0"

NP.33.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = NP.33, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(NP.33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 33 @ -0.6 MPA")
summary(NP.33.mod.mpa6) 
NP33.mpa6 <- as.data.frame(tidy(NP.33.mod.mpa6))
NP33.mpa6$pop <- "Ninepipe"
NP33.mpa6$temp <- "33"
NP33.mpa6$MPA <- "-0.6"

NP.33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = NP.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(NP.33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 33 @ -1.2 MPA")
summary(NP.33.mod.mpa12)
NP33.mpa12 <- as.data.frame(tidy(NP.33.mod.mpa12))
NP33.mpa12$pop <- "Ninepipe"
NP33.mpa12$temp <- "33"
NP33.mpa12$MPA <- "-1.2"

### Temp 36 ----
NP.36 <- filter(NP.long, temp == "36")
NP.36.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = NP.36, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(NP.36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 36 @ 0 MPA")
summary(NP.36.mod.mpa0)
NP36.mpa0 <- as.data.frame(tidy(NP.36.mod.mpa0))
NP36.mpa0$pop <- "Ninepipe"
NP36.mpa0$temp <- "36"
NP36.mpa0$MPA <- "0"

NP.36.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = NP.36, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(NP.36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 36 @ -0.6 MPA")
summary(NP.36.mod.mpa6) 
NP36.mpa6 <- as.data.frame(tidy(NP.36.mod.mpa6))
NP36.mpa6$pop <- "Ninepipe"
NP36.mpa6$temp <- "36"
NP36.mpa6$MPA <- "-0.6"

NP.36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = NP.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(NP.36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="NP 36 @ -1.2 MPA")
summary(NP.36.mod.mpa12)
NP36.mpa12 <- as.data.frame(tidy(NP.36.mod.mpa12))
NP36.mpa12$pop <- "Ninepipe"
NP36.mpa12$temp <- "36"
NP36.mpa12$MPA <- "-1.2"

## Pahranagat models ----
### Temp 23 ----
P.23 <- filter(P.long, temp == "23")
P.23.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = P.23, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(P.23.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 23 @ 0 MPA")
summary(P.23.mod.mpa0)
P23.mpa0 <- as.data.frame(tidy(P.23.mod.mpa0))
P23.mpa0$pop <- "Pahranagat"
P23.mpa0$temp <- "23"
P23.mpa0$MPA <- "0"

P.23.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = P.23, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(P.23.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 23 @ -0.6 MPA")
summary(P.23.mod.mpa6) 
P23.mpa6 <- as.data.frame(tidy(P.23.mod.mpa6))
P23.mpa6$pop <- "Pahranagat"
P23.mpa6$temp <- "23"
P23.mpa6$MPA <- "-0.6"

P.23.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = P.23, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(P.23.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 23 @ -1.2 MPA")
summary(P.23.mod.mpa12)
P23.mpa12 <- as.data.frame(tidy(P.23.mod.mpa12))
P23.mpa12$pop <- "Pahranagat"
P23.mpa12$temp <- "23"
P23.mpa12$MPA <- "-1.2"

### Temp 28 ----
P.28 <- filter(P.long, temp == "28")
P.28.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = P.28, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(P.28.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 28 @ 0 MPA")
summary(P.28.mod.mpa0)
P28.mpa0 <- as.data.frame(tidy(P.28.mod.mpa0))
P28.mpa0$pop <- "Pahranagat"
P28.mpa0$temp <- "28"
P28.mpa0$MPA <- "0"

P.28.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = P.28, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(P.28.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 28 @ -0.6 MPA")
summary(P.28.mod.mpa6) 
P28.mpa6 <- as.data.frame(tidy(P.28.mod.mpa6))
P28.mpa6$pop <- "Pahranagat"
P28.mpa6$temp <- "28"
P28.mpa6$MPA <- "-0.6"

P.28.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = P.28, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(P.28.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 28 @ -1.2 MPA")
summary(P.28.mod.mpa12)
P28.mpa12 <- as.data.frame(tidy(P.28.mod.mpa12))
P28.mpa12$pop <- "Pahranagat"
P28.mpa12$temp <- "28"
P28.mpa12$MPA <- "-1.2"

### Temp 33 ----
P.33 <- filter(P.long, temp == "33")
P.33.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = P.33, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(P.33.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 33 @ 0 MPA")
summary(P.33.mod.mpa0)
P33.mpa0 <- as.data.frame(tidy(P.33.mod.mpa0))
P33.mpa0$pop <- "Pahranagat"
P33.mpa0$temp <- "33"
P33.mpa0$MPA <- "0"

P.33.mod.mpa6 <- drm(Germ ~ Start + End, fct=W2.4(),
                      data = P.33, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(P.33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 33 @ -0.6 MPA")
summary(P.33.mod.mpa6) 

# rep 3 total germ is higher than N, causing issues with model. Remove rep 3
P.33.6 <- filter(P.33, MPA == "-0.6" & rep != "3")
P.33.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(), 
                     data = P.33.6, type = "event", 
                     curveid = rep)
plot(P.33.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 33 @ -0.6 MPA")
summary(P.33.mod.mpa6) 
P33.mpa6 <- as.data.frame(tidy(P.33.mod.mpa6))
P33.mpa6$pop <- "Pahranagat"
P33.mpa6$temp <- "33"
P33.mpa6$MPA <- "-0.6"

P.33.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = P.33, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(P.33.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 33 @ -1.2 MPA")
summary(P.33.mod.mpa12)
P33.mpa12 <- as.data.frame(tidy(P.33.mod.mpa12))
P33.mpa12$pop <- "Pahranagat"
P33.mpa12$temp <- "33"
P33.mpa12$MPA <- "-1.2"

### Temp 36 ----
P.36 <- filter(P.long, temp == "36")
P.36.mod.mpa0 <- drm(Germ ~ Start + End, fct=LL.3(), 
                      data = P.36, type = "event", 
                      curveid = rep, subset = c(MPA == "0"))
plot(P.36.mod.mpa0, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 36 @ 0 MPA")
summary(P.36.mod.mpa0)
P36.mpa0 <- as.data.frame(tidy(P.36.mod.mpa0))
P36.mpa0$pop <- "Pahranagat"
P36.mpa0$temp <- "36"
P36.mpa0$MPA <- "0"

P.36.mod.mpa6 <- drm(Germ ~ Start + End, fct=LL.3(),
                      data = P.36, type = "event", 
                      curveid = rep, subset = c(MPA == "-0.6"))
plot(P.36.mod.mpa6, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 36 @ -0.6 MPA")
summary(P.36.mod.mpa6) 
P36.mpa6 <- as.data.frame(tidy(P.36.mod.mpa6))
P36.mpa6$pop <- "Pahranagat"
P36.mpa6$temp <- "36"
P36.mpa6$MPA <- "-0.6"

P.36.mod.mpa12 <- drm(Germ ~ Start + End, fct = LL.3(),
                       data = P.36, type = "event", 
                       curveid = rep, subset = c(MPA == "-1.2"))
plot(P.36.mod.mpa12, log = "", xlab = "Time", 
     ylab = "Proportion of germinated seeds",
     xlim = c(0, 30), ylim = c(0,1), main="P 36 @ -1.2 MPA")
summary(P.36.mod.mpa12)
P36.mpa12 <- as.data.frame(tidy(P.36.mod.mpa12))
P36.mpa12$pop <- "Pahranagat"
P36.mpa12$temp <- "36"
P36.mpa12$MPA <- "-1.2"

# Combine all data sets ----
germ.met.all <- rbind(BL23.mpa0, BL23.mpa6, BL23.mpa12, BL28.mpa0, BL28.mpa6, BL28.mpa12, 
                      BL33.mpa0, BL33.mpa6, BL33.mpa12, BL36.mpa0, BL36.mpa6, BL36.mpa12, 
                      GSL23.mpa0, GSL23.mpa6, GSL23.mpa12, GSL28.mpa0, GSL28.mpa6, GSL28.mpa12, 
                      GSL33.mpa0, GSL33.mpa6, GSL33.mpa12, GSL36.mpa0, GSL36.mpa6, GSL36.mpa12, 
                      K23.mpa0, K23.mpa6, K23.mpa12, K28.mpa0, K28.mpa6, K28.mpa12, 
                      K33.mpa0, K33.mpa6, K33.mpa12, K36.mpa0, K36.mpa6, K36.mpa12, 
                      NP23.mpa0, NP23.mpa6, NP23.mpa12, NP28.mpa0, NP28.mpa6, NP28.mpa12, 
                      NP33.mpa0, NP33.mpa6, NP33.mpa12, NP36.mpa0, NP36.mpa6, NP36.mpa12, 
                      P23.mpa0, P23.mpa6, P23.mpa12, P28.mpa0, P28.mpa6, P28.mpa12, 
                      P33.mpa0, P33.mpa6, P33.mpa12, P36.mpa0, P36.mpa6, P36.mpa12) 
                      
                      

# Remove negative for 'b' term (see reference below): disregard negative, comes from the drc package being rooted in bioassays
# https://www.statforbiology.com/2021/stat_drcte_2-methods/
for (i in 1:nrow(germ.met.all)){
  if (germ.met.all$term[i] == "b") (germ.met.all$estimate[i] <- abs(germ.met.all$estimate[i]))
}

# Define & relevel factors
str(germ.met.all)
germ.met.all$term <- as.factor(germ.met.all$term)
germ.met.all$curve <- as.factor(germ.met.all$curve)
germ.met.all$pop <- as.factor(germ.met.all$pop)
germ.met.all$temp <- as.factor(germ.met.all$temp)
germ.met.all$MPA <- as.factor(germ.met.all$MPA)
germ.met.all$MPA <- factor(germ.met.all$MPA, levels = c("0", "-0.6", "-1.2"))

ggplot(filter(germ.met.all, term == "d"), aes(x=estimate, y=pop)) +
  geom_point() + xlim(0,1)

# Calculate average values
size <- function(X) sum(!is.na(X))
gmet.stats <- germ.met.all %>%
  group_by(term, pop, temp, MPA) %>%
  dplyr::summarize_at(vars(estimate), 
                      funs(mean, sd, size))
gmet.stats$SE <- gmet.stats$sd/sqrt(gmet.stats$size)

# Plot germination metrics ----
# d: max germination
# e: time to 50% germination (time to 50% of max germination)
## d plots; max germ ----
d.g1 <- ggplot(filter(germ.met.all, term == "d"), aes(x=MPA, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha=0.2, position = position_jitter(width = 0.2, height = 0.1)) +
  facet_wrap(~pop*temp) + ylim(0,1); d.g1

## e plots; t50 ----
e.g1 <- ggplot(filter(germ.met.all, term == "e"), aes(x=MPA, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(alpha=0.2, position = position_jitter(width = 0.2, height = 0.1)) +
  facet_wrap(~pop*temp); e.g1

# Statistical models ----
## d; max germ ----
d.estimate <- filter(germ.met.all, term == "d")
hist(qlogis(d.estimate$estimate))

for (i in 1:nrow(d.estimate)){
  if (d.estimate$estimate[i] == 0.00000000) d.estimate$estimate[i] <-  0.0005
}

# Run model!
d.mod1 <- lmer(qlogis(estimate) ~ temp * MPA +
                 (1|pop), 
               data=d.estimate, control=lmerControl(optimizer="bobyqa"))

summary(d.mod1) 
Anova(d.mod1)
plot(d.mod1) # not great, but better than other transformations (log & sqrt looks worse)
simulationOutput.d.mod1 <- DHARMa::simulateResiduals(d.mod1, plot=T) # ehh, not great; better than sqrt
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$pop) # wow!
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$temp) 
emmip(d.mod1, temp ~ MPA)
emmip(d.mod1, MPA ~ temp)

# Test the significance (and visualize) source
d.mod2 <- lmer(qlogis(estimate) ~ temp * MPA +
                 (1|temp) + (1|temp:MPA:pop), 
               data=d.estimate, control=lmerControl(optimizer="Nelder_Mead"))
lmtest::lrtest(d.mod1, d.mod2) # source is significant
sjPlot::plot_model(d.mod1, type="re", show.p=T) # see plot [[2]]
ranef(d.mod2) 
# https://cran.microsoft.com/snapshot/2017-04-21/web/packages/sjPlot/vignettes/sjplmer.html

####
# Try model with populations as fixed effect
d.mod3 <- lm(qlogis(estimate) ~ temp * MPA * pop,
               data=d.estimate) # look at 3-way here
# (1|temp) = chamber level variation; (1|temp:MPA) = variation among sets of cups across chambers; 
# (1|temp:MPA:pop) = variation among sets of cups within chambers; 
# residual = variation in germination within cup (i.e., response)
summary(d.mod3) 
Anova(d.mod3)
plot(d.mod3) 
emmip(d.mod3, MPA ~ temp|pop) 

simulationOutput.d.mod3 <- DHARMa::simulateResiduals(d.mod3, plot=T)
DHARMa::plotResiduals(simulationOutput.d.mod3, form = d.estimate$pop) 
DHARMa::plotResiduals(simulationOutput.d.mod3, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod3, form = d.estimate$Temp)

# Beta distributed model
# Try beta distributed model
d.mod2 = glmmTMB::glmmTMB(estimate ~ temp * MPA * pop, 
                          data = d.estimate, beta_family(link="logit"))
summary(d.mod2)
Anova(d.mod2)
simulationOutput.d.mod2 <- DHARMa::simulateResiduals(d.mod2, plot=T) 
DHARMa::plotResiduals(simulationOutput.d.mod2, form = d.estimate$species) 
DHARMa::plotResiduals(simulationOutput.d.mod2, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod2, form = d.estimate$Temp)

## Visualize significant interactions
emmip(d.mod2, pop ~ MPA|temp) 
emmip(d.mod2, pop ~ temp|MPA) 
emmip(d.mod2, temp ~ MPA|pop) 

d.emm1.1 <- emmeans(d.mod2, pairwise ~ MPA|temp|pop, adjust='tukey', type = "response") 
d.emm1.2 <- emmeans(d.mod2, pairwise ~ temp|pop|MPA, adjust='tukey', type = "response") 
d.emm1.1$contrasts %>%
  summary(infer = TRUE)
d.emm1.2$contrasts %>%
  summary(infer = TRUE)

d.mod3.rg <- update(ref_grid(d.mod2), tran = "logit") # back-transform from logit scale
emmeans(d.mod3.rg, ~ MPA + temp + pop, type="response") # By default, 95% CIs; what's up with these df??
dmod3.RG.temp <- confint(emmeans(d.mod3.rg, ~ MPA + temp + pop, type="response"), adjust = "none", level = 0.95) 
dmod3.RG.temp$adj <- dmod3.RG.temp$response * 100
dmod3.RG.temp$adjLCL <- (dmod3.RG.temp$lower.CL)*100
dmod3.RG.temp$adjUCL <- (dmod3.RG.temp$upper.CL)*100
dmod3.RG.temp$adjSE <- (dmod3.RG.temp$SE)*100 

col.pal <- c("tomato3", "cadetblue", "goldenrod", "seagreen4", "slateblue4")
bw.pal <- c("gray80", "gray10", "gray40")
germ.plot <- ggplot(data=dmod3.RG.temp, aes(x = MPA, y = adj, fill = pop)) + 
  geom_errorbar(aes(ymin=adjLCL, ymax=adjUCL), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = pop, color = pop), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=pop), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Water potential (MPA)") + facet_wrap(~ temp) +
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw() 
germ.plot <- germ.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
germ.plot

germ.plot2 <- ggplot(data=dmod3.RG.temp, aes(x = MPA, y = adj, fill = pop)) + 
  geom_errorbar(aes(ymin=adjLCL, ymax=adjUCL), width = 0.3, size=0.5, position=position_dodge(width=0.5)) +
  geom_line(aes(group = pop, color = pop), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=pop), pch=21, color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Water potential (MPA)") + facet_wrap(~ temp) +
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw() 
germ.plot2 <- germ.plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
germ.plot2

# pull in observed data to overlay
d.estimate$adj <- d.estimate$estimate * 100
obs.data.pc <- ggplot(data=d.estimate, aes(x = temp, y = adj, fill=pop)) + 
  geom_point(aes(fill=pop), alpha=0.1, pch = 21, size = 3.5, position = position_dodge(width = 0.5)) +
  ylab("Maximum germination (%)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Temperature") + 
  scale_fill_manual(values=col.pal) + theme_bw()
obs.data.pc <- obs.data.pc + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
obs.data.pc <- obs.data.pc + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none")
obs.data.pc <- obs.data.pc + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(), strip.text.x = element_text(size = 8))
obs.data.pc

# Merge graphs
g1 <- ggplotGrob(germ.plot)
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
#ggsave("C:/Users/egtar/Desktop/wat-pt2.png", g, width=5, height=3.5, unit="in", dpi=800)


## e; t50 ---- 
e.estimate <- filter(germ.met.all, term == "e")
hist(log(e.estimate$estimate)) # big dip in middle...

ggplot(e.estimate, aes(x=pop, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter() +
  facet_wrap (~temp*MPA)

for (i in 1:nrow(e.estimate)){
  if (e.estimate$estimate[i] == 0.00000000) e.estimate$estimate[i] <-  0.0005
}

# Run model!
e.mod1 <- glm(log(estimate) ~ temp * MPA * pop,
             data=e.estimate, family=Gamma)
summary(e.mod1)
Anova(e.mod1)
simulationOutput.e.mod1 <- DHARMa::simulateResiduals(e.mod1, plot=T) # pretty good
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate$pop) 
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate$Temp) 

# Test the significance (and visualize) source
e.mod2 <- lmer(log(estimate) ~ temp * MPA +
                 (1|temp) + (1|temp:MPA:pop), 
               data=e.estimate, control=lmerControl(optimizer="bobyqa"))
lmtest::lrtest(e.mod1, e.mod2) # source is significant
sjPlot::plot_model(e.mod1, type="re", show.p=T) # see plot [[2]]
ranef(d.mod2) 
# https://cran.microsoft.com/snapshot/2017-04-21/web/packages/sjPlot/vignettes/sjplmer.html

# Try a random slopes model
e.mod.rs1 <- lmer(log(estimate) ~ temp * MPA +
                    (1|temp) + (1|temp:MPA:pop) + (temp|pop), 
                  data=e.estimate, control=lmerControl(optimizer="bobyqa"))
summary(e.mod.rs1)
sjPlot::plot_model(e.mod.rs1, type="re", show.p=T) # this does not work when you turn off covariance estimates with || in model
sjPlot::plot_model(e.mod.rs1,type="pred",
                  terms=c("temp","pop"),pred.type="re") # not really an interaction here..
# follow-up with Susan, can you do temp:MPA interaction with pop as random slopes?
# How should I visualize (and test 'significance') of the random effect of population in my model?

e.mod.rs3 <- lmer(log(estimate) ~ temp * MPA +
                    (1|temp) + (1|temp:MPA:pop) + (MPA|pop), 
                  data=e.estimate, control=lmerControl(optimizer="bobyqa"))
summary(e.mod.rs3)
sjPlot::plot_model(e.mod.rs3, type="re", show.p=T) # this does not work when you turn off covariance estimates with || in model
sjPlot::plot_model(e.mod.rs3,type="pred",
                   terms=c("MPA","pop"),pred.type="re") # definitely does not seem to be an interaction here

# Try model with populations as fixed effect
e.mod3 <- lmer(log(estimate) ~ temp * MPA * pop +
                 (1|temp) + (1|temp:MPA:pop),
               data=e.estimate, control=lmerControl(optimizer="bobyqa"))
# (1|temp) = chamber level variation; (1|temp:MPA) = variation among sets of cups across chambers; 
# (1|temp:MPA:pop) = variation among sets of cups within chambers; 
# residual = variation in germination within cup (i.e., response)
summary(e.mod3) 
Anova(e.mod3) # strange, not the results I would have expected; there's variation in pop germination, but it's not attributed to temp or MPA (or temp*MPA)??
plot(e.mod3) # decent
simulationOutput.e.mod3 <- DHARMa::simulateResiduals(d.mod3, plot=T)
DHARMa::plotResiduals(simulationOutput.e.mod3, form = e.estimate$pop) 
DHARMa::plotResiduals(simulationOutput.e.mod3, form = e.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.e.mod3, form = e.estimate$Temp) 

## Visualize significant interactions
emmip(e.mod3, ~ pop) 
emmip(e.mod3, ~ temp)
emmip(e.mod3, ~ MPA)

