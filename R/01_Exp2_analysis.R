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
                 (1|temp) + (1|temp:MPA) + (1|pop), 
               data=d.estimate, control=lmerControl(optimizer="bobyqa"))

summary(d.mod1) # Check with Susan - the repetition in random effect estimates looks fishy..
Anova(d.mod1)
plot(d.mod1) # outlier..
simulationOutput.d.mod1 <- DHARMa::simulateResiduals(d.mod1, plot=T) # ehh, not great; better than sqrt
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$pop) 
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$Temp) 
# there's quite a bit of underdispersion; overparameterizing model?

# Test the significance (and visualize) source
d.mod2 <- lmer(qlogis(estimate) ~ temp * MPA +
                 (1|temp) + (1|temp:MPA), 
               data=d.estimate, control=lmerControl(optimizer="bobyqa"))
lmtest::lrtest(d.mod1, d.mod2) # source is significant
sjPlot::plot_model(d.mod1, type="re", show.p=T)

ranef(d.mod2)
# https://cran.microsoft.com/snapshot/2017-04-21/web/packages/sjPlot/vignettes/sjplmer.html

lattice::dotplot(ranef(d.mod1,condVar=TRUE),
        lattice.options=list(layout=c(3,1)))

# Try a random slopes model
d.mod.rs1 <- lmer(qlogis(estimate) ~ temp * MPA +
                 (1|temp) + (1|temp:MPA) + (temp-1|pop), # should I suppress the intercept here?
               data=d.estimate, control=lmerControl(optimizer="bobyqa"))
summary(d.mod.rs1)
sjPlot::plot_model(d.mod.rs1, type="re", show.p=T)
# follow-up with Susan, can you do temp:MPA interaction with pop as random slopes?
# How should I visualize (and test 'significance') of the random effect of population in my model?

# Try model with populations as fixed effect
d.mod3 <- lmer(qlogis(estimate) ~ temp * MPA * pop +
                 (1|temp) + (1|temp:MPA) + (1|temp:MPA:pop), # I'm not sure if this last RE is necessary
               data=d.estimate, control=lmerControl(optimizer="bobyqa"))
summary(d.mod3) # Check with Susan - the repetition in random effect estimates looks fishy..
Anova(d.mod3)
plot(d.mod1) # outlier..
simulationOutput.d.mod1 <- DHARMa::simulateResiduals(d.mod1, plot=T) # ehh, not great; better than sqrt
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$pop) 
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$MPA)
DHARMa::plotResiduals(simulationOutput.d.mod1, form = d.estimate$Temp) 
# there's quite a bit of underdispersion; overparameterizing model?

##Pairwise comparisons for significant interactions\
emmip(d.mod3, pop ~ MPA) 
emmip(d.mod3, pop ~ temp)

## e; t50 ---- left off here 3/16/22
e.estimate <- filter(germ.met.all, term == "e")
hist(log(e.estimate$estimate))

ggplot(e.estimate, aes(x=species, y=estimate)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter() +
  facet_wrap (~temp*MPA)

for (i in 1:nrow(e.estimate)){
  if (e.estimate$estimate[i] == 0.00000000) e.estimate$estimate[i] <-  0.0005
}

# Look at paired t-test between 0 MPA at time 1 & 2; if not different, drop reps 4-6 of 0 MPA for models
# calculate t-test for MPA0 across data collection periods
time1 <- filter(e.estimate, MPA == 0 & curve != 4 & curve != 5 & curve != 6)
time2 <- filter(e.estimate, MPA == 0 & curve != 1 & curve != 2 & curve != 3)
all.time <- rbind(time1,time2)
ggplot(all.time, aes(x=curve, y=estimate, color=species)) +
  geom_point()

phau.time1 <- filter(time1, species == "PHAU")
phau.time2 <- filter(time2, species == "PHAU")
phau.t <- t.test(phau.time1$estimate, phau.time2$estimate, paired = TRUE); phau.t # accept null: no difference between time periods for PHAU

scac.time1 <- filter(time1, species == "SCAC")
scac.time2 <- filter(time2, species == "SCAC")
scac.t <- t.test(scac.time1$estimate, scac.time2$estimate, paired = TRUE); scac.t # reject null: there is a difference between time periods for SCAM

scam.time1 <- filter(time1, species == "SCAM")
scam.time2 <- filter(time2, species == "SCAM")
scam.t <- t.test(scam.time1$estimate, scam.time2$estimate, paired = TRUE); scam.t # reject null: there is a difference between time periods for SCAM

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

# drop reps 4-6 (will account for time as a random effect)
e.estimate2 <- filter(e.estimate, curve != 4 & curve != 5 & curve != 6)

# Run model!
e.mod1 <- lmer(log(estimate) ~ temp * MPA * species +
                 (1|temp) + (1|temp:MPA) + (1|temp:MPA:species) + (1|time), data=e.estimate2,
               control=lmerControl(optimizer="Nelder_Mead")) 
# (1|temp) = chamber level variation; (1|temp:MPA) = variation among sets of cups across chambers; 
# (1|temp:MPA:species) = variation among sets of cups within chambers; 
# (1|time) accounts for potential differences across the two time collections; 
# residual = variation in germination within cup (our response!)
summary(e.mod1) # Check with Susan - the repetition in random effect estimates looks fishy..
Anova(e.mod1)
plot(e.mod1)
simulationOutput.e.mod1 <- DHARMa::simulateResiduals(e.mod1, plot=T) # ehh, not great; better than sqrt
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate2$species) 
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate2$MPA)
DHARMa::plotResiduals(simulationOutput.e.mod1, form = e.estimate2$Temp)

# compare model with and without time
e.mod2 <- lmer(log(estimate) ~ temp * MPA * species +
                 (1|temp) + (1|temp:MPA) + (1|temp:MPA:species), data=e.estimate2,
               control=lmerControl(optimizer="Nelder_Mead")) 

lmtest::lrtest(e.mod1,e.mod2) # time is not siginficant

##Pairwise comparisons for significant interactions (temp*MPA*species)
emmip(e.mod1, species ~ MPA|temp) 
emmip(e.mod1, species ~ temp|MPA) 
# https://stackoverflow.com/questions/67338560/zeroinfl-model-warning-message-in-sqrtdiagobjectvcov-nans-produced

emmeans(e.mod1, pairwise ~ MPA|temp|species, adjust='tukey') 

e.mod1.rg <- update(ref_grid(e.mod1), tran = "log") # back-transform from log scale
emmeans(e.mod1.rg, ~ MPA + temp + species, type="response") # By default, 95% CIs
emm.e.mod1.rg <- confint(emmeans(e.mod1.rg, ~ MPA + temp + species, type="response"), adjust = "none", level = 0.68) # specify 68% CIs (roughly equivalent to +/- 1 SE)

# something is not working here - I think the model optimizer is preventing estimate of SE

col.pal <- c("tomato3", "cadetblue", "goldenrod")
bw.pal <- c("gray80", "gray10", "gray40")
t50.plot <- ggplot(data=emm.e.mod1.rg, aes(x = MPA, y = adj, shape = species, fill = species)) + 
  geom_errorbar(aes(ymin=adj-adjSE, ymax=adj+adjSE), width = 0.3, size=0.5, position=position_dodge(width=0.5), alpha=0.3) +
  geom_line(aes(group = species, color = species, alpha=0.05), size = 1, position=position_dodge(width = 0.5)) +
  geom_point(aes(fill=species, shape=species), color="black", size=3.5, position=position_dodge(width=0.5)) +
  ylab("Time to 50% germination (days)") + scale_y_continuous(breaks = seq(0,100,20), limits = c(0,100)) + 
  xlab("Water potential") + facet_wrap(~temp) + scale_shape_manual(values=c(21, 23, 24)) +
  scale_fill_manual(values=col.pal, name="fill") +
  scale_color_manual(values=col.pal) + 
  theme_bw()
t50.plot <- t50.plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
t50.plot
