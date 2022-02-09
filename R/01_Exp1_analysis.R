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
library(reshape2)
library(dplyr)
library(agricolae)
library(gtools)
library(lme4)
library(glmmTMB)
#devtools::install_github("onofriAndreaPG/drcSeedGerm")
library(drcSeedGerm)

# Load cleaned (non-cumulative) germination data ----
phau.exp1 <- read.csv("./data/exp1_PHAU_noncum.csv", as.is=T)
scac.exp1 <- read.csv("./data/exp1_SCAC_noncum.csv", as.is=T)
scam.exp1 <- read.csv("./data/exp1_SCAM_noncum.csv", as.is=T)

# Merge species data into one dataframe ----
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

# Time to event germination curve ----
## Flip dataframe ----
n.vec <- as.vector(rep(all.spp$N, each=23))
TTE.all.spp <- makeDrm(all.spp[,6:27], all.spp[,1:3], n.vec, 1:22)

plot(propCum ~ timeAf, data = TTE.all.spp, subset=c(is.finite(TTE.all.spp)==T), pch = 20,
     xlab= "Time (days)", ylab = "Proportion of germinated seeds", xlim=c(0,25),
     ylim = c(0, 1))
for(i in 2:15){
  polygon(c(TTE.all.spp$timeBef[i], TTE.all.spp$timeAf[i], TTE.all.spp$timeAf[i], TTE.all.spp$timeBef[i]), 
          c(TTE.all.spp$propCum[i-1], TTE.all.spp$propCum[i-1], TTE.all.spp$propCum[i], TTE.all.spp$propCum[i] ),
          density=100, col="grey", border=NA)
}
polygon(c(15, 17, 17, 15), 
        c(0.86, 0.86, 1, 1 ),
        density=100, col="grey")
lines(c(propCum, 0.86) ~ c(timeAf, 17), type="s", data = TTE.all.spp,
      col="blue", subset=is.finite(timeAf)==T)
lines(c(propCum, 1) ~ c(timeBef, 17), type="s", data = TTE.all.spp,
      col="red")


# https://www.statforbiology.com/seedgermination/censoring








# Calculate germination indices for response variables ----
## Germination proportion ----
all.spp$germsum <- rowSums(all.spp[,6:27]) # Calculate total germination per row
all.spp$germprop <- (all.spp$germsum / all.spp$N) # germination proportion

## Germination synchrony ----
# need to address interval length here! Data collection intervals varied across
# this experiment...
# calculated from a fitted curve of germination data vs from the raw data itself!

for (i in 1:nrow(all.spp)){
  x <- dod.noncum[i, 6:27] # pull out germination counts for each row
  germsum <- (sum(x))
  maxgerm <- dod.noncum$max.germ[i] # use max germination across all cold strat treatments per species * source
  germ50 <- maxgerm/2
  c.germ <- cumsum(as.numeric(x)) # create a vector of cumulative germination for each row
  interval <- as.matrix(seq(from=1, to=((21-7)*2)+1, by=2), dim(15,1)) # create data collection interval; days on which data was collected
  # Check if length of germination data and data collection intervals are of equal length
  if (length(x) != length(interval)) {
    stop("'germ counts' and 'interval' lengths differ")
  }
  
  if (germsum < germ50) {
    dod.noncum$d50[i] <- interval[nrow(interval),1] # Set d50 to upper limit of trial if there was no germination; double check that this comes through on germ dataframe
  }else {
    if (x[1] > germ50){
      dod.noncum$d50[i] <- 0 # Set d50 to 0 if 50 was hit on first day of trial
    } else{
      if (x[1] < germ50) {
        nearest <- c(match(max(c.germ[c.germ <= germ50]), c.germ), match(min(c.germ[c.germ >= germ50]), c.germ)) 
        if (nearest[2] == nearest[1]) {
          dod.noncum$d50[i] <- as.numeric(interval[nearest[1]]) 
        } else {
          if (nearest[2] > nearest[1]) {  
            dod.noncum$d50[i] <- interval[nearest[1]] + ((germ50 - c.germ[nearest[1]])*(interval[nearest[2]] - interval[nearest[1]]))/(c.germ[nearest[2]] - c.germ[nearest[1]])
          } else {
            dod.noncum$d50[i] <- "NA"
          }
        }
      }
    } 
  }
}


## t50 germination ----
for (i in 1:nrow(all.spp)){
  x <- dod.noncum[i, 6:27] # pull out germination counts for each row
  germsum <- (sum(x))
  maxgerm <- dod.noncum$max.germ[i] # use max germination across all cold strat treatments per species * source
  germ50 <- maxgerm/2
  c.germ <- cumsum(as.numeric(x)) # create a vector of cumulative germination for each row
  interval <- as.matrix(seq(from=1, to=((21-7)*2)+1, by=2), dim(15,1)) # create data collection interval; days on which data was collected
  # Check if length of germination data and data collection intervals are of equal length
  if (length(x) != length(interval)) {
    stop("'germ counts' and 'interval' lengths differ")
  }
  
  if (germsum < germ50) {
    dod.noncum$d50[i] <- interval[nrow(interval),1] # Set d50 to upper limit of trial if there was no germination; double check that this comes through on germ dataframe
  }else {
    if (x[1] > germ50){
      dod.noncum$d50[i] <- 0 # Set d50 to 0 if 50 was hit on first day of trial
    } else{
      if (x[1] < germ50) {
        nearest <- c(match(max(c.germ[c.germ <= germ50]), c.germ), match(min(c.germ[c.germ >= germ50]), c.germ)) 
        if (nearest[2] == nearest[1]) {
          dod.noncum$d50[i] <- as.numeric(interval[nearest[1]]) 
        } else {
          if (nearest[2] > nearest[1]) {  
            dod.noncum$d50[i] <- interval[nearest[1]] + ((germ50 - c.germ[nearest[1]])*(interval[nearest[2]] - interval[nearest[1]]))/(c.germ[nearest[2]] - c.germ[nearest[1]])
          } else {
            dod.noncum$d50[i] <- "NA"
          }
        }
      }
    } 
  }
}




# Build models ----
## Germ Proprtion Model ----
hist(all.spp$germprop)
fitdistrplus::descdist(all.spp$germprop) # beta distribution seems appropriate here

for (i in 1:nrow(all.spp)){
  if (all.spp$germprop[i] == 0.00000000) all.spp$germprop[i] <- 0.00000001
}
fitgrp.beta <- fitdistrplus::fitdist(all.spp$germprop, "beta")
fitgrp.beta <- fitdistrplus::fitdist(all.spp$germprop, "beta")
plot(fitgrp)

grp.mod <- glmmTMB()