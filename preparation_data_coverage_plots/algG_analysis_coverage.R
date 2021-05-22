# Arbeitsort spezifizieren
# Results_* kommen vom Skript des Auszählens und alle im gleichen Ordner
# Adapt number of isolates per pool if necessary
# check location of gene on fwd or reverse strand --> Adapt header in section input
# Check flag #### if some positions are SNPs in all isolates compared to PA14 --> correction needs to be adapted manually
# CF16, CF17 muss noch rein je nachdem wie es sequenziert wurde

setwd("/mnt/sfb900nfs/groups/tuemmler/sebastian/Berichte/pathoadaptive_mutationen/algG")
library(WriteXLS)
library(openxlsx)
library(tibble)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)

# Import Vorlage

#input <- read.xlsx("algG_Vorlage.xlsx")

# Set number of isolates per pooL
isolates_CF_1 <- 20
isolates_CF_2 <- 20
isolates_CF_3 <- 20
isolates_CF_4 <- 20
isolates_CF_5 <- 20
isolates_CF_6 <- 20
isolates_CF_7 <- 20
isolates_CF_8 <- 20
isolates_CF_9 <- 20
isolates_CF_10 <- 20
isolates_CF_11 <- 20
isolates_CF_12 <- 20
isolates_CF_13 <- 20
isolates_CF_14 <- 23
isolates_CF_15 <- 22
isolates_CF_16 <- 20
isolates_CF_17 <- 20
isolates_acute_1 <- 21
isolates_acute_2 <- 23
isolates_acute_3 <- 23
isolates_env_1 <- 20
isolates_env_2  <- 20
isolates_env_3 <- 20
isolates_env_4 <- 19
isolates_env_5 <- 18
isolates_COPD <- 25


# Import bases for coordinates

PA14_bases <- read.table("algG_bases.txt",header = TRUE)
names(PA14_bases) <- c("PA14-Koordinate","GenBase_PA14")
PA14_bases_5UTR <- read.table("algG_5UTR_bases.txt",header = TRUE)
names(PA14_bases_5UTR) <- c("PA14-Koordinate","GenBase_PA14")
PA14_bases_3UTR <- read.table("algG_3UTR_bases.txt",header = TRUE)
names(PA14_bases_3UTR) <- c("PA14-Koordinate","GenBase_PA14")


# Gene is on reverse strand so the As are Ts and the Cs Gs and vice versa in the genome. Column names are changed to have the correct data

### CF_1

# Input 
data_algG_CF_1 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T1/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_1) <- c("PA14-Koordinate","Ts_CF_1","As_CF_1","Gs_CF_1","Cs_CF_1")
data_algG_CF_1 <- merge(PA14_bases, data_algG_CF_1, "PA14-Koordinate")
data_algG_CF_1_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T1/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_1_5UTR) <- c("PA14-Koordinate","As_CF_1_5UTR","Ts_CF_1_5UTR","Cs_CF_1_5UTR","Gs_CF_1_5UTR")
data_algG_CF_1_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_1_5UTR, "PA14-Koordinate")
data_algG_CF_1_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T1/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_1_3UTR) <- c("PA14-Koordinate","As_CF_1_3UTR","Ts_CF_1_3UTR","Cs_CF_1_3UTR","Gs_CF_1_3UTR")
data_algG_CF_1_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_1_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_1$coverage_CF_1 <- data_algG_CF_1$Ts_CF_1 + data_algG_CF_1$As_CF_1 + data_algG_CF_1$Gs_CF_1 + data_algG_CF_1$Cs_CF_1
data_algG_CF_1_5UTR$coverage_CF_1_5UTR <- data_algG_CF_1_5UTR$Ts_CF_1_5UTR + data_algG_CF_1_5UTR$As_CF_1_5UTR + data_algG_CF_1_5UTR$Gs_CF_1_5UTR + data_algG_CF_1_5UTR$Cs_CF_1_5UTR
data_algG_CF_1_3UTR$coverage_CF_1_3UTR <- data_algG_CF_1_3UTR$Ts_CF_1_3UTR + data_algG_CF_1_3UTR$As_CF_1_3UTR + data_algG_CF_1_3UTR$Gs_CF_1_3UTR + data_algG_CF_1_3UTR$Cs_CF_1_3UTR


### CF_2

# Input 
data_algG_CF_2 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T2/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_2) <- c("PA14-Koordinate","Ts_CF_2","As_CF_2","Gs_CF_2","Cs_CF_2")
data_algG_CF_2 <- merge(PA14_bases, data_algG_CF_2, "PA14-Koordinate")
data_algG_CF_2_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T2/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_2_5UTR) <- c("PA14-Koordinate","As_CF_2_5UTR","Ts_CF_2_5UTR","Cs_CF_2_5UTR","Gs_CF_2_5UTR")
data_algG_CF_2_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_2_5UTR, "PA14-Koordinate")
data_algG_CF_2_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T2/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_2_3UTR) <- c("PA14-Koordinate","As_CF_2_3UTR","Ts_CF_2_3UTR","Cs_CF_2_3UTR","Gs_CF_2_3UTR")
data_algG_CF_2_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_2_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_2$coverage_CF_2 <- data_algG_CF_2$Ts_CF_2 + data_algG_CF_2$As_CF_2 + data_algG_CF_2$Gs_CF_2 + data_algG_CF_2$Cs_CF_2
data_algG_CF_2_5UTR$coverage_CF_2_5UTR <- data_algG_CF_2_5UTR$Ts_CF_2_5UTR + data_algG_CF_2_5UTR$As_CF_2_5UTR + data_algG_CF_2_5UTR$Gs_CF_2_5UTR + data_algG_CF_2_5UTR$Cs_CF_2_5UTR
data_algG_CF_2_3UTR$coverage_CF_2_3UTR <- data_algG_CF_2_3UTR$Ts_CF_2_3UTR + data_algG_CF_2_3UTR$As_CF_2_3UTR + data_algG_CF_2_3UTR$Gs_CF_2_3UTR + data_algG_CF_2_3UTR$Cs_CF_2_3UTR


### CF_3

# Input 
data_algG_CF_3 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T3/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_3) <- c("PA14-Koordinate","Ts_CF_3","As_CF_3","Gs_CF_3","Cs_CF_3")
data_algG_CF_3 <- merge(PA14_bases, data_algG_CF_3, "PA14-Koordinate")
data_algG_CF_3_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T3/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_3_5UTR) <- c("PA14-Koordinate","As_CF_3_5UTR","Ts_CF_3_5UTR","Cs_CF_3_5UTR","Gs_CF_3_5UTR")
data_algG_CF_3_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_3_5UTR, "PA14-Koordinate")
data_algG_CF_3_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T3/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_3_3UTR) <- c("PA14-Koordinate","As_CF_3_3UTR","Ts_CF_3_3UTR","Cs_CF_3_3UTR","Gs_CF_3_3UTR")
data_algG_CF_3_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_3_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_3$coverage_CF_3 <- data_algG_CF_3$Ts_CF_3 + data_algG_CF_3$As_CF_3 + data_algG_CF_3$Gs_CF_3 + data_algG_CF_3$Cs_CF_3
data_algG_CF_3_5UTR$coverage_CF_3_5UTR <- data_algG_CF_3_5UTR$Ts_CF_3_5UTR + data_algG_CF_3_5UTR$As_CF_3_5UTR + data_algG_CF_3_5UTR$Gs_CF_3_5UTR + data_algG_CF_3_5UTR$Cs_CF_3_5UTR
data_algG_CF_3_3UTR$coverage_CF_3_3UTR <- data_algG_CF_3_3UTR$Ts_CF_3_3UTR + data_algG_CF_3_3UTR$As_CF_3_3UTR + data_algG_CF_3_3UTR$Gs_CF_3_3UTR + data_algG_CF_3_3UTR$Cs_CF_3_3UTR

### CF_4

# Input 
data_algG_CF_4 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T4/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_4) <- c("PA14-Koordinate","Ts_CF_4","As_CF_4","Gs_CF_4","Cs_CF_4")
data_algG_CF_4 <- merge(PA14_bases, data_algG_CF_4, "PA14-Koordinate")
data_algG_CF_4_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T4/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_4_5UTR) <- c("PA14-Koordinate","As_CF_4_5UTR","Ts_CF_4_5UTR","Cs_CF_4_5UTR","Gs_CF_4_5UTR")
data_algG_CF_4_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_4_5UTR, "PA14-Koordinate")
data_algG_CF_4_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T4/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_4_3UTR) <- c("PA14-Koordinate","As_CF_4_3UTR","Ts_CF_4_3UTR","Cs_CF_4_3UTR","Gs_CF_4_3UTR")
data_algG_CF_4_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_4_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_4$coverage_CF_4 <- data_algG_CF_4$Ts_CF_4 + data_algG_CF_4$As_CF_4 + data_algG_CF_4$Gs_CF_4 + data_algG_CF_4$Cs_CF_4
data_algG_CF_4_5UTR$coverage_CF_4_5UTR <- data_algG_CF_4_5UTR$Ts_CF_4_5UTR + data_algG_CF_4_5UTR$As_CF_4_5UTR + data_algG_CF_4_5UTR$Gs_CF_4_5UTR + data_algG_CF_4_5UTR$Cs_CF_4_5UTR
data_algG_CF_4_3UTR$coverage_CF_4_3UTR <- data_algG_CF_4_3UTR$Ts_CF_4_3UTR + data_algG_CF_4_3UTR$As_CF_4_3UTR + data_algG_CF_4_3UTR$Gs_CF_4_3UTR + data_algG_CF_4_3UTR$Cs_CF_4_3UTR

### CF_5

# Input 
data_algG_CF_5 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T5/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_5) <- c("PA14-Koordinate","Ts_CF_5","As_CF_5","Gs_CF_5","Cs_CF_5")
data_algG_CF_5 <- merge(PA14_bases, data_algG_CF_5, "PA14-Koordinate")
data_algG_CF_5_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T5/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_5_5UTR) <- c("PA14-Koordinate","As_CF_5_5UTR","Ts_CF_5_5UTR","Cs_CF_5_5UTR","Gs_CF_5_5UTR")
data_algG_CF_5_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_5_5UTR, "PA14-Koordinate")
data_algG_CF_5_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T5/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_5_3UTR) <- c("PA14-Koordinate","As_CF_5_3UTR","Ts_CF_5_3UTR","Cs_CF_5_3UTR","Gs_CF_5_3UTR")
data_algG_CF_5_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_5_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_5$coverage_CF_5 <- data_algG_CF_5$Ts_CF_5 + data_algG_CF_5$As_CF_5 + data_algG_CF_5$Gs_CF_5 + data_algG_CF_5$Cs_CF_5
data_algG_CF_5_5UTR$coverage_CF_5_5UTR <- data_algG_CF_5_5UTR$Ts_CF_5_5UTR + data_algG_CF_5_5UTR$As_CF_5_5UTR + data_algG_CF_5_5UTR$Gs_CF_5_5UTR + data_algG_CF_5_5UTR$Cs_CF_5_5UTR
data_algG_CF_5_3UTR$coverage_CF_5_3UTR <- data_algG_CF_5_3UTR$Ts_CF_5_3UTR + data_algG_CF_5_3UTR$As_CF_5_3UTR + data_algG_CF_5_3UTR$Gs_CF_5_3UTR + data_algG_CF_5_3UTR$Cs_CF_5_3UTR

### CF_6

# Input 
data_algG_CF_6 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T6/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_6) <- c("PA14-Koordinate","Ts_CF_6","As_CF_6","Gs_CF_6","Cs_CF_6")
data_algG_CF_6 <- merge(PA14_bases, data_algG_CF_6, "PA14-Koordinate")
data_algG_CF_6_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T6/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_6_5UTR) <- c("PA14-Koordinate","As_CF_6_5UTR","Ts_CF_6_5UTR","Cs_CF_6_5UTR","Gs_CF_6_5UTR")
data_algG_CF_6_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_6_5UTR, "PA14-Koordinate")
data_algG_CF_6_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T6/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_6_3UTR) <- c("PA14-Koordinate","As_CF_6_3UTR","Ts_CF_6_3UTR","Cs_CF_6_3UTR","Gs_CF_6_3UTR")
data_algG_CF_6_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_6_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_6$coverage_CF_6 <- data_algG_CF_6$Ts_CF_6 + data_algG_CF_6$As_CF_6 + data_algG_CF_6$Gs_CF_6 + data_algG_CF_6$Cs_CF_6
data_algG_CF_6_5UTR$coverage_CF_6_5UTR <- data_algG_CF_6_5UTR$Ts_CF_6_5UTR + data_algG_CF_6_5UTR$As_CF_6_5UTR + data_algG_CF_6_5UTR$Gs_CF_6_5UTR + data_algG_CF_6_5UTR$Cs_CF_6_5UTR
data_algG_CF_6_3UTR$coverage_CF_6_3UTR <- data_algG_CF_6_3UTR$Ts_CF_6_3UTR + data_algG_CF_6_3UTR$As_CF_6_3UTR + data_algG_CF_6_3UTR$Gs_CF_6_3UTR + data_algG_CF_6_3UTR$Cs_CF_6_3UTR

### CF_7

# Input 
data_algG_CF_7 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T7/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_7) <- c("PA14-Koordinate","Ts_CF_7","As_CF_7","Gs_CF_7","Cs_CF_7")
data_algG_CF_7 <- merge(PA14_bases, data_algG_CF_7, "PA14-Koordinate")
data_algG_CF_7_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T7/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_7_5UTR) <- c("PA14-Koordinate","As_CF_7_5UTR","Ts_CF_7_5UTR","Cs_CF_7_5UTR","Gs_CF_7_5UTR")
data_algG_CF_7_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_7_5UTR, "PA14-Koordinate")
data_algG_CF_7_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T7/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_7_3UTR) <- c("PA14-Koordinate","As_CF_7_3UTR","Ts_CF_7_3UTR","Cs_CF_7_3UTR","Gs_CF_7_3UTR")
data_algG_CF_7_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_7_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_7$coverage_CF_7 <- data_algG_CF_7$Ts_CF_7 + data_algG_CF_7$As_CF_7 + data_algG_CF_7$Gs_CF_7 + data_algG_CF_7$Cs_CF_7
data_algG_CF_7_5UTR$coverage_CF_7_5UTR <- data_algG_CF_7_5UTR$Ts_CF_7_5UTR + data_algG_CF_7_5UTR$As_CF_7_5UTR + data_algG_CF_7_5UTR$Gs_CF_7_5UTR + data_algG_CF_7_5UTR$Cs_CF_7_5UTR
data_algG_CF_7_3UTR$coverage_CF_7_3UTR <- data_algG_CF_7_3UTR$Ts_CF_7_3UTR + data_algG_CF_7_3UTR$As_CF_7_3UTR + data_algG_CF_7_3UTR$Gs_CF_7_3UTR + data_algG_CF_7_3UTR$Cs_CF_7_3UTR

### CF_8

# Input 
data_algG_CF_8 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T8/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_8) <- c("PA14-Koordinate","Ts_CF_8","As_CF_8","Gs_CF_8","Cs_CF_8")
data_algG_CF_8 <- merge(PA14_bases, data_algG_CF_8, "PA14-Koordinate")
data_algG_CF_8_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T8/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_8_5UTR) <- c("PA14-Koordinate","As_CF_8_5UTR","Ts_CF_8_5UTR","Cs_CF_8_5UTR","Gs_CF_8_5UTR")
data_algG_CF_8_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_8_5UTR, "PA14-Koordinate")
data_algG_CF_8_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T8/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_8_3UTR) <- c("PA14-Koordinate","As_CF_8_3UTR","Ts_CF_8_3UTR","Cs_CF_8_3UTR","Gs_CF_8_3UTR")
data_algG_CF_8_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_8_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_8$coverage_CF_8 <- data_algG_CF_8$Ts_CF_8 + data_algG_CF_8$As_CF_8 + data_algG_CF_8$Gs_CF_8 + data_algG_CF_8$Cs_CF_8
data_algG_CF_8_5UTR$coverage_CF_8_5UTR <- data_algG_CF_8_5UTR$Ts_CF_8_5UTR + data_algG_CF_8_5UTR$As_CF_8_5UTR + data_algG_CF_8_5UTR$Gs_CF_8_5UTR + data_algG_CF_8_5UTR$Cs_CF_8_5UTR
data_algG_CF_8_3UTR$coverage_CF_8_3UTR <- data_algG_CF_8_3UTR$Ts_CF_8_3UTR + data_algG_CF_8_3UTR$As_CF_8_3UTR + data_algG_CF_8_3UTR$Gs_CF_8_3UTR + data_algG_CF_8_3UTR$Cs_CF_8_3UTR

### CF_9

# Input 
data_algG_CF_9 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T9/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_9) <- c("PA14-Koordinate","Ts_CF_9","As_CF_9","Gs_CF_9","Cs_CF_9")
data_algG_CF_9 <- merge(PA14_bases, data_algG_CF_9, "PA14-Koordinate")
data_algG_CF_9_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T9/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_9_5UTR) <- c("PA14-Koordinate","As_CF_9_5UTR","Ts_CF_9_5UTR","Cs_CF_9_5UTR","Gs_CF_9_5UTR")
data_algG_CF_9_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_9_5UTR, "PA14-Koordinate")
data_algG_CF_9_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T9/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_9_3UTR) <- c("PA14-Koordinate","As_CF_9_3UTR","Ts_CF_9_3UTR","Cs_CF_9_3UTR","Gs_CF_9_3UTR")
data_algG_CF_9_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_9_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_9$coverage_CF_9 <- data_algG_CF_9$Ts_CF_9 + data_algG_CF_9$As_CF_9 + data_algG_CF_9$Gs_CF_9 + data_algG_CF_9$Cs_CF_9
data_algG_CF_9_5UTR$coverage_CF_9_5UTR <- data_algG_CF_9_5UTR$Ts_CF_9_5UTR + data_algG_CF_9_5UTR$As_CF_9_5UTR + data_algG_CF_9_5UTR$Gs_CF_9_5UTR + data_algG_CF_9_5UTR$Cs_CF_9_5UTR
data_algG_CF_9_3UTR$coverage_CF_9_3UTR <- data_algG_CF_9_3UTR$Ts_CF_9_3UTR + data_algG_CF_9_3UTR$As_CF_9_3UTR + data_algG_CF_9_3UTR$Gs_CF_9_3UTR + data_algG_CF_9_3UTR$Cs_CF_9_3UTR

### CF_10

# Input 
data_algG_CF_10 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T10/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_10) <- c("PA14-Koordinate","Ts_CF_10","As_CF_10","Gs_CF_10","Cs_CF_10")
data_algG_CF_10 <- merge(PA14_bases, data_algG_CF_10, "PA14-Koordinate")
data_algG_CF_10_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T10/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_10_5UTR) <- c("PA14-Koordinate","As_CF_10_5UTR","Ts_CF_10_5UTR","Cs_CF_10_5UTR","Gs_CF_10_5UTR")
data_algG_CF_10_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_10_5UTR, "PA14-Koordinate")
data_algG_CF_10_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T10/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_10_3UTR) <- c("PA14-Koordinate","As_CF_10_3UTR","Ts_CF_10_3UTR","Cs_CF_10_3UTR","Gs_CF_10_3UTR")
data_algG_CF_10_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_10_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_10$coverage_CF_10 <- data_algG_CF_10$Ts_CF_10 + data_algG_CF_10$As_CF_10 + data_algG_CF_10$Gs_CF_10 + data_algG_CF_10$Cs_CF_10
data_algG_CF_10_5UTR$coverage_CF_10_5UTR <- data_algG_CF_10_5UTR$Ts_CF_10_5UTR + data_algG_CF_10_5UTR$As_CF_10_5UTR + data_algG_CF_10_5UTR$Gs_CF_10_5UTR + data_algG_CF_10_5UTR$Cs_CF_10_5UTR
data_algG_CF_10_3UTR$coverage_CF_10_3UTR <- data_algG_CF_10_3UTR$Ts_CF_10_3UTR + data_algG_CF_10_3UTR$As_CF_10_3UTR + data_algG_CF_10_3UTR$Gs_CF_10_3UTR + data_algG_CF_10_3UTR$Cs_CF_10_3UTR

### CF_11

# Input 
data_algG_CF_11 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T11/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_11) <- c("PA14-Koordinate","Ts_CF_11","As_CF_11","Gs_CF_11","Cs_CF_11")
data_algG_CF_11 <- merge(PA14_bases, data_algG_CF_11, "PA14-Koordinate")
data_algG_CF_11_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T11/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_11_5UTR) <- c("PA14-Koordinate","As_CF_11_5UTR","Ts_CF_11_5UTR","Cs_CF_11_5UTR","Gs_CF_11_5UTR")
data_algG_CF_11_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_11_5UTR, "PA14-Koordinate")
data_algG_CF_11_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T11/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_11_3UTR) <- c("PA14-Koordinate","As_CF_11_3UTR","Ts_CF_11_3UTR","Cs_CF_11_3UTR","Gs_CF_11_3UTR")
data_algG_CF_11_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_11_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_11$coverage_CF_11 <- data_algG_CF_11$Ts_CF_11 + data_algG_CF_11$As_CF_11 + data_algG_CF_11$Gs_CF_11 + data_algG_CF_11$Cs_CF_11
data_algG_CF_11_5UTR$coverage_CF_11_5UTR <- data_algG_CF_11_5UTR$Ts_CF_11_5UTR + data_algG_CF_11_5UTR$As_CF_11_5UTR + data_algG_CF_11_5UTR$Gs_CF_11_5UTR + data_algG_CF_11_5UTR$Cs_CF_11_5UTR
data_algG_CF_11_3UTR$coverage_CF_11_3UTR <- data_algG_CF_11_3UTR$Ts_CF_11_3UTR + data_algG_CF_11_3UTR$As_CF_11_3UTR + data_algG_CF_11_3UTR$Gs_CF_11_3UTR + data_algG_CF_11_3UTR$Cs_CF_11_3UTR

### CF_12

# Input 
data_algG_CF_12 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T12/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_12) <- c("PA14-Koordinate","Ts_CF_12","As_CF_12","Gs_CF_12","Cs_CF_12")
data_algG_CF_12 <- merge(PA14_bases, data_algG_CF_12, "PA14-Koordinate")
data_algG_CF_12_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T12/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_12_5UTR) <- c("PA14-Koordinate","As_CF_12_5UTR","Ts_CF_12_5UTR","Cs_CF_12_5UTR","Gs_CF_12_5UTR")
data_algG_CF_12_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_12_5UTR, "PA14-Koordinate")
data_algG_CF_12_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T12/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_12_3UTR) <- c("PA14-Koordinate","As_CF_12_3UTR","Ts_CF_12_3UTR","Cs_CF_12_3UTR","Gs_CF_12_3UTR")
data_algG_CF_12_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_12_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_12$coverage_CF_12 <- data_algG_CF_12$Ts_CF_12 + data_algG_CF_12$As_CF_12 + data_algG_CF_12$Gs_CF_12 + data_algG_CF_12$Cs_CF_12
data_algG_CF_12_5UTR$coverage_CF_12_5UTR <- data_algG_CF_12_5UTR$Ts_CF_12_5UTR + data_algG_CF_12_5UTR$As_CF_12_5UTR + data_algG_CF_12_5UTR$Gs_CF_12_5UTR + data_algG_CF_12_5UTR$Cs_CF_12_5UTR
data_algG_CF_12_3UTR$coverage_CF_12_3UTR <- data_algG_CF_12_3UTR$Ts_CF_12_3UTR + data_algG_CF_12_3UTR$As_CF_12_3UTR + data_algG_CF_12_3UTR$Gs_CF_12_3UTR + data_algG_CF_12_3UTR$Cs_CF_12_3UTR

### CF_13

# Input 
data_algG_CF_13 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T13/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_13) <- c("PA14-Koordinate","Ts_CF_13","As_CF_13","Gs_CF_13","Cs_CF_13")
data_algG_CF_13 <- merge(PA14_bases, data_algG_CF_13, "PA14-Koordinate")
data_algG_CF_13_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T13/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_13_5UTR) <- c("PA14-Koordinate","As_CF_13_5UTR","Ts_CF_13_5UTR","Cs_CF_13_5UTR","Gs_CF_13_5UTR")
data_algG_CF_13_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_13_5UTR, "PA14-Koordinate")
data_algG_CF_13_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T13/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_13_3UTR) <- c("PA14-Koordinate","As_CF_13_3UTR","Ts_CF_13_3UTR","Cs_CF_13_3UTR","Gs_CF_13_3UTR")
data_algG_CF_13_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_13_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_13$coverage_CF_13 <- data_algG_CF_13$Ts_CF_13 + data_algG_CF_13$As_CF_13 + data_algG_CF_13$Gs_CF_13 + data_algG_CF_13$Cs_CF_13
data_algG_CF_13_5UTR$coverage_CF_13_5UTR <- data_algG_CF_13_5UTR$Ts_CF_13_5UTR + data_algG_CF_13_5UTR$As_CF_13_5UTR + data_algG_CF_13_5UTR$Gs_CF_13_5UTR + data_algG_CF_13_5UTR$Cs_CF_13_5UTR
data_algG_CF_13_3UTR$coverage_CF_13_3UTR <- data_algG_CF_13_3UTR$Ts_CF_13_3UTR + data_algG_CF_13_3UTR$As_CF_13_3UTR + data_algG_CF_13_3UTR$Gs_CF_13_3UTR + data_algG_CF_13_3UTR$Cs_CF_13_3UTR

### CF_14

# Input 
data_algG_CF_14 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T14/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_14) <- c("PA14-Koordinate","Ts_CF_14","As_CF_14","Gs_CF_14","Cs_CF_14")
data_algG_CF_14 <- merge(PA14_bases, data_algG_CF_14, "PA14-Koordinate")
data_algG_CF_14_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T14/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_14_5UTR) <- c("PA14-Koordinate","As_CF_14_5UTR","Ts_CF_14_5UTR","Cs_CF_14_5UTR","Gs_CF_14_5UTR")
data_algG_CF_14_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_14_5UTR, "PA14-Koordinate")
data_algG_CF_14_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T14/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_14_3UTR) <- c("PA14-Koordinate","As_CF_14_3UTR","Ts_CF_14_3UTR","Cs_CF_14_3UTR","Gs_CF_14_3UTR")
data_algG_CF_14_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_14_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_14$coverage_CF_14 <- data_algG_CF_14$Ts_CF_14 + data_algG_CF_14$As_CF_14 + data_algG_CF_14$Gs_CF_14 + data_algG_CF_14$Cs_CF_14
data_algG_CF_14_5UTR$coverage_CF_14_5UTR <- data_algG_CF_14_5UTR$Ts_CF_14_5UTR + data_algG_CF_14_5UTR$As_CF_14_5UTR + data_algG_CF_14_5UTR$Gs_CF_14_5UTR + data_algG_CF_14_5UTR$Cs_CF_14_5UTR
data_algG_CF_14_3UTR$coverage_CF_14_3UTR <- data_algG_CF_14_3UTR$Ts_CF_14_3UTR + data_algG_CF_14_3UTR$As_CF_14_3UTR + data_algG_CF_14_3UTR$Gs_CF_14_3UTR + data_algG_CF_14_3UTR$Cs_CF_14_3UTR

### CF_15

# Input 
data_algG_CF_15 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T15/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_15) <- c("PA14-Koordinate","Ts_CF_15","As_CF_15","Gs_CF_15","Cs_CF_15")
data_algG_CF_15 <- merge(PA14_bases, data_algG_CF_15, "PA14-Koordinate")
data_algG_CF_15_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T15/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_15_5UTR) <- c("PA14-Koordinate","As_CF_15_5UTR","Ts_CF_15_5UTR","Cs_CF_15_5UTR","Gs_CF_15_5UTR")
data_algG_CF_15_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_15_5UTR, "PA14-Koordinate")
data_algG_CF_15_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T15/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_15_3UTR) <- c("PA14-Koordinate","As_CF_15_3UTR","Ts_CF_15_3UTR","Cs_CF_15_3UTR","Gs_CF_15_3UTR")
data_algG_CF_15_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_15_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_15$coverage_CF_15 <- data_algG_CF_15$Ts_CF_15 + data_algG_CF_15$As_CF_15 + data_algG_CF_15$Gs_CF_15 + data_algG_CF_15$Cs_CF_15
data_algG_CF_15_5UTR$coverage_CF_15_5UTR <- data_algG_CF_15_5UTR$Ts_CF_15_5UTR + data_algG_CF_15_5UTR$As_CF_15_5UTR + data_algG_CF_15_5UTR$Gs_CF_15_5UTR + data_algG_CF_15_5UTR$Cs_CF_15_5UTR
data_algG_CF_15_3UTR$coverage_CF_15_3UTR <- data_algG_CF_15_3UTR$Ts_CF_15_3UTR + data_algG_CF_15_3UTR$As_CF_15_3UTR + data_algG_CF_15_3UTR$Gs_CF_15_3UTR + data_algG_CF_15_3UTR$Cs_CF_15_3UTR

### CF_16

# Input 
data_algG_CF_16 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T16/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_16) <- c("PA14-Koordinate","Ts_CF_16","As_CF_16","Gs_CF_16","Cs_CF_16")
data_algG_CF_16 <- merge(PA14_bases, data_algG_CF_16, "PA14-Koordinate")
data_algG_CF_16_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T16/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_16_5UTR) <- c("PA14-Koordinate","As_CF_16_5UTR","Ts_CF_16_5UTR","Cs_CF_16_5UTR","Gs_CF_16_5UTR")
data_algG_CF_16_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_16_5UTR, "PA14-Koordinate")
data_algG_CF_16_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T16/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_16_3UTR) <- c("PA14-Koordinate","As_CF_16_3UTR","Ts_CF_16_3UTR","Cs_CF_16_3UTR","Gs_CF_16_3UTR")
data_algG_CF_16_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_16_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_16$coverage_CF_16 <- data_algG_CF_16$Ts_CF_16 + data_algG_CF_16$As_CF_16 + data_algG_CF_16$Gs_CF_16 + data_algG_CF_16$Cs_CF_16
data_algG_CF_16_5UTR$coverage_CF_16_5UTR <- data_algG_CF_16_5UTR$Ts_CF_16_5UTR + data_algG_CF_16_5UTR$As_CF_16_5UTR + data_algG_CF_16_5UTR$Gs_CF_16_5UTR + data_algG_CF_16_5UTR$Cs_CF_16_5UTR
data_algG_CF_16_3UTR$coverage_CF_16_3UTR <- data_algG_CF_16_3UTR$Ts_CF_16_3UTR + data_algG_CF_16_3UTR$As_CF_16_3UTR + data_algG_CF_16_3UTR$Gs_CF_16_3UTR + data_algG_CF_16_3UTR$Cs_CF_16_3UTR

### CF_17

# Input 
data_algG_CF_17 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T17/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_CF_17) <- c("PA14-Koordinate","Ts_CF_17","As_CF_17","Gs_CF_17","Cs_CF_17")
data_algG_CF_17 <- merge(PA14_bases, data_algG_CF_17, "PA14-Koordinate")
data_algG_CF_17_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T17/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_CF_17_5UTR) <- c("PA14-Koordinate","As_CF_17_5UTR","Ts_CF_17_5UTR","Cs_CF_17_5UTR","Gs_CF_17_5UTR")
data_algG_CF_17_5UTR <- merge(PA14_bases_5UTR, data_algG_CF_17_5UTR, "PA14-Koordinate")
data_algG_CF_17_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T17/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_CF_17_3UTR) <- c("PA14-Koordinate","As_CF_17_3UTR","Ts_CF_17_3UTR","Cs_CF_17_3UTR","Gs_CF_17_3UTR")
data_algG_CF_17_3UTR <- merge(PA14_bases_3UTR, data_algG_CF_17_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_CF_17$coverage_CF_17 <- data_algG_CF_17$Ts_CF_17 + data_algG_CF_17$As_CF_17 + data_algG_CF_17$Gs_CF_17 + data_algG_CF_17$Cs_CF_17
data_algG_CF_17_5UTR$coverage_CF_17_5UTR <- data_algG_CF_17_5UTR$Ts_CF_17_5UTR + data_algG_CF_17_5UTR$As_CF_17_5UTR + data_algG_CF_17_5UTR$Gs_CF_17_5UTR + data_algG_CF_17_5UTR$Cs_CF_17_5UTR
data_algG_CF_17_3UTR$coverage_CF_17_3UTR <- data_algG_CF_17_3UTR$Ts_CF_17_3UTR + data_algG_CF_17_3UTR$As_CF_17_3UTR + data_algG_CF_17_3UTR$Gs_CF_17_3UTR + data_algG_CF_17_3UTR$Cs_CF_17_3UTR

### acute_1

# Input 
data_algG_acute_1 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T18/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_acute_1) <- c("PA14-Koordinate","Ts_acute_1","As_acute_1","Gs_acute_1","Cs_acute_1")
data_algG_acute_1 <- merge(PA14_bases, data_algG_acute_1, "PA14-Koordinate")
data_algG_acute_1_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T18/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_acute_1_5UTR) <- c("PA14-Koordinate","As_acute_1_5UTR","Ts_acute_1_5UTR","Cs_acute_1_5UTR","Gs_acute_1_5UTR")
data_algG_acute_1_5UTR <- merge(PA14_bases_5UTR, data_algG_acute_1_5UTR, "PA14-Koordinate")
data_algG_acute_1_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T18/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_acute_1_3UTR) <- c("PA14-Koordinate","As_acute_1_3UTR","Ts_acute_1_3UTR","Cs_acute_1_3UTR","Gs_acute_1_3UTR")
data_algG_acute_1_3UTR <- merge(PA14_bases_3UTR, data_algG_acute_1_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_acute_1$coverage_acute_1 <- data_algG_acute_1$Ts_acute_1 + data_algG_acute_1$As_acute_1 + data_algG_acute_1$Gs_acute_1 + data_algG_acute_1$Cs_acute_1
data_algG_acute_1_5UTR$coverage_acute_1_5UTR <- data_algG_acute_1_5UTR$Ts_acute_1_5UTR + data_algG_acute_1_5UTR$As_acute_1_5UTR + data_algG_acute_1_5UTR$Gs_acute_1_5UTR + data_algG_acute_1_5UTR$Cs_acute_1_5UTR
data_algG_acute_1_3UTR$coverage_acute_1_3UTR <- data_algG_acute_1_3UTR$Ts_acute_1_3UTR + data_algG_acute_1_3UTR$As_acute_1_3UTR + data_algG_acute_1_3UTR$Gs_acute_1_3UTR + data_algG_acute_1_3UTR$Cs_acute_1_3UTR

### acute_2

# Input 
data_algG_acute_2 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T19/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_acute_2) <- c("PA14-Koordinate","Ts_acute_2","As_acute_2","Gs_acute_2","Cs_acute_2")
data_algG_acute_2 <- merge(PA14_bases, data_algG_acute_2, "PA14-Koordinate")
data_algG_acute_2_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T19/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_acute_2_5UTR) <- c("PA14-Koordinate","As_acute_2_5UTR","Ts_acute_2_5UTR","Cs_acute_2_5UTR","Gs_acute_2_5UTR")
data_algG_acute_2_5UTR <- merge(PA14_bases_5UTR, data_algG_acute_2_5UTR, "PA14-Koordinate")
data_algG_acute_2_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T19/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_acute_2_3UTR) <- c("PA14-Koordinate","As_acute_2_3UTR","Ts_acute_2_3UTR","Cs_acute_2_3UTR","Gs_acute_2_3UTR")
data_algG_acute_2_3UTR <- merge(PA14_bases_3UTR, data_algG_acute_2_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_acute_2$coverage_acute_2 <- data_algG_acute_2$Ts_acute_2 + data_algG_acute_2$As_acute_2 + data_algG_acute_2$Gs_acute_2 + data_algG_acute_2$Cs_acute_2
data_algG_acute_2_5UTR$coverage_acute_2_5UTR <- data_algG_acute_2_5UTR$Ts_acute_2_5UTR + data_algG_acute_2_5UTR$As_acute_2_5UTR + data_algG_acute_2_5UTR$Gs_acute_2_5UTR + data_algG_acute_2_5UTR$Cs_acute_2_5UTR
data_algG_acute_2_3UTR$coverage_acute_2_3UTR <- data_algG_acute_2_3UTR$Ts_acute_2_3UTR + data_algG_acute_2_3UTR$As_acute_2_3UTR + data_algG_acute_2_3UTR$Gs_acute_2_3UTR + data_algG_acute_2_3UTR$Cs_acute_2_3UTR

### acute_3

# Input 
data_algG_acute_3 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T20/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_acute_3) <- c("PA14-Koordinate","Ts_acute_3","As_acute_3","Gs_acute_3","Cs_acute_3")
data_algG_acute_3 <- merge(PA14_bases, data_algG_acute_3, "PA14-Koordinate")
data_algG_acute_3_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T20/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_acute_3_5UTR) <- c("PA14-Koordinate","As_acute_3_5UTR","Ts_acute_3_5UTR","Cs_acute_3_5UTR","Gs_acute_3_5UTR")
data_algG_acute_3_5UTR <- merge(PA14_bases_5UTR, data_algG_acute_3_5UTR, "PA14-Koordinate")
data_algG_acute_3_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T20/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_acute_3_3UTR) <- c("PA14-Koordinate","As_acute_3_3UTR","Ts_acute_3_3UTR","Cs_acute_3_3UTR","Gs_acute_3_3UTR")
data_algG_acute_3_3UTR <- merge(PA14_bases_3UTR, data_algG_acute_3_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_acute_3$coverage_acute_3 <- data_algG_acute_3$Ts_acute_3 + data_algG_acute_3$As_acute_3 + data_algG_acute_3$Gs_acute_3 + data_algG_acute_3$Cs_acute_3
data_algG_acute_3_5UTR$coverage_acute_3_5UTR <- data_algG_acute_3_5UTR$Ts_acute_3_5UTR + data_algG_acute_3_5UTR$As_acute_3_5UTR + data_algG_acute_3_5UTR$Gs_acute_3_5UTR + data_algG_acute_3_5UTR$Cs_acute_3_5UTR
data_algG_acute_3_3UTR$coverage_acute_3_3UTR <- data_algG_acute_3_3UTR$Ts_acute_3_3UTR + data_algG_acute_3_3UTR$As_acute_3_3UTR + data_algG_acute_3_3UTR$Gs_acute_3_3UTR + data_algG_acute_3_3UTR$Cs_acute_3_3UTR

### env_1

# Input 
data_algG_env_1 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T21/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_env_1) <- c("PA14-Koordinate","Ts_env_1","As_env_1","Gs_env_1","Cs_env_1")
data_algG_env_1 <- merge(PA14_bases, data_algG_env_1, "PA14-Koordinate")
data_algG_env_1_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T21/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_env_1_5UTR) <- c("PA14-Koordinate","As_env_1_5UTR","Ts_env_1_5UTR","Cs_env_1_5UTR","Gs_env_1_5UTR")
data_algG_env_1_5UTR <- merge(PA14_bases_5UTR, data_algG_env_1_5UTR, "PA14-Koordinate")
data_algG_env_1_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T21/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_env_1_3UTR) <- c("PA14-Koordinate","As_env_1_3UTR","Ts_env_1_3UTR","Cs_env_1_3UTR","Gs_env_1_3UTR")
data_algG_env_1_3UTR <- merge(PA14_bases_3UTR, data_algG_env_1_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_env_1$coverage_env_1 <- data_algG_env_1$Ts_env_1 + data_algG_env_1$As_env_1 + data_algG_env_1$Gs_env_1 + data_algG_env_1$Cs_env_1
data_algG_env_1_5UTR$coverage_env_1_5UTR <- data_algG_env_1_5UTR$Ts_env_1_5UTR + data_algG_env_1_5UTR$As_env_1_5UTR + data_algG_env_1_5UTR$Gs_env_1_5UTR + data_algG_env_1_5UTR$Cs_env_1_5UTR
data_algG_env_1_3UTR$coverage_env_1_3UTR <- data_algG_env_1_3UTR$Ts_env_1_3UTR + data_algG_env_1_3UTR$As_env_1_3UTR + data_algG_env_1_3UTR$Gs_env_1_3UTR + data_algG_env_1_3UTR$Cs_env_1_3UTR

### env_2

# Input 
data_algG_env_2 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T22/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_env_2) <- c("PA14-Koordinate","Ts_env_2","As_env_2","Gs_env_2","Cs_env_2")
data_algG_env_2 <- merge(PA14_bases, data_algG_env_2, "PA14-Koordinate")
data_algG_env_2_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T22/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_env_2_5UTR) <- c("PA14-Koordinate","As_env_2_5UTR","Ts_env_2_5UTR","Cs_env_2_5UTR","Gs_env_2_5UTR")
data_algG_env_2_5UTR <- merge(PA14_bases_5UTR, data_algG_env_2_5UTR, "PA14-Koordinate")
data_algG_env_2_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T22/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_env_2_3UTR) <- c("PA14-Koordinate","As_env_2_3UTR","Ts_env_2_3UTR","Cs_env_2_3UTR","Gs_env_2_3UTR")
data_algG_env_2_3UTR <- merge(PA14_bases_3UTR, data_algG_env_2_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_env_2$coverage_env_2 <- data_algG_env_2$Ts_env_2 + data_algG_env_2$As_env_2 + data_algG_env_2$Gs_env_2 + data_algG_env_2$Cs_env_2
data_algG_env_2_5UTR$coverage_env_2_5UTR <- data_algG_env_2_5UTR$Ts_env_2_5UTR + data_algG_env_2_5UTR$As_env_2_5UTR + data_algG_env_2_5UTR$Gs_env_2_5UTR + data_algG_env_2_5UTR$Cs_env_2_5UTR
data_algG_env_2_3UTR$coverage_env_2_3UTR <- data_algG_env_2_3UTR$Ts_env_2_3UTR + data_algG_env_2_3UTR$As_env_2_3UTR + data_algG_env_2_3UTR$Gs_env_2_3UTR + data_algG_env_2_3UTR$Cs_env_2_3UTR

### env_3

# Input 
data_algG_env_3 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T23/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_env_3) <- c("PA14-Koordinate","Ts_env_3","As_env_3","Gs_env_3","Cs_env_3")
data_algG_env_3 <- merge(PA14_bases, data_algG_env_3, "PA14-Koordinate")
data_algG_env_3_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T23/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_env_3_5UTR) <- c("PA14-Koordinate","As_env_3_5UTR","Ts_env_3_5UTR","Cs_env_3_5UTR","Gs_env_3_5UTR")
data_algG_env_3_5UTR <- merge(PA14_bases_5UTR, data_algG_env_3_5UTR, "PA14-Koordinate")
data_algG_env_3_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T23/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_env_3_3UTR) <- c("PA14-Koordinate","As_env_3_3UTR","Ts_env_3_3UTR","Cs_env_3_3UTR","Gs_env_3_3UTR")
data_algG_env_3_3UTR <- merge(PA14_bases_3UTR, data_algG_env_3_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_env_3$coverage_env_3 <- data_algG_env_3$Ts_env_3 + data_algG_env_3$As_env_3 + data_algG_env_3$Gs_env_3 + data_algG_env_3$Cs_env_3
data_algG_env_3_5UTR$coverage_env_3_5UTR <- data_algG_env_3_5UTR$Ts_env_3_5UTR + data_algG_env_3_5UTR$As_env_3_5UTR + data_algG_env_3_5UTR$Gs_env_3_5UTR + data_algG_env_3_5UTR$Cs_env_3_5UTR
data_algG_env_3_3UTR$coverage_env_3_3UTR <- data_algG_env_3_3UTR$Ts_env_3_3UTR + data_algG_env_3_3UTR$As_env_3_3UTR + data_algG_env_3_3UTR$Gs_env_3_3UTR + data_algG_env_3_3UTR$Cs_env_3_3UTR

### env_4

# Input 
data_algG_env_4 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T24/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_env_4) <- c("PA14-Koordinate","Ts_env_4","As_env_4","Gs_env_4","Cs_env_4")
data_algG_env_4 <- merge(PA14_bases, data_algG_env_4, "PA14-Koordinate")
data_algG_env_4_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T24/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_env_4_5UTR) <- c("PA14-Koordinate","As_env_4_5UTR","Ts_env_4_5UTR","Cs_env_4_5UTR","Gs_env_4_5UTR")
data_algG_env_4_5UTR <- merge(PA14_bases_5UTR, data_algG_env_4_5UTR, "PA14-Koordinate")
data_algG_env_4_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T24/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_env_4_3UTR) <- c("PA14-Koordinate","As_env_4_3UTR","Ts_env_4_3UTR","Cs_env_4_3UTR","Gs_env_4_3UTR")
data_algG_env_4_3UTR <- merge(PA14_bases_3UTR, data_algG_env_4_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_env_4$coverage_env_4 <- data_algG_env_4$Ts_env_4 + data_algG_env_4$As_env_4 + data_algG_env_4$Gs_env_4 + data_algG_env_4$Cs_env_4
data_algG_env_4_5UTR$coverage_env_4_5UTR <- data_algG_env_4_5UTR$Ts_env_4_5UTR + data_algG_env_4_5UTR$As_env_4_5UTR + data_algG_env_4_5UTR$Gs_env_4_5UTR + data_algG_env_4_5UTR$Cs_env_4_5UTR
data_algG_env_4_3UTR$coverage_env_4_3UTR <- data_algG_env_4_3UTR$Ts_env_4_3UTR + data_algG_env_4_3UTR$As_env_4_3UTR + data_algG_env_4_3UTR$Gs_env_4_3UTR + data_algG_env_4_3UTR$Cs_env_4_3UTR

### env_5

# Input 
data_algG_env_5 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T25/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_env_5) <- c("PA14-Koordinate","Ts_env_5","As_env_5","Gs_env_5","Cs_env_5")
data_algG_env_5 <- merge(PA14_bases, data_algG_env_5, "PA14-Koordinate")
data_algG_env_5_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T25/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_env_5_5UTR) <- c("PA14-Koordinate","As_env_5_5UTR","Ts_env_5_5UTR","Cs_env_5_5UTR","Gs_env_5_5UTR")
data_algG_env_5_5UTR <- merge(PA14_bases_5UTR, data_algG_env_5_5UTR, "PA14-Koordinate")
data_algG_env_5_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T25/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_env_5_3UTR) <- c("PA14-Koordinate","As_env_5_3UTR","Ts_env_5_3UTR","Cs_env_5_3UTR","Gs_env_5_3UTR")
data_algG_env_5_3UTR <- merge(PA14_bases_3UTR, data_algG_env_5_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_env_5$coverage_env_5 <- data_algG_env_5$Ts_env_5 + data_algG_env_5$As_env_5 + data_algG_env_5$Gs_env_5 + data_algG_env_5$Cs_env_5
data_algG_env_5_5UTR$coverage_env_5_5UTR <- data_algG_env_5_5UTR$Ts_env_5_5UTR + data_algG_env_5_5UTR$As_env_5_5UTR + data_algG_env_5_5UTR$Gs_env_5_5UTR + data_algG_env_5_5UTR$Cs_env_5_5UTR
data_algG_env_5_3UTR$coverage_env_5_3UTR <- data_algG_env_5_3UTR$Ts_env_5_3UTR + data_algG_env_5_3UTR$As_env_5_3UTR + data_algG_env_5_3UTR$Gs_env_5_3UTR + data_algG_env_5_3UTR$Cs_env_5_3UTR

### COPD

# Input 
data_algG_COPD <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T26/Results/algG/Results_algG.txt", header = TRUE)
names(data_algG_COPD) <- c("PA14-Koordinate","Ts_COPD","As_COPD","Gs_COPD","Cs_COPD")
data_algG_COPD <- merge(PA14_bases, data_algG_COPD, "PA14-Koordinate")
data_algG_COPD_5UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T26/Results/algG_5UTR/Results_algG_5UTR.txt", header = TRUE)
names(data_algG_COPD_5UTR) <- c("PA14-Koordinate","As_COPD_5UTR","Ts_COPD_5UTR","Cs_COPD_5UTR","Gs_COPD_5UTR")
data_algG_COPD_5UTR <- merge(PA14_bases_5UTR, data_algG_COPD_5UTR, "PA14-Koordinate")
data_algG_COPD_3UTR <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T26/Results/algG_3UTR/Results_algG_3UTR.txt", header = TRUE)
names(data_algG_COPD_3UTR) <- c("PA14-Koordinate","As_COPD_3UTR","Ts_COPD_3UTR","Cs_COPD_3UTR","Gs_COPD_3UTR")
data_algG_COPD_3UTR <- merge(PA14_bases_3UTR, data_algG_COPD_3UTR, "PA14-Koordinate")

# Add coverage

data_algG_COPD$coverage_COPD <- data_algG_COPD$Ts_COPD + data_algG_COPD$As_COPD + data_algG_COPD$Gs_COPD + data_algG_COPD$Cs_COPD
data_algG_COPD_5UTR$coverage_COPD_5UTR <- data_algG_COPD_5UTR$Ts_COPD_5UTR + data_algG_COPD_5UTR$As_COPD_5UTR + data_algG_COPD_5UTR$Gs_COPD_5UTR + data_algG_COPD_5UTR$Cs_COPD_5UTR
data_algG_COPD_3UTR$coverage_COPD_3UTR <- data_algG_COPD_3UTR$Ts_COPD_3UTR + data_algG_COPD_3UTR$As_COPD_3UTR + data_algG_COPD_3UTR$Gs_COPD_3UTR + data_algG_COPD_3UTR$Cs_COPD_3UTR
















# Coverage plot


a <- ggplot(data = data_algG_CF_1, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_1)) +
  geom_point()+
  geom_point(data = data_algG_CF_1_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_1_5UTR),color="blue")+
  geom_point(data = data_algG_CF_1_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_1_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

b <- ggplot(data = data_algG_CF_2, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_2)) +
  geom_point()+
  geom_point(data = data_algG_CF_2_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_2_5UTR),color="blue")+
  geom_point(data = data_algG_CF_2_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_2_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

c <- ggplot(data = data_algG_CF_3, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_3)) +
  geom_point()+
  geom_point(data = data_algG_CF_3_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_3_5UTR),color="blue")+
  geom_point(data = data_algG_CF_3_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_3_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

d <- ggplot(data = data_algG_CF_4, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_4)) +
  geom_point()+
  geom_point(data = data_algG_CF_4_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_4_5UTR),color="blue")+
  geom_point(data = data_algG_CF_4_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_4_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

e <- ggplot(data = data_algG_CF_5, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_5)) +
  geom_point()+
  geom_point(data = data_algG_CF_5_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_5_5UTR),color="blue")+
  geom_point(data = data_algG_CF_5_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_5_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

f <- ggplot(data = data_algG_CF_6, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_6)) +
  geom_point()+
  geom_point(data = data_algG_CF_6_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_6_5UTR),color="blue")+
  geom_point(data = data_algG_CF_6_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_6_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

g <- ggplot(data = data_algG_CF_7, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_7)) +
  geom_point()+
  geom_point(data = data_algG_CF_7_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_7_5UTR),color="blue")+
  geom_point(data = data_algG_CF_7_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_7_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

h <- ggplot(data = data_algG_CF_8, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_8)) +
  geom_point()+
  geom_point(data = data_algG_CF_8_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_8_5UTR),color="blue")+
  geom_point(data = data_algG_CF_8_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_8_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

i <- ggplot(data = data_algG_CF_9, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_9)) +
  geom_point()+
  geom_point(data = data_algG_CF_9_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_9_5UTR),color="blue")+
  geom_point(data = data_algG_CF_9_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_9_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

j <- ggplot(data = data_algG_CF_10, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_10)) +
  geom_point()+
  geom_point(data = data_algG_CF_10_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_10_5UTR),color="blue")+
  geom_point(data = data_algG_CF_10_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_10_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

k <- ggplot(data = data_algG_CF_11, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_11)) +
  geom_point()+
  geom_point(data = data_algG_CF_11_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_11_5UTR),color="blue")+
  geom_point(data = data_algG_CF_11_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_11_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

l <- ggplot(data = data_algG_CF_12, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_12)) +
  geom_point()+
  geom_point(data = data_algG_CF_12_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_12_5UTR),color="blue")+
  geom_point(data = data_algG_CF_12_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_12_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

m <- ggplot(data = data_algG_CF_13, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_13)) +
  geom_point()+
  geom_point(data = data_algG_CF_13_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_13_5UTR),color="blue")+
  geom_point(data = data_algG_CF_13_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_13_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

n <- ggplot(data = data_algG_CF_14, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_14)) +
  geom_point()+
  geom_point(data = data_algG_CF_14_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_14_5UTR),color="blue")+
  geom_point(data = data_algG_CF_14_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_14_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

o <- ggplot(data = data_algG_CF_15, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_15)) +
  geom_point()+
  geom_point(data = data_algG_CF_15_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_15_5UTR),color="blue")+
  geom_point(data = data_algG_CF_15_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_15_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

p <- ggplot(data = data_algG_CF_16, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_16)) +
  geom_point()+
  geom_point(data = data_algG_CF_16_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_16_5UTR),color="blue")+
  geom_point(data = data_algG_CF_16_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_16_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

q <- ggplot(data = data_algG_CF_17, mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_17)) +
  geom_point()+
  geom_point(data = data_algG_CF_17_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_17_5UTR),color="blue")+
  geom_point(data = data_algG_CF_17_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_CF_17_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

r <- ggplot(data = data_algG_acute_1, mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_1)) +
  geom_point()+
  geom_point(data = data_algG_acute_1_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_1_5UTR),color="blue")+
  geom_point(data = data_algG_acute_1_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_1_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

s <- ggplot(data = data_algG_acute_2, mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_2)) +
  geom_point()+
  geom_point(data = data_algG_acute_2_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_2_5UTR),color="blue")+
  geom_point(data = data_algG_acute_2_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_2_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

t <- ggplot(data = data_algG_acute_3, mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_3)) +
  geom_point()+
  geom_point(data = data_algG_acute_3_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_3_5UTR),color="blue")+
  geom_point(data = data_algG_acute_3_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_acute_3_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

u <- ggplot(data = data_algG_env_1, mapping = aes(x = `PA14-Koordinate`, y = coverage_env_1)) +
  geom_point()+
  geom_point(data = data_algG_env_1_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_1_5UTR),color="blue")+
  geom_point(data = data_algG_env_1_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_1_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

v <- ggplot(data = data_algG_env_2, mapping = aes(x = `PA14-Koordinate`, y = coverage_env_2)) +
  geom_point()+
  geom_point(data = data_algG_env_2_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_2_5UTR),color="blue")+
  geom_point(data = data_algG_env_2_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_2_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

w <- ggplot(data = data_algG_env_3, mapping = aes(x = `PA14-Koordinate`, y = coverage_env_3)) +
  geom_point()+
  geom_point(data = data_algG_env_3_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_3_5UTR),color="blue")+
  geom_point(data = data_algG_env_3_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_3_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

x <- ggplot(data = data_algG_env_4, mapping = aes(x = `PA14-Koordinate`, y = coverage_env_4)) +
  geom_point()+
  geom_point(data = data_algG_env_4_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_4_5UTR),color="blue")+
  geom_point(data = data_algG_env_4_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_4_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

y <- ggplot(data = data_algG_env_5, mapping = aes(x = `PA14-Koordinate`, y = coverage_env_5)) +
  geom_point()+
  geom_point(data = data_algG_env_5_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_5_5UTR),color="blue")+
  geom_point(data = data_algG_env_5_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_env_5_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

z <- ggplot(data = data_algG_COPD, mapping = aes(x = `PA14-Koordinate`, y = coverage_COPD)) +
  geom_point()+
  geom_point(data = data_algG_COPD_5UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_COPD_5UTR),color="blue")+
  geom_point(data = data_algG_COPD_3UTR , mapping = aes(x = `PA14-Koordinate`, y = coverage_COPD_3UTR),color="red")+
  ylim(0,12000)+
  geom_hline(yintercept=isolates_COPD*10,color="red")+
  geom_hline(yintercept=1000,color="blue")

png("algG_cov.png", width=2400, height=2880)
grid.arrange(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,ncol=5)
dev.off()


########### CF_1   ##############

# Calculate percentage for each base

data_algG_CF_1$percent_Ts_CF_1 <- data_algG_CF_1$Ts_CF_1/data_algG_CF_1$coverage_CF_1
data_algG_CF_1$percent_As_CF_1 <- data_algG_CF_1$As_CF_1/data_algG_CF_1$coverage_CF_1
data_algG_CF_1$percent_Gs_CF_1 <- data_algG_CF_1$Gs_CF_1/data_algG_CF_1$coverage_CF_1
data_algG_CF_1$percent_Cs_CF_1 <- data_algG_CF_1$Cs_CF_1/data_algG_CF_1$coverage_CF_1


# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_1$highest_percentage <- 0
data_algG_CF_1$highest_percentage <- apply(data_algG_CF_1[, 8:11], 1, max)
data_algG_CF_1$second_percentage <- apply(data_algG_CF_1[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_1[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_1_filtered <- subset(data_algG_CF_1, second_percentage >= 0.025) 
data_algG_CF_1_filtered$main_base_CF_1 <- colnames(data_algG_CF_1_filtered[, 8:11])[apply(data_algG_CF_1_filtered[, 8:11],1,which.max)]

data_algG_CF_1_filtered_4percent <- subset(data_algG_CF_1, second_percentage >= 0.04) 
data_algG_CF_1_filtered_4percent$main_base_CF_1 <- colnames(data_algG_CF_1_filtered_4percent[, 8:11])[apply(data_algG_CF_1_filtered_4percent[, 8:11],1,which.max)]


########### CF_2   ##############

# Calculate percentage for each base

data_algG_CF_2$percent_Ts_CF_2 <- data_algG_CF_2$Ts_CF_2/data_algG_CF_2$coverage_CF_2
data_algG_CF_2$percent_As_CF_2 <- data_algG_CF_2$As_CF_2/data_algG_CF_2$coverage_CF_2
data_algG_CF_2$percent_Gs_CF_2 <- data_algG_CF_2$Gs_CF_2/data_algG_CF_2$coverage_CF_2
data_algG_CF_2$percent_Cs_CF_2 <- data_algG_CF_2$Cs_CF_2/data_algG_CF_2$coverage_CF_2




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_2$highest_percentage <- 0
data_algG_CF_2$highest_percentage <- apply(data_algG_CF_2[, 8:11], 1, max)
data_algG_CF_2$second_percentage <- apply(data_algG_CF_2[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_2[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_2_filtered <- subset(data_algG_CF_2, second_percentage >= 0.025) 
data_algG_CF_2_filtered$main_base_CF_2 <- colnames(data_algG_CF_2_filtered[, 8:11])[apply(data_algG_CF_2_filtered[, 8:11],1,which.max)]

data_algG_CF_2_filtered_4percent <- subset(data_algG_CF_2, second_percentage >= 0.04) 
data_algG_CF_2_filtered_4percent$main_base_CF_2 <- colnames(data_algG_CF_2_filtered_4percent[, 8:11])[apply(data_algG_CF_2_filtered_4percent[, 8:11],1,which.max)]

########### CF_3   ##############

# Calculate percentage for each base

data_algG_CF_3$percent_Ts_CF_3 <- data_algG_CF_3$Ts_CF_3/data_algG_CF_3$coverage_CF_3
data_algG_CF_3$percent_As_CF_3 <- data_algG_CF_3$As_CF_3/data_algG_CF_3$coverage_CF_3
data_algG_CF_3$percent_Gs_CF_3 <- data_algG_CF_3$Gs_CF_3/data_algG_CF_3$coverage_CF_3
data_algG_CF_3$percent_Cs_CF_3 <- data_algG_CF_3$Cs_CF_3/data_algG_CF_3$coverage_CF_3




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_3$highest_percentage <- 0
data_algG_CF_3$highest_percentage <- apply(data_algG_CF_3[, 8:11], 1, max)
data_algG_CF_3$second_percentage <- apply(data_algG_CF_3[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_3[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_3_filtered <- subset(data_algG_CF_3, second_percentage >= 0.025) 
data_algG_CF_3_filtered$main_base_CF_3 <- colnames(data_algG_CF_3_filtered[, 8:11])[apply(data_algG_CF_3_filtered[, 8:11],1,which.max)]

data_algG_CF_3_filtered_4percent <- subset(data_algG_CF_3, second_percentage >= 0.04) 
data_algG_CF_3_filtered_4percent$main_base_CF_3 <- colnames(data_algG_CF_3_filtered_4percent[, 8:11])[apply(data_algG_CF_3_filtered_4percent[, 8:11],1,which.max)]


########### CF_4   ##############

# Calculate percentage for each base

data_algG_CF_4$percent_Ts_CF_4 <- data_algG_CF_4$Ts_CF_4/data_algG_CF_4$coverage_CF_4
data_algG_CF_4$percent_As_CF_4 <- data_algG_CF_4$As_CF_4/data_algG_CF_4$coverage_CF_4
data_algG_CF_4$percent_Gs_CF_4 <- data_algG_CF_4$Gs_CF_4/data_algG_CF_4$coverage_CF_4
data_algG_CF_4$percent_Cs_CF_4 <- data_algG_CF_4$Cs_CF_4/data_algG_CF_4$coverage_CF_4




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_4$highest_percentage <- 0
data_algG_CF_4$highest_percentage <- apply(data_algG_CF_4[, 8:11], 1, max)
data_algG_CF_4$second_percentage <- apply(data_algG_CF_4[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_4[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_4_filtered <- subset(data_algG_CF_4, second_percentage >= 0.025) 
data_algG_CF_4_filtered$main_base_CF_4 <- colnames(data_algG_CF_4_filtered[, 8:11])[apply(data_algG_CF_4_filtered[, 8:11],1,which.max)]

data_algG_CF_4_filtered_4percent <- subset(data_algG_CF_4, second_percentage >= 0.04) 
data_algG_CF_4_filtered_4percent$main_base_CF_4 <- colnames(data_algG_CF_4_filtered_4percent[, 8:11])[apply(data_algG_CF_4_filtered_4percent[, 8:11],1,which.max)]



save(data_algG_CF_4_filtered,data_algG_CF_4_filtered_4percent,file="algG.Rdata")


########### CF_5   ##############

# Calculate percentage for each base

data_algG_CF_5$percent_Ts_CF_5 <- data_algG_CF_5$Ts_CF_5/data_algG_CF_5$coverage_CF_5
data_algG_CF_5$percent_As_CF_5 <- data_algG_CF_5$As_CF_5/data_algG_CF_5$coverage_CF_5
data_algG_CF_5$percent_Gs_CF_5 <- data_algG_CF_5$Gs_CF_5/data_algG_CF_5$coverage_CF_5
data_algG_CF_5$percent_Cs_CF_5 <- data_algG_CF_5$Cs_CF_5/data_algG_CF_5$coverage_CF_5




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_5$highest_percentage <- 0
data_algG_CF_5$highest_percentage <- apply(data_algG_CF_5[, 8:11], 1, max)
data_algG_CF_5$second_percentage <- apply(data_algG_CF_5[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_5[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_5_filtered <- subset(data_algG_CF_5, second_percentage >= 0.025) 
data_algG_CF_5_filtered$main_base_CF_5 <- colnames(data_algG_CF_5_filtered[, 8:11])[apply(data_algG_CF_5_filtered[, 8:11],1,which.max)]

data_algG_CF_5_filtered_4percent <- subset(data_algG_CF_5, second_percentage >= 0.04) 
data_algG_CF_5_filtered_4percent$main_base_CF_5 <- colnames(data_algG_CF_5_filtered_4percent[, 8:11])[apply(data_algG_CF_5_filtered_4percent[, 8:11],1,which.max)]

########### CF_6   ##############

# Calculate percentage for each base

data_algG_CF_6$percent_Ts_CF_6 <- data_algG_CF_6$Ts_CF_6/data_algG_CF_6$coverage_CF_6
data_algG_CF_6$percent_As_CF_6 <- data_algG_CF_6$As_CF_6/data_algG_CF_6$coverage_CF_6
data_algG_CF_6$percent_Gs_CF_6 <- data_algG_CF_6$Gs_CF_6/data_algG_CF_6$coverage_CF_6
data_algG_CF_6$percent_Cs_CF_6 <- data_algG_CF_6$Cs_CF_6/data_algG_CF_6$coverage_CF_6




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_6$highest_percentage <- 0
data_algG_CF_6$highest_percentage <- apply(data_algG_CF_6[, 8:11], 1, max)
data_algG_CF_6$second_percentage <- apply(data_algG_CF_6[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_6[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_6_filtered <- subset(data_algG_CF_6, second_percentage >= 0.025) 
data_algG_CF_6_filtered$main_base_CF_6 <- colnames(data_algG_CF_6_filtered[, 8:11])[apply(data_algG_CF_6_filtered[, 8:11],1,which.max)]

data_algG_CF_6_filtered_4percent <- subset(data_algG_CF_6, second_percentage >= 0.04) 
data_algG_CF_6_filtered_4percent$main_base_CF_6 <- colnames(data_algG_CF_6_filtered_4percent[, 8:11])[apply(data_algG_CF_6_filtered_4percent[, 8:11],1,which.max)]


########### CF_7   ##############

# Calculate percentage for each base

data_algG_CF_7$percent_Ts_CF_7 <- data_algG_CF_7$Ts_CF_7/data_algG_CF_7$coverage_CF_7
data_algG_CF_7$percent_As_CF_7 <- data_algG_CF_7$As_CF_7/data_algG_CF_7$coverage_CF_7
data_algG_CF_7$percent_Gs_CF_7 <- data_algG_CF_7$Gs_CF_7/data_algG_CF_7$coverage_CF_7
data_algG_CF_7$percent_Cs_CF_7 <- data_algG_CF_7$Cs_CF_7/data_algG_CF_7$coverage_CF_7




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_7$highest_percentage <- 0
data_algG_CF_7$highest_percentage <- apply(data_algG_CF_7[, 8:11], 1, max)
data_algG_CF_7$second_percentage <- apply(data_algG_CF_7[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_7[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_7_filtered <- subset(data_algG_CF_7, second_percentage >= 0.025) 
data_algG_CF_7_filtered$main_base_CF_7 <- colnames(data_algG_CF_7_filtered[, 8:11])[apply(data_algG_CF_7_filtered[, 8:11],1,which.max)]

data_algG_CF_7_filtered_4percent <- subset(data_algG_CF_7, second_percentage >= 0.04) 
data_algG_CF_7_filtered_4percent$main_base_CF_7 <- colnames(data_algG_CF_7_filtered_4percent[, 8:11])[apply(data_algG_CF_7_filtered_4percent[, 8:11],1,which.max)]


########### CF_8   ##############

# Calculate percentage for each base

data_algG_CF_8$percent_Ts_CF_8 <- data_algG_CF_8$Ts_CF_8/data_algG_CF_8$coverage_CF_8
data_algG_CF_8$percent_As_CF_8 <- data_algG_CF_8$As_CF_8/data_algG_CF_8$coverage_CF_8
data_algG_CF_8$percent_Gs_CF_8 <- data_algG_CF_8$Gs_CF_8/data_algG_CF_8$coverage_CF_8
data_algG_CF_8$percent_Cs_CF_8 <- data_algG_CF_8$Cs_CF_8/data_algG_CF_8$coverage_CF_8




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_8$highest_percentage <- 0
data_algG_CF_8$highest_percentage <- apply(data_algG_CF_8[, 8:11], 1, max)
data_algG_CF_8$second_percentage <- apply(data_algG_CF_8[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_8[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_8_filtered <- subset(data_algG_CF_8, second_percentage >= 0.025) 
data_algG_CF_8_filtered$main_base_CF_8 <- colnames(data_algG_CF_8_filtered[, 8:11])[apply(data_algG_CF_8_filtered[, 8:11],1,which.max)]

data_algG_CF_8_filtered_4percent <- subset(data_algG_CF_8, second_percentage >= 0.04) 
data_algG_CF_8_filtered_4percent$main_base_CF_8 <- colnames(data_algG_CF_8_filtered_4percent[, 8:11])[apply(data_algG_CF_8_filtered_4percent[, 8:11],1,which.max)]

########### CF_9   ##############

# Calculate percentage for each base

data_algG_CF_9$percent_Ts_CF_9 <- data_algG_CF_9$Ts_CF_9/data_algG_CF_9$coverage_CF_9
data_algG_CF_9$percent_As_CF_9 <- data_algG_CF_9$As_CF_9/data_algG_CF_9$coverage_CF_9
data_algG_CF_9$percent_Gs_CF_9 <- data_algG_CF_9$Gs_CF_9/data_algG_CF_9$coverage_CF_9
data_algG_CF_9$percent_Cs_CF_9 <- data_algG_CF_9$Cs_CF_9/data_algG_CF_9$coverage_CF_9




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_9$highest_percentage <- 0
data_algG_CF_9$highest_percentage <- apply(data_algG_CF_9[, 8:11], 1, max)
data_algG_CF_9$second_percentage <- apply(data_algG_CF_9[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_9[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_9_filtered <- subset(data_algG_CF_9, second_percentage >= 0.025) 
data_algG_CF_9_filtered$main_base_CF_9 <- colnames(data_algG_CF_9_filtered[, 8:11])[apply(data_algG_CF_9_filtered[, 8:11],1,which.max)]

data_algG_CF_9_filtered_4percent <- subset(data_algG_CF_9, second_percentage >= 0.04) 
data_algG_CF_9_filtered_4percent$main_base_CF_9 <- colnames(data_algG_CF_9_filtered_4percent[, 8:11])[apply(data_algG_CF_9_filtered_4percent[, 8:11],1,which.max)]


########### CF_10   ##############

# Calculate percentage for each base

data_algG_CF_10$percent_Ts_CF_10 <- data_algG_CF_10$Ts_CF_10/data_algG_CF_10$coverage_CF_10
data_algG_CF_10$percent_As_CF_10 <- data_algG_CF_10$As_CF_10/data_algG_CF_10$coverage_CF_10
data_algG_CF_10$percent_Gs_CF_10 <- data_algG_CF_10$Gs_CF_10/data_algG_CF_10$coverage_CF_10
data_algG_CF_10$percent_Cs_CF_10 <- data_algG_CF_10$Cs_CF_10/data_algG_CF_10$coverage_CF_10




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_10$highest_percentage <- 0
data_algG_CF_10$highest_percentage <- apply(data_algG_CF_10[, 8:11], 1, max)
data_algG_CF_10$second_percentage <- apply(data_algG_CF_10[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_10[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_10_filtered <- subset(data_algG_CF_10, second_percentage >= 0.025) 
data_algG_CF_10_filtered$main_base_CF_10 <- colnames(data_algG_CF_10_filtered[, 8:11])[apply(data_algG_CF_10_filtered[, 8:11],1,which.max)]

data_algG_CF_10_filtered_4percent <- subset(data_algG_CF_10, second_percentage >= 0.04) 
data_algG_CF_10_filtered_4percent$main_base_CF_10 <- colnames(data_algG_CF_10_filtered_4percent[, 8:11])[apply(data_algG_CF_10_filtered_4percent[, 8:11],1,which.max)]

########### CF_11   ##############

# Calculate percentage for each base

data_algG_CF_11$percent_Ts_CF_11 <- data_algG_CF_11$Ts_CF_11/data_algG_CF_11$coverage_CF_11
data_algG_CF_11$percent_As_CF_11 <- data_algG_CF_11$As_CF_11/data_algG_CF_11$coverage_CF_11
data_algG_CF_11$percent_Gs_CF_11 <- data_algG_CF_11$Gs_CF_11/data_algG_CF_11$coverage_CF_11
data_algG_CF_11$percent_Cs_CF_11 <- data_algG_CF_11$Cs_CF_11/data_algG_CF_11$coverage_CF_11




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_11$highest_percentage <- 0
data_algG_CF_11$highest_percentage <- apply(data_algG_CF_11[, 8:11], 1, max)
data_algG_CF_11$second_percentage <- apply(data_algG_CF_11[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_11[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_11_filtered <- subset(data_algG_CF_11, second_percentage >= 0.025) 
data_algG_CF_11_filtered$main_base_CF_11 <- colnames(data_algG_CF_11_filtered[, 8:11])[apply(data_algG_CF_11_filtered[, 8:11],1,which.max)]

data_algG_CF_11_filtered_4percent <- subset(data_algG_CF_11, second_percentage >= 0.04) 
data_algG_CF_11_filtered_4percent$main_base_CF_11 <- colnames(data_algG_CF_11_filtered_4percent[, 8:11])[apply(data_algG_CF_11_filtered_4percent[, 8:11],1,which.max)]

########### CF_12   ##############

# Calculate percentage for each base

data_algG_CF_12$percent_Ts_CF_12 <- data_algG_CF_12$Ts_CF_12/data_algG_CF_12$coverage_CF_12
data_algG_CF_12$percent_As_CF_12 <- data_algG_CF_12$As_CF_12/data_algG_CF_12$coverage_CF_12
data_algG_CF_12$percent_Gs_CF_12 <- data_algG_CF_12$Gs_CF_12/data_algG_CF_12$coverage_CF_12
data_algG_CF_12$percent_Cs_CF_12 <- data_algG_CF_12$Cs_CF_12/data_algG_CF_12$coverage_CF_12




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_12$highest_percentage <- 0
data_algG_CF_12$highest_percentage <- apply(data_algG_CF_12[, 8:11], 1, max)
data_algG_CF_12$second_percentage <- apply(data_algG_CF_12[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_12[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_12_filtered <- subset(data_algG_CF_12, second_percentage >= 0.025) 
data_algG_CF_12_filtered$main_base_CF_12 <- colnames(data_algG_CF_12_filtered[, 8:11])[apply(data_algG_CF_12_filtered[, 8:11],1,which.max)]

data_algG_CF_12_filtered_4percent <- subset(data_algG_CF_12, second_percentage >= 0.04) 
data_algG_CF_12_filtered_4percent$main_base_CF_12 <- colnames(data_algG_CF_12_filtered_4percent[, 8:11])[apply(data_algG_CF_12_filtered_4percent[, 8:11],1,which.max)]

########### CF_13   ##############

# Calculate percentage for each base

data_algG_CF_13$percent_Ts_CF_13 <- data_algG_CF_13$Ts_CF_13/data_algG_CF_13$coverage_CF_13
data_algG_CF_13$percent_As_CF_13 <- data_algG_CF_13$As_CF_13/data_algG_CF_13$coverage_CF_13
data_algG_CF_13$percent_Gs_CF_13 <- data_algG_CF_13$Gs_CF_13/data_algG_CF_13$coverage_CF_13
data_algG_CF_13$percent_Cs_CF_13 <- data_algG_CF_13$Cs_CF_13/data_algG_CF_13$coverage_CF_13




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_13$highest_percentage <- 0
data_algG_CF_13$highest_percentage <- apply(data_algG_CF_13[, 8:11], 1, max)
data_algG_CF_13$second_percentage <- apply(data_algG_CF_13[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_13[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_13_filtered <- subset(data_algG_CF_13, second_percentage >= 0.025) 
data_algG_CF_13_filtered$main_base_CF_13 <- colnames(data_algG_CF_13_filtered[, 8:11])[apply(data_algG_CF_13_filtered[, 8:11],1,which.max)]

data_algG_CF_13_filtered_4percent <- subset(data_algG_CF_13, second_percentage >= 0.04) 
data_algG_CF_13_filtered_4percent$main_base_CF_13 <- colnames(data_algG_CF_13_filtered_4percent[, 8:11])[apply(data_algG_CF_13_filtered_4percent[, 8:11],1,which.max)]

########### CF_14   ##############

# Calculate percentage for each base

data_algG_CF_14$percent_Ts_CF_14 <- data_algG_CF_14$Ts_CF_14/data_algG_CF_14$coverage_CF_14
data_algG_CF_14$percent_As_CF_14 <- data_algG_CF_14$As_CF_14/data_algG_CF_14$coverage_CF_14
data_algG_CF_14$percent_Gs_CF_14 <- data_algG_CF_14$Gs_CF_14/data_algG_CF_14$coverage_CF_14
data_algG_CF_14$percent_Cs_CF_14 <- data_algG_CF_14$Cs_CF_14/data_algG_CF_14$coverage_CF_14




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_14$highest_percentage <- 0
data_algG_CF_14$highest_percentage <- apply(data_algG_CF_14[, 8:11], 1, max)
data_algG_CF_14$second_percentage <- apply(data_algG_CF_14[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_14[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_14_filtered <- subset(data_algG_CF_14, second_percentage >= 0.025) 
data_algG_CF_14_filtered$main_base_CF_14 <- colnames(data_algG_CF_14_filtered[, 8:11])[apply(data_algG_CF_14_filtered[, 8:11],1,which.max)]

data_algG_CF_14_filtered_4percent <- subset(data_algG_CF_14, second_percentage >= 0.04) 
data_algG_CF_14_filtered_4percent$main_base_CF_14 <- colnames(data_algG_CF_14_filtered_4percent[, 8:11])[apply(data_algG_CF_14_filtered_4percent[, 8:11],1,which.max)]

########### CF_15   ##############

# Calculate percentage for each base

data_algG_CF_15$percent_Ts_CF_15 <- data_algG_CF_15$Ts_CF_15/data_algG_CF_15$coverage_CF_15
data_algG_CF_15$percent_As_CF_15 <- data_algG_CF_15$As_CF_15/data_algG_CF_15$coverage_CF_15
data_algG_CF_15$percent_Gs_CF_15 <- data_algG_CF_15$Gs_CF_15/data_algG_CF_15$coverage_CF_15
data_algG_CF_15$percent_Cs_CF_15 <- data_algG_CF_15$Cs_CF_15/data_algG_CF_15$coverage_CF_15




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_15$highest_percentage <- 0
data_algG_CF_15$highest_percentage <- apply(data_algG_CF_15[, 8:11], 1, max)
data_algG_CF_15$second_percentage <- apply(data_algG_CF_15[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_15[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_15_filtered <- subset(data_algG_CF_15, second_percentage >= 0.025) 
data_algG_CF_15_filtered$main_base_CF_15 <- colnames(data_algG_CF_15_filtered[, 8:11])[apply(data_algG_CF_15_filtered[, 8:11],1,which.max)]

data_algG_CF_15_filtered_4percent <- subset(data_algG_CF_15, second_percentage >= 0.04) 
data_algG_CF_15_filtered_4percent$main_base_CF_15 <- colnames(data_algG_CF_15_filtered_4percent[, 8:11])[apply(data_algG_CF_15_filtered_4percent[, 8:11],1,which.max)]

########### CF_16   ##############

# Calculate percentage for each base

data_algG_CF_16$percent_Ts_CF_16 <- data_algG_CF_16$Ts_CF_16/data_algG_CF_16$coverage_CF_16
data_algG_CF_16$percent_As_CF_16 <- data_algG_CF_16$As_CF_16/data_algG_CF_16$coverage_CF_16
data_algG_CF_16$percent_Gs_CF_16 <- data_algG_CF_16$Gs_CF_16/data_algG_CF_16$coverage_CF_16
data_algG_CF_16$percent_Cs_CF_16 <- data_algG_CF_16$Cs_CF_16/data_algG_CF_16$coverage_CF_16




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_16$highest_percentage <- 0
data_algG_CF_16$highest_percentage <- apply(data_algG_CF_16[, 8:11], 1, max)
data_algG_CF_16$second_percentage <- apply(data_algG_CF_16[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_16[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_16_filtered <- subset(data_algG_CF_16, second_percentage >= 0.025) 
data_algG_CF_16_filtered$main_base_CF_16 <- colnames(data_algG_CF_16_filtered[, 8:11])[apply(data_algG_CF_16_filtered[, 8:11],1,which.max)]

data_algG_CF_16_filtered_4percent <- subset(data_algG_CF_16, second_percentage >= 0.04) 
data_algG_CF_16_filtered_4percent$main_base_CF_16 <- colnames(data_algG_CF_16_filtered_4percent[, 8:11])[apply(data_algG_CF_16_filtered_4percent[, 8:11],1,which.max)]

########### CF_17   ##############

# Calculate percentage for each base

data_algG_CF_17$percent_Ts_CF_17 <- data_algG_CF_17$Ts_CF_17/data_algG_CF_17$coverage_CF_17
data_algG_CF_17$percent_As_CF_17 <- data_algG_CF_17$As_CF_17/data_algG_CF_17$coverage_CF_17
data_algG_CF_17$percent_Gs_CF_17 <- data_algG_CF_17$Gs_CF_17/data_algG_CF_17$coverage_CF_17
data_algG_CF_17$percent_Cs_CF_17 <- data_algG_CF_17$Cs_CF_17/data_algG_CF_17$coverage_CF_17




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_CF_17$highest_percentage <- 0
data_algG_CF_17$highest_percentage <- apply(data_algG_CF_17[, 8:11], 1, max)
data_algG_CF_17$second_percentage <- apply(data_algG_CF_17[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_17[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_17_filtered <- subset(data_algG_CF_17, second_percentage >= 0.025) 
data_algG_CF_17_filtered$main_base_CF_17 <- colnames(data_algG_CF_17_filtered[, 8:11])[apply(data_algG_CF_17_filtered[, 8:11],1,which.max)]

data_algG_CF_17_filtered_4percent <- subset(data_algG_CF_17, second_percentage >= 0.04) 
data_algG_CF_17_filtered_4percent$main_base_CF_17 <- colnames(data_algG_CF_17_filtered_4percent[, 8:11])[apply(data_algG_CF_17_filtered_4percent[, 8:11],1,which.max)]


########### acute_1   ##############

# Calculate percentage for each base

data_algG_acute_1$percent_Ts_acute_1 <- data_algG_acute_1$Ts_acute_1/data_algG_acute_1$coverage_acute_1
data_algG_acute_1$percent_As_acute_1 <- data_algG_acute_1$As_acute_1/data_algG_acute_1$coverage_acute_1
data_algG_acute_1$percent_Gs_acute_1 <- data_algG_acute_1$Gs_acute_1/data_algG_acute_1$coverage_acute_1
data_algG_acute_1$percent_Cs_acute_1 <- data_algG_acute_1$Cs_acute_1/data_algG_acute_1$coverage_acute_1




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_acute_1$highest_percentage <- 0
data_algG_acute_1$highest_percentage <- apply(data_algG_acute_1[, 8:11], 1, max)
data_algG_acute_1$second_percentage <- apply(data_algG_acute_1[, 8:11], 1, function(i) sort(i)[ dim(data_algG_acute_1[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_acute_1_filtered <- subset(data_algG_acute_1, second_percentage >= 0.025) 
data_algG_acute_1_filtered$main_base_acute_1 <- colnames(data_algG_acute_1_filtered[, 8:11])[apply(data_algG_acute_1_filtered[, 8:11],1,which.max)]

data_algG_acute_1_filtered_4percent <- subset(data_algG_acute_1, second_percentage >= 0.04) 
data_algG_acute_1_filtered_4percent$main_base_acute_1 <- colnames(data_algG_acute_1_filtered_4percent[, 8:11])[apply(data_algG_acute_1_filtered_4percent[, 8:11],1,which.max)]

########### acute_2   ##############

# Calculate percentage for each base

data_algG_acute_2$percent_Ts_acute_2 <- data_algG_acute_2$Ts_acute_2/data_algG_acute_2$coverage_acute_2
data_algG_acute_2$percent_As_acute_2 <- data_algG_acute_2$As_acute_2/data_algG_acute_2$coverage_acute_2
data_algG_acute_2$percent_Gs_acute_2 <- data_algG_acute_2$Gs_acute_2/data_algG_acute_2$coverage_acute_2
data_algG_acute_2$percent_Cs_acute_2 <- data_algG_acute_2$Cs_acute_2/data_algG_acute_2$coverage_acute_2




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_acute_2$highest_percentage <- 0
data_algG_acute_2$highest_percentage <- apply(data_algG_acute_2[, 8:11], 1, max)
data_algG_acute_2$second_percentage <- apply(data_algG_acute_2[, 8:11], 1, function(i) sort(i)[ dim(data_algG_acute_2[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_acute_2_filtered <- subset(data_algG_acute_2, second_percentage >= 0.025) 
data_algG_acute_2_filtered$main_base_acute_2 <- colnames(data_algG_acute_2_filtered[, 8:11])[apply(data_algG_acute_2_filtered[, 8:11],1,which.max)]

data_algG_acute_2_filtered_4percent <- subset(data_algG_acute_2, second_percentage >= 0.04) 
data_algG_acute_2_filtered_4percent$main_base_acute_2 <- colnames(data_algG_acute_2_filtered_4percent[, 8:11])[apply(data_algG_acute_2_filtered_4percent[, 8:11],1,which.max)]

########### acute_3   ##############

# Calculate percentage for each base

data_algG_acute_3$percent_Ts_acute_3 <- data_algG_acute_3$Ts_acute_3/data_algG_acute_3$coverage_acute_3
data_algG_acute_3$percent_As_acute_3 <- data_algG_acute_3$As_acute_3/data_algG_acute_3$coverage_acute_3
data_algG_acute_3$percent_Gs_acute_3 <- data_algG_acute_3$Gs_acute_3/data_algG_acute_3$coverage_acute_3
data_algG_acute_3$percent_Cs_acute_3 <- data_algG_acute_3$Cs_acute_3/data_algG_acute_3$coverage_acute_3




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_acute_3$highest_percentage <- 0
data_algG_acute_3$highest_percentage <- apply(data_algG_acute_3[, 8:11], 1, max)
data_algG_acute_3$second_percentage <- apply(data_algG_acute_3[, 8:11], 1, function(i) sort(i)[ dim(data_algG_acute_3[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_acute_3_filtered <- subset(data_algG_acute_3, second_percentage >= 0.025) 
data_algG_acute_3_filtered$main_base_acute_3 <- colnames(data_algG_acute_3_filtered[, 8:11])[apply(data_algG_acute_3_filtered[, 8:11],1,which.max)]

data_algG_acute_3_filtered_4percent <- subset(data_algG_acute_3, second_percentage >= 0.04) 
data_algG_acute_3_filtered_4percent$main_base_acute_3 <- colnames(data_algG_acute_3_filtered_4percent[, 8:11])[apply(data_algG_acute_3_filtered_4percent[, 8:11],1,which.max)]

########### env_1   ##############

# Calculate percentage for each base

data_algG_env_1$percent_Ts_env_1 <- data_algG_env_1$Ts_env_1/data_algG_env_1$coverage_env_1
data_algG_env_1$percent_As_env_1 <- data_algG_env_1$As_env_1/data_algG_env_1$coverage_env_1
data_algG_env_1$percent_Gs_env_1 <- data_algG_env_1$Gs_env_1/data_algG_env_1$coverage_env_1
data_algG_env_1$percent_Cs_env_1 <- data_algG_env_1$Cs_env_1/data_algG_env_1$coverage_env_1




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_env_1$highest_percentage <- 0
data_algG_env_1$highest_percentage <- apply(data_algG_env_1[, 8:11], 1, max)
data_algG_env_1$second_percentage <- apply(data_algG_env_1[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_1[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_env_1_filtered <- subset(data_algG_env_1, second_percentage >= 0.025) 
data_algG_env_1_filtered$main_base_env_1 <- colnames(data_algG_env_1_filtered[, 8:11])[apply(data_algG_env_1_filtered[, 8:11],1,which.max)]

data_algG_env_1_filtered_4percent <- subset(data_algG_env_1, second_percentage >= 0.04) 
data_algG_env_1_filtered_4percent$main_base_env_1 <- colnames(data_algG_env_1_filtered_4percent[, 8:11])[apply(data_algG_env_1_filtered_4percent[, 8:11],1,which.max)]

########### env_2   ##############

# Calculate percentage for each base

data_algG_env_2$percent_Ts_env_2 <- data_algG_env_2$Ts_env_2/data_algG_env_2$coverage_env_2
data_algG_env_2$percent_As_env_2 <- data_algG_env_2$As_env_2/data_algG_env_2$coverage_env_2
data_algG_env_2$percent_Gs_env_2 <- data_algG_env_2$Gs_env_2/data_algG_env_2$coverage_env_2
data_algG_env_2$percent_Cs_env_2 <- data_algG_env_2$Cs_env_2/data_algG_env_2$coverage_env_2




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_env_2$highest_percentage <- 0
data_algG_env_2$highest_percentage <- apply(data_algG_env_2[, 8:11], 1, max)
data_algG_env_2$second_percentage <- apply(data_algG_env_2[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_2[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_env_2_filtered <- subset(data_algG_env_2, second_percentage >= 0.025) 
data_algG_env_2_filtered$main_base_env_2 <- colnames(data_algG_env_2_filtered[, 8:11])[apply(data_algG_env_2_filtered[, 8:11],1,which.max)]

data_algG_env_2_filtered_4percent <- subset(data_algG_env_2, second_percentage >= 0.04) 
data_algG_env_2_filtered_4percent$main_base_env_2 <- colnames(data_algG_env_2_filtered_4percent[, 8:11])[apply(data_algG_env_2_filtered_4percent[, 8:11],1,which.max)]

########### env_3   ##############

# Calculate percentage for each base

data_algG_env_3$percent_Ts_env_3 <- data_algG_env_3$Ts_env_3/data_algG_env_3$coverage_env_3
data_algG_env_3$percent_As_env_3 <- data_algG_env_3$As_env_3/data_algG_env_3$coverage_env_3
data_algG_env_3$percent_Gs_env_3 <- data_algG_env_3$Gs_env_3/data_algG_env_3$coverage_env_3
data_algG_env_3$percent_Cs_env_3 <- data_algG_env_3$Cs_env_3/data_algG_env_3$coverage_env_3




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_env_3$highest_percentage <- 0
data_algG_env_3$highest_percentage <- apply(data_algG_env_3[, 8:11], 1, max)
data_algG_env_3$second_percentage <- apply(data_algG_env_3[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_3[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_env_3_filtered <- subset(data_algG_env_3, second_percentage >= 0.025) 
data_algG_env_3_filtered$main_base_env_3 <- colnames(data_algG_env_3_filtered[, 8:11])[apply(data_algG_env_3_filtered[, 8:11],1,which.max)]

data_algG_env_3_filtered_4percent <- subset(data_algG_env_3, second_percentage >= 0.04) 
data_algG_env_3_filtered_4percent$main_base_env_3 <- colnames(data_algG_env_3_filtered_4percent[, 8:11])[apply(data_algG_env_3_filtered_4percent[, 8:11],1,which.max)]

########### env_4   ##############

# Calculate percentage for each base

data_algG_env_4$percent_Ts_env_4 <- data_algG_env_4$Ts_env_4/data_algG_env_4$coverage_env_4
data_algG_env_4$percent_As_env_4 <- data_algG_env_4$As_env_4/data_algG_env_4$coverage_env_4
data_algG_env_4$percent_Gs_env_4 <- data_algG_env_4$Gs_env_4/data_algG_env_4$coverage_env_4
data_algG_env_4$percent_Cs_env_4 <- data_algG_env_4$Cs_env_4/data_algG_env_4$coverage_env_4




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_env_4$highest_percentage <- 0
data_algG_env_4$highest_percentage <- apply(data_algG_env_4[, 8:11], 1, max)
data_algG_env_4$second_percentage <- apply(data_algG_env_4[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_4[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_env_4_filtered <- subset(data_algG_env_4, second_percentage >= 0.025) 
data_algG_env_4_filtered$main_base_env_4 <- colnames(data_algG_env_4_filtered[, 8:11])[apply(data_algG_env_4_filtered[, 8:11],1,which.max)]

data_algG_env_4_filtered_4percent <- subset(data_algG_env_4, second_percentage >= 0.04) 
data_algG_env_4_filtered_4percent$main_base_env_4 <- colnames(data_algG_env_4_filtered_4percent[, 8:11])[apply(data_algG_env_4_filtered_4percent[, 8:11],1,which.max)]


########### env_5   ##############

# Calculate percentage for each base

data_algG_env_5$percent_Ts_env_5 <- data_algG_env_5$Ts_env_5/data_algG_env_5$coverage_env_5
data_algG_env_5$percent_As_env_5 <- data_algG_env_5$As_env_5/data_algG_env_5$coverage_env_5
data_algG_env_5$percent_Gs_env_5 <- data_algG_env_5$Gs_env_5/data_algG_env_5$coverage_env_5
data_algG_env_5$percent_Cs_env_5 <- data_algG_env_5$Cs_env_5/data_algG_env_5$coverage_env_5




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_env_5$highest_percentage <- 0
data_algG_env_5$highest_percentage <- apply(data_algG_env_5[, 8:11], 1, max)
data_algG_env_5$second_percentage <- apply(data_algG_env_5[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_5[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_env_5_filtered <- subset(data_algG_env_5, second_percentage >= 0.025) 
data_algG_env_5_filtered$main_base_env_5 <- colnames(data_algG_env_5_filtered[, 8:11])[apply(data_algG_env_5_filtered[, 8:11],1,which.max)]

data_algG_env_5_filtered_4percent <- subset(data_algG_env_5, second_percentage >= 0.04) 
data_algG_env_5_filtered_4percent$main_base_env_5 <- colnames(data_algG_env_5_filtered_4percent[, 8:11])[apply(data_algG_env_5_filtered_4percent[, 8:11],1,which.max)]

########### COPD   ##############

# Calculate percentage for each base

data_algG_COPD$percent_Ts_COPD <- data_algG_COPD$Ts_COPD/data_algG_COPD$coverage_COPD
data_algG_COPD$percent_As_COPD <- data_algG_COPD$As_COPD/data_algG_COPD$coverage_COPD
data_algG_COPD$percent_Gs_COPD <- data_algG_COPD$Gs_COPD/data_algG_COPD$coverage_COPD
data_algG_COPD$percent_Cs_COPD <- data_algG_COPD$Cs_COPD/data_algG_COPD$coverage_COPD




# Get highest and second highest. This is important because at least 2.5 % of the Reads need to be a SNP.

data_algG_COPD$highest_percentage <- 0
data_algG_COPD$highest_percentage <- apply(data_algG_COPD[, 8:11], 1, max)
data_algG_COPD$second_percentage <- apply(data_algG_COPD[, 8:11], 1, function(i) sort(i)[ dim(data_algG_COPD[, 8:11])[2]-1])


# Subset only those with the second highest minimum 2.5 %

data_algG_COPD_filtered <- subset(data_algG_COPD, second_percentage >= 0.025) 
data_algG_COPD_filtered$main_base_COPD <- colnames(data_algG_COPD_filtered[, 8:11])[apply(data_algG_COPD_filtered[, 8:11],1,which.max)]

data_algG_COPD_filtered_4percent <- subset(data_algG_COPD, second_percentage >= 0.04) 
data_algG_COPD_filtered_4percent$main_base_COPD <- colnames(data_algG_COPD_filtered_4percent[, 8:11])[apply(data_algG_COPD_filtered_4percent[, 8:11],1,which.max)]


