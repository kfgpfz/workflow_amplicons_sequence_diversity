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
library(Rfast)
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
strand <- "rev"
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
primer_correction_algG_CF_1_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T1/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_1_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T1/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_1 <- data_algG_CF_1[!(data_algG_CF_1$Position %in% primer_correction_algG_CF_1_fwd$Position),]
data_algG_CF_1 <- data_algG_CF_1[!(data_algG_CF_1$Position %in% primer_correction_algG_CF_1_rev$Position),]
data_algG_CF_1 <- rbind(data_algG_CF_1,primer_correction_algG_CF_1_fwd,primer_correction_algG_CF_1_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_1 <- subset(data_algG_CF_1, coverage_CF_1 >199)
data_algG_CF_1_5UTR <- subset(data_algG_CF_1_5UTR, coverage_CF_1_5UTR >199)
data_algG_CF_1_3UTR <- subset(data_algG_CF_1_3UTR, coverage_CF_1_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_1_only_SNP <- subset(data_algG_CF_1, second_percentage <= 0.025) 
data_algG_CF_1_only_SNP$main_base_CF_1 <- colnames(data_algG_CF_1_only_SNP[, 8:11])[apply(data_algG_CF_1_only_SNP[, 8:11],1,which.max)]
data_algG_CF_1_only_SNP$main_base_CF_1 <- ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_As_CF_1","A", ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_Cs_CF_1","C", ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_Gs_CF_1","G", ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_Ts_CF_1","T",z<-0))))
data_algG_CF_1_only_SNP <- subset(data_algG_CF_1_only_SNP,GenBase_PA14 != main_base_CF_1)
ifelse(nrow(data_algG_CF_1_only_SNP) > 0,data_algG_CF_1_only_SNP$Isolates_CF_1 <- isolates_CF_1, z <- 0)


ifelse(nrow(data_algG_CF_1_only_SNP) > 0,data_algG_CF_1_only_SNP_clean <- data_algG_CF_1_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_1_only_SNP) > 0,names(data_algG_CF_1_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_1","Isolates_CF_1"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_1_filtered <- subset(data_algG_CF_1, second_percentage >= 0.025) 
data_algG_CF_1_filtered$main_base_CF_1 <- colnames(data_algG_CF_1_filtered[, 8:11])[apply(data_algG_CF_1_filtered[, 8:11],1,which.max)]
data_algG_CF_1_filtered$second_base_CF_1 <- colnames(data_algG_CF_1_filtered[, 8:11])[apply(data_algG_CF_1_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_1_filtered$third_base_CF_1 <- ifelse(apply(data_algG_CF_1_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_1_filtered$`PA14-Koordinate`)
  if(data_algG_CF_1_filtered$main_base_CF_1[a] == "percent_As_CF_1"){
    data_algG_CF_1_filtered$main_base_CF_1[a] <- "A"                     
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] == "percent_Cs_CF_1"){
    data_algG_CF_1_filtered$main_base_CF_1[a] <- "C"                     
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] == "percent_Gs_CF_1"){
    data_algG_CF_1_filtered$main_base_CF_1[a] <- "G"                     
  } else{
    data_algG_CF_1_filtered$main_base_CF_1[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_1_filtered$`PA14-Koordinate`)
  if(data_algG_CF_1_filtered$second_base_CF_1[a] == "percent_As_CF_1"){
    data_algG_CF_1_filtered$second_base_CF_1[a] <- "A"                     
  } else if(data_algG_CF_1_filtered$second_base_CF_1[a] == "percent_Cs_CF_1"){
    data_algG_CF_1_filtered$second_base_CF_1[a] <- "C"                     
  } else if(data_algG_CF_1_filtered$second_base_CF_1[a] == "percent_Gs_CF_1"){
    data_algG_CF_1_filtered$second_base_CF_1[a] <- "G"                     
  } else{
    data_algG_CF_1_filtered$second_base_CF_1[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_1_filtered$`PA14-Koordinate`)
  if(data_algG_CF_1_filtered$main_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$third_base_CF_1[a] == 0){
    data_algG_CF_1_filtered$third_base_CF_1[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_1_filtered$`PA14-Koordinate`)
  if(data_algG_CF_1_filtered$third_base_CF_1[a] == 1){
    data_algG_CF_1_filtered_third_variant <- data_algG_CF_1_filtered[a,]
    data_algG_CF_1_filtered$third_base_CF_1[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_1_filtered$`PA14-Koordinate`)
  if(data_algG_CF_1_filtered$third_base_CF_1[a] == 2){
    data_algG_CF_1_filtered_third_variant_no_ref <- data_algG_CF_1_filtered[a,]
    data_algG_CF_1_filtered$third_base_CF_1[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_1_filtered$third_base_CF_1) == 0, data_algG_CF_1_filtered <- subset(data_algG_CF_1_filtered, select = -c(third_base_CF_1)),data_algG_CF_1_filtered <- data_algG_CF_1_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_filtered_third_variant$third_percentage <- apply(data_algG_CF_1_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_1_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_filtered_third_variant$third_base_CF_1 <- colnames(data_algG_CF_1_filtered_third_variant[, 8:11])[apply(data_algG_CF_1_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_As_CF_1","A", ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_Cs_CF_1","C", ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_Gs_CF_1","G", ifelse(data_algG_CF_1_only_SNP$main_base_CF_1 == "percent_Ts_CF_1","T",z<-0))))

ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,
       data_algG_CF_1_filtered_third_variant$third_base_CF_1 <- 
         ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "percent_As_CF_1","A", 
                ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "percent_Cs_CF_1","C",
                       ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "percent_Gs_CF_1","G",
                              ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "percent_Ts_CF_1","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,
       data_algG_CF_1_filtered_third_variant$third_base_CF_1 <- ifelse(
         data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14,data_algG_CF_1_filtered_third_variant$third_base_CF_1 <- data_algG_CF_1_filtered_third_variant$second_base_CF_1,data_algG_CF_1_filtered_third_variant$third_base_CF_1
       ),
       z<-0
)

ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,
       data_algG_CF_1_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14, data_algG_CF_1_filtered_third_variant$second_percentage,data_algG_CF_1_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_filtered_third_variant <- data_algG_CF_1_filtered_third_variant[,-which(names(data_algG_CF_1_filtered_third_variant) %in% c("second_percentage","second_base_CF_1"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_1_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_1_filtered_third_variant_no_ref)[names(data_algG_CF_1_filtered_third_variant_no_ref) == "second_base_CF_1"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_1_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_1_filtered_third_variant_no_ref)[names(data_algG_CF_1_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_1",z<-0)




# Define SNP percentage & Add it

data_algG_CF_1_filtered$SNP_base <- 0
data_algG_CF_1_filtered$SNP_percentage_CF_1 <- 0


for(i in data_algG_CF_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_1_filtered$`PA14-Koordinate`)
  if(data_algG_CF_1_filtered$main_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] == "A"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_As_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] == "C"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Cs_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] == "G"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Gs_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] == "T"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Ts_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_1_filtered$second_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "A"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_As_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_1_filtered$second_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "C"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Cs_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_1_filtered$second_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "G"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Gs_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_1_filtered$second_base_CF_1[a] == data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "T"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Ts_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "A"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_As_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "C"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Cs_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "G"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Gs_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_1_filtered$main_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$second_base_CF_1[a] != data_algG_CF_1_filtered$GenBase_PA14[a] & data_algG_CF_1_filtered$main_base_CF_1[a] == "T"){
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- data_algG_CF_1_filtered$percent_Ts_CF_1[a]
    data_algG_CF_1_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_1_filtered$SNP_percentage_CF_1[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_filtered_third_variant$SNP_percentage_CF_1 <- 0,z<-0)


ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,
       data_algG_CF_1_filtered_third_variant$SNP_percentage_CF_1 <- 
         ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "A",
                data_algG_CF_1_filtered_third_variant$percent_As_CF_1, 
                ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "C",
                       data_algG_CF_1_filtered_third_variant$percent_Cs_CF_1,
                       ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "G",
                              data_algG_CF_1_filtered_third_variant$percent_Gs_CF_1,
                              ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "T",
                                     data_algG_CF_1_filtered_third_variant$percent_Ts_CF_1,
                                     ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "A",
                                            data_algG_CF_1_filtered_third_variant$percent_As_CF_1,
                                            ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "C",
                                                   data_algG_CF_1_filtered_third_variant$percent_Cs_CF_1,
                                                   ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "G",
                                                          data_algG_CF_1_filtered_third_variant$percent_Gs_CF_1,
                                                          ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "T",
                                                                 data_algG_CF_1_filtered_third_variant$percent_Ts_CF_1,
                                                                 ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "A",
                                                                        data_algG_CF_1_filtered_third_variant$percent_As_CF_1,
                                                                        ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "C",
                                                                               data_algG_CF_1_filtered_third_variant$percent_Cs_CF_1,
                                                                               ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "G",
                                                                                      data_algG_CF_1_filtered_third_variant$percent_Gs_CF_1,
                                                                                      ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "T",
                                                                                             data_algG_CF_1_filtered_third_variant$percent_Ts_CF_1,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,
       data_algG_CF_1_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "A",
                "A", 
                ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "C",
                       "C",
                       ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "G",
                              "G",
                              ifelse(data_algG_CF_1_filtered_third_variant$main_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "T",
                                     "T",
                                     ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "A",
                                            "A",
                                            ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 == data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$main_base_CF_1 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_1_filtered_third_variant$third_base_CF_1 != data_algG_CF_1_filtered_third_variant$GenBase_PA14 & data_algG_CF_1_filtered_third_variant$third_base_CF_1 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_1_filtered$Isolates_CF_1 <- round((data_algG_CF_1_filtered$SNP_percentage_CF_1*100)/(100/isolates_CF_1))
ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_filtered_third_variant$Isolates_CF_1 <- round((data_algG_CF_1_filtered_third_variant$SNP_percentage_CF_1*100)/(100/isolates_CF_1)),z<-0)
ifelse(exists("data_algG_CF_1_filtered_third_variant_no_ref") == TRUE,data_algG_CF_1_filtered_third_variant_no_ref$Isolates_CF_1 <- round((data_algG_CF_1_filtered_third_variant_no_ref$SNP_percentage_CF_1*100)/(100/isolates_CF_1)),z<-0)



# Clean data

data_algG_CF_1_clean <- data_algG_CF_1_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_clean_third_base <- data_algG_CF_1_filtered_third_variant[,which(names(data_algG_CF_1_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_1","Isolates_CF_1"))],z<-0)
ifelse(exists("data_algG_CF_1_filtered_third_variant") == TRUE,data_algG_CF_1_clean <- rbind(data_algG_CF_1_clean,data_algG_CF_1_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_1_filtered_third_variant_no_ref") == TRUE,data_algG_CF_1_clean_third_base_no_ref <- data_algG_CF_1_filtered_third_variant_no_ref[,which(names(data_algG_CF_1_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_1","Isolates_CF_1"))],z<-0)
ifelse(exists("data_algG_CF_1_filtered_third_variant_no_ref") == TRUE,data_algG_CF_1_clean <- rbind(data_algG_CF_1_clean,data_algG_CF_1_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_1_clean_20_filter <- merge(data_algG_CF_1_clean, data_algG_CF_1_filtered[,c(1,3:6)],by="PA14-Koordinate")
#data_algG_CF_1_clean_20_filter$helper <- paste(data_algG_CF_1_clean_20_filter$`PA14-Koordinate`,data_algG_CF_1_clean_20_filter$SNP_base,sep="")
#data_algG_CF_1_filtered$helper <- paste(data_algG_CF_1_filtered$`PA14-Koordinate`,data_algG_CF_1_filtered$SNP_base,sep="")
data_algG_CF_1_clean_cleared <- data_algG_CF_1_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_1_clean_20_filter)){
    if(data_algG_CF_1_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_1_clean_20_filter$As_CF_1[i] >= 20){
    data_algG_CF_1_clean_cleared <- rbind(data_algG_CF_1_clean_cleared,data_algG_CF_1_clean_20_filter[i,])
  } else if(data_algG_CF_1_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_1_clean_20_filter$Cs_CF_1[i] >= 20){
    data_algG_CF_1_clean_cleared <- rbind(data_algG_CF_1_clean_cleared,data_algG_CF_1_clean_20_filter[i,])
  } else if(data_algG_CF_1_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_1_clean_20_filter$Gs_CF_1[i] >= 20){
    data_algG_CF_1_clean_cleared <- rbind(data_algG_CF_1_clean_cleared,data_algG_CF_1_clean_20_filter[i,])
  }else if(data_algG_CF_1_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_1_clean_20_filter$Ts_CF_1[i] >= 20){
    data_algG_CF_1_clean_cleared <- rbind(data_algG_CF_1_clean_cleared,data_algG_CF_1_clean_20_filter[i,])
  }
}


data_algG_CF_1_clean <- data_algG_CF_1_clean_cleared[,c(1:5)]

### CF_2

# Input 
data_algG_CF_2 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T2/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_2_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T2/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_2_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T2/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_2 <- data_algG_CF_2[!(data_algG_CF_2$Position %in% primer_correction_algG_CF_2_fwd$Position),]
data_algG_CF_2 <- data_algG_CF_2[!(data_algG_CF_2$Position %in% primer_correction_algG_CF_2_rev$Position),]
data_algG_CF_2 <- rbind(data_algG_CF_2,primer_correction_algG_CF_2_fwd,primer_correction_algG_CF_2_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_2 <- subset(data_algG_CF_2, coverage_CF_2 >199)
data_algG_CF_2_5UTR <- subset(data_algG_CF_2_5UTR, coverage_CF_2_5UTR >199)
data_algG_CF_2_3UTR <- subset(data_algG_CF_2_3UTR, coverage_CF_2_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_2_only_SNP <- subset(data_algG_CF_2, second_percentage <= 0.025) 
data_algG_CF_2_only_SNP$main_base_CF_2 <- colnames(data_algG_CF_2_only_SNP[, 8:11])[apply(data_algG_CF_2_only_SNP[, 8:11],1,which.max)]
data_algG_CF_2_only_SNP$main_base_CF_2 <- ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_As_CF_2","A", ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_Cs_CF_2","C", ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_Gs_CF_2","G", ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_Ts_CF_2","T",z<-0))))
data_algG_CF_2_only_SNP <- subset(data_algG_CF_2_only_SNP,GenBase_PA14 != main_base_CF_2)
ifelse(nrow(data_algG_CF_2_only_SNP) > 0,data_algG_CF_2_only_SNP$Isolates_CF_2 <- isolates_CF_2, z <- 0)


ifelse(nrow(data_algG_CF_2_only_SNP) > 0,data_algG_CF_2_only_SNP_clean <- data_algG_CF_2_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_2_only_SNP) > 0,names(data_algG_CF_2_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_2","Isolates_CF_2"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_2_filtered <- subset(data_algG_CF_2, second_percentage >= 0.025) 
data_algG_CF_2_filtered$main_base_CF_2 <- colnames(data_algG_CF_2_filtered[, 8:11])[apply(data_algG_CF_2_filtered[, 8:11],1,which.max)]
data_algG_CF_2_filtered$second_base_CF_2 <- colnames(data_algG_CF_2_filtered[, 8:11])[apply(data_algG_CF_2_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_2_filtered$third_base_CF_2 <- ifelse(apply(data_algG_CF_2_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_2_filtered$`PA14-Koordinate`)
  if(data_algG_CF_2_filtered$main_base_CF_2[a] == "percent_As_CF_2"){
    data_algG_CF_2_filtered$main_base_CF_2[a] <- "A"                     
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] == "percent_Cs_CF_2"){
    data_algG_CF_2_filtered$main_base_CF_2[a] <- "C"                     
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] == "percent_Gs_CF_2"){
    data_algG_CF_2_filtered$main_base_CF_2[a] <- "G"                     
  } else{
    data_algG_CF_2_filtered$main_base_CF_2[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_2_filtered$`PA14-Koordinate`)
  if(data_algG_CF_2_filtered$second_base_CF_2[a] == "percent_As_CF_2"){
    data_algG_CF_2_filtered$second_base_CF_2[a] <- "A"                     
  } else if(data_algG_CF_2_filtered$second_base_CF_2[a] == "percent_Cs_CF_2"){
    data_algG_CF_2_filtered$second_base_CF_2[a] <- "C"                     
  } else if(data_algG_CF_2_filtered$second_base_CF_2[a] == "percent_Gs_CF_2"){
    data_algG_CF_2_filtered$second_base_CF_2[a] <- "G"                     
  } else{
    data_algG_CF_2_filtered$second_base_CF_2[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_2_filtered$`PA14-Koordinate`)
  if(data_algG_CF_2_filtered$main_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$third_base_CF_2[a] == 0){
    data_algG_CF_2_filtered$third_base_CF_2[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_2_filtered$`PA14-Koordinate`)
  if(data_algG_CF_2_filtered$third_base_CF_2[a] == 1){
    data_algG_CF_2_filtered_third_variant <- data_algG_CF_2_filtered[a,]
    data_algG_CF_2_filtered$third_base_CF_2[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_2_filtered$`PA14-Koordinate`)
  if(data_algG_CF_2_filtered$third_base_CF_2[a] == 2){
    data_algG_CF_2_filtered_third_variant_no_ref <- data_algG_CF_2_filtered[a,]
    data_algG_CF_2_filtered$third_base_CF_2[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_2_filtered$third_base_CF_2) == 0, data_algG_CF_2_filtered <- subset(data_algG_CF_2_filtered, select = -c(third_base_CF_2)),data_algG_CF_2_filtered <- data_algG_CF_2_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_filtered_third_variant$third_percentage <- apply(data_algG_CF_2_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_2_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_filtered_third_variant$third_base_CF_2 <- colnames(data_algG_CF_2_filtered_third_variant[, 8:11])[apply(data_algG_CF_2_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_As_CF_2","A", ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_Cs_CF_2","C", ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_Gs_CF_2","G", ifelse(data_algG_CF_2_only_SNP$main_base_CF_2 == "percent_Ts_CF_2","T",z<-0))))

ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,
       data_algG_CF_2_filtered_third_variant$third_base_CF_2 <- 
         ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "percent_As_CF_2","A", 
                ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "percent_Cs_CF_2","C",
                       ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "percent_Gs_CF_2","G",
                              ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "percent_Ts_CF_2","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,
       data_algG_CF_2_filtered_third_variant$third_base_CF_2 <- ifelse(
         data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14,data_algG_CF_2_filtered_third_variant$third_base_CF_2 <- data_algG_CF_2_filtered_third_variant$second_base_CF_2,data_algG_CF_2_filtered_third_variant$third_base_CF_2
       ),
       z<-0
)

ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,
       data_algG_CF_2_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14, data_algG_CF_2_filtered_third_variant$second_percentage,data_algG_CF_2_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_filtered_third_variant <- data_algG_CF_2_filtered_third_variant[,-which(names(data_algG_CF_2_filtered_third_variant) %in% c("second_percentage","second_base_CF_2"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_2_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_2_filtered_third_variant_no_ref)[names(data_algG_CF_2_filtered_third_variant_no_ref) == "second_base_CF_2"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_2_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_2_filtered_third_variant_no_ref)[names(data_algG_CF_2_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_2",z<-0)




# Define SNP percentage & Add it

data_algG_CF_2_filtered$SNP_base <- 0
data_algG_CF_2_filtered$SNP_percentage_CF_2 <- 0


for(i in data_algG_CF_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_2_filtered$`PA14-Koordinate`)
  if(data_algG_CF_2_filtered$main_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] == "A"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_As_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] == "C"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Cs_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] == "G"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Gs_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] == "T"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Ts_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_2_filtered$second_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "A"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_As_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_2_filtered$second_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "C"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Cs_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_2_filtered$second_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "G"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Gs_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_2_filtered$second_base_CF_2[a] == data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "T"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Ts_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "A"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_As_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "C"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Cs_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "G"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Gs_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_2_filtered$main_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$second_base_CF_2[a] != data_algG_CF_2_filtered$GenBase_PA14[a] & data_algG_CF_2_filtered$main_base_CF_2[a] == "T"){
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- data_algG_CF_2_filtered$percent_Ts_CF_2[a]
    data_algG_CF_2_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_2_filtered$SNP_percentage_CF_2[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_filtered_third_variant$SNP_percentage_CF_2 <- 0,z<-0)


ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,
       data_algG_CF_2_filtered_third_variant$SNP_percentage_CF_2 <- 
         ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "A",
                data_algG_CF_2_filtered_third_variant$percent_As_CF_2, 
                ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "C",
                       data_algG_CF_2_filtered_third_variant$percent_Cs_CF_2,
                       ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "G",
                              data_algG_CF_2_filtered_third_variant$percent_Gs_CF_2,
                              ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "T",
                                     data_algG_CF_2_filtered_third_variant$percent_Ts_CF_2,
                                     ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "A",
                                            data_algG_CF_2_filtered_third_variant$percent_As_CF_2,
                                            ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "C",
                                                   data_algG_CF_2_filtered_third_variant$percent_Cs_CF_2,
                                                   ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "G",
                                                          data_algG_CF_2_filtered_third_variant$percent_Gs_CF_2,
                                                          ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "T",
                                                                 data_algG_CF_2_filtered_third_variant$percent_Ts_CF_2,
                                                                 ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "A",
                                                                        data_algG_CF_2_filtered_third_variant$percent_As_CF_2,
                                                                        ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "C",
                                                                               data_algG_CF_2_filtered_third_variant$percent_Cs_CF_2,
                                                                               ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "G",
                                                                                      data_algG_CF_2_filtered_third_variant$percent_Gs_CF_2,
                                                                                      ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "T",
                                                                                             data_algG_CF_2_filtered_third_variant$percent_Ts_CF_2,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,
       data_algG_CF_2_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "A",
                "A", 
                ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "C",
                       "C",
                       ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "G",
                              "G",
                              ifelse(data_algG_CF_2_filtered_third_variant$main_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "T",
                                     "T",
                                     ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "A",
                                            "A",
                                            ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 == data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$main_base_CF_2 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_2_filtered_third_variant$third_base_CF_2 != data_algG_CF_2_filtered_third_variant$GenBase_PA14 & data_algG_CF_2_filtered_third_variant$third_base_CF_2 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_2_filtered$Isolates_CF_2 <- round((data_algG_CF_2_filtered$SNP_percentage_CF_2*100)/(100/isolates_CF_2))
ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_filtered_third_variant$Isolates_CF_2 <- round((data_algG_CF_2_filtered_third_variant$SNP_percentage_CF_2*100)/(100/isolates_CF_2)),z<-0)
ifelse(exists("data_algG_CF_2_filtered_third_variant_no_ref") == TRUE,data_algG_CF_2_filtered_third_variant_no_ref$Isolates_CF_2 <- round((data_algG_CF_2_filtered_third_variant_no_ref$SNP_percentage_CF_2*100)/(100/isolates_CF_2)),z<-0)



# Clean data

data_algG_CF_2_clean <- data_algG_CF_2_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_clean_third_base <- data_algG_CF_2_filtered_third_variant[,which(names(data_algG_CF_2_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_2","Isolates_CF_2"))],z<-0)
ifelse(exists("data_algG_CF_2_filtered_third_variant") == TRUE,data_algG_CF_2_clean <- rbind(data_algG_CF_2_clean,data_algG_CF_2_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_2_filtered_third_variant_no_ref") == TRUE,data_algG_CF_2_clean_third_base_no_ref <- data_algG_CF_2_filtered_third_variant_no_ref[,which(names(data_algG_CF_2_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_2","Isolates_CF_2"))],z<-0)
ifelse(exists("data_algG_CF_2_filtered_third_variant_no_ref") == TRUE,data_algG_CF_2_clean <- rbind(data_algG_CF_2_clean,data_algG_CF_2_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_2_clean_20_filter <- merge(data_algG_CF_2_clean, data_algG_CF_2_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_2_clean_cleared <- data_algG_CF_2_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_2_clean_20_filter)){
  if(data_algG_CF_2_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_2_clean_20_filter$As_CF_2[i] >= 20){
    data_algG_CF_2_clean_cleared <- rbind(data_algG_CF_2_clean_cleared,data_algG_CF_2_clean_20_filter[i,])
  } else if(data_algG_CF_2_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_2_clean_20_filter$Cs_CF_2[i] >= 20){
    data_algG_CF_2_clean_cleared <- rbind(data_algG_CF_2_clean_cleared,data_algG_CF_2_clean_20_filter[i,])
  } else if(data_algG_CF_2_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_2_clean_20_filter$Gs_CF_2[i] >= 20){
    data_algG_CF_2_clean_cleared <- rbind(data_algG_CF_2_clean_cleared,data_algG_CF_2_clean_20_filter[i,])
  }else if(data_algG_CF_2_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_2_clean_20_filter$Ts_CF_2[i] >= 20){
    data_algG_CF_2_clean_cleared <- rbind(data_algG_CF_2_clean_cleared,data_algG_CF_2_clean_20_filter[i,])
  }
}

data_algG_CF_2_clean <- data_algG_CF_2_clean_cleared[,c(1:5)]

### CF_3

# Input 
data_algG_CF_3 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T3/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_3_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T3/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_3_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T3/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_3 <- data_algG_CF_3[!(data_algG_CF_3$Position %in% primer_correction_algG_CF_3_fwd$Position),]
data_algG_CF_3 <- data_algG_CF_3[!(data_algG_CF_3$Position %in% primer_correction_algG_CF_3_rev$Position),]
data_algG_CF_3 <- rbind(data_algG_CF_3,primer_correction_algG_CF_3_fwd,primer_correction_algG_CF_3_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_3 <- subset(data_algG_CF_3, coverage_CF_3 >199)
data_algG_CF_3_5UTR <- subset(data_algG_CF_3_5UTR, coverage_CF_3_5UTR >199)
data_algG_CF_3_3UTR <- subset(data_algG_CF_3_3UTR, coverage_CF_3_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_3_only_SNP <- subset(data_algG_CF_3, second_percentage <= 0.025) 
data_algG_CF_3_only_SNP$main_base_CF_3 <- colnames(data_algG_CF_3_only_SNP[, 8:11])[apply(data_algG_CF_3_only_SNP[, 8:11],1,which.max)]
data_algG_CF_3_only_SNP$main_base_CF_3 <- ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_As_CF_3","A", ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_Cs_CF_3","C", ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_Gs_CF_3","G", ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_Ts_CF_3","T",z<-0))))
data_algG_CF_3_only_SNP <- subset(data_algG_CF_3_only_SNP,GenBase_PA14 != main_base_CF_3)
ifelse(nrow(data_algG_CF_3_only_SNP) > 0,data_algG_CF_3_only_SNP$Isolates_CF_3 <- isolates_CF_3, z <- 0)


ifelse(nrow(data_algG_CF_3_only_SNP) > 0,data_algG_CF_3_only_SNP_clean <- data_algG_CF_3_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_3_only_SNP) > 0,names(data_algG_CF_3_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_3","Isolates_CF_3"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_3_filtered <- subset(data_algG_CF_3, second_percentage >= 0.025) 
data_algG_CF_3_filtered$main_base_CF_3 <- colnames(data_algG_CF_3_filtered[, 8:11])[apply(data_algG_CF_3_filtered[, 8:11],1,which.max)]
data_algG_CF_3_filtered$second_base_CF_3 <- colnames(data_algG_CF_3_filtered[, 8:11])[apply(data_algG_CF_3_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_3_filtered$third_base_CF_3 <- ifelse(apply(data_algG_CF_3_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_3_filtered$`PA14-Koordinate`)
  if(data_algG_CF_3_filtered$main_base_CF_3[a] == "percent_As_CF_3"){
    data_algG_CF_3_filtered$main_base_CF_3[a] <- "A"                     
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] == "percent_Cs_CF_3"){
    data_algG_CF_3_filtered$main_base_CF_3[a] <- "C"                     
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] == "percent_Gs_CF_3"){
    data_algG_CF_3_filtered$main_base_CF_3[a] <- "G"                     
  } else{
    data_algG_CF_3_filtered$main_base_CF_3[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_3_filtered$`PA14-Koordinate`)
  if(data_algG_CF_3_filtered$second_base_CF_3[a] == "percent_As_CF_3"){
    data_algG_CF_3_filtered$second_base_CF_3[a] <- "A"                     
  } else if(data_algG_CF_3_filtered$second_base_CF_3[a] == "percent_Cs_CF_3"){
    data_algG_CF_3_filtered$second_base_CF_3[a] <- "C"                     
  } else if(data_algG_CF_3_filtered$second_base_CF_3[a] == "percent_Gs_CF_3"){
    data_algG_CF_3_filtered$second_base_CF_3[a] <- "G"                     
  } else{
    data_algG_CF_3_filtered$second_base_CF_3[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_3_filtered$`PA14-Koordinate`)
  if(data_algG_CF_3_filtered$main_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$third_base_CF_3[a] == 0){
    data_algG_CF_3_filtered$third_base_CF_3[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_3_filtered$`PA14-Koordinate`)
  if(data_algG_CF_3_filtered$third_base_CF_3[a] == 1){
    data_algG_CF_3_filtered_third_variant <- data_algG_CF_3_filtered[a,]
    data_algG_CF_3_filtered$third_base_CF_3[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_3_filtered$`PA14-Koordinate`)
  if(data_algG_CF_3_filtered$third_base_CF_3[a] == 2){
    data_algG_CF_3_filtered_third_variant_no_ref <- data_algG_CF_3_filtered[a,]
    data_algG_CF_3_filtered$third_base_CF_3[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_3_filtered$third_base_CF_3) == 0, data_algG_CF_3_filtered <- subset(data_algG_CF_3_filtered, select = -c(third_base_CF_3)),data_algG_CF_3_filtered <- data_algG_CF_3_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_filtered_third_variant$third_percentage <- apply(data_algG_CF_3_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_3_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_filtered_third_variant$third_base_CF_3 <- colnames(data_algG_CF_3_filtered_third_variant[, 8:11])[apply(data_algG_CF_3_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_As_CF_3","A", ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_Cs_CF_3","C", ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_Gs_CF_3","G", ifelse(data_algG_CF_3_only_SNP$main_base_CF_3 == "percent_Ts_CF_3","T",z<-0))))

ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,
       data_algG_CF_3_filtered_third_variant$third_base_CF_3 <- 
         ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "percent_As_CF_3","A", 
                ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "percent_Cs_CF_3","C",
                       ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "percent_Gs_CF_3","G",
                              ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "percent_Ts_CF_3","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,
       data_algG_CF_3_filtered_third_variant$third_base_CF_3 <- ifelse(
         data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14,data_algG_CF_3_filtered_third_variant$third_base_CF_3 <- data_algG_CF_3_filtered_third_variant$second_base_CF_3,data_algG_CF_3_filtered_third_variant$third_base_CF_3
       ),
       z<-0
)

ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,
       data_algG_CF_3_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14, data_algG_CF_3_filtered_third_variant$second_percentage,data_algG_CF_3_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_filtered_third_variant <- data_algG_CF_3_filtered_third_variant[,-which(names(data_algG_CF_3_filtered_third_variant) %in% c("second_percentage","second_base_CF_3"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_3_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_3_filtered_third_variant_no_ref)[names(data_algG_CF_3_filtered_third_variant_no_ref) == "second_base_CF_3"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_3_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_3_filtered_third_variant_no_ref)[names(data_algG_CF_3_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_3",z<-0)




# Define SNP percentage & Add it

data_algG_CF_3_filtered$SNP_base <- 0
data_algG_CF_3_filtered$SNP_percentage_CF_3 <- 0


for(i in data_algG_CF_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_3_filtered$`PA14-Koordinate`)
  if(data_algG_CF_3_filtered$main_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] == "A"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_As_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] == "C"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Cs_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] == "G"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Gs_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] == "T"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Ts_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_3_filtered$second_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "A"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_As_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_3_filtered$second_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "C"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Cs_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_3_filtered$second_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "G"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Gs_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_3_filtered$second_base_CF_3[a] == data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "T"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Ts_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "A"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_As_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "C"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Cs_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "G"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Gs_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_3_filtered$main_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$second_base_CF_3[a] != data_algG_CF_3_filtered$GenBase_PA14[a] & data_algG_CF_3_filtered$main_base_CF_3[a] == "T"){
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- data_algG_CF_3_filtered$percent_Ts_CF_3[a]
    data_algG_CF_3_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_3_filtered$SNP_percentage_CF_3[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_filtered_third_variant$SNP_percentage_CF_3 <- 0,z<-0)


ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,
       data_algG_CF_3_filtered_third_variant$SNP_percentage_CF_3 <- 
         ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "A",
                data_algG_CF_3_filtered_third_variant$percent_As_CF_3, 
                ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "C",
                       data_algG_CF_3_filtered_third_variant$percent_Cs_CF_3,
                       ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "G",
                              data_algG_CF_3_filtered_third_variant$percent_Gs_CF_3,
                              ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "T",
                                     data_algG_CF_3_filtered_third_variant$percent_Ts_CF_3,
                                     ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "A",
                                            data_algG_CF_3_filtered_third_variant$percent_As_CF_3,
                                            ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "C",
                                                   data_algG_CF_3_filtered_third_variant$percent_Cs_CF_3,
                                                   ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "G",
                                                          data_algG_CF_3_filtered_third_variant$percent_Gs_CF_3,
                                                          ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "T",
                                                                 data_algG_CF_3_filtered_third_variant$percent_Ts_CF_3,
                                                                 ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "A",
                                                                        data_algG_CF_3_filtered_third_variant$percent_As_CF_3,
                                                                        ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "C",
                                                                               data_algG_CF_3_filtered_third_variant$percent_Cs_CF_3,
                                                                               ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "G",
                                                                                      data_algG_CF_3_filtered_third_variant$percent_Gs_CF_3,
                                                                                      ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "T",
                                                                                             data_algG_CF_3_filtered_third_variant$percent_Ts_CF_3,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,
       data_algG_CF_3_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "A",
                "A", 
                ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "C",
                       "C",
                       ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "G",
                              "G",
                              ifelse(data_algG_CF_3_filtered_third_variant$main_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "T",
                                     "T",
                                     ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "A",
                                            "A",
                                            ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 == data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$main_base_CF_3 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_3_filtered_third_variant$third_base_CF_3 != data_algG_CF_3_filtered_third_variant$GenBase_PA14 & data_algG_CF_3_filtered_third_variant$third_base_CF_3 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_3_filtered$Isolates_CF_3 <- round((data_algG_CF_3_filtered$SNP_percentage_CF_3*100)/(100/isolates_CF_3))
ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_filtered_third_variant$Isolates_CF_3 <- round((data_algG_CF_3_filtered_third_variant$SNP_percentage_CF_3*100)/(100/isolates_CF_3)),z<-0)
ifelse(exists("data_algG_CF_3_filtered_third_variant_no_ref") == TRUE,data_algG_CF_3_filtered_third_variant_no_ref$Isolates_CF_3 <- round((data_algG_CF_3_filtered_third_variant_no_ref$SNP_percentage_CF_3*100)/(100/isolates_CF_3)),z<-0)



# Clean data

data_algG_CF_3_clean <- data_algG_CF_3_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_clean_third_base <- data_algG_CF_3_filtered_third_variant[,which(names(data_algG_CF_3_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_3","Isolates_CF_3"))],z<-0)
ifelse(exists("data_algG_CF_3_filtered_third_variant") == TRUE,data_algG_CF_3_clean <- rbind(data_algG_CF_3_clean,data_algG_CF_3_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_3_filtered_third_variant_no_ref") == TRUE,data_algG_CF_3_clean_third_base_no_ref <- data_algG_CF_3_filtered_third_variant_no_ref[,which(names(data_algG_CF_3_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_3","Isolates_CF_3"))],z<-0)
ifelse(exists("data_algG_CF_3_filtered_third_variant_no_ref") == TRUE,data_algG_CF_3_clean <- rbind(data_algG_CF_3_clean,data_algG_CF_3_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_3_clean_20_filter <- merge(data_algG_CF_3_clean, data_algG_CF_3_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_3_clean_cleared <- data_algG_CF_3_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_3_clean_20_filter)){
  if(data_algG_CF_3_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_3_clean_20_filter$As_CF_3[i] >= 20){
    data_algG_CF_3_clean_cleared <- rbind(data_algG_CF_3_clean_cleared,data_algG_CF_3_clean_20_filter[i,])
  } else if(data_algG_CF_3_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_3_clean_20_filter$Cs_CF_3[i] >= 20){
    data_algG_CF_3_clean_cleared <- rbind(data_algG_CF_3_clean_cleared,data_algG_CF_3_clean_20_filter[i,])
  } else if(data_algG_CF_3_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_3_clean_20_filter$Gs_CF_3[i] >= 20){
    data_algG_CF_3_clean_cleared <- rbind(data_algG_CF_3_clean_cleared,data_algG_CF_3_clean_20_filter[i,])
  }else if(data_algG_CF_3_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_3_clean_20_filter$Ts_CF_3[i] >= 20){
    data_algG_CF_3_clean_cleared <- rbind(data_algG_CF_3_clean_cleared,data_algG_CF_3_clean_20_filter[i,])
  }
}

data_algG_CF_3_clean <- data_algG_CF_3_clean_cleared[,c(1:5)]


### CF_4

# Input 
data_algG_CF_4 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T4/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_4_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T4/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_4_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T4/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_4 <- data_algG_CF_4[!(data_algG_CF_4$Position %in% primer_correction_algG_CF_4_fwd$Position),]
data_algG_CF_4 <- data_algG_CF_4[!(data_algG_CF_4$Position %in% primer_correction_algG_CF_4_rev$Position),]
data_algG_CF_4 <- rbind(data_algG_CF_4,primer_correction_algG_CF_4_fwd,primer_correction_algG_CF_4_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_4 <- subset(data_algG_CF_4, coverage_CF_4 >199)
data_algG_CF_4_5UTR <- subset(data_algG_CF_4_5UTR, coverage_CF_4_5UTR >199)
data_algG_CF_4_3UTR <- subset(data_algG_CF_4_3UTR, coverage_CF_4_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_4_only_SNP <- subset(data_algG_CF_4, second_percentage <= 0.025) 
data_algG_CF_4_only_SNP$main_base_CF_4 <- colnames(data_algG_CF_4_only_SNP[, 8:11])[apply(data_algG_CF_4_only_SNP[, 8:11],1,which.max)]
data_algG_CF_4_only_SNP$main_base_CF_4 <- ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_As_CF_4","A", ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_Cs_CF_4","C", ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_Gs_CF_4","G", ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_Ts_CF_4","T",z<-0))))
data_algG_CF_4_only_SNP <- subset(data_algG_CF_4_only_SNP,GenBase_PA14 != main_base_CF_4)
ifelse(nrow(data_algG_CF_4_only_SNP) > 0,data_algG_CF_4_only_SNP$Isolates_CF_4 <- isolates_CF_4, z <- 0)


ifelse(nrow(data_algG_CF_4_only_SNP) > 0,data_algG_CF_4_only_SNP_clean <- data_algG_CF_4_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_4_only_SNP) > 0,names(data_algG_CF_4_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_4","Isolates_CF_4"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_4_filtered <- subset(data_algG_CF_4, second_percentage >= 0.025) 
data_algG_CF_4_filtered$main_base_CF_4 <- colnames(data_algG_CF_4_filtered[, 8:11])[apply(data_algG_CF_4_filtered[, 8:11],1,which.max)]
data_algG_CF_4_filtered$second_base_CF_4 <- colnames(data_algG_CF_4_filtered[, 8:11])[apply(data_algG_CF_4_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_4_filtered$third_base_CF_4 <- ifelse(apply(data_algG_CF_4_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_4_filtered$`PA14-Koordinate`)
  if(data_algG_CF_4_filtered$main_base_CF_4[a] == "percent_As_CF_4"){
    data_algG_CF_4_filtered$main_base_CF_4[a] <- "A"                     
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] == "percent_Cs_CF_4"){
    data_algG_CF_4_filtered$main_base_CF_4[a] <- "C"                     
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] == "percent_Gs_CF_4"){
    data_algG_CF_4_filtered$main_base_CF_4[a] <- "G"                     
  } else{
    data_algG_CF_4_filtered$main_base_CF_4[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_4_filtered$`PA14-Koordinate`)
  if(data_algG_CF_4_filtered$second_base_CF_4[a] == "percent_As_CF_4"){
    data_algG_CF_4_filtered$second_base_CF_4[a] <- "A"                     
  } else if(data_algG_CF_4_filtered$second_base_CF_4[a] == "percent_Cs_CF_4"){
    data_algG_CF_4_filtered$second_base_CF_4[a] <- "C"                     
  } else if(data_algG_CF_4_filtered$second_base_CF_4[a] == "percent_Gs_CF_4"){
    data_algG_CF_4_filtered$second_base_CF_4[a] <- "G"                     
  } else{
    data_algG_CF_4_filtered$second_base_CF_4[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_4_filtered$`PA14-Koordinate`)
  if(data_algG_CF_4_filtered$main_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$third_base_CF_4[a] == 0){
    data_algG_CF_4_filtered$third_base_CF_4[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_4_filtered$`PA14-Koordinate`)
  if(data_algG_CF_4_filtered$third_base_CF_4[a] == 1){
    data_algG_CF_4_filtered_third_variant <- data_algG_CF_4_filtered[a,]
    data_algG_CF_4_filtered$third_base_CF_4[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_4_filtered$`PA14-Koordinate`)
  if(data_algG_CF_4_filtered$third_base_CF_4[a] == 2){
    data_algG_CF_4_filtered_third_variant_no_ref <- data_algG_CF_4_filtered[a,]
    data_algG_CF_4_filtered$third_base_CF_4[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_4_filtered$third_base_CF_4) == 0, data_algG_CF_4_filtered <- subset(data_algG_CF_4_filtered, select = -c(third_base_CF_4)),data_algG_CF_4_filtered <- data_algG_CF_4_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_filtered_third_variant$third_percentage <- apply(data_algG_CF_4_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_4_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_filtered_third_variant$third_base_CF_4 <- colnames(data_algG_CF_4_filtered_third_variant[, 8:11])[apply(data_algG_CF_4_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_As_CF_4","A", ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_Cs_CF_4","C", ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_Gs_CF_4","G", ifelse(data_algG_CF_4_only_SNP$main_base_CF_4 == "percent_Ts_CF_4","T",z<-0))))

ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,
       data_algG_CF_4_filtered_third_variant$third_base_CF_4 <- 
         ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "percent_As_CF_4","A", 
                ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "percent_Cs_CF_4","C",
                       ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "percent_Gs_CF_4","G",
                              ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "percent_Ts_CF_4","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,
       data_algG_CF_4_filtered_third_variant$third_base_CF_4 <- ifelse(
         data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14,data_algG_CF_4_filtered_third_variant$third_base_CF_4 <- data_algG_CF_4_filtered_third_variant$second_base_CF_4,data_algG_CF_4_filtered_third_variant$third_base_CF_4
       ),
       z<-0
)

ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,
       data_algG_CF_4_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14, data_algG_CF_4_filtered_third_variant$second_percentage,data_algG_CF_4_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_filtered_third_variant <- data_algG_CF_4_filtered_third_variant[,-which(names(data_algG_CF_4_filtered_third_variant) %in% c("second_percentage","second_base_CF_4"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_4_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_4_filtered_third_variant_no_ref)[names(data_algG_CF_4_filtered_third_variant_no_ref) == "second_base_CF_4"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_4_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_4_filtered_third_variant_no_ref)[names(data_algG_CF_4_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_4",z<-0)




# Define SNP percentage & Add it

data_algG_CF_4_filtered$SNP_base <- 0
data_algG_CF_4_filtered$SNP_percentage_CF_4 <- 0


for(i in data_algG_CF_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_4_filtered$`PA14-Koordinate`)
  if(data_algG_CF_4_filtered$main_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] == "A"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_As_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] == "C"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Cs_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] == "G"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Gs_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] == "T"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Ts_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_4_filtered$second_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "A"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_As_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_4_filtered$second_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "C"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Cs_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_4_filtered$second_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "G"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Gs_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_4_filtered$second_base_CF_4[a] == data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "T"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Ts_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "A"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_As_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "C"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Cs_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "G"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Gs_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_4_filtered$main_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$second_base_CF_4[a] != data_algG_CF_4_filtered$GenBase_PA14[a] & data_algG_CF_4_filtered$main_base_CF_4[a] == "T"){
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- data_algG_CF_4_filtered$percent_Ts_CF_4[a]
    data_algG_CF_4_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_4_filtered$SNP_percentage_CF_4[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_filtered_third_variant$SNP_percentage_CF_4 <- 0,z<-0)


ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,
       data_algG_CF_4_filtered_third_variant$SNP_percentage_CF_4 <- 
         ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "A",
                data_algG_CF_4_filtered_third_variant$percent_As_CF_4, 
                ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "C",
                       data_algG_CF_4_filtered_third_variant$percent_Cs_CF_4,
                       ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "G",
                              data_algG_CF_4_filtered_third_variant$percent_Gs_CF_4,
                              ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "T",
                                     data_algG_CF_4_filtered_third_variant$percent_Ts_CF_4,
                                     ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "A",
                                            data_algG_CF_4_filtered_third_variant$percent_As_CF_4,
                                            ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "C",
                                                   data_algG_CF_4_filtered_third_variant$percent_Cs_CF_4,
                                                   ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "G",
                                                          data_algG_CF_4_filtered_third_variant$percent_Gs_CF_4,
                                                          ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "T",
                                                                 data_algG_CF_4_filtered_third_variant$percent_Ts_CF_4,
                                                                 ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "A",
                                                                        data_algG_CF_4_filtered_third_variant$percent_As_CF_4,
                                                                        ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "C",
                                                                               data_algG_CF_4_filtered_third_variant$percent_Cs_CF_4,
                                                                               ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "G",
                                                                                      data_algG_CF_4_filtered_third_variant$percent_Gs_CF_4,
                                                                                      ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "T",
                                                                                             data_algG_CF_4_filtered_third_variant$percent_Ts_CF_4,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,
       data_algG_CF_4_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "A",
                "A", 
                ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "C",
                       "C",
                       ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "G",
                              "G",
                              ifelse(data_algG_CF_4_filtered_third_variant$main_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "T",
                                     "T",
                                     ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "A",
                                            "A",
                                            ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 == data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$main_base_CF_4 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_4_filtered_third_variant$third_base_CF_4 != data_algG_CF_4_filtered_third_variant$GenBase_PA14 & data_algG_CF_4_filtered_third_variant$third_base_CF_4 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_4_filtered$Isolates_CF_4 <- round((data_algG_CF_4_filtered$SNP_percentage_CF_4*100)/(100/isolates_CF_4))
ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_filtered_third_variant$Isolates_CF_4 <- round((data_algG_CF_4_filtered_third_variant$SNP_percentage_CF_4*100)/(100/isolates_CF_4)),z<-0)
ifelse(exists("data_algG_CF_4_filtered_third_variant_no_ref") == TRUE,data_algG_CF_4_filtered_third_variant_no_ref$Isolates_CF_4 <- round((data_algG_CF_4_filtered_third_variant_no_ref$SNP_percentage_CF_4*100)/(100/isolates_CF_4)),z<-0)



# Clean data

data_algG_CF_4_clean <- data_algG_CF_4_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_clean_third_base <- data_algG_CF_4_filtered_third_variant[,which(names(data_algG_CF_4_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_4","Isolates_CF_4"))],z<-0)
ifelse(exists("data_algG_CF_4_filtered_third_variant") == TRUE,data_algG_CF_4_clean <- rbind(data_algG_CF_4_clean,data_algG_CF_4_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_4_filtered_third_variant_no_ref") == TRUE,data_algG_CF_4_clean_third_base_no_ref <- data_algG_CF_4_filtered_third_variant_no_ref[,which(names(data_algG_CF_4_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_4","Isolates_CF_4"))],z<-0)
ifelse(exists("data_algG_CF_4_filtered_third_variant_no_ref") == TRUE,data_algG_CF_4_clean <- rbind(data_algG_CF_4_clean,data_algG_CF_4_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_4_clean_20_filter <- merge(data_algG_CF_4_clean, data_algG_CF_4_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_4_clean_cleared <- data_algG_CF_4_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_4_clean_20_filter)){
  if(data_algG_CF_4_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_4_clean_20_filter$As_CF_4[i] >= 20){
    data_algG_CF_4_clean_cleared <- rbind(data_algG_CF_4_clean_cleared,data_algG_CF_4_clean_20_filter[i,])
  } else if(data_algG_CF_4_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_4_clean_20_filter$Cs_CF_4[i] >= 20){
    data_algG_CF_4_clean_cleared <- rbind(data_algG_CF_4_clean_cleared,data_algG_CF_4_clean_20_filter[i,])
  } else if(data_algG_CF_4_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_4_clean_20_filter$Gs_CF_4[i] >= 20){
    data_algG_CF_4_clean_cleared <- rbind(data_algG_CF_4_clean_cleared,data_algG_CF_4_clean_20_filter[i,])
  }else if(data_algG_CF_4_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_4_clean_20_filter$Ts_CF_4[i] >= 20){
    data_algG_CF_4_clean_cleared <- rbind(data_algG_CF_4_clean_cleared,data_algG_CF_4_clean_20_filter[i,])
  }
}
data_algG_CF_4_clean <- data_algG_CF_4_clean_cleared[,c(1:5)]


### CF_5

# Input 
data_algG_CF_5 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T5/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_5_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T5/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_5_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T5/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_5 <- data_algG_CF_5[!(data_algG_CF_5$Position %in% primer_correction_algG_CF_5_fwd$Position),]
data_algG_CF_5 <- data_algG_CF_5[!(data_algG_CF_5$Position %in% primer_correction_algG_CF_5_rev$Position),]
data_algG_CF_5 <- rbind(data_algG_CF_5,primer_correction_algG_CF_5_fwd,primer_correction_algG_CF_5_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_5 <- subset(data_algG_CF_5, coverage_CF_5 >199)
data_algG_CF_5_5UTR <- subset(data_algG_CF_5_5UTR, coverage_CF_5_5UTR >199)
data_algG_CF_5_3UTR <- subset(data_algG_CF_5_3UTR, coverage_CF_5_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_5_only_SNP <- subset(data_algG_CF_5, second_percentage <= 0.025) 
data_algG_CF_5_only_SNP$main_base_CF_5 <- colnames(data_algG_CF_5_only_SNP[, 8:11])[apply(data_algG_CF_5_only_SNP[, 8:11],1,which.max)]
data_algG_CF_5_only_SNP$main_base_CF_5 <- ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_As_CF_5","A", ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_Cs_CF_5","C", ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_Gs_CF_5","G", ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_Ts_CF_5","T",z<-0))))
data_algG_CF_5_only_SNP <- subset(data_algG_CF_5_only_SNP,GenBase_PA14 != main_base_CF_5)
ifelse(nrow(data_algG_CF_5_only_SNP) > 0,data_algG_CF_5_only_SNP$Isolates_CF_5 <- isolates_CF_5, z <- 0)


ifelse(nrow(data_algG_CF_5_only_SNP) > 0,data_algG_CF_5_only_SNP_clean <- data_algG_CF_5_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_5_only_SNP) > 0,names(data_algG_CF_5_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_5","Isolates_CF_5"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_5_filtered <- subset(data_algG_CF_5, second_percentage >= 0.025) 
data_algG_CF_5_filtered$main_base_CF_5 <- colnames(data_algG_CF_5_filtered[, 8:11])[apply(data_algG_CF_5_filtered[, 8:11],1,which.max)]
data_algG_CF_5_filtered$second_base_CF_5 <- colnames(data_algG_CF_5_filtered[, 8:11])[apply(data_algG_CF_5_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_5_filtered$third_base_CF_5 <- ifelse(apply(data_algG_CF_5_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_5_filtered$`PA14-Koordinate`)
  if(data_algG_CF_5_filtered$main_base_CF_5[a] == "percent_As_CF_5"){
    data_algG_CF_5_filtered$main_base_CF_5[a] <- "A"                     
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] == "percent_Cs_CF_5"){
    data_algG_CF_5_filtered$main_base_CF_5[a] <- "C"                     
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] == "percent_Gs_CF_5"){
    data_algG_CF_5_filtered$main_base_CF_5[a] <- "G"                     
  } else{
    data_algG_CF_5_filtered$main_base_CF_5[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_5_filtered$`PA14-Koordinate`)
  if(data_algG_CF_5_filtered$second_base_CF_5[a] == "percent_As_CF_5"){
    data_algG_CF_5_filtered$second_base_CF_5[a] <- "A"                     
  } else if(data_algG_CF_5_filtered$second_base_CF_5[a] == "percent_Cs_CF_5"){
    data_algG_CF_5_filtered$second_base_CF_5[a] <- "C"                     
  } else if(data_algG_CF_5_filtered$second_base_CF_5[a] == "percent_Gs_CF_5"){
    data_algG_CF_5_filtered$second_base_CF_5[a] <- "G"                     
  } else{
    data_algG_CF_5_filtered$second_base_CF_5[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_5_filtered$`PA14-Koordinate`)
  if(data_algG_CF_5_filtered$main_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$third_base_CF_5[a] == 0){
    data_algG_CF_5_filtered$third_base_CF_5[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_5_filtered$`PA14-Koordinate`)
  if(data_algG_CF_5_filtered$third_base_CF_5[a] == 1){
    data_algG_CF_5_filtered_third_variant <- data_algG_CF_5_filtered[a,]
    data_algG_CF_5_filtered$third_base_CF_5[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_5_filtered$`PA14-Koordinate`)
  if(data_algG_CF_5_filtered$third_base_CF_5[a] == 2){
    data_algG_CF_5_filtered_third_variant_no_ref <- data_algG_CF_5_filtered[a,]
    data_algG_CF_5_filtered$third_base_CF_5[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_5_filtered$third_base_CF_5) == 0, data_algG_CF_5_filtered <- subset(data_algG_CF_5_filtered, select = -c(third_base_CF_5)),data_algG_CF_5_filtered <- data_algG_CF_5_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_filtered_third_variant$third_percentage <- apply(data_algG_CF_5_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_5_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_filtered_third_variant$third_base_CF_5 <- colnames(data_algG_CF_5_filtered_third_variant[, 8:11])[apply(data_algG_CF_5_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_As_CF_5","A", ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_Cs_CF_5","C", ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_Gs_CF_5","G", ifelse(data_algG_CF_5_only_SNP$main_base_CF_5 == "percent_Ts_CF_5","T",z<-0))))

ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,
       data_algG_CF_5_filtered_third_variant$third_base_CF_5 <- 
         ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "percent_As_CF_5","A", 
                ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "percent_Cs_CF_5","C",
                       ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "percent_Gs_CF_5","G",
                              ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "percent_Ts_CF_5","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,
       data_algG_CF_5_filtered_third_variant$third_base_CF_5 <- ifelse(
         data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14,data_algG_CF_5_filtered_third_variant$third_base_CF_5 <- data_algG_CF_5_filtered_third_variant$second_base_CF_5,data_algG_CF_5_filtered_third_variant$third_base_CF_5
       ),
       z<-0
)

ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,
       data_algG_CF_5_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14, data_algG_CF_5_filtered_third_variant$second_percentage,data_algG_CF_5_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_filtered_third_variant <- data_algG_CF_5_filtered_third_variant[,-which(names(data_algG_CF_5_filtered_third_variant) %in% c("second_percentage","second_base_CF_5"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_5_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_5_filtered_third_variant_no_ref)[names(data_algG_CF_5_filtered_third_variant_no_ref) == "second_base_CF_5"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_5_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_5_filtered_third_variant_no_ref)[names(data_algG_CF_5_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_5",z<-0)




# Define SNP percentage & Add it

data_algG_CF_5_filtered$SNP_base <- 0
data_algG_CF_5_filtered$SNP_percentage_CF_5 <- 0


for(i in data_algG_CF_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_5_filtered$`PA14-Koordinate`)
  if(data_algG_CF_5_filtered$main_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] == "A"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_As_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] == "C"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Cs_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] == "G"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Gs_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] == "T"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Ts_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_5_filtered$second_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "A"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_As_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_5_filtered$second_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "C"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Cs_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_5_filtered$second_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "G"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Gs_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_5_filtered$second_base_CF_5[a] == data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "T"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Ts_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "A"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_As_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "C"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Cs_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "G"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Gs_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_5_filtered$main_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$second_base_CF_5[a] != data_algG_CF_5_filtered$GenBase_PA14[a] & data_algG_CF_5_filtered$main_base_CF_5[a] == "T"){
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- data_algG_CF_5_filtered$percent_Ts_CF_5[a]
    data_algG_CF_5_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_5_filtered$SNP_percentage_CF_5[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_filtered_third_variant$SNP_percentage_CF_5 <- 0,z<-0)


ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,
       data_algG_CF_5_filtered_third_variant$SNP_percentage_CF_5 <- 
         ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "A",
                data_algG_CF_5_filtered_third_variant$percent_As_CF_5, 
                ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "C",
                       data_algG_CF_5_filtered_third_variant$percent_Cs_CF_5,
                       ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "G",
                              data_algG_CF_5_filtered_third_variant$percent_Gs_CF_5,
                              ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "T",
                                     data_algG_CF_5_filtered_third_variant$percent_Ts_CF_5,
                                     ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "A",
                                            data_algG_CF_5_filtered_third_variant$percent_As_CF_5,
                                            ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "C",
                                                   data_algG_CF_5_filtered_third_variant$percent_Cs_CF_5,
                                                   ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "G",
                                                          data_algG_CF_5_filtered_third_variant$percent_Gs_CF_5,
                                                          ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "T",
                                                                 data_algG_CF_5_filtered_third_variant$percent_Ts_CF_5,
                                                                 ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "A",
                                                                        data_algG_CF_5_filtered_third_variant$percent_As_CF_5,
                                                                        ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "C",
                                                                               data_algG_CF_5_filtered_third_variant$percent_Cs_CF_5,
                                                                               ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "G",
                                                                                      data_algG_CF_5_filtered_third_variant$percent_Gs_CF_5,
                                                                                      ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "T",
                                                                                             data_algG_CF_5_filtered_third_variant$percent_Ts_CF_5,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,
       data_algG_CF_5_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "A",
                "A", 
                ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "C",
                       "C",
                       ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "G",
                              "G",
                              ifelse(data_algG_CF_5_filtered_third_variant$main_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "T",
                                     "T",
                                     ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "A",
                                            "A",
                                            ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 == data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$main_base_CF_5 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_5_filtered_third_variant$third_base_CF_5 != data_algG_CF_5_filtered_third_variant$GenBase_PA14 & data_algG_CF_5_filtered_third_variant$third_base_CF_5 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_5_filtered$Isolates_CF_5 <- round((data_algG_CF_5_filtered$SNP_percentage_CF_5*100)/(100/isolates_CF_5))
ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_filtered_third_variant$Isolates_CF_5 <- round((data_algG_CF_5_filtered_third_variant$SNP_percentage_CF_5*100)/(100/isolates_CF_5)),z<-0)
ifelse(exists("data_algG_CF_5_filtered_third_variant_no_ref") == TRUE,data_algG_CF_5_filtered_third_variant_no_ref$Isolates_CF_5 <- round((data_algG_CF_5_filtered_third_variant_no_ref$SNP_percentage_CF_5*100)/(100/isolates_CF_5)),z<-0)



# Clean data

data_algG_CF_5_clean <- data_algG_CF_5_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_clean_third_base <- data_algG_CF_5_filtered_third_variant[,which(names(data_algG_CF_5_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_5","Isolates_CF_5"))],z<-0)
ifelse(exists("data_algG_CF_5_filtered_third_variant") == TRUE,data_algG_CF_5_clean <- rbind(data_algG_CF_5_clean,data_algG_CF_5_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_5_filtered_third_variant_no_ref") == TRUE,data_algG_CF_5_clean_third_base_no_ref <- data_algG_CF_5_filtered_third_variant_no_ref[,which(names(data_algG_CF_5_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_5","Isolates_CF_5"))],z<-0)
ifelse(exists("data_algG_CF_5_filtered_third_variant_no_ref") == TRUE,data_algG_CF_5_clean <- rbind(data_algG_CF_5_clean,data_algG_CF_5_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_5_clean_20_filter <- merge(data_algG_CF_5_clean, data_algG_CF_5_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_5_clean_cleared <- data_algG_CF_5_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_5_clean_20_filter)){
  if(data_algG_CF_5_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_5_clean_20_filter$As_CF_5[i] >= 20){
    data_algG_CF_5_clean_cleared <- rbind(data_algG_CF_5_clean_cleared,data_algG_CF_5_clean_20_filter[i,])
  } else if(data_algG_CF_5_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_5_clean_20_filter$Cs_CF_5[i] >= 20){
    data_algG_CF_5_clean_cleared <- rbind(data_algG_CF_5_clean_cleared,data_algG_CF_5_clean_20_filter[i,])
  } else if(data_algG_CF_5_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_5_clean_20_filter$Gs_CF_5[i] >= 20){
    data_algG_CF_5_clean_cleared <- rbind(data_algG_CF_5_clean_cleared,data_algG_CF_5_clean_20_filter[i,])
  }else if(data_algG_CF_5_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_5_clean_20_filter$Ts_CF_5[i] >= 20){
    data_algG_CF_5_clean_cleared <- rbind(data_algG_CF_5_clean_cleared,data_algG_CF_5_clean_20_filter[i,])
  }
}

data_algG_CF_5_clean <- data_algG_CF_5_clean_cleared[,c(1:5)]


### CF_6

# Input 
data_algG_CF_6 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T6/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_6_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T6/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_6_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T6/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_6 <- data_algG_CF_6[!(data_algG_CF_6$Position %in% primer_correction_algG_CF_6_fwd$Position),]
data_algG_CF_6 <- data_algG_CF_6[!(data_algG_CF_6$Position %in% primer_correction_algG_CF_6_rev$Position),]
data_algG_CF_6 <- rbind(data_algG_CF_6,primer_correction_algG_CF_6_fwd,primer_correction_algG_CF_6_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_6 <- subset(data_algG_CF_6, coverage_CF_6 >199)
data_algG_CF_6_5UTR <- subset(data_algG_CF_6_5UTR, coverage_CF_6_5UTR >199)
data_algG_CF_6_3UTR <- subset(data_algG_CF_6_3UTR, coverage_CF_6_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_6_only_SNP <- subset(data_algG_CF_6, second_percentage <= 0.025) 
data_algG_CF_6_only_SNP$main_base_CF_6 <- colnames(data_algG_CF_6_only_SNP[, 8:11])[apply(data_algG_CF_6_only_SNP[, 8:11],1,which.max)]
data_algG_CF_6_only_SNP$main_base_CF_6 <- ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_As_CF_6","A", ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_Cs_CF_6","C", ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_Gs_CF_6","G", ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_Ts_CF_6","T",z<-0))))
data_algG_CF_6_only_SNP <- subset(data_algG_CF_6_only_SNP,GenBase_PA14 != main_base_CF_6)
ifelse(nrow(data_algG_CF_6_only_SNP) > 0,data_algG_CF_6_only_SNP$Isolates_CF_6 <- isolates_CF_6, z <- 0)


ifelse(nrow(data_algG_CF_6_only_SNP) > 0,data_algG_CF_6_only_SNP_clean <- data_algG_CF_6_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_6_only_SNP) > 0,names(data_algG_CF_6_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_6","Isolates_CF_6"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_6_filtered <- subset(data_algG_CF_6, second_percentage >= 0.025) 
data_algG_CF_6_filtered$main_base_CF_6 <- colnames(data_algG_CF_6_filtered[, 8:11])[apply(data_algG_CF_6_filtered[, 8:11],1,which.max)]
data_algG_CF_6_filtered$second_base_CF_6 <- colnames(data_algG_CF_6_filtered[, 8:11])[apply(data_algG_CF_6_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_6_filtered$third_base_CF_6 <- ifelse(apply(data_algG_CF_6_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_6_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_6_filtered$`PA14-Koordinate`)
  if(data_algG_CF_6_filtered$main_base_CF_6[a] == "percent_As_CF_6"){
    data_algG_CF_6_filtered$main_base_CF_6[a] <- "A"                     
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] == "percent_Cs_CF_6"){
    data_algG_CF_6_filtered$main_base_CF_6[a] <- "C"                     
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] == "percent_Gs_CF_6"){
    data_algG_CF_6_filtered$main_base_CF_6[a] <- "G"                     
  } else{
    data_algG_CF_6_filtered$main_base_CF_6[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_6_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_6_filtered$`PA14-Koordinate`)
  if(data_algG_CF_6_filtered$second_base_CF_6[a] == "percent_As_CF_6"){
    data_algG_CF_6_filtered$second_base_CF_6[a] <- "A"                     
  } else if(data_algG_CF_6_filtered$second_base_CF_6[a] == "percent_Cs_CF_6"){
    data_algG_CF_6_filtered$second_base_CF_6[a] <- "C"                     
  } else if(data_algG_CF_6_filtered$second_base_CF_6[a] == "percent_Gs_CF_6"){
    data_algG_CF_6_filtered$second_base_CF_6[a] <- "G"                     
  } else{
    data_algG_CF_6_filtered$second_base_CF_6[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_6_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_6_filtered$`PA14-Koordinate`)
  if(data_algG_CF_6_filtered$main_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$third_base_CF_6[a] == 0){
    data_algG_CF_6_filtered$third_base_CF_6[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_6_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_6_filtered$`PA14-Koordinate`)
  if(data_algG_CF_6_filtered$third_base_CF_6[a] == 1){
    data_algG_CF_6_filtered_third_variant <- data_algG_CF_6_filtered[a,]
    data_algG_CF_6_filtered$third_base_CF_6[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_6_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_6_filtered$`PA14-Koordinate`)
  if(data_algG_CF_6_filtered$third_base_CF_6[a] == 2){
    data_algG_CF_6_filtered_third_variant_no_ref <- data_algG_CF_6_filtered[a,]
    data_algG_CF_6_filtered$third_base_CF_6[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_6_filtered$third_base_CF_6) == 0, data_algG_CF_6_filtered <- subset(data_algG_CF_6_filtered, select = -c(third_base_CF_6)),data_algG_CF_6_filtered <- data_algG_CF_6_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_filtered_third_variant$third_percentage <- apply(data_algG_CF_6_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_6_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_filtered_third_variant$third_base_CF_6 <- colnames(data_algG_CF_6_filtered_third_variant[, 8:11])[apply(data_algG_CF_6_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_As_CF_6","A", ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_Cs_CF_6","C", ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_Gs_CF_6","G", ifelse(data_algG_CF_6_only_SNP$main_base_CF_6 == "percent_Ts_CF_6","T",z<-0))))

ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,
       data_algG_CF_6_filtered_third_variant$third_base_CF_6 <- 
         ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "percent_As_CF_6","A", 
                ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "percent_Cs_CF_6","C",
                       ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "percent_Gs_CF_6","G",
                              ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "percent_Ts_CF_6","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,
       data_algG_CF_6_filtered_third_variant$third_base_CF_6 <- ifelse(
         data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14,data_algG_CF_6_filtered_third_variant$third_base_CF_6 <- data_algG_CF_6_filtered_third_variant$second_base_CF_6,data_algG_CF_6_filtered_third_variant$third_base_CF_6
       ),
       z<-0
)

ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,
       data_algG_CF_6_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14, data_algG_CF_6_filtered_third_variant$second_percentage,data_algG_CF_6_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_filtered_third_variant <- data_algG_CF_6_filtered_third_variant[,-which(names(data_algG_CF_6_filtered_third_variant) %in% c("second_percentage","second_base_CF_6"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_6_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_6_filtered_third_variant_no_ref)[names(data_algG_CF_6_filtered_third_variant_no_ref) == "second_base_CF_6"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_6_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_6_filtered_third_variant_no_ref)[names(data_algG_CF_6_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_6",z<-0)




# Define SNP percentage & Add it

data_algG_CF_6_filtered$SNP_base <- 0
data_algG_CF_6_filtered$SNP_percentage_CF_6 <- 0


for(i in data_algG_CF_6_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_6_filtered$`PA14-Koordinate`)
  if(data_algG_CF_6_filtered$main_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] == "A"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_As_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] == "C"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Cs_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] == "G"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Gs_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] == "T"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Ts_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_6_filtered$second_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "A"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_As_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_6_filtered$second_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "C"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Cs_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_6_filtered$second_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "G"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Gs_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_6_filtered$second_base_CF_6[a] == data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "T"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Ts_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "A"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_As_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "C"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Cs_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "G"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Gs_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_6_filtered$main_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$second_base_CF_6[a] != data_algG_CF_6_filtered$GenBase_PA14[a] & data_algG_CF_6_filtered$main_base_CF_6[a] == "T"){
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- data_algG_CF_6_filtered$percent_Ts_CF_6[a]
    data_algG_CF_6_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_6_filtered$SNP_percentage_CF_6[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_filtered_third_variant$SNP_percentage_CF_6 <- 0,z<-0)


ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,
       data_algG_CF_6_filtered_third_variant$SNP_percentage_CF_6 <- 
         ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "A",
                data_algG_CF_6_filtered_third_variant$percent_As_CF_6, 
                ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "C",
                       data_algG_CF_6_filtered_third_variant$percent_Cs_CF_6,
                       ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "G",
                              data_algG_CF_6_filtered_third_variant$percent_Gs_CF_6,
                              ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "T",
                                     data_algG_CF_6_filtered_third_variant$percent_Ts_CF_6,
                                     ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "A",
                                            data_algG_CF_6_filtered_third_variant$percent_As_CF_6,
                                            ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "C",
                                                   data_algG_CF_6_filtered_third_variant$percent_Cs_CF_6,
                                                   ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "G",
                                                          data_algG_CF_6_filtered_third_variant$percent_Gs_CF_6,
                                                          ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "T",
                                                                 data_algG_CF_6_filtered_third_variant$percent_Ts_CF_6,
                                                                 ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "A",
                                                                        data_algG_CF_6_filtered_third_variant$percent_As_CF_6,
                                                                        ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "C",
                                                                               data_algG_CF_6_filtered_third_variant$percent_Cs_CF_6,
                                                                               ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "G",
                                                                                      data_algG_CF_6_filtered_third_variant$percent_Gs_CF_6,
                                                                                      ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "T",
                                                                                             data_algG_CF_6_filtered_third_variant$percent_Ts_CF_6,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,
       data_algG_CF_6_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "A",
                "A", 
                ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "C",
                       "C",
                       ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "G",
                              "G",
                              ifelse(data_algG_CF_6_filtered_third_variant$main_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "T",
                                     "T",
                                     ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "A",
                                            "A",
                                            ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 == data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$main_base_CF_6 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_6_filtered_third_variant$third_base_CF_6 != data_algG_CF_6_filtered_third_variant$GenBase_PA14 & data_algG_CF_6_filtered_third_variant$third_base_CF_6 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_6_filtered$Isolates_CF_6 <- round((data_algG_CF_6_filtered$SNP_percentage_CF_6*100)/(100/isolates_CF_6))
ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_filtered_third_variant$Isolates_CF_6 <- round((data_algG_CF_6_filtered_third_variant$SNP_percentage_CF_6*100)/(100/isolates_CF_6)),z<-0)
ifelse(exists("data_algG_CF_6_filtered_third_variant_no_ref") == TRUE,data_algG_CF_6_filtered_third_variant_no_ref$Isolates_CF_6 <- round((data_algG_CF_6_filtered_third_variant_no_ref$SNP_percentage_CF_6*100)/(100/isolates_CF_6)),z<-0)



# Clean data

data_algG_CF_6_clean <- data_algG_CF_6_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_clean_third_base <- data_algG_CF_6_filtered_third_variant[,which(names(data_algG_CF_6_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_6","Isolates_CF_6"))],z<-0)
ifelse(exists("data_algG_CF_6_filtered_third_variant") == TRUE,data_algG_CF_6_clean <- rbind(data_algG_CF_6_clean,data_algG_CF_6_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_6_filtered_third_variant_no_ref") == TRUE,data_algG_CF_6_clean_third_base_no_ref <- data_algG_CF_6_filtered_third_variant_no_ref[,which(names(data_algG_CF_6_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_6","Isolates_CF_6"))],z<-0)
ifelse(exists("data_algG_CF_6_filtered_third_variant_no_ref") == TRUE,data_algG_CF_6_clean <- rbind(data_algG_CF_6_clean,data_algG_CF_6_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_6_clean_20_filter <- merge(data_algG_CF_6_clean, data_algG_CF_6_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_6_clean_cleared <- data_algG_CF_6_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_6_clean_20_filter)){
  if(data_algG_CF_6_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_6_clean_20_filter$As_CF_6[i] >= 20){
    data_algG_CF_6_clean_cleared <- rbind(data_algG_CF_6_clean_cleared,data_algG_CF_6_clean_20_filter[i,])
  } else if(data_algG_CF_6_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_6_clean_20_filter$Cs_CF_6[i] >= 20){
    data_algG_CF_6_clean_cleared <- rbind(data_algG_CF_6_clean_cleared,data_algG_CF_6_clean_20_filter[i,])
  } else if(data_algG_CF_6_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_6_clean_20_filter$Gs_CF_6[i] >= 20){
    data_algG_CF_6_clean_cleared <- rbind(data_algG_CF_6_clean_cleared,data_algG_CF_6_clean_20_filter[i,])
  }else if(data_algG_CF_6_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_6_clean_20_filter$Ts_CF_6[i] >= 20){
    data_algG_CF_6_clean_cleared <- rbind(data_algG_CF_6_clean_cleared,data_algG_CF_6_clean_20_filter[i,])
  }
}

data_algG_CF_6_clean <- data_algG_CF_6_clean_cleared[,c(1:5)]



### CF_7

# Input 
data_algG_CF_7 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T7/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_7_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T7/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_7_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T7/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_7 <- data_algG_CF_7[!(data_algG_CF_7$Position %in% primer_correction_algG_CF_7_fwd$Position),]
data_algG_CF_7 <- data_algG_CF_7[!(data_algG_CF_7$Position %in% primer_correction_algG_CF_7_rev$Position),]
data_algG_CF_7 <- rbind(data_algG_CF_7,primer_correction_algG_CF_7_fwd,primer_correction_algG_CF_7_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_7 <- subset(data_algG_CF_7, coverage_CF_7 >199)
data_algG_CF_7_5UTR <- subset(data_algG_CF_7_5UTR, coverage_CF_7_5UTR >199)
data_algG_CF_7_3UTR <- subset(data_algG_CF_7_3UTR, coverage_CF_7_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_7_only_SNP <- subset(data_algG_CF_7, second_percentage <= 0.025) 
data_algG_CF_7_only_SNP$main_base_CF_7 <- colnames(data_algG_CF_7_only_SNP[, 8:11])[apply(data_algG_CF_7_only_SNP[, 8:11],1,which.max)]
data_algG_CF_7_only_SNP$main_base_CF_7 <- ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_As_CF_7","A", ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_Cs_CF_7","C", ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_Gs_CF_7","G", ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_Ts_CF_7","T",z<-0))))
data_algG_CF_7_only_SNP <- subset(data_algG_CF_7_only_SNP,GenBase_PA14 != main_base_CF_7)
ifelse(nrow(data_algG_CF_7_only_SNP) > 0,data_algG_CF_7_only_SNP$Isolates_CF_7 <- isolates_CF_7, z <- 0)


ifelse(nrow(data_algG_CF_7_only_SNP) > 0,data_algG_CF_7_only_SNP_clean <- data_algG_CF_7_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_7_only_SNP) > 0,names(data_algG_CF_7_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_7","Isolates_CF_7"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_7_filtered <- subset(data_algG_CF_7, second_percentage >= 0.025) 
data_algG_CF_7_filtered$main_base_CF_7 <- colnames(data_algG_CF_7_filtered[, 8:11])[apply(data_algG_CF_7_filtered[, 8:11],1,which.max)]
data_algG_CF_7_filtered$second_base_CF_7 <- colnames(data_algG_CF_7_filtered[, 8:11])[apply(data_algG_CF_7_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_7_filtered$third_base_CF_7 <- ifelse(apply(data_algG_CF_7_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_7_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_7_filtered$`PA14-Koordinate`)
  if(data_algG_CF_7_filtered$main_base_CF_7[a] == "percent_As_CF_7"){
    data_algG_CF_7_filtered$main_base_CF_7[a] <- "A"                     
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] == "percent_Cs_CF_7"){
    data_algG_CF_7_filtered$main_base_CF_7[a] <- "C"                     
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] == "percent_Gs_CF_7"){
    data_algG_CF_7_filtered$main_base_CF_7[a] <- "G"                     
  } else{
    data_algG_CF_7_filtered$main_base_CF_7[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_7_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_7_filtered$`PA14-Koordinate`)
  if(data_algG_CF_7_filtered$second_base_CF_7[a] == "percent_As_CF_7"){
    data_algG_CF_7_filtered$second_base_CF_7[a] <- "A"                     
  } else if(data_algG_CF_7_filtered$second_base_CF_7[a] == "percent_Cs_CF_7"){
    data_algG_CF_7_filtered$second_base_CF_7[a] <- "C"                     
  } else if(data_algG_CF_7_filtered$second_base_CF_7[a] == "percent_Gs_CF_7"){
    data_algG_CF_7_filtered$second_base_CF_7[a] <- "G"                     
  } else{
    data_algG_CF_7_filtered$second_base_CF_7[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_7_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_7_filtered$`PA14-Koordinate`)
  if(data_algG_CF_7_filtered$main_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$third_base_CF_7[a] == 0){
    data_algG_CF_7_filtered$third_base_CF_7[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_7_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_7_filtered$`PA14-Koordinate`)
  if(data_algG_CF_7_filtered$third_base_CF_7[a] == 1){
    data_algG_CF_7_filtered_third_variant <- data_algG_CF_7_filtered[a,]
    data_algG_CF_7_filtered$third_base_CF_7[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_7_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_7_filtered$`PA14-Koordinate`)
  if(data_algG_CF_7_filtered$third_base_CF_7[a] == 2){
    data_algG_CF_7_filtered_third_variant_no_ref <- data_algG_CF_7_filtered[a,]
    data_algG_CF_7_filtered$third_base_CF_7[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_7_filtered$third_base_CF_7) == 0, data_algG_CF_7_filtered <- subset(data_algG_CF_7_filtered, select = -c(third_base_CF_7)),data_algG_CF_7_filtered <- data_algG_CF_7_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_filtered_third_variant$third_percentage <- apply(data_algG_CF_7_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_7_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_filtered_third_variant$third_base_CF_7 <- colnames(data_algG_CF_7_filtered_third_variant[, 8:11])[apply(data_algG_CF_7_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_As_CF_7","A", ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_Cs_CF_7","C", ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_Gs_CF_7","G", ifelse(data_algG_CF_7_only_SNP$main_base_CF_7 == "percent_Ts_CF_7","T",z<-0))))

ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,
       data_algG_CF_7_filtered_third_variant$third_base_CF_7 <- 
         ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "percent_As_CF_7","A", 
                ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "percent_Cs_CF_7","C",
                       ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "percent_Gs_CF_7","G",
                              ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "percent_Ts_CF_7","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,
       data_algG_CF_7_filtered_third_variant$third_base_CF_7 <- ifelse(
         data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14,data_algG_CF_7_filtered_third_variant$third_base_CF_7 <- data_algG_CF_7_filtered_third_variant$second_base_CF_7,data_algG_CF_7_filtered_third_variant$third_base_CF_7
       ),
       z<-0
)

ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,
       data_algG_CF_7_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14, data_algG_CF_7_filtered_third_variant$second_percentage,data_algG_CF_7_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_filtered_third_variant <- data_algG_CF_7_filtered_third_variant[,-which(names(data_algG_CF_7_filtered_third_variant) %in% c("second_percentage","second_base_CF_7"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_7_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_7_filtered_third_variant_no_ref)[names(data_algG_CF_7_filtered_third_variant_no_ref) == "second_base_CF_7"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_7_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_7_filtered_third_variant_no_ref)[names(data_algG_CF_7_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_7",z<-0)




# Define SNP percentage & Add it

data_algG_CF_7_filtered$SNP_base <- 0
data_algG_CF_7_filtered$SNP_percentage_CF_7 <- 0


for(i in data_algG_CF_7_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_7_filtered$`PA14-Koordinate`)
  if(data_algG_CF_7_filtered$main_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] == "A"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_As_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] == "C"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Cs_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] == "G"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Gs_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] == "T"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Ts_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_7_filtered$second_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "A"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_As_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_7_filtered$second_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "C"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Cs_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_7_filtered$second_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "G"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Gs_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_7_filtered$second_base_CF_7[a] == data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "T"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Ts_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "A"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_As_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "C"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Cs_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "G"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Gs_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_7_filtered$main_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$second_base_CF_7[a] != data_algG_CF_7_filtered$GenBase_PA14[a] & data_algG_CF_7_filtered$main_base_CF_7[a] == "T"){
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- data_algG_CF_7_filtered$percent_Ts_CF_7[a]
    data_algG_CF_7_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_7_filtered$SNP_percentage_CF_7[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_filtered_third_variant$SNP_percentage_CF_7 <- 0,z<-0)


ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,
       data_algG_CF_7_filtered_third_variant$SNP_percentage_CF_7 <- 
         ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "A",
                data_algG_CF_7_filtered_third_variant$percent_As_CF_7, 
                ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "C",
                       data_algG_CF_7_filtered_third_variant$percent_Cs_CF_7,
                       ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "G",
                              data_algG_CF_7_filtered_third_variant$percent_Gs_CF_7,
                              ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "T",
                                     data_algG_CF_7_filtered_third_variant$percent_Ts_CF_7,
                                     ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "A",
                                            data_algG_CF_7_filtered_third_variant$percent_As_CF_7,
                                            ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "C",
                                                   data_algG_CF_7_filtered_third_variant$percent_Cs_CF_7,
                                                   ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "G",
                                                          data_algG_CF_7_filtered_third_variant$percent_Gs_CF_7,
                                                          ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "T",
                                                                 data_algG_CF_7_filtered_third_variant$percent_Ts_CF_7,
                                                                 ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "A",
                                                                        data_algG_CF_7_filtered_third_variant$percent_As_CF_7,
                                                                        ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "C",
                                                                               data_algG_CF_7_filtered_third_variant$percent_Cs_CF_7,
                                                                               ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "G",
                                                                                      data_algG_CF_7_filtered_third_variant$percent_Gs_CF_7,
                                                                                      ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "T",
                                                                                             data_algG_CF_7_filtered_third_variant$percent_Ts_CF_7,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,
       data_algG_CF_7_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "A",
                "A", 
                ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "C",
                       "C",
                       ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "G",
                              "G",
                              ifelse(data_algG_CF_7_filtered_third_variant$main_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "T",
                                     "T",
                                     ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "A",
                                            "A",
                                            ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 == data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$main_base_CF_7 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_7_filtered_third_variant$third_base_CF_7 != data_algG_CF_7_filtered_third_variant$GenBase_PA14 & data_algG_CF_7_filtered_third_variant$third_base_CF_7 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_7_filtered$Isolates_CF_7 <- round((data_algG_CF_7_filtered$SNP_percentage_CF_7*100)/(100/isolates_CF_7))
ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_filtered_third_variant$Isolates_CF_7 <- round((data_algG_CF_7_filtered_third_variant$SNP_percentage_CF_7*100)/(100/isolates_CF_7)),z<-0)
ifelse(exists("data_algG_CF_7_filtered_third_variant_no_ref") == TRUE,data_algG_CF_7_filtered_third_variant_no_ref$Isolates_CF_7 <- round((data_algG_CF_7_filtered_third_variant_no_ref$SNP_percentage_CF_7*100)/(100/isolates_CF_7)),z<-0)



# Clean data

data_algG_CF_7_clean <- data_algG_CF_7_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_clean_third_base <- data_algG_CF_7_filtered_third_variant[,which(names(data_algG_CF_7_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_7","Isolates_CF_7"))],z<-0)
ifelse(exists("data_algG_CF_7_filtered_third_variant") == TRUE,data_algG_CF_7_clean <- rbind(data_algG_CF_7_clean,data_algG_CF_7_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_7_filtered_third_variant_no_ref") == TRUE,data_algG_CF_7_clean_third_base_no_ref <- data_algG_CF_7_filtered_third_variant_no_ref[,which(names(data_algG_CF_7_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_7","Isolates_CF_7"))],z<-0)
ifelse(exists("data_algG_CF_7_filtered_third_variant_no_ref") == TRUE,data_algG_CF_7_clean <- rbind(data_algG_CF_7_clean,data_algG_CF_7_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_7_clean_20_filter <- merge(data_algG_CF_7_clean, data_algG_CF_7_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_7_clean_cleared <- data_algG_CF_7_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_7_clean_20_filter)){
  if(data_algG_CF_7_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_7_clean_20_filter$As_CF_7[i] >= 20){
    data_algG_CF_7_clean_cleared <- rbind(data_algG_CF_7_clean_cleared,data_algG_CF_7_clean_20_filter[i,])
  } else if(data_algG_CF_7_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_7_clean_20_filter$Cs_CF_7[i] >= 20){
    data_algG_CF_7_clean_cleared <- rbind(data_algG_CF_7_clean_cleared,data_algG_CF_7_clean_20_filter[i,])
  } else if(data_algG_CF_7_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_7_clean_20_filter$Gs_CF_7[i] >= 20){
    data_algG_CF_7_clean_cleared <- rbind(data_algG_CF_7_clean_cleared,data_algG_CF_7_clean_20_filter[i,])
  }else if(data_algG_CF_7_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_7_clean_20_filter$Ts_CF_7[i] >= 20){
    data_algG_CF_7_clean_cleared <- rbind(data_algG_CF_7_clean_cleared,data_algG_CF_7_clean_20_filter[i,])
  }
}

data_algG_CF_7_clean <- data_algG_CF_7_clean_cleared[,c(1:5)]


### CF_8

# Input 
data_algG_CF_8 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T8/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_8_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T8/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_8_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T8/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_8 <- data_algG_CF_8[!(data_algG_CF_8$Position %in% primer_correction_algG_CF_8_fwd$Position),]
data_algG_CF_8 <- data_algG_CF_8[!(data_algG_CF_8$Position %in% primer_correction_algG_CF_8_rev$Position),]
data_algG_CF_8 <- rbind(data_algG_CF_8,primer_correction_algG_CF_8_fwd,primer_correction_algG_CF_8_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_8 <- subset(data_algG_CF_8, coverage_CF_8 >199)
data_algG_CF_8_5UTR <- subset(data_algG_CF_8_5UTR, coverage_CF_8_5UTR >199)
data_algG_CF_8_3UTR <- subset(data_algG_CF_8_3UTR, coverage_CF_8_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_8_only_SNP <- subset(data_algG_CF_8, second_percentage <= 0.025) 
data_algG_CF_8_only_SNP$main_base_CF_8 <- colnames(data_algG_CF_8_only_SNP[, 8:11])[apply(data_algG_CF_8_only_SNP[, 8:11],1,which.max)]
data_algG_CF_8_only_SNP$main_base_CF_8 <- ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_As_CF_8","A", ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_Cs_CF_8","C", ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_Gs_CF_8","G", ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_Ts_CF_8","T",z<-0))))
data_algG_CF_8_only_SNP <- subset(data_algG_CF_8_only_SNP,GenBase_PA14 != main_base_CF_8)
ifelse(nrow(data_algG_CF_8_only_SNP) > 0,data_algG_CF_8_only_SNP$Isolates_CF_8 <- isolates_CF_8, z <- 0)


ifelse(nrow(data_algG_CF_8_only_SNP) > 0,data_algG_CF_8_only_SNP_clean <- data_algG_CF_8_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_8_only_SNP) > 0,names(data_algG_CF_8_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_8","Isolates_CF_8"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_8_filtered <- subset(data_algG_CF_8, second_percentage >= 0.025) 
data_algG_CF_8_filtered$main_base_CF_8 <- colnames(data_algG_CF_8_filtered[, 8:11])[apply(data_algG_CF_8_filtered[, 8:11],1,which.max)]
data_algG_CF_8_filtered$second_base_CF_8 <- colnames(data_algG_CF_8_filtered[, 8:11])[apply(data_algG_CF_8_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_8_filtered$third_base_CF_8 <- ifelse(apply(data_algG_CF_8_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_8_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_8_filtered$`PA14-Koordinate`)
  if(data_algG_CF_8_filtered$main_base_CF_8[a] == "percent_As_CF_8"){
    data_algG_CF_8_filtered$main_base_CF_8[a] <- "A"                     
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] == "percent_Cs_CF_8"){
    data_algG_CF_8_filtered$main_base_CF_8[a] <- "C"                     
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] == "percent_Gs_CF_8"){
    data_algG_CF_8_filtered$main_base_CF_8[a] <- "G"                     
  } else{
    data_algG_CF_8_filtered$main_base_CF_8[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_8_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_8_filtered$`PA14-Koordinate`)
  if(data_algG_CF_8_filtered$second_base_CF_8[a] == "percent_As_CF_8"){
    data_algG_CF_8_filtered$second_base_CF_8[a] <- "A"                     
  } else if(data_algG_CF_8_filtered$second_base_CF_8[a] == "percent_Cs_CF_8"){
    data_algG_CF_8_filtered$second_base_CF_8[a] <- "C"                     
  } else if(data_algG_CF_8_filtered$second_base_CF_8[a] == "percent_Gs_CF_8"){
    data_algG_CF_8_filtered$second_base_CF_8[a] <- "G"                     
  } else{
    data_algG_CF_8_filtered$second_base_CF_8[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_8_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_8_filtered$`PA14-Koordinate`)
  if(data_algG_CF_8_filtered$main_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$third_base_CF_8[a] == 0){
    data_algG_CF_8_filtered$third_base_CF_8[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_8_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_8_filtered$`PA14-Koordinate`)
  if(data_algG_CF_8_filtered$third_base_CF_8[a] == 1){
    data_algG_CF_8_filtered_third_variant <- data_algG_CF_8_filtered[a,]
    data_algG_CF_8_filtered$third_base_CF_8[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_8_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_8_filtered$`PA14-Koordinate`)
  if(data_algG_CF_8_filtered$third_base_CF_8[a] == 2){
    data_algG_CF_8_filtered_third_variant_no_ref <- data_algG_CF_8_filtered[a,]
    data_algG_CF_8_filtered$third_base_CF_8[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_8_filtered$third_base_CF_8) == 0, data_algG_CF_8_filtered <- subset(data_algG_CF_8_filtered, select = -c(third_base_CF_8)),data_algG_CF_8_filtered <- data_algG_CF_8_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_filtered_third_variant$third_percentage <- apply(data_algG_CF_8_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_8_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_filtered_third_variant$third_base_CF_8 <- colnames(data_algG_CF_8_filtered_third_variant[, 8:11])[apply(data_algG_CF_8_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_As_CF_8","A", ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_Cs_CF_8","C", ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_Gs_CF_8","G", ifelse(data_algG_CF_8_only_SNP$main_base_CF_8 == "percent_Ts_CF_8","T",z<-0))))

ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,
       data_algG_CF_8_filtered_third_variant$third_base_CF_8 <- 
         ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "percent_As_CF_8","A", 
                ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "percent_Cs_CF_8","C",
                       ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "percent_Gs_CF_8","G",
                              ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "percent_Ts_CF_8","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,
       data_algG_CF_8_filtered_third_variant$third_base_CF_8 <- ifelse(
         data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14,data_algG_CF_8_filtered_third_variant$third_base_CF_8 <- data_algG_CF_8_filtered_third_variant$second_base_CF_8,data_algG_CF_8_filtered_third_variant$third_base_CF_8
       ),
       z<-0
)

ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,
       data_algG_CF_8_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14, data_algG_CF_8_filtered_third_variant$second_percentage,data_algG_CF_8_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_filtered_third_variant <- data_algG_CF_8_filtered_third_variant[,-which(names(data_algG_CF_8_filtered_third_variant) %in% c("second_percentage","second_base_CF_8"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_8_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_8_filtered_third_variant_no_ref)[names(data_algG_CF_8_filtered_third_variant_no_ref) == "second_base_CF_8"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_8_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_8_filtered_third_variant_no_ref)[names(data_algG_CF_8_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_8",z<-0)




# Define SNP percentage & Add it

data_algG_CF_8_filtered$SNP_base <- 0
data_algG_CF_8_filtered$SNP_percentage_CF_8 <- 0


for(i in data_algG_CF_8_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_8_filtered$`PA14-Koordinate`)
  if(data_algG_CF_8_filtered$main_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] == "A"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_As_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] == "C"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Cs_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] == "G"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Gs_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] == "T"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Ts_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_8_filtered$second_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "A"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_As_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_8_filtered$second_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "C"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Cs_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_8_filtered$second_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "G"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Gs_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_8_filtered$second_base_CF_8[a] == data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "T"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Ts_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "A"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_As_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "C"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Cs_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "G"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Gs_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_8_filtered$main_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$second_base_CF_8[a] != data_algG_CF_8_filtered$GenBase_PA14[a] & data_algG_CF_8_filtered$main_base_CF_8[a] == "T"){
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- data_algG_CF_8_filtered$percent_Ts_CF_8[a]
    data_algG_CF_8_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_8_filtered$SNP_percentage_CF_8[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_filtered_third_variant$SNP_percentage_CF_8 <- 0,z<-0)


ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,
       data_algG_CF_8_filtered_third_variant$SNP_percentage_CF_8 <- 
         ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "A",
                data_algG_CF_8_filtered_third_variant$percent_As_CF_8, 
                ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "C",
                       data_algG_CF_8_filtered_third_variant$percent_Cs_CF_8,
                       ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "G",
                              data_algG_CF_8_filtered_third_variant$percent_Gs_CF_8,
                              ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "T",
                                     data_algG_CF_8_filtered_third_variant$percent_Ts_CF_8,
                                     ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "A",
                                            data_algG_CF_8_filtered_third_variant$percent_As_CF_8,
                                            ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "C",
                                                   data_algG_CF_8_filtered_third_variant$percent_Cs_CF_8,
                                                   ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "G",
                                                          data_algG_CF_8_filtered_third_variant$percent_Gs_CF_8,
                                                          ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "T",
                                                                 data_algG_CF_8_filtered_third_variant$percent_Ts_CF_8,
                                                                 ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "A",
                                                                        data_algG_CF_8_filtered_third_variant$percent_As_CF_8,
                                                                        ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "C",
                                                                               data_algG_CF_8_filtered_third_variant$percent_Cs_CF_8,
                                                                               ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "G",
                                                                                      data_algG_CF_8_filtered_third_variant$percent_Gs_CF_8,
                                                                                      ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "T",
                                                                                             data_algG_CF_8_filtered_third_variant$percent_Ts_CF_8,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,
       data_algG_CF_8_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "A",
                "A", 
                ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "C",
                       "C",
                       ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "G",
                              "G",
                              ifelse(data_algG_CF_8_filtered_third_variant$main_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "T",
                                     "T",
                                     ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "A",
                                            "A",
                                            ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 == data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$main_base_CF_8 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_8_filtered_third_variant$third_base_CF_8 != data_algG_CF_8_filtered_third_variant$GenBase_PA14 & data_algG_CF_8_filtered_third_variant$third_base_CF_8 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_8_filtered$Isolates_CF_8 <- round((data_algG_CF_8_filtered$SNP_percentage_CF_8*100)/(100/isolates_CF_8))
ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_filtered_third_variant$Isolates_CF_8 <- round((data_algG_CF_8_filtered_third_variant$SNP_percentage_CF_8*100)/(100/isolates_CF_8)),z<-0)
ifelse(exists("data_algG_CF_8_filtered_third_variant_no_ref") == TRUE,data_algG_CF_8_filtered_third_variant_no_ref$Isolates_CF_8 <- round((data_algG_CF_8_filtered_third_variant_no_ref$SNP_percentage_CF_8*100)/(100/isolates_CF_8)),z<-0)



# Clean data

data_algG_CF_8_clean <- data_algG_CF_8_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_clean_third_base <- data_algG_CF_8_filtered_third_variant[,which(names(data_algG_CF_8_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_8","Isolates_CF_8"))],z<-0)
ifelse(exists("data_algG_CF_8_filtered_third_variant") == TRUE,data_algG_CF_8_clean <- rbind(data_algG_CF_8_clean,data_algG_CF_8_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_8_filtered_third_variant_no_ref") == TRUE,data_algG_CF_8_clean_third_base_no_ref <- data_algG_CF_8_filtered_third_variant_no_ref[,which(names(data_algG_CF_8_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_8","Isolates_CF_8"))],z<-0)
ifelse(exists("data_algG_CF_8_filtered_third_variant_no_ref") == TRUE,data_algG_CF_8_clean <- rbind(data_algG_CF_8_clean,data_algG_CF_8_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_8_clean_20_filter <- merge(data_algG_CF_8_clean, data_algG_CF_8_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_8_clean_cleared <- data_algG_CF_8_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_8_clean_20_filter)){
  if(data_algG_CF_8_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_8_clean_20_filter$As_CF_8[i] >= 20){
    data_algG_CF_8_clean_cleared <- rbind(data_algG_CF_8_clean_cleared,data_algG_CF_8_clean_20_filter[i,])
  } else if(data_algG_CF_8_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_8_clean_20_filter$Cs_CF_8[i] >= 20){
    data_algG_CF_8_clean_cleared <- rbind(data_algG_CF_8_clean_cleared,data_algG_CF_8_clean_20_filter[i,])
  } else if(data_algG_CF_8_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_8_clean_20_filter$Gs_CF_8[i] >= 20){
    data_algG_CF_8_clean_cleared <- rbind(data_algG_CF_8_clean_cleared,data_algG_CF_8_clean_20_filter[i,])
  }else if(data_algG_CF_8_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_8_clean_20_filter$Ts_CF_8[i] >= 20){
    data_algG_CF_8_clean_cleared <- rbind(data_algG_CF_8_clean_cleared,data_algG_CF_8_clean_20_filter[i,])
  }
}

data_algG_CF_8_clean <- data_algG_CF_8_clean_cleared[,c(1:5)]


### CF_9

# Input 
data_algG_CF_9 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T9/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_9_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T9/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_9_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T9/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_9 <- data_algG_CF_9[!(data_algG_CF_9$Position %in% primer_correction_algG_CF_9_fwd$Position),]
data_algG_CF_9 <- data_algG_CF_9[!(data_algG_CF_9$Position %in% primer_correction_algG_CF_9_rev$Position),]
data_algG_CF_9 <- rbind(data_algG_CF_9,primer_correction_algG_CF_9_fwd,primer_correction_algG_CF_9_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_9 <- subset(data_algG_CF_9, coverage_CF_9 >199)
data_algG_CF_9_5UTR <- subset(data_algG_CF_9_5UTR, coverage_CF_9_5UTR >199)
data_algG_CF_9_3UTR <- subset(data_algG_CF_9_3UTR, coverage_CF_9_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_9_only_SNP <- subset(data_algG_CF_9, second_percentage <= 0.025) 
data_algG_CF_9_only_SNP$main_base_CF_9 <- colnames(data_algG_CF_9_only_SNP[, 8:11])[apply(data_algG_CF_9_only_SNP[, 8:11],1,which.max)]
data_algG_CF_9_only_SNP$main_base_CF_9 <- ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_As_CF_9","A", ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_Cs_CF_9","C", ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_Gs_CF_9","G", ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_Ts_CF_9","T",z<-0))))
data_algG_CF_9_only_SNP <- subset(data_algG_CF_9_only_SNP,GenBase_PA14 != main_base_CF_9)
ifelse(nrow(data_algG_CF_9_only_SNP) > 0,data_algG_CF_9_only_SNP$Isolates_CF_9 <- isolates_CF_9, z <- 0)


ifelse(nrow(data_algG_CF_9_only_SNP) > 0,data_algG_CF_9_only_SNP_clean <- data_algG_CF_9_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_9_only_SNP) > 0,names(data_algG_CF_9_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_9","Isolates_CF_9"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_9_filtered <- subset(data_algG_CF_9, second_percentage >= 0.025) 
data_algG_CF_9_filtered$main_base_CF_9 <- colnames(data_algG_CF_9_filtered[, 8:11])[apply(data_algG_CF_9_filtered[, 8:11],1,which.max)]
data_algG_CF_9_filtered$second_base_CF_9 <- colnames(data_algG_CF_9_filtered[, 8:11])[apply(data_algG_CF_9_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_9_filtered$third_base_CF_9 <- ifelse(apply(data_algG_CF_9_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_9_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_9_filtered$`PA14-Koordinate`)
  if(data_algG_CF_9_filtered$main_base_CF_9[a] == "percent_As_CF_9"){
    data_algG_CF_9_filtered$main_base_CF_9[a] <- "A"                     
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] == "percent_Cs_CF_9"){
    data_algG_CF_9_filtered$main_base_CF_9[a] <- "C"                     
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] == "percent_Gs_CF_9"){
    data_algG_CF_9_filtered$main_base_CF_9[a] <- "G"                     
  } else{
    data_algG_CF_9_filtered$main_base_CF_9[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_9_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_9_filtered$`PA14-Koordinate`)
  if(data_algG_CF_9_filtered$second_base_CF_9[a] == "percent_As_CF_9"){
    data_algG_CF_9_filtered$second_base_CF_9[a] <- "A"                     
  } else if(data_algG_CF_9_filtered$second_base_CF_9[a] == "percent_Cs_CF_9"){
    data_algG_CF_9_filtered$second_base_CF_9[a] <- "C"                     
  } else if(data_algG_CF_9_filtered$second_base_CF_9[a] == "percent_Gs_CF_9"){
    data_algG_CF_9_filtered$second_base_CF_9[a] <- "G"                     
  } else{
    data_algG_CF_9_filtered$second_base_CF_9[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_9_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_9_filtered$`PA14-Koordinate`)
  if(data_algG_CF_9_filtered$main_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$third_base_CF_9[a] == 0){
    data_algG_CF_9_filtered$third_base_CF_9[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_9_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_9_filtered$`PA14-Koordinate`)
  if(data_algG_CF_9_filtered$third_base_CF_9[a] == 1){
    data_algG_CF_9_filtered_third_variant <- data_algG_CF_9_filtered[a,]
    data_algG_CF_9_filtered$third_base_CF_9[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_9_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_9_filtered$`PA14-Koordinate`)
  if(data_algG_CF_9_filtered$third_base_CF_9[a] == 2){
    data_algG_CF_9_filtered_third_variant_no_ref <- data_algG_CF_9_filtered[a,]
    data_algG_CF_9_filtered$third_base_CF_9[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_9_filtered$third_base_CF_9) == 0, data_algG_CF_9_filtered <- subset(data_algG_CF_9_filtered, select = -c(third_base_CF_9)),data_algG_CF_9_filtered <- data_algG_CF_9_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_filtered_third_variant$third_percentage <- apply(data_algG_CF_9_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_9_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_filtered_third_variant$third_base_CF_9 <- colnames(data_algG_CF_9_filtered_third_variant[, 8:11])[apply(data_algG_CF_9_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_As_CF_9","A", ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_Cs_CF_9","C", ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_Gs_CF_9","G", ifelse(data_algG_CF_9_only_SNP$main_base_CF_9 == "percent_Ts_CF_9","T",z<-0))))

ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,
       data_algG_CF_9_filtered_third_variant$third_base_CF_9 <- 
         ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "percent_As_CF_9","A", 
                ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "percent_Cs_CF_9","C",
                       ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "percent_Gs_CF_9","G",
                              ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "percent_Ts_CF_9","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,
       data_algG_CF_9_filtered_third_variant$third_base_CF_9 <- ifelse(
         data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14,data_algG_CF_9_filtered_third_variant$third_base_CF_9 <- data_algG_CF_9_filtered_third_variant$second_base_CF_9,data_algG_CF_9_filtered_third_variant$third_base_CF_9
       ),
       z<-0
)

ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,
       data_algG_CF_9_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14, data_algG_CF_9_filtered_third_variant$second_percentage,data_algG_CF_9_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_filtered_third_variant <- data_algG_CF_9_filtered_third_variant[,-which(names(data_algG_CF_9_filtered_third_variant) %in% c("second_percentage","second_base_CF_9"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_9_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_9_filtered_third_variant_no_ref)[names(data_algG_CF_9_filtered_third_variant_no_ref) == "second_base_CF_9"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_9_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_9_filtered_third_variant_no_ref)[names(data_algG_CF_9_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_9",z<-0)




# Define SNP percentage & Add it

data_algG_CF_9_filtered$SNP_base <- 0
data_algG_CF_9_filtered$SNP_percentage_CF_9 <- 0


for(i in data_algG_CF_9_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_9_filtered$`PA14-Koordinate`)
  if(data_algG_CF_9_filtered$main_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] == "A"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_As_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] == "C"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Cs_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] == "G"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Gs_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] == "T"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Ts_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_9_filtered$second_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "A"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_As_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_9_filtered$second_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "C"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Cs_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_9_filtered$second_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "G"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Gs_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_9_filtered$second_base_CF_9[a] == data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "T"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Ts_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "A"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_As_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "C"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Cs_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "G"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Gs_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_9_filtered$main_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$second_base_CF_9[a] != data_algG_CF_9_filtered$GenBase_PA14[a] & data_algG_CF_9_filtered$main_base_CF_9[a] == "T"){
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- data_algG_CF_9_filtered$percent_Ts_CF_9[a]
    data_algG_CF_9_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_9_filtered$SNP_percentage_CF_9[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_filtered_third_variant$SNP_percentage_CF_9 <- 0,z<-0)


ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,
       data_algG_CF_9_filtered_third_variant$SNP_percentage_CF_9 <- 
         ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "A",
                data_algG_CF_9_filtered_third_variant$percent_As_CF_9, 
                ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "C",
                       data_algG_CF_9_filtered_third_variant$percent_Cs_CF_9,
                       ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "G",
                              data_algG_CF_9_filtered_third_variant$percent_Gs_CF_9,
                              ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "T",
                                     data_algG_CF_9_filtered_third_variant$percent_Ts_CF_9,
                                     ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "A",
                                            data_algG_CF_9_filtered_third_variant$percent_As_CF_9,
                                            ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "C",
                                                   data_algG_CF_9_filtered_third_variant$percent_Cs_CF_9,
                                                   ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "G",
                                                          data_algG_CF_9_filtered_third_variant$percent_Gs_CF_9,
                                                          ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "T",
                                                                 data_algG_CF_9_filtered_third_variant$percent_Ts_CF_9,
                                                                 ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "A",
                                                                        data_algG_CF_9_filtered_third_variant$percent_As_CF_9,
                                                                        ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "C",
                                                                               data_algG_CF_9_filtered_third_variant$percent_Cs_CF_9,
                                                                               ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "G",
                                                                                      data_algG_CF_9_filtered_third_variant$percent_Gs_CF_9,
                                                                                      ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "T",
                                                                                             data_algG_CF_9_filtered_third_variant$percent_Ts_CF_9,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,
       data_algG_CF_9_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "A",
                "A", 
                ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "C",
                       "C",
                       ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "G",
                              "G",
                              ifelse(data_algG_CF_9_filtered_third_variant$main_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "T",
                                     "T",
                                     ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "A",
                                            "A",
                                            ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 == data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$main_base_CF_9 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_9_filtered_third_variant$third_base_CF_9 != data_algG_CF_9_filtered_third_variant$GenBase_PA14 & data_algG_CF_9_filtered_third_variant$third_base_CF_9 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_9_filtered$Isolates_CF_9 <- round((data_algG_CF_9_filtered$SNP_percentage_CF_9*100)/(100/isolates_CF_9))
ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_filtered_third_variant$Isolates_CF_9 <- round((data_algG_CF_9_filtered_third_variant$SNP_percentage_CF_9*100)/(100/isolates_CF_9)),z<-0)
ifelse(exists("data_algG_CF_9_filtered_third_variant_no_ref") == TRUE,data_algG_CF_9_filtered_third_variant_no_ref$Isolates_CF_9 <- round((data_algG_CF_9_filtered_third_variant_no_ref$SNP_percentage_CF_9*100)/(100/isolates_CF_9)),z<-0)



# Clean data

data_algG_CF_9_clean <- data_algG_CF_9_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_clean_third_base <- data_algG_CF_9_filtered_third_variant[,which(names(data_algG_CF_9_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_9","Isolates_CF_9"))],z<-0)
ifelse(exists("data_algG_CF_9_filtered_third_variant") == TRUE,data_algG_CF_9_clean <- rbind(data_algG_CF_9_clean,data_algG_CF_9_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_9_filtered_third_variant_no_ref") == TRUE,data_algG_CF_9_clean_third_base_no_ref <- data_algG_CF_9_filtered_third_variant_no_ref[,which(names(data_algG_CF_9_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_9","Isolates_CF_9"))],z<-0)
ifelse(exists("data_algG_CF_9_filtered_third_variant_no_ref") == TRUE,data_algG_CF_9_clean <- rbind(data_algG_CF_9_clean,data_algG_CF_9_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_9_clean_20_filter <- merge(data_algG_CF_9_clean, data_algG_CF_9_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_9_clean_cleared <- data_algG_CF_9_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_9_clean_20_filter)){
  if(data_algG_CF_9_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_9_clean_20_filter$As_CF_9[i] >= 20){
    data_algG_CF_9_clean_cleared <- rbind(data_algG_CF_9_clean_cleared,data_algG_CF_9_clean_20_filter[i,])
  } else if(data_algG_CF_9_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_9_clean_20_filter$Cs_CF_9[i] >= 20){
    data_algG_CF_9_clean_cleared <- rbind(data_algG_CF_9_clean_cleared,data_algG_CF_9_clean_20_filter[i,])
  } else if(data_algG_CF_9_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_9_clean_20_filter$Gs_CF_9[i] >= 20){
    data_algG_CF_9_clean_cleared <- rbind(data_algG_CF_9_clean_cleared,data_algG_CF_9_clean_20_filter[i,])
  }else if(data_algG_CF_9_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_9_clean_20_filter$Ts_CF_9[i] >= 20){
    data_algG_CF_9_clean_cleared <- rbind(data_algG_CF_9_clean_cleared,data_algG_CF_9_clean_20_filter[i,])
  }
}

data_algG_CF_9_clean <- data_algG_CF_9_clean_cleared[,c(1:5)]


### CF_10

# Input 
data_algG_CF_10 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T10/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_10_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T10/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_10_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T10/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_10 <- data_algG_CF_10[!(data_algG_CF_10$Position %in% primer_correction_algG_CF_10_fwd$Position),]
data_algG_CF_10 <- data_algG_CF_10[!(data_algG_CF_10$Position %in% primer_correction_algG_CF_10_rev$Position),]
data_algG_CF_10 <- rbind(data_algG_CF_10,primer_correction_algG_CF_10_fwd,primer_correction_algG_CF_10_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_10 <- subset(data_algG_CF_10, coverage_CF_10 >199)
data_algG_CF_10_5UTR <- subset(data_algG_CF_10_5UTR, coverage_CF_10_5UTR >199)
data_algG_CF_10_3UTR <- subset(data_algG_CF_10_3UTR, coverage_CF_10_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_10_only_SNP <- subset(data_algG_CF_10, second_percentage <= 0.025) 
data_algG_CF_10_only_SNP$main_base_CF_10 <- colnames(data_algG_CF_10_only_SNP[, 8:11])[apply(data_algG_CF_10_only_SNP[, 8:11],1,which.max)]
data_algG_CF_10_only_SNP$main_base_CF_10 <- ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_As_CF_10","A", ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_Cs_CF_10","C", ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_Gs_CF_10","G", ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_Ts_CF_10","T",z<-0))))
data_algG_CF_10_only_SNP <- subset(data_algG_CF_10_only_SNP,GenBase_PA14 != main_base_CF_10)
ifelse(nrow(data_algG_CF_10_only_SNP) > 0,data_algG_CF_10_only_SNP$Isolates_CF_10 <- isolates_CF_10, z <- 0)


ifelse(nrow(data_algG_CF_10_only_SNP) > 0,data_algG_CF_10_only_SNP_clean <- data_algG_CF_10_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_10_only_SNP) > 0,names(data_algG_CF_10_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_10","Isolates_CF_10"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_10_filtered <- subset(data_algG_CF_10, second_percentage >= 0.025) 
data_algG_CF_10_filtered$main_base_CF_10 <- colnames(data_algG_CF_10_filtered[, 8:11])[apply(data_algG_CF_10_filtered[, 8:11],1,which.max)]
data_algG_CF_10_filtered$second_base_CF_10 <- colnames(data_algG_CF_10_filtered[, 8:11])[apply(data_algG_CF_10_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_10_filtered$third_base_CF_10 <- ifelse(apply(data_algG_CF_10_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_10_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_10_filtered$`PA14-Koordinate`)
  if(data_algG_CF_10_filtered$main_base_CF_10[a] == "percent_As_CF_10"){
    data_algG_CF_10_filtered$main_base_CF_10[a] <- "A"                     
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] == "percent_Cs_CF_10"){
    data_algG_CF_10_filtered$main_base_CF_10[a] <- "C"                     
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] == "percent_Gs_CF_10"){
    data_algG_CF_10_filtered$main_base_CF_10[a] <- "G"                     
  } else{
    data_algG_CF_10_filtered$main_base_CF_10[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_10_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_10_filtered$`PA14-Koordinate`)
  if(data_algG_CF_10_filtered$second_base_CF_10[a] == "percent_As_CF_10"){
    data_algG_CF_10_filtered$second_base_CF_10[a] <- "A"                     
  } else if(data_algG_CF_10_filtered$second_base_CF_10[a] == "percent_Cs_CF_10"){
    data_algG_CF_10_filtered$second_base_CF_10[a] <- "C"                     
  } else if(data_algG_CF_10_filtered$second_base_CF_10[a] == "percent_Gs_CF_10"){
    data_algG_CF_10_filtered$second_base_CF_10[a] <- "G"                     
  } else{
    data_algG_CF_10_filtered$second_base_CF_10[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_10_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_10_filtered$`PA14-Koordinate`)
  if(data_algG_CF_10_filtered$main_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$third_base_CF_10[a] == 0){
    data_algG_CF_10_filtered$third_base_CF_10[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_10_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_10_filtered$`PA14-Koordinate`)
  if(data_algG_CF_10_filtered$third_base_CF_10[a] == 1){
    data_algG_CF_10_filtered_third_variant <- data_algG_CF_10_filtered[a,]
    data_algG_CF_10_filtered$third_base_CF_10[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_10_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_10_filtered$`PA14-Koordinate`)
  if(data_algG_CF_10_filtered$third_base_CF_10[a] == 2){
    data_algG_CF_10_filtered_third_variant_no_ref <- data_algG_CF_10_filtered[a,]
    data_algG_CF_10_filtered$third_base_CF_10[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_10_filtered$third_base_CF_10) == 0, data_algG_CF_10_filtered <- subset(data_algG_CF_10_filtered, select = -c(third_base_CF_10)),data_algG_CF_10_filtered <- data_algG_CF_10_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_filtered_third_variant$third_percentage <- apply(data_algG_CF_10_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_10_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_filtered_third_variant$third_base_CF_10 <- colnames(data_algG_CF_10_filtered_third_variant[, 8:11])[apply(data_algG_CF_10_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_As_CF_10","A", ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_Cs_CF_10","C", ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_Gs_CF_10","G", ifelse(data_algG_CF_10_only_SNP$main_base_CF_10 == "percent_Ts_CF_10","T",z<-0))))

ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,
       data_algG_CF_10_filtered_third_variant$third_base_CF_10 <- 
         ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "percent_As_CF_10","A", 
                ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "percent_Cs_CF_10","C",
                       ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "percent_Gs_CF_10","G",
                              ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "percent_Ts_CF_10","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,
       data_algG_CF_10_filtered_third_variant$third_base_CF_10 <- ifelse(
         data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14,data_algG_CF_10_filtered_third_variant$third_base_CF_10 <- data_algG_CF_10_filtered_third_variant$second_base_CF_10,data_algG_CF_10_filtered_third_variant$third_base_CF_10
       ),
       z<-0
)

ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,
       data_algG_CF_10_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14, data_algG_CF_10_filtered_third_variant$second_percentage,data_algG_CF_10_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_filtered_third_variant <- data_algG_CF_10_filtered_third_variant[,-which(names(data_algG_CF_10_filtered_third_variant) %in% c("second_percentage","second_base_CF_10"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_10_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_10_filtered_third_variant_no_ref)[names(data_algG_CF_10_filtered_third_variant_no_ref) == "second_base_CF_10"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_10_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_10_filtered_third_variant_no_ref)[names(data_algG_CF_10_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_10",z<-0)




# Define SNP percentage & Add it

data_algG_CF_10_filtered$SNP_base <- 0
data_algG_CF_10_filtered$SNP_percentage_CF_10 <- 0


for(i in data_algG_CF_10_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_10_filtered$`PA14-Koordinate`)
  if(data_algG_CF_10_filtered$main_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] == "A"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_As_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] == "C"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Cs_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] == "G"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Gs_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] == "T"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Ts_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_10_filtered$second_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "A"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_As_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_10_filtered$second_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "C"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Cs_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_10_filtered$second_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "G"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Gs_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_10_filtered$second_base_CF_10[a] == data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "T"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Ts_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "A"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_As_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "C"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Cs_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "G"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Gs_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_10_filtered$main_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$second_base_CF_10[a] != data_algG_CF_10_filtered$GenBase_PA14[a] & data_algG_CF_10_filtered$main_base_CF_10[a] == "T"){
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- data_algG_CF_10_filtered$percent_Ts_CF_10[a]
    data_algG_CF_10_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_10_filtered$SNP_percentage_CF_10[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_filtered_third_variant$SNP_percentage_CF_10 <- 0,z<-0)


ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,
       data_algG_CF_10_filtered_third_variant$SNP_percentage_CF_10 <- 
         ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "A",
                data_algG_CF_10_filtered_third_variant$percent_As_CF_10, 
                ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "C",
                       data_algG_CF_10_filtered_third_variant$percent_Cs_CF_10,
                       ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "G",
                              data_algG_CF_10_filtered_third_variant$percent_Gs_CF_10,
                              ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "T",
                                     data_algG_CF_10_filtered_third_variant$percent_Ts_CF_10,
                                     ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "A",
                                            data_algG_CF_10_filtered_third_variant$percent_As_CF_10,
                                            ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "C",
                                                   data_algG_CF_10_filtered_third_variant$percent_Cs_CF_10,
                                                   ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "G",
                                                          data_algG_CF_10_filtered_third_variant$percent_Gs_CF_10,
                                                          ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "T",
                                                                 data_algG_CF_10_filtered_third_variant$percent_Ts_CF_10,
                                                                 ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "A",
                                                                        data_algG_CF_10_filtered_third_variant$percent_As_CF_10,
                                                                        ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "C",
                                                                               data_algG_CF_10_filtered_third_variant$percent_Cs_CF_10,
                                                                               ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "G",
                                                                                      data_algG_CF_10_filtered_third_variant$percent_Gs_CF_10,
                                                                                      ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "T",
                                                                                             data_algG_CF_10_filtered_third_variant$percent_Ts_CF_10,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,
       data_algG_CF_10_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "A",
                "A", 
                ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "C",
                       "C",
                       ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "G",
                              "G",
                              ifelse(data_algG_CF_10_filtered_third_variant$main_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "T",
                                     "T",
                                     ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "A",
                                            "A",
                                            ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 == data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$main_base_CF_10 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_10_filtered_third_variant$third_base_CF_10 != data_algG_CF_10_filtered_third_variant$GenBase_PA14 & data_algG_CF_10_filtered_third_variant$third_base_CF_10 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_10_filtered$Isolates_CF_10 <- round((data_algG_CF_10_filtered$SNP_percentage_CF_10*100)/(100/isolates_CF_10))
ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_filtered_third_variant$Isolates_CF_10 <- round((data_algG_CF_10_filtered_third_variant$SNP_percentage_CF_10*100)/(100/isolates_CF_10)),z<-0)
ifelse(exists("data_algG_CF_10_filtered_third_variant_no_ref") == TRUE,data_algG_CF_10_filtered_third_variant_no_ref$Isolates_CF_10 <- round((data_algG_CF_10_filtered_third_variant_no_ref$SNP_percentage_CF_10*100)/(100/isolates_CF_10)),z<-0)



# Clean data

data_algG_CF_10_clean <- data_algG_CF_10_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_clean_third_base <- data_algG_CF_10_filtered_third_variant[,which(names(data_algG_CF_10_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_10","Isolates_CF_10"))],z<-0)
ifelse(exists("data_algG_CF_10_filtered_third_variant") == TRUE,data_algG_CF_10_clean <- rbind(data_algG_CF_10_clean,data_algG_CF_10_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_10_filtered_third_variant_no_ref") == TRUE,data_algG_CF_10_clean_third_base_no_ref <- data_algG_CF_10_filtered_third_variant_no_ref[,which(names(data_algG_CF_10_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_10","Isolates_CF_10"))],z<-0)
ifelse(exists("data_algG_CF_10_filtered_third_variant_no_ref") == TRUE,data_algG_CF_10_clean <- rbind(data_algG_CF_10_clean,data_algG_CF_10_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_10_clean_20_filter <- merge(data_algG_CF_10_clean, data_algG_CF_10_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_10_clean_cleared <- data_algG_CF_10_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_10_clean_20_filter)){
  if(data_algG_CF_10_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_10_clean_20_filter$As_CF_10[i] >= 20){
    data_algG_CF_10_clean_cleared <- rbind(data_algG_CF_10_clean_cleared,data_algG_CF_10_clean_20_filter[i,])
  } else if(data_algG_CF_10_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_10_clean_20_filter$Cs_CF_10[i] >= 20){
    data_algG_CF_10_clean_cleared <- rbind(data_algG_CF_10_clean_cleared,data_algG_CF_10_clean_20_filter[i,])
  } else if(data_algG_CF_10_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_10_clean_20_filter$Gs_CF_10[i] >= 20){
    data_algG_CF_10_clean_cleared <- rbind(data_algG_CF_10_clean_cleared,data_algG_CF_10_clean_20_filter[i,])
  }else if(data_algG_CF_10_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_10_clean_20_filter$Ts_CF_10[i] >= 20){
    data_algG_CF_10_clean_cleared <- rbind(data_algG_CF_10_clean_cleared,data_algG_CF_10_clean_20_filter[i,])
  }
}

data_algG_CF_10_clean <- data_algG_CF_10_clean_cleared[,c(1:5)]


### CF_11

# Input 
data_algG_CF_11 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T11/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_11_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T11/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_11_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T11/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_11 <- data_algG_CF_11[!(data_algG_CF_11$Position %in% primer_correction_algG_CF_11_fwd$Position),]
data_algG_CF_11 <- data_algG_CF_11[!(data_algG_CF_11$Position %in% primer_correction_algG_CF_11_rev$Position),]
data_algG_CF_11 <- rbind(data_algG_CF_11,primer_correction_algG_CF_11_fwd,primer_correction_algG_CF_11_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_11 <- subset(data_algG_CF_11, coverage_CF_11 >199)
data_algG_CF_11_5UTR <- subset(data_algG_CF_11_5UTR, coverage_CF_11_5UTR >199)
data_algG_CF_11_3UTR <- subset(data_algG_CF_11_3UTR, coverage_CF_11_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_11_only_SNP <- subset(data_algG_CF_11, second_percentage <= 0.025) 
data_algG_CF_11_only_SNP$main_base_CF_11 <- colnames(data_algG_CF_11_only_SNP[, 8:11])[apply(data_algG_CF_11_only_SNP[, 8:11],1,which.max)]
data_algG_CF_11_only_SNP$main_base_CF_11 <- ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_As_CF_11","A", ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_Cs_CF_11","C", ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_Gs_CF_11","G", ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_Ts_CF_11","T",z<-0))))
data_algG_CF_11_only_SNP <- subset(data_algG_CF_11_only_SNP,GenBase_PA14 != main_base_CF_11)
ifelse(nrow(data_algG_CF_11_only_SNP) > 0,data_algG_CF_11_only_SNP$Isolates_CF_11 <- isolates_CF_11, z <- 0)


ifelse(nrow(data_algG_CF_11_only_SNP) > 0,data_algG_CF_11_only_SNP_clean <- data_algG_CF_11_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_11_only_SNP) > 0,names(data_algG_CF_11_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_11","Isolates_CF_11"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_11_filtered <- subset(data_algG_CF_11, second_percentage >= 0.025) 
data_algG_CF_11_filtered$main_base_CF_11 <- colnames(data_algG_CF_11_filtered[, 8:11])[apply(data_algG_CF_11_filtered[, 8:11],1,which.max)]
data_algG_CF_11_filtered$second_base_CF_11 <- colnames(data_algG_CF_11_filtered[, 8:11])[apply(data_algG_CF_11_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_11_filtered$third_base_CF_11 <- ifelse(apply(data_algG_CF_11_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_11_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_11_filtered$`PA14-Koordinate`)
  if(data_algG_CF_11_filtered$main_base_CF_11[a] == "percent_As_CF_11"){
    data_algG_CF_11_filtered$main_base_CF_11[a] <- "A"                     
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] == "percent_Cs_CF_11"){
    data_algG_CF_11_filtered$main_base_CF_11[a] <- "C"                     
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] == "percent_Gs_CF_11"){
    data_algG_CF_11_filtered$main_base_CF_11[a] <- "G"                     
  } else{
    data_algG_CF_11_filtered$main_base_CF_11[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_11_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_11_filtered$`PA14-Koordinate`)
  if(data_algG_CF_11_filtered$second_base_CF_11[a] == "percent_As_CF_11"){
    data_algG_CF_11_filtered$second_base_CF_11[a] <- "A"                     
  } else if(data_algG_CF_11_filtered$second_base_CF_11[a] == "percent_Cs_CF_11"){
    data_algG_CF_11_filtered$second_base_CF_11[a] <- "C"                     
  } else if(data_algG_CF_11_filtered$second_base_CF_11[a] == "percent_Gs_CF_11"){
    data_algG_CF_11_filtered$second_base_CF_11[a] <- "G"                     
  } else{
    data_algG_CF_11_filtered$second_base_CF_11[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_11_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_11_filtered$`PA14-Koordinate`)
  if(data_algG_CF_11_filtered$main_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$third_base_CF_11[a] == 0){
    data_algG_CF_11_filtered$third_base_CF_11[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_11_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_11_filtered$`PA14-Koordinate`)
  if(data_algG_CF_11_filtered$third_base_CF_11[a] == 1){
    data_algG_CF_11_filtered_third_variant <- data_algG_CF_11_filtered[a,]
    data_algG_CF_11_filtered$third_base_CF_11[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_11_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_11_filtered$`PA14-Koordinate`)
  if(data_algG_CF_11_filtered$third_base_CF_11[a] == 2){
    data_algG_CF_11_filtered_third_variant_no_ref <- data_algG_CF_11_filtered[a,]
    data_algG_CF_11_filtered$third_base_CF_11[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_11_filtered$third_base_CF_11) == 0, data_algG_CF_11_filtered <- subset(data_algG_CF_11_filtered, select = -c(third_base_CF_11)),data_algG_CF_11_filtered <- data_algG_CF_11_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_filtered_third_variant$third_percentage <- apply(data_algG_CF_11_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_11_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_filtered_third_variant$third_base_CF_11 <- colnames(data_algG_CF_11_filtered_third_variant[, 8:11])[apply(data_algG_CF_11_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_As_CF_11","A", ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_Cs_CF_11","C", ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_Gs_CF_11","G", ifelse(data_algG_CF_11_only_SNP$main_base_CF_11 == "percent_Ts_CF_11","T",z<-0))))

ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,
       data_algG_CF_11_filtered_third_variant$third_base_CF_11 <- 
         ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "percent_As_CF_11","A", 
                ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "percent_Cs_CF_11","C",
                       ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "percent_Gs_CF_11","G",
                              ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "percent_Ts_CF_11","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,
       data_algG_CF_11_filtered_third_variant$third_base_CF_11 <- ifelse(
         data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14,data_algG_CF_11_filtered_third_variant$third_base_CF_11 <- data_algG_CF_11_filtered_third_variant$second_base_CF_11,data_algG_CF_11_filtered_third_variant$third_base_CF_11
       ),
       z<-0
)

ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,
       data_algG_CF_11_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14, data_algG_CF_11_filtered_third_variant$second_percentage,data_algG_CF_11_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_filtered_third_variant <- data_algG_CF_11_filtered_third_variant[,-which(names(data_algG_CF_11_filtered_third_variant) %in% c("second_percentage","second_base_CF_11"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_11_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_11_filtered_third_variant_no_ref)[names(data_algG_CF_11_filtered_third_variant_no_ref) == "second_base_CF_11"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_11_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_11_filtered_third_variant_no_ref)[names(data_algG_CF_11_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_11",z<-0)




# Define SNP percentage & Add it

data_algG_CF_11_filtered$SNP_base <- 0
data_algG_CF_11_filtered$SNP_percentage_CF_11 <- 0


for(i in data_algG_CF_11_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_11_filtered$`PA14-Koordinate`)
  if(data_algG_CF_11_filtered$main_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] == "A"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_As_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] == "C"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Cs_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] == "G"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Gs_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] == "T"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Ts_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_11_filtered$second_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "A"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_As_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_11_filtered$second_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "C"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Cs_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_11_filtered$second_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "G"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Gs_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_11_filtered$second_base_CF_11[a] == data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "T"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Ts_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "A"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_As_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "C"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Cs_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "G"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Gs_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_11_filtered$main_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$second_base_CF_11[a] != data_algG_CF_11_filtered$GenBase_PA14[a] & data_algG_CF_11_filtered$main_base_CF_11[a] == "T"){
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- data_algG_CF_11_filtered$percent_Ts_CF_11[a]
    data_algG_CF_11_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_11_filtered$SNP_percentage_CF_11[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_filtered_third_variant$SNP_percentage_CF_11 <- 0,z<-0)


ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,
       data_algG_CF_11_filtered_third_variant$SNP_percentage_CF_11 <- 
         ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "A",
                data_algG_CF_11_filtered_third_variant$percent_As_CF_11, 
                ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "C",
                       data_algG_CF_11_filtered_third_variant$percent_Cs_CF_11,
                       ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "G",
                              data_algG_CF_11_filtered_third_variant$percent_Gs_CF_11,
                              ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "T",
                                     data_algG_CF_11_filtered_third_variant$percent_Ts_CF_11,
                                     ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "A",
                                            data_algG_CF_11_filtered_third_variant$percent_As_CF_11,
                                            ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "C",
                                                   data_algG_CF_11_filtered_third_variant$percent_Cs_CF_11,
                                                   ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "G",
                                                          data_algG_CF_11_filtered_third_variant$percent_Gs_CF_11,
                                                          ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "T",
                                                                 data_algG_CF_11_filtered_third_variant$percent_Ts_CF_11,
                                                                 ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "A",
                                                                        data_algG_CF_11_filtered_third_variant$percent_As_CF_11,
                                                                        ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "C",
                                                                               data_algG_CF_11_filtered_third_variant$percent_Cs_CF_11,
                                                                               ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "G",
                                                                                      data_algG_CF_11_filtered_third_variant$percent_Gs_CF_11,
                                                                                      ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "T",
                                                                                             data_algG_CF_11_filtered_third_variant$percent_Ts_CF_11,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,
       data_algG_CF_11_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "A",
                "A", 
                ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "C",
                       "C",
                       ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "G",
                              "G",
                              ifelse(data_algG_CF_11_filtered_third_variant$main_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "T",
                                     "T",
                                     ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "A",
                                            "A",
                                            ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 == data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$main_base_CF_11 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_11_filtered_third_variant$third_base_CF_11 != data_algG_CF_11_filtered_third_variant$GenBase_PA14 & data_algG_CF_11_filtered_third_variant$third_base_CF_11 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_11_filtered$Isolates_CF_11 <- round((data_algG_CF_11_filtered$SNP_percentage_CF_11*100)/(100/isolates_CF_11))
ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_filtered_third_variant$Isolates_CF_11 <- round((data_algG_CF_11_filtered_third_variant$SNP_percentage_CF_11*100)/(100/isolates_CF_11)),z<-0)
ifelse(exists("data_algG_CF_11_filtered_third_variant_no_ref") == TRUE,data_algG_CF_11_filtered_third_variant_no_ref$Isolates_CF_11 <- round((data_algG_CF_11_filtered_third_variant_no_ref$SNP_percentage_CF_11*100)/(100/isolates_CF_11)),z<-0)



# Clean data

data_algG_CF_11_clean <- data_algG_CF_11_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_clean_third_base <- data_algG_CF_11_filtered_third_variant[,which(names(data_algG_CF_11_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_11","Isolates_CF_11"))],z<-0)
ifelse(exists("data_algG_CF_11_filtered_third_variant") == TRUE,data_algG_CF_11_clean <- rbind(data_algG_CF_11_clean,data_algG_CF_11_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_11_filtered_third_variant_no_ref") == TRUE,data_algG_CF_11_clean_third_base_no_ref <- data_algG_CF_11_filtered_third_variant_no_ref[,which(names(data_algG_CF_11_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_11","Isolates_CF_11"))],z<-0)
ifelse(exists("data_algG_CF_11_filtered_third_variant_no_ref") == TRUE,data_algG_CF_11_clean <- rbind(data_algG_CF_11_clean,data_algG_CF_11_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_11_clean_20_filter <- merge(data_algG_CF_11_clean, data_algG_CF_11_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_11_clean_cleared <- data_algG_CF_11_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_11_clean_20_filter)){
  if(data_algG_CF_11_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_11_clean_20_filter$As_CF_11[i] >= 20){
    data_algG_CF_11_clean_cleared <- rbind(data_algG_CF_11_clean_cleared,data_algG_CF_11_clean_20_filter[i,])
  } else if(data_algG_CF_11_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_11_clean_20_filter$Cs_CF_11[i] >= 20){
    data_algG_CF_11_clean_cleared <- rbind(data_algG_CF_11_clean_cleared,data_algG_CF_11_clean_20_filter[i,])
  } else if(data_algG_CF_11_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_11_clean_20_filter$Gs_CF_11[i] >= 20){
    data_algG_CF_11_clean_cleared <- rbind(data_algG_CF_11_clean_cleared,data_algG_CF_11_clean_20_filter[i,])
  }else if(data_algG_CF_11_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_11_clean_20_filter$Ts_CF_11[i] >= 20){
    data_algG_CF_11_clean_cleared <- rbind(data_algG_CF_11_clean_cleared,data_algG_CF_11_clean_20_filter[i,])
  }
}

data_algG_CF_11_clean <- data_algG_CF_11_clean_cleared[,c(1:5)]

### CF_12

# Input 
data_algG_CF_12 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T12/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_12_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T12/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_12_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T12/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_12 <- data_algG_CF_12[!(data_algG_CF_12$Position %in% primer_correction_algG_CF_12_fwd$Position),]
data_algG_CF_12 <- data_algG_CF_12[!(data_algG_CF_12$Position %in% primer_correction_algG_CF_12_rev$Position),]
data_algG_CF_12 <- rbind(data_algG_CF_12,primer_correction_algG_CF_12_fwd,primer_correction_algG_CF_12_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_12 <- subset(data_algG_CF_12, coverage_CF_12 >199)
data_algG_CF_12_5UTR <- subset(data_algG_CF_12_5UTR, coverage_CF_12_5UTR >199)
data_algG_CF_12_3UTR <- subset(data_algG_CF_12_3UTR, coverage_CF_12_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_12_only_SNP <- subset(data_algG_CF_12, second_percentage <= 0.025) 
data_algG_CF_12_only_SNP$main_base_CF_12 <- colnames(data_algG_CF_12_only_SNP[, 8:11])[apply(data_algG_CF_12_only_SNP[, 8:11],1,which.max)]
data_algG_CF_12_only_SNP$main_base_CF_12 <- ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_As_CF_12","A", ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_Cs_CF_12","C", ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_Gs_CF_12","G", ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_Ts_CF_12","T",z<-0))))
data_algG_CF_12_only_SNP <- subset(data_algG_CF_12_only_SNP,GenBase_PA14 != main_base_CF_12)
ifelse(nrow(data_algG_CF_12_only_SNP) > 0,data_algG_CF_12_only_SNP$Isolates_CF_12 <- isolates_CF_12, z <- 0)


ifelse(nrow(data_algG_CF_12_only_SNP) > 0,data_algG_CF_12_only_SNP_clean <- data_algG_CF_12_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_12_only_SNP) > 0,names(data_algG_CF_12_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_12","Isolates_CF_12"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_12_filtered <- subset(data_algG_CF_12, second_percentage >= 0.025) 
data_algG_CF_12_filtered$main_base_CF_12 <- colnames(data_algG_CF_12_filtered[, 8:11])[apply(data_algG_CF_12_filtered[, 8:11],1,which.max)]
data_algG_CF_12_filtered$second_base_CF_12 <- colnames(data_algG_CF_12_filtered[, 8:11])[apply(data_algG_CF_12_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_12_filtered$third_base_CF_12 <- ifelse(apply(data_algG_CF_12_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_12_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_12_filtered$`PA14-Koordinate`)
  if(data_algG_CF_12_filtered$main_base_CF_12[a] == "percent_As_CF_12"){
    data_algG_CF_12_filtered$main_base_CF_12[a] <- "A"                     
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] == "percent_Cs_CF_12"){
    data_algG_CF_12_filtered$main_base_CF_12[a] <- "C"                     
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] == "percent_Gs_CF_12"){
    data_algG_CF_12_filtered$main_base_CF_12[a] <- "G"                     
  } else{
    data_algG_CF_12_filtered$main_base_CF_12[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_12_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_12_filtered$`PA14-Koordinate`)
  if(data_algG_CF_12_filtered$second_base_CF_12[a] == "percent_As_CF_12"){
    data_algG_CF_12_filtered$second_base_CF_12[a] <- "A"                     
  } else if(data_algG_CF_12_filtered$second_base_CF_12[a] == "percent_Cs_CF_12"){
    data_algG_CF_12_filtered$second_base_CF_12[a] <- "C"                     
  } else if(data_algG_CF_12_filtered$second_base_CF_12[a] == "percent_Gs_CF_12"){
    data_algG_CF_12_filtered$second_base_CF_12[a] <- "G"                     
  } else{
    data_algG_CF_12_filtered$second_base_CF_12[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_12_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_12_filtered$`PA14-Koordinate`)
  if(data_algG_CF_12_filtered$main_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$third_base_CF_12[a] == 0){
    data_algG_CF_12_filtered$third_base_CF_12[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_12_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_12_filtered$`PA14-Koordinate`)
  if(data_algG_CF_12_filtered$third_base_CF_12[a] == 1){
    data_algG_CF_12_filtered_third_variant <- data_algG_CF_12_filtered[a,]
    data_algG_CF_12_filtered$third_base_CF_12[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_12_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_12_filtered$`PA14-Koordinate`)
  if(data_algG_CF_12_filtered$third_base_CF_12[a] == 2){
    data_algG_CF_12_filtered_third_variant_no_ref <- data_algG_CF_12_filtered[a,]
    data_algG_CF_12_filtered$third_base_CF_12[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_12_filtered$third_base_CF_12) == 0, data_algG_CF_12_filtered <- subset(data_algG_CF_12_filtered, select = -c(third_base_CF_12)),data_algG_CF_12_filtered <- data_algG_CF_12_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_filtered_third_variant$third_percentage <- apply(data_algG_CF_12_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_12_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_filtered_third_variant$third_base_CF_12 <- colnames(data_algG_CF_12_filtered_third_variant[, 8:11])[apply(data_algG_CF_12_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_As_CF_12","A", ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_Cs_CF_12","C", ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_Gs_CF_12","G", ifelse(data_algG_CF_12_only_SNP$main_base_CF_12 == "percent_Ts_CF_12","T",z<-0))))

ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,
       data_algG_CF_12_filtered_third_variant$third_base_CF_12 <- 
         ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "percent_As_CF_12","A", 
                ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "percent_Cs_CF_12","C",
                       ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "percent_Gs_CF_12","G",
                              ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "percent_Ts_CF_12","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,
       data_algG_CF_12_filtered_third_variant$third_base_CF_12 <- ifelse(
         data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14,data_algG_CF_12_filtered_third_variant$third_base_CF_12 <- data_algG_CF_12_filtered_third_variant$second_base_CF_12,data_algG_CF_12_filtered_third_variant$third_base_CF_12
       ),
       z<-0
)

ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,
       data_algG_CF_12_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14, data_algG_CF_12_filtered_third_variant$second_percentage,data_algG_CF_12_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_filtered_third_variant <- data_algG_CF_12_filtered_third_variant[,-which(names(data_algG_CF_12_filtered_third_variant) %in% c("second_percentage","second_base_CF_12"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_12_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_12_filtered_third_variant_no_ref)[names(data_algG_CF_12_filtered_third_variant_no_ref) == "second_base_CF_12"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_12_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_12_filtered_third_variant_no_ref)[names(data_algG_CF_12_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_12",z<-0)




# Define SNP percentage & Add it

data_algG_CF_12_filtered$SNP_base <- 0
data_algG_CF_12_filtered$SNP_percentage_CF_12 <- 0


for(i in data_algG_CF_12_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_12_filtered$`PA14-Koordinate`)
  if(data_algG_CF_12_filtered$main_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] == "A"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_As_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] == "C"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Cs_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] == "G"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Gs_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] == "T"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Ts_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_12_filtered$second_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "A"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_As_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_12_filtered$second_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "C"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Cs_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_12_filtered$second_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "G"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Gs_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_12_filtered$second_base_CF_12[a] == data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "T"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Ts_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "A"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_As_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "C"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Cs_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "G"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Gs_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_12_filtered$main_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$second_base_CF_12[a] != data_algG_CF_12_filtered$GenBase_PA14[a] & data_algG_CF_12_filtered$main_base_CF_12[a] == "T"){
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- data_algG_CF_12_filtered$percent_Ts_CF_12[a]
    data_algG_CF_12_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_12_filtered$SNP_percentage_CF_12[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_filtered_third_variant$SNP_percentage_CF_12 <- 0,z<-0)


ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,
       data_algG_CF_12_filtered_third_variant$SNP_percentage_CF_12 <- 
         ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "A",
                data_algG_CF_12_filtered_third_variant$percent_As_CF_12, 
                ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "C",
                       data_algG_CF_12_filtered_third_variant$percent_Cs_CF_12,
                       ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "G",
                              data_algG_CF_12_filtered_third_variant$percent_Gs_CF_12,
                              ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "T",
                                     data_algG_CF_12_filtered_third_variant$percent_Ts_CF_12,
                                     ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "A",
                                            data_algG_CF_12_filtered_third_variant$percent_As_CF_12,
                                            ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "C",
                                                   data_algG_CF_12_filtered_third_variant$percent_Cs_CF_12,
                                                   ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "G",
                                                          data_algG_CF_12_filtered_third_variant$percent_Gs_CF_12,
                                                          ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "T",
                                                                 data_algG_CF_12_filtered_third_variant$percent_Ts_CF_12,
                                                                 ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "A",
                                                                        data_algG_CF_12_filtered_third_variant$percent_As_CF_12,
                                                                        ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "C",
                                                                               data_algG_CF_12_filtered_third_variant$percent_Cs_CF_12,
                                                                               ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "G",
                                                                                      data_algG_CF_12_filtered_third_variant$percent_Gs_CF_12,
                                                                                      ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "T",
                                                                                             data_algG_CF_12_filtered_third_variant$percent_Ts_CF_12,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,
       data_algG_CF_12_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "A",
                "A", 
                ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "C",
                       "C",
                       ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "G",
                              "G",
                              ifelse(data_algG_CF_12_filtered_third_variant$main_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "T",
                                     "T",
                                     ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "A",
                                            "A",
                                            ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 == data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$main_base_CF_12 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_12_filtered_third_variant$third_base_CF_12 != data_algG_CF_12_filtered_third_variant$GenBase_PA14 & data_algG_CF_12_filtered_third_variant$third_base_CF_12 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_12_filtered$Isolates_CF_12 <- round((data_algG_CF_12_filtered$SNP_percentage_CF_12*100)/(100/isolates_CF_12))
ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_filtered_third_variant$Isolates_CF_12 <- round((data_algG_CF_12_filtered_third_variant$SNP_percentage_CF_12*100)/(100/isolates_CF_12)),z<-0)
ifelse(exists("data_algG_CF_12_filtered_third_variant_no_ref") == TRUE,data_algG_CF_12_filtered_third_variant_no_ref$Isolates_CF_12 <- round((data_algG_CF_12_filtered_third_variant_no_ref$SNP_percentage_CF_12*100)/(100/isolates_CF_12)),z<-0)



# Clean data

data_algG_CF_12_clean <- data_algG_CF_12_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_clean_third_base <- data_algG_CF_12_filtered_third_variant[,which(names(data_algG_CF_12_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_12","Isolates_CF_12"))],z<-0)
ifelse(exists("data_algG_CF_12_filtered_third_variant") == TRUE,data_algG_CF_12_clean <- rbind(data_algG_CF_12_clean,data_algG_CF_12_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_12_filtered_third_variant_no_ref") == TRUE,data_algG_CF_12_clean_third_base_no_ref <- data_algG_CF_12_filtered_third_variant_no_ref[,which(names(data_algG_CF_12_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_12","Isolates_CF_12"))],z<-0)
ifelse(exists("data_algG_CF_12_filtered_third_variant_no_ref") == TRUE,data_algG_CF_12_clean <- rbind(data_algG_CF_12_clean,data_algG_CF_12_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_12_clean_20_filter <- merge(data_algG_CF_12_clean, data_algG_CF_12_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_12_clean_cleared <- data_algG_CF_12_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_12_clean_20_filter)){
  if(data_algG_CF_12_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_12_clean_20_filter$As_CF_12[i] >= 20){
    data_algG_CF_12_clean_cleared <- rbind(data_algG_CF_12_clean_cleared,data_algG_CF_12_clean_20_filter[i,])
  } else if(data_algG_CF_12_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_12_clean_20_filter$Cs_CF_12[i] >= 20){
    data_algG_CF_12_clean_cleared <- rbind(data_algG_CF_12_clean_cleared,data_algG_CF_12_clean_20_filter[i,])
  } else if(data_algG_CF_12_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_12_clean_20_filter$Gs_CF_12[i] >= 20){
    data_algG_CF_12_clean_cleared <- rbind(data_algG_CF_12_clean_cleared,data_algG_CF_12_clean_20_filter[i,])
  }else if(data_algG_CF_12_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_12_clean_20_filter$Ts_CF_12[i] >= 20){
    data_algG_CF_12_clean_cleared <- rbind(data_algG_CF_12_clean_cleared,data_algG_CF_12_clean_20_filter[i,])
  }
}

data_algG_CF_12_clean <- data_algG_CF_12_clean_cleared[,c(1:5)]


### CF_13

# Input 
data_algG_CF_13 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T13/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_13_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T13/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_13_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T13/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_13 <- data_algG_CF_13[!(data_algG_CF_13$Position %in% primer_correction_algG_CF_13_fwd$Position),]
data_algG_CF_13 <- data_algG_CF_13[!(data_algG_CF_13$Position %in% primer_correction_algG_CF_13_rev$Position),]
data_algG_CF_13 <- rbind(data_algG_CF_13,primer_correction_algG_CF_13_fwd,primer_correction_algG_CF_13_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_13 <- subset(data_algG_CF_13, coverage_CF_13 >199)
data_algG_CF_13_5UTR <- subset(data_algG_CF_13_5UTR, coverage_CF_13_5UTR >199)
data_algG_CF_13_3UTR <- subset(data_algG_CF_13_3UTR, coverage_CF_13_3UTR >199)


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

# Get positions with only SNP

data_algG_CF_13_only_SNP <- subset(data_algG_CF_13, second_percentage <= 0.025) 
data_algG_CF_13_only_SNP$main_base_CF_13 <- colnames(data_algG_CF_13_only_SNP[, 8:11])[apply(data_algG_CF_13_only_SNP[, 8:11],1,which.max)]
data_algG_CF_13_only_SNP$main_base_CF_13 <- ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_As_CF_13","A", ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_Cs_CF_13","C", ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_Gs_CF_13","G", ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_Ts_CF_13","T",z<-0))))
data_algG_CF_13_only_SNP <- subset(data_algG_CF_13_only_SNP,GenBase_PA14 != main_base_CF_13)
ifelse(nrow(data_algG_CF_13_only_SNP) > 0,data_algG_CF_13_only_SNP$Isolates_CF_13 <- isolates_CF_13, z <- 0)


ifelse(nrow(data_algG_CF_13_only_SNP) > 0,data_algG_CF_13_only_SNP_clean <- data_algG_CF_13_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_13_only_SNP) > 0,names(data_algG_CF_13_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_13","Isolates_CF_13"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_13_filtered <- subset(data_algG_CF_13, second_percentage >= 0.025) 
data_algG_CF_13_filtered$main_base_CF_13 <- colnames(data_algG_CF_13_filtered[, 8:11])[apply(data_algG_CF_13_filtered[, 8:11],1,which.max)]
data_algG_CF_13_filtered$second_base_CF_13 <- colnames(data_algG_CF_13_filtered[, 8:11])[apply(data_algG_CF_13_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_13_filtered$third_base_CF_13 <- ifelse(apply(data_algG_CF_13_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_13_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_13_filtered$`PA14-Koordinate`)
  if(data_algG_CF_13_filtered$main_base_CF_13[a] == "percent_As_CF_13"){
    data_algG_CF_13_filtered$main_base_CF_13[a] <- "A"                     
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] == "percent_Cs_CF_13"){
    data_algG_CF_13_filtered$main_base_CF_13[a] <- "C"                     
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] == "percent_Gs_CF_13"){
    data_algG_CF_13_filtered$main_base_CF_13[a] <- "G"                     
  } else{
    data_algG_CF_13_filtered$main_base_CF_13[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_13_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_13_filtered$`PA14-Koordinate`)
  if(data_algG_CF_13_filtered$second_base_CF_13[a] == "percent_As_CF_13"){
    data_algG_CF_13_filtered$second_base_CF_13[a] <- "A"                     
  } else if(data_algG_CF_13_filtered$second_base_CF_13[a] == "percent_Cs_CF_13"){
    data_algG_CF_13_filtered$second_base_CF_13[a] <- "C"                     
  } else if(data_algG_CF_13_filtered$second_base_CF_13[a] == "percent_Gs_CF_13"){
    data_algG_CF_13_filtered$second_base_CF_13[a] <- "G"                     
  } else{
    data_algG_CF_13_filtered$second_base_CF_13[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_13_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_13_filtered$`PA14-Koordinate`)
  if(data_algG_CF_13_filtered$main_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$third_base_CF_13[a] == 0){
    data_algG_CF_13_filtered$third_base_CF_13[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_13_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_13_filtered$`PA14-Koordinate`)
  if(data_algG_CF_13_filtered$third_base_CF_13[a] == 1){
    data_algG_CF_13_filtered_third_variant <- data_algG_CF_13_filtered[a,]
    data_algG_CF_13_filtered$third_base_CF_13[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_13_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_13_filtered$`PA14-Koordinate`)
  if(data_algG_CF_13_filtered$third_base_CF_13[a] == 2){
    data_algG_CF_13_filtered_third_variant_no_ref <- data_algG_CF_13_filtered[a,]
    data_algG_CF_13_filtered$third_base_CF_13[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_13_filtered$third_base_CF_13) == 0, data_algG_CF_13_filtered <- subset(data_algG_CF_13_filtered, select = -c(third_base_CF_13)),data_algG_CF_13_filtered <- data_algG_CF_13_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_filtered_third_variant$third_percentage <- apply(data_algG_CF_13_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_13_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_filtered_third_variant$third_base_CF_13 <- colnames(data_algG_CF_13_filtered_third_variant[, 8:11])[apply(data_algG_CF_13_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_As_CF_13","A", ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_Cs_CF_13","C", ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_Gs_CF_13","G", ifelse(data_algG_CF_13_only_SNP$main_base_CF_13 == "percent_Ts_CF_13","T",z<-0))))

ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,
       data_algG_CF_13_filtered_third_variant$third_base_CF_13 <- 
         ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "percent_As_CF_13","A", 
                ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "percent_Cs_CF_13","C",
                       ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "percent_Gs_CF_13","G",
                              ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "percent_Ts_CF_13","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,
       data_algG_CF_13_filtered_third_variant$third_base_CF_13 <- ifelse(
         data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14,data_algG_CF_13_filtered_third_variant$third_base_CF_13 <- data_algG_CF_13_filtered_third_variant$second_base_CF_13,data_algG_CF_13_filtered_third_variant$third_base_CF_13
       ),
       z<-0
)

ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,
       data_algG_CF_13_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14, data_algG_CF_13_filtered_third_variant$second_percentage,data_algG_CF_13_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_filtered_third_variant <- data_algG_CF_13_filtered_third_variant[,-which(names(data_algG_CF_13_filtered_third_variant) %in% c("second_percentage","second_base_CF_13"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_13_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_13_filtered_third_variant_no_ref)[names(data_algG_CF_13_filtered_third_variant_no_ref) == "second_base_CF_13"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_13_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_13_filtered_third_variant_no_ref)[names(data_algG_CF_13_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_13",z<-0)




# Define SNP percentage & Add it

data_algG_CF_13_filtered$SNP_base <- 0
data_algG_CF_13_filtered$SNP_percentage_CF_13 <- 0


for(i in data_algG_CF_13_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_13_filtered$`PA14-Koordinate`)
  if(data_algG_CF_13_filtered$main_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] == "A"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_As_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] == "C"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Cs_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] == "G"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Gs_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] == "T"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Ts_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_13_filtered$second_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "A"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_As_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_13_filtered$second_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "C"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Cs_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_13_filtered$second_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "G"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Gs_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_13_filtered$second_base_CF_13[a] == data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "T"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Ts_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "A"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_As_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "C"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Cs_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "G"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Gs_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_13_filtered$main_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$second_base_CF_13[a] != data_algG_CF_13_filtered$GenBase_PA14[a] & data_algG_CF_13_filtered$main_base_CF_13[a] == "T"){
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- data_algG_CF_13_filtered$percent_Ts_CF_13[a]
    data_algG_CF_13_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_13_filtered$SNP_percentage_CF_13[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_filtered_third_variant$SNP_percentage_CF_13 <- 0,z<-0)


ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,
       data_algG_CF_13_filtered_third_variant$SNP_percentage_CF_13 <- 
         ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "A",
                data_algG_CF_13_filtered_third_variant$percent_As_CF_13, 
                ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "C",
                       data_algG_CF_13_filtered_third_variant$percent_Cs_CF_13,
                       ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "G",
                              data_algG_CF_13_filtered_third_variant$percent_Gs_CF_13,
                              ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "T",
                                     data_algG_CF_13_filtered_third_variant$percent_Ts_CF_13,
                                     ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "A",
                                            data_algG_CF_13_filtered_third_variant$percent_As_CF_13,
                                            ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "C",
                                                   data_algG_CF_13_filtered_third_variant$percent_Cs_CF_13,
                                                   ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "G",
                                                          data_algG_CF_13_filtered_third_variant$percent_Gs_CF_13,
                                                          ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "T",
                                                                 data_algG_CF_13_filtered_third_variant$percent_Ts_CF_13,
                                                                 ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "A",
                                                                        data_algG_CF_13_filtered_third_variant$percent_As_CF_13,
                                                                        ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "C",
                                                                               data_algG_CF_13_filtered_third_variant$percent_Cs_CF_13,
                                                                               ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "G",
                                                                                      data_algG_CF_13_filtered_third_variant$percent_Gs_CF_13,
                                                                                      ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "T",
                                                                                             data_algG_CF_13_filtered_third_variant$percent_Ts_CF_13,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,
       data_algG_CF_13_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "A",
                "A", 
                ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "C",
                       "C",
                       ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "G",
                              "G",
                              ifelse(data_algG_CF_13_filtered_third_variant$main_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "T",
                                     "T",
                                     ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "A",
                                            "A",
                                            ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 == data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$main_base_CF_13 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_13_filtered_third_variant$third_base_CF_13 != data_algG_CF_13_filtered_third_variant$GenBase_PA14 & data_algG_CF_13_filtered_third_variant$third_base_CF_13 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_13_filtered$Isolates_CF_13 <- round((data_algG_CF_13_filtered$SNP_percentage_CF_13*100)/(100/isolates_CF_13))
ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_filtered_third_variant$Isolates_CF_13 <- round((data_algG_CF_13_filtered_third_variant$SNP_percentage_CF_13*100)/(100/isolates_CF_13)),z<-0)
ifelse(exists("data_algG_CF_13_filtered_third_variant_no_ref") == TRUE,data_algG_CF_13_filtered_third_variant_no_ref$Isolates_CF_13 <- round((data_algG_CF_13_filtered_third_variant_no_ref$SNP_percentage_CF_13*100)/(100/isolates_CF_13)),z<-0)



# Clean data

data_algG_CF_13_clean <- data_algG_CF_13_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_clean_third_base <- data_algG_CF_13_filtered_third_variant[,which(names(data_algG_CF_13_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_13","Isolates_CF_13"))],z<-0)
ifelse(exists("data_algG_CF_13_filtered_third_variant") == TRUE,data_algG_CF_13_clean <- rbind(data_algG_CF_13_clean,data_algG_CF_13_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_13_filtered_third_variant_no_ref") == TRUE,data_algG_CF_13_clean_third_base_no_ref <- data_algG_CF_13_filtered_third_variant_no_ref[,which(names(data_algG_CF_13_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_13","Isolates_CF_13"))],z<-0)
ifelse(exists("data_algG_CF_13_filtered_third_variant_no_ref") == TRUE,data_algG_CF_13_clean <- rbind(data_algG_CF_13_clean,data_algG_CF_13_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_13_clean_20_filter <- merge(data_algG_CF_13_clean, data_algG_CF_13_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_13_clean_cleared <- data_algG_CF_13_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_13_clean_20_filter)){
  if(data_algG_CF_13_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_13_clean_20_filter$As_CF_13[i] >= 20){
    data_algG_CF_13_clean_cleared <- rbind(data_algG_CF_13_clean_cleared,data_algG_CF_13_clean_20_filter[i,])
  } else if(data_algG_CF_13_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_13_clean_20_filter$Cs_CF_13[i] >= 20){
    data_algG_CF_13_clean_cleared <- rbind(data_algG_CF_13_clean_cleared,data_algG_CF_13_clean_20_filter[i,])
  } else if(data_algG_CF_13_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_13_clean_20_filter$Gs_CF_13[i] >= 20){
    data_algG_CF_13_clean_cleared <- rbind(data_algG_CF_13_clean_cleared,data_algG_CF_13_clean_20_filter[i,])
  }else if(data_algG_CF_13_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_13_clean_20_filter$Ts_CF_13[i] >= 20){
    data_algG_CF_13_clean_cleared <- rbind(data_algG_CF_13_clean_cleared,data_algG_CF_13_clean_20_filter[i,])
  }
}

data_algG_CF_13_clean <- data_algG_CF_13_clean_cleared[,c(1:5)]

### CF_14

# Input 
data_algG_CF_14 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T14/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_14_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T14/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_14_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T14/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_14 <- data_algG_CF_14[!(data_algG_CF_14$Position %in% primer_correction_algG_CF_14_fwd$Position),]
data_algG_CF_14 <- data_algG_CF_14[!(data_algG_CF_14$Position %in% primer_correction_algG_CF_14_rev$Position),]
data_algG_CF_14 <- rbind(data_algG_CF_14,primer_correction_algG_CF_14_fwd,primer_correction_algG_CF_14_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_14 <- subset(data_algG_CF_14, coverage_CF_14 >199)
data_algG_CF_14_5UTR <- subset(data_algG_CF_14_5UTR, coverage_CF_14_5UTR >199)
data_algG_CF_14_3UTR <- subset(data_algG_CF_14_3UTR, coverage_CF_14_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_14_only_SNP <- subset(data_algG_CF_14, second_percentage <= 0.025) 
data_algG_CF_14_only_SNP$main_base_CF_14 <- colnames(data_algG_CF_14_only_SNP[, 8:11])[apply(data_algG_CF_14_only_SNP[, 8:11],1,which.max)]
data_algG_CF_14_only_SNP$main_base_CF_14 <- ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_As_CF_14","A", ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_Cs_CF_14","C", ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_Gs_CF_14","G", ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_Ts_CF_14","T",z<-0))))
data_algG_CF_14_only_SNP <- subset(data_algG_CF_14_only_SNP,GenBase_PA14 != main_base_CF_14)
ifelse(nrow(data_algG_CF_14_only_SNP) > 0,data_algG_CF_14_only_SNP$Isolates_CF_14 <- isolates_CF_14, z <- 0)


ifelse(nrow(data_algG_CF_14_only_SNP) > 0,data_algG_CF_14_only_SNP_clean <- data_algG_CF_14_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_14_only_SNP) > 0,names(data_algG_CF_14_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_14","Isolates_CF_14"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_14_filtered <- subset(data_algG_CF_14, second_percentage >= 0.025) 
data_algG_CF_14_filtered$main_base_CF_14 <- colnames(data_algG_CF_14_filtered[, 8:11])[apply(data_algG_CF_14_filtered[, 8:11],1,which.max)]
data_algG_CF_14_filtered$second_base_CF_14 <- colnames(data_algG_CF_14_filtered[, 8:11])[apply(data_algG_CF_14_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_14_filtered$third_base_CF_14 <- ifelse(apply(data_algG_CF_14_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_14_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_14_filtered$`PA14-Koordinate`)
  if(data_algG_CF_14_filtered$main_base_CF_14[a] == "percent_As_CF_14"){
    data_algG_CF_14_filtered$main_base_CF_14[a] <- "A"                     
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] == "percent_Cs_CF_14"){
    data_algG_CF_14_filtered$main_base_CF_14[a] <- "C"                     
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] == "percent_Gs_CF_14"){
    data_algG_CF_14_filtered$main_base_CF_14[a] <- "G"                     
  } else{
    data_algG_CF_14_filtered$main_base_CF_14[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_14_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_14_filtered$`PA14-Koordinate`)
  if(data_algG_CF_14_filtered$second_base_CF_14[a] == "percent_As_CF_14"){
    data_algG_CF_14_filtered$second_base_CF_14[a] <- "A"                     
  } else if(data_algG_CF_14_filtered$second_base_CF_14[a] == "percent_Cs_CF_14"){
    data_algG_CF_14_filtered$second_base_CF_14[a] <- "C"                     
  } else if(data_algG_CF_14_filtered$second_base_CF_14[a] == "percent_Gs_CF_14"){
    data_algG_CF_14_filtered$second_base_CF_14[a] <- "G"                     
  } else{
    data_algG_CF_14_filtered$second_base_CF_14[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_14_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_14_filtered$`PA14-Koordinate`)
  if(data_algG_CF_14_filtered$main_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$third_base_CF_14[a] == 0){
    data_algG_CF_14_filtered$third_base_CF_14[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_14_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_14_filtered$`PA14-Koordinate`)
  if(data_algG_CF_14_filtered$third_base_CF_14[a] == 1){
    data_algG_CF_14_filtered_third_variant <- data_algG_CF_14_filtered[a,]
    data_algG_CF_14_filtered$third_base_CF_14[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_14_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_14_filtered$`PA14-Koordinate`)
  if(data_algG_CF_14_filtered$third_base_CF_14[a] == 2){
    data_algG_CF_14_filtered_third_variant_no_ref <- data_algG_CF_14_filtered[a,]
    data_algG_CF_14_filtered$third_base_CF_14[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_14_filtered$third_base_CF_14) == 0, data_algG_CF_14_filtered <- subset(data_algG_CF_14_filtered, select = -c(third_base_CF_14)),data_algG_CF_14_filtered <- data_algG_CF_14_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_filtered_third_variant$third_percentage <- apply(data_algG_CF_14_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_14_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_filtered_third_variant$third_base_CF_14 <- colnames(data_algG_CF_14_filtered_third_variant[, 8:11])[apply(data_algG_CF_14_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_As_CF_14","A", ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_Cs_CF_14","C", ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_Gs_CF_14","G", ifelse(data_algG_CF_14_only_SNP$main_base_CF_14 == "percent_Ts_CF_14","T",z<-0))))

ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,
       data_algG_CF_14_filtered_third_variant$third_base_CF_14 <- 
         ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "percent_As_CF_14","A", 
                ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "percent_Cs_CF_14","C",
                       ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "percent_Gs_CF_14","G",
                              ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "percent_Ts_CF_14","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,
       data_algG_CF_14_filtered_third_variant$third_base_CF_14 <- ifelse(
         data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14,data_algG_CF_14_filtered_third_variant$third_base_CF_14 <- data_algG_CF_14_filtered_third_variant$second_base_CF_14,data_algG_CF_14_filtered_third_variant$third_base_CF_14
       ),
       z<-0
)

ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,
       data_algG_CF_14_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14, data_algG_CF_14_filtered_third_variant$second_percentage,data_algG_CF_14_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_filtered_third_variant <- data_algG_CF_14_filtered_third_variant[,-which(names(data_algG_CF_14_filtered_third_variant) %in% c("second_percentage","second_base_CF_14"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_14_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_14_filtered_third_variant_no_ref)[names(data_algG_CF_14_filtered_third_variant_no_ref) == "second_base_CF_14"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_14_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_14_filtered_third_variant_no_ref)[names(data_algG_CF_14_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_14",z<-0)




# Define SNP percentage & Add it

data_algG_CF_14_filtered$SNP_base <- 0
data_algG_CF_14_filtered$SNP_percentage_CF_14 <- 0


for(i in data_algG_CF_14_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_14_filtered$`PA14-Koordinate`)
  if(data_algG_CF_14_filtered$main_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] == "A"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_As_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] == "C"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Cs_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] == "G"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Gs_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] == "T"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Ts_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_14_filtered$second_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "A"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_As_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_14_filtered$second_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "C"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Cs_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_14_filtered$second_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "G"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Gs_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_14_filtered$second_base_CF_14[a] == data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "T"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Ts_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "A"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_As_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "C"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Cs_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "G"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Gs_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_14_filtered$main_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$second_base_CF_14[a] != data_algG_CF_14_filtered$GenBase_PA14[a] & data_algG_CF_14_filtered$main_base_CF_14[a] == "T"){
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- data_algG_CF_14_filtered$percent_Ts_CF_14[a]
    data_algG_CF_14_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_14_filtered$SNP_percentage_CF_14[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_filtered_third_variant$SNP_percentage_CF_14 <- 0,z<-0)


ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,
       data_algG_CF_14_filtered_third_variant$SNP_percentage_CF_14 <- 
         ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "A",
                data_algG_CF_14_filtered_third_variant$percent_As_CF_14, 
                ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "C",
                       data_algG_CF_14_filtered_third_variant$percent_Cs_CF_14,
                       ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "G",
                              data_algG_CF_14_filtered_third_variant$percent_Gs_CF_14,
                              ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "T",
                                     data_algG_CF_14_filtered_third_variant$percent_Ts_CF_14,
                                     ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "A",
                                            data_algG_CF_14_filtered_third_variant$percent_As_CF_14,
                                            ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "C",
                                                   data_algG_CF_14_filtered_third_variant$percent_Cs_CF_14,
                                                   ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "G",
                                                          data_algG_CF_14_filtered_third_variant$percent_Gs_CF_14,
                                                          ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "T",
                                                                 data_algG_CF_14_filtered_third_variant$percent_Ts_CF_14,
                                                                 ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "A",
                                                                        data_algG_CF_14_filtered_third_variant$percent_As_CF_14,
                                                                        ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "C",
                                                                               data_algG_CF_14_filtered_third_variant$percent_Cs_CF_14,
                                                                               ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "G",
                                                                                      data_algG_CF_14_filtered_third_variant$percent_Gs_CF_14,
                                                                                      ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "T",
                                                                                             data_algG_CF_14_filtered_third_variant$percent_Ts_CF_14,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,
       data_algG_CF_14_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "A",
                "A", 
                ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "C",
                       "C",
                       ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "G",
                              "G",
                              ifelse(data_algG_CF_14_filtered_third_variant$main_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "T",
                                     "T",
                                     ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "A",
                                            "A",
                                            ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 == data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$main_base_CF_14 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_14_filtered_third_variant$third_base_CF_14 != data_algG_CF_14_filtered_third_variant$GenBase_PA14 & data_algG_CF_14_filtered_third_variant$third_base_CF_14 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_14_filtered$Isolates_CF_14 <- round((data_algG_CF_14_filtered$SNP_percentage_CF_14*100)/(100/isolates_CF_14))
ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_filtered_third_variant$Isolates_CF_14 <- round((data_algG_CF_14_filtered_third_variant$SNP_percentage_CF_14*100)/(100/isolates_CF_14)),z<-0)
ifelse(exists("data_algG_CF_14_filtered_third_variant_no_ref") == TRUE,data_algG_CF_14_filtered_third_variant_no_ref$Isolates_CF_14 <- round((data_algG_CF_14_filtered_third_variant_no_ref$SNP_percentage_CF_14*100)/(100/isolates_CF_14)),z<-0)



# Clean data

data_algG_CF_14_clean <- data_algG_CF_14_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_clean_third_base <- data_algG_CF_14_filtered_third_variant[,which(names(data_algG_CF_14_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_14","Isolates_CF_14"))],z<-0)
ifelse(exists("data_algG_CF_14_filtered_third_variant") == TRUE,data_algG_CF_14_clean <- rbind(data_algG_CF_14_clean,data_algG_CF_14_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_14_filtered_third_variant_no_ref") == TRUE,data_algG_CF_14_clean_third_base_no_ref <- data_algG_CF_14_filtered_third_variant_no_ref[,which(names(data_algG_CF_14_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_14","Isolates_CF_14"))],z<-0)
ifelse(exists("data_algG_CF_14_filtered_third_variant_no_ref") == TRUE,data_algG_CF_14_clean <- rbind(data_algG_CF_14_clean,data_algG_CF_14_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_14_clean_20_filter <- merge(data_algG_CF_14_clean, data_algG_CF_14_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_14_clean_cleared <- data_algG_CF_14_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_14_clean_20_filter)){
  if(data_algG_CF_14_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_14_clean_20_filter$As_CF_14[i] >= 20){
    data_algG_CF_14_clean_cleared <- rbind(data_algG_CF_14_clean_cleared,data_algG_CF_14_clean_20_filter[i,])
  } else if(data_algG_CF_14_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_14_clean_20_filter$Cs_CF_14[i] >= 20){
    data_algG_CF_14_clean_cleared <- rbind(data_algG_CF_14_clean_cleared,data_algG_CF_14_clean_20_filter[i,])
  } else if(data_algG_CF_14_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_14_clean_20_filter$Gs_CF_14[i] >= 20){
    data_algG_CF_14_clean_cleared <- rbind(data_algG_CF_14_clean_cleared,data_algG_CF_14_clean_20_filter[i,])
  }else if(data_algG_CF_14_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_14_clean_20_filter$Ts_CF_14[i] >= 20){
    data_algG_CF_14_clean_cleared <- rbind(data_algG_CF_14_clean_cleared,data_algG_CF_14_clean_20_filter[i,])
  }
}

data_algG_CF_14_clean <- data_algG_CF_14_clean_cleared[,c(1:5)]


### CF_15

# Input 
data_algG_CF_15 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T15/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_15_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T15/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_15_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T15/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_15 <- data_algG_CF_15[!(data_algG_CF_15$Position %in% primer_correction_algG_CF_15_fwd$Position),]
data_algG_CF_15 <- data_algG_CF_15[!(data_algG_CF_15$Position %in% primer_correction_algG_CF_15_rev$Position),]
data_algG_CF_15 <- rbind(data_algG_CF_15,primer_correction_algG_CF_15_fwd,primer_correction_algG_CF_15_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_15 <- subset(data_algG_CF_15, coverage_CF_15 >199)
data_algG_CF_15_5UTR <- subset(data_algG_CF_15_5UTR, coverage_CF_15_5UTR >199)
data_algG_CF_15_3UTR <- subset(data_algG_CF_15_3UTR, coverage_CF_15_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_15_only_SNP <- subset(data_algG_CF_15, second_percentage <= 0.025) 
data_algG_CF_15_only_SNP$main_base_CF_15 <- colnames(data_algG_CF_15_only_SNP[, 8:11])[apply(data_algG_CF_15_only_SNP[, 8:11],1,which.max)]
data_algG_CF_15_only_SNP$main_base_CF_15 <- ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_As_CF_15","A", ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_Cs_CF_15","C", ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_Gs_CF_15","G", ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_Ts_CF_15","T",z<-0))))
data_algG_CF_15_only_SNP <- subset(data_algG_CF_15_only_SNP,GenBase_PA14 != main_base_CF_15)
ifelse(nrow(data_algG_CF_15_only_SNP) > 0,data_algG_CF_15_only_SNP$Isolates_CF_15 <- isolates_CF_15, z <- 0)


ifelse(nrow(data_algG_CF_15_only_SNP) > 0,data_algG_CF_15_only_SNP_clean <- data_algG_CF_15_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_15_only_SNP) > 0,names(data_algG_CF_15_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_15","Isolates_CF_15"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_15_filtered <- subset(data_algG_CF_15, second_percentage >= 0.025) 
data_algG_CF_15_filtered$main_base_CF_15 <- colnames(data_algG_CF_15_filtered[, 8:11])[apply(data_algG_CF_15_filtered[, 8:11],1,which.max)]
data_algG_CF_15_filtered$second_base_CF_15 <- colnames(data_algG_CF_15_filtered[, 8:11])[apply(data_algG_CF_15_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_15_filtered$third_base_CF_15 <- ifelse(apply(data_algG_CF_15_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_15_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_15_filtered$`PA14-Koordinate`)
  if(data_algG_CF_15_filtered$main_base_CF_15[a] == "percent_As_CF_15"){
    data_algG_CF_15_filtered$main_base_CF_15[a] <- "A"                     
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] == "percent_Cs_CF_15"){
    data_algG_CF_15_filtered$main_base_CF_15[a] <- "C"                     
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] == "percent_Gs_CF_15"){
    data_algG_CF_15_filtered$main_base_CF_15[a] <- "G"                     
  } else{
    data_algG_CF_15_filtered$main_base_CF_15[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_15_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_15_filtered$`PA14-Koordinate`)
  if(data_algG_CF_15_filtered$second_base_CF_15[a] == "percent_As_CF_15"){
    data_algG_CF_15_filtered$second_base_CF_15[a] <- "A"                     
  } else if(data_algG_CF_15_filtered$second_base_CF_15[a] == "percent_Cs_CF_15"){
    data_algG_CF_15_filtered$second_base_CF_15[a] <- "C"                     
  } else if(data_algG_CF_15_filtered$second_base_CF_15[a] == "percent_Gs_CF_15"){
    data_algG_CF_15_filtered$second_base_CF_15[a] <- "G"                     
  } else{
    data_algG_CF_15_filtered$second_base_CF_15[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_15_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_15_filtered$`PA14-Koordinate`)
  if(data_algG_CF_15_filtered$main_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$third_base_CF_15[a] == 0){
    data_algG_CF_15_filtered$third_base_CF_15[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_15_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_15_filtered$`PA14-Koordinate`)
  if(data_algG_CF_15_filtered$third_base_CF_15[a] == 1){
    data_algG_CF_15_filtered_third_variant <- data_algG_CF_15_filtered[a,]
    data_algG_CF_15_filtered$third_base_CF_15[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_15_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_15_filtered$`PA14-Koordinate`)
  if(data_algG_CF_15_filtered$third_base_CF_15[a] == 2){
    data_algG_CF_15_filtered_third_variant_no_ref <- data_algG_CF_15_filtered[a,]
    data_algG_CF_15_filtered$third_base_CF_15[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_15_filtered$third_base_CF_15) == 0, data_algG_CF_15_filtered <- subset(data_algG_CF_15_filtered, select = -c(third_base_CF_15)),data_algG_CF_15_filtered <- data_algG_CF_15_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_filtered_third_variant$third_percentage <- apply(data_algG_CF_15_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_15_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_filtered_third_variant$third_base_CF_15 <- colnames(data_algG_CF_15_filtered_third_variant[, 8:11])[apply(data_algG_CF_15_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_As_CF_15","A", ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_Cs_CF_15","C", ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_Gs_CF_15","G", ifelse(data_algG_CF_15_only_SNP$main_base_CF_15 == "percent_Ts_CF_15","T",z<-0))))

ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,
       data_algG_CF_15_filtered_third_variant$third_base_CF_15 <- 
         ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "percent_As_CF_15","A", 
                ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "percent_Cs_CF_15","C",
                       ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "percent_Gs_CF_15","G",
                              ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "percent_Ts_CF_15","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,
       data_algG_CF_15_filtered_third_variant$third_base_CF_15 <- ifelse(
         data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14,data_algG_CF_15_filtered_third_variant$third_base_CF_15 <- data_algG_CF_15_filtered_third_variant$second_base_CF_15,data_algG_CF_15_filtered_third_variant$third_base_CF_15
       ),
       z<-0
)

ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,
       data_algG_CF_15_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14, data_algG_CF_15_filtered_third_variant$second_percentage,data_algG_CF_15_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_filtered_third_variant <- data_algG_CF_15_filtered_third_variant[,-which(names(data_algG_CF_15_filtered_third_variant) %in% c("second_percentage","second_base_CF_15"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_15_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_15_filtered_third_variant_no_ref)[names(data_algG_CF_15_filtered_third_variant_no_ref) == "second_base_CF_15"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_15_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_15_filtered_third_variant_no_ref)[names(data_algG_CF_15_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_15",z<-0)




# Define SNP percentage & Add it

data_algG_CF_15_filtered$SNP_base <- 0
data_algG_CF_15_filtered$SNP_percentage_CF_15 <- 0


for(i in data_algG_CF_15_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_15_filtered$`PA14-Koordinate`)
  if(data_algG_CF_15_filtered$main_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] == "A"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_As_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] == "C"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Cs_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] == "G"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Gs_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] == "T"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Ts_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_15_filtered$second_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "A"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_As_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_15_filtered$second_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "C"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Cs_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_15_filtered$second_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "G"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Gs_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_15_filtered$second_base_CF_15[a] == data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "T"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Ts_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "A"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_As_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "C"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Cs_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "G"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Gs_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_15_filtered$main_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$second_base_CF_15[a] != data_algG_CF_15_filtered$GenBase_PA14[a] & data_algG_CF_15_filtered$main_base_CF_15[a] == "T"){
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- data_algG_CF_15_filtered$percent_Ts_CF_15[a]
    data_algG_CF_15_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_15_filtered$SNP_percentage_CF_15[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_filtered_third_variant$SNP_percentage_CF_15 <- 0,z<-0)


ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,
       data_algG_CF_15_filtered_third_variant$SNP_percentage_CF_15 <- 
         ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "A",
                data_algG_CF_15_filtered_third_variant$percent_As_CF_15, 
                ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "C",
                       data_algG_CF_15_filtered_third_variant$percent_Cs_CF_15,
                       ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "G",
                              data_algG_CF_15_filtered_third_variant$percent_Gs_CF_15,
                              ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "T",
                                     data_algG_CF_15_filtered_third_variant$percent_Ts_CF_15,
                                     ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "A",
                                            data_algG_CF_15_filtered_third_variant$percent_As_CF_15,
                                            ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "C",
                                                   data_algG_CF_15_filtered_third_variant$percent_Cs_CF_15,
                                                   ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "G",
                                                          data_algG_CF_15_filtered_third_variant$percent_Gs_CF_15,
                                                          ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "T",
                                                                 data_algG_CF_15_filtered_third_variant$percent_Ts_CF_15,
                                                                 ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "A",
                                                                        data_algG_CF_15_filtered_third_variant$percent_As_CF_15,
                                                                        ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "C",
                                                                               data_algG_CF_15_filtered_third_variant$percent_Cs_CF_15,
                                                                               ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "G",
                                                                                      data_algG_CF_15_filtered_third_variant$percent_Gs_CF_15,
                                                                                      ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "T",
                                                                                             data_algG_CF_15_filtered_third_variant$percent_Ts_CF_15,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,
       data_algG_CF_15_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "A",
                "A", 
                ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "C",
                       "C",
                       ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "G",
                              "G",
                              ifelse(data_algG_CF_15_filtered_third_variant$main_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "T",
                                     "T",
                                     ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "A",
                                            "A",
                                            ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 == data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$main_base_CF_15 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_15_filtered_third_variant$third_base_CF_15 != data_algG_CF_15_filtered_third_variant$GenBase_PA14 & data_algG_CF_15_filtered_third_variant$third_base_CF_15 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_15_filtered$Isolates_CF_15 <- round((data_algG_CF_15_filtered$SNP_percentage_CF_15*100)/(100/isolates_CF_15))
ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_filtered_third_variant$Isolates_CF_15 <- round((data_algG_CF_15_filtered_third_variant$SNP_percentage_CF_15*100)/(100/isolates_CF_15)),z<-0)
ifelse(exists("data_algG_CF_15_filtered_third_variant_no_ref") == TRUE,data_algG_CF_15_filtered_third_variant_no_ref$Isolates_CF_15 <- round((data_algG_CF_15_filtered_third_variant_no_ref$SNP_percentage_CF_15*100)/(100/isolates_CF_15)),z<-0)



# Clean data

data_algG_CF_15_clean <- data_algG_CF_15_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_clean_third_base <- data_algG_CF_15_filtered_third_variant[,which(names(data_algG_CF_15_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_15","Isolates_CF_15"))],z<-0)
ifelse(exists("data_algG_CF_15_filtered_third_variant") == TRUE,data_algG_CF_15_clean <- rbind(data_algG_CF_15_clean,data_algG_CF_15_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_15_filtered_third_variant_no_ref") == TRUE,data_algG_CF_15_clean_third_base_no_ref <- data_algG_CF_15_filtered_third_variant_no_ref[,which(names(data_algG_CF_15_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_15","Isolates_CF_15"))],z<-0)
ifelse(exists("data_algG_CF_15_filtered_third_variant_no_ref") == TRUE,data_algG_CF_15_clean <- rbind(data_algG_CF_15_clean,data_algG_CF_15_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_15_clean_20_filter <- merge(data_algG_CF_15_clean, data_algG_CF_15_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_15_clean_cleared <- data_algG_CF_15_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_15_clean_20_filter)){
  if(data_algG_CF_15_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_15_clean_20_filter$As_CF_15[i] >= 20){
    data_algG_CF_15_clean_cleared <- rbind(data_algG_CF_15_clean_cleared,data_algG_CF_15_clean_20_filter[i,])
  } else if(data_algG_CF_15_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_15_clean_20_filter$Cs_CF_15[i] >= 20){
    data_algG_CF_15_clean_cleared <- rbind(data_algG_CF_15_clean_cleared,data_algG_CF_15_clean_20_filter[i,])
  } else if(data_algG_CF_15_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_15_clean_20_filter$Gs_CF_15[i] >= 20){
    data_algG_CF_15_clean_cleared <- rbind(data_algG_CF_15_clean_cleared,data_algG_CF_15_clean_20_filter[i,])
  }else if(data_algG_CF_15_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_15_clean_20_filter$Ts_CF_15[i] >= 20){
    data_algG_CF_15_clean_cleared <- rbind(data_algG_CF_15_clean_cleared,data_algG_CF_15_clean_20_filter[i,])
  }
}

data_algG_CF_15_clean <- data_algG_CF_15_clean_cleared[,c(1:5)]


### CF_16

# Input 
data_algG_CF_16 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T16/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_16_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T16/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_16_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T16/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_16 <- data_algG_CF_16[!(data_algG_CF_16$Position %in% primer_correction_algG_CF_16_fwd$Position),]
data_algG_CF_16 <- data_algG_CF_16[!(data_algG_CF_16$Position %in% primer_correction_algG_CF_16_rev$Position),]
data_algG_CF_16 <- rbind(data_algG_CF_16,primer_correction_algG_CF_16_fwd,primer_correction_algG_CF_16_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_16 <- subset(data_algG_CF_16, coverage_CF_16 >199)
data_algG_CF_16_5UTR <- subset(data_algG_CF_16_5UTR, coverage_CF_16_5UTR >199)
data_algG_CF_16_3UTR <- subset(data_algG_CF_16_3UTR, coverage_CF_16_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_16_only_SNP <- subset(data_algG_CF_16, second_percentage <= 0.025) 
data_algG_CF_16_only_SNP$main_base_CF_16 <- colnames(data_algG_CF_16_only_SNP[, 8:11])[apply(data_algG_CF_16_only_SNP[, 8:11],1,which.max)]
data_algG_CF_16_only_SNP$main_base_CF_16 <- ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_As_CF_16","A", ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_Cs_CF_16","C", ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_Gs_CF_16","G", ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_Ts_CF_16","T",z<-0))))
data_algG_CF_16_only_SNP <- subset(data_algG_CF_16_only_SNP,GenBase_PA14 != main_base_CF_16)
ifelse(nrow(data_algG_CF_16_only_SNP) > 0,data_algG_CF_16_only_SNP$Isolates_CF_16 <- isolates_CF_16, z <- 0)


ifelse(nrow(data_algG_CF_16_only_SNP) > 0,data_algG_CF_16_only_SNP_clean <- data_algG_CF_16_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_16_only_SNP) > 0,names(data_algG_CF_16_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_16","Isolates_CF_16"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_16_filtered <- subset(data_algG_CF_16, second_percentage >= 0.025) 
data_algG_CF_16_filtered$main_base_CF_16 <- colnames(data_algG_CF_16_filtered[, 8:11])[apply(data_algG_CF_16_filtered[, 8:11],1,which.max)]
data_algG_CF_16_filtered$second_base_CF_16 <- colnames(data_algG_CF_16_filtered[, 8:11])[apply(data_algG_CF_16_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_16_filtered$third_base_CF_16 <- ifelse(apply(data_algG_CF_16_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_16_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_16_filtered$`PA14-Koordinate`)
  if(data_algG_CF_16_filtered$main_base_CF_16[a] == "percent_As_CF_16"){
    data_algG_CF_16_filtered$main_base_CF_16[a] <- "A"                     
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] == "percent_Cs_CF_16"){
    data_algG_CF_16_filtered$main_base_CF_16[a] <- "C"                     
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] == "percent_Gs_CF_16"){
    data_algG_CF_16_filtered$main_base_CF_16[a] <- "G"                     
  } else{
    data_algG_CF_16_filtered$main_base_CF_16[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_16_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_16_filtered$`PA14-Koordinate`)
  if(data_algG_CF_16_filtered$second_base_CF_16[a] == "percent_As_CF_16"){
    data_algG_CF_16_filtered$second_base_CF_16[a] <- "A"                     
  } else if(data_algG_CF_16_filtered$second_base_CF_16[a] == "percent_Cs_CF_16"){
    data_algG_CF_16_filtered$second_base_CF_16[a] <- "C"                     
  } else if(data_algG_CF_16_filtered$second_base_CF_16[a] == "percent_Gs_CF_16"){
    data_algG_CF_16_filtered$second_base_CF_16[a] <- "G"                     
  } else{
    data_algG_CF_16_filtered$second_base_CF_16[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_16_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_16_filtered$`PA14-Koordinate`)
  if(data_algG_CF_16_filtered$main_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$third_base_CF_16[a] == 0){
    data_algG_CF_16_filtered$third_base_CF_16[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_16_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_16_filtered$`PA14-Koordinate`)
  if(data_algG_CF_16_filtered$third_base_CF_16[a] == 1){
    data_algG_CF_16_filtered_third_variant <- data_algG_CF_16_filtered[a,]
    data_algG_CF_16_filtered$third_base_CF_16[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_16_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_16_filtered$`PA14-Koordinate`)
  if(data_algG_CF_16_filtered$third_base_CF_16[a] == 2){
    data_algG_CF_16_filtered_third_variant_no_ref <- data_algG_CF_16_filtered[a,]
    data_algG_CF_16_filtered$third_base_CF_16[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_16_filtered$third_base_CF_16) == 0, data_algG_CF_16_filtered <- subset(data_algG_CF_16_filtered, select = -c(third_base_CF_16)),data_algG_CF_16_filtered <- data_algG_CF_16_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_filtered_third_variant$third_percentage <- apply(data_algG_CF_16_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_16_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_filtered_third_variant$third_base_CF_16 <- colnames(data_algG_CF_16_filtered_third_variant[, 8:11])[apply(data_algG_CF_16_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_As_CF_16","A", ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_Cs_CF_16","C", ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_Gs_CF_16","G", ifelse(data_algG_CF_16_only_SNP$main_base_CF_16 == "percent_Ts_CF_16","T",z<-0))))

ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,
       data_algG_CF_16_filtered_third_variant$third_base_CF_16 <- 
         ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "percent_As_CF_16","A", 
                ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "percent_Cs_CF_16","C",
                       ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "percent_Gs_CF_16","G",
                              ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "percent_Ts_CF_16","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,
       data_algG_CF_16_filtered_third_variant$third_base_CF_16 <- ifelse(
         data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14,data_algG_CF_16_filtered_third_variant$third_base_CF_16 <- data_algG_CF_16_filtered_third_variant$second_base_CF_16,data_algG_CF_16_filtered_third_variant$third_base_CF_16
       ),
       z<-0
)

ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,
       data_algG_CF_16_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14, data_algG_CF_16_filtered_third_variant$second_percentage,data_algG_CF_16_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_filtered_third_variant <- data_algG_CF_16_filtered_third_variant[,-which(names(data_algG_CF_16_filtered_third_variant) %in% c("second_percentage","second_base_CF_16"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_16_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_16_filtered_third_variant_no_ref)[names(data_algG_CF_16_filtered_third_variant_no_ref) == "second_base_CF_16"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_16_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_16_filtered_third_variant_no_ref)[names(data_algG_CF_16_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_16",z<-0)




# Define SNP percentage & Add it

data_algG_CF_16_filtered$SNP_base <- 0
data_algG_CF_16_filtered$SNP_percentage_CF_16 <- 0


for(i in data_algG_CF_16_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_16_filtered$`PA14-Koordinate`)
  if(data_algG_CF_16_filtered$main_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] == "A"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_As_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] == "C"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Cs_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] == "G"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Gs_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] == "T"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Ts_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_16_filtered$second_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "A"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_As_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_16_filtered$second_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "C"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Cs_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_16_filtered$second_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "G"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Gs_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_16_filtered$second_base_CF_16[a] == data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "T"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Ts_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "A"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_As_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "C"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Cs_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "G"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Gs_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_16_filtered$main_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$second_base_CF_16[a] != data_algG_CF_16_filtered$GenBase_PA14[a] & data_algG_CF_16_filtered$main_base_CF_16[a] == "T"){
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- data_algG_CF_16_filtered$percent_Ts_CF_16[a]
    data_algG_CF_16_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_16_filtered$SNP_percentage_CF_16[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_filtered_third_variant$SNP_percentage_CF_16 <- 0,z<-0)


ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,
       data_algG_CF_16_filtered_third_variant$SNP_percentage_CF_16 <- 
         ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "A",
                data_algG_CF_16_filtered_third_variant$percent_As_CF_16, 
                ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "C",
                       data_algG_CF_16_filtered_third_variant$percent_Cs_CF_16,
                       ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "G",
                              data_algG_CF_16_filtered_third_variant$percent_Gs_CF_16,
                              ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "T",
                                     data_algG_CF_16_filtered_third_variant$percent_Ts_CF_16,
                                     ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "A",
                                            data_algG_CF_16_filtered_third_variant$percent_As_CF_16,
                                            ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "C",
                                                   data_algG_CF_16_filtered_third_variant$percent_Cs_CF_16,
                                                   ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "G",
                                                          data_algG_CF_16_filtered_third_variant$percent_Gs_CF_16,
                                                          ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "T",
                                                                 data_algG_CF_16_filtered_third_variant$percent_Ts_CF_16,
                                                                 ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "A",
                                                                        data_algG_CF_16_filtered_third_variant$percent_As_CF_16,
                                                                        ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "C",
                                                                               data_algG_CF_16_filtered_third_variant$percent_Cs_CF_16,
                                                                               ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "G",
                                                                                      data_algG_CF_16_filtered_third_variant$percent_Gs_CF_16,
                                                                                      ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "T",
                                                                                             data_algG_CF_16_filtered_third_variant$percent_Ts_CF_16,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,
       data_algG_CF_16_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "A",
                "A", 
                ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "C",
                       "C",
                       ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "G",
                              "G",
                              ifelse(data_algG_CF_16_filtered_third_variant$main_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "T",
                                     "T",
                                     ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "A",
                                            "A",
                                            ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 == data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$main_base_CF_16 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_16_filtered_third_variant$third_base_CF_16 != data_algG_CF_16_filtered_third_variant$GenBase_PA14 & data_algG_CF_16_filtered_third_variant$third_base_CF_16 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_16_filtered$Isolates_CF_16 <- round((data_algG_CF_16_filtered$SNP_percentage_CF_16*100)/(100/isolates_CF_16))
ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_filtered_third_variant$Isolates_CF_16 <- round((data_algG_CF_16_filtered_third_variant$SNP_percentage_CF_16*100)/(100/isolates_CF_16)),z<-0)
ifelse(exists("data_algG_CF_16_filtered_third_variant_no_ref") == TRUE,data_algG_CF_16_filtered_third_variant_no_ref$Isolates_CF_16 <- round((data_algG_CF_16_filtered_third_variant_no_ref$SNP_percentage_CF_16*100)/(100/isolates_CF_16)),z<-0)



# Clean data

data_algG_CF_16_clean <- data_algG_CF_16_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_clean_third_base <- data_algG_CF_16_filtered_third_variant[,which(names(data_algG_CF_16_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_16","Isolates_CF_16"))],z<-0)
ifelse(exists("data_algG_CF_16_filtered_third_variant") == TRUE,data_algG_CF_16_clean <- rbind(data_algG_CF_16_clean,data_algG_CF_16_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_16_filtered_third_variant_no_ref") == TRUE,data_algG_CF_16_clean_third_base_no_ref <- data_algG_CF_16_filtered_third_variant_no_ref[,which(names(data_algG_CF_16_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_16","Isolates_CF_16"))],z<-0)
ifelse(exists("data_algG_CF_16_filtered_third_variant_no_ref") == TRUE,data_algG_CF_16_clean <- rbind(data_algG_CF_16_clean,data_algG_CF_16_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_16_clean_20_filter <- merge(data_algG_CF_16_clean, data_algG_CF_16_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_16_clean_cleared <- data_algG_CF_16_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_16_clean_20_filter)){
  if(data_algG_CF_16_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_16_clean_20_filter$As_CF_16[i] >= 20){
    data_algG_CF_16_clean_cleared <- rbind(data_algG_CF_16_clean_cleared,data_algG_CF_16_clean_20_filter[i,])
  } else if(data_algG_CF_16_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_16_clean_20_filter$Cs_CF_16[i] >= 20){
    data_algG_CF_16_clean_cleared <- rbind(data_algG_CF_16_clean_cleared,data_algG_CF_16_clean_20_filter[i,])
  } else if(data_algG_CF_16_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_16_clean_20_filter$Gs_CF_16[i] >= 20){
    data_algG_CF_16_clean_cleared <- rbind(data_algG_CF_16_clean_cleared,data_algG_CF_16_clean_20_filter[i,])
  }else if(data_algG_CF_16_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_16_clean_20_filter$Ts_CF_16[i] >= 20){
    data_algG_CF_16_clean_cleared <- rbind(data_algG_CF_16_clean_cleared,data_algG_CF_16_clean_20_filter[i,])
  }
}

data_algG_CF_16_clean <- data_algG_CF_16_clean_cleared[,c(1:5)]


### CF_17

# Input 
data_algG_CF_17 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T17/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_CF_17_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T17/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_CF_17_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T17/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_CF_17 <- data_algG_CF_17[!(data_algG_CF_17$Position %in% primer_correction_algG_CF_17_fwd$Position),]
data_algG_CF_17 <- data_algG_CF_17[!(data_algG_CF_17$Position %in% primer_correction_algG_CF_17_rev$Position),]
data_algG_CF_17 <- rbind(data_algG_CF_17,primer_correction_algG_CF_17_fwd,primer_correction_algG_CF_17_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_CF_17 <- subset(data_algG_CF_17, coverage_CF_17 >199)
data_algG_CF_17_5UTR <- subset(data_algG_CF_17_5UTR, coverage_CF_17_5UTR >199)
data_algG_CF_17_3UTR <- subset(data_algG_CF_17_3UTR, coverage_CF_17_3UTR >199)

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

# Get positions with only SNP

data_algG_CF_17_only_SNP <- subset(data_algG_CF_17, second_percentage <= 0.025) 
data_algG_CF_17_only_SNP$main_base_CF_17 <- colnames(data_algG_CF_17_only_SNP[, 8:11])[apply(data_algG_CF_17_only_SNP[, 8:11],1,which.max)]
data_algG_CF_17_only_SNP$main_base_CF_17 <- ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_As_CF_17","A", ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_Cs_CF_17","C", ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_Gs_CF_17","G", ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_Ts_CF_17","T",z<-0))))
data_algG_CF_17_only_SNP <- subset(data_algG_CF_17_only_SNP,GenBase_PA14 != main_base_CF_17)
ifelse(nrow(data_algG_CF_17_only_SNP) > 0,data_algG_CF_17_only_SNP$Isolates_CF_17 <- isolates_CF_17, z <- 0)


ifelse(nrow(data_algG_CF_17_only_SNP) > 0,data_algG_CF_17_only_SNP_clean <- data_algG_CF_17_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_CF_17_only_SNP) > 0,names(data_algG_CF_17_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_17","Isolates_CF_17"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_CF_17_filtered <- subset(data_algG_CF_17, second_percentage >= 0.025) 
data_algG_CF_17_filtered$main_base_CF_17 <- colnames(data_algG_CF_17_filtered[, 8:11])[apply(data_algG_CF_17_filtered[, 8:11],1,which.max)]
data_algG_CF_17_filtered$second_base_CF_17 <- colnames(data_algG_CF_17_filtered[, 8:11])[apply(data_algG_CF_17_filtered[, 8:11], 1, maxn(2))]
data_algG_CF_17_filtered$third_base_CF_17 <- ifelse(apply(data_algG_CF_17_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_CF_17_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_17_filtered$`PA14-Koordinate`)
  if(data_algG_CF_17_filtered$main_base_CF_17[a] == "percent_As_CF_17"){
    data_algG_CF_17_filtered$main_base_CF_17[a] <- "A"                     
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] == "percent_Cs_CF_17"){
    data_algG_CF_17_filtered$main_base_CF_17[a] <- "C"                     
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] == "percent_Gs_CF_17"){
    data_algG_CF_17_filtered$main_base_CF_17[a] <- "G"                     
  } else{
    data_algG_CF_17_filtered$main_base_CF_17[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_CF_17_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_17_filtered$`PA14-Koordinate`)
  if(data_algG_CF_17_filtered$second_base_CF_17[a] == "percent_As_CF_17"){
    data_algG_CF_17_filtered$second_base_CF_17[a] <- "A"                     
  } else if(data_algG_CF_17_filtered$second_base_CF_17[a] == "percent_Cs_CF_17"){
    data_algG_CF_17_filtered$second_base_CF_17[a] <- "C"                     
  } else if(data_algG_CF_17_filtered$second_base_CF_17[a] == "percent_Gs_CF_17"){
    data_algG_CF_17_filtered$second_base_CF_17[a] <- "G"                     
  } else{
    data_algG_CF_17_filtered$second_base_CF_17[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_CF_17_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_17_filtered$`PA14-Koordinate`)
  if(data_algG_CF_17_filtered$main_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$third_base_CF_17[a] == 0){
    data_algG_CF_17_filtered$third_base_CF_17[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_CF_17_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_17_filtered$`PA14-Koordinate`)
  if(data_algG_CF_17_filtered$third_base_CF_17[a] == 1){
    data_algG_CF_17_filtered_third_variant <- data_algG_CF_17_filtered[a,]
    data_algG_CF_17_filtered$third_base_CF_17[a] <- 0
  }
  rm(a)
}


for(i in data_algG_CF_17_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_17_filtered$`PA14-Koordinate`)
  if(data_algG_CF_17_filtered$third_base_CF_17[a] == 2){
    data_algG_CF_17_filtered_third_variant_no_ref <- data_algG_CF_17_filtered[a,]
    data_algG_CF_17_filtered$third_base_CF_17[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_CF_17_filtered$third_base_CF_17) == 0, data_algG_CF_17_filtered <- subset(data_algG_CF_17_filtered, select = -c(third_base_CF_17)),data_algG_CF_17_filtered <- data_algG_CF_17_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_filtered_third_variant$third_percentage <- apply(data_algG_CF_17_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_CF_17_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_filtered_third_variant$third_base_CF_17 <- colnames(data_algG_CF_17_filtered_third_variant[, 8:11])[apply(data_algG_CF_17_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_As_CF_17","A", ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_Cs_CF_17","C", ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_Gs_CF_17","G", ifelse(data_algG_CF_17_only_SNP$main_base_CF_17 == "percent_Ts_CF_17","T",z<-0))))

ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,
       data_algG_CF_17_filtered_third_variant$third_base_CF_17 <- 
         ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "percent_As_CF_17","A", 
                ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "percent_Cs_CF_17","C",
                       ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "percent_Gs_CF_17","G",
                              ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "percent_Ts_CF_17","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,
       data_algG_CF_17_filtered_third_variant$third_base_CF_17 <- ifelse(
         data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14,data_algG_CF_17_filtered_third_variant$third_base_CF_17 <- data_algG_CF_17_filtered_third_variant$second_base_CF_17,data_algG_CF_17_filtered_third_variant$third_base_CF_17
       ),
       z<-0
)

ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,
       data_algG_CF_17_filtered_third_variant$third_percentage <- ifelse(
         data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14, data_algG_CF_17_filtered_third_variant$second_percentage,data_algG_CF_17_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_filtered_third_variant <- data_algG_CF_17_filtered_third_variant[,-which(names(data_algG_CF_17_filtered_third_variant) %in% c("second_percentage","second_base_CF_17"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_CF_17_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_17_filtered_third_variant_no_ref)[names(data_algG_CF_17_filtered_third_variant_no_ref) == "second_base_CF_17"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_CF_17_filtered_third_variant_no_ref") == TRUE, names(data_algG_CF_17_filtered_third_variant_no_ref)[names(data_algG_CF_17_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_CF_17",z<-0)




# Define SNP percentage & Add it

data_algG_CF_17_filtered$SNP_base <- 0
data_algG_CF_17_filtered$SNP_percentage_CF_17 <- 0


for(i in data_algG_CF_17_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_CF_17_filtered$`PA14-Koordinate`)
  if(data_algG_CF_17_filtered$main_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] == "A"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_As_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] == "C"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Cs_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] == "G"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Gs_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] == "T"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Ts_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_17_filtered$second_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "A"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_As_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_17_filtered$second_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "C"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Cs_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_17_filtered$second_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "G"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Gs_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_17_filtered$second_base_CF_17[a] == data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "T"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Ts_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "T"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "A"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_As_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "A"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "C"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Cs_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "C"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "G"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Gs_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "G"
  } else if(data_algG_CF_17_filtered$main_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$second_base_CF_17[a] != data_algG_CF_17_filtered$GenBase_PA14[a] & data_algG_CF_17_filtered$main_base_CF_17[a] == "T"){
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- data_algG_CF_17_filtered$percent_Ts_CF_17[a]
    data_algG_CF_17_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_CF_17_filtered$SNP_percentage_CF_17[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_filtered_third_variant$SNP_percentage_CF_17 <- 0,z<-0)


ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,
       data_algG_CF_17_filtered_third_variant$SNP_percentage_CF_17 <- 
         ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "A",
                data_algG_CF_17_filtered_third_variant$percent_As_CF_17, 
                ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "C",
                       data_algG_CF_17_filtered_third_variant$percent_Cs_CF_17,
                       ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "G",
                              data_algG_CF_17_filtered_third_variant$percent_Gs_CF_17,
                              ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "T",
                                     data_algG_CF_17_filtered_third_variant$percent_Ts_CF_17,
                                     ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "A",
                                            data_algG_CF_17_filtered_third_variant$percent_As_CF_17,
                                            ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "C",
                                                   data_algG_CF_17_filtered_third_variant$percent_Cs_CF_17,
                                                   ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "G",
                                                          data_algG_CF_17_filtered_third_variant$percent_Gs_CF_17,
                                                          ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "T",
                                                                 data_algG_CF_17_filtered_third_variant$percent_Ts_CF_17,
                                                                 ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "A",
                                                                        data_algG_CF_17_filtered_third_variant$percent_As_CF_17,
                                                                        ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "C",
                                                                               data_algG_CF_17_filtered_third_variant$percent_Cs_CF_17,
                                                                               ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "G",
                                                                                      data_algG_CF_17_filtered_third_variant$percent_Gs_CF_17,
                                                                                      ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "T",
                                                                                             data_algG_CF_17_filtered_third_variant$percent_Ts_CF_17,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,
       data_algG_CF_17_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "A",
                "A", 
                ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "C",
                       "C",
                       ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "G",
                              "G",
                              ifelse(data_algG_CF_17_filtered_third_variant$main_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "T",
                                     "T",
                                     ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "A",
                                            "A",
                                            ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "C",
                                                   "C",
                                                   ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "G",
                                                          "G",
                                                          ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 == data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$main_base_CF_17 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_CF_17_filtered_third_variant$third_base_CF_17 != data_algG_CF_17_filtered_third_variant$GenBase_PA14 & data_algG_CF_17_filtered_third_variant$third_base_CF_17 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_CF_17_filtered$Isolates_CF_17 <- round((data_algG_CF_17_filtered$SNP_percentage_CF_17*100)/(100/isolates_CF_17))
ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_filtered_third_variant$Isolates_CF_17 <- round((data_algG_CF_17_filtered_third_variant$SNP_percentage_CF_17*100)/(100/isolates_CF_17)),z<-0)
ifelse(exists("data_algG_CF_17_filtered_third_variant_no_ref") == TRUE,data_algG_CF_17_filtered_third_variant_no_ref$Isolates_CF_17 <- round((data_algG_CF_17_filtered_third_variant_no_ref$SNP_percentage_CF_17*100)/(100/isolates_CF_17)),z<-0)



# Clean data

data_algG_CF_17_clean <- data_algG_CF_17_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_clean_third_base <- data_algG_CF_17_filtered_third_variant[,which(names(data_algG_CF_17_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_17","Isolates_CF_17"))],z<-0)
ifelse(exists("data_algG_CF_17_filtered_third_variant") == TRUE,data_algG_CF_17_clean <- rbind(data_algG_CF_17_clean,data_algG_CF_17_clean_third_base),z<-0)


ifelse(exists("data_algG_CF_17_filtered_third_variant_no_ref") == TRUE,data_algG_CF_17_clean_third_base_no_ref <- data_algG_CF_17_filtered_third_variant_no_ref[,which(names(data_algG_CF_17_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_17","Isolates_CF_17"))],z<-0)
ifelse(exists("data_algG_CF_17_filtered_third_variant_no_ref") == TRUE,data_algG_CF_17_clean <- rbind(data_algG_CF_17_clean,data_algG_CF_17_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_CF_17_clean_20_filter <- merge(data_algG_CF_17_clean, data_algG_CF_17_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_CF_17_clean_cleared <- data_algG_CF_17_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_CF_17_clean_20_filter)){
  if(data_algG_CF_17_clean_20_filter$SNP_base[i] == "A" & data_algG_CF_17_clean_20_filter$As_CF_17[i] >= 20){
    data_algG_CF_17_clean_cleared <- rbind(data_algG_CF_17_clean_cleared,data_algG_CF_17_clean_20_filter[i,])
  } else if(data_algG_CF_17_clean_20_filter$SNP_base[i] == "C" & data_algG_CF_17_clean_20_filter$Cs_CF_17[i] >= 20){
    data_algG_CF_17_clean_cleared <- rbind(data_algG_CF_17_clean_cleared,data_algG_CF_17_clean_20_filter[i,])
  } else if(data_algG_CF_17_clean_20_filter$SNP_base[i] == "G" & data_algG_CF_17_clean_20_filter$Gs_CF_17[i] >= 20){
    data_algG_CF_17_clean_cleared <- rbind(data_algG_CF_17_clean_cleared,data_algG_CF_17_clean_20_filter[i,])
  }else if(data_algG_CF_17_clean_20_filter$SNP_base[i] == "T" & data_algG_CF_17_clean_20_filter$Ts_CF_17[i] >= 20){
    data_algG_CF_17_clean_cleared <- rbind(data_algG_CF_17_clean_cleared,data_algG_CF_17_clean_20_filter[i,])
  }
}

data_algG_CF_17_clean <- data_algG_CF_17_clean_cleared[,c(1:5)]

### acute_1

# Input 
data_algG_acute_1 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T18/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_acute_1_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T18/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_acute_1_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T18/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_acute_1 <- data_algG_acute_1[!(data_algG_acute_1$Position %in% primer_correction_algG_acute_1_fwd$Position),]
data_algG_acute_1 <- data_algG_acute_1[!(data_algG_acute_1$Position %in% primer_correction_algG_acute_1_rev$Position),]
data_algG_acute_1 <- rbind(data_algG_acute_1,primer_correction_algG_acute_1_fwd,primer_correction_algG_acute_1_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_acute_1 <- subset(data_algG_acute_1, coverage_acute_1 >199)
data_algG_acute_1_5UTR <- subset(data_algG_acute_1_5UTR, coverage_acute_1_5UTR >199)
data_algG_acute_1_3UTR <- subset(data_algG_acute_1_3UTR, coverage_acute_1_3UTR >199)


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

# Get positions with only SNP

data_algG_acute_1_only_SNP <- subset(data_algG_acute_1, second_percentage <= 0.025) 
data_algG_acute_1_only_SNP$main_base_acute_1 <- colnames(data_algG_acute_1_only_SNP[, 8:11])[apply(data_algG_acute_1_only_SNP[, 8:11],1,which.max)]
data_algG_acute_1_only_SNP$main_base_acute_1 <- ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_As_acute_1","A", ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_Cs_acute_1","C", ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_Gs_acute_1","G", ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_Ts_acute_1","T",z<-0))))
data_algG_acute_1_only_SNP <- subset(data_algG_acute_1_only_SNP,GenBase_PA14 != main_base_acute_1)
ifelse(nrow(data_algG_acute_1_only_SNP) > 0,data_algG_acute_1_only_SNP$Isolates_acute_1 <- isolates_acute_1, z <- 0)


ifelse(nrow(data_algG_acute_1_only_SNP) > 0,data_algG_acute_1_only_SNP_clean <- data_algG_acute_1_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_acute_1_only_SNP) > 0,names(data_algG_acute_1_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_1","Isolates_acute_1"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_acute_1_filtered <- subset(data_algG_acute_1, second_percentage >= 0.025) 
data_algG_acute_1_filtered$main_base_acute_1 <- colnames(data_algG_acute_1_filtered[, 8:11])[apply(data_algG_acute_1_filtered[, 8:11],1,which.max)]
data_algG_acute_1_filtered$second_base_acute_1 <- colnames(data_algG_acute_1_filtered[, 8:11])[apply(data_algG_acute_1_filtered[, 8:11], 1, maxn(2))]
data_algG_acute_1_filtered$third_base_acute_1 <- ifelse(apply(data_algG_acute_1_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_acute_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_1_filtered$`PA14-Koordinate`)
  if(data_algG_acute_1_filtered$main_base_acute_1[a] == "percent_As_acute_1"){
    data_algG_acute_1_filtered$main_base_acute_1[a] <- "A"                     
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] == "percent_Cs_acute_1"){
    data_algG_acute_1_filtered$main_base_acute_1[a] <- "C"                     
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] == "percent_Gs_acute_1"){
    data_algG_acute_1_filtered$main_base_acute_1[a] <- "G"                     
  } else{
    data_algG_acute_1_filtered$main_base_acute_1[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_acute_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_1_filtered$`PA14-Koordinate`)
  if(data_algG_acute_1_filtered$second_base_acute_1[a] == "percent_As_acute_1"){
    data_algG_acute_1_filtered$second_base_acute_1[a] <- "A"                     
  } else if(data_algG_acute_1_filtered$second_base_acute_1[a] == "percent_Cs_acute_1"){
    data_algG_acute_1_filtered$second_base_acute_1[a] <- "C"                     
  } else if(data_algG_acute_1_filtered$second_base_acute_1[a] == "percent_Gs_acute_1"){
    data_algG_acute_1_filtered$second_base_acute_1[a] <- "G"                     
  } else{
    data_algG_acute_1_filtered$second_base_acute_1[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_acute_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_1_filtered$`PA14-Koordinate`)
  if(data_algG_acute_1_filtered$main_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$third_base_acute_1[a] == 0){
    data_algG_acute_1_filtered$third_base_acute_1[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_acute_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_1_filtered$`PA14-Koordinate`)
  if(data_algG_acute_1_filtered$third_base_acute_1[a] == 1){
    data_algG_acute_1_filtered_third_variant <- data_algG_acute_1_filtered[a,]
    data_algG_acute_1_filtered$third_base_acute_1[a] <- 0
  }
  rm(a)
}


for(i in data_algG_acute_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_1_filtered$`PA14-Koordinate`)
  if(data_algG_acute_1_filtered$third_base_acute_1[a] == 2){
    data_algG_acute_1_filtered_third_variant_no_ref <- data_algG_acute_1_filtered[a,]
    data_algG_acute_1_filtered$third_base_acute_1[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_acute_1_filtered$third_base_acute_1) == 0, data_algG_acute_1_filtered <- subset(data_algG_acute_1_filtered, select = -c(third_base_acute_1)),data_algG_acute_1_filtered <- data_algG_acute_1_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_filtered_third_variant$third_percentage <- apply(data_algG_acute_1_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_acute_1_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_filtered_third_variant$third_base_acute_1 <- colnames(data_algG_acute_1_filtered_third_variant[, 8:11])[apply(data_algG_acute_1_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_As_acute_1","A", ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_Cs_acute_1","C", ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_Gs_acute_1","G", ifelse(data_algG_acute_1_only_SNP$main_base_acute_1 == "percent_Ts_acute_1","T",z<-0))))

ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,
       data_algG_acute_1_filtered_third_variant$third_base_acute_1 <- 
         ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "percent_As_acute_1","A", 
                ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "percent_Cs_acute_1","C",
                       ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "percent_Gs_acute_1","G",
                              ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "percent_Ts_acute_1","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,
       data_algG_acute_1_filtered_third_variant$third_base_acute_1 <- ifelse(
         data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14,data_algG_acute_1_filtered_third_variant$third_base_acute_1 <- data_algG_acute_1_filtered_third_variant$second_base_acute_1,data_algG_acute_1_filtered_third_variant$third_base_acute_1
       ),
       z<-0
)

ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,
       data_algG_acute_1_filtered_third_variant$third_percentage <- ifelse(
         data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14, data_algG_acute_1_filtered_third_variant$second_percentage,data_algG_acute_1_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_filtered_third_variant <- data_algG_acute_1_filtered_third_variant[,-which(names(data_algG_acute_1_filtered_third_variant) %in% c("second_percentage","second_base_acute_1"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_acute_1_filtered_third_variant_no_ref") == TRUE, names(data_algG_acute_1_filtered_third_variant_no_ref)[names(data_algG_acute_1_filtered_third_variant_no_ref) == "second_base_acute_1"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_acute_1_filtered_third_variant_no_ref") == TRUE, names(data_algG_acute_1_filtered_third_variant_no_ref)[names(data_algG_acute_1_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_acute_1",z<-0)




# Define SNP percentage & Add it

data_algG_acute_1_filtered$SNP_base <- 0
data_algG_acute_1_filtered$SNP_percentage_acute_1 <- 0


for(i in data_algG_acute_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_1_filtered$`PA14-Koordinate`)
  if(data_algG_acute_1_filtered$main_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] == "A"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_As_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] == "C"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Cs_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] == "G"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Gs_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] == "T"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Ts_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "T"
  } else if(data_algG_acute_1_filtered$second_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "A"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_As_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_1_filtered$second_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "C"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Cs_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_1_filtered$second_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "G"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Gs_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_1_filtered$second_base_acute_1[a] == data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "T"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Ts_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "T"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "A"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_As_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "C"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Cs_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "G"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Gs_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_1_filtered$main_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$second_base_acute_1[a] != data_algG_acute_1_filtered$GenBase_PA14[a] & data_algG_acute_1_filtered$main_base_acute_1[a] == "T"){
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- data_algG_acute_1_filtered$percent_Ts_acute_1[a]
    data_algG_acute_1_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_acute_1_filtered$SNP_percentage_acute_1[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_filtered_third_variant$SNP_percentage_acute_1 <- 0,z<-0)


ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,
       data_algG_acute_1_filtered_third_variant$SNP_percentage_acute_1 <- 
         ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "A",
                data_algG_acute_1_filtered_third_variant$percent_As_acute_1, 
                ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "C",
                       data_algG_acute_1_filtered_third_variant$percent_Cs_acute_1,
                       ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "G",
                              data_algG_acute_1_filtered_third_variant$percent_Gs_acute_1,
                              ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "T",
                                     data_algG_acute_1_filtered_third_variant$percent_Ts_acute_1,
                                     ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "A",
                                            data_algG_acute_1_filtered_third_variant$percent_As_acute_1,
                                            ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "C",
                                                   data_algG_acute_1_filtered_third_variant$percent_Cs_acute_1,
                                                   ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "G",
                                                          data_algG_acute_1_filtered_third_variant$percent_Gs_acute_1,
                                                          ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "T",
                                                                 data_algG_acute_1_filtered_third_variant$percent_Ts_acute_1,
                                                                 ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "A",
                                                                        data_algG_acute_1_filtered_third_variant$percent_As_acute_1,
                                                                        ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "C",
                                                                               data_algG_acute_1_filtered_third_variant$percent_Cs_acute_1,
                                                                               ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "G",
                                                                                      data_algG_acute_1_filtered_third_variant$percent_Gs_acute_1,
                                                                                      ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "T",
                                                                                             data_algG_acute_1_filtered_third_variant$percent_Ts_acute_1,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,
       data_algG_acute_1_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "A",
                "A", 
                ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "C",
                       "C",
                       ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "G",
                              "G",
                              ifelse(data_algG_acute_1_filtered_third_variant$main_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "T",
                                     "T",
                                     ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "A",
                                            "A",
                                            ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "C",
                                                   "C",
                                                   ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "G",
                                                          "G",
                                                          ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 == data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$main_base_acute_1 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_acute_1_filtered_third_variant$third_base_acute_1 != data_algG_acute_1_filtered_third_variant$GenBase_PA14 & data_algG_acute_1_filtered_third_variant$third_base_acute_1 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_acute_1_filtered$Isolates_acute_1 <- round((data_algG_acute_1_filtered$SNP_percentage_acute_1*100)/(100/isolates_acute_1))
ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_filtered_third_variant$Isolates_acute_1 <- round((data_algG_acute_1_filtered_third_variant$SNP_percentage_acute_1*100)/(100/isolates_acute_1)),z<-0)
ifelse(exists("data_algG_acute_1_filtered_third_variant_no_ref") == TRUE,data_algG_acute_1_filtered_third_variant_no_ref$Isolates_acute_1 <- round((data_algG_acute_1_filtered_third_variant_no_ref$SNP_percentage_acute_1*100)/(100/isolates_acute_1)),z<-0)



# Clean data

data_algG_acute_1_clean <- data_algG_acute_1_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_clean_third_base <- data_algG_acute_1_filtered_third_variant[,which(names(data_algG_acute_1_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_1","Isolates_acute_1"))],z<-0)
ifelse(exists("data_algG_acute_1_filtered_third_variant") == TRUE,data_algG_acute_1_clean <- rbind(data_algG_acute_1_clean,data_algG_acute_1_clean_third_base),z<-0)


ifelse(exists("data_algG_acute_1_filtered_third_variant_no_ref") == TRUE,data_algG_acute_1_clean_third_base_no_ref <- data_algG_acute_1_filtered_third_variant_no_ref[,which(names(data_algG_acute_1_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_1","Isolates_acute_1"))],z<-0)
ifelse(exists("data_algG_acute_1_filtered_third_variant_no_ref") == TRUE,data_algG_acute_1_clean <- rbind(data_algG_acute_1_clean,data_algG_acute_1_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_acute_1_clean_20_filter <- merge(data_algG_acute_1_clean, data_algG_acute_1_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_acute_1_clean_cleared <- data_algG_acute_1_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_acute_1_clean_20_filter)){
  if(data_algG_acute_1_clean_20_filter$SNP_base[i] == "A" & data_algG_acute_1_clean_20_filter$As_acute_1[i] >= 20){
    data_algG_acute_1_clean_cleared <- rbind(data_algG_acute_1_clean_cleared,data_algG_acute_1_clean_20_filter[i,])
  } else if(data_algG_acute_1_clean_20_filter$SNP_base[i] == "C" & data_algG_acute_1_clean_20_filter$Cs_acute_1[i] >= 20){
    data_algG_acute_1_clean_cleared <- rbind(data_algG_acute_1_clean_cleared,data_algG_acute_1_clean_20_filter[i,])
  } else if(data_algG_acute_1_clean_20_filter$SNP_base[i] == "G" & data_algG_acute_1_clean_20_filter$Gs_acute_1[i] >= 20){
    data_algG_acute_1_clean_cleared <- rbind(data_algG_acute_1_clean_cleared,data_algG_acute_1_clean_20_filter[i,])
  }else if(data_algG_acute_1_clean_20_filter$SNP_base[i] == "T" & data_algG_acute_1_clean_20_filter$Ts_acute_1[i] >= 20){
    data_algG_acute_1_clean_cleared <- rbind(data_algG_acute_1_clean_cleared,data_algG_acute_1_clean_20_filter[i,])
  }
}

data_algG_acute_1_clean <- data_algG_acute_1_clean_cleared[,c(1:5)]


### acute_2

# Input 
data_algG_acute_2 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T19/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_acute_2_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T19/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_acute_2_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T19/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_acute_2 <- data_algG_acute_2[!(data_algG_acute_2$Position %in% primer_correction_algG_acute_2_fwd$Position),]
data_algG_acute_2 <- data_algG_acute_2[!(data_algG_acute_2$Position %in% primer_correction_algG_acute_2_rev$Position),]
data_algG_acute_2 <- rbind(data_algG_acute_2,primer_correction_algG_acute_2_fwd,primer_correction_algG_acute_2_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_acute_2 <- subset(data_algG_acute_2, coverage_acute_2 >199)
data_algG_acute_2_5UTR <- subset(data_algG_acute_2_5UTR, coverage_acute_2_5UTR >199)
data_algG_acute_2_3UTR <- subset(data_algG_acute_2_3UTR, coverage_acute_2_3UTR >199)

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

# Get positions with only SNP

data_algG_acute_2_only_SNP <- subset(data_algG_acute_2, second_percentage <= 0.025) 
data_algG_acute_2_only_SNP$main_base_acute_2 <- colnames(data_algG_acute_2_only_SNP[, 8:11])[apply(data_algG_acute_2_only_SNP[, 8:11],1,which.max)]
data_algG_acute_2_only_SNP$main_base_acute_2 <- ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_As_acute_2","A", ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_Cs_acute_2","C", ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_Gs_acute_2","G", ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_Ts_acute_2","T",z<-0))))
data_algG_acute_2_only_SNP <- subset(data_algG_acute_2_only_SNP,GenBase_PA14 != main_base_acute_2)
ifelse(nrow(data_algG_acute_2_only_SNP) > 0,data_algG_acute_2_only_SNP$Isolates_acute_2 <- isolates_acute_2, z <- 0)


ifelse(nrow(data_algG_acute_2_only_SNP) > 0,data_algG_acute_2_only_SNP_clean <- data_algG_acute_2_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_acute_2_only_SNP) > 0,names(data_algG_acute_2_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_2","Isolates_acute_2"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_acute_2_filtered <- subset(data_algG_acute_2, second_percentage >= 0.025) 
data_algG_acute_2_filtered$main_base_acute_2 <- colnames(data_algG_acute_2_filtered[, 8:11])[apply(data_algG_acute_2_filtered[, 8:11],1,which.max)]
data_algG_acute_2_filtered$second_base_acute_2 <- colnames(data_algG_acute_2_filtered[, 8:11])[apply(data_algG_acute_2_filtered[, 8:11], 1, maxn(2))]
data_algG_acute_2_filtered$third_base_acute_2 <- ifelse(apply(data_algG_acute_2_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_acute_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_2_filtered$`PA14-Koordinate`)
  if(data_algG_acute_2_filtered$main_base_acute_2[a] == "percent_As_acute_2"){
    data_algG_acute_2_filtered$main_base_acute_2[a] <- "A"                     
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] == "percent_Cs_acute_2"){
    data_algG_acute_2_filtered$main_base_acute_2[a] <- "C"                     
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] == "percent_Gs_acute_2"){
    data_algG_acute_2_filtered$main_base_acute_2[a] <- "G"                     
  } else{
    data_algG_acute_2_filtered$main_base_acute_2[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_acute_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_2_filtered$`PA14-Koordinate`)
  if(data_algG_acute_2_filtered$second_base_acute_2[a] == "percent_As_acute_2"){
    data_algG_acute_2_filtered$second_base_acute_2[a] <- "A"                     
  } else if(data_algG_acute_2_filtered$second_base_acute_2[a] == "percent_Cs_acute_2"){
    data_algG_acute_2_filtered$second_base_acute_2[a] <- "C"                     
  } else if(data_algG_acute_2_filtered$second_base_acute_2[a] == "percent_Gs_acute_2"){
    data_algG_acute_2_filtered$second_base_acute_2[a] <- "G"                     
  } else{
    data_algG_acute_2_filtered$second_base_acute_2[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_acute_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_2_filtered$`PA14-Koordinate`)
  if(data_algG_acute_2_filtered$main_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$third_base_acute_2[a] == 0){
    data_algG_acute_2_filtered$third_base_acute_2[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_acute_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_2_filtered$`PA14-Koordinate`)
  if(data_algG_acute_2_filtered$third_base_acute_2[a] == 1){
    data_algG_acute_2_filtered_third_variant <- data_algG_acute_2_filtered[a,]
    data_algG_acute_2_filtered$third_base_acute_2[a] <- 0
  }
  rm(a)
}


for(i in data_algG_acute_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_2_filtered$`PA14-Koordinate`)
  if(data_algG_acute_2_filtered$third_base_acute_2[a] == 2){
    data_algG_acute_2_filtered_third_variant_no_ref <- data_algG_acute_2_filtered[a,]
    data_algG_acute_2_filtered$third_base_acute_2[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_acute_2_filtered$third_base_acute_2) == 0, data_algG_acute_2_filtered <- subset(data_algG_acute_2_filtered, select = -c(third_base_acute_2)),data_algG_acute_2_filtered <- data_algG_acute_2_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_filtered_third_variant$third_percentage <- apply(data_algG_acute_2_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_acute_2_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_filtered_third_variant$third_base_acute_2 <- colnames(data_algG_acute_2_filtered_third_variant[, 8:11])[apply(data_algG_acute_2_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_As_acute_2","A", ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_Cs_acute_2","C", ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_Gs_acute_2","G", ifelse(data_algG_acute_2_only_SNP$main_base_acute_2 == "percent_Ts_acute_2","T",z<-0))))

ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,
       data_algG_acute_2_filtered_third_variant$third_base_acute_2 <- 
         ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "percent_As_acute_2","A", 
                ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "percent_Cs_acute_2","C",
                       ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "percent_Gs_acute_2","G",
                              ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "percent_Ts_acute_2","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,
       data_algG_acute_2_filtered_third_variant$third_base_acute_2 <- ifelse(
         data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14,data_algG_acute_2_filtered_third_variant$third_base_acute_2 <- data_algG_acute_2_filtered_third_variant$second_base_acute_2,data_algG_acute_2_filtered_third_variant$third_base_acute_2
       ),
       z<-0
)

ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,
       data_algG_acute_2_filtered_third_variant$third_percentage <- ifelse(
         data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14, data_algG_acute_2_filtered_third_variant$second_percentage,data_algG_acute_2_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_filtered_third_variant <- data_algG_acute_2_filtered_third_variant[,-which(names(data_algG_acute_2_filtered_third_variant) %in% c("second_percentage","second_base_acute_2"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_acute_2_filtered_third_variant_no_ref") == TRUE, names(data_algG_acute_2_filtered_third_variant_no_ref)[names(data_algG_acute_2_filtered_third_variant_no_ref) == "second_base_acute_2"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_acute_2_filtered_third_variant_no_ref") == TRUE, names(data_algG_acute_2_filtered_third_variant_no_ref)[names(data_algG_acute_2_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_acute_2",z<-0)




# Define SNP percentage & Add it

data_algG_acute_2_filtered$SNP_base <- 0
data_algG_acute_2_filtered$SNP_percentage_acute_2 <- 0


for(i in data_algG_acute_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_2_filtered$`PA14-Koordinate`)
  if(data_algG_acute_2_filtered$main_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] == "A"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_As_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] == "C"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Cs_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] == "G"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Gs_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] == "T"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Ts_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "T"
  } else if(data_algG_acute_2_filtered$second_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "A"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_As_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_2_filtered$second_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "C"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Cs_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_2_filtered$second_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "G"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Gs_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_2_filtered$second_base_acute_2[a] == data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "T"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Ts_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "T"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "A"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_As_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "C"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Cs_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "G"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Gs_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_2_filtered$main_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$second_base_acute_2[a] != data_algG_acute_2_filtered$GenBase_PA14[a] & data_algG_acute_2_filtered$main_base_acute_2[a] == "T"){
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- data_algG_acute_2_filtered$percent_Ts_acute_2[a]
    data_algG_acute_2_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_acute_2_filtered$SNP_percentage_acute_2[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_filtered_third_variant$SNP_percentage_acute_2 <- 0,z<-0)


ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,
       data_algG_acute_2_filtered_third_variant$SNP_percentage_acute_2 <- 
         ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "A",
                data_algG_acute_2_filtered_third_variant$percent_As_acute_2, 
                ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "C",
                       data_algG_acute_2_filtered_third_variant$percent_Cs_acute_2,
                       ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "G",
                              data_algG_acute_2_filtered_third_variant$percent_Gs_acute_2,
                              ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "T",
                                     data_algG_acute_2_filtered_third_variant$percent_Ts_acute_2,
                                     ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "A",
                                            data_algG_acute_2_filtered_third_variant$percent_As_acute_2,
                                            ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "C",
                                                   data_algG_acute_2_filtered_third_variant$percent_Cs_acute_2,
                                                   ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "G",
                                                          data_algG_acute_2_filtered_third_variant$percent_Gs_acute_2,
                                                          ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "T",
                                                                 data_algG_acute_2_filtered_third_variant$percent_Ts_acute_2,
                                                                 ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "A",
                                                                        data_algG_acute_2_filtered_third_variant$percent_As_acute_2,
                                                                        ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "C",
                                                                               data_algG_acute_2_filtered_third_variant$percent_Cs_acute_2,
                                                                               ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "G",
                                                                                      data_algG_acute_2_filtered_third_variant$percent_Gs_acute_2,
                                                                                      ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "T",
                                                                                             data_algG_acute_2_filtered_third_variant$percent_Ts_acute_2,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,
       data_algG_acute_2_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "A",
                "A", 
                ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "C",
                       "C",
                       ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "G",
                              "G",
                              ifelse(data_algG_acute_2_filtered_third_variant$main_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "T",
                                     "T",
                                     ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "A",
                                            "A",
                                            ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "C",
                                                   "C",
                                                   ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "G",
                                                          "G",
                                                          ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 == data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$main_base_acute_2 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_acute_2_filtered_third_variant$third_base_acute_2 != data_algG_acute_2_filtered_third_variant$GenBase_PA14 & data_algG_acute_2_filtered_third_variant$third_base_acute_2 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_acute_2_filtered$Isolates_acute_2 <- round((data_algG_acute_2_filtered$SNP_percentage_acute_2*100)/(100/isolates_acute_2))
ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_filtered_third_variant$Isolates_acute_2 <- round((data_algG_acute_2_filtered_third_variant$SNP_percentage_acute_2*100)/(100/isolates_acute_2)),z<-0)
ifelse(exists("data_algG_acute_2_filtered_third_variant_no_ref") == TRUE,data_algG_acute_2_filtered_third_variant_no_ref$Isolates_acute_2 <- round((data_algG_acute_2_filtered_third_variant_no_ref$SNP_percentage_acute_2*100)/(100/isolates_acute_2)),z<-0)



# Clean data

data_algG_acute_2_clean <- data_algG_acute_2_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_clean_third_base <- data_algG_acute_2_filtered_third_variant[,which(names(data_algG_acute_2_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_2","Isolates_acute_2"))],z<-0)
ifelse(exists("data_algG_acute_2_filtered_third_variant") == TRUE,data_algG_acute_2_clean <- rbind(data_algG_acute_2_clean,data_algG_acute_2_clean_third_base),z<-0)


ifelse(exists("data_algG_acute_2_filtered_third_variant_no_ref") == TRUE,data_algG_acute_2_clean_third_base_no_ref <- data_algG_acute_2_filtered_third_variant_no_ref[,which(names(data_algG_acute_2_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_2","Isolates_acute_2"))],z<-0)
ifelse(exists("data_algG_acute_2_filtered_third_variant_no_ref") == TRUE,data_algG_acute_2_clean <- rbind(data_algG_acute_2_clean,data_algG_acute_2_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_acute_2_clean_20_filter <- merge(data_algG_acute_2_clean, data_algG_acute_2_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_acute_2_clean_cleared <- data_algG_acute_2_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_acute_2_clean_20_filter)){
  if(data_algG_acute_2_clean_20_filter$SNP_base[i] == "A" & data_algG_acute_2_clean_20_filter$As_acute_2[i] >= 20){
    data_algG_acute_2_clean_cleared <- rbind(data_algG_acute_2_clean_cleared,data_algG_acute_2_clean_20_filter[i,])
  } else if(data_algG_acute_2_clean_20_filter$SNP_base[i] == "C" & data_algG_acute_2_clean_20_filter$Cs_acute_2[i] >= 20){
    data_algG_acute_2_clean_cleared <- rbind(data_algG_acute_2_clean_cleared,data_algG_acute_2_clean_20_filter[i,])
  } else if(data_algG_acute_2_clean_20_filter$SNP_base[i] == "G" & data_algG_acute_2_clean_20_filter$Gs_acute_2[i] >= 20){
    data_algG_acute_2_clean_cleared <- rbind(data_algG_acute_2_clean_cleared,data_algG_acute_2_clean_20_filter[i,])
  }else if(data_algG_acute_2_clean_20_filter$SNP_base[i] == "T" & data_algG_acute_2_clean_20_filter$Ts_acute_2[i] >= 20){
    data_algG_acute_2_clean_cleared <- rbind(data_algG_acute_2_clean_cleared,data_algG_acute_2_clean_20_filter[i,])
  }
}

data_algG_acute_2_clean <- data_algG_acute_2_clean_cleared[,c(1:5)]


### acute_3

# Input 
data_algG_acute_3 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T20/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_acute_3_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T20/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_acute_3_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T20/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_acute_3 <- data_algG_acute_3[!(data_algG_acute_3$Position %in% primer_correction_algG_acute_3_fwd$Position),]
data_algG_acute_3 <- data_algG_acute_3[!(data_algG_acute_3$Position %in% primer_correction_algG_acute_3_rev$Position),]
data_algG_acute_3 <- rbind(data_algG_acute_3,primer_correction_algG_acute_3_fwd,primer_correction_algG_acute_3_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_acute_3 <- subset(data_algG_acute_3, coverage_acute_3 >199)
data_algG_acute_3_5UTR <- subset(data_algG_acute_3_5UTR, coverage_acute_3_5UTR >199)
data_algG_acute_3_3UTR <- subset(data_algG_acute_3_3UTR, coverage_acute_3_3UTR >199)

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

# Get positions with only SNP

data_algG_acute_3_only_SNP <- subset(data_algG_acute_3, second_percentage <= 0.025) 
data_algG_acute_3_only_SNP$main_base_acute_3 <- colnames(data_algG_acute_3_only_SNP[, 8:11])[apply(data_algG_acute_3_only_SNP[, 8:11],1,which.max)]
data_algG_acute_3_only_SNP$main_base_acute_3 <- ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_As_acute_3","A", ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_Cs_acute_3","C", ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_Gs_acute_3","G", ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_Ts_acute_3","T",z<-0))))
data_algG_acute_3_only_SNP <- subset(data_algG_acute_3_only_SNP,GenBase_PA14 != main_base_acute_3)
ifelse(nrow(data_algG_acute_3_only_SNP) > 0,data_algG_acute_3_only_SNP$Isolates_acute_3 <- isolates_acute_3, z <- 0)


ifelse(nrow(data_algG_acute_3_only_SNP) > 0,data_algG_acute_3_only_SNP_clean <- data_algG_acute_3_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_acute_3_only_SNP) > 0,names(data_algG_acute_3_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_3","Isolates_acute_3"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_acute_3_filtered <- subset(data_algG_acute_3, second_percentage >= 0.025) 
data_algG_acute_3_filtered$main_base_acute_3 <- colnames(data_algG_acute_3_filtered[, 8:11])[apply(data_algG_acute_3_filtered[, 8:11],1,which.max)]
data_algG_acute_3_filtered$second_base_acute_3 <- colnames(data_algG_acute_3_filtered[, 8:11])[apply(data_algG_acute_3_filtered[, 8:11], 1, maxn(2))]
data_algG_acute_3_filtered$third_base_acute_3 <- ifelse(apply(data_algG_acute_3_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_acute_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_3_filtered$`PA14-Koordinate`)
  if(data_algG_acute_3_filtered$main_base_acute_3[a] == "percent_As_acute_3"){
    data_algG_acute_3_filtered$main_base_acute_3[a] <- "A"                     
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] == "percent_Cs_acute_3"){
    data_algG_acute_3_filtered$main_base_acute_3[a] <- "C"                     
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] == "percent_Gs_acute_3"){
    data_algG_acute_3_filtered$main_base_acute_3[a] <- "G"                     
  } else{
    data_algG_acute_3_filtered$main_base_acute_3[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_acute_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_3_filtered$`PA14-Koordinate`)
  if(data_algG_acute_3_filtered$second_base_acute_3[a] == "percent_As_acute_3"){
    data_algG_acute_3_filtered$second_base_acute_3[a] <- "A"                     
  } else if(data_algG_acute_3_filtered$second_base_acute_3[a] == "percent_Cs_acute_3"){
    data_algG_acute_3_filtered$second_base_acute_3[a] <- "C"                     
  } else if(data_algG_acute_3_filtered$second_base_acute_3[a] == "percent_Gs_acute_3"){
    data_algG_acute_3_filtered$second_base_acute_3[a] <- "G"                     
  } else{
    data_algG_acute_3_filtered$second_base_acute_3[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_acute_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_3_filtered$`PA14-Koordinate`)
  if(data_algG_acute_3_filtered$main_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$third_base_acute_3[a] == 0){
    data_algG_acute_3_filtered$third_base_acute_3[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_acute_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_3_filtered$`PA14-Koordinate`)
  if(data_algG_acute_3_filtered$third_base_acute_3[a] == 1){
    data_algG_acute_3_filtered_third_variant <- data_algG_acute_3_filtered[a,]
    data_algG_acute_3_filtered$third_base_acute_3[a] <- 0
  }
  rm(a)
}


for(i in data_algG_acute_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_3_filtered$`PA14-Koordinate`)
  if(data_algG_acute_3_filtered$third_base_acute_3[a] == 2){
    data_algG_acute_3_filtered_third_variant_no_ref <- data_algG_acute_3_filtered[a,]
    data_algG_acute_3_filtered$third_base_acute_3[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_acute_3_filtered$third_base_acute_3) == 0, data_algG_acute_3_filtered <- subset(data_algG_acute_3_filtered, select = -c(third_base_acute_3)),data_algG_acute_3_filtered <- data_algG_acute_3_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_filtered_third_variant$third_percentage <- apply(data_algG_acute_3_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_acute_3_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_filtered_third_variant$third_base_acute_3 <- colnames(data_algG_acute_3_filtered_third_variant[, 8:11])[apply(data_algG_acute_3_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_As_acute_3","A", ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_Cs_acute_3","C", ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_Gs_acute_3","G", ifelse(data_algG_acute_3_only_SNP$main_base_acute_3 == "percent_Ts_acute_3","T",z<-0))))

ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,
       data_algG_acute_3_filtered_third_variant$third_base_acute_3 <- 
         ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "percent_As_acute_3","A", 
                ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "percent_Cs_acute_3","C",
                       ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "percent_Gs_acute_3","G",
                              ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "percent_Ts_acute_3","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,
       data_algG_acute_3_filtered_third_variant$third_base_acute_3 <- ifelse(
         data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14,data_algG_acute_3_filtered_third_variant$third_base_acute_3 <- data_algG_acute_3_filtered_third_variant$second_base_acute_3,data_algG_acute_3_filtered_third_variant$third_base_acute_3
       ),
       z<-0
)

ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,
       data_algG_acute_3_filtered_third_variant$third_percentage <- ifelse(
         data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14, data_algG_acute_3_filtered_third_variant$second_percentage,data_algG_acute_3_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_filtered_third_variant <- data_algG_acute_3_filtered_third_variant[,-which(names(data_algG_acute_3_filtered_third_variant) %in% c("second_percentage","second_base_acute_3"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_acute_3_filtered_third_variant_no_ref") == TRUE, names(data_algG_acute_3_filtered_third_variant_no_ref)[names(data_algG_acute_3_filtered_third_variant_no_ref) == "second_base_acute_3"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_acute_3_filtered_third_variant_no_ref") == TRUE, names(data_algG_acute_3_filtered_third_variant_no_ref)[names(data_algG_acute_3_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_acute_3",z<-0)




# Define SNP percentage & Add it

data_algG_acute_3_filtered$SNP_base <- 0
data_algG_acute_3_filtered$SNP_percentage_acute_3 <- 0


for(i in data_algG_acute_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_acute_3_filtered$`PA14-Koordinate`)
  if(data_algG_acute_3_filtered$main_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] == "A"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_As_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] == "C"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Cs_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] == "G"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Gs_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] == "T"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Ts_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "T"
  } else if(data_algG_acute_3_filtered$second_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "A"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_As_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_3_filtered$second_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "C"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Cs_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_3_filtered$second_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "G"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Gs_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_3_filtered$second_base_acute_3[a] == data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "T"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Ts_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "T"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "A"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_As_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "C"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Cs_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "G"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Gs_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_acute_3_filtered$main_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$second_base_acute_3[a] != data_algG_acute_3_filtered$GenBase_PA14[a] & data_algG_acute_3_filtered$main_base_acute_3[a] == "T"){
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- data_algG_acute_3_filtered$percent_Ts_acute_3[a]
    data_algG_acute_3_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_acute_3_filtered$SNP_percentage_acute_3[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_filtered_third_variant$SNP_percentage_acute_3 <- 0,z<-0)


ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,
       data_algG_acute_3_filtered_third_variant$SNP_percentage_acute_3 <- 
         ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "A",
                data_algG_acute_3_filtered_third_variant$percent_As_acute_3, 
                ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "C",
                       data_algG_acute_3_filtered_third_variant$percent_Cs_acute_3,
                       ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "G",
                              data_algG_acute_3_filtered_third_variant$percent_Gs_acute_3,
                              ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "T",
                                     data_algG_acute_3_filtered_third_variant$percent_Ts_acute_3,
                                     ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "A",
                                            data_algG_acute_3_filtered_third_variant$percent_As_acute_3,
                                            ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "C",
                                                   data_algG_acute_3_filtered_third_variant$percent_Cs_acute_3,
                                                   ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "G",
                                                          data_algG_acute_3_filtered_third_variant$percent_Gs_acute_3,
                                                          ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "T",
                                                                 data_algG_acute_3_filtered_third_variant$percent_Ts_acute_3,
                                                                 ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "A",
                                                                        data_algG_acute_3_filtered_third_variant$percent_As_acute_3,
                                                                        ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "C",
                                                                               data_algG_acute_3_filtered_third_variant$percent_Cs_acute_3,
                                                                               ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "G",
                                                                                      data_algG_acute_3_filtered_third_variant$percent_Gs_acute_3,
                                                                                      ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "T",
                                                                                             data_algG_acute_3_filtered_third_variant$percent_Ts_acute_3,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,
       data_algG_acute_3_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "A",
                "A", 
                ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "C",
                       "C",
                       ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "G",
                              "G",
                              ifelse(data_algG_acute_3_filtered_third_variant$main_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "T",
                                     "T",
                                     ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "A",
                                            "A",
                                            ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "C",
                                                   "C",
                                                   ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "G",
                                                          "G",
                                                          ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 == data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$main_base_acute_3 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_acute_3_filtered_third_variant$third_base_acute_3 != data_algG_acute_3_filtered_third_variant$GenBase_PA14 & data_algG_acute_3_filtered_third_variant$third_base_acute_3 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_acute_3_filtered$Isolates_acute_3 <- round((data_algG_acute_3_filtered$SNP_percentage_acute_3*100)/(100/isolates_acute_3))
ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_filtered_third_variant$Isolates_acute_3 <- round((data_algG_acute_3_filtered_third_variant$SNP_percentage_acute_3*100)/(100/isolates_acute_3)),z<-0)
ifelse(exists("data_algG_acute_3_filtered_third_variant_no_ref") == TRUE,data_algG_acute_3_filtered_third_variant_no_ref$Isolates_acute_3 <- round((data_algG_acute_3_filtered_third_variant_no_ref$SNP_percentage_acute_3*100)/(100/isolates_acute_3)),z<-0)



# Clean data

data_algG_acute_3_clean <- data_algG_acute_3_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_clean_third_base <- data_algG_acute_3_filtered_third_variant[,which(names(data_algG_acute_3_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_3","Isolates_acute_3"))],z<-0)
ifelse(exists("data_algG_acute_3_filtered_third_variant") == TRUE,data_algG_acute_3_clean <- rbind(data_algG_acute_3_clean,data_algG_acute_3_clean_third_base),z<-0)


ifelse(exists("data_algG_acute_3_filtered_third_variant_no_ref") == TRUE,data_algG_acute_3_clean_third_base_no_ref <- data_algG_acute_3_filtered_third_variant_no_ref[,which(names(data_algG_acute_3_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_3","Isolates_acute_3"))],z<-0)
ifelse(exists("data_algG_acute_3_filtered_third_variant_no_ref") == TRUE,data_algG_acute_3_clean <- rbind(data_algG_acute_3_clean,data_algG_acute_3_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_acute_3_clean_20_filter <- merge(data_algG_acute_3_clean, data_algG_acute_3_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_acute_3_clean_cleared <- data_algG_acute_3_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_acute_3_clean_20_filter)){
  if(data_algG_acute_3_clean_20_filter$SNP_base[i] == "A" & data_algG_acute_3_clean_20_filter$As_acute_3[i] >= 20){
    data_algG_acute_3_clean_cleared <- rbind(data_algG_acute_3_clean_cleared,data_algG_acute_3_clean_20_filter[i,])
  } else if(data_algG_acute_3_clean_20_filter$SNP_base[i] == "C" & data_algG_acute_3_clean_20_filter$Cs_acute_3[i] >= 20){
    data_algG_acute_3_clean_cleared <- rbind(data_algG_acute_3_clean_cleared,data_algG_acute_3_clean_20_filter[i,])
  } else if(data_algG_acute_3_clean_20_filter$SNP_base[i] == "G" & data_algG_acute_3_clean_20_filter$Gs_acute_3[i] >= 20){
    data_algG_acute_3_clean_cleared <- rbind(data_algG_acute_3_clean_cleared,data_algG_acute_3_clean_20_filter[i,])
  }else if(data_algG_acute_3_clean_20_filter$SNP_base[i] == "T" & data_algG_acute_3_clean_20_filter$Ts_acute_3[i] >= 20){
    data_algG_acute_3_clean_cleared <- rbind(data_algG_acute_3_clean_cleared,data_algG_acute_3_clean_20_filter[i,])
  }
}

data_algG_acute_3_clean <- data_algG_acute_3_clean_cleared[,c(1:5)]

### env_1

# Input 
data_algG_env_1 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T21/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_env_1_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T21/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_env_1_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T21/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_env_1 <- data_algG_env_1[!(data_algG_env_1$Position %in% primer_correction_algG_env_1_fwd$Position),]
data_algG_env_1 <- data_algG_env_1[!(data_algG_env_1$Position %in% primer_correction_algG_env_1_rev$Position),]
data_algG_env_1 <- rbind(data_algG_env_1,primer_correction_algG_env_1_fwd,primer_correction_algG_env_1_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_env_1 <- subset(data_algG_env_1, coverage_env_1 >199)
data_algG_env_1_5UTR <- subset(data_algG_env_1_5UTR, coverage_env_1_5UTR >199)
data_algG_env_1_3UTR <- subset(data_algG_env_1_3UTR, coverage_env_1_3UTR >199)


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

# Get positions with only SNP

data_algG_env_1_only_SNP <- subset(data_algG_env_1, second_percentage <= 0.025) 
data_algG_env_1_only_SNP$main_base_env_1 <- colnames(data_algG_env_1_only_SNP[, 8:11])[apply(data_algG_env_1_only_SNP[, 8:11],1,which.max)]
data_algG_env_1_only_SNP$main_base_env_1 <- ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_As_env_1","A", ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_Cs_env_1","C", ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_Gs_env_1","G", ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_Ts_env_1","T",z<-0))))
data_algG_env_1_only_SNP <- subset(data_algG_env_1_only_SNP,GenBase_PA14 != main_base_env_1)
ifelse(nrow(data_algG_env_1_only_SNP) > 0,data_algG_env_1_only_SNP$Isolates_env_1 <- isolates_env_1, z <- 0)


ifelse(nrow(data_algG_env_1_only_SNP) > 0,data_algG_env_1_only_SNP_clean <- data_algG_env_1_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_env_1_only_SNP) > 0,names(data_algG_env_1_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_1","Isolates_env_1"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_env_1_filtered <- subset(data_algG_env_1, second_percentage >= 0.025) 
data_algG_env_1_filtered$main_base_env_1 <- colnames(data_algG_env_1_filtered[, 8:11])[apply(data_algG_env_1_filtered[, 8:11],1,which.max)]
data_algG_env_1_filtered$second_base_env_1 <- colnames(data_algG_env_1_filtered[, 8:11])[apply(data_algG_env_1_filtered[, 8:11], 1, maxn(2))]
data_algG_env_1_filtered$third_base_env_1 <- ifelse(apply(data_algG_env_1_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_env_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_1_filtered$`PA14-Koordinate`)
  if(data_algG_env_1_filtered$main_base_env_1[a] == "percent_As_env_1"){
    data_algG_env_1_filtered$main_base_env_1[a] <- "A"                     
  } else if(data_algG_env_1_filtered$main_base_env_1[a] == "percent_Cs_env_1"){
    data_algG_env_1_filtered$main_base_env_1[a] <- "C"                     
  } else if(data_algG_env_1_filtered$main_base_env_1[a] == "percent_Gs_env_1"){
    data_algG_env_1_filtered$main_base_env_1[a] <- "G"                     
  } else{
    data_algG_env_1_filtered$main_base_env_1[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_env_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_1_filtered$`PA14-Koordinate`)
  if(data_algG_env_1_filtered$second_base_env_1[a] == "percent_As_env_1"){
    data_algG_env_1_filtered$second_base_env_1[a] <- "A"                     
  } else if(data_algG_env_1_filtered$second_base_env_1[a] == "percent_Cs_env_1"){
    data_algG_env_1_filtered$second_base_env_1[a] <- "C"                     
  } else if(data_algG_env_1_filtered$second_base_env_1[a] == "percent_Gs_env_1"){
    data_algG_env_1_filtered$second_base_env_1[a] <- "G"                     
  } else{
    data_algG_env_1_filtered$second_base_env_1[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_env_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_1_filtered$`PA14-Koordinate`)
  if(data_algG_env_1_filtered$main_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$third_base_env_1[a] == 0){
    data_algG_env_1_filtered$third_base_env_1[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_env_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_1_filtered$`PA14-Koordinate`)
  if(data_algG_env_1_filtered$third_base_env_1[a] == 1){
    data_algG_env_1_filtered_third_variant <- data_algG_env_1_filtered[a,]
    data_algG_env_1_filtered$third_base_env_1[a] <- 0
  }
  rm(a)
}


for(i in data_algG_env_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_1_filtered$`PA14-Koordinate`)
  if(data_algG_env_1_filtered$third_base_env_1[a] == 2){
    data_algG_env_1_filtered_third_variant_no_ref <- data_algG_env_1_filtered[a,]
    data_algG_env_1_filtered$third_base_env_1[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_env_1_filtered$third_base_env_1) == 0, data_algG_env_1_filtered <- subset(data_algG_env_1_filtered, select = -c(third_base_env_1)),data_algG_env_1_filtered <- data_algG_env_1_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_filtered_third_variant$third_percentage <- apply(data_algG_env_1_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_1_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_filtered_third_variant$third_base_env_1 <- colnames(data_algG_env_1_filtered_third_variant[, 8:11])[apply(data_algG_env_1_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_As_env_1","A", ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_Cs_env_1","C", ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_Gs_env_1","G", ifelse(data_algG_env_1_only_SNP$main_base_env_1 == "percent_Ts_env_1","T",z<-0))))

ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,
       data_algG_env_1_filtered_third_variant$third_base_env_1 <- 
         ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == "percent_As_env_1","A", 
                ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == "percent_Cs_env_1","C",
                       ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == "percent_Gs_env_1","G",
                              ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == "percent_Ts_env_1","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,
       data_algG_env_1_filtered_third_variant$third_base_env_1 <- ifelse(
         data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14,data_algG_env_1_filtered_third_variant$third_base_env_1 <- data_algG_env_1_filtered_third_variant$second_base_env_1,data_algG_env_1_filtered_third_variant$third_base_env_1
       ),
       z<-0
)

ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,
       data_algG_env_1_filtered_third_variant$third_percentage <- ifelse(
         data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14, data_algG_env_1_filtered_third_variant$second_percentage,data_algG_env_1_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_filtered_third_variant <- data_algG_env_1_filtered_third_variant[,-which(names(data_algG_env_1_filtered_third_variant) %in% c("second_percentage","second_base_env_1"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_env_1_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_1_filtered_third_variant_no_ref)[names(data_algG_env_1_filtered_third_variant_no_ref) == "second_base_env_1"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_env_1_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_1_filtered_third_variant_no_ref)[names(data_algG_env_1_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_env_1",z<-0)




# Define SNP percentage & Add it

data_algG_env_1_filtered$SNP_base <- 0
data_algG_env_1_filtered$SNP_percentage_env_1 <- 0


for(i in data_algG_env_1_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_1_filtered$`PA14-Koordinate`)
  if(data_algG_env_1_filtered$main_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] == "A"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_As_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] == "C"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Cs_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] == "G"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Gs_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] == "T"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Ts_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_1_filtered$second_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "A"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_As_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_1_filtered$second_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "C"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Cs_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_1_filtered$second_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "G"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Gs_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_1_filtered$second_base_env_1[a] == data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "T"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Ts_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "A"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_As_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "C"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Cs_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "G"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Gs_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_1_filtered$main_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$second_base_env_1[a] != data_algG_env_1_filtered$GenBase_PA14[a] & data_algG_env_1_filtered$main_base_env_1[a] == "T"){
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- data_algG_env_1_filtered$percent_Ts_env_1[a]
    data_algG_env_1_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_env_1_filtered$SNP_percentage_env_1[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_filtered_third_variant$SNP_percentage_env_1 <- 0,z<-0)


ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,
       data_algG_env_1_filtered_third_variant$SNP_percentage_env_1 <- 
         ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "A",
                data_algG_env_1_filtered_third_variant$percent_As_env_1, 
                ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "C",
                       data_algG_env_1_filtered_third_variant$percent_Cs_env_1,
                       ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "G",
                              data_algG_env_1_filtered_third_variant$percent_Gs_env_1,
                              ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "T",
                                     data_algG_env_1_filtered_third_variant$percent_Ts_env_1,
                                     ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "A",
                                            data_algG_env_1_filtered_third_variant$percent_As_env_1,
                                            ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "C",
                                                   data_algG_env_1_filtered_third_variant$percent_Cs_env_1,
                                                   ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "G",
                                                          data_algG_env_1_filtered_third_variant$percent_Gs_env_1,
                                                          ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "T",
                                                                 data_algG_env_1_filtered_third_variant$percent_Ts_env_1,
                                                                 ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "A",
                                                                        data_algG_env_1_filtered_third_variant$percent_As_env_1,
                                                                        ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "C",
                                                                               data_algG_env_1_filtered_third_variant$percent_Cs_env_1,
                                                                               ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "G",
                                                                                      data_algG_env_1_filtered_third_variant$percent_Gs_env_1,
                                                                                      ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "T",
                                                                                             data_algG_env_1_filtered_third_variant$percent_Ts_env_1,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,
       data_algG_env_1_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "A",
                "A", 
                ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "C",
                       "C",
                       ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "G",
                              "G",
                              ifelse(data_algG_env_1_filtered_third_variant$main_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "T",
                                     "T",
                                     ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "A",
                                            "A",
                                            ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "C",
                                                   "C",
                                                   ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "G",
                                                          "G",
                                                          ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 == data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$main_base_env_1 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_env_1_filtered_third_variant$third_base_env_1 != data_algG_env_1_filtered_third_variant$GenBase_PA14 & data_algG_env_1_filtered_third_variant$third_base_env_1 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_env_1_filtered$Isolates_env_1 <- round((data_algG_env_1_filtered$SNP_percentage_env_1*100)/(100/isolates_env_1))
ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_filtered_third_variant$Isolates_env_1 <- round((data_algG_env_1_filtered_third_variant$SNP_percentage_env_1*100)/(100/isolates_env_1)),z<-0)
ifelse(exists("data_algG_env_1_filtered_third_variant_no_ref") == TRUE,data_algG_env_1_filtered_third_variant_no_ref$Isolates_env_1 <- round((data_algG_env_1_filtered_third_variant_no_ref$SNP_percentage_env_1*100)/(100/isolates_env_1)),z<-0)



# Clean data

data_algG_env_1_clean <- data_algG_env_1_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_clean_third_base <- data_algG_env_1_filtered_third_variant[,which(names(data_algG_env_1_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_1","Isolates_env_1"))],z<-0)
ifelse(exists("data_algG_env_1_filtered_third_variant") == TRUE,data_algG_env_1_clean <- rbind(data_algG_env_1_clean,data_algG_env_1_clean_third_base),z<-0)


ifelse(exists("data_algG_env_1_filtered_third_variant_no_ref") == TRUE,data_algG_env_1_clean_third_base_no_ref <- data_algG_env_1_filtered_third_variant_no_ref[,which(names(data_algG_env_1_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_1","Isolates_env_1"))],z<-0)
ifelse(exists("data_algG_env_1_filtered_third_variant_no_ref") == TRUE,data_algG_env_1_clean <- rbind(data_algG_env_1_clean,data_algG_env_1_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_env_1_clean_20_filter <- merge(data_algG_env_1_clean, data_algG_env_1_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_env_1_clean_cleared <- data_algG_env_1_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_env_1_clean_20_filter)){
  if(data_algG_env_1_clean_20_filter$SNP_base[i] == "A" & data_algG_env_1_clean_20_filter$As_env_1[i] >= 20){
    data_algG_env_1_clean_cleared <- rbind(data_algG_env_1_clean_cleared,data_algG_env_1_clean_20_filter[i,])
  } else if(data_algG_env_1_clean_20_filter$SNP_base[i] == "C" & data_algG_env_1_clean_20_filter$Cs_env_1[i] >= 20){
    data_algG_env_1_clean_cleared <- rbind(data_algG_env_1_clean_cleared,data_algG_env_1_clean_20_filter[i,])
  } else if(data_algG_env_1_clean_20_filter$SNP_base[i] == "G" & data_algG_env_1_clean_20_filter$Gs_env_1[i] >= 20){
    data_algG_env_1_clean_cleared <- rbind(data_algG_env_1_clean_cleared,data_algG_env_1_clean_20_filter[i,])
  }else if(data_algG_env_1_clean_20_filter$SNP_base[i] == "T" & data_algG_env_1_clean_20_filter$Ts_env_1[i] >= 20){
    data_algG_env_1_clean_cleared <- rbind(data_algG_env_1_clean_cleared,data_algG_env_1_clean_20_filter[i,])
  }
}

data_algG_env_1_clean <- data_algG_env_1_clean_cleared[,c(1:5)]


### env_2

# Input 
data_algG_env_2 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T22/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_env_2_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T22/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_env_2_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T22/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_env_2 <- data_algG_env_2[!(data_algG_env_2$Position %in% primer_correction_algG_env_2_fwd$Position),]
data_algG_env_2 <- data_algG_env_2[!(data_algG_env_2$Position %in% primer_correction_algG_env_2_rev$Position),]
data_algG_env_2 <- rbind(data_algG_env_2,primer_correction_algG_env_2_fwd,primer_correction_algG_env_2_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_env_2 <- subset(data_algG_env_2, coverage_env_2 >199)
data_algG_env_2_5UTR <- subset(data_algG_env_2_5UTR, coverage_env_2_5UTR >199)
data_algG_env_2_3UTR <- subset(data_algG_env_2_3UTR, coverage_env_2_3UTR >199)

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

# Get positions with only SNP

data_algG_env_2_only_SNP <- subset(data_algG_env_2, second_percentage <= 0.025) 
data_algG_env_2_only_SNP$main_base_env_2 <- colnames(data_algG_env_2_only_SNP[, 8:11])[apply(data_algG_env_2_only_SNP[, 8:11],1,which.max)]
data_algG_env_2_only_SNP$main_base_env_2 <- ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_As_env_2","A", ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_Cs_env_2","C", ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_Gs_env_2","G", ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_Ts_env_2","T",z<-0))))
data_algG_env_2_only_SNP <- subset(data_algG_env_2_only_SNP,GenBase_PA14 != main_base_env_2)
ifelse(nrow(data_algG_env_2_only_SNP) > 0,data_algG_env_2_only_SNP$Isolates_env_2 <- isolates_env_2, z <- 0)


ifelse(nrow(data_algG_env_2_only_SNP) > 0,data_algG_env_2_only_SNP_clean <- data_algG_env_2_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_env_2_only_SNP) > 0,names(data_algG_env_2_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_2","Isolates_env_2"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_env_2_filtered <- subset(data_algG_env_2, second_percentage >= 0.025) 
data_algG_env_2_filtered$main_base_env_2 <- colnames(data_algG_env_2_filtered[, 8:11])[apply(data_algG_env_2_filtered[, 8:11],1,which.max)]
data_algG_env_2_filtered$second_base_env_2 <- colnames(data_algG_env_2_filtered[, 8:11])[apply(data_algG_env_2_filtered[, 8:11], 1, maxn(2))]
data_algG_env_2_filtered$third_base_env_2 <- ifelse(apply(data_algG_env_2_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_env_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_2_filtered$`PA14-Koordinate`)
  if(data_algG_env_2_filtered$main_base_env_2[a] == "percent_As_env_2"){
    data_algG_env_2_filtered$main_base_env_2[a] <- "A"                     
  } else if(data_algG_env_2_filtered$main_base_env_2[a] == "percent_Cs_env_2"){
    data_algG_env_2_filtered$main_base_env_2[a] <- "C"                     
  } else if(data_algG_env_2_filtered$main_base_env_2[a] == "percent_Gs_env_2"){
    data_algG_env_2_filtered$main_base_env_2[a] <- "G"                     
  } else{
    data_algG_env_2_filtered$main_base_env_2[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_env_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_2_filtered$`PA14-Koordinate`)
  if(data_algG_env_2_filtered$second_base_env_2[a] == "percent_As_env_2"){
    data_algG_env_2_filtered$second_base_env_2[a] <- "A"                     
  } else if(data_algG_env_2_filtered$second_base_env_2[a] == "percent_Cs_env_2"){
    data_algG_env_2_filtered$second_base_env_2[a] <- "C"                     
  } else if(data_algG_env_2_filtered$second_base_env_2[a] == "percent_Gs_env_2"){
    data_algG_env_2_filtered$second_base_env_2[a] <- "G"                     
  } else{
    data_algG_env_2_filtered$second_base_env_2[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_env_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_2_filtered$`PA14-Koordinate`)
  if(data_algG_env_2_filtered$main_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$third_base_env_2[a] == 0){
    data_algG_env_2_filtered$third_base_env_2[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_env_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_2_filtered$`PA14-Koordinate`)
  if(data_algG_env_2_filtered$third_base_env_2[a] == 1){
    data_algG_env_2_filtered_third_variant <- data_algG_env_2_filtered[a,]
    data_algG_env_2_filtered$third_base_env_2[a] <- 0
  }
  rm(a)
}


for(i in data_algG_env_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_2_filtered$`PA14-Koordinate`)
  if(data_algG_env_2_filtered$third_base_env_2[a] == 2){
    data_algG_env_2_filtered_third_variant_no_ref <- data_algG_env_2_filtered[a,]
    data_algG_env_2_filtered$third_base_env_2[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_env_2_filtered$third_base_env_2) == 0, data_algG_env_2_filtered <- subset(data_algG_env_2_filtered, select = -c(third_base_env_2)),data_algG_env_2_filtered <- data_algG_env_2_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_filtered_third_variant$third_percentage <- apply(data_algG_env_2_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_2_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_filtered_third_variant$third_base_env_2 <- colnames(data_algG_env_2_filtered_third_variant[, 8:11])[apply(data_algG_env_2_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_As_env_2","A", ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_Cs_env_2","C", ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_Gs_env_2","G", ifelse(data_algG_env_2_only_SNP$main_base_env_2 == "percent_Ts_env_2","T",z<-0))))

ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,
       data_algG_env_2_filtered_third_variant$third_base_env_2 <- 
         ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == "percent_As_env_2","A", 
                ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == "percent_Cs_env_2","C",
                       ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == "percent_Gs_env_2","G",
                              ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == "percent_Ts_env_2","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,
       data_algG_env_2_filtered_third_variant$third_base_env_2 <- ifelse(
         data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14,data_algG_env_2_filtered_third_variant$third_base_env_2 <- data_algG_env_2_filtered_third_variant$second_base_env_2,data_algG_env_2_filtered_third_variant$third_base_env_2
       ),
       z<-0
)

ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,
       data_algG_env_2_filtered_third_variant$third_percentage <- ifelse(
         data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14, data_algG_env_2_filtered_third_variant$second_percentage,data_algG_env_2_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_filtered_third_variant <- data_algG_env_2_filtered_third_variant[,-which(names(data_algG_env_2_filtered_third_variant) %in% c("second_percentage","second_base_env_2"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_env_2_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_2_filtered_third_variant_no_ref)[names(data_algG_env_2_filtered_third_variant_no_ref) == "second_base_env_2"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_env_2_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_2_filtered_third_variant_no_ref)[names(data_algG_env_2_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_env_2",z<-0)




# Define SNP percentage & Add it

data_algG_env_2_filtered$SNP_base <- 0
data_algG_env_2_filtered$SNP_percentage_env_2 <- 0


for(i in data_algG_env_2_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_2_filtered$`PA14-Koordinate`)
  if(data_algG_env_2_filtered$main_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] == "A"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_As_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] == "C"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Cs_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] == "G"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Gs_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] == "T"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Ts_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_2_filtered$second_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "A"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_As_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_2_filtered$second_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "C"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Cs_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_2_filtered$second_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "G"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Gs_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_2_filtered$second_base_env_2[a] == data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "T"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Ts_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "A"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_As_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "C"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Cs_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "G"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Gs_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_2_filtered$main_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$second_base_env_2[a] != data_algG_env_2_filtered$GenBase_PA14[a] & data_algG_env_2_filtered$main_base_env_2[a] == "T"){
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- data_algG_env_2_filtered$percent_Ts_env_2[a]
    data_algG_env_2_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_env_2_filtered$SNP_percentage_env_2[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_filtered_third_variant$SNP_percentage_env_2 <- 0,z<-0)


ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,
       data_algG_env_2_filtered_third_variant$SNP_percentage_env_2 <- 
         ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "A",
                data_algG_env_2_filtered_third_variant$percent_As_env_2, 
                ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "C",
                       data_algG_env_2_filtered_third_variant$percent_Cs_env_2,
                       ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "G",
                              data_algG_env_2_filtered_third_variant$percent_Gs_env_2,
                              ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "T",
                                     data_algG_env_2_filtered_third_variant$percent_Ts_env_2,
                                     ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "A",
                                            data_algG_env_2_filtered_third_variant$percent_As_env_2,
                                            ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "C",
                                                   data_algG_env_2_filtered_third_variant$percent_Cs_env_2,
                                                   ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "G",
                                                          data_algG_env_2_filtered_third_variant$percent_Gs_env_2,
                                                          ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "T",
                                                                 data_algG_env_2_filtered_third_variant$percent_Ts_env_2,
                                                                 ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "A",
                                                                        data_algG_env_2_filtered_third_variant$percent_As_env_2,
                                                                        ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "C",
                                                                               data_algG_env_2_filtered_third_variant$percent_Cs_env_2,
                                                                               ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "G",
                                                                                      data_algG_env_2_filtered_third_variant$percent_Gs_env_2,
                                                                                      ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "T",
                                                                                             data_algG_env_2_filtered_third_variant$percent_Ts_env_2,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,
       data_algG_env_2_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "A",
                "A", 
                ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "C",
                       "C",
                       ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "G",
                              "G",
                              ifelse(data_algG_env_2_filtered_third_variant$main_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "T",
                                     "T",
                                     ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "A",
                                            "A",
                                            ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "C",
                                                   "C",
                                                   ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "G",
                                                          "G",
                                                          ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 == data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$main_base_env_2 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_env_2_filtered_third_variant$third_base_env_2 != data_algG_env_2_filtered_third_variant$GenBase_PA14 & data_algG_env_2_filtered_third_variant$third_base_env_2 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_env_2_filtered$Isolates_env_2 <- round((data_algG_env_2_filtered$SNP_percentage_env_2*100)/(100/isolates_env_2))
ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_filtered_third_variant$Isolates_env_2 <- round((data_algG_env_2_filtered_third_variant$SNP_percentage_env_2*100)/(100/isolates_env_2)),z<-0)
ifelse(exists("data_algG_env_2_filtered_third_variant_no_ref") == TRUE,data_algG_env_2_filtered_third_variant_no_ref$Isolates_env_2 <- round((data_algG_env_2_filtered_third_variant_no_ref$SNP_percentage_env_2*100)/(100/isolates_env_2)),z<-0)



# Clean data

data_algG_env_2_clean <- data_algG_env_2_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_clean_third_base <- data_algG_env_2_filtered_third_variant[,which(names(data_algG_env_2_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_2","Isolates_env_2"))],z<-0)
ifelse(exists("data_algG_env_2_filtered_third_variant") == TRUE,data_algG_env_2_clean <- rbind(data_algG_env_2_clean,data_algG_env_2_clean_third_base),z<-0)


ifelse(exists("data_algG_env_2_filtered_third_variant_no_ref") == TRUE,data_algG_env_2_clean_third_base_no_ref <- data_algG_env_2_filtered_third_variant_no_ref[,which(names(data_algG_env_2_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_2","Isolates_env_2"))],z<-0)
ifelse(exists("data_algG_env_2_filtered_third_variant_no_ref") == TRUE,data_algG_env_2_clean <- rbind(data_algG_env_2_clean,data_algG_env_2_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_env_2_clean_20_filter <- merge(data_algG_env_2_clean, data_algG_env_2_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_env_2_clean_cleared <- data_algG_env_2_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_env_2_clean_20_filter)){
  if(data_algG_env_2_clean_20_filter$SNP_base[i] == "A" & data_algG_env_2_clean_20_filter$As_env_2[i] >= 20){
    data_algG_env_2_clean_cleared <- rbind(data_algG_env_2_clean_cleared,data_algG_env_2_clean_20_filter[i,])
  } else if(data_algG_env_2_clean_20_filter$SNP_base[i] == "C" & data_algG_env_2_clean_20_filter$Cs_env_2[i] >= 20){
    data_algG_env_2_clean_cleared <- rbind(data_algG_env_2_clean_cleared,data_algG_env_2_clean_20_filter[i,])
  } else if(data_algG_env_2_clean_20_filter$SNP_base[i] == "G" & data_algG_env_2_clean_20_filter$Gs_env_2[i] >= 20){
    data_algG_env_2_clean_cleared <- rbind(data_algG_env_2_clean_cleared,data_algG_env_2_clean_20_filter[i,])
  }else if(data_algG_env_2_clean_20_filter$SNP_base[i] == "T" & data_algG_env_2_clean_20_filter$Ts_env_2[i] >= 20){
    data_algG_env_2_clean_cleared <- rbind(data_algG_env_2_clean_cleared,data_algG_env_2_clean_20_filter[i,])
  }
}

data_algG_env_2_clean <- data_algG_env_2_clean_cleared[,c(1:5)]


### env_3

# Input 
data_algG_env_3 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T23/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_env_3_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T23/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_env_3_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T23/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_env_3 <- data_algG_env_3[!(data_algG_env_3$Position %in% primer_correction_algG_env_3_fwd$Position),]
data_algG_env_3 <- data_algG_env_3[!(data_algG_env_3$Position %in% primer_correction_algG_env_3_rev$Position),]
data_algG_env_3 <- rbind(data_algG_env_3,primer_correction_algG_env_3_fwd,primer_correction_algG_env_3_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_env_3 <- subset(data_algG_env_3, coverage_env_3 >199)
data_algG_env_3_5UTR <- subset(data_algG_env_3_5UTR, coverage_env_3_5UTR >199)
data_algG_env_3_3UTR <- subset(data_algG_env_3_3UTR, coverage_env_3_3UTR >199)

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

# Get positions with only SNP

data_algG_env_3_only_SNP <- subset(data_algG_env_3, second_percentage <= 0.025) 
data_algG_env_3_only_SNP$main_base_env_3 <- colnames(data_algG_env_3_only_SNP[, 8:11])[apply(data_algG_env_3_only_SNP[, 8:11],1,which.max)]
data_algG_env_3_only_SNP$main_base_env_3 <- ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_As_env_3","A", ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_Cs_env_3","C", ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_Gs_env_3","G", ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_Ts_env_3","T",z<-0))))
data_algG_env_3_only_SNP <- subset(data_algG_env_3_only_SNP,GenBase_PA14 != main_base_env_3)
ifelse(nrow(data_algG_env_3_only_SNP) > 0,data_algG_env_3_only_SNP$Isolates_env_3 <- isolates_env_3, z <- 0)


ifelse(nrow(data_algG_env_3_only_SNP) > 0,data_algG_env_3_only_SNP_clean <- data_algG_env_3_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_env_3_only_SNP) > 0,names(data_algG_env_3_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_3","Isolates_env_3"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_env_3_filtered <- subset(data_algG_env_3, second_percentage >= 0.025) 
data_algG_env_3_filtered$main_base_env_3 <- colnames(data_algG_env_3_filtered[, 8:11])[apply(data_algG_env_3_filtered[, 8:11],1,which.max)]
data_algG_env_3_filtered$second_base_env_3 <- colnames(data_algG_env_3_filtered[, 8:11])[apply(data_algG_env_3_filtered[, 8:11], 1, maxn(2))]
data_algG_env_3_filtered$third_base_env_3 <- ifelse(apply(data_algG_env_3_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_env_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_3_filtered$`PA14-Koordinate`)
  if(data_algG_env_3_filtered$main_base_env_3[a] == "percent_As_env_3"){
    data_algG_env_3_filtered$main_base_env_3[a] <- "A"                     
  } else if(data_algG_env_3_filtered$main_base_env_3[a] == "percent_Cs_env_3"){
    data_algG_env_3_filtered$main_base_env_3[a] <- "C"                     
  } else if(data_algG_env_3_filtered$main_base_env_3[a] == "percent_Gs_env_3"){
    data_algG_env_3_filtered$main_base_env_3[a] <- "G"                     
  } else{
    data_algG_env_3_filtered$main_base_env_3[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_env_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_3_filtered$`PA14-Koordinate`)
  if(data_algG_env_3_filtered$second_base_env_3[a] == "percent_As_env_3"){
    data_algG_env_3_filtered$second_base_env_3[a] <- "A"                     
  } else if(data_algG_env_3_filtered$second_base_env_3[a] == "percent_Cs_env_3"){
    data_algG_env_3_filtered$second_base_env_3[a] <- "C"                     
  } else if(data_algG_env_3_filtered$second_base_env_3[a] == "percent_Gs_env_3"){
    data_algG_env_3_filtered$second_base_env_3[a] <- "G"                     
  } else{
    data_algG_env_3_filtered$second_base_env_3[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_env_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_3_filtered$`PA14-Koordinate`)
  if(data_algG_env_3_filtered$main_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$third_base_env_3[a] == 0){
    data_algG_env_3_filtered$third_base_env_3[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_env_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_3_filtered$`PA14-Koordinate`)
  if(data_algG_env_3_filtered$third_base_env_3[a] == 1){
    data_algG_env_3_filtered_third_variant <- data_algG_env_3_filtered[a,]
    data_algG_env_3_filtered$third_base_env_3[a] <- 0
  }
  rm(a)
}


for(i in data_algG_env_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_3_filtered$`PA14-Koordinate`)
  if(data_algG_env_3_filtered$third_base_env_3[a] == 2){
    data_algG_env_3_filtered_third_variant_no_ref <- data_algG_env_3_filtered[a,]
    data_algG_env_3_filtered$third_base_env_3[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_env_3_filtered$third_base_env_3) == 0, data_algG_env_3_filtered <- subset(data_algG_env_3_filtered, select = -c(third_base_env_3)),data_algG_env_3_filtered <- data_algG_env_3_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_filtered_third_variant$third_percentage <- apply(data_algG_env_3_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_3_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_filtered_third_variant$third_base_env_3 <- colnames(data_algG_env_3_filtered_third_variant[, 8:11])[apply(data_algG_env_3_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_As_env_3","A", ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_Cs_env_3","C", ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_Gs_env_3","G", ifelse(data_algG_env_3_only_SNP$main_base_env_3 == "percent_Ts_env_3","T",z<-0))))

ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,
       data_algG_env_3_filtered_third_variant$third_base_env_3 <- 
         ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == "percent_As_env_3","A", 
                ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == "percent_Cs_env_3","C",
                       ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == "percent_Gs_env_3","G",
                              ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == "percent_Ts_env_3","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,
       data_algG_env_3_filtered_third_variant$third_base_env_3 <- ifelse(
         data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14,data_algG_env_3_filtered_third_variant$third_base_env_3 <- data_algG_env_3_filtered_third_variant$second_base_env_3,data_algG_env_3_filtered_third_variant$third_base_env_3
       ),
       z<-0
)

ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,
       data_algG_env_3_filtered_third_variant$third_percentage <- ifelse(
         data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14, data_algG_env_3_filtered_third_variant$second_percentage,data_algG_env_3_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_filtered_third_variant <- data_algG_env_3_filtered_third_variant[,-which(names(data_algG_env_3_filtered_third_variant) %in% c("second_percentage","second_base_env_3"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_env_3_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_3_filtered_third_variant_no_ref)[names(data_algG_env_3_filtered_third_variant_no_ref) == "second_base_env_3"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_env_3_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_3_filtered_third_variant_no_ref)[names(data_algG_env_3_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_env_3",z<-0)




# Define SNP percentage & Add it

data_algG_env_3_filtered$SNP_base <- 0
data_algG_env_3_filtered$SNP_percentage_env_3 <- 0


for(i in data_algG_env_3_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_3_filtered$`PA14-Koordinate`)
  if(data_algG_env_3_filtered$main_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] == "A"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_As_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] == "C"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Cs_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] == "G"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Gs_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] == "T"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Ts_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_3_filtered$second_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "A"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_As_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_3_filtered$second_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "C"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Cs_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_3_filtered$second_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "G"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Gs_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_3_filtered$second_base_env_3[a] == data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "T"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Ts_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "A"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_As_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "C"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Cs_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "G"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Gs_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_3_filtered$main_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$second_base_env_3[a] != data_algG_env_3_filtered$GenBase_PA14[a] & data_algG_env_3_filtered$main_base_env_3[a] == "T"){
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- data_algG_env_3_filtered$percent_Ts_env_3[a]
    data_algG_env_3_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_env_3_filtered$SNP_percentage_env_3[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_filtered_third_variant$SNP_percentage_env_3 <- 0,z<-0)


ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,
       data_algG_env_3_filtered_third_variant$SNP_percentage_env_3 <- 
         ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "A",
                data_algG_env_3_filtered_third_variant$percent_As_env_3, 
                ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "C",
                       data_algG_env_3_filtered_third_variant$percent_Cs_env_3,
                       ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "G",
                              data_algG_env_3_filtered_third_variant$percent_Gs_env_3,
                              ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "T",
                                     data_algG_env_3_filtered_third_variant$percent_Ts_env_3,
                                     ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "A",
                                            data_algG_env_3_filtered_third_variant$percent_As_env_3,
                                            ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "C",
                                                   data_algG_env_3_filtered_third_variant$percent_Cs_env_3,
                                                   ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "G",
                                                          data_algG_env_3_filtered_third_variant$percent_Gs_env_3,
                                                          ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "T",
                                                                 data_algG_env_3_filtered_third_variant$percent_Ts_env_3,
                                                                 ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "A",
                                                                        data_algG_env_3_filtered_third_variant$percent_As_env_3,
                                                                        ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "C",
                                                                               data_algG_env_3_filtered_third_variant$percent_Cs_env_3,
                                                                               ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "G",
                                                                                      data_algG_env_3_filtered_third_variant$percent_Gs_env_3,
                                                                                      ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "T",
                                                                                             data_algG_env_3_filtered_third_variant$percent_Ts_env_3,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,
       data_algG_env_3_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "A",
                "A", 
                ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "C",
                       "C",
                       ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "G",
                              "G",
                              ifelse(data_algG_env_3_filtered_third_variant$main_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "T",
                                     "T",
                                     ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "A",
                                            "A",
                                            ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "C",
                                                   "C",
                                                   ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "G",
                                                          "G",
                                                          ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 == data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$main_base_env_3 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_env_3_filtered_third_variant$third_base_env_3 != data_algG_env_3_filtered_third_variant$GenBase_PA14 & data_algG_env_3_filtered_third_variant$third_base_env_3 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_env_3_filtered$Isolates_env_3 <- round((data_algG_env_3_filtered$SNP_percentage_env_3*100)/(100/isolates_env_3))
ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_filtered_third_variant$Isolates_env_3 <- round((data_algG_env_3_filtered_third_variant$SNP_percentage_env_3*100)/(100/isolates_env_3)),z<-0)
ifelse(exists("data_algG_env_3_filtered_third_variant_no_ref") == TRUE,data_algG_env_3_filtered_third_variant_no_ref$Isolates_env_3 <- round((data_algG_env_3_filtered_third_variant_no_ref$SNP_percentage_env_3*100)/(100/isolates_env_3)),z<-0)



# Clean data

data_algG_env_3_clean <- data_algG_env_3_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_clean_third_base <- data_algG_env_3_filtered_third_variant[,which(names(data_algG_env_3_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_3","Isolates_env_3"))],z<-0)
ifelse(exists("data_algG_env_3_filtered_third_variant") == TRUE,data_algG_env_3_clean <- rbind(data_algG_env_3_clean,data_algG_env_3_clean_third_base),z<-0)


ifelse(exists("data_algG_env_3_filtered_third_variant_no_ref") == TRUE,data_algG_env_3_clean_third_base_no_ref <- data_algG_env_3_filtered_third_variant_no_ref[,which(names(data_algG_env_3_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_3","Isolates_env_3"))],z<-0)
ifelse(exists("data_algG_env_3_filtered_third_variant_no_ref") == TRUE,data_algG_env_3_clean <- rbind(data_algG_env_3_clean,data_algG_env_3_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_env_3_clean_20_filter <- merge(data_algG_env_3_clean, data_algG_env_3_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_env_3_clean_cleared <- data_algG_env_3_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_env_3_clean_20_filter)){
  if(data_algG_env_3_clean_20_filter$SNP_base[i] == "A" & data_algG_env_3_clean_20_filter$As_env_3[i] >= 20){
    data_algG_env_3_clean_cleared <- rbind(data_algG_env_3_clean_cleared,data_algG_env_3_clean_20_filter[i,])
  } else if(data_algG_env_3_clean_20_filter$SNP_base[i] == "C" & data_algG_env_3_clean_20_filter$Cs_env_3[i] >= 20){
    data_algG_env_3_clean_cleared <- rbind(data_algG_env_3_clean_cleared,data_algG_env_3_clean_20_filter[i,])
  } else if(data_algG_env_3_clean_20_filter$SNP_base[i] == "G" & data_algG_env_3_clean_20_filter$Gs_env_3[i] >= 20){
    data_algG_env_3_clean_cleared <- rbind(data_algG_env_3_clean_cleared,data_algG_env_3_clean_20_filter[i,])
  }else if(data_algG_env_3_clean_20_filter$SNP_base[i] == "T" & data_algG_env_3_clean_20_filter$Ts_env_3[i] >= 20){
    data_algG_env_3_clean_cleared <- rbind(data_algG_env_3_clean_cleared,data_algG_env_3_clean_20_filter[i,])
  }
}

data_algG_env_3_clean <- data_algG_env_3_clean_cleared[,c(1:5)]

### env_4

# Input 
data_algG_env_4 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T24/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_env_4_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T24/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_env_4_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T24/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_env_4 <- data_algG_env_4[!(data_algG_env_4$Position %in% primer_correction_algG_env_4_fwd$Position),]
data_algG_env_4 <- data_algG_env_4[!(data_algG_env_4$Position %in% primer_correction_algG_env_4_rev$Position),]
data_algG_env_4 <- rbind(data_algG_env_4,primer_correction_algG_env_4_fwd,primer_correction_algG_env_4_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_env_4 <- subset(data_algG_env_4, coverage_env_4 >199)
data_algG_env_4_5UTR <- subset(data_algG_env_4_5UTR, coverage_env_4_5UTR >199)
data_algG_env_4_3UTR <- subset(data_algG_env_4_3UTR, coverage_env_4_3UTR >199)


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

# Get positions with only SNP

data_algG_env_4_only_SNP <- subset(data_algG_env_4, second_percentage <= 0.025) 
data_algG_env_4_only_SNP$main_base_env_4 <- colnames(data_algG_env_4_only_SNP[, 8:11])[apply(data_algG_env_4_only_SNP[, 8:11],1,which.max)]
data_algG_env_4_only_SNP$main_base_env_4 <- ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_As_env_4","A", ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_Cs_env_4","C", ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_Gs_env_4","G", ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_Ts_env_4","T",z<-0))))
data_algG_env_4_only_SNP <- subset(data_algG_env_4_only_SNP,GenBase_PA14 != main_base_env_4)
ifelse(nrow(data_algG_env_4_only_SNP) > 0,data_algG_env_4_only_SNP$Isolates_env_4 <- isolates_env_4, z <- 0)


ifelse(nrow(data_algG_env_4_only_SNP) > 0,data_algG_env_4_only_SNP_clean <- data_algG_env_4_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_env_4_only_SNP) > 0,names(data_algG_env_4_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_4","Isolates_env_4"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_env_4_filtered <- subset(data_algG_env_4, second_percentage >= 0.025) 
data_algG_env_4_filtered$main_base_env_4 <- colnames(data_algG_env_4_filtered[, 8:11])[apply(data_algG_env_4_filtered[, 8:11],1,which.max)]
data_algG_env_4_filtered$second_base_env_4 <- colnames(data_algG_env_4_filtered[, 8:11])[apply(data_algG_env_4_filtered[, 8:11], 1, maxn(2))]
data_algG_env_4_filtered$third_base_env_4 <- ifelse(apply(data_algG_env_4_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_env_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_4_filtered$`PA14-Koordinate`)
  if(data_algG_env_4_filtered$main_base_env_4[a] == "percent_As_env_4"){
    data_algG_env_4_filtered$main_base_env_4[a] <- "A"                     
  } else if(data_algG_env_4_filtered$main_base_env_4[a] == "percent_Cs_env_4"){
    data_algG_env_4_filtered$main_base_env_4[a] <- "C"                     
  } else if(data_algG_env_4_filtered$main_base_env_4[a] == "percent_Gs_env_4"){
    data_algG_env_4_filtered$main_base_env_4[a] <- "G"                     
  } else{
    data_algG_env_4_filtered$main_base_env_4[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_env_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_4_filtered$`PA14-Koordinate`)
  if(data_algG_env_4_filtered$second_base_env_4[a] == "percent_As_env_4"){
    data_algG_env_4_filtered$second_base_env_4[a] <- "A"                     
  } else if(data_algG_env_4_filtered$second_base_env_4[a] == "percent_Cs_env_4"){
    data_algG_env_4_filtered$second_base_env_4[a] <- "C"                     
  } else if(data_algG_env_4_filtered$second_base_env_4[a] == "percent_Gs_env_4"){
    data_algG_env_4_filtered$second_base_env_4[a] <- "G"                     
  } else{
    data_algG_env_4_filtered$second_base_env_4[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_env_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_4_filtered$`PA14-Koordinate`)
  if(data_algG_env_4_filtered$main_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$third_base_env_4[a] == 0){
    data_algG_env_4_filtered$third_base_env_4[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_env_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_4_filtered$`PA14-Koordinate`)
  if(data_algG_env_4_filtered$third_base_env_4[a] == 1){
    data_algG_env_4_filtered_third_variant <- data_algG_env_4_filtered[a,]
    data_algG_env_4_filtered$third_base_env_4[a] <- 0
  }
  rm(a)
}


for(i in data_algG_env_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_4_filtered$`PA14-Koordinate`)
  if(data_algG_env_4_filtered$third_base_env_4[a] == 2){
    data_algG_env_4_filtered_third_variant_no_ref <- data_algG_env_4_filtered[a,]
    data_algG_env_4_filtered$third_base_env_4[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_env_4_filtered$third_base_env_4) == 0, data_algG_env_4_filtered <- subset(data_algG_env_4_filtered, select = -c(third_base_env_4)),data_algG_env_4_filtered <- data_algG_env_4_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_filtered_third_variant$third_percentage <- apply(data_algG_env_4_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_4_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_filtered_third_variant$third_base_env_4 <- colnames(data_algG_env_4_filtered_third_variant[, 8:11])[apply(data_algG_env_4_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_As_env_4","A", ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_Cs_env_4","C", ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_Gs_env_4","G", ifelse(data_algG_env_4_only_SNP$main_base_env_4 == "percent_Ts_env_4","T",z<-0))))

ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,
       data_algG_env_4_filtered_third_variant$third_base_env_4 <- 
         ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == "percent_As_env_4","A", 
                ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == "percent_Cs_env_4","C",
                       ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == "percent_Gs_env_4","G",
                              ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == "percent_Ts_env_4","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,
       data_algG_env_4_filtered_third_variant$third_base_env_4 <- ifelse(
         data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14,data_algG_env_4_filtered_third_variant$third_base_env_4 <- data_algG_env_4_filtered_third_variant$second_base_env_4,data_algG_env_4_filtered_third_variant$third_base_env_4
       ),
       z<-0
)

ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,
       data_algG_env_4_filtered_third_variant$third_percentage <- ifelse(
         data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14, data_algG_env_4_filtered_third_variant$second_percentage,data_algG_env_4_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_filtered_third_variant <- data_algG_env_4_filtered_third_variant[,-which(names(data_algG_env_4_filtered_third_variant) %in% c("second_percentage","second_base_env_4"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_env_4_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_4_filtered_third_variant_no_ref)[names(data_algG_env_4_filtered_third_variant_no_ref) == "second_base_env_4"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_env_4_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_4_filtered_third_variant_no_ref)[names(data_algG_env_4_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_env_4",z<-0)




# Define SNP percentage & Add it

data_algG_env_4_filtered$SNP_base <- 0
data_algG_env_4_filtered$SNP_percentage_env_4 <- 0


for(i in data_algG_env_4_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_4_filtered$`PA14-Koordinate`)
  if(data_algG_env_4_filtered$main_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] == "A"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_As_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] == "C"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Cs_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] == "G"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Gs_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] == "T"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Ts_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_4_filtered$second_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "A"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_As_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_4_filtered$second_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "C"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Cs_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_4_filtered$second_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "G"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Gs_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_4_filtered$second_base_env_4[a] == data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "T"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Ts_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "A"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_As_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "C"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Cs_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "G"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Gs_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_4_filtered$main_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$second_base_env_4[a] != data_algG_env_4_filtered$GenBase_PA14[a] & data_algG_env_4_filtered$main_base_env_4[a] == "T"){
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- data_algG_env_4_filtered$percent_Ts_env_4[a]
    data_algG_env_4_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_env_4_filtered$SNP_percentage_env_4[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_filtered_third_variant$SNP_percentage_env_4 <- 0,z<-0)


ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,
       data_algG_env_4_filtered_third_variant$SNP_percentage_env_4 <- 
         ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "A",
                data_algG_env_4_filtered_third_variant$percent_As_env_4, 
                ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "C",
                       data_algG_env_4_filtered_third_variant$percent_Cs_env_4,
                       ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "G",
                              data_algG_env_4_filtered_third_variant$percent_Gs_env_4,
                              ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "T",
                                     data_algG_env_4_filtered_third_variant$percent_Ts_env_4,
                                     ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "A",
                                            data_algG_env_4_filtered_third_variant$percent_As_env_4,
                                            ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "C",
                                                   data_algG_env_4_filtered_third_variant$percent_Cs_env_4,
                                                   ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "G",
                                                          data_algG_env_4_filtered_third_variant$percent_Gs_env_4,
                                                          ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "T",
                                                                 data_algG_env_4_filtered_third_variant$percent_Ts_env_4,
                                                                 ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "A",
                                                                        data_algG_env_4_filtered_third_variant$percent_As_env_4,
                                                                        ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "C",
                                                                               data_algG_env_4_filtered_third_variant$percent_Cs_env_4,
                                                                               ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "G",
                                                                                      data_algG_env_4_filtered_third_variant$percent_Gs_env_4,
                                                                                      ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "T",
                                                                                             data_algG_env_4_filtered_third_variant$percent_Ts_env_4,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,
       data_algG_env_4_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "A",
                "A", 
                ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "C",
                       "C",
                       ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "G",
                              "G",
                              ifelse(data_algG_env_4_filtered_third_variant$main_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "T",
                                     "T",
                                     ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "A",
                                            "A",
                                            ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "C",
                                                   "C",
                                                   ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "G",
                                                          "G",
                                                          ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 == data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$main_base_env_4 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_env_4_filtered_third_variant$third_base_env_4 != data_algG_env_4_filtered_third_variant$GenBase_PA14 & data_algG_env_4_filtered_third_variant$third_base_env_4 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_env_4_filtered$Isolates_env_4 <- round((data_algG_env_4_filtered$SNP_percentage_env_4*100)/(100/isolates_env_4))
ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_filtered_third_variant$Isolates_env_4 <- round((data_algG_env_4_filtered_third_variant$SNP_percentage_env_4*100)/(100/isolates_env_4)),z<-0)
ifelse(exists("data_algG_env_4_filtered_third_variant_no_ref") == TRUE,data_algG_env_4_filtered_third_variant_no_ref$Isolates_env_4 <- round((data_algG_env_4_filtered_third_variant_no_ref$SNP_percentage_env_4*100)/(100/isolates_env_4)),z<-0)



# Clean data

data_algG_env_4_clean <- data_algG_env_4_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_clean_third_base <- data_algG_env_4_filtered_third_variant[,which(names(data_algG_env_4_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_4","Isolates_env_4"))],z<-0)
ifelse(exists("data_algG_env_4_filtered_third_variant") == TRUE,data_algG_env_4_clean <- rbind(data_algG_env_4_clean,data_algG_env_4_clean_third_base),z<-0)


ifelse(exists("data_algG_env_4_filtered_third_variant_no_ref") == TRUE,data_algG_env_4_clean_third_base_no_ref <- data_algG_env_4_filtered_third_variant_no_ref[,which(names(data_algG_env_4_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_4","Isolates_env_4"))],z<-0)
ifelse(exists("data_algG_env_4_filtered_third_variant_no_ref") == TRUE,data_algG_env_4_clean <- rbind(data_algG_env_4_clean,data_algG_env_4_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_env_4_clean_20_filter <- merge(data_algG_env_4_clean, data_algG_env_4_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_env_4_clean_cleared <- data_algG_env_4_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_env_4_clean_20_filter)){
  if(data_algG_env_4_clean_20_filter$SNP_base[i] == "A" & data_algG_env_4_clean_20_filter$As_env_4[i] >= 20){
    data_algG_env_4_clean_cleared <- rbind(data_algG_env_4_clean_cleared,data_algG_env_4_clean_20_filter[i,])
  } else if(data_algG_env_4_clean_20_filter$SNP_base[i] == "C" & data_algG_env_4_clean_20_filter$Cs_env_4[i] >= 20){
    data_algG_env_4_clean_cleared <- rbind(data_algG_env_4_clean_cleared,data_algG_env_4_clean_20_filter[i,])
  } else if(data_algG_env_4_clean_20_filter$SNP_base[i] == "G" & data_algG_env_4_clean_20_filter$Gs_env_4[i] >= 20){
    data_algG_env_4_clean_cleared <- rbind(data_algG_env_4_clean_cleared,data_algG_env_4_clean_20_filter[i,])
  }else if(data_algG_env_4_clean_20_filter$SNP_base[i] == "T" & data_algG_env_4_clean_20_filter$Ts_env_4[i] >= 20){
    data_algG_env_4_clean_cleared <- rbind(data_algG_env_4_clean_cleared,data_algG_env_4_clean_20_filter[i,])
  }
}

data_algG_env_4_clean <- data_algG_env_4_clean_cleared[,c(1:5)]


### env_5

# Input 
data_algG_env_5 <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T25/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_env_5_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T25/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_env_5_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T25/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_env_5 <- data_algG_env_5[!(data_algG_env_5$Position %in% primer_correction_algG_env_5_fwd$Position),]
data_algG_env_5 <- data_algG_env_5[!(data_algG_env_5$Position %in% primer_correction_algG_env_5_rev$Position),]
data_algG_env_5 <- rbind(data_algG_env_5,primer_correction_algG_env_5_fwd,primer_correction_algG_env_5_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_env_5 <- subset(data_algG_env_5, coverage_env_5 >199)
data_algG_env_5_5UTR <- subset(data_algG_env_5_5UTR, coverage_env_5_5UTR >199)
data_algG_env_5_3UTR <- subset(data_algG_env_5_3UTR, coverage_env_5_3UTR >199)

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

# Get positions with only SNP

data_algG_env_5_only_SNP <- subset(data_algG_env_5, second_percentage <= 0.025) 
data_algG_env_5_only_SNP$main_base_env_5 <- colnames(data_algG_env_5_only_SNP[, 8:11])[apply(data_algG_env_5_only_SNP[, 8:11],1,which.max)]
data_algG_env_5_only_SNP$main_base_env_5 <- ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_As_env_5","A", ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_Cs_env_5","C", ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_Gs_env_5","G", ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_Ts_env_5","T",z<-0))))
data_algG_env_5_only_SNP <- subset(data_algG_env_5_only_SNP,GenBase_PA14 != main_base_env_5)
ifelse(nrow(data_algG_env_5_only_SNP) > 0,data_algG_env_5_only_SNP$Isolates_env_5 <- isolates_env_5, z <- 0)


ifelse(nrow(data_algG_env_5_only_SNP) > 0,data_algG_env_5_only_SNP_clean <- data_algG_env_5_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_env_5_only_SNP) > 0,names(data_algG_env_5_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_5","Isolates_env_5"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_env_5_filtered <- subset(data_algG_env_5, second_percentage >= 0.025) 
data_algG_env_5_filtered$main_base_env_5 <- colnames(data_algG_env_5_filtered[, 8:11])[apply(data_algG_env_5_filtered[, 8:11],1,which.max)]
data_algG_env_5_filtered$second_base_env_5 <- colnames(data_algG_env_5_filtered[, 8:11])[apply(data_algG_env_5_filtered[, 8:11], 1, maxn(2))]
data_algG_env_5_filtered$third_base_env_5 <- ifelse(apply(data_algG_env_5_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_env_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_5_filtered$`PA14-Koordinate`)
  if(data_algG_env_5_filtered$main_base_env_5[a] == "percent_As_env_5"){
    data_algG_env_5_filtered$main_base_env_5[a] <- "A"                     
  } else if(data_algG_env_5_filtered$main_base_env_5[a] == "percent_Cs_env_5"){
    data_algG_env_5_filtered$main_base_env_5[a] <- "C"                     
  } else if(data_algG_env_5_filtered$main_base_env_5[a] == "percent_Gs_env_5"){
    data_algG_env_5_filtered$main_base_env_5[a] <- "G"                     
  } else{
    data_algG_env_5_filtered$main_base_env_5[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_env_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_5_filtered$`PA14-Koordinate`)
  if(data_algG_env_5_filtered$second_base_env_5[a] == "percent_As_env_5"){
    data_algG_env_5_filtered$second_base_env_5[a] <- "A"                     
  } else if(data_algG_env_5_filtered$second_base_env_5[a] == "percent_Cs_env_5"){
    data_algG_env_5_filtered$second_base_env_5[a] <- "C"                     
  } else if(data_algG_env_5_filtered$second_base_env_5[a] == "percent_Gs_env_5"){
    data_algG_env_5_filtered$second_base_env_5[a] <- "G"                     
  } else{
    data_algG_env_5_filtered$second_base_env_5[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_env_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_5_filtered$`PA14-Koordinate`)
  if(data_algG_env_5_filtered$main_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$third_base_env_5[a] == 0){
    data_algG_env_5_filtered$third_base_env_5[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_env_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_5_filtered$`PA14-Koordinate`)
  if(data_algG_env_5_filtered$third_base_env_5[a] == 1){
    data_algG_env_5_filtered_third_variant <- data_algG_env_5_filtered[a,]
    data_algG_env_5_filtered$third_base_env_5[a] <- 0
  }
  rm(a)
}


for(i in data_algG_env_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_5_filtered$`PA14-Koordinate`)
  if(data_algG_env_5_filtered$third_base_env_5[a] == 2){
    data_algG_env_5_filtered_third_variant_no_ref <- data_algG_env_5_filtered[a,]
    data_algG_env_5_filtered$third_base_env_5[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_env_5_filtered$third_base_env_5) == 0, data_algG_env_5_filtered <- subset(data_algG_env_5_filtered, select = -c(third_base_env_5)),data_algG_env_5_filtered <- data_algG_env_5_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_filtered_third_variant$third_percentage <- apply(data_algG_env_5_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_env_5_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_filtered_third_variant$third_base_env_5 <- colnames(data_algG_env_5_filtered_third_variant[, 8:11])[apply(data_algG_env_5_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_As_env_5","A", ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_Cs_env_5","C", ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_Gs_env_5","G", ifelse(data_algG_env_5_only_SNP$main_base_env_5 == "percent_Ts_env_5","T",z<-0))))

ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,
       data_algG_env_5_filtered_third_variant$third_base_env_5 <- 
         ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == "percent_As_env_5","A", 
                ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == "percent_Cs_env_5","C",
                       ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == "percent_Gs_env_5","G",
                              ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == "percent_Ts_env_5","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,
       data_algG_env_5_filtered_third_variant$third_base_env_5 <- ifelse(
         data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14,data_algG_env_5_filtered_third_variant$third_base_env_5 <- data_algG_env_5_filtered_third_variant$second_base_env_5,data_algG_env_5_filtered_third_variant$third_base_env_5
       ),
       z<-0
)

ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,
       data_algG_env_5_filtered_third_variant$third_percentage <- ifelse(
         data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14, data_algG_env_5_filtered_third_variant$second_percentage,data_algG_env_5_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_filtered_third_variant <- data_algG_env_5_filtered_third_variant[,-which(names(data_algG_env_5_filtered_third_variant) %in% c("second_percentage","second_base_env_5"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_env_5_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_5_filtered_third_variant_no_ref)[names(data_algG_env_5_filtered_third_variant_no_ref) == "second_base_env_5"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_env_5_filtered_third_variant_no_ref") == TRUE, names(data_algG_env_5_filtered_third_variant_no_ref)[names(data_algG_env_5_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_env_5",z<-0)




# Define SNP percentage & Add it

data_algG_env_5_filtered$SNP_base <- 0
data_algG_env_5_filtered$SNP_percentage_env_5 <- 0


for(i in data_algG_env_5_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_env_5_filtered$`PA14-Koordinate`)
  if(data_algG_env_5_filtered$main_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] == "A"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_As_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] == "C"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Cs_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] == "G"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Gs_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] == "T"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Ts_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_5_filtered$second_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "A"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_As_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_5_filtered$second_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "C"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Cs_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_5_filtered$second_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "G"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Gs_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_5_filtered$second_base_env_5[a] == data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "T"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Ts_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "T"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "A"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_As_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "A"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "C"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Cs_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "C"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "G"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Gs_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "G"
  } else if(data_algG_env_5_filtered$main_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$second_base_env_5[a] != data_algG_env_5_filtered$GenBase_PA14[a] & data_algG_env_5_filtered$main_base_env_5[a] == "T"){
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- data_algG_env_5_filtered$percent_Ts_env_5[a]
    data_algG_env_5_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_env_5_filtered$SNP_percentage_env_5[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_filtered_third_variant$SNP_percentage_env_5 <- 0,z<-0)


ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,
       data_algG_env_5_filtered_third_variant$SNP_percentage_env_5 <- 
         ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "A",
                data_algG_env_5_filtered_third_variant$percent_As_env_5, 
                ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "C",
                       data_algG_env_5_filtered_third_variant$percent_Cs_env_5,
                       ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "G",
                              data_algG_env_5_filtered_third_variant$percent_Gs_env_5,
                              ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "T",
                                     data_algG_env_5_filtered_third_variant$percent_Ts_env_5,
                                     ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "A",
                                            data_algG_env_5_filtered_third_variant$percent_As_env_5,
                                            ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "C",
                                                   data_algG_env_5_filtered_third_variant$percent_Cs_env_5,
                                                   ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "G",
                                                          data_algG_env_5_filtered_third_variant$percent_Gs_env_5,
                                                          ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "T",
                                                                 data_algG_env_5_filtered_third_variant$percent_Ts_env_5,
                                                                 ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "A",
                                                                        data_algG_env_5_filtered_third_variant$percent_As_env_5,
                                                                        ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "C",
                                                                               data_algG_env_5_filtered_third_variant$percent_Cs_env_5,
                                                                               ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "G",
                                                                                      data_algG_env_5_filtered_third_variant$percent_Gs_env_5,
                                                                                      ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "T",
                                                                                             data_algG_env_5_filtered_third_variant$percent_Ts_env_5,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,
       data_algG_env_5_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "A",
                "A", 
                ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "C",
                       "C",
                       ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "G",
                              "G",
                              ifelse(data_algG_env_5_filtered_third_variant$main_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "T",
                                     "T",
                                     ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "A",
                                            "A",
                                            ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "C",
                                                   "C",
                                                   ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "G",
                                                          "G",
                                                          ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 == data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$main_base_env_5 == "T",
                                                                 "T",
                                                                 ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "A",
                                                                        "A",
                                                                        ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "C",
                                                                               "C",
                                                                               ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_env_5_filtered_third_variant$third_base_env_5 != data_algG_env_5_filtered_third_variant$GenBase_PA14 & data_algG_env_5_filtered_third_variant$third_base_env_5 == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_env_5_filtered$Isolates_env_5 <- round((data_algG_env_5_filtered$SNP_percentage_env_5*100)/(100/isolates_env_5))
ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_filtered_third_variant$Isolates_env_5 <- round((data_algG_env_5_filtered_third_variant$SNP_percentage_env_5*100)/(100/isolates_env_5)),z<-0)
ifelse(exists("data_algG_env_5_filtered_third_variant_no_ref") == TRUE,data_algG_env_5_filtered_third_variant_no_ref$Isolates_env_5 <- round((data_algG_env_5_filtered_third_variant_no_ref$SNP_percentage_env_5*100)/(100/isolates_env_5)),z<-0)



# Clean data

data_algG_env_5_clean <- data_algG_env_5_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_clean_third_base <- data_algG_env_5_filtered_third_variant[,which(names(data_algG_env_5_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_5","Isolates_env_5"))],z<-0)
ifelse(exists("data_algG_env_5_filtered_third_variant") == TRUE,data_algG_env_5_clean <- rbind(data_algG_env_5_clean,data_algG_env_5_clean_third_base),z<-0)


ifelse(exists("data_algG_env_5_filtered_third_variant_no_ref") == TRUE,data_algG_env_5_clean_third_base_no_ref <- data_algG_env_5_filtered_third_variant_no_ref[,which(names(data_algG_env_5_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_5","Isolates_env_5"))],z<-0)
ifelse(exists("data_algG_env_5_filtered_third_variant_no_ref") == TRUE,data_algG_env_5_clean <- rbind(data_algG_env_5_clean,data_algG_env_5_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_env_5_clean_20_filter <- merge(data_algG_env_5_clean, data_algG_env_5_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_env_5_clean_cleared <- data_algG_env_5_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_env_5_clean_20_filter)){
  if(data_algG_env_5_clean_20_filter$SNP_base[i] == "A" & data_algG_env_5_clean_20_filter$As_env_5[i] >= 20){
    data_algG_env_5_clean_cleared <- rbind(data_algG_env_5_clean_cleared,data_algG_env_5_clean_20_filter[i,])
  } else if(data_algG_env_5_clean_20_filter$SNP_base[i] == "C" & data_algG_env_5_clean_20_filter$Cs_env_5[i] >= 20){
    data_algG_env_5_clean_cleared <- rbind(data_algG_env_5_clean_cleared,data_algG_env_5_clean_20_filter[i,])
  } else if(data_algG_env_5_clean_20_filter$SNP_base[i] == "G" & data_algG_env_5_clean_20_filter$Gs_env_5[i] >= 20){
    data_algG_env_5_clean_cleared <- rbind(data_algG_env_5_clean_cleared,data_algG_env_5_clean_20_filter[i,])
  }else if(data_algG_env_5_clean_20_filter$SNP_base[i] == "T" & data_algG_env_5_clean_20_filter$Ts_env_5[i] >= 20){
    data_algG_env_5_clean_cleared <- rbind(data_algG_env_5_clean_cleared,data_algG_env_5_clean_20_filter[i,])
  }
}

data_algG_env_5_clean <- data_algG_env_5_clean_cleared[,c(1:5)]


### COPD

# Input 
data_algG_COPD <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/T26/Results/algG/Results_algG.txt", header = TRUE)
primer_correction_algG_COPD_fwd <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T26/Results/algG_only_fwd/Results_algG_only_fwd.txt", header = TRUE)
primer_correction_algG_COPD_rev <- read.table("/mnt/sfb900nfs/groups/tuemmler/sebastian/pathoadaptive_mutationen/primer_elimination/T26/Results/algG_only_rev/Results_algG_only_rev.txt", header = TRUE)

# Primer correction
data_algG_COPD <- data_algG_COPD[!(data_algG_COPD$Position %in% primer_correction_algG_COPD_fwd$Position),]
data_algG_COPD <- data_algG_COPD[!(data_algG_COPD$Position %in% primer_correction_algG_COPD_rev$Position),]
data_algG_COPD <- rbind(data_algG_COPD,primer_correction_algG_COPD_fwd,primer_correction_algG_COPD_rev)

# Reformatting
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

#Subset coverage > 200

data_algG_COPD <- subset(data_algG_COPD, coverage_COPD >199)
data_algG_COPD_5UTR <- subset(data_algG_COPD_5UTR, coverage_COPD_5UTR >199)
data_algG_COPD_3UTR <- subset(data_algG_COPD_3UTR, coverage_COPD_3UTR >199)

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

# Get positions with only SNP

data_algG_COPD_only_SNP <- subset(data_algG_COPD, second_percentage <= 0.025) 
data_algG_COPD_only_SNP$main_base_COPD <- colnames(data_algG_COPD_only_SNP[, 8:11])[apply(data_algG_COPD_only_SNP[, 8:11],1,which.max)]
data_algG_COPD_only_SNP$main_base_COPD <- ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_As_COPD","A", ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_Cs_COPD","C", ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_Gs_COPD","G", ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_Ts_COPD","T",z<-0))))
data_algG_COPD_only_SNP <- subset(data_algG_COPD_only_SNP,GenBase_PA14 != main_base_COPD)
ifelse(nrow(data_algG_COPD_only_SNP) > 0,data_algG_COPD_only_SNP$Isolates_COPD <- isolates_COPD, z <- 0)


ifelse(nrow(data_algG_COPD_only_SNP) > 0,data_algG_COPD_only_SNP_clean <- data_algG_COPD_only_SNP[,c(1:2,14,12,15)], z <- 0)
ifelse(nrow(data_algG_COPD_only_SNP) > 0,names(data_algG_COPD_only_SNP_clean) <- c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_COPD","Isolates_COPD"),z <-0)


# Subset only those with the second highest minimum 2.5 %

data_algG_COPD_filtered <- subset(data_algG_COPD, second_percentage >= 0.025) 
data_algG_COPD_filtered$main_base_COPD <- colnames(data_algG_COPD_filtered[, 8:11])[apply(data_algG_COPD_filtered[, 8:11],1,which.max)]
data_algG_COPD_filtered$second_base_COPD <- colnames(data_algG_COPD_filtered[, 8:11])[apply(data_algG_COPD_filtered[, 8:11], 1, maxn(2))]
data_algG_COPD_filtered$third_base_COPD <- ifelse(apply(data_algG_COPD_filtered[, 8:11], 1, function(x)x[maxn(3)(x)])>=0.025,1,0)

# Reformat 

for(i in data_algG_COPD_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_COPD_filtered$`PA14-Koordinate`)
  if(data_algG_COPD_filtered$main_base_COPD[a] == "percent_As_COPD"){
    data_algG_COPD_filtered$main_base_COPD[a] <- "A"                     
  } else if(data_algG_COPD_filtered$main_base_COPD[a] == "percent_Cs_COPD"){
    data_algG_COPD_filtered$main_base_COPD[a] <- "C"                     
  } else if(data_algG_COPD_filtered$main_base_COPD[a] == "percent_Gs_COPD"){
    data_algG_COPD_filtered$main_base_COPD[a] <- "G"                     
  } else{
    data_algG_COPD_filtered$main_base_COPD[a] <- "T"
  }
  rm(a)
}


for(i in data_algG_COPD_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_COPD_filtered$`PA14-Koordinate`)
  if(data_algG_COPD_filtered$second_base_COPD[a] == "percent_As_COPD"){
    data_algG_COPD_filtered$second_base_COPD[a] <- "A"                     
  } else if(data_algG_COPD_filtered$second_base_COPD[a] == "percent_Cs_COPD"){
    data_algG_COPD_filtered$second_base_COPD[a] <- "C"                     
  } else if(data_algG_COPD_filtered$second_base_COPD[a] == "percent_Gs_COPD"){
    data_algG_COPD_filtered$second_base_COPD[a] <- "G"                     
  } else{
    data_algG_COPD_filtered$second_base_COPD[a] <- "T"
  }
  rm(a)
}


## Add third variant if no reference base is present



for(i in data_algG_COPD_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_COPD_filtered$`PA14-Koordinate`)
  if(data_algG_COPD_filtered$main_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$third_base_COPD[a] == 0){
    data_algG_COPD_filtered$third_base_COPD[a] <- 2
  } else{
    z<-0
  }
  rm(a)
}


for(i in data_algG_COPD_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_COPD_filtered$`PA14-Koordinate`)
  if(data_algG_COPD_filtered$third_base_COPD[a] == 1){
    data_algG_COPD_filtered_third_variant <- data_algG_COPD_filtered[a,]
    data_algG_COPD_filtered$third_base_COPD[a] <- 0
  }
  rm(a)
}


for(i in data_algG_COPD_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_COPD_filtered$`PA14-Koordinate`)
  if(data_algG_COPD_filtered$third_base_COPD[a] == 2){
    data_algG_COPD_filtered_third_variant_no_ref <- data_algG_COPD_filtered[a,]
    data_algG_COPD_filtered$third_base_COPD[a] <- 0
  }
  rm(a)
}

ifelse(sum(data_algG_COPD_filtered$third_base_COPD) == 0, data_algG_COPD_filtered <- subset(data_algG_COPD_filtered, select = -c(third_base_COPD)),data_algG_COPD_filtered <- data_algG_COPD_filtered)

## Processing third base if ref base is present

ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_filtered_third_variant$third_percentage <- apply(data_algG_COPD_filtered_third_variant[, 8:11], 1, function(i) sort(i)[ dim(data_algG_COPD_filtered_third_variant[, 8:11])[2]-2]),z<-0)
ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_filtered_third_variant$third_base_COPD <- colnames(data_algG_COPD_filtered_third_variant[, 8:11])[apply(data_algG_COPD_filtered_third_variant[, 8:11], 1, maxn(3))],z<-0)


ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_As_COPD","A", ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_Cs_COPD","C", ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_Gs_COPD","G", ifelse(data_algG_COPD_only_SNP$main_base_COPD == "percent_Ts_COPD","T",z<-0))))

ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,
       data_algG_COPD_filtered_third_variant$third_base_COPD <- 
         ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == "percent_As_COPD","A", 
                ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == "percent_Cs_COPD","C",
                       ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == "percent_Gs_COPD","G",
                              ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == "percent_Ts_COPD","T",
                                     z <- 0)))),
       z <- 0
)         


## Define second base as third base if reference == third base

ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,
       data_algG_COPD_filtered_third_variant$third_base_COPD <- ifelse(
         data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14,data_algG_COPD_filtered_third_variant$third_base_COPD <- data_algG_COPD_filtered_third_variant$second_base_COPD,data_algG_COPD_filtered_third_variant$third_base_COPD
       ),
       z<-0
)

ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,
       data_algG_COPD_filtered_third_variant$third_percentage <- ifelse(
         data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14, data_algG_COPD_filtered_third_variant$second_percentage,data_algG_COPD_filtered_third_variant$third_percentage
       ),
       z<-0
)


ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_filtered_third_variant <- data_algG_COPD_filtered_third_variant[,-which(names(data_algG_COPD_filtered_third_variant) %in% c("second_percentage","second_base_COPD"))],z<-0)


#### Processing third base if reference is not present

ifelse(exists("data_algG_COPD_filtered_third_variant_no_ref") == TRUE, names(data_algG_COPD_filtered_third_variant_no_ref)[names(data_algG_COPD_filtered_third_variant_no_ref) == "second_base_COPD"] <- "SNP_base",z<-0)
ifelse(exists("data_algG_COPD_filtered_third_variant_no_ref") == TRUE, names(data_algG_COPD_filtered_third_variant_no_ref)[names(data_algG_COPD_filtered_third_variant_no_ref) == "second_percentage"] <- "SNP_percentage_COPD",z<-0)




# Define SNP percentage & Add it

data_algG_COPD_filtered$SNP_base <- 0
data_algG_COPD_filtered$SNP_percentage_COPD <- 0


for(i in data_algG_COPD_filtered$`PA14-Koordinate`){
  a <- which(i == data_algG_COPD_filtered$`PA14-Koordinate`)
  if(data_algG_COPD_filtered$main_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] == "A"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_As_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "A"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] == "C"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Cs_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "C"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] == "G"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Gs_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "G"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] == "T"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Ts_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "T"
  } else if(data_algG_COPD_filtered$second_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "A"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_As_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "A"
  } else if(data_algG_COPD_filtered$second_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "C"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Cs_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "C"
  } else if(data_algG_COPD_filtered$second_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "G"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Gs_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "G"
  } else if(data_algG_COPD_filtered$second_base_COPD[a] == data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "T"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Ts_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "T"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "A"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_As_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "A"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "C"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Cs_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "C"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "G"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Gs_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "G"
  } else if(data_algG_COPD_filtered$main_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$second_base_COPD[a] != data_algG_COPD_filtered$GenBase_PA14[a] & data_algG_COPD_filtered$main_base_COPD[a] == "T"){
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- data_algG_COPD_filtered$percent_Ts_COPD[a]
    data_algG_COPD_filtered$SNP_base[a] <- "T"
  } else{
    data_algG_COPD_filtered$SNP_percentage_COPD[a] <- "Fehler"                
    
  }
  rm(a)
}


ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_filtered_third_variant$SNP_base <- 0,z<-0)
ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_filtered_third_variant$SNP_percentage_COPD <- 0,z<-0)


ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,
       data_algG_COPD_filtered_third_variant$SNP_percentage_COPD <- 
         ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "A",
                data_algG_COPD_filtered_third_variant$percent_As_COPD, 
                ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "C",
                       data_algG_COPD_filtered_third_variant$percent_Cs_COPD,
                       ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "G",
                              data_algG_COPD_filtered_third_variant$percent_Gs_COPD,
                              ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "T",
                                     data_algG_COPD_filtered_third_variant$percent_Ts_COPD,
                                     ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "A",
                                            data_algG_COPD_filtered_third_variant$percent_As_COPD,
                                            ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "C",
                                                   data_algG_COPD_filtered_third_variant$percent_Cs_COPD,
                                                   ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "G",
                                                          data_algG_COPD_filtered_third_variant$percent_Gs_COPD,
                                                          ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "T",
                                                                 data_algG_COPD_filtered_third_variant$percent_Ts_COPD,
                                                                 ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "A",
                                                                        data_algG_COPD_filtered_third_variant$percent_As_COPD,
                                                                        ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "C",
                                                                               data_algG_COPD_filtered_third_variant$percent_Cs_COPD,
                                                                               ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "G",
                                                                                      data_algG_COPD_filtered_third_variant$percent_Gs_COPD,
                                                                                      ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "T",
                                                                                             data_algG_COPD_filtered_third_variant$percent_Ts_COPD,
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,
       data_algG_COPD_filtered_third_variant$SNP_base <- 
         ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "A",
                "A", 
                ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "C",
                       "C",
                       ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "G",
                              "G",
                              ifelse(data_algG_COPD_filtered_third_variant$main_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "T",
                                     "T",
                                     ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "A",
                                            "A",
                                            ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "C",
                                                   "C",
                                                   ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "G",
                                                          "G",
                                                          ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD == data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$main_base_COPD == "T",
                                                                 "T",
                                                                 ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "A",
                                                                        "A",
                                                                        ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "C",
                                                                               "C",
                                                                               ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "G",
                                                                                      "G",
                                                                                      ifelse(data_algG_COPD_filtered_third_variant$third_base_COPD != data_algG_COPD_filtered_third_variant$GenBase_PA14 & data_algG_COPD_filtered_third_variant$third_base_COPD == "T",
                                                                                             "T",
                                                                                             "Fehler"
                                                                                      )))))))))))
         ),z<-0
)


# calculate Isolates

data_algG_COPD_filtered$Isolates_COPD <- round((data_algG_COPD_filtered$SNP_percentage_COPD*100)/(100/isolates_COPD))
ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_filtered_third_variant$Isolates_COPD <- round((data_algG_COPD_filtered_third_variant$SNP_percentage_COPD*100)/(100/isolates_COPD)),z<-0)
ifelse(exists("data_algG_COPD_filtered_third_variant_no_ref") == TRUE,data_algG_COPD_filtered_third_variant_no_ref$Isolates_COPD <- round((data_algG_COPD_filtered_third_variant_no_ref$SNP_percentage_COPD*100)/(100/isolates_COPD)),z<-0)



# Clean data

data_algG_COPD_clean <- data_algG_COPD_filtered[,c(1:2,16:18)]

ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_clean_third_base <- data_algG_COPD_filtered_third_variant[,which(names(data_algG_COPD_filtered_third_variant) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_COPD","Isolates_COPD"))],z<-0)
ifelse(exists("data_algG_COPD_filtered_third_variant") == TRUE,data_algG_COPD_clean <- rbind(data_algG_COPD_clean,data_algG_COPD_clean_third_base),z<-0)


ifelse(exists("data_algG_COPD_filtered_third_variant_no_ref") == TRUE,data_algG_COPD_clean_third_base_no_ref <- data_algG_COPD_filtered_third_variant_no_ref[,which(names(data_algG_COPD_filtered_third_variant_no_ref) %in% c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_COPD","Isolates_COPD"))],z<-0)
ifelse(exists("data_algG_COPD_filtered_third_variant_no_ref") == TRUE,data_algG_COPD_clean <- rbind(data_algG_COPD_clean,data_algG_COPD_clean_third_base_no_ref),z<-0)

# filter 20 SNP Reads


data_algG_COPD_clean_20_filter <- merge(data_algG_COPD_clean, data_algG_COPD_filtered[,c(1,3:6)],by="PA14-Koordinate")
data_algG_COPD_clean_cleared <- data_algG_COPD_clean_20_filter[FALSE,]


for(i in 1:nrow(data_algG_COPD_clean_20_filter)){
  if(data_algG_COPD_clean_20_filter$SNP_base[i] == "A" & data_algG_COPD_clean_20_filter$As_COPD[i] >= 20){
    data_algG_COPD_clean_cleared <- rbind(data_algG_COPD_clean_cleared,data_algG_COPD_clean_20_filter[i,])
  } else if(data_algG_COPD_clean_20_filter$SNP_base[i] == "C" & data_algG_COPD_clean_20_filter$Cs_COPD[i] >= 20){
    data_algG_COPD_clean_cleared <- rbind(data_algG_COPD_clean_cleared,data_algG_COPD_clean_20_filter[i,])
  } else if(data_algG_COPD_clean_20_filter$SNP_base[i] == "G" & data_algG_COPD_clean_20_filter$Gs_COPD[i] >= 20){
    data_algG_COPD_clean_cleared <- rbind(data_algG_COPD_clean_cleared,data_algG_COPD_clean_20_filter[i,])
  }else if(data_algG_COPD_clean_20_filter$SNP_base[i] == "T" & data_algG_COPD_clean_20_filter$Ts_COPD[i] >= 20){
    data_algG_COPD_clean_cleared <- rbind(data_algG_COPD_clean_cleared,data_algG_COPD_clean_20_filter[i,])
  }
}

data_algG_COPD_clean <- data_algG_COPD_clean_cleared[,c(1:5)]

# Merge clean SNPs and only SNPs

ifelse(nrow(data_algG_CF_1_only_SNP) > 0,data_algG_CF_1_clean <- merge(data_algG_CF_1_clean,data_algG_CF_1_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_1","Isolates_CF_1"),all=TRUE),data_algG_CF_1_clean <- data_algG_CF_1_clean)
ifelse(nrow(data_algG_CF_2_only_SNP) > 0,data_algG_CF_2_clean <- merge(data_algG_CF_2_clean,data_algG_CF_2_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_2","Isolates_CF_2"),all=TRUE),data_algG_CF_2_clean <- data_algG_CF_2_clean)
ifelse(nrow(data_algG_CF_3_only_SNP) > 0,data_algG_CF_3_clean <- merge(data_algG_CF_3_clean,data_algG_CF_3_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_3","Isolates_CF_3"),all=TRUE),data_algG_CF_3_clean <- data_algG_CF_3_clean)
ifelse(nrow(data_algG_CF_4_only_SNP) > 0,data_algG_CF_4_clean <- merge(data_algG_CF_4_clean,data_algG_CF_4_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_4","Isolates_CF_4"),all=TRUE),data_algG_CF_4_clean <- data_algG_CF_4_clean)
ifelse(nrow(data_algG_CF_5_only_SNP) > 0,data_algG_CF_5_clean <- merge(data_algG_CF_5_clean,data_algG_CF_5_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_5","Isolates_CF_5"),all=TRUE),data_algG_CF_5_clean <- data_algG_CF_5_clean)
ifelse(nrow(data_algG_CF_6_only_SNP) > 0,data_algG_CF_6_clean <- merge(data_algG_CF_6_clean,data_algG_CF_6_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_6","Isolates_CF_6"),all=TRUE),data_algG_CF_6_clean <- data_algG_CF_6_clean)
ifelse(nrow(data_algG_CF_7_only_SNP) > 0,data_algG_CF_7_clean <- merge(data_algG_CF_7_clean,data_algG_CF_7_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_7","Isolates_CF_7"),all=TRUE),data_algG_CF_7_clean <- data_algG_CF_7_clean)
ifelse(nrow(data_algG_CF_8_only_SNP) > 0,data_algG_CF_8_clean <- merge(data_algG_CF_8_clean,data_algG_CF_8_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_8","Isolates_CF_8"),all=TRUE),data_algG_CF_8_clean <- data_algG_CF_8_clean)
ifelse(nrow(data_algG_CF_9_only_SNP) > 0,data_algG_CF_9_clean <- merge(data_algG_CF_9_clean,data_algG_CF_9_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_9","Isolates_CF_9"),all=TRUE),data_algG_CF_9_clean <- data_algG_CF_9_clean)
ifelse(nrow(data_algG_CF_10_only_SNP) > 0,data_algG_CF_10_clean <- merge(data_algG_CF_10_clean,data_algG_CF_10_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_10","Isolates_CF_10"),all=TRUE),data_algG_CF_10_clean <- data_algG_CF_10_clean)
ifelse(nrow(data_algG_CF_11_only_SNP) > 0,data_algG_CF_11_clean <- merge(data_algG_CF_11_clean,data_algG_CF_11_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_11","Isolates_CF_11"),all=TRUE),data_algG_CF_11_clean <- data_algG_CF_11_clean)
ifelse(nrow(data_algG_CF_12_only_SNP) > 0,data_algG_CF_12_clean <- merge(data_algG_CF_12_clean,data_algG_CF_12_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_12","Isolates_CF_12"),all=TRUE),data_algG_CF_12_clean <- data_algG_CF_12_clean)
ifelse(nrow(data_algG_CF_13_only_SNP) > 0,data_algG_CF_13_clean <- merge(data_algG_CF_13_clean,data_algG_CF_13_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_13","Isolates_CF_13"),all=TRUE),data_algG_CF_13_clean <- data_algG_CF_13_clean)
ifelse(nrow(data_algG_CF_14_only_SNP) > 0,data_algG_CF_14_clean <- merge(data_algG_CF_14_clean,data_algG_CF_14_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_14","Isolates_CF_14"),all=TRUE),data_algG_CF_14_clean <- data_algG_CF_14_clean)
ifelse(nrow(data_algG_CF_15_only_SNP) > 0,data_algG_CF_15_clean <- merge(data_algG_CF_15_clean,data_algG_CF_15_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_15","Isolates_CF_15"),all=TRUE),data_algG_CF_15_clean <- data_algG_CF_15_clean)
ifelse(nrow(data_algG_CF_16_only_SNP) > 0,data_algG_CF_16_clean <- merge(data_algG_CF_16_clean,data_algG_CF_16_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_16","Isolates_CF_16"),all=TRUE),data_algG_CF_16_clean <- data_algG_CF_16_clean)
ifelse(nrow(data_algG_CF_17_only_SNP) > 0,data_algG_CF_17_clean <- merge(data_algG_CF_17_clean,data_algG_CF_17_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_CF_17","Isolates_CF_17"),all=TRUE),data_algG_CF_17_clean <- data_algG_CF_17_clean)
ifelse(nrow(data_algG_acute_1_only_SNP) > 0,data_algG_acute_1_clean <- merge(data_algG_acute_1_clean,data_algG_acute_1_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_1","Isolates_acute_1"),all=TRUE),data_algG_acute_1_clean <- data_algG_acute_1_clean)
ifelse(nrow(data_algG_acute_2_only_SNP) > 0,data_algG_acute_2_clean <- merge(data_algG_acute_2_clean,data_algG_acute_2_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_2","Isolates_acute_2"),all=TRUE),data_algG_acute_2_clean <- data_algG_acute_2_clean)
ifelse(nrow(data_algG_acute_3_only_SNP) > 0,data_algG_acute_3_clean <- merge(data_algG_acute_3_clean,data_algG_acute_3_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_acute_3","Isolates_acute_3"),all=TRUE),data_algG_acute_3_clean <- data_algG_acute_3_clean)
ifelse(nrow(data_algG_env_1_only_SNP) > 0,data_algG_env_1_clean <- merge(data_algG_env_1_clean,data_algG_env_1_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_1","Isolates_env_1"),all=TRUE),data_algG_env_1_clean <- data_algG_env_1_clean)
ifelse(nrow(data_algG_env_2_only_SNP) > 0,data_algG_env_2_clean <- merge(data_algG_env_2_clean,data_algG_env_2_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_2","Isolates_env_2"),all=TRUE),data_algG_env_2_clean <- data_algG_env_2_clean)
ifelse(nrow(data_algG_env_3_only_SNP) > 0,data_algG_env_3_clean <- merge(data_algG_env_3_clean,data_algG_env_3_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_3","Isolates_env_3"),all=TRUE),data_algG_env_3_clean <- data_algG_env_3_clean)
ifelse(nrow(data_algG_env_4_only_SNP) > 0,data_algG_env_4_clean <- merge(data_algG_env_4_clean,data_algG_env_4_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_4","Isolates_env_4"),all=TRUE),data_algG_env_4_clean <- data_algG_env_4_clean)
ifelse(nrow(data_algG_env_5_only_SNP) > 0,data_algG_env_5_clean <- merge(data_algG_env_5_clean,data_algG_env_5_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_env_5","Isolates_env_5"),all=TRUE),data_algG_env_5_clean <- data_algG_env_5_clean)
ifelse(nrow(data_algG_COPD_only_SNP) > 0,data_algG_COPD_clean <- merge(data_algG_COPD_clean,data_algG_COPD_only_SNP_clean, by=c("PA14-Koordinate","GenBase_PA14","SNP_base","SNP_percentage_COPD","Isolates_COPD"),all=TRUE),data_algG_COPD_clean <- data_algG_COPD_clean)


# Merge data & Reformat I

algG_merged <- merge(data_algG_CF_1_clean,data_algG_CF_2_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat II

algG_merged <- merge(algG_merged,data_algG_CF_3_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0


# Merge data & Reformat III

algG_merged <- merge(algG_merged,data_algG_CF_4_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat IV

algG_merged <- merge(algG_merged,data_algG_CF_5_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat V

algG_merged <- merge(algG_merged,data_algG_CF_6_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat VI

algG_merged <- merge(algG_merged,data_algG_CF_7_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat VII

algG_merged <- merge(algG_merged,data_algG_CF_8_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat VIII

algG_merged <- merge(algG_merged,data_algG_CF_9_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat IX

algG_merged <- merge(algG_merged,data_algG_CF_10_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat X

algG_merged <- merge(algG_merged,data_algG_CF_11_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XI

algG_merged <- merge(algG_merged,data_algG_CF_12_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XII

algG_merged <- merge(algG_merged,data_algG_CF_13_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XIII

algG_merged <- merge(algG_merged,data_algG_CF_14_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XIV

algG_merged <- merge(algG_merged,data_algG_CF_15_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XV

algG_merged <- merge(algG_merged,data_algG_CF_16_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XVI

algG_merged <- merge(algG_merged,data_algG_CF_17_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XVII

algG_merged <- merge(algG_merged,data_algG_acute_1_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XVIII

algG_merged <- merge(algG_merged,data_algG_acute_2_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XIX

algG_merged <- merge(algG_merged,data_algG_acute_3_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XX

algG_merged <- merge(algG_merged,data_algG_env_1_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XXI

algG_merged <- merge(algG_merged,data_algG_env_2_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XXII

algG_merged <- merge(algG_merged,data_algG_env_3_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XXIII

algG_merged <- merge(algG_merged,data_algG_env_4_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XXIV

algG_merged <- merge(algG_merged,data_algG_env_5_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Merge data & Reformat XXV

algG_merged <- merge(algG_merged,data_algG_COPD_clean,by=c("PA14-Koordinate","GenBase_PA14","SNP_base"),all=TRUE)
colnames(algG_merged)[2] <- "GenBase_PA14"
colnames(algG_merged)[3] <- "SNP_base"
algG_merged[is.na(algG_merged)] <- 0

# Sum Isolates for CF and non-CF

algG_merged <- add_column(algG_merged, isolates_CF = 0, .after = "SNP_base")
algG_merged <- add_column(algG_merged, isolates_non_CF = 0, .after = "isolates_CF")


algG_merged$isolates_CF <- algG_merged$Isolates_CF_1+algG_merged$Isolates_CF_2+algG_merged$Isolates_CF_3+algG_merged$Isolates_CF_4+algG_merged$Isolates_CF_5+algG_merged$Isolates_CF_6+algG_merged$Isolates_CF_7+algG_merged$Isolates_CF_8+algG_merged$Isolates_CF_9+algG_merged$Isolates_CF_10+algG_merged$Isolates_CF_11+algG_merged$Isolates_CF_12+algG_merged$Isolates_CF_13+algG_merged$Isolates_CF_14+algG_merged$Isolates_CF_15+algG_merged$Isolates_CF_16+algG_merged$Isolates_CF_17
algG_merged$isolates_non_CF <- +algG_merged$Isolates_acute_1+algG_merged$Isolates_acute_2+algG_merged$Isolates_acute_3+algG_merged$Isolates_env_1+algG_merged$Isolates_env_2+algG_merged$Isolates_env_3+algG_merged$Isolates_env_4+algG_merged$Isolates_env_5+algG_merged$Isolates_COPD

algG_merged_switched <- algG_merged


algG_merged_switched$GenBase_PA14 <- apply(algG_merged_switched[,2,drop=F], 1,function(i) if(i == "A") i <- "T" else if(i == "C") i <- "G" else if(i == "G") i <- "C" else if(i == "T") i <- "A")
algG_merged_switched$SNP_base <- apply(algG_merged_switched[,3,drop=F], 1,function(i) if(i == "A") i <- "T" else if(i == "C") i <- "G" else if(i == "G") i <- "C" else if(i == "T") i <- "A")


save(algG_merged_switched,algG_merged,file="algG_200_20.Rdata")

# Write Output

write.xlsx(algG_merged,"algG_all_200_20.xlsx")
write.table(algG_merged$`PA14-Koordinate`,"SNP_positions_200_20.txt",row.names = FALSE,col.names = FALSE)

if(strand == "rev"){
       write.table(algG_merged_switched$GenBase_PA14,"SNP_reference_base_200_20.txt",row.names = FALSE,col.names = FALSE)
  } else{
       write.table(algG_merged$GenBase_PA14,"SNP_reference_base_200_20.txt",row.names = FALSE,col.names = FALSE)
  }
      

if(strand=="rev"){
       write.table(algG_merged_switched$SNP_base,"SNP_SNP_base_200_20.txt",row.names = FALSE,col.names = FALSE)
  }else{
       write.table(algG_merged$SNP_base,"SNP_SNP_base_200_20.txt",row.names = FALSE,col.names = FALSE)
}
