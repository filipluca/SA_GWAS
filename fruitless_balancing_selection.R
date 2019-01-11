#Making a change#
#fru analyses####

asNumeric <- function(x) as.numeric(as.character(x))
factorsNumeric <- function(x) modifyList(x, lapply(x[, sapply(x, is.factor)], asNumeric))
library(PopGenome)
library(ggplot2)

#Coordinates####

#Coordinate start: 18,414,273..18,545,586 (r6)
#http://flybase.org/reports/FBgn0004652
#What coordinates in r5?
#r5 coordinates: 14,239,995..14,371,308 

#DGRP stats####
#2L
vcf.2L.DGRP <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dgrp_sequences/dgrp_Chr2L/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.2L.DGRP <- sliding.window.transform(vcf.2L.DGRP,1000,500,type=2)
sliding.2L.DGRP <- neutrality.stats(sliding.2L.DGRP)
neutrality.2L.DGRP <- get.neutrality(sliding.2L.DGRP)[[1]]
sliding.2L.DGRP <- linkage.stats(sliding.2L.DGRP,do.ZnS = T)
linkage.2L.DGRP <- get.linkage(sliding.2L.DGRP)[[1]]
sliding.2L.DGRP <- diversity.stats(sliding.2L.DGRP)
diversity.2L.DGRP <- get.diversity(sliding.2L.DGRP)[[1]]
stats.2L.DGRP <- cbind(neutrality.2L.DGRP,linkage.2L.DGRP,diversity.2L.DGRP)
#clean up stats dataframe
stats.2L.DGRP <- cbind(Row.names=rownames(stats.2L.DGRP),stats.2L.DGRP)
rownames(stats.2L.DGRP) <- NULL
stats.2L.DGRP <- as.data.frame(stats.2L.DGRP)
stats.2L.DGRP$Start_End <- do.call(rbind,strsplit(as.character(stats.2L.DGRP$Row.names),"-"))
stats.2L.DGRP$Start <- as.numeric(as.character(stats.2L.DGRP$Start_End[,1]))
stats.2L.DGRP$End <- stats.2L.DGRP$Start+999
stats.2L.DGRP <- factorsNumeric(stats.2L.DGRP)
stats.2L.DGRP <- stats.2L.DGRP[,c(22:23,2:3,5:6,11:17)]
stats.2L.DGRP$Chrom <- "2L"
write.table(stats.2L.DGRP,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2L.DGRP.txt",quote=F,row.names=F)

#2R
vcf.2R.DGRP <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dgrp_sequences/dgrp_Chr2R/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.2R.DGRP <- sliding.window.transform(vcf.2R.DGRP,1000,500,type=2)
sliding.2R.DGRP <- neutrality.stats(sliding.2R.DGRP)
neutrality.2R.DGRP <- get.neutrality(sliding.2R.DGRP)[[1]]
sliding.2R.DGRP <- linkage.stats(sliding.2R.DGRP,do.ZnS = T)
linkage.2R.DGRP <- get.linkage(sliding.2R.DGRP)[[1]]
sliding.2R.DGRP <- diversity.stats(sliding.2R.DGRP)
diversity.2R.DGRP <- get.diversity(sliding.2R.DGRP)[[1]]
stats.2R.DGRP <- cbind(neutrality.2R.DGRP,linkage.2R.DGRP,diversity.2R.DGRP)
#clean up stats dataframe
stats.2R.DGRP <- cbind(Row.names=rownames(stats.2R.DGRP),stats.2R.DGRP)
rownames(stats.2R.DGRP) <- NULL
stats.2R.DGRP <- as.data.frame(stats.2R.DGRP)
stats.2R.DGRP$Start_End <- do.call(rbind,strsplit(as.character(stats.2R.DGRP$Row.names),"-"))
stats.2R.DGRP$Start <- as.numeric(as.character(stats.2R.DGRP$Start_End[,1]))
stats.2R.DGRP$End <- stats.2R.DGRP$Start+999
stats.2R.DGRP <- factorsNumeric(stats.2R.DGRP)
stats.2R.DGRP <- stats.2R.DGRP[,c(22:23,2:3,5:6,11:17)]
stats.2R.DGRP$Chrom <- "2R"
write.table(stats.2R.DGRP,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2R.DGRP.txt",quote=F,row.names=F)

#3L
vcf.3L.DGRP <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dgrp_sequences/dgrp_Chr3L/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.3L.DGRP <- sliding.window.transform(vcf.3L.DGRP,1000,500,type=2)
sliding.3L.DGRP <- neutrality.stats(sliding.3L.DGRP)
neutrality.3L.DGRP <- get.neutrality(sliding.3L.DGRP)[[1]]
sliding.3L.DGRP <- linkage.stats(sliding.3L.DGRP,do.ZnS = T)
linkage.3L.DGRP <- get.linkage(sliding.3L.DGRP)[[1]]
sliding.3L.DGRP <- diversity.stats(sliding.3L.DGRP)
diversity.3L.DGRP <- get.diversity(sliding.3L.DGRP)[[1]]
stats.3L.DGRP <- cbind(neutrality.3L.DGRP,linkage.3L.DGRP,diversity.3L.DGRP)
#clean up stats dataframe
stats.3L.DGRP <- cbind(Row.names=rownames(stats.3L.DGRP),stats.3L.DGRP)
rownames(stats.3L.DGRP) <- NULL
stats.3L.DGRP <- as.data.frame(stats.3L.DGRP)
stats.3L.DGRP$Start_End <- do.call(rbind,strsplit(as.character(stats.3L.DGRP$Row.names),"-"))
stats.3L.DGRP$Start <- as.numeric(as.character(stats.3L.DGRP$Start_End[,1]))
stats.3L.DGRP$End <- stats.3L.DGRP$Start+999
stats.3L.DGRP <- factorsNumeric(stats.3L.DGRP)
stats.3L.DGRP <- stats.3L.DGRP[,c(22:23,2:3,5:6,11:17)]
stats.3L.DGRP$Chrom <- "3L"
write.table(stats.3L.DGRP,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3L.DGRP.txt",quote=F,row.names=F)

#3R
vcf.3R.DGRP <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dgrp_sequences/dgrp_Chr3R/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.3R.DGRP <- sliding.window.transform(vcf.3R.DGRP,1000,500,type=2)
sliding.3R.DGRP <- neutrality.stats(sliding.3R.DGRP)
neutrality.3R.DGRP <- get.neutrality(sliding.3R.DGRP)[[1]]
sliding.3R.DGRP <- linkage.stats(sliding.3R.DGRP)
linkage.3R.DGRP <- get.linkage(sliding.3R.DGRP)[[1]]
sliding.3R.DGRP <- diversity.stats(sliding.3R.DGRP)
diversity.3R.DGRP <- get.diversity(sliding.3R.DGRP)[[1]]
stats.3R.DGRP <- cbind(neutrality.3R.DGRP,linkage.3R.DGRP,diversity.3R.DGRP)
#clean up stats dataframe
stats.3R.DGRP <- cbind(Row.names=rownames(stats.3R.DGRP),stats.3R.DGRP)
rownames(stats.3R.DGRP) <- NULL
stats.3R.DGRP <- as.data.frame(stats.3R.DGRP)
stats.3R.DGRP$Start_End <- do.call(rbind,strsplit(as.character(stats.3R.DGRP$Row.names),"-"))
stats.3R.DGRP$Start <- as.numeric(as.character(stats.3R.DGRP$Start_End[,1]))
stats.3R.DGRP$End <- stats.3R.DGRP$Start+999
stats.3R.DGRP <- factorsNumeric(stats.3R.DGRP)
stats.3R.DGRP <- stats.3R.DGRP[,c(22:23,2:3,5:6,11:17)]
stats.3R.DGRP$Chrom <- "3R"
write.table(stats.3R.DGRP,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3R.DGRP.txt",quote=F,row.names=F)

#X
vcf.X.DGRP <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dgrp_sequences/dgrp_ChrX/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.X.DGRP <- sliding.window.transform(vcf.X.DGRP,1000,500,type=2)
sliding.X.DGRP <- neutrality.stats(sliding.X.DGRP)
neutrality.X.DGRP <- get.neutrality(sliding.X.DGRP)[[1]]
sliding.X.DGRP <- linkage.stats(sliding.X.DGRP,do.ZnS = T)
linkage.X.DGRP <- get.linkage(sliding.X.DGRP)[[1]]
sliding.X.DGRP <- diversity.stats(sliding.X.DGRP)
diversity.X.DGRP <- get.diversity(sliding.X.DGRP)[[1]]
stats.X.DGRP <- cbind(neutrality.X.DGRP,linkage.X.DGRP,diversity.X.DGRP)
#clean up stats dataframe
stats.X.DGRP <- cbind(Row.names=rownames(stats.X.DGRP),stats.X.DGRP)
rownames(stats.X.DGRP) <- NULL
stats.X.DGRP <- as.data.frame(stats.X.DGRP)
stats.X.DGRP$Start_End <- do.call(rbind,strsplit(as.character(stats.X.DGRP$Row.names),"-"))
stats.X.DGRP$Start <- as.numeric(as.character(stats.X.DGRP$Start_End[,1]))
stats.X.DGRP$End <- stats.X.DGRP$Start+999
stats.X.DGRP <- factorsNumeric(stats.X.DGRP)
stats.X.DGRP <- stats.X.DGRP[,c(22:23,2:3,5:6,11:17)]
stats.X.DGRP$Chrom <- "X"
write.table(stats.X.DGRP,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.X.DGRP.txt",quote=F,row.names=F)


stats.2L.DGRP <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2L.DGRP.txt",head=T)
stats.2R.DGRP <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2R.DGRP.txt",head=T)
stats.3L.DGRP <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3L.DGRP.txt",head=T)
stats.3R.DGRP <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3R.DGRP.txt",head=T)
stats.X.DGRP <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.X.DGRP.txt",head=T)

stats.all.DGRP <- rbind(stats.2L.DGRP,stats.2R.DGRP,stats.3L.DGRP,stats.3R.DGRP,stats.X.DGRP)

#ZI stats####
#2L
vcf.2L.ZI <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dpgp3_sequences/dpgp3_Chr2L/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.2L.ZI <- sliding.window.transform(vcf.2L.ZI,1000,500,type=2)
sliding.2L.ZI <- neutrality.stats(sliding.2L.ZI)
neutrality.2L.ZI <- get.neutrality(sliding.2L.ZI)[[1]]
sliding.2L.ZI <- linkage.stats(sliding.2L.ZI,do.ZnS = T)
linkage.2L.ZI <- get.linkage(sliding.2L.ZI)[[1]]
sliding.2L.ZI <- diversity.stats(sliding.2L.ZI)
diversity.2L.ZI <- get.diversity(sliding.2L.ZI)[[1]]
stats.2L.ZI <- cbind(neutrality.2L.ZI,linkage.2L.ZI,diversity.2L.ZI)
#clean up stats dataframe
stats.2L.ZI <- cbind(Row.names=rownames(stats.2L.ZI),stats.2L.ZI)
rownames(stats.2L.ZI) <- NULL
stats.2L.ZI <- as.data.frame(stats.2L.ZI)
stats.2L.ZI$Start_End <- do.call(rbind,strsplit(as.character(stats.2L.ZI$Row.names),"-"))
stats.2L.ZI$Start <- as.numeric(as.character(stats.2L.ZI$Start_End[,1]))
stats.2L.ZI$End <- stats.2L.ZI$Start+999
stats.2L.ZI <- factorsNumeric(stats.2L.ZI)
stats.2L.ZI <- stats.2L.ZI[,c(22:23,2:3,5:6,11:17)]
stats.2L.ZI$Chrom <- "2L"
write.table(stats.2L.ZI,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2L.ZI.txt",quote=F,row.names=F)
rm(sliding.2L.ZI)

#2R
vcf.2R.ZI <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dpgp3_sequences/dpgp3_Chr2R/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.2R.ZI <- sliding.window.transform(vcf.2R.ZI,1000,500,type=2)
sliding.2R.ZI <- neutrality.stats(sliding.2R.ZI)
neutrality.2R.ZI <- get.neutrality(sliding.2R.ZI)[[1]]
sliding.2R.ZI <- linkage.stats(sliding.2R.ZI,do.ZnS = T)
linkage.2R.ZI <- get.linkage(sliding.2R.ZI)[[1]]
sliding.2R.ZI <- diversity.stats(sliding.2R.ZI)
diversity.2R.ZI <- get.diversity(sliding.2R.ZI)[[1]]
stats.2R.ZI <- cbind(neutrality.2R.ZI,linkage.2R.ZI,diversity.2R.ZI)
#clean up stats dataframe
stats.2R.ZI <- cbind(Row.names=rownames(stats.2R.ZI),stats.2R.ZI)
rownames(stats.2R.ZI) <- NULL
stats.2R.ZI <- as.data.frame(stats.2R.ZI)
stats.2R.ZI$Start_End <- do.call(rbind,strsplit(as.character(stats.2R.ZI$Row.names),"-"))
stats.2R.ZI$Start <- as.numeric(as.character(stats.2R.ZI$Start_End[,1]))
stats.2R.ZI$End <- stats.2R.ZI$Start+999
stats.2R.ZI <- factorsNumeric(stats.2R.ZI)
stats.2R.ZI <- stats.2R.ZI[,c(22:23,2:3,5:6,11:17)]
stats.2R.ZI$Chrom <- "2R"
write.table(stats.2R.ZI,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2R.ZI.txt",quote=F,row.names=F)
rm(sliding.2R.ZI)

#3L
vcf.3L.ZI <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dpgp3_sequences/dpgp3_Chr3L/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.3L.ZI <- sliding.window.transform(vcf.3L.ZI,1000,500,type=2)
sliding.3L.ZI <- neutrality.stats(sliding.3L.ZI)
neutrality.3L.ZI <- get.neutrality(sliding.3L.ZI)[[1]]
sliding.3L.ZI <- linkage.stats(sliding.3L.ZI,do.ZnS = T)
linkage.3L.ZI <- get.linkage(sliding.3L.ZI)[[1]]
sliding.3L.ZI <- diversity.stats(sliding.3L.ZI)
diversity.3L.ZI <- get.diversity(sliding.3L.ZI)[[1]]
stats.3L.ZI <- cbind(neutrality.3L.ZI,linkage.3L.ZI,diversity.3L.ZI)
#clean up stats dataframe
stats.3L.ZI <- cbind(Row.names=rownames(stats.3L.ZI),stats.3L.ZI)
rownames(stats.3L.ZI) <- NULL
stats.3L.ZI <- as.data.frame(stats.3L.ZI)
stats.3L.ZI$Start_End <- do.call(rbind,strsplit(as.character(stats.3L.ZI$Row.names),"-"))
stats.3L.ZI$Start <- as.numeric(as.character(stats.3L.ZI$Start_End[,1]))
stats.3L.ZI$End <- stats.3L.ZI$Start+999
stats.3L.ZI <- factorsNumeric(stats.3L.ZI)
stats.3L.ZI <- stats.3L.ZI[,c(22:23,2:3,5:6,11:17)]
stats.3L.ZI$Chrom <- "3L"
write.table(stats.3L.ZI,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3L.ZI.txt",quote=F,row.names=F)
rm(sliding.3L.ZI)

#3R
vcf.3R.ZI <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dpgp3_sequences/dpgp3_Chr3R/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.3R.ZI <- sliding.window.transform(vcf.3R.ZI,1000,500,type=2)
sliding.3R.ZI <- neutrality.stats(sliding.3R.ZI)
neutrality.3R.ZI <- get.neutrality(sliding.3R.ZI)[[1]]
sliding.3R.ZI <- linkage.stats(sliding.3R.ZI)
linkage.3R.ZI <- get.linkage(sliding.3R.ZI)[[1]]
sliding.3R.ZI <- diversity.stats(sliding.3R.ZI)
diversity.3R.ZI <- get.diversity(sliding.3R.ZI)[[1]]
stats.3R.ZI <- cbind(neutrality.3R.ZI,linkage.3R.ZI,diversity.3R.ZI)
#clean up stats dataframe
stats.3R.ZI <- cbind(Row.names=rownames(stats.3R.ZI),stats.3R.ZI)
rownames(stats.3R.ZI) <- NULL
stats.3R.ZI <- as.data.frame(stats.3R.ZI)
stats.3R.ZI$Start_End <- do.call(rbind,strsplit(as.character(stats.3R.ZI$Row.names),"-"))
stats.3R.ZI$Start <- as.numeric(as.character(stats.3R.ZI$Start_End[,1]))
stats.3R.ZI$End <- stats.3R.ZI$Start+999
stats.3R.ZI <- factorsNumeric(stats.3R.ZI)
stats.3R.ZI <- stats.3R.ZI[,c(22:23,2:3,5:6,11:17)]
stats.3R.ZI$Chrom <- "3R"
write.table(stats.3R.ZI,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3R.ZI.txt",quote=F,row.names=F)
rm(sliding.3R.ZI)

#Fru
vcf.3R.ZI <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dpgp3_sequences/dpgp3_Chrfru/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.fru.ZI <- splitting.data(vcf.3R.ZI,positions=seq(18414200,18545600,250),type=2)
sliding.fru.ZI <- neutrality.stats(sliding.fru.ZI)
neutrality.fru.ZI <- get.neutrality(sliding.fru.ZI)[[1]]
sliding.fru.ZI <- linkage.stats(sliding.fru.ZI)
linkage.fru.ZI <- get.linkage(sliding.fru.ZI)[[1]]
sliding.fru.ZI <- diversity.stats(sliding.fru.ZI)
diversity.fru.ZI <- get.diversity(sliding.fru.ZI)[[1]]
stats.fru.ZI <- cbind(neutrality.fru.ZI,linkage.fru.ZI,diversity.fru.ZI)
#clean up stats dataframe
stats.fru.ZI <- cbind(Row.names=rownames(stats.fru.ZI),stats.fru.ZI)
rownames(stats.fru.ZI) <- NULL
stats.fru.ZI <- as.data.frame(stats.fru.ZI)
stats.fru.ZI$Start_End <- do.call(rbind,strsplit(as.character(stats.fru.ZI$Row.names),"-"))
stats.fru.ZI$Start <- as.numeric(as.character(stats.fru.ZI$Start_End[,1]))
stats.fru.ZI$End <- stats.fru.ZI$Start+999
stats.fru.ZI <- factorsNumeric(stats.fru.ZI)
stats.fru.ZI <- stats.fru.ZI[,c(22:23,2:3,5:6,11:17)]
stats.fru.ZI$Chrom <- "fru"
write.table(stats.fru.ZI,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.fru.ZI.txt",quote=F,row.names=F)
rm(sliding.fru.ZI)


#X
vcf.X.ZI <- readData("/Volumes/Time_Machine_Backups/nexus_originals/dpgp3_sequences/dpgp3_ChrX/vcf_snpsites/vcf_without_reference/r6/original/",include.unknown=T,format="VCF",FAST=T)
sliding.X.ZI <- sliding.window.transform(vcf.X.ZI,1000,500,type=2)
sliding.X.ZI <- neutrality.stats(sliding.X.ZI)
neutrality.X.ZI <- get.neutrality(sliding.X.ZI)[[1]]
sliding.X.ZI <- linkage.stats(sliding.X.ZI,do.ZnS = T)
linkage.X.ZI <- get.linkage(sliding.X.ZI)[[1]]
sliding.X.ZI <- diversity.stats(sliding.X.ZI)
diversity.X.ZI <- get.diversity(sliding.X.ZI)[[1]]
stats.X.ZI <- cbind(neutrality.X.ZI,linkage.X.ZI,diversity.X.ZI)
#clean up stats dataframe
stats.X.ZI <- cbind(Row.names=rownames(stats.X.ZI),stats.X.ZI)
rownames(stats.X.ZI) <- NULL
stats.X.ZI <- as.data.frame(stats.X.ZI)
stats.X.ZI$Start_End <- do.call(rbind,strsplit(as.character(stats.X.ZI$Row.names),"-"))
stats.X.ZI$Start <- as.numeric(as.character(stats.X.ZI$Start_End[,1]))
stats.X.ZI$End <- stats.X.ZI$Start+999
stats.X.ZI <- factorsNumeric(stats.X.ZI)
stats.X.ZI <- stats.X.ZI[,c(22:23,2:3,5:6,11:17)]
stats.X.ZI$Chrom <- "X"
write.table(stats.X.ZI,"~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.X.ZI.txt",quote=F,row.names=F)
rm(sliding.X.ZI)


stats.2L.ZI <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2L.ZI.txt",head=T)
stats.2R.ZI <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.2R.ZI.txt",head=T)
stats.3L.ZI <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3L.ZI.txt",head=T)
stats.3R.ZI <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.3R.ZI.txt",head=T)
stats.X.ZI <- read.table("~/Documents/data/pop_gen/polymorphism/sliding_windows/r6/original/stats.X.ZI.txt",head=T)

stats.all.ZI <- rbind(stats.2L.ZI,stats.2R.ZI,stats.3L.ZI,stats.3R.ZI,stats.X.ZI)

#D. simulans stats####




#Plots####
#DGRP####
#Tajima's D plots
stats.all.DGRP$is.fru <- ifelse(stats.all.DGRP$Start>=18414273 & stats.all.DGRP$End<=18545586 & stats.all.DGRP$Chrom=="3R",1,0)

stats.all.DGRP$is.region.of.interest <- ifelse(stats.all.DGRP$Start>=18520500 & stats.all.DGRP$End<=18521500 & stats.all.DGRP$Chrom=="3R",1,0)

fru_mean <- mean(stats.all.DGRP$Tajima.D[stats.all.DGRP$is.fru==1])
roi_mean <- mean(stats.all.DGRP$Tajima.D[stats.all.DGRP$is.region.of.interest==1])

#Histogram of Tajima's D values across the whole genome
ggplot(stats.all.DGRP,aes(Tajima.D))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=fru_mean,xend=fru_mean,y=0,yend=100000,col="black",size=1)+
  geom_segment(x=roi_mean,xend=roi_mean,y=0,yend=100000,col="red",size=1,linetype="dashed")
  
#Histogram of Tajima's D/ Kelly's Zns for fru
ggplot(subset(stats.all.DGRP,is.fru==1),aes(Tajima.D))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=roi_mean,xend=roi_mean,y=0,yend=100000,col="red",size=1,linetype="dashed")

ggplot(subset(stats.all.DGRP,is.fru==1),aes(Kelly.Z_nS))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=roi_mean,xend=roi_mean,y=0,yend=100000,col="red",size=1,linetype="dashed")


#Line chart of Tajima's D / Kelly's ZnS for fru
ggplot(subset(stats.all.DGRP,is.fru==1),aes(x=Start))+
  geom_line(aes(y=Tajima.D),col="forestgreen")+
  geom_line(aes(y=Kelly.Z_nS),col="salmon")+
  geom_line(aes(y=Tajima.D),col="darkblue")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=18521000,xend=18521000,y=-100000,yend=100000,col="red",size=1,linetype="dashed")


#ZI####
#Tajima's D plots
stats.all.ZI$is.fru <- ifelse(stats.all.ZI$Start>=18414273 & stats.all.ZI$End<=18545586 & stats.all.ZI$Chrom=="3R",1,0)

stats.all.ZI$is.region.of.interest <- ifelse(stats.all.ZI$Start>=18520500 & stats.all.ZI$End<=18521500 & stats.all.ZI$Chrom=="3R",1,0)

fru_mean <- mean(stats.all.ZI$Tajima.D[stats.all.ZI$is.fru==1])
roi_mean <- mean(stats.all.ZI$Tajima.D[stats.all.ZI$is.region.of.interest==1])
fru_mean_KellyZnS <- mean(stats.all.ZI$Kelly.Z_nS[stats.all.ZI$is.fru==1])
roi_mean_KellyZnS <- mean(stats.all.ZI$Kelly.Z_nS[stats.all.ZI$is.region.of.interest==1])

#Histogram of Tajima's D values across the whole genome
ggplot(stats.all.ZI,aes(Tajima.D))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=fru_mean,xend=fru_mean,y=0,yend=100000,col="black",size=1)+
  geom_segment(x=roi_mean,xend=roi_mean,y=0,yend=100000,col="red",size=1,linetype="dashed")

ggplot(stats.all.ZI,aes(Kelly.Z_nS))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=fru_mean_KellyZnS,xend=fru_mean_KellyZnS,y=0,yend=100000,col="black",size=1)+
  geom_segment(x=roi_mean_KellyZnS,xend=roi_mean_KellyZnS,y=0,yend=100000,col="red",size=1,linetype="dashed")

#Histogram of Tajima's D/ Kelly's Zns for fru
ggplot(subset(stats.all.ZI,is.fru==1),aes(Tajima.D))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=roi_mean,xend=roi_mean,y=0,yend=100000,col="red",size=1,linetype="dashed")

ggplot(subset(stats.all.ZI,is.fru==1),aes(Kelly.Z_nS))+
  geom_histogram(fill="grey75",col="white",bins=50)+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=fru_mean_KellyZnS,xend=roi_mean_KellyZnS,y=0,yend=100000,col="red",size=1,linetype="dashed")


#Line chart of Tajima's D / Kelly's ZnS for fru
ggplot(subset(stats.all.ZI,is.fru==1),aes(x=Start))+
  geom_line(aes(y=Tajima.D),col="forestgreen")+
  geom_line(aes(y=Kelly.Z_nS),col="salmon")+
  geom_line(aes(y=Tajima.D),col="darkblue")+
  theme_bw()+
  theme(legend.title=element_blank(),axis.text = element_text(size=25),axis.title = element_text(size=30),legend.text = element_text(size=30),legend.key.height=unit(3,"line"),legend.key.width=unit(3,"line"),legend.position="none",legend.background =element_blank())+
  geom_segment(x=18521000,xend=18521000,y=-100000,yend=100000,col="red",size=1,linetype="dashed")

