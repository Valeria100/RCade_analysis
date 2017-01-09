#----------------------------------------------------------------------------------------------------------------------------------
# Integrate CHIPSEQ DATA - apply_rcade --------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1 - Read the differential expression of the microarray data ----
VCap_DE <- read.table(paste(dirname(getwd()),"/Bayesian_acf/VCaPoutput.tsv",sep=""), header=TRUE)[,c(1,4,5,6)]

# 2 - Get the corresponding EnsemblID ----
library(DBI)
library("illuminaHumanv3.db")
ExpAnnotation_VCap<-illuminaHumanv3fullReannotation()
annotation_vcap1 <- ExpAnnotation_VCap[which(ExpAnnotation_VCap$IlluminaID%in%VCap_DE$Probe),]
annotation_vcap2 <- annotation_vcap1[order(match(annotation_vcap1$IlluminaID,VCap_DE$Probe)),]
VCaP_DE1 <- cbind(VCap_DE$Probe,annotation_vcap2$EnsemblReannotated,VCap_DE[,-1],annotation_vcap2$GenomicLocation)
colnames(VCaP_DE1) <- c("Probe","EnsemblID","Symbol","logFC","B","GenomicLocation")
rm(annotation_vcap1);rm(annotation_vcap2)
head(VCaP_DE1)
dim(VCaP_DE1)
#     Separate the genomic location into multiple rows
library(tidyr)
VCaP_DE<-separate_rows(VCaP_DE1, GenomicLocation, sep=",")
dim(VCaP_DE)
# 3 - DElookup table ----
DElookup <- list(GeneID="Ensembl.Gene.ID", logFC="logFC", B="B")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 4- Targets information - A matrix containing information about the .Bam files to be used in the analysis ----

# samples_file <- read.csv(paste(dirname(getwd()),"/PEAKS_analysis/samples_file.csv",sep=""), 
#                          header=TRUE, colClasses = c(rep(NA,6),rep("NULL",11))) 

samples_file <- read.csv("/Volumes/Valeria/CRUK/AR_all_folder/AR_final_code/Dora_AR/samples.csv", 
                         header=TRUE, colClasses = c(NA,"NULL","NULL","NULL","NULL","NULL",NA,NA,NA,"NULL","NULL",NA,NA)) 

fileID <- as.character(samples_file$SRA)
sampleID <- samples_file$Experiment
factorID <- samples_file$Sample.Input
levels(factorID) <- c("Input", "S")

# 4a- Order the files' name in the same order as in the table ----
bam_bai <- dir("/Volumes/Valeria/CRUK/AR_all_folder/AR_final_code/Dora_AR/bam", pattern = ".bam")
file_dir <- bam_bai[-c(grep("bai", bam_bai),grep("samstat", bam_bai))]

file_dir_name <- sapply(strsplit(file_dir,"_"), function(i) i[1])

#x nell'ordine di y --> x[order(match(x,y))]
filepath <- file_dir[order(match(file_dir_name,fileID))]

targets <- data.frame(fileID,sampleID,factorID,filepath)

# 4b- Use only the samples you are interested in analysing ----
interest <- c("SRR039769", "SRR039774", "SRR039773", "SRR039775")
targets_int <- targets[which(targets$fileID%in%interest),]
# write.table(targets_int, file="/Volumes/Valeria/CRUK/AR_downstream_analysis/EncodeData/targets_chipseq_peaks_id.txt", 
#             quote=FALSE, row.names=FALSE, col.names=FALSE)
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 5- Annotation information ----
EnsemblID <- as.character(VCaP_DE$EnsemblID)
chr <- sapply(strsplit(VCaP_DE$GenomicLocation,":"), function(i) i[1])
start <- sapply(strsplit(VCaP_DE$GenomicLocation,":"), function(i) i[2])
end <- sapply(strsplit(VCaP_DE$GenomicLocation,":"), function(i) i[3])
# str <- sapply(strsplit(VCaP_DE$GenomicLocation,":"), function(i) i[4])
genomic_location <- cbind(EnsemblID, chr, start, end, str)
head(genomic_location)
library(Rcade)
ChIPannoZones <- defineBins(genomic_location, zone=c(-5000, 5000), geneID="EnsemblID")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 6- Prior specification ----
DE.prior <- 0.01
prior.mode <- "keepChIP"
prior <- c("D|C"=0.05, "D|notC"=0.005)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 7- Analysis function ----
library(Rcade)
Dir <- "/Volumes/Valeria/CRUK/AR_all_folder/AR_final_code/Dora_AR/bam"
cl <- NULL
for(i in 1:length(DE)){
  Rcade <- RcadeAnalysis(DE[[i]], ChIPannoZones[[i]], annoZoneGeneidName = "Ensembl.Gene.ID", 
                         ChIPtargets = targets_int, ChIPfileDir = Dir, cl=cl, 
                         DE.prior=DE.prior, prior.mode=prior.mode, prior=prior, DElookup=DElookup)
  # dir.create("Rcade_files_B")
  exportRcade(Rcade, directory = paste("Rcade_files_B_",i,sep=""),  cutoffMode = "B", cutoffArg = -5, removeDuplicates = "beforeCutoff")
  # exportRcade(Rcade, directory = paste("Rcade_files_FDR",i,sep=""),  cutoffMode = "FDR", removeDuplicates = "beforeCutoff")
  # rcade_output <- getRcade(Rcade)
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read the output files
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
rcade_results1 <- vector("list",length(DE))
for(i in 1:length(DE)){
  rcade_results1[[i]] <- read.csv(paste("Rcade_files_B_",i,"/DEandChIP.csv",sep=""), header=TRUE, colClasses = c(NA,rep("NULL",11),NA)) 
}
