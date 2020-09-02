suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(SmartSVA))
library(isva)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(openxlsx)
library(minfiData)
library(CpGFilter)
library(structSSI)
library(GEOquery)
library(readxl)
library(stringr)

Args <- commandArgs()
print(Args)
experiment <- Args[6]

id_conversion <- read_xlsx("/home/bailing/projects/ewas/doc/id_conversion.xlsx",
                           sheet = "Sheet1")
ID <- id_conversion$New_ID[match(experiment, id_conversion$Raw_ID)]

#### load ref data (xReactiveProbes, ann450k, iccref) =====================
load("/home/bailing/projects/ewas/analysis/reffile/ref.RData")
ann450k$refgene_pos <- sub(";.*","",ann450k$UCSC_RefGene_Group)
ann450k$refgene_pos <- ifelse(ann450k$refgene_pos=="",
                              "Non_gene",ann450k$refgene_pos)

#### set data directory and working directory =============================
outputDir <- "/home/bailing/projects/ewas/analysis/TSFDR-data"
cmd <- paste0("mkdir -p ", outputDir, "/", experiment)
system(cmd)
baseDir <- "/home/bailing/projects/ewas/data"
setwd(paste0(outputDir, "/", experiment))
dataset <- unlist(str_split(experiment, pattern = "_"))[1]
datadir <- paste0(baseDir, "/", dataset)


beta2m <- function(beta.values){
  log2(beta.values/(1 - beta.values))
}
m2beta <- function(m.values){
  2^m.values/(2^m.values + 1)
}

#### read in data and QC for samples ======================================
if(file.exists(paste0(datadir, "/SampleSheet.csv"))){
  targets <- read.metharray.sheet(datadir, pattern="SampleSheet.csv")
  rgSet <- read.metharray.exp(targets = targets)
  ## general information
  platform <- rgSet@annotation["array"]
  num_probes <- nrow(rgSet)
  col_name <- c()
  for(i in 1:ncol(rgSet)){
    samplename <- unlist(str_split(colnames(rgSet)[i], pattern = "_"))[1]
    col_name <- c(col_name, samplename)
  }
  colnames(rgSet) <- col_name
  ## read in phenotype information
  targets <- read.xlsx(paste0(baseDir, "/targets/", experiment, ".xlsx"))
  ## match phenotype info
  rgSet <- rgSet[, targets$Sample_Name]
  ## QC for samples
  detP <- detectionP(rgSet)
  keep <- colMeans(detP) < 0.01
  if(sum(keep) != nrow(targets)){
    cat("==========", sum(!keep)," sample was filted==========\n")
    rgSet <- rgSet[,keep]
    targets <- targets[keep,]
    detP <- detP[,keep]
  }
  sample_size <- ncol(rgSet)
  # Normalization
  mSetSq <- preprocessQuantile(rgSet)
  detP <- detP[rownames(mSetSq), ]
}else{
  filename1 <- paste0(datadir, "/", dataset, "_family.soft.gz")
  filename2 <- paste0(datadir, "/", dataset, "_family.soft")
  if(file.exists(filename1)){
    gset <- getGEO(filename = filename1)
  }else if(file.exists(filename2)){
    gset <- getGEO(filename = filename2)
  }else{
    stop("This dataset doesn't have soft file!")
  }
  
  
  ## general information
  platform <- gset@gpls[[1]]@header$title
  num_probes <- nrow(Table(gset@gsms[[1]]))
  ## beta values and m values
  bVals <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                  ncol = length(gset@header$sample_id))
  rownames(bVals) <- Table(gset@gsms[[1]])$ID_REF
  colnames(bVals) <- gset@header$sample_id
  for (i in 1:length(gset@gsms)) {
    bVals[, i] <- Table(gset@gsms[[i]])[[2]]
  }
  if (all(is.na(bVals) | (bVals>=0 & bVals<=1))){
    mVals <- beta2m(bVals)
  } else {
    mVals <- bVals
    bVals <- m2beta(mVals)
  }
  
  ## P values (if exists)
  if (dim(Table(gset@gsms[[1]]))[2]>=3){
    detP <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                   ncol = length(gset@header$sample_id))
    rownames(detP) <- rownames(bVals)
    colnames(detP) <- colnames(bVals)
    for (i in 1:length(gset@gsms)) {
      detP[, i] <- Table(gset@gsms[[i]])[[3]]
    }
  }
  
  ## read in phenotype information
  targets <- read.xlsx(paste0(baseDir, "/targets/", experiment, ".xlsx"))
  ## remove duplicated samples
  index <- duplicated(targets$Rep_Design)
  targets <- targets[!index, ]
  ## match phenotype info with values(beta, M, P) and probe info in ann450k
  bVals <- bVals[is.element(rownames(bVals), ann450k$Name), targets$Sample_Name]
  mVals <- mVals[is.element(rownames(mVals), ann450k$Name), targets$Sample_Name]
  if (exists("detP")){
    detP <- detP[is.element(rownames(detP), ann450k$Name), targets$Sample_Name]
    ## QC for samples
    keep <- colMeans(detP,na.rm=TRUE) < 0.01
    if(sum(keep) != nrow(targets)){
      cat("==========", sum(!keep)," sample was filtered==========\n")
      targets <- targets[keep,]
      detP <- detP[,keep]
      bVals <- bVals[, keep]
      mVals <- mVals[, keep]
    }
  }
  sample_size <- ncol(bVals)
  ## convert to GenomicRatioSet class
  mSetSq <- makeGenomicRatioSetFromMatrix(mat = bVals, pData = targets, what = "Beta")
  if(exists("detP")){detP <- detP[rownames(mSetSq), ]}
}

#### QC for probes ========================================================
if (exists("detP")){
  keep <- rowSums(detP < 0.01,na.rm=TRUE) >= ncol(detP) * 0.5
  cat(paste(rep("=",10),collapse = ""),sum(keep==FALSE)," probes have bad quality ",
      paste(rep("=",10),collapse = ""),"\n")
  mSetSqFlt <- mSetSq[keep,]
} else {
  mSetSqFlt <- mSetSq
}

#filter out the probes from the X and Y chromosomes
keep <- !(featureNames(mSetSqFlt) %in%
            ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
#filter out probes that have shown to be cross-reactive
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]
#filter out the probes where common SNPs may affect the CpG or single base extension (SBE) site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))
# the probes' number of experiment after filter
num_probes_filter <- nrow(mSetSqFlt)


#### get beta and M values =========================================
bVals <- getBeta(mSetSqFlt)

#Remove the lines with NA values
bVals <- na.omit(bVals)

#Adjust the beta values of 0/1
if(sum(bVals == 0) > 0){
  bVals[bVals == 0] <- min(bVals[bVals != 0])
}
if(sum(bVals == 1) > 0){
  bVals[bVals == 1] <- max(bVals[bVals != 1])
}

# M values
mVals <- beta2m(bVals)
index <- duplicated(targets$Rep_Design)
targets <- targets[!index, ]
mVals <- mVals[, targets$Sample_Name]

#### SVA ==================================================================
mod <- model.matrix( ~ targets$Sample_Group)

#for m value
Y.r <- t(resid(lm(t(mVals) ~ targets$Sample_Group)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
sva.res <- smartsva.cpp(mVals, mod, mod0=NULL, n.sv=n.sv)


#### save objects =============================================
sv <- sva.res$sv
saveRDS(sv, file = paste0(ID, ".sv.RDS"))
saveRDS(mVals, file = paste0(ID, ".mVals.RDS"))
saveRDS(targets, file = paste0(ID, ".targets.RDS"))

#### tdfdr ====================================================
library(tdfdr)

if((!is.numeric(targets$Sample_Group)) & nlevels(as.factor(targets$Sample_Group)) > 2){
  cat('This data set has multiple-level groups !')
}else{
  # y - matrix， Methylation data (M value), row: sample, column: CpG
  # x - vector,  Numerical values (Two-level factor, please convert into 0 and 1)
  # currently,  x can not be multi-level factor
  # z - matrix,  Surrogate variables.
  y <- t(mVals)
  z <- sva.res$sv
  if(is.numeric(targets$Sample_Group)){
    x <- targets$Sample_Group
  }else{
    x <- model.matrix(~targets$Sample_Group)[, 2]
  }
  
  cat('Perform TS FDR ...\n')
  cat(date())
  obj1 <- tdfdr(y, x, z, alpha = 0.20)
  cat(date())
  obj2 <- fastLM(t(y), cbind(1, x, z))
  cat(date())
  cat('Finished!\n')
  save(obj1, file = paste0(ID, '.obj1.RDS'))
  save(obj2, file = paste0(ID, ".obj2.RDS"))
  tdfdr.plot(obj1, file.name = paste0(ID, '.tdfdr.plot'))
}
