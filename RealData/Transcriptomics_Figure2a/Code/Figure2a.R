library(qvalue)
library(tdfdr)
library(edgeR)
setwd('~/Dropbox/YangLu/2020_03_08_TSFDR/Submit/RealDataAnalysis/RNASeq_Figure2a/')

dataMat <- read.delim('Data/dataMat_hbv_paired.txt',header=T,stringsAsFactors=F,check.names=F)
sampleinfo <- read.delim('Data/sampleinfo_hbv_paired.txt',header=T,stringsAsFactors=T)
normcount <- rpkm(dataMat[,-c(1:7)]+1,gene.length=dataMat$Length)
normcount_tpm <- sweep(normcount,2,apply(normcount,2,sum),'/')*10^6
idx <- which(apply(normcount_tpm,1,mean)>=1)
normcount_tpm_filter <- normcount_tpm[idx,]
normcount_tpm_filter <- preprocessCore::normalize.quantiles(normcount_tpm_filter)
normcount_tpm_filter <- log2(normcount_tpm_filter)
rownames(normcount_tpm_filter) <- dataMat$GeneName[idx]
colnames(normcount_tpm_filter) <- colnames(dataMat)[-c(1:7)]
dim(normcount_tpm_filter);normcount_tpm_filter[1:5,1:5]
# One stage procedure
obj1 <- summary(lm(t(normcount_tpm_filter) ~ sampleinfo$diagnosis + sampleinfo$race + sampleinfo$gender ))
pv1 <- sapply(obj1, function(x) x$coefficient[2, 4])
qv1 <- qvalue(pv1)$qvalue
pos1 <- qv1 <= 0.05
sum(pos1)

# Two stage procedure
## increase grid size
obj2 <- tdfdr(t(normcount_tpm_filter), sampleinfo$diagnosis, sampleinfo[,3:4], alpha = 0.05, ngrid = 100)
tsfdr.plot(obj2,file.name = 'tsfdr_ngrid100_0.1_070620',fdr.level = seq(0.01, 0.1, len = 10))
pos2 <- obj2$pos
sum(pos2)
