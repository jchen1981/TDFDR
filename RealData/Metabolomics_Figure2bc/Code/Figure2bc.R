library(tdfdr); library(dplyr)

setwd('~/Dropbox/YangLu/2020_03_08_TSFDR/Submit/RealDataAnalysis/Metabolics_Figure2bc/')
# loading data
pheno <- read.table('Data/phenotypes.tab', sep = '\t', head = TRUE, row.names = 1);head(pheno);dim(pheno)
metabolic <- read.table('Data/metabolomic.tab', sep = '\t', head = TRUE, row.names = 1) 
lipidomic <- read.table('Data/lipidomic.tab', sep = '\t', head = TRUE, row.names = 1)
metabolome <- inner_join(metabolic%>% rownames_to_column('id'), lipidomic%>% rownames_to_column('id')) %>% column_to_rownames('id')
dim(metabolome);dim(lipidomic);dim(metabolic)

# filter all NAs, most occurred in Homa.IR
idx <- rownames(pheno[!is.na(pheno$Homa.IR),]);length(idx)
metabolic1 <- metabolic[(rownames(metabolic)%in% idx),] %>% as.matrix();dim(metabolic1);dim(metabolic)
lipidomic1 <- lipidomic[(rownames(lipidomic)%in% idx),] %>% as.matrix();dim(lipidomic1);dim(lipidomic)
metabolome1 <- metabolome[(rownames(metabolome)%in% idx),] %>% as.matrix() 

# distribution plots
par(mfrow=c(1,3))
hist(t(metabolome1));hist(metabolic1);hist(lipidomic1)

# replace 0 with (0.5 *minimal non-zero value)
replace_halfmin <- function(df){
  for(i in 1:nrow(df)){
    idx = which(df[i,] == 0)
    df[i, idx] = min(df[i, df[i,] > 0] * 0.5) 
  }
  return(df)
}

metabolic1 <- replace_halfmin(metabolic1)
lipidomic1 <- replace_halfmin(lipidomic1)
metabolome1 <- replace_halfmin(metabolome1)

# log2 transformation
metabolic1 <- log2(metabolic1);dim(metabolome1)
lipidomic1 <- log2(lipidomic1);dim(metabolome1)
metabolome1 <- log2(metabolome1);dim(metabolome1)
pheno1 <- pheno[!is.na(pheno$Homa.IR),];dim(pheno1)
table(pheno$Diabetes);table(pheno1$Diabetes)
cor.test(pheno1$Homa.IR, pheno1$BMI.kg.m2, method = 'spearman')

par(mfrow=c(1,3))
hist(t(metabolome1));hist(metabolic1);hist(lipidomic1)

#  tdfdr analysis for Insulin resistance metabolomics dataset pooling polar metabolites and molecular lipids
tdfdr.obj <- tdfdr(y = metabolome1, x = pheno1$Homa.IR, z = pheno1$BMI.kg.m2, alpha = 0.05)
TSFDR <- tdfdr.select(tdfdr.obj, fdr.level = 0.05)$pos %>% as.data.frame() %>% rownames_to_column('metabolic')
tdfdr.plot(tdfdr.obj,fdr.level = seq(0.01, 0.1, len=10),file.name = "Figure2b")

# tdfdr analysis for Insulin resistance metabolomics dataset using polar metabolites only
tdfdr.met <- tdfdr(y = metabolic1, x = pheno1$Homa.IR, z = pheno1$BMI.kg.m2, alpha = 0.05)
TDFDR.met <- tdfdr.select(tdfdr.met, fdr.level = 0.05)$pos %>% as.data.frame() %>% rownames_to_column('metabolic')
tdfdr.plot(tdfdr.met,fdr.level = seq(0.01, 0.1, len=10),file.name = "Figure2c")

