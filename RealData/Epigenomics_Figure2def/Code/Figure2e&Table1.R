library(reshape2)
library(ggpubr)
library(dplyr)
setwd('~/Dropbox/YangLu/2020_03_08_TSFDR/Submit/RealDataAnalysis/EWAS_Figure2def/')

files = list.files(path='Data/Step1/.')
file0 = list.files(path='Data/Step1/.',pattern = '*.obj1.RDS')
file0 = gsub('.obj1.RDS','',file0)
file0 = file0[grep('^EWAS',file0)]

adj.RSquare = n.svs =  RSquare = list()
for (i in 1: length(file0)){
  tryCatch({
  name = file0[i]
  ewas.sv = readRDS(paste0('Data/Step1/',name,'.sv.RDS'))
  n.sv = ncol(ewas.sv)
  n.svs[name] = n.sv
  ewas.targets = readRDS(paste0('Data/Step1/',name,'.targets.RDS'))
  if(length(unique(ewas.targets$Sample_Group)) > 2){
    ewas = cbind(ewas.targets$Sample_Group, ewas.sv) %>% as.data.frame()
  }else{
    ewas = cbind(as.factor(ewas.targets$Sample_Group), ewas.sv) %>% as.data.frame()
    ewas[,1] = ewas[,1] -1
  }
  ewas.sum = summary(lm(V1 ~ ., data = ewas))
  adj.RSquare[name] = ewas.sum$adj.r.squared
  RSquare[name] = ewas.sum$r.squared
  }, error = function(e){cat("ERROR :",name,conditionMessage(e), "\n")})
}

adj.RSquares = melt(adj.RSquare)
colnames(adj.RSquares) = c('adj.RSquare','New_ID')
RSquares = melt(RSquare)
colnames(RSquares) = c('RSquare','New_ID')
svs = melt(n.svs)
colnames(svs) = c('nSVs','New_ID')


#### create DataSetInfo
getwd()
meta = read.csv('Data/Step2/DatasetInfo.csv')
findings = read.csv('Data/Step2/comp.csv')
colnames(findings)[1] = colnames(meta)[1]
DataSetInfo = inner_join(findings, svs) %>% 
  inner_join(adj.RSquares) %>% 
  inner_join(RSquares) %>% 
  inner_join(meta) %>% mutate(improveP = ((X2dFDR-`X1dFDR.A`)/(`X1dFDR.A`+1)))
DataSetInfo_prep = DataSetInfo %>% dplyr::select(-c('Replicates','Sample_Group'))
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'New_ID'] = 'ID'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'X2dFDR'] = '# TDFDR hits'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'X1dFDR.A'] = '# OSFDR hits'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'nSVs'] = '# SVs'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'improveP'] = '% improvement *'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'RSquare'] = 'R2'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'GEO_Datasets'] = 'GEO accesion'
colnames(DataSetInfo_prep)[colnames(DataSetInfo_prep) == 'adj.RSquare'] = 'Adjusted R2'

DataSetInfo_prep  = DataSetInfo_prep %>% dplyr::select(c('ID','# TDFDR hits','# OSFDR hits','% improvement *',
                                                         '# SVs','Adjusted R2','R2','GEO accesion','Phenotype',
                                                         'PMID','Platform','Tissue','Size'))

write.csv(DataSetInfo_prep,'SupplementaryTable1.csv',row.names = F)

# Figure 2e
DataSetInfo1 = DataSetInfo %>% mutate(improveP = 100*((X2dFDR-`X1dFDR.A`)/(`X1dFDR.A`+1)))
DataSetInfo1[DataSetInfo1$improveP > 100,'improveP'] = 100
hist_improvement = ggplot(DataSetInfo1, aes(improveP)) +
  geom_histogram(fill ='#36989a', color = 'white', binwidth = 10) +
  theme_bw()+
  geom_vline(xintercept = -1, color = 'red', size=0.8) +
  xlab('Improvement, %') +
  ylab('Datasets') +
  theme(axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 20, color = 'black'),
        plot.margin = margin(1, 1, 1, 1, "cm"))
hist_improvement
# ggsave('Fig2e.pdf', width = 5, height = 5, dpi = 100)
