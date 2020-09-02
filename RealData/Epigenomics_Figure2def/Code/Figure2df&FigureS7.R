pkg = c('dplyr','tidyr','tibble','ggpubr','ggplot2','qvalue','reshape2','RColorBrewer','plyr')
lapply(pkg, require, character=T)
setwd('~/Dropbox/YangLu/2020_03_08_TSFDR/Submit/RealDataAnalysis/EWAS_Figure2def/')
########################  Data arrangement ########################
# =========== loading EWAS results of tsfdr and classic methods
files = list.files('Data/Step1', pattern = '*obj1.RDS')
qval.sum = pos.sum = data.frame()
for (file in files){
  ewas = readRDS(paste0('Data/Step1/',file))
  pos <- tdfdr.select(ewas, fdr.level = 0.05)$pos %>% as.data.frame() %>% rownames_to_column('id')
  qval <- qvalue(ewas$p.value,fdr.level=0.05)$qvalue %>% as.data.frame() %>% rownames_to_column('id')

  colnames(pos)[2] <- paste0('tsfdr_',gsub('.obj1.RDS','',file))
  colnames(qval)[2] <- paste0('qval_',gsub('.obj1.RDS','',file))
  
  if (ncol(pos.sum)==0){
    pos.sum <- pos
  }else{
    pos.sum <- full_join(pos.sum, pos)
  }

  if (ncol(qval.sum)==0){
    qval.sum <- qval
  }else{
    qval.sum <- full_join(qval.sum, qval)
  }
}
saveRDS(qval.sum,'Data/Step2/qval.sum.Rdata')
saveRDS(pos.sum,'Data/Step2/pos.sum.Rdata')

# =========== loading saved RDS result of tsfdr and classic methods 
pos.sum = readRDS('Data/Step2/pos.sum.Rdata')
qval.sum = readRDS('Data/Step2/qval.sum.Rdata')
colnames(qval.sum);colnames(pos.sum)
# one case example of EWAS1 dataset
ewas1 = qval.sum[,c('id',"qval_EWAS1")] %>% full_join(pos.sum[,c('id',"tsfdr_EWAS1")])
tsfdr.sig = ewas1[(ewas1$qval_EWAS1 > 0.05 & ewas1$tsfdr_EWAS1 == TRUE),]
classic.sig = ewas1[(ewas1$qval_EWAS1 <= 0.05 & ewas1$tsfdr_EWAS1 == FALSE),]
tsfdr.classic.sig = ewas1[(ewas1$qval_EWAS1 <= 0.05 & ewas1$tsfdr_EWAS1 == TRUE),]
tsfdr.classic.nonsig = ewas1[(ewas1$qval_EWAS1 > 0.05 & ewas1$tsfdr_EWAS1 == FALSE),]
nrow(tsfdr.sig) + nrow(classic.sig) + nrow(tsfdr.classic.sig) + nrow(tsfdr.classic.nonsig); nrow(ewas1)

########################  boxplot compare performance of 2 methods ########################
# =========== 1. compare number of positives on log scale with boxplot
qval.sum[is.na(qval.sum)] = 1
classic = tsfdr = list()
for (i in 2:ncol(qval.sum)){
  name1 = colnames(qval.sum)[i] 
  name2 = colnames(pos.sum)[i] 
  classic[[name1]] = sum(ifelse(qval.sum[,i] <= 0.05, TRUE, FALSE))
  tsfdr[[name2]] = sum(pos.sum[,i][!is.na(pos.sum[,i])])
}
classic = melt(classic)
tsfdr = melt(tsfdr)
head(tsfdr);head(classic)
tsfdr$L1 = gsub('tsfdr_\\.*','',tsfdr$L1)
classic$L1 = gsub('qval_\\.*','',classic$L1)
colnames(tsfdr)[1] = '2dFDR'
colnames(classic)[1] = '1dFDR-A'
comp = inner_join(tsfdr, classic);nrow(classic); nrow(tsfdr);nrow(comp)
comp = comp[,c(2,1,3)]
colnames(comp)[1] = 'DataSet'
# write.csv(comp,'Data/Step2/comp.csv',row.names = F)
median = plyr::ddply(melt(comp), c("variable"), summarise,
                     median = median(value, na.rm = TRUE))

diff = comp %>% mutate(diff = `2dFDR` - `1dFDR-A`) 
head(diff);str(diff)
nrow(diff[diff$diff>0,])/nrow(diff)
worse = diff %>% filter(diff <= 0) # classic outperforms/equal to tsfdr
better = diff %>% filter(diff > 0) # tsfdr otuperforms classic
comp[,c(2:3)] = log(comp[,c(2:3)] +0.5) # add sudo ct to 0s
comp1 = melt(comp)
comp1$variable = factor(comp1$variable, levels = c('1dFDR-A','2dFDR'),ordered = TRUE)
head(comp1)
p1 = ggplot(comp1, aes(x = variable, y = value, fill = factor(variable),color = factor(variable)))+
  stat_boxplot(aes(variable, value), geom='errorbar', linetype=1, width=0.2, lwd = 1)+
  geom_boxplot(width = 0.5, aes(color = factor(variable)), alpha = 0.7, lwd = 1)+
  geom_point(size=2, alpha=0.4, color = 'grey20') +
  geom_line(aes(group=DataSet), colour="grey60", alpha=0.5) +
  theme_bw()+
  scale_fill_manual(values = c(brewer.pal(8,'Dark2')[c(1,2)])) +
  scale_color_manual(values = c(brewer.pal(8,'Dark2')[c(1,2)])) +
  theme(axis.text.x = element_text(color="black", size = 18,vjust = 0.5),
        axis.text.y = element_text(color="black", size = 18),
        axis.title = element_text(color="black", size = 18),
        legend.position = 'none')+
  labs(x='',y='log (Number of DMPs)')+
  guides(fill=FALSE) 
p1
getwd()
ggsave(path = 'plots','Figure2d.pdf',width=6, height = 4, dpi = 200)
  

########################  Figure 2f: age associated ########################
meta = read.csv('Data/Step2/DatasetInfo.csv')
age = meta %>% filter(Phenotype =='Age')
age = as.vector(age$New_ID)
list.files('Data/Step1/')
age = paste0(age,'.obj1.RDS')
sub = age[age!='EWAS14.obj1.RDS']
ps = pos = pval = data.frame()
for (i in sub){
  ewas = readRDS(paste0('Data/Step1/',i))
  p = as.data.frame(ewas$p.value) %>% rownames_to_column('id')
  colnames(p)[2] = gsub('.obj1.RDS','',i)
  
  t <- tdfdr.select(ewas, fdr.level = 0.05)$pos %>% as.data.frame() %>% rownames_to_column('id') 
  colnames(t)[2] = paste0(gsub('.obj1.RDS','',i),'.tsfdr')
  
  c <- qvalue(ewas$p.value,fdr.level=0.05)$qvalue %>% as.data.frame() %>% rownames_to_column('id') 
  colnames(c)[2] = paste0(gsub('.obj1.RDS','',i),'.classic')
  
  if (nrow(ps)==0){
    ps <- p
  }else{
    ps <- full_join(ps, p)
  }
  
  if (nrow(pos)==0){
    pos <- t
  }else{
    pos <- full_join(pos, t)
  }
  
  if (nrow(pval)==0){
    pval <- c
  }else{
    pval <- full_join(pval, c)
  }
}

bench = full_join(pos, pval)

# (1) more findings in tdfdr of each EWAS dataset
EWAS26 = bench %>% filter(EWAS26.tsfdr ==TRUE & EWAS26.classic > 0.05)
EWAS27 = bench %>% filter(EWAS27.tsfdr ==TRUE & EWAS27.classic > 0.05)
EWAS30 = bench %>% filter(EWAS30.tsfdr ==TRUE & EWAS30.classic > 0.05)
EWAS39 = bench %>% filter(EWAS39.tsfdr ==TRUE & EWAS39.classic > 0.05)
EWAS45 = bench %>% filter(EWAS45.tsfdr ==TRUE & EWAS45.classic > 0.05)

# (2) EWAS26 found in other
g = list(EWAS26,EWAS27,EWAS30,EWAS39,EWAS45)
name = c('EWAS26','EWAS27','EWAS30','EWAS39','EWAS45')
target = name[1]
name1 = name[!(name %in% target)]
EWAS26.p = ps[ps$id %in% (EWAS26[,'id']),] %>% 
  dplyr::select(c('id',name1)) %>% 
  melt() %>% dplyr::select(-id) %>% mutate(target = target) # group means ewas26 is the base, while others siginificant findings in ewas26

target = name[2]
name1 = name[!(name %in% target)]
EWAS27.p = ps[ps$id %in% (EWAS27[,'id']),] %>% 
  dplyr::select(c('id',name1)) %>% 
  melt() %>% dplyr::select(-id) %>% mutate(target = target)

target = name[3]
name1 = name[!(name %in% target)]
EWAS30.p = ps[ps$id %in% (EWAS30[,'id']),] %>% 
  dplyr::select(c('id',name1)) %>% 
  melt() %>% dplyr::select(-id) %>% mutate(target = target)

target = name[4]
name1 = name[!(name %in% target)]
EWAS39.p = ps[ps$id %in% (EWAS39[,'id']),] %>% 
  dplyr::select(c('id',name1)) %>% 
  melt() %>% dplyr::select(-id) %>% mutate(target = target)

target = name[5]
name1 = name[!(name %in% target)]
EWAS45.p = ps[ps$id %in% (EWAS45[,'id']),] %>% 
  dplyr::select(c('id',name1)) %>% 
  melt() %>% dplyr::select(-id) %>% mutate(target = target)

EWAS = rbind(EWAS26.p,EWAS27.p, EWAS30.p,EWAS39.p,EWAS45.p)

colors = c('EWAS26' = brewer.pal(8,'Set1')[c(2:5,7)][1], 'EWAS27' = brewer.pal(8,'Set1')[c(2:5,7)][2],
           'EWAS30' = brewer.pal(8,'Set1')[c(2:5,7)][3], 'EWAS39' = brewer.pal(8,'Set1')[c(2:5,7)][4],
           'EWAS45' = brewer.pal(8,'Set1')[c(2:5,7)][5])
hist = ggplot(EWAS, aes(value, fill = variable)) +
  geom_histogram(binwidth = 0.05)+
  scale_x_continuous(limit = c(0, 1), oob = function(x, limits) x) +
  facet_wrap(target ~ variable, scales="free_y", nrow = 5, ncol = 4) +
  scale_fill_manual(values = colors) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color="black", size = 12),
        axis.title = element_text(color="black", size = 12, margin = margin(t = 0, r = 0, b = 0, l = 10)),
        plot.title = element_text(size = 10),
        legend.position = 'none',
        plot.margin = unit(c(0.5,0,0,0), "cm"),
        axis.ticks.length=unit(-0.07, "cm")) +
  labs(x='',y='Counts') +
  theme(strip.background = element_blank(),
        strip.text = element_blank())

ps1 = melt(ps) %>% filter(variable!= 'EWAS14')
density = ggplot(ps1) +
  geom_density(aes(x=value, fill = variable,color = variable)) +
  theme_bw() +
  facet_grid(variable~ ., scales="free_y") +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  theme(panel.background = element_rect(fill = 'grey95', color = 'grey95'),
        axis.text.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color="black", size = 12),
        axis.title = element_text(color="black", size = 12,margin = unit(c(0, 0, 0, 0), "cm")),
        legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0,0), "cm"),
        axis.ticks.length=unit(-0.07, "cm"))+
  labs(x='',y='Density')

p = ggarrange(hist,density, widths = c(4,1), heights = c(1,1))
p1= annotate_figure(p,bottom = text_grob("p-value", color = "black",hjust = 1.2, x =0.55,  size = 12))
getwd()
# ggsave('Figure2f.pdf',width = 10, height = 8,dpi = 200)


# ========== 2. SLE associated : Figure S7
meta = read.csv('Data/Step2/DatasetInfo.csv')
sub = table(meta$Phenotype) %>% as.data.frame() %>% filter(Freq > 1)
n = sub$Var1
select = n[8]
pheno = meta %>% filter(Phenotype == select)
pheno = as.vector(pheno$New_ID)
list.files('data')
sub = paste0(pheno,'.obj1.RDS')

ewas = readRDS(paste0("Data/Step1/",sub[1]))
p = as.data.frame(ewas$p.value) %>% rownames_to_column('id')
colnames(p)[2] = gsub('.obj1.RDS','',sub[1])

t <- tdfdr.select(ewas, fdr.level = 0.05)$pos %>% as.data.frame() %>% rownames_to_column('id') 
colnames(t)[2] = paste0(gsub('.obj1.RDS','',sub[1]),'.tsfdr')

c <- qvalue(ewas$p.value,fdr.level=0.05)$qvalue %>% as.data.frame() %>% rownames_to_column('id') 
colnames(c)[2] = paste0(gsub('.obj1.RDS','',sub[1]),'.classic')
ewasA = inner_join(p,c) %>% inner_join(t)


ewas = readRDS(paste0("Data/Step1/",sub[2]))
p = as.data.frame(ewas$p.value) %>% rownames_to_column('id')
colnames(p)[2] = gsub('.obj1.RDS','',sub[2])

t <- tdfdr.select(ewas, fdr.level = 0.05)$pos %>% as.data.frame() %>% rownames_to_column('id') 
colnames(t)[2] = paste0(gsub('.obj1.RDS','',sub[2]),'.tsfdr')

c <- qvalue(ewas$p.value,fdr.level=0.05)$qvalue %>% as.data.frame() %>% rownames_to_column('id') 
colnames(c)[2] = paste0(gsub('.obj1.RDS','',sub[2]),'.classic')
ewasB = inner_join(p,c) %>% inner_join(t)
bench = full_join(ewasA, ewasB)

# (1) more findings in tsfdr
EWASA = bench[bench[,4]==TRUE & bench[,3] > 0.05,] %>% na.omit() #more finidngs in tsfdr in A
EWASB = bench[bench[,7]==TRUE & bench[,6] > 0.05,] %>% na.omit() #more finidngs in tsfdr in B

# (2) tsfdr more findings in EWAS28's distribution in EWAS29 
name = colnames(EWASA)
histAinB = ggplot(EWASA, aes(!!as.name(name[5]))) +
  geom_histogram(color = brewer.pal(8,'Set2')[2],fill = brewer.pal(8,'Set2')[2]) +
  scale_x_continuous(limit = c(0, 1), oob = function(x, limits) x) +
  theme_bw()+
  theme(#axis.text.x=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color="black", size = 12),
    axis.title = element_text(color="black", size = 12, margin = margin(t = 0, r = 0, b = 0, l = 10)),
    plot.title = element_text(size = 10),
    legend.position = 'none',
    plot.margin = unit(c(0.5,0,0,0), "cm"),
    axis.ticks.length=unit(-0.07, "cm"))+
  labs(x='',y='Counts')+
  theme(strip.background = element_blank(),
        strip.text = element_blank())
histBinA = ggplot(EWASB, aes(!!as.name(name[2]))) +
  geom_histogram(color = brewer.pal(8,'Set2')[1],fill = brewer.pal(8,'Set2')[1]) +
  scale_x_continuous(limit = c(0, 1), oob = function(x, limits) x) +
  theme_bw() +
  theme(#axis.text.x=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color="black", size = 12),
    axis.title = element_text(color="black", size = 12, margin = margin(t = 0, r = 0, b = 0, l = 10)),
    plot.title = element_text(size = 10),
    legend.position = 'none',
    plot.margin = unit(c(0.5,0,0,0), "cm"),
    axis.ticks.length=unit(-0.07, "cm"))+
  labs(x='',y='Counts')+
  theme(strip.background = element_blank(),
        strip.text = element_blank())
hist = ggarrange(histAinB, histBinA, nrow = 2)
densityA_B = bench %>% dplyr::select(name[c(1:2,5)]) %>% melt()
density = ggplot(densityA_B) +
  geom_density(aes(x=value, fill = variable,color = variable))+
  theme_bw()+
  facet_grid(variable~ ., scales="free_y")+
  scale_fill_manual(values = (brewer.pal(8,'Set2')[c(1:2)]))+
  scale_color_manual(values = (brewer.pal(8,'Set2')[c(1:2)]))+
  theme(axis.text.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color="black", size = 12),
        axis.title = element_text(color="black", size = 12,
                                  margin = unit(c(0, 0, 0, 0), "cm")),
        legend.position = 'none',
        plot.margin = unit(c(0.5,0.5,0,0), "cm"),
        axis.ticks.length=unit(-0.07, "cm"))+
  labs(x='',y='Density')
SLE_grid = ggarrange(hist, density, nrow = 1)
SLE_grid1 = annotate_figure(SLE_grid,bottom = text_grob("p-value", color = "black",hjust = 1.2, x =0.55,  size = 12))

ggsave('FigureS7.pdf',width = 6, height =6,dpi = 200)
  




