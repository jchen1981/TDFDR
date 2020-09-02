pkg = c('tdfdr','qvalue','reshape2','ggplot2')
sapply(pkg, require, character = TRUE)

setwd('~/Documents/Mayo_Research/2020_03_08_TSFDR/TSFDR_submit/EWAS/')
ewas22.mVals = readRDS('Data/Step3/EWAS22.mVals.RDS')
ewas22.targets = readRDS('Data/Step3/EWAS22.targets.RDS')
ewas22.sv = readRDS('Data/Step3/EWAS22.sv.RDS')
grp1 = (nrow(ewas22.targets[ewas22.targets$Sample_Group ==unique(ewas22.targets$Sample_Group)[1],]))/nrow(ewas22.targets)
grp2 = (nrow(ewas22.targets[ewas22.targets$Sample_Group ==unique(ewas22.targets$Sample_Group)[2],]))/nrow(ewas22.targets)
ratio = grp1/grp2

#### tdfdr 
y <- t(ewas22.mVals)
z <- ewas22.sv;is.data.frame(z)
x <- model.matrix(~ ewas22.targets$Sample_Group)[, 2]
obj0 <- tdfdr(t(ewas22.mVals), x, z, alpha = 0.05)
# saveRDS(obj0, 'ewas22_full.RData')
ewas22_select = tdfdr.select(obj0, fdr.level = 0.05)
P.value = obj0$p.value
ewas22_full.bonferroni = sum(p.adjust(obj0$p.value, method = 'bonferroni') <= 0.05)
GoldStandard = obj0$p.value[p.adjust(obj0$p.value, method = 'bonferroni') <= 0.05]
saveRDS(GoldStandard,'GoldStandard.RData')
write.csv(GoldStandard,'GoldStandard.csv')
GoldStandard = readRDS('GoldStandard.RData')

#### downsampling at different depth
sums = list()
sum = as.data.frame(matrix(0, ncol = 3, nrow = 100))
colnames(sum) = c('tdfdr', 'Bonferroni', 'qvalue')
for (sample_size in c(20,40,60,80,100)){
  for (iter in 1:100){
    tryCatch({
      cat(date())
      group1 = ewas22.targets[ewas22.targets$Sample_Group == unique(ewas22.targets$Sample_Group)[1],]
      group2 = ewas22.targets[ewas22.targets$Sample_Group == unique(ewas22.targets$Sample_Group)[2],]
      
      group1_index = which(ewas22.targets$Sample_Name %in% group1$Sample_Name)
      group2_index = which(ewas22.targets$Sample_Name %in% group2$Sample_Name)
      
      sample1 = sample(group1_index, round(sample_size * grp1))
      sample2 = sample(group2_index, round(sample_size * grp2))
      
      targets_subset <- ewas22.targets[c(sample1, sample2), ]
      mVals_subset <- ewas22.mVals[, targets_subset$Sample_Name]
      sv_subset <- ewas22.sv[c(sample1,sample2),]
      
      y <- t(mVals_subset);y[1:5,1:5]
      z <- sv_subset;is.data.frame(z);z[1:5,1:5]
      x <- model.matrix(~ targets_subset$Sample_Group)[, 2]
      obj <- tdfdr(y, x, z, alpha = 0.05)
      # tdfdr_select = tdfdr.select(obj, fdr.level = 0.05)
      
      # saveRDS(obj, paste0('iter',iter,'sample_size',sample_size,'.RData'))
      tdfdr = sum(names(obj$pos[obj$pos == TRUE]) %in% names(GoldStandard))
      bonferroni = sum(names(obj$p.value[p.adjust(obj$p.value, method = 'bonferroni') <= 0.05]) %in% names(GoldStandard))
      qvalue = sum(names(obj$p.value[qvalue::qvalue(obj$p.value)$qvalue <= 0.05]) %in% names(GoldStandard)) 
      
      sum[iter,] = c(tdfdr,bonferroni,qvalue)
      rownames(sum)[iter] = paste0('res',iter)

    }, error = function(e){cat("ERROR :",'iter:',iter,'sample_size:',sample_size,conditionMessage(e), "\n")})
  }
  sums[[paste0(sample_size)]] = sum
  saveRDS(sums, paste0('sample_size_',sample_size,'_ewas22.RData'))
  ewas22.res = melt(sums)
  write.csv(ewas22.res, paste0('sample_size_',sample_size,'_ewas22.csv'))
}

# saveRDS(sums, 'stratified_res_ewas22.RData')

##### making boxplot
size20 = read.csv('Data/Step3/sample_size_20_ewas22.csv',row.names = 1) %>% filter(variable %in% c('tdfdr','qvalue')) 
size40 = read.csv('Data/Step3/sample_size_40_ewas22.csv',row.names = 1) %>% filter(variable %in% c('tdfdr','qvalue')) 
size60 = read.csv('Data/Step3/sample_size_60_ewas22.csv',row.names = 1) %>% filter(variable %in% c('tdfdr','qvalue')) 
size80 = read.csv('Data/Step3/sample_size_80_ewas22.csv',row.names = 1) %>% filter(variable %in% c('tdfdr','qvalue')) 
size100 = read.csv('Data/Step3/sample_size_100_ewas22.csv',row.names = 1) %>% filter(variable %in% c('tdfdr','qvalue'))
size = rbind(size20, size40, size60, size80, size100)
size$variable = as.factor(size$variable)
size$L1 <- as.factor(size$L1)
colnames(size) = c('measure','Power','sample_size')
size$Power <- size$Power/10
cdata <- plyr::ddply(size, c("measure", "sample_size"), summarise,
                     N    = sum(!is.na(Power)),
                     mean = mean(Power, na.rm = TRUE),
                     sd   = sd(Power, na.rm = TRUE),
                     se   = sd / sqrt(N))
cdata$measure = gsub(unique(cdata$measure)[1],'1dFDR-A',cdata$measure)
cdata$measure = gsub(unique(cdata$measure)[2],'2dFDR',cdata$measure)

boxplt = ggplot(data = cdata, aes(x = sample_size, y = mean, color = measure)) +
  geom_point(position = position_dodge(width = 0.3), size = 3, shape = 18) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width = 0.3,
                position = position_dodge(width = 0.3),
                lwd = 0.7) +
  scale_color_brewer(palette = 'Dark2') +
  xlab("Sample size") + ylab("Percentage of gDMPs recovered") +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
ggsave('SuplementaryNote5.pdf', width = 7, height = 5, dpi = 100)
