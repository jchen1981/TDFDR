# Process the results to generate the figures
# set to the folder containing the results.
setwd('Result/')

require(reshape)
require(ggplot2)
require(cowplot)
###################################################################################################
# Combine the results

sample.nos <- c(25, 50, 100)
feature.nos <- c(100, 500, 10000)
sig.densities <- c('Low', 'Medium', 'High')
sig.strengths <- c('Weak', 'Moderate', 'Strong')
conf.densities <- c('Low', 'Medium', 'High')
conf.strengths <- c('Weak', 'Moderate', 'Strong')
conf.sig.cors <- c('Low', 'Medium', 'High')
conf.sig.locs <- c('Random', 'NonCoLoc', 'CoLoc')
cor.structs <- c('Indep', 'Block1',  'AR1')
methods <- c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
measures <- c('FDR', 'Power')

nIter <- 100

res.a <- array(NA, c(length(sample.nos), length(feature.nos), length(sig.densities), length(sig.strengths), 
				length(conf.densities), length(conf.strengths), length(conf.sig.cors), length(conf.sig.locs), length(cor.structs),
				length(methods), length(measures), nIter), 
		dimnames=list(SampleNumber=paste(sample.nos), FeatureNumber=paste(feature.nos),  SignalDensity=sig.densities, 
				SignalStrength=sig.strengths, ConfounderDensity=conf.densities,
				ConfounderStrength=conf.strengths, ConfounderSignalCor=conf.sig.cors, ConfounderSignalLoc=conf.sig.locs, 
				CorStruct = cor.structs, Method=methods, Measure=measures, Iteration=paste(1:nIter)))	

load('Base1_res.Rdata')

for (i in 1:length(res)) {
	res.a['100', '10000', , , c('Medium'),  c('Moderate'), , 'Random', 'Indep', , , i] <- res[[i]]
}

load('Base2_res.Rdata')

for (i in 1:length(res)) {
	res.a['100', '10000',  c('Medium'),  c('Moderate'), , , , 'Random', 'Indep', , , i] <- res[[i]]
}

load('Loc_res.Rdata')

for (i in 1:length(res)) {
	res.a['100', '10000' , c('Medium'), 'Moderate', , , , c('NonCoLoc', 'CoLoc'), 'Indep', , , i] <- res[[i]]
}

load('P_res.Rdata')

for (i in 1:length(res)) {
	res.a['100', c('100', '500') , c('Medium'), 'Moderate', , , , c('Random'), 'Indep', , , i] <- res[[i]]
}

load('N_res.Rdata')

for (i in 1:length(res)) {
	res.a[c('25', '50'), c('10000') , c('Medium'), 'Moderate', , , , c('Random'), 'Indep', , , i] <- res[[i]]
}

load('Cor1_res.Rdata')

for (i in 1:length(res)) {
	res.a[c('100'), c('10000') , c('Medium'), 'Moderate', c('Low', 'Medium', 'High') , , , c('Random'), c('Block1'), , , i] <- res[[i]]
}

load('Cor2_res.Rdata')

for (i in 1:length(res)) {
	res.a[c('100'), c('10000') , c('Medium'), 'Moderate', c('Low', 'Medium', 'High') , , , c('Random'), c('AR1'), , , i] <- res[[i]]
}


###################################################################################################
# Calculate the standard error
res.df <- melt(res.a)
colnames(res.df)[ncol(res.df)] <- 'Value'
# Error counts
error.info <- aggregate(Value ~ SampleNumber + FeatureNumber + SignalDensity + SignalStrength + ConfounderDensity +
				ConfounderStrength + ConfounderSignalCor + ConfounderSignalLoc + CorStruct +  Method + Measure, res.df, function(x) mean(is.na(x)))
write.csv(error.info, 'NumberOfAlgFailure.csv')

m  <- aggregate(Value ~ SampleNumber + FeatureNumber + SignalDensity + SignalStrength + ConfounderDensity +
				ConfounderStrength + ConfounderSignalCor + ConfounderSignalLoc + CorStruct + Method + Measure, res.df, function(x) mean(x[!is.na(x)]))
se <- aggregate(Value ~ SampleNumber + FeatureNumber + SignalDensity + SignalStrength + ConfounderDensity +
				ConfounderStrength + ConfounderSignalCor + ConfounderSignalLoc + CorStruct + Method + Measure, res.df, function(x) {
			x <- x[!is.na(x)]
			n <- length(x)
			sd(x) / sqrt(n)
		})
			#sd(x / sqrt(length(x))})
sd <- aggregate(Value ~ SampleNumber + FeatureNumber + SignalDensity + SignalStrength + ConfounderDensity +
				ConfounderStrength + ConfounderSignalCor + ConfounderSignalLoc + CorStruct + Method + Measure, res.df, function(x) {
			ind <- !is.na(x)
			sd(x[ind])})

res.df2 <- cbind(m, SD = sd[, ncol(sd)], ymax=m[, ncol(m)] + 1.96 * se[, ncol(m)], ymin=m[, ncol(m)] - 1.96 * se[, ncol(m)],  
		SE=se[, ncol(se)], ErrRate=error.info[, ncol(m)])

res.df2$SampleNumber <- factor(res.df2$SampleNumber, levels=sample.nos)
res.df2$FeatureNumber <- factor(res.df2$FeatureNumber, levels=feature.nos)
res.df2$SignalDensity <- factor(res.df2$SignalDensity, levels=sig.densities)
res.df2$SignalStrength <- factor(res.df2$SignalStrength, levels=sig.strengths)
res.df2$ConfounderDensity <- factor(res.df2$ConfounderDensity, levels=conf.densities)
res.df2$ConfounderStrength <- factor(res.df2$ConfounderStrength, levels=conf.strengths)
res.df2$ConfounderSignalCor <- factor(res.df2$ConfounderSignalCor, levels=conf.sig.cors)
res.df2$ConfounderSignalLoc <- factor(res.df2$ConfounderSignalLoc, levels=conf.sig.locs )
res.df2$CorStruct <- factor(res.df2$CorStruct, levels=cor.structs)
res.df2$Method <- factor(res.df2$Method, levels=methods)

###################################################################################################
# Plot - Annonation, theme
algs <- c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
algs.ann <- c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
names(algs.ann) <- algs

sig.density.ann <- c('Low'='Low density', 'Medium'='Medium density', 'High'='High density')
sig.strength.ann <- c('Weak'='Weak effect', 'Moderate'='Moderate effect', 'Strong'='Strong effect')

conf.density.ann <- c('Low'='Low density', 'Medium'='Medium density', 'High'='High density')
conf.strength.ann <- c('Weak'='Weak effect', 'Moderate'='Moderate effect', 'Strong'='Strong effect')

conf.sig.cor.ann <- c('Low'='+', 'Medium'='++', 'High'='++')
conf.sig.loc.ann <- c('Random'='Random', 'NonCoLoc'='NonCoLoc', 'CoLoc'='NonCoLoc')


conf.sig.cor.ann <- c('Low'='+', 'Medium'='++', 'High'='+++')
conf.sig.loc.ann <- c('Random'='Random', 'NonCoLoc'='NonCoLoc', 'CoLoc'='NonCoLoc')


cor.struct.ann <- c('Indep'='Indep', 'Block1'='Block1', 'Block2'='Block2', 'AR1' = 'AR1')

nMeth <- length(algs)
cols <- scales::hue_pal()(nMeth)
cols <- c("#999999", "#56B4E9", "#E69F00",  'red', "#009E73",  "#0072B2", "#D55E00", "#CC79A7", "#F0E442", "orange",
		"steelblue", 'darkseagreen', 'firebrick2', 'gold2')[1:nMeth]  # Colorblind friendly
shapes <- c(rep(21:25, 3))[1:nMeth]
ltys <- rep(c(1, 2), ceiling(nMeth/2))[1:nMeth]
names(cols) <- names(shapes) <- names(ltys) <- algs.ann

d1 <- subset(res3, Method == '1dFDR-A')
d2 <- subset(res3, Method == '2dFDR')

temp <- data.frame(d2$SignalDensity, d2$SignalStrength, d2$ConfounderSignalCor, improve=(d2$Value - d1$Value) / d1$Value * 100)

###################################################################################################
library(grid)
library(gtable)
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {

		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
		if (Mea == 'Power') algs2 <-  c( '1dFDR-A',  '2dFDR')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & ConfounderDensity %in% 'Medium' & ConfounderStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% 'Indep', drop=TRUE)	
		
		levels(res3$Method) <- algs.ann[levels(res3$Method)]
		levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
		levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
		levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
		levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
		levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
		levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]


		dodge <- position_dodge(width=0.9)
		obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
				geom_bar(stat='identity', position=dodge) +
				geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
		
		if (Mea == 'FDR') {
			obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
					ylab('False Discovery Proportion') +
					facet_grid(SignalDensity ~ SignalStrength)
		} else {
			obj <- obj + 
					ylab('True Positive Rate') +
					facet_grid(SignalDensity ~ SignalStrength, scales='free')
		}
		
		obj <- obj +
				xlab("Confounding") +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom") +
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 
}

pdf('Figure1d_e.pdf', width = 15, height = 8)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


###################################################################################################
# Figure S1

obj.list <- list()

i <- 1

for (Mea in c('FDR', 'Power')) {
#	for (Strength in c( 'Weak', 'Moderate')) {
	cat(".")
	
	if (Mea == 'FDR') algs2 <-  c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
	if (Mea == 'Power') algs2 <-  c( '1dFDR-A',  '2dFDR')
	res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
					SampleNumber %in% '100' & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% 'Indep', drop=TRUE)	
	
	levels(res3$Method) <- algs.ann[levels(res3$Method)]
	levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
	levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
	levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
	levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
	levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
	levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]
	
#	pdf(paste0(Dist, "_", Mea, "_grid_bar.pdf"), width=12, height=8)
	dodge <- position_dodge(width=0.9)
	obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
			geom_bar(stat='identity', position=dodge) +
			geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
	
	if (Mea == 'FDR') {
		obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
				ylab('False Discovery Proportion') +
				facet_grid(ConfounderDensity ~ ConfounderStrength)
	} else {
		obj <- obj + 
				ylab('True Positive Rate') +
				facet_grid(ConfounderDensity ~ ConfounderStrength, scales='free')
	}
	
	obj <- obj +
			xlab("Confounding") +
			scale_fill_manual(values=cols[algs.ann[algs2]]) +
			theme_bw(base_size = 22) +
			theme(legend.position="bottom")+
			theme(legend.title = element_blank())
	

	obj.list[[i]] <- obj
	i <- i + 1
	print(obj) 
	# dev.off()
#	}
}

pdf('FigureS1.pdf', width = 15, height = 8)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()

############################################################################################################
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (C in c('NonCoLoc', 'CoLoc')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
		if (Mea == 'Power') algs2 <-  c( '1dFDR-A',  '2dFDR')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% C & CorStruct %in% 'Indep', drop=TRUE)	
		
		levels(res3$Method) <- algs.ann[levels(res3$Method)]
		levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
		levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
		levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
		levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
		levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
		levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]
		
#	pdf(paste0(Dist, "_", Mea, "_grid_bar.pdf"), width=12, height=8)
		dodge <- position_dodge(width=0.9)
		obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
				geom_bar(stat='identity', position=dodge) +
				geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
		
		if (Mea == 'FDR') {
			obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
					ylab('False Discovery Proportion') +
					facet_grid(ConfounderDensity ~ ConfounderStrength)
		} else {
			obj <- obj + 
					ylab('True Positive Rate') +
					facet_grid(ConfounderDensity ~ ConfounderStrength, scales='free')
		}
		
		obj <- obj +
				xlab("Confounding") +
				ggtitle(paste0('Relation = ', C)) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom") +
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 
		# dev.off()
	}
}

pdf('FigureS2.pdf', width = 20, height = 16)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], obj.list[[3]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()

############################################################################################################
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (C in c('Block1', 'AR1')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
		if (Mea == 'Power') algs2 <-  c( '1dFDR-A',  '2dFDR')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% C, drop=TRUE)	
		
		levels(res3$Method) <- algs.ann[levels(res3$Method)]
		levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
		levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
		levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
		levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
		levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
		levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]
		
#	pdf(paste0(Dist, "_", Mea, "_grid_bar.pdf"), width=12, height=8)
		dodge <- position_dodge(width=0.9)
		obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
				geom_bar(stat='identity', position=dodge) +
				geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
		
		if (Mea == 'FDR') {
			obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
					ylab('False Discovery Proportion') +
					facet_grid(ConfounderDensity ~ ConfounderStrength)
		} else {
			obj <- obj + 
					ylab('True Positive Rate') +
					facet_grid(ConfounderDensity ~ ConfounderStrength, scales='free')
		}
		
		obj <- obj +
				xlab("Confounding") +
				ggtitle(paste0('CorStruct = ', C)) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom") +
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 
		# dev.off()
	}
}


pdf('FigureS3.pdf', width = 15, height = 8)
obj <- plot_grid(obj.list[[1]],  obj.list[[3]],  labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()

pdf('FigureS4.pdf', width = 15, height = 8)
obj <- plot_grid( obj.list[[2]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


############################################################################################################
# The effect of sample size
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (N in c( '50', '25')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
		if (Mea == 'Power') algs2 <-  c( '1dFDR-A',  '2dFDR')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% N & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% 'Indep', drop=TRUE)	
		
		levels(res3$Method) <- algs.ann[levels(res3$Method)]
		levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
		levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
		levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
		levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
		levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
		levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]
		
#	pdf(paste0(Dist, "_", Mea, "_grid_bar.pdf"), width=12, height=8)
		dodge <- position_dodge(width=0.9)
		obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
				geom_bar(stat='identity', position=dodge) +
				geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
		
		if (Mea == 'FDR') {
			obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
					ylab('False Discovery Proportion') +
					facet_grid(ConfounderDensity ~ ConfounderStrength)
		} else {
			obj <- obj + 
					ylab('True Positive Rate') +
					facet_grid(ConfounderDensity ~ ConfounderStrength, scales='free')
		}
		
		obj <- obj +
				xlab("Confounding") +
				ggtitle(paste0('N=', N)) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom")+
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 
		# dev.off()
	}
}

pdf('FigureS5.pdf', width = 20, height = 16)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], obj.list[[3]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()

############################################################################################################
# The effect of feature size
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (P in c( '500', '100')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
		if (Mea == 'Power') algs2 <-  c( '1dFDR-A',  '2dFDR')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% P & ConfounderSignalLoc %in% 'Random' & CorStruct %in% 'Indep', drop=TRUE)	
		
		levels(res3$Method) <- algs.ann[levels(res3$Method)]
		levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
		levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
		levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
		levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
		levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
		levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]
		
#	pdf(paste0(Dist, "_", Mea, "_grid_bar.pdf"), width=12, height=8)
		dodge <- position_dodge(width=0.9)
		obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
				geom_bar(stat='identity', position=dodge) +
				geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
		
		if (Mea == 'FDR') {
			obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
					ylab('False Discovery Proportion') +
					facet_grid(ConfounderDensity ~ ConfounderStrength)
		} else {
			obj <- obj + 
					ylab('True Positive Rate') +
					facet_grid(ConfounderDensity ~ ConfounderStrength, scales='free')
		}
		
		obj <- obj +
				xlab("Confounding") +
				ggtitle(paste0('P=', P)) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom")+
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 
		# dev.off()
	}
}

pdf('FigureS6.pdf', width = 20, height = 16)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], obj.list[[3]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=30, label_fontface = "plain")
print(obj)
dev.off()





