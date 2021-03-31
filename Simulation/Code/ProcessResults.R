# Process the results to generate the figures
# set to the folder containing the results.
setwd('Result/')

require(reshape)
require(ggplot2)
require(cowplot)
###################################################################################################
# Combine the results

sample.nos <- c(25, 50, 100, 250)
feature.nos <- c(100, 500, 10000, 50000)
sig.densities <- c('Low', 'Medium', 'High')
sig.strengths <- c('Weak', 'Moderate', 'Strong')
conf.densities <- c('Low', 'Medium', 'High')
conf.strengths <- c('Weak', 'Moderate', 'Strong')
conf.sig.cors <- c('Low', 'Medium', 'High')
conf.sig.locs <- c('Random', 'NonCoLoc', 'CoLoc')
cor.structs <- c('Indep', 'Block1',  'AR1')
methods <- c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
measures <- c('FDR', 'Power')

nIter <- 100
nIter2 <- 1000
load('DimZ_res.Rdata')
res.DimZ.a <- array(NA, c(1, 1, 1, 1, 
				length(conf.densities), length(conf.strengths), length(conf.sig.cors), 1, 1,
				length(methods), length(measures), nIter), 
		dimnames=list(SampleNumber='100', FeatureNumber='10000',  SignalDensity='Medium', 
				SignalStrength='Moderate', ConfounderDensity=conf.densities,
				ConfounderStrength=conf.strengths, ConfounderSignalCor=conf.sig.cors, ConfounderSignalLoc='Random', 
				CorStruct = 'Indep', Method=methods, Measure=measures, Iteration=paste(1:nIter)))

for (i in 1:length(res)) {
	res.DimZ.a['100', '10000',  c('Medium'),  c('Moderate'), , , ,'Random', 'Indep', , , i] <- res[[i]]
}


load('Base3_res.Rdata')
res.Base3.a <- array(NA, c(1, 1, 1, 1, 1, length(conf.strengths), length(conf.sig.cors), 1, 1,
				length(methods), length(measures), nIter), 
		dimnames=list(SampleNumber='100', FeatureNumber='10000',  SignalDensity='Medium', 
				SignalStrength='Moderate', ConfounderDensity='High',
				ConfounderStrength=conf.strengths, ConfounderSignalCor=conf.sig.cors, ConfounderSignalLoc='Random', 
				CorStruct = 'Indep', Method=methods, Measure=measures, Iteration=paste(1:nIter)))
for (i in 1:length(res)) {
	res.Base3.a['100', '10000',  c('Medium'),  c('Moderate'), 'High', , , 'Random', 'Indep', , , i] <- res[[i]]
}

load('P_res.Rdata')
res.P.a <- array(NA, c(1, 2, 1, 1, length(conf.densities), length(conf.strengths), length(conf.sig.cors), 1, 1,
				length(methods), length(measures), nIter2), 
		dimnames=list(SampleNumber='100', FeatureNumber= c('100', '500'),  SignalDensity='Medium', 
				SignalStrength='Moderate', ConfounderDensity=conf.densities,
				ConfounderStrength=conf.strengths, ConfounderSignalCor=conf.sig.cors, ConfounderSignalLoc='Random', 
				CorStruct = 'Indep', Method=methods, Measure=measures, Iteration=paste(1:nIter2)))

for (i in 1:length(res)) {
	res.P.a['100', c('100', '500') , c('Medium'), 'Moderate', , , , c('Random'), 'Indep', , , i] <- res[[i]]
}


load('NP_res.Rdata')
res.NP.a <- array(NA, c(3, 3, 1, 1, 
				1, 1, length(conf.sig.cors), 1, 1,
				length(methods), length(measures), nIter2), 
		dimnames=list(SampleNumber=c('25', '100', '250'), FeatureNumber=c('100', '10000', '50000'),  SignalDensity='Medium', 
				SignalStrength='Moderate', ConfounderDensity='Medium',
				ConfounderStrength='Moderate', ConfounderSignalCor=conf.sig.cors, ConfounderSignalLoc='Random', 
				CorStruct = 'Indep', Method=methods, Measure=measures, Iteration=paste(1:nIter2)))
for (i in 1:length(res)) {
	res.NP.a[,  , c('Medium'), 'Moderate', 'Medium', 'Moderate',  , c('Random'), 'Indep', , , i] <- res[[i]]
}


load('Cor1_res.Rdata')
res.Cor1.a <- array(NA, c(1, 1, 1, 1, 
				length(conf.densities), length(conf.strengths), length(conf.sig.cors), 1, 1,
				length(methods), length(measures), nIter), 
		dimnames=list(SampleNumber='100', FeatureNumber='10000',  SignalDensity='Medium', 
				SignalStrength='Moderate', ConfounderDensity=conf.densities,
				ConfounderStrength=conf.strengths, ConfounderSignalCor=conf.sig.cors, ConfounderSignalLoc='Random', 
				CorStruct = 'Block1', Method=methods, Measure=measures, Iteration=paste(1:nIter)))

for (i in 1:length(res)) {
	res.Cor1.a[c('100'), c('10000') , c('Medium'), 'Moderate', c('Low', 'Medium', 'High') , , , c('Random'), c('Block1'), , , i] <- res[[i]]
}



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
create.data <- function (res.a) {
	# Calculate the standard error
	res.df <- melt(res.a)
	colnames(res.df)[ncol(res.df)] <- 'Value'
# Error counts
	error.info <- aggregate(Value ~ SampleNumber + FeatureNumber + SignalDensity + SignalStrength + ConfounderDensity +
					ConfounderStrength + ConfounderSignalCor + ConfounderSignalLoc + CorStruct +  Method + Measure, res.df, function(x) mean(is.na(x)))
	
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
	res.df2$SampleNumber <- factor(res.df2$SampleNumber)
	res.df2$FeatureNumber <- factor(res.df2$FeatureNumber)
	res.df2$SampleNumber <- factor(res.df2$SampleNumber, levels=sample.nos[sample.nos %in% levels(res.df2$SampleNumber)])
	res.df2$FeatureNumber <- factor(res.df2$FeatureNumber, levels=feature.nos[feature.nos %in% levels(res.df2$FeatureNumber)])
	res.df2$SignalDensity <- factor(res.df2$SignalDensity, levels=sig.densities[sig.densities %in% levels(res.df2$SignalDensity)])
	res.df2$SignalStrength <- factor(res.df2$SignalStrength, levels=sig.strengths[sig.strengths %in% levels(res.df2$SignalStrength)])
	res.df2$ConfounderDensity <- factor(res.df2$ConfounderDensity, levels=conf.densities[conf.densities %in% levels(res.df2$ConfounderDensity)])
	res.df2$ConfounderStrength <- factor(res.df2$ConfounderStrength, levels=conf.strengths[conf.strengths %in% levels(res.df2$ConfounderStrength)])
	res.df2$ConfounderSignalCor <- factor(res.df2$ConfounderSignalCor, levels=conf.sig.cors[conf.sig.cors %in% levels(res.df2$ConfounderSignalCor)])
	res.df2$ConfounderSignalLoc <- factor(res.df2$ConfounderSignalLoc, levels=conf.sig.locs[conf.sig.locs %in% levels(res.df2$ConfounderSignalLoc)] )
	res.df2$CorStruct <- factor(res.df2$CorStruct, levels=cor.structs[cor.structs %in% levels(res.df2$CorStruct)])
	res.df2$Method <- factor(res.df2$Method, levels=methods[methods %in% levels(res.df2$Method)])
	return(res.df2)
}


###################################################################################################
# Plot - Annonation, theme
algs <- c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
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


###################################################################################################

res.df2 <- create.data(res.a)
# Methods to be plotted and their orders
library(grid)
library(gtable)
obj.list <- list()
i <- 1
for (Mea in c('FDR', 'Power')) {
	cat(".")
	
	if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
	if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
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
	
	if (Mea == 'Power') {
		dat1 <- res3[res3$Method == '2dFDR', ]
		dat2 <- res3[res3$Method == '1dFDR-A', ]
		labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
		Value <- pmax(dat1$Value, dat2$Value)
		dat1$labels <- labels
		dat1$Value <- Value
		
		obj <- obj + 	geom_text(data = dat1,
				aes(x = ConfounderSignalCor,  y = Value + 0.075, 
						label = labels), 
				color="black", position=position_dodge(.9), hjust=.5)
	}
	
	
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
			xlab("Confounding level") +
			scale_fill_manual(values=cols[algs.ann[algs2]]) +
			theme_bw(base_size = 22) +
			theme(legend.position="bottom") +
			theme(legend.title = element_blank())
	

	
	obj.list[[i]] <- obj
	i <- i + 1
	print(obj) 

}

pdf('Figure2.pdf', width = 15, height = 9)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


###################################################################################################
# Methods to be plotted and their orders

obj.list <- list()

i <- 1

for (Mea in c('FDR', 'Power')) {

	cat(".")
	
	if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
	if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
	res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
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
	
	if (Mea == 'Power') {
		dat1 <- res3[res3$Method == '2dFDR', ]
		dat2 <- res3[res3$Method == '1dFDR-A', ]
		labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
		Value <- pmax(dat1$Value, dat2$Value)
		dat1$labels <- labels
		dat1$Value <- Value
		
		obj <- obj + 	geom_text(data = dat1,
				aes(x = ConfounderSignalCor,  y = Value + 0.075, 
						label = labels), 
				color="black", position=position_dodge(.9), hjust=.5)
	}
	
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
			xlab("Confounding level") +
			scale_fill_manual(values=cols[algs.ann[algs2]]) +
			theme_bw(base_size = 22) +
			theme(legend.position="bottom")+
			theme(legend.title = element_blank())
	
	
	obj.list[[i]] <- obj
	i <- i + 1
	print(obj) 

}

pdf('Figure3.pdf', width = 15, height = 9)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()

############################################################################################################
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (N in c( '50', '25')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
		if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
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
		
		if (Mea == 'Power') {
			dat1 <- res3[res3$Method == '2dFDR', ]
			dat2 <- res3[res3$Method == '1dFDR-A', ]
			labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
			Value <- pmax(dat1$Value, dat2$Value)
			dat1$labels <- labels
			dat1$Value <- Value
			
			obj <- obj + 	geom_text(data = dat1,
					aes(x = ConfounderSignalCor,  y = Value + 0.075, 
							label = labels), 
					color="black", position=position_dodge(.9), hjust=.5)
		}
		
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
				xlab("Confounding level") +
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

pdf('FigureS6.pdf', width = 20, height = 18)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], obj.list[[3]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


############################################################################################################

obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (C in c('NonCoLoc', 'CoLoc')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
		if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
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
		
		if (Mea == 'Power') {
			dat1 <- res3[res3$Method == '2dFDR', ]
			dat2 <- res3[res3$Method == '1dFDR-A', ]
			labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
			Value <- pmax(dat1$Value, dat2$Value)
			dat1$labels <- labels
			dat1$Value <- Value
			
			obj <- obj + 	geom_text(data = dat1,
					aes(x = ConfounderSignalCor,  y = Value + 0.075, 
							label = labels), 
					color="black", position=position_dodge(.9), hjust=.5)
		}
		
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
				xlab("Confounding level") +
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

pdf('FigureS3.pdf', width = 20, height = 18)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], obj.list[[3]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()



############################################################################################################
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (C in c( 'Block1')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
		if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% C, drop=TRUE)	
		
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
		
		if (Mea == 'Power') {
			dat1 <- res3[res3$Method == '2dFDR', ]
			dat2 <- res3[res3$Method == '1dFDR-A', ]
			labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
			Value <- pmax(dat1$Value, dat2$Value)
			dat1$labels <- labels
			dat1$Value <- Value
			
			obj <- obj + 	geom_text(data = dat1,
					aes(x = ConfounderSignalCor,  y = Value + 0.075, 
							label = labels), 
					color="black", position=position_dodge(.9), hjust=.5)
		}
		
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
				xlab("Confounding level") +
				ggtitle(paste0('CorStruct = ', C)) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom") +
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 

	}
}

pdf('FigureS4.pdf', width = 15, height = 9)
obj <- plot_grid( obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()
############################################################################################################
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (C in c( 'AR1')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
		if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% '10000' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% C, drop=TRUE)	
		
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
		
		if (Mea == 'Power') {
			dat1 <- res3[res3$Method == '2dFDR', ]
			dat2 <- res3[res3$Method == '1dFDR-A', ]
			labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
			Value <- pmax(dat1$Value, dat2$Value)
			dat1$labels <- labels
			dat1$Value <- Value
			
			obj <- obj + 	geom_text(data = dat1,
					aes(x = ConfounderSignalCor,  y = Value + 0.075, 
							label = labels), 
					color="black", position=position_dodge(.9), hjust=.5)
		}
		
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
				xlab("Confounding level") +
				ggtitle(paste0('CorStruct = ', C)) +
				scale_fill_manual(values=cols[algs.ann[algs2]]) +
				theme_bw(base_size = 22) +
				theme(legend.position="bottom") +
				theme(legend.title = element_blank())
		
		obj.list[[i]] <- obj
		i <- i + 1
		print(obj) 

	}
}

pdf('FigureS5.pdf', width = 15, height = 9)
obj <- plot_grid( obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


############################################################################################################
res.df2 <- create.data(res.NP.a)
i <- 1
for (Mea in c('FDR', 'Power')) {

	cat(".")
	
	if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
	if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
	res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & ConfounderDensity %in% 'Medium' & ConfounderStrength %in% 'Moderate' &
					SignalStrength %in% 'Moderate' & SignalDensity %in% 'Medium' & ConfounderSignalLoc %in% 'Random' & CorStruct %in% 'Indep', drop=TRUE)	
	
	levels(res3$Method) <- algs.ann[levels(res3$Method)]
	levels(res3$SignalDensity) <- sig.density.ann[levels(res3$SignalDensity)]
	levels(res3$SignalStrength) <- sig.strength.ann[levels(res3$SignalStrength)]
	levels(res3$ConfounderDensity) <- conf.density.ann[levels(res3$ConfounderDensity)]
	levels(res3$ConfounderStrength) <- conf.strength.ann[levels(res3$ConfounderStrength)]
	levels(res3$ConfounderSignalCor) <- conf.sig.cor.ann[levels(res3$ConfounderSignalCor)]
	levels(res3$ConfounderSignalLoc) <- conf.sig.loc.ann[levels(res3$ConfounderSignalLoc)]
	
	levels(res3$SampleNumber) <- paste0('n=', levels(res3$SampleNumber))
	levels(res3$FeatureNumber) <- paste0('m=', levels(res3$FeatureNumber))
	

	dodge <- position_dodge(width=0.9)
	obj <- ggplot(res3, aes(x=ConfounderSignalCor, y=Value,  fill=Method)) +
			geom_bar(stat='identity', position=dodge) +
			geom_errorbar(aes(ymax=ymax, ymin=ymin), position=dodge, width=0.25)
	
	if (Mea == 'Power') {
		dat1 <- res3[res3$Method == '2dFDR', ]
		dat2 <- res3[res3$Method == '1dFDR-A', ]
		labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
		Value <- pmax(dat1$Value, dat2$Value)
		dat1$labels <- labels
		dat1$Value <- Value
		
		obj <- obj + 	geom_text(data = dat1,
				aes(x = ConfounderSignalCor,  y = Value + 0.075, 
						label = labels), 
				color="black", position=position_dodge(.9), hjust=.5)
	}
	
	
	if (Mea == 'FDR') {
		obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
				ylab('False Discovery Proportion') +
				facet_grid(SampleNumber ~ FeatureNumber)
	} else {
		obj <- obj + 
				ylab('True Positive Rate') +
				facet_grid(SampleNumber ~ FeatureNumber, scales='free')
	}
	
	obj <- obj +
			xlab("Confounding level") +
			scale_fill_manual(values=cols[algs.ann[algs2]]) +
			theme_bw(base_size = 22) +
			theme(legend.position="bottom") +
			theme(legend.title = element_blank())
	
	obj.list[[i]] <- obj
	i <- i + 1
	print(obj) 

}

pdf('FigureS8.pdf', width = 15, height = 9)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()



############################################################################################################
res.df2 <- create.data(res.P.a)
obj.list <- list()

i <- 1
for (Mea in c('FDR', 'Power')) {
	for (P in c( '500', '100')) {
		cat(".")
		
		if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
		if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
		res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
						SampleNumber %in% '100' & FeatureNumber %in% P & ConfounderSignalLoc %in% 'Random' & CorStruct %in% 'Indep', drop=TRUE)	
		
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
		
		if (Mea == 'Power') {
			dat1 <- res3[res3$Method == '2dFDR', ]
			dat2 <- res3[res3$Method == '1dFDR-A', ]
			labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
			Value <- pmax(dat1$Value, dat2$Value)
			dat1$labels <- labels
			dat1$Value <- Value
			
			obj <- obj + 	geom_text(data = dat1,
					aes(x = ConfounderSignalCor,  y = Value + 0.075, 
							label = labels), 
					color="black", position=position_dodge(.9), hjust=.5)
		}
		
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
				xlab("Confounding level") +
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

pdf('FigureS7.pdf', width = 20, height = 18)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], obj.list[[3]], obj.list[[4]], labels = "AUTO", ncol = 2, align = 'h', label_size=30, label_fontface = "plain")
print(obj)
dev.off()



###################################################################################################
# Methods to be plotted and their orders
res.df2 <- create.data(res.Base3.a)
obj.list <- list()

i <- 1

for (Mea in c('FDR', 'Power')) {
	cat(".")
	
	if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
	if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
	res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
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
	
	if (Mea == 'Power') {
		dat1 <- res3[res3$Method == '2dFDR', ]
		dat2 <- res3[res3$Method == '1dFDR-A', ]
		labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
		Value <- pmax(dat1$Value, dat2$Value)
		dat1$labels <- labels
		dat1$Value <- Value
		
		obj <- obj + 	geom_text(data = dat1,
				aes(x = ConfounderSignalCor,  y = Value + 0.075, 
						label = labels), 
				color="black", position=position_dodge(.9), hjust=.5)
	}
	
	if (Mea == 'FDR') {
		obj <- obj + geom_hline(aes(yintercept=0.05), linetype=2) + 
				ylab('False Discovery Proportion') +
				facet_grid( ~ ConfounderStrength)
	} else {
		obj <- obj + 
				ylab('True Positive Rate') +
				facet_grid( ~ ConfounderStrength, scales='free')
	}
	
	obj <- obj +
			xlab("Confounding level") +
			scale_fill_manual(values=cols[algs.ann[algs2]]) +
			theme_bw(base_size = 22) +
			theme(legend.position="bottom")+
			theme(legend.title = element_blank())
	
	
	obj.list[[i]] <- obj
	i <- i + 1
	print(obj) 

}

pdf('FigureS1.pdf', width = 16, height = 6)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()


###################################################################################################
# Methods to be plotted and their orders
res.df2 <- create.data(res.DimZ.a)
obj.list <- list()

i <- 1

for (Mea in c('FDR', 'Power')) {
	cat(".")
	
	if (Mea == 'FDR') algs2 <-  c('OneStage-U', 'OneStage-A', 'TwoStage-N', 'TwoStage-T')
	if (Mea == 'Power') algs2 <-  c('OneStage-A',  'TwoStage-T')
	res3 <- subset(res.df2, Measure %in% Mea & Method %in% algs2 & SignalDensity %in% 'Medium' & SignalStrength %in% 'Moderate' &
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
	
	if (Mea == 'Power') {
		dat1 <- res3[res3$Method == '2dFDR', ]
		dat2 <- res3[res3$Method == '1dFDR-A', ]
		labels <- paste0(round((dat1$Value - dat2$Value) / dat2$Value * 100), '%')
		Value <- pmax(dat1$Value, dat2$Value)
		dat1$labels <- labels
		dat1$Value <- Value
		
		obj <- obj + 	geom_text(data = dat1,
				aes(x = ConfounderSignalCor,  y = Value + 0.075, 
						label = labels), 
				color="black", position=position_dodge(.9), hjust=.5)
	}
	
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
			xlab("Confounding level") +
			scale_fill_manual(values=cols[algs.ann[algs2]]) +
			theme_bw(base_size = 22) +
			theme(legend.position="bottom")+
			theme(legend.title = element_blank())
	
	
	obj.list[[i]] <- obj
	i <- i + 1
	print(obj) 

}

pdf('FigureS2.pdf', width = 15, height = 9)
obj <- plot_grid(obj.list[[1]], obj.list[[2]], labels = "AUTO", ncol = 2, align = 'h', label_size=40, label_fontface = "plain")
print(obj)
dev.off()




