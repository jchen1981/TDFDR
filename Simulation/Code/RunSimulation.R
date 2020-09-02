# Run simulations of all settings on the cluster
# All the simulations are contained in this file
func <- function(part, paras) {
	require(qvalue)
	require(pbivnorm)
	require(REBayes)
	require(limma)
	require(MASS)
	require(tdfdr)
	
	set.seed(part)
	sample.nos <- paras$sample.nos
	feature.nos <- paras$feature.nos
	
	sig.densities <- paras$sig.densities
	sig.strengths <- paras$sig.strengths
	
	conf.densities <- paras$conf.densities
	conf.strengths <- paras$conf.strengths
	
	conf.sig.cors <- paras$conf.sig.cors
	conf.sig.locs <- paras$conf.sig.locs
	coloc.prob <- paras$coloc.prob
	
	cor.structs <- paras$cor.structs
	nblock <- paras$nblock
	rho <- paras$rho
	
	fdr.cutoff <- paras$fdr.cutoff
	pval.cutoff <- paras$pval.cutoff
	
	resdir <- paras$resdir
	prefix <- paras$prefix
	prefix2 <- paras$prefix2
	
	resdir <- gsub('/$', '', resdir)
	source(file.path(resdir, "Simulation.R"))
	
	sink(file.path(resdir, paste(prefix2, "_",  part, ".log", sep="")))
	cat(date(), '\n')
	
	methods <- c('1dFDR-U', '1dFDR-A', '1dFDR-H', '2dFDR')
	measures <- c('FDR', 'Power')
	
	res <- array(NA, c(length(sample.nos), length(feature.nos), length(sig.densities), length(sig.strengths), 
					length(conf.densities), length(conf.strengths), length(conf.sig.cors), length(conf.sig.locs),
					length(cor.structs), length(methods), length(measures)), 
			dimnames=list(paste(sample.nos), paste(feature.nos),  sig.densities, sig.strengths, conf.densities,
					conf.strengths, conf.sig.cors, conf.sig.locs, cor.structs, methods, measures))	
	
	for (sample.no in sample.nos) {
		
		for (feature.no in feature.nos) {
			
			for (sig.density in sig.densities) {
				
				for (sig.strength in sig.strengths) {
					
					for (conf.density in conf.densities) {
						cat('*')
						for (conf.strength in conf.strengths) {
							cat('!')
							for (conf.sig.cor in conf.sig.cors) {
								cat('%')
								for (conf.sig.loc in conf.sig.locs) {
									cat('.')
									for (cor.struct in cor.structs) {
										
										sig.density0 <- switch(sig.density, Low = 0.05, Medium = 0.1, High = 0.2)
										sig.strength0 <- switch(sig.strength, Weak = 0.2, Moderate = 0.3, Strong = 0.4) * sqrt(100 / sample.no) * 0.5
										conf.density0 <- switch(conf.density, Low = 0.05, Medium = 0.1, High = 0.2)
										conf.strength0 <- switch(conf.strength, Weak = 0.2, Moderate = 0.3, Strong = 0.4) * sqrt(100 / sample.no)
										conf.sig.cor0 <- switch(conf.sig.cor, Low = 0.5, Medium = 1.25, High = 2)
										
										data <- simulate.data(n = sample.no, p = feature.no, 
												sig.density = sig.density0, sig.strength.m = sig.strength0,
												conf.density = conf.density0, conf.strength.m = conf.strength0, 
												conf.sig.cor = conf.sig.cor0, conf.sig.loc = conf.sig.loc, coloc.prob = coloc.prob,
												cor.struct = cor.struct, rho = rho, nblock = nblock)
										x <- data$x
										z <- data$z
										y <- data$y
										
										truth <- data$truth
										
										# Call methods
									    qvalue2 <- function (pvalue) {
											try.obj <- try(obj <- qvalue(pvalue, pi0.method = 'smoother'))
											if (class(try.obj) == 'try-error') {
												try.obj <- try(obj <- qvalue(pvalue, pi0.method = 'bootstrap'))
												if (class(try.obj) == 'try-error') {
													obj <- list()
													obj$qvalue <- p.adjust(pvalue, 'fdr')
												}
											} 
											return(obj)
										}
										
										res.obj <- list()
										
										obj1 <- summary(lm(y ~ x))
										pv1 <- sapply(obj1, function(x) x$coefficient[2, 4])
										res.obj[['1dFDR-U']] <- list(fdr = qvalue2(pv1)$qvalue)
										
										obj2 <- summary(lm(y ~ x + z))
										pv2 <- sapply(obj2, function(x) x$coefficient[2, 4])
										pv3 <- sapply(obj2, function(x) x$coefficient[3, 4])
										res.obj[['1dFDR-A']] <- list(fdr = qvalue2(pv2)$qvalue)
										
										pv4 <- ifelse(pv3 >= pval.cutoff, pv1, pv2)
										res.obj[['1dFDR-H']] <- list(fdr = qvalue2(pv4)$qvalue)
										
										fdr <- rep(1, feature.no)
										obj3 <- tdfdr(y, x, z, alpha = fdr.cutoff, ngrid = 50)
										
										cat(obj3$t1, ' ', obj3$t2, '\n')
										fdr[obj3$pos]  <- fdr.cutoff / 2
										res.obj[['2dFDR-T']] <- list(fdr = fdr)
										
										for (method in methods) {
											
											# Evaluation
											if (is.null(res.obj[[method]])) {
												res[paste(sample.no), paste(feature.no), sig.density, sig.strength, conf.density, conf.strength,
														conf.sig.cor, conf.sig.loc, cor.struct, method, c('FDR', 'Power')] <- NA
												next
											} 
											
											fdr <- res.obj[[method]]$fdr
											
											if (sum(fdr <= fdr.cutoff) == 0) {
												res[paste(sample.no), paste(feature.no), sig.density, sig.strength, conf.density, conf.strength,
														conf.sig.cor, conf.sig.loc, cor.struct, method, 'FDR'] <- 0
												res[paste(sample.no), paste(feature.no), sig.density, sig.strength, conf.density, conf.strength,
														conf.sig.cor, conf.sig.loc, cor.struct, method, 'Power'] <- 0
											} else {
												res[paste(sample.no), paste(feature.no), sig.density, sig.strength, conf.density, conf.strength,
														conf.sig.cor, conf.sig.loc, cor.struct, method, 'FDR'] <- mean(truth[fdr <= fdr.cutoff] == 0)
												res[paste(sample.no), paste(feature.no), sig.density, sig.strength, conf.density, conf.strength,
														conf.sig.cor, conf.sig.loc, cor.struct, method, 'Power'] <- sum(truth[fdr <= fdr.cutoff] != 0) / sum(truth != 0)
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	cat('\n', date(), '\n')
	sink()
	save(res, file=file.path(resdir, paste(prefix2, "_res",  part, ".Rdata", sep="")))
	return(res)
}

prefix <- 'Sim1'
resdir <- paste0("~/project/2dfdr/", prefix)
source(paste0("~/project/2dfdr/", prefix, "/ClusterFunc.R"))

paras.list <- list()
tempdir.list <- list()
dat.list <- list()
func.list <- list()
dat.list[[1]] <- 1:100
func.list[[1]] <- func

# Varying density and strength of the true signals, and varying confounding strength
# Figure 1d, e
paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'Base1'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- 10000
paras$sample.nos <- 100
paras$conf.densities <- c('Medium')
paras$conf.strengths <- c('Moderate')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('Random')
paras$sig.densities <- c('Low', 'Medium', 'High')
paras$sig.strengths <- c('Weak', 'Moderate', 'Strong')
paras$coloc.prob <- 0.5
paras$cor.structs <- c('Indep')
paras$nblock <- 100
paras$rho <- 0.6

paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)

# Varying density and strength of the confounding signals, and varying confounding strength
# Figure S1
paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'Base2'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- 10000
paras$sample.nos <- 100
paras$conf.densities <- c('Low', 'Medium', 'High')
paras$conf.strengths <- c('Weak', 'Moderate', 'Strong')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('Random')
paras$sig.densities <- c('Medium')
paras$sig.strengths <- c('Moderate')
paras$coloc.prob <- 0.5
paras$cor.structs <- c('Indep')
paras$nblock <- 100
paras$rho <- 0.6
paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)

# Varying density and strength of the confounding signals, and varying confounding strength
# Different degrees of colocation of the true and confounding signals
# Figure S2

paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'Loc'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- 10000
paras$sample.nos <- 100
paras$conf.densities <- c('Low', 'Medium', 'High')
paras$conf.strengths <- c('Weak', 'Moderate', 'Strong')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('NonCoLoc', 'CoLoc')
paras$sig.densities <- c('Medium')
paras$sig.strengths <- c('Moderate')
paras$coloc.prob <- 0.5
paras$cor.structs <- c('Indep')
paras$nblock <- 100
paras$rho <- 0.6

paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)

# Varying density and strength of the confounding signals, and varying confounding strength
# Block correlation structure
# Figure S3
paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'Cor1'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- c(10000)
paras$sample.nos <- c(100)
paras$conf.densities <- c('Low', 'Medium', 'High')
paras$conf.strengths <- c('Weak', 'Moderate', 'Strong')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('Random')
paras$sig.densities <- c('Medium')
paras$sig.strengths <- c('Moderate')
paras$coloc.prob <- 0.5
paras$cor.structs <- c('Block1')
paras$nblock <- 100
paras$rho <- 0.6
paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)

# Varying density and strength of the confounding signals, and varying confounding strength
# AR(1) correlation structure
# Figure S4
paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'Cor2'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- c(10000)
paras$sample.nos <- c(100)
paras$conf.densities <- c('Low', 'Medium', 'High')
paras$conf.strengths <- c('Weak', 'Moderate', 'Strong')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('Random')
paras$sig.densities <- c('Medium')
paras$sig.strengths <- c('Moderate')
paras$coloc.prob <- 0.5
paras$cor.structs <- c( 'AR1')
paras$nblock <- 1000
paras$rho <- 0.6

paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)

# Varying density and strength of the confounding signals, and varying confoudning strength
# Two smaller sample sizes
# Figure S5
paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'N'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- c(10000)
paras$sample.nos <- c(25, 50)
paras$conf.densities <- c('Low', 'Medium', 'High')
paras$conf.strengths <- c('Weak', 'Moderate', 'Strong')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('Random')
paras$sig.densities <- c('Medium')
paras$sig.strengths <- c('Moderate')
paras$coloc.prob <- 0.5
paras$cor.structs <- c('Indep')
paras$nblock <- 100
paras$rho <- 0.6
paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)

# Varying density and strength of the confounding signals, and varying confounding strength
# Two smaller feature sizes
# Figure S6
paras <- list()
paras$prefix <- prefix
paras$prefix2 <- 'P'
paras$resdir <- resdir
paras$fdr.cutoff <- 0.05
paras$pval.cutoff <- 0.05
paras$feature.nos <- c(100, 500)
paras$sample.nos <- 100
paras$conf.densities <- c('Low', 'Medium', 'High')
paras$conf.strengths <- c('Weak', 'Moderate', 'Strong')
paras$conf.sig.cors <- c('Low', 'Medium', 'High')
paras$conf.sig.locs <- c('Random')
paras$sig.densities <- c('Medium')
paras$sig.strengths <- c('Moderate')
paras$coloc.prob <- 0.5
paras$cor.structs <- c('Indep')
paras$nblock <- 100
paras$rho <- 0.6

paras.list[[paras$prefix2]] <- paras
tempdir.list[[paras$prefix2]] <- file.path(resdir, paras$prefix2)





# Run on the cluster
clsapply.plus(dat.list, func.list, paras.list, tempdir.list, resdir = resdir)
