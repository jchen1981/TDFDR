

fastLM <- function(Y, M) {
	Y <- as.matrix(Y)
	XXI <- solve(t(M) %*% M)
	dof <- ncol(Y) - ncol(M)
	est <- XXI %*% t(M) %*% t(Y)
	resid <- t(Y) - M %*% est
	sigma <- sqrt(colSums(resid^2)/dof)
	pval <- 2 * pt(-abs(t(est/(sqrt(diag(XXI)))) / sigma), dof)
	return(list(pval = pval, sigma = sigma, dof = dof, est = est))
}


#' Two-dimensional false discovery rate control for powerful confounder adjustment in omics association studies
#'
#' The function implements the two-dimensional false discovery rate control for powerful confounder adjustment in omics-wide association studies.
#' The method is based on the idea that the confounder(s) usually affect part of the omics features, and thus adjusting the confounders(s)
#' for ALL genomic features will be over-adjustment, leading to reduced statistical power.  The proposed procedure starts with performing
#' the unadjusted analysis (first dimension - filtering) to narrow down the list of omics features which are more likely to be affected by either the confounder
#' or the variable of interest or both. In the second dimension, we conduct adjusted analysis on these 'top' candidates to reduce multiple testing burden.
#' In other words, the unadjusted p-values tell us about the probability of the null hypotheses being false, and the multiple testing can be focused on those
#' promising hypotheses. The procedure is theoretically guaranteed to control the false discovery rate while maximizing the power of discovery. 
#'  
#' @param y a matrix of numeric values for the high-dimensional data (sample size by number of features).
#' @param x a vector of numeric values for the variable of interest or a factor of two levels.
#' @param z a factor, data.frame, or a matrix of numeric values, representing the confounding variables. 
#' @param est.pi0 a logical value indicating whether John-Storey's correction (pi0 estimation) should be applied. Default TRUE. 
#' @param alpha a numeric value, the target FDR level.
#' @param etype a character string indicating the type of error control ("FDR" or "FWER"). Currently only supports "FDR".
#' @param t1a,t2a a numeric vector giving the search grid for the z-score of adjusted and unadjusted statistics, respectively. Deafult is NULL, and the grid
#' will be determined automatically.
#' @param ngrid an integer number indicating the number of grids. Default is 50. For more refined results, increase the number.
#' @param parallel a logical value indicating whether to search in parallel. Default is FALSE.
#' @param cores an integer number indicating the number of cores needed if \code{parallel} is TRUE.
#' @param verbose a logical value indicating whether to print out progress messages during computation.
#'
#' @return A list with the elements
#' \item{pos}{a vector of logical values indicating the significant findings (positives) at the given type I error level.}
#' \item{t1, t2}{a numerical value indicating the best point on the search grid.}
#' \item{t1t2}{a numeric vector containing the values of the search grid.}
#' \item{Zu, Za}{a numeric vector of unadjusted or adjusted z-scores.}
#' \item{NP}{a numeric vector of the nubmer of postives at all the grid points.}
#' \item{FDP}{a numeric vector of the false discovery proportions at all the grid points.}
#' \item{p.value}{a numeric vector of the p-values of the adjusted analyses.}
#' \item{mm}{the object from GLmix, containing the NPEB estimate of the location parameters.}
#' \item{pi0}{A numeric value of the estimated null proportion.}
#' 
#' @author ***
#' @references Two-dimensional false discovery rate control for powerful confounder adjustment in omics-wide association studies. Submitted.
#' @keywords FDR  confounder 
#' @import doMC
#' @importFrom stats model.matrix pt pnorm lm 
#' @importFrom qvalue qvalue
#' @examples
#'
#' require(qvalue)
#' truth <- c(rep(1, 50), rep(0, 50), rep(1, 50), rep(0, 850))
#' x <- rnorm(50)
#' z <- x + rnorm(50)
#' z <- scale(z)
#'
#' y1 <- x %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y2 <- z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y3 <- x %*% t(rep(0.5, 50)) + z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y <- cbind(y1, y2, y3, matrix(rnorm(50 * 850), nrow = 50))
#'
#' # The traditional procedure
#' obj1 <- summary(lm(y ~ x + z))
#' pv1 <- sapply(obj1, function(x) x$coefficient[2, 4])
#' qv1 <- qvalue(pv1)$qvalue
#' pos1 <- qv1 <= 0.05
#'
#' # The proposed procedure
#' obj2 <- tdfdr(y, x, z, alpha = 0.05)
#' pos2 <- obj2$pos
#'
#' sum(pos1 & !truth) / max(sum(pos1), 1)
#' sum(pos2 & !truth) / max(sum(pos2), 1)
#' sum(pos1 & truth)
#' sum(pos2 & truth)
#' 
#' @rdname tdfdr
#' @export

tdfdr <- function (y,  x,  z,  
		est.pi0 = TRUE, lambda = 0.5, alpha = 0.05, etype = c('FDR', 'FWER'), 
		t1a = NULL, t2a = NULL, ngrid = 50, npeb.grid = 500,
		parallel = FALSE, cores = NULL, verbose = TRUE) {	

	etype <- match.arg(etype)

	
	if (is.factor(x)) {
		if (nlevels(x) > 2) {
			stop('x should be a vector of two levels!\n')
		} else {
			x <- as.numeric(x)
		}
	}
	
	if (is.vector(z)) {
		z <- as.matrix(z)
	}
	
	if (is.factor(z)) {
		z <- model.matrix(~ z)[, -1, drop = FALSE]
	}
	
	if (is.data.frame(z)) {
		z <- model.matrix(~., data = z)[, -1, drop = FALSE]
	}
	
	if (sum(z[, 1] != 1) == 0) {
		z <- z[, -1, drop = FALSE]
	}
	
	# scaling to improve the stability
	y <- scale(y)
	z <- scale(z)
	x <- as.vector(scale(x))
	
	qrZ <- qr(z, tol = 1e-07)
	z <- qr.Q(qrZ)
	z <- z[, 1:qrZ$rank, drop = FALSE]
	
	n <- length(x)
	p <- ncol(y)
	
	Pz <- diag(n) - z %*% solve(t(z) %*% z) %*% t(z) 
	Px <- diag(n) - x %*% t(x) / sum(x^2) 
	Oz <- as.numeric(t(x) %*% Pz %*% x / n)
	
	O <- sum(x^2) / n
	G <- t(x) %*% z / n
	
	rho <- sqrt(Oz / O)
	
	cat(date(), '\n')
	# We still include intercept for better performance for extremely small sample size
    if (verbose) cat('Perform linear regression ...\n')
	
	M1 <- cbind(1, x)
	M2 <- cbind(1, x, z)
	
	y <- t(y)
	lm.obj1 <- fastLM(y, M1)
	lm.obj2 <- fastLM(y, M2)
	
	coef.x.u <- lm.obj1$est[2, ]
	coef.x.a <- lm.obj2$est[2, ]
	coef.z.a <- lm.obj2$est[-(1:2), , drop = FALSE]
	sigma.est <- lm.obj2$sigma
	
	# need to adjust for multivariate confounder, and collinearity between Z and X
    # if (verbose) cat('Perform EB estimation of sigma ...\n')
	out <- squeezeVar(sigma.est^2, lm.obj2$dof)
	sigma.est <- sqrt(out$var.post)
	
	Zu <- (sqrt(n * O) * coef.x.u / sigma.est)
	Za <- (sqrt(n * Oz) * coef.x.a / sigma.est)
	
	A <- as.numeric(sqrt(1 / O) * G %*% solve(t(z) %*% Px %*% z / n) %*% t(G) * sqrt(1 / O))
	eta <- as.numeric((1 / sqrt(A)) * sqrt(n / O) * G %*% coef.z.a / sigma.est)
	
	if (verbose) cat('Perform NPEB estimation of eta...\n')
	
	mm <- REBayes::GLmix(x = eta)
	normalized.prob <- mm$y / sum(mm$y)
	
#	mm <- cnm(npnorm(eta), grid = npeb.grid)   
#	normalized.prob <- mm$mix$pr
#	mmx <- mm$mix$pt
	
	# Here pi0 is a vector, and could be a global estimate pi0 <- mean(pi0)
	if (est.pi0 == TRUE) {
		pi0 <- (abs(Za) <= lambda) / (1 - 2 * pnorm(-lambda))
		
	} else {
		pi0 <- rep(1, p)
	}
	
	# Search algorithmn to be optimized in the future
	fdr.est <- function (t1, t2, etype) {
		
		# Number of positives
		NP <- sum(abs(Zu) >= t1 & abs(Za) >= t2)
		
		if (NP != 0) {
#			x1 <- -t1 - mmx * sqrt(A)
#			x2 <- t1 - mmx * sqrt(A)
			x1 <- -t1 - mm$x * sqrt(A)
    		x2 <- t1 - mm$x * sqrt(A)
			y1 <- -t2
			y2 <- t2
			
			A1 <- pbivnorm(x = x1, y = y1, rho = rho)
			A2 <- pbivnorm(x = x2, y = y1, rho = rho)
			A3 <- pbivnorm(x = x2, y = y2, rho = rho)
			A4 <- pbivnorm(x = x1, y = y2, rho = rho)
			
			B1 <- pnorm(x1)
			B2 <- pnorm(x2)
			C1 <- pnorm(y1)
			C2 <- pnorm(y2)
			
			if (etype == 'FDR') {
				FD <- sum(pi0) * sum(normalized.prob * (1 + A1 + B1 + C1 + A3 - A2 - B2 - C2 - A4))
				FDP <- FD / NP
			}
			
			if (etype == 'FWER') {
				FDP <- 1 - exp(sum(log(-(A1 + B1 + C1 + A3 - A2 - B2 - C2 - A4)) + log(normalized.prob * length(normalized.prob))
						))
			}
		} else {
			FDP <- 0
		}

		return(list(FDP = FDP, NP = NP))	
	}
	
	if (verbose) cat('Search the best combination of thresholds ...\n')
	
	if (is.null(t1a)) {
		t1a <- seq(min(abs(Zu)), max(abs(Zu)), len = ngrid)
	}
	
	if (is.null(t2a)) {
		t2a <- seq(min(abs(Za)), max(abs(Za)), len = ngrid)
	}
	
	
	if(isTRUE(parallel)){
		
		registerDoMC(cores = cores)
		rr <- foreach(t1 = rep(t1a, each = length(t2a)), t2 = rep(t2a, length(t1a)), .combine = cbind) %dopar% fdr.est(t1, t2, etype)   
		FDP <- unname(unlist(rr[1, ]))
		NP <- unname(unlist(rr[2, ]))
		t1t2 <- paste(rep(t1a, each = length(t2a)), t2a)
		
	} else{
		
		FDP <- NULL
		NP <- NULL
		t1t2 <- NULL
		
		for (t1 in t1a) {
			for (t2 in t2a) {
				obj <- fdr.est(t1, t2, etype)
				NP <- c(NP, obj$NP)
				FDP <- c(FDP, obj$FDP)
				t1t2 <- c(t1t2, paste(t1, t2))
			}
		}
		
	}
	
	
	ind <- which(FDP <= alpha)
	
	if (length(ind) == 0) {
		pos <- rep(FALSE, p)
	} else {
		NP2 <- NP[ind]
		FDP2 <- FDP[ind]
		names(NP2) <- names(FDP2) <- paste(ind)

#		ind2 <- as.numeric(names(which.min(FDP2[NP2 == max(NP2)])))
		
		ind2 <- as.numeric(names(which.max(NP2)))
		temp <- unlist(strsplit(t1t2[ind2], ' '))
		t1 <- as.numeric(temp[1])
		t2 <- as.numeric(temp[2])
		
		pos <- abs(Zu) >= t1 & abs(Za) >= t2
	}
	

	if (verbose) cat('Completed!\n')
	
	return(list(pos = pos, t1 = t1, t2 = t2, t1t2 = t1t2,
					Zu = Zu, Za = Za, NP = NP, FDP = FDP,
					p.value = lm.obj2$pval[, 2], mm = mm, pi0 = mean(pi0)))
}

#' Select the features at a target FDR level
#'
#' The function outputs the features at a target FDR level.
#' 
#' @param tdfdr.obj a list returned from running 'tdfdr'.
#' @param fdr.level a numeric for the target FDR level.
#'
#' @return A list with the elements
#' \item{pos}{a vector of logical values indicating the significant findings (positives) at the target FDR level.}
#' \item{t1, t2}{a numerical value indicating the t1, t2 cutoff values, which achieve the target FDR level.}

#' 
#' 
#' @author ***
#' @references Two-dimensional false discovery rate control for powerful confounder adjustment in omics association studies. Submitted.
#' @keywords FDR, multiple testing, confounder

#' @examples
#' x <- rnorm(50)
#' z <- x + rnorm(50)
#' z <- scale(z)
#' y1 <- x %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y2 <- z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y3 <- x %*% t(rep(0.5, 50)) + z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y <- cbind(y1, y2, y3, matrix(rnorm(50 * 850), nrow = 50))
#' obj2 <- tdfdr(y, x, z, alpha = 0.05)
#' tdfdr.select(obj2, fdr.level = 0.1)
#' @rdname tdfdr.select
#' @export

tdfdr.select <- function (tdfdr.obj, fdr.level = 0.05) {
	FDP <- tdfdr.obj$FDP
	NP <- tdfdr.obj$NP
	Zu <- tdfdr.obj$Zu
	Za <- tdfdr.obj$Za
	t1t2 <- tdfdr.obj$t1t2
	ind <- which(FDP <= fdr.level)
	
	if (length(ind) == 0) {
		pos <- rep(FALSE, length(FDP))
		t1 <- NA
		t2 <- NA
	} else {
		NP2 <- NP[ind]
		names(NP2) <- paste(ind)
		ind <- as.numeric(names(which.max(NP2)))
		temp <- unlist(strsplit(t1t2[ind], ' '))
		t1 <- as.numeric(temp[1])
		t2 <- as.numeric(temp[2])
		
		pos <- abs(Zu) >= t1 & abs(Za) >= t2
	}
	return(list(pos = pos, t1 = t1, t2 = t2))
	
}

#' Plot the tdfdr results
#'
#' The function plots the number of hits (discoveries) against different FDR levels and plots the unadjusted z-score against adjusted z-score with highlighted discoveries. 
#' 
#' @param tdfdr.obj a list returned from running 'tdfdr'.
#' @param fdr.level a numeric vector of the target FDR levels to be plotted.
#' @param nlimit an integer, the number of insignificant z-scores to be sampled to reduce the visualization complexity.
#' @param file.name a character string for the name of the generated pdf file.
#'
#' @return A list of \code{'ggplot2'} objects.
#' 
#' @author ***
#' @references Two-dimensional false discovery rate control for powerful confounder adjustment in omics association studies. Submitted.
#' @keywords FDR, multiple testing, confounder
#' @import ggplot2
#' @import reshape2
#' @importFrom qvalue qvalue
#' @examples
#' x <- rnorm(50)
#' z <- x + rnorm(50)
#' z <- scale(z)
#' y1 <- x %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y2 <- z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y3 <- x %*% t(rep(0.5, 50)) + z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
#' y <- cbind(y1, y2, y3, matrix(rnorm(50 * 850), nrow = 50))
#' obj2 <- tdfdr(y, x, z, alpha = 0.05)
#' tdfdr.plot(obj2)
#' @rdname tdfdr.plot 
#' @export

tdfdr.plot <- function (tdfdr.obj, fdr.level = seq(0.01, 0.2, len = 20), nlimit = 10000, file.name = 'tdfdr.plot') {
	
	qval <- qvalue(tdfdr.obj$p.value)$qvalue
	
	data <- matrix(NA, 2, length(fdr.level))
	dimnames(data) <- list(Method = c('Traditional', 'TDFDR'), Cutoff = paste(fdr.level))
	
	fdr <- rep(NA, length(qval))
	
	for (alpha in rev(fdr.level)) {
		pos <- tdfdr.select(tdfdr.obj, alpha)$pos
		data['TDFDR', paste(alpha)] <- sum(pos)
		fdr[pos] <- alpha
		data['Traditional', paste(alpha)] <- sum(qval <= alpha)
	}
	
	data <- melt(data)
	data$Method <- factor(data$Method, levels = c('TDFDR', 'Traditional'))
	plot.list <- NULL
	pdf(paste0(file.name, '.hits.pdf'), height = 5, width = 6)
	obj <- ggplot(data, aes(x = Cutoff, y = value, col = Method, shape = Method)) + 
			ylab('Number of hits') +
			xlab('FDR cutoff') +
			geom_line() + 
			geom_point(size = 2) + 
			theme_bw()
	print(obj)
	dev.off()
	plot.list[['Hits']] <- obj
	
	data <- data.frame(Zu = abs(tdfdr.obj$Zu), Za = abs(tdfdr.obj$Za),  fdr.cutoff = fdr)
	
	if (nrow(data) > nlimit) {
		data <- data[union(sample(1 : nrow(data), nlimit), which(data$fdr <= max(fdr.level))), ]
	}
	
	cex.pt <- min(1 * nlimit / nrow(data), 2)
			
	pdf(paste0(file.name, '.Zscore.pdf'), height = 5, width = 6)
	obj <- ggplot(data, aes(x = Zu, y = Za, col = fdr.cutoff)) + 
			geom_point(alpha = 0.75, size = cex.pt) +
			ylab('Adjusted absolute Z-score') +
			xlab('Unadjusted absolute Z-score') +
			scale_color_gradient(low="darkred", high="orange", na.value = "grey50") +
			theme_bw()
	print(obj)
	dev.off()
	plot.list[['Zscore1']] <- obj

	obj <- tdfdr.select(tdfdr.obj, 0.05)
	
	if(sum(obj$pos) != 0) {
		
		fdr <- rep('None', length(qval))
		pos1 <- obj$pos
		
		fdr[pos1] <- 'TDFDR'
	
		fdr <- factor(fdr, levels = c('None', 'TDFDR'))
		
		data <- data.frame(Zu = abs(tdfdr.obj$Zu), Za = abs(tdfdr.obj$Za),  Significant = fdr)
		if (nrow(data) > nlimit) {
			data <- data[union(sample(1 : nrow(data), nlimit), which(data$Significant != 'None')), ]
		}
		
		pdf(paste0(file.name, '.Zscore(0.05FDR).pdf'), height = 5, width = 5)
		obj <- ggplot(data, aes(x = Zu, y = Za, col = Significant)) + 
				geom_point(alpha = 0.75, size = cex.pt) +
				ylab('Adjusted absolute Z-score') +
				xlab('Unadjusted absolute Z-score') +
				geom_hline(yintercept = obj$t2, color = 'red') +
				geom_vline(xintercept = obj$t1, color = 'red')
#		if (sum(pos2) >= 1) {
#			obj <- obj + geom_hline(yintercept = min(abs(tdfdr.obj$Za)[pos2]), color = 'blue')
#		}
		obj <- obj +
				scale_color_manual(values = c("grey50", 'darkred')) +
				theme_bw() +
				theme(legend.position = "none")
		print(obj)
		dev.off()
		plot.list[['Zscore2']] <- obj
	}
	
	return(invisible(plot.list))
}

