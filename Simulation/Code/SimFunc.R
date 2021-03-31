#' The main simulation function
#'
#' The simulation function allows varying density and strength of the true and confounding signals, varying confounding strength, 
#' varying degree of co-location between the true and confounding signals, and two common correlation structures (block and AR(1)).
#'  
#' @param n a numeric value, sample size.
#' @param p a numeric value, feature size.
#' @param conf.sig.cor a numeric value controlling the correlation between the variable of interest and the confounder. 
#' @param sig.density a numeric value (0-1), the density of true signals.
#' @param sig.strength.m a numeric value, minimum strength (regression coefficient) for the true signals.
#' @param sig.strength.sd a numeric value, the maximum increment over the minimum strength for the true signals.
#' @param conf.density a numeric value (0-1), the density of confounding signals.
#' @param conf.strength.m a numeric value, minimum strength (regression coefficient) for the confounding signals.
#' @param conf.strength.sd a numeric value, the maximum increment over the minimum strength for the confounding signals.
#' @param conf.sig.loc a character string indicating the type of co-location between the true and confounding signals.
#' @param coloc.prob a numeric vector, the proability of co-location.
#' @param cor.struct a character string indicating the type of correlation structure. 
#' @param rho a numeric value, the correlation between features
#' @param nblock an integer number, the number of blocks 
#'
#' @return A list with the elements
#' \item{x}{a vector of numeric values, the variable of interest.}
#' \item{z}{a vector of numeric values, the confounding variable.}
#' \item{y}{a numeric matrix (n by p) of the measurements of genomic features.}
#' \item{coeff.x}{a numeric vector of the coefficients for the variable of interest.}
#' \item{coeff.z}{a numeric vector of the coefficients for the confounding variable.}
#' \item{truth}{a numeric vector of 0 and 1 indicating the true signals. }
#' 
#' @author ***
#' @references Two-dimensional false discovery rate control for powerful confounder adjustment in omics association studies.
#' @keywords Simulation
#' @importFrom stats arma.sim
#' 
#' @rdname simulate.data
#' @export

simulate.data <- function (n = 100, p = 10000, conf.sig.cor = 1.25, dimZ = 1,
		sig.density = 0.1, sig.strength.m = 0.4, sig.strength.sd = 0.2,
		conf.density = 0.1, conf.strength.m = 0.4, conf.strength.sd = 0.2,
		conf.sig.loc = c('Random', 'NonCoLoc', 'CoLoc'), coloc.prob = 0.5, 
		cor.struct = c('Indep', 'Block1', 'AR1'), rho = 0, nblock = 100) {
	
	conf.sig.loc <- match.arg(conf.sig.loc)
	cor.struct <- match.arg(cor.struct)
	# Generate correlated x and z
	x0 <- rnorm(n)
	x <- scale(conf.sig.cor * x0 + rnorm(n))
	z <- scale(conf.sig.cor * x0 + matrix(rnorm(n * dimZ), n, dimZ))
	
	p.sig <- round(p * sig.density)
	p.conf <- round(p * conf.density)
	
	if (cor.struct %in% c('Block1')) {
		# Signal conforms to the block structure
		coef.x <- c(runif(round(p.sig / 2), sig.strength.m, sig.strength.m + sig.strength.sd), 
						runif(p.sig - round(p.sig / 2),  -(sig.strength.m + sig.strength.sd), -sig.strength.m), rep(0, p - p.sig))
	} else {
		coef.x <- sample(c(runif(round(p.sig / 2), sig.strength.m, sig.strength.m + sig.strength.sd), 
						runif(p.sig - round(p.sig / 2),  -(sig.strength.m + sig.strength.sd), -sig.strength.m), rep(0, p - p.sig)))
	}

	coef.z <- sample(c(runif(round(p.conf / 2), conf.strength.m, conf.strength.m + conf.strength.sd), 
					runif(p.conf - round(p.conf / 2), -(conf.strength.m + conf.strength.sd), -conf.strength.m), rep(0, p - p.conf)))
	
	if (conf.sig.loc == 'NonCoLoc') {
		coef.x0 <- rep(0, p)
		if (sum(coef.z == 0) > p.sig) {
			coef.x0[sample(which(coef.z == 0), p.sig)] <- coef.x[coef.x != 0]
		} else {
			stop("Couldn't achieve non-colocation!\n")
		}
		coef.x <- coef.x0
	}
	
	if (conf.sig.loc == 'CoLoc') {
		if (sig.density < conf.density) {
			x.n1 <- round(p.sig * coloc.prob)
			x.n3 <- p.sig - x.n1
			x.n2 <- p - x.n1 - x.n3
			x.vec <- c(rep(TRUE, x.n1), rep(FALSE, x.n2), rep(TRUE, x.n3))
			z.vec <- c(rep(TRUE, p.conf), rep(FALSE, p - p.conf))
			
		} else {
			z.n1 <- round(p.conf * coloc.prob)
			z.n3 <- p.conf - z.n1
			z.n2 <- p - z.n1 - z.n3
			z.vec <- c(rep(TRUE, z.n1), rep(FALSE, z.n2), rep(TRUE, z.n3))
			x.vec <- c(rep(TRUE, p.sig), rep(FALSE, p - p.sig))
		}
		coef.x0 <- coef.z0 <- rep(0, p)
		coef.x0[x.vec] <- coef.x[coef.x != 0]
		coef.z0[z.vec] <- coef.z[coef.z != 0]
		ind <- sample(p)
		coef.x <- coef.x0[ind]
		coef.z <- coef.z0[ind]
	}
	
	
	if (cor.struct == 'Indep') {
		epsilon <- matrix(rnorm(n * p), n, p)
	}
	
	if (cor.struct == 'Block1') {
		
		Tmat <- function (feature.no = 10000, nb = 100, rho = 0) {
			# Simulate special correlation structure
			bs <- feature.no / nb
			mat <- diag(bs)
			mat[, ] <- rho
			diag(mat) <- 1
			obj <- eigen(mat)
			T1 <- obj$vectors %*% diag(sqrt(obj$values)) %*% t(obj$vectors)
			
			mat <- diag(bs)
			mat[, ] <- -rho
			mat[1:(bs/2), 1:(bs/2)] <- rho
			mat[(bs/2+1):bs, (bs/2+1):bs] <- rho
			diag(mat) <- 1
			obj <- eigen(mat)
			T2 <- obj$vectors %*% diag(sqrt(obj$values)) %*% t(obj$vectors)
			
			return(list(T1=T1, T2=T2))
		}
		
		obj <- Tmat(p, nblock, rho)
		bs <- p / nblock
		epsilon <- t(sapply(1:n, function (i) as.vector(obj$T1 %*% matrix(rnorm(p), bs, nblock))))
	} 

	if (cor.struct == 'AR1') {
		epsilon <- t(matrix(arima.sim(n = p * n * 2, list(ar = c(rho)), n.start = 100000) * sqrt(1 - rho^2), p, 2 * n))
		epsilon <- epsilon[2 * (1 : n), ]

	}
	
	y <- x %*% t(coef.x) + z %*% matrix(coef.z, nrow = dimZ, ncol = p, byrow = TRUE) + epsilon
	
	truth <- as.numeric(coef.x != 0)
	
	return(list(x = x, z = z, y = y, coef.x = coef.x, coef.z = coef.z, truth = truth))
	
}

