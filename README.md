# tdfdr
Two-dimensional false discovery rate control for powerful confounder adjustment in omics association studies

## Overview
The function implements the two-dimensional false discovery rate control for powerful confounder adjustment in omics association analysis. The method is based on the idea that the confounder(s) usually affect part of the omics features, and thus adjusting the confounder(s) for ALL omics features will be over-adjustment, leading to reduced statistical power.  The proposed procedure starts with performing the unadjusted analysis (first dimension - filtering) to narrow down the list of omics features which are more likely to be affected by either the confounder or the variable of interest or both. In the second dimension, we conduct confounder-adjusted analysis on these 'top' candidates, which are enriched in signals, to reduce multiple testing burden and increase the power. The method belongs to the general topic of using auxiliary data to increase the power of multiple testing, which has recently received tremendous research interest. In our case, the auxiliary data are the the unadjusted statistics, which could inform the probability of the null hypotheses being true.  The difficulty here is to take into account the correlation between the auxiliary data (unadjusted statistics) and the main data (adjusted statistics). We provide a procedure that is theoretically guaranteed to control the false discovery rate while maximizing the power.

## Installation 
### Please install mosek first 

** Following instruction is for macOS. For linux, change 'osx' to 'linux'.**


```
The following instruction was modified from
https://gist.github.com/mikelove/67ea44d5be5a053e599257fe357483dc

1) Download mosek from here (** please download version 8 **):
https://www.mosek.com/downloads/list/8/
(I downloaded this to my ~/bin)

cd ~/bin
tar -xvf mosektoolsosx64x86.tar.bz2

2) Add this to your ~/.bashrc
export PATH=$PATH:~/bin/mosek/8/tools/platform/osx64x86/bin

3) Get academic license:
https://www.mosek.com/products/academic-licenses/
Check email, put licsense file in ~/mosek
(You need to create ~/mosek directory)

4) Install:

export PKG_MOSEKHOME=~/bin/mosek/8/tools/platform/osx64x86
export PKG_MOSEKLIB=mosek64

5) For macOS Catalina, run "xattr -dr com.apple.quarantine mosek" to 
prevent security exceptions.

Then in R:
install.packages("Rmosek", type="source", INSTALL_opts="--no-multiarch", 
repos="http://download.mosek.com/R/8")

We also provide an alternative function 'install.mosek.R' to install mosek if the procedure above fails.
```
source('installmosek.R') 
installmosek()

```
### Install dependent packages 

"ggplot2", "reshape2", "doMC", "pbivnorm", "REBayes", "limma", "qvalue"

```
install.packages(c("ggplot2", "reshape2", "doMC", "pbivnorm", "REBayes"))
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
BiocManager::install(c("limma", "qvalue"))
install.packages("devtools")
devtools::install_github("jchen1981/TDFDR/R_package")
```



### An Example
We illustrate the usage of tdfdr package using simulated data.

```
     require(tdfdr)
     require(qvalue)
     
     # Generate simulated data with 100 true positives out of 1000
     truth <- c(rep(1, 50), rep(0, 50), rep(1, 50), rep(0, 850))
     x <- rnorm(50)
     z <- x + rnorm(50)
     z <- scale(z)

     y1 <- x %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
     y2 <- z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
     y3 <- x %*% t(rep(0.5, 50)) + z %*% t(rep(0.5, 50)) + matrix(rnorm(50 * 50), 50, 50)
     y <- cbind(y1, y2, y3, matrix(rnorm(50 * 850), nrow = 50))

     # One-dimenisonal procedure - classic adjusted procedure
     obj1 <- summary(lm(y ~ x + z))
     pv1 <- sapply(obj1, function(x) x$coefficient[2, 4])
     qv1 <- qvalue(pv1)$qvalue
     pos1 <- qv1 <= 0.05

     # Two-dimensional procedure
     obj2 <- tdfdr(y, x, z, alpha = 0.05)
     pos2 <- obj2$pos

     # Compare the false discovery proportions
     sum(pos1 & !truth) / max(sum(pos1), 1)
     sum(pos2 & !truth) / max(sum(pos2), 1)
     
     # Compare the number of hits
     sum(pos1 & truth)
     sum(pos2 & truth)
  
```

### Reproducibility
All the codes to reproduce the results in the manuscript are contained in 'Simulation' and 'RealData' folders.

