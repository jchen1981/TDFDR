# tdfdr
Two-dimensional false discovery rate control for powerful confounder adjustment in omics association studies

## Overview
The package implements the two-dimensional false discovery rate control for powerful confounder adjustment in omics association analysis. The method is based on the idea that the confounder(s) usually affect part of the omics features, and thus adjusting the confounder(s) for ALL omics features will be over-adjustment, leading to reduced statistical power.  The proposed procedure starts with performing the unadjusted analysis (first dimension - filtering) to narrow down the list of omics features that are more likely to be affected by either the confounder or the variable of interest or both. In the second dimension, we conduct confounder-adjusted analysis on these 'top' candidates, which are enriched in signals, to reduce multiple testing burden and increase the power. The method belongs to the general topic of using auxiliary data to increase the power of multiple testing, which has recently received tremendous research interest. In our case, the auxiliary data are the the unadjusted statistics, which could inform the probability of the null hypotheses being true.  The difficulty here is to take into account the correlation between the auxiliary data (unadjusted statistics) and the main data (adjusted statistics). We provide a procedure that is theoretically guaranteed to control the false discovery rate while maximizing the power.

## License
GPL version 2 or any later version

## Installation 

### Install dependent packages "ggplot2", "reshape2", "doMC", "pbivnorm", "nspmix", "limma", "qvalue"

```
install.packages(c("ggplot2", "reshape2", "doMC", "pbivnorm", "nspmix"))
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
 source('https://raw.githubusercontent.com/jchen1981/TDFDR/master/Simulation/Code/SimFunc.R')

 
 # Simulate data
 data <- simulate.data(n = 100, p = 1000, conf.sig.cor = 1.25, dimZ = 2,
	    sig.density = 0.1, sig.strength.m = 0.4, sig.strength.sd = 0.2,
	    conf.density = 0.1, conf.strength.m = 0.4, conf.strength.sd = 0.2)
 x <- data$x
 z <- data$z
 y <- data$y
 truth <- data$truth
		
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
All the codes to reproduce the results in the manuscript are contained in 'Simulation' and 'RealData' folders. Due to file limit on github, we provide a dropbox link to download the relevant data.


