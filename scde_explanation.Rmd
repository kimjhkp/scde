---
title: "R Notebook"
output: html_notebook
---

A personal note on SCDE package developed by the Kharchenko lab  
http://hms-dbmi.github.io/scde/diffexp.html

Let's load the data and associated functions

```{r}
library(flexmix)

load(file.path(proj.dir, "scde-master/data/es.mef.small.rda"))
source(file.path(proj.dir, "scde-master/R/functions.R"))

```
We can group the cells based on prior information 

``` {r}
# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
# the group factor should be named accordingly
names(sg) <- colnames(es.mef.small)  
table(sg)

```

Clean the data

``` {r}

cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)

```

**scde.error.models**

This is the main model for fitting the errors. This can be divided into two parts:  
1) Crossfit error model  
2) Individual error model  

``` {r}

scde.error.models <- function(counts, groups = NULL, min.nonfailed = 3, threshold.segmentation = TRUE, min.count.threshold = 4, zero.count.threshold = min.count.threshold, zero.lambda = 0.1, save.crossfit.plots = FALSE, save.model.plots = TRUE, n.cores = 12, min.size.entries = 2e3, max.pairs = 5000, min.pairs.per.cell = 10, verbose = 0, linear.fit = TRUE, local.theta.fit = linear.fit, theta.fit.range = c(1e-2, 1e2)) {
    # default same group
    if(is.null(groups)) {
        groups <- as.factor(rep("cell", ncol(counts)))
    }
    # check for integer counts
    if(any(!unlist(lapply(counts,is.integer)))) {
      stop("Some of the supplied counts are not integer values (or stored as non-integer types). Aborting!\nThe method is designed to work on read counts - do not pass normalized read counts (e.g. FPKM values). If matrix contains read counts, but they are stored as numeric values, use counts<-apply(counts,2,function(x) {storage.mode(x) <- 'integer'; x}) to recast.");
    }

    # crossfit
    if(verbose) {
        cat("cross-fitting cells.\n")
    }
    cfm <- calculate.crossfit.models(counts, groups, n.cores = n.cores, threshold.segmentation = threshold.segmentation, min.count.threshold = min.count.threshold, zero.lambda = zero.lambda, max.pairs = max.pairs, save.plots = save.crossfit.plots, min.pairs.per.cell = min.pairs.per.cell, verbose = verbose)
    # error model for each cell
    if(verbose) {
        cat("building individual error models.\n")
    }
    ifm <- calculate.individual.models(counts, groups, cfm, min.nonfailed = min.nonfailed, zero.count.threshold = zero.count.threshold, n.cores = n.cores, save.plots = save.model.plots, linear.fit = linear.fit, return.compressed.models = TRUE, verbose = verbose, min.size.entries = min.size.entries, local.theta.fit = local.theta.fit, theta.fit.range = theta.fit.range)
    rm(cfm)
    gc()
    return(ifm)
}
```

**1. Calculating the crossfit error**

This function uses a 3 components mixture model to determine correlated and expressed genes across different cell pairs

``` {r}

cfm <- calculate.crossfit.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, min.count.threshold = 4, zero.lambda = 0.1,
                                 max.pairs = 5000, save.plots = FALSE, min.pairs.per.cell = 10)
```

stepFlexmix determines the maximum likelihood of different mixture models - eg. which of the three out of four models defined here give the best ML estimate ?

``` {r}

mo1 <- FLXMRglmCf(c1~1, family = "poisson", components = c(1), mu = log(zero.lambda))
mo2 <- FLXMRnb2glmC(c1~1+I(log(c2+1)), components = c(2))
mo3 <- FLXMRnb2glmC(c2~1+I(log(c1+1)), components = c(2))
mo4 <- FLXMRglmCf(c2~1, family = "poisson", components = c(3), mu = log(zero.lambda))
m1 <- mc.stepFlexmix(c1~1, data = df[vi, ], k = 3, model = list(mo1, mo2, mo3, mo4), control = list(verbose = verbose, minprior = min.prior),
                     concomitant = FLXPmultinom(~I((log(c1+1)+log(c2+1))/2)+1), 
                     cluster = cbind(df$c1[vi]<= min.count.threshold, df$c1[vi] > min.count.threshold & df$c2[vi] > min.count.threshold, df$c2[vi]<= min.count.threshold), 
                     nrep = nrep)

```


**2. Calculate the individual fit**

This function uses the prior information calculated from the crossfit model to:

1) Select genes that are well correlated between different cell pairs
2) Uses the median expression of the genes (expected value) and failure rate (due to dropout) and determines  
parameters for negative binomial distribution (dispersion and mu)

``` {r}

o.ifm <- calculate.individual.models(counts = cd, groups = sg, cfm = cfm, min.nonfailed = 3, zero.count.threshold = 4, n.cores = 1, save.plots = FALSE,
                                   linear.fit = TRUE, return.compressed.models = TRUE, verbose = 0, min.size.entries = 2e3, local.theta.fit = linear.fit, theta.fit.range = c(1e-2,
                                                                                                                                                                              1e2))

```

corr.a and corr.b are slope and intercept of the correlated component fit
corr.theta is the over-dispersion parameter
fail.r is the background Poisson rate

``` {r}

load(file.path(proj.dir,"scde-master/data/o.ifm.rda"))
o.ifm

```
**3. Gene expression prior**    

Calculate the gene expression prior using the parameters defined in the previous step (dispersion and mu)

``` {r}

o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
o.prior

```

**4. Differential expression analysis**

``` {r}

groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels  =  c("ESC", "MEF"))
names(groups) <- row.names(o.ifm)

diff <- scde.expression.difference(o.ifm, cd, o.prior, 
                                   groups = groups, n.randomizations = 100,
                                   n.cores = 1, verbose = 1)

```










