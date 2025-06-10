# mremaR
An R package for implementing gene set analysis tests using a mixture of random-effect meta-analysis approach. 

## Installation 
 To install this package in R, use 
 
 ```
    library("devtools");
    install_github("osedo/mremaR")
 ```
 
 ## References
 
 Markrooni et al., Random-effects meta-analysis of effect sizes as a unified framework for gene set analysis. PLOS Comp Biol (2022). https://doi.org/10.1371/journal.pcbi.1010278
 
 ## Example Pipeline
 
 ```
library(mremaR)
```

Here we will show a simple workflow for carrying out cell type-specific gene set analysis (CT-GSA). First let's load the CAR-seq package which is designed for carrying out cell type-specific differential expression analysis. 

```
if(!require(CARseq)){
    devtools::install_github("chongjin/CARseq")
    library(CARseq)
}
```

Now we will simulate gene expression for samples which are a mixture of six different cell types. This defaults to 1,000 genes in 100 cases and 100 controls. In each cell type 5% of genes are assigned a log-fold-change drawn from a normal distribution with a mean of plus/minus 4 and a standard deviation of 0.2. 

```
sim <- CTexpSimulation()
sim$bulk.expression[c(1:5),c(1:5)] # bulk expression counts
sim$ols.mixture.estimate[c(1:5),] # cell type proportion estimates 
table(sim$trait) # sample traits
```

These outputs can then be passed to the CAR-seq package. The important results for us are the shrunken_lfc and shrunken_lfcSE matrices

```
# this will take a minute or two to run
res <- run_CARseq(
  count_matrix = sim$bulk.expression,
  cellular_proportions = sim$ols.mixture.estimate,
  groups = sim$trait,
  cores = 1
)
```

Let's simulate a few gene sets
```
# a few power gene sets, gs1 enriched for DE genes in cellType1, same for gs2 and so on
# we will give sets of 100 approx 20% DE gene
gs.power <- lapply(c(1:6), function(x){
  paste0("gene.", 
         # sample DE genes
         c(sample(which(sim$simulated.lfc[,x] != 0), 20),
         # sample non DE genes  
           sample(which(sim$simulated.lfc[,x] == 0), 80)))
})

# random gene sets
gs.null <- lapply(c(7:20), function(x){
  paste0("gene.", sample.int(1000, 100))
})

gs <- c(gs.power, gs.null)
names(gs) <- paste0("gs", 1:20)
```


Let's run mrema() on all cell-types

```
mrema.res <- lapply(c(1:6), function(ct){
  postdata <- data.frame("Ensembl" = rownames(res$shrunken_lfcSE), "effect" = res$shrunken_lfc[,ct], "variance" = res$shrunken_lfcSE[,ct]^2, "pval" = res$p[,ct])
  mrema_res <- mrema(postdata = postdata, raw.gs = gs, threshold = 1.5, DF = 4, ncores = 1)
  mrema_res
  })
names(mrema.res) <- colnames(res$p)
saveRDS(mrema.res, file = "sample.RDS")
```

We can now launch the shiny app to investigate our results

```
mremApp()
```

Test

