library(dplyr)
source("code/R/plot_fnc.R")


## run CCD on subsamples
ccd_subsample_dep <- subsamples %>% 
  map(~ccdKP(df=.x, dataType = "continuous", alpha=alpha))
# create an adjacency matrix for PAG
mat_subsample_dep <- ccd_subsample_dep %>% 
  map(~CreateAdjMat(.x, length(.x$nodes)))

# save(mat_subsample_dep, file="data/empirical/mat_subsample_dep.RData")

## run FCI on subsamples
fci_subsample_dep <- subsamples %>%
  map(~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
           alpha = alpha, doPdsep = TRUE, selectionBias= FALSE,
           labels = colnames(.x)) %>% 
        .@amat
  )
# save(fci_subsample_dep, file="data/empirical/fci_subsample_dep.RData")

## run CCI on subsamples
cci_subsample_dep <- subsamples %>%
  map(~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha,
           labels = colnames(.x), p = ncol(.x)) %>% 
        .$maag
  )
# save(cci_subsample_dep, file="data/empirical/cci_subsample_dep.RData")

load("data/empirical/mat_subsample_dep.RData")
load("data/empirical/fci_subsample_dep.RData")
load("data/empirical/cci_subsample_dep.RData")

# frequency table per cell
combmat_CCD <- simplify2array(mat_subsample_dep) %>% apply(c(1, 2), table)
combmat_FCI <- simplify2array(fci_subsample_dep) %>% apply(c(1, 2), table)
combmat_CCI <- simplify2array(cci_subsample_dep) %>% apply(c(1, 2), table)


# get the highest freq value .... what do i want with this? 
# Lourens: what occurs more than 80% --> keep, otherwise discard: place circle instead?
freqlistCCD <- combmat_CCD %>% 
  purrr::map_dfr(
    ~as.data.frame(.x) %>% 
      summarise(val = Var1[which.max(Freq)], 
                freq = max(Freq))
  )

freqlistFCI <- combmat_FCI %>% 
  purrr::map_dfr(
    ~as.data.frame(.x) %>% 
      summarise(val = Var1[which.max(Freq)], 
                freq = max(Freq))
  )

freqlistCCI <- combmat_CCI %>% 
  purrr::map_dfr(
    ~as.data.frame(.x) %>% 
      summarise(val = Var1[which.max(Freq)], 
                freq = max(Freq))
  )


# reconstruct matrix
adjmat_CCD <- freqlistCCD %>% 
  rowwise() %>% 
  summarize(val = case_when(
    freq > 800 ~ val,
    .default = "0"
  )) %>% 
  unlist %>% 
  as.numeric %>% 
  matrix(nrow=16, ncol=16) %>% 
  `dimnames<-`(dimnames(combmat_CCD))

adjmat_FCI <- freqlistFCI %>% 
  rowwise() %>% 
  summarize(val = case_when(
    freq > 800 ~ val,
    .default = "0"
  )) %>% 
  unlist %>% 
  as.numeric %>% 
  matrix(nrow=16, ncol=16) %>% 
  `dimnames<-`(dimnames(combmat_FCI))

adjmat_CCI <- freqlistCCI %>% 
  rowwise() %>% 
  summarize(val = case_when(
    freq > 800 ~ val,
    .default = "0"
  )) %>% 
  unlist %>% 
  as.numeric %>% 
  matrix(nrow=16, ncol=16) %>% 
  `dimnames<-`(dimnames(combmat_CCI))

# plot the models
plotAG(adjmat_CCD)
plotAG(adjmat_FCI)
plotAG(adjmat_CCI)
