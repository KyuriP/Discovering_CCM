library(dplyr)

mats <- mat_subsample_dep

# frequency table per cell
comb_mat <- simplify2array(mats) %>% apply(c(1, 2), table)

# get the highest freq value .... what do i want with this? 
# Lourens: what occurs more than 80% --> keep, otherwise discard: place circle instead?
freqlist <- comb_mat %>% 
  purrr::map_dfr(
    ~as.data.frame(.x) %>% 
      summarise(val = Var1[which.max(Freq)], 
                freq = max(Freq))
  )


# reconstruct matrix
subsampling_adjmat <- freqlist %>% 
  rowwise() %>% 
  summarize(val = case_when(
    freq > 800 ~ val,
    .default = "0"
  )) %>% 
  unlist %>% 
  as.numeric %>% 
  matrix(nrow=16, ncol=16) %>% 
  `dimnames<-`(dimnames(comb_mat))


pcalg::plotAG(subsampling_adjmat)
