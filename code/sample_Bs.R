B5dense = matrix(c(0, 0, 0, 0, 0,
                   1.4, 0, 0.8, 0, 0,
                   0, 0, 0, 0.9, 0,
                   0, 0.7, 0, 0, 1,
                   1, 0, 0, 0, 0), 5, 5, byrow = T)

qgraph(cor(simdata_5pdense[[11]][[3]]), layout = layout5, theme="colorblind", graph = "pcor")


### conditional indep test 

alpha <- .05
data <- simdata_5pdense[[11]][[3]]

cor.test(data[,"X2"], data[,"X5"])
ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X1")])

count <- matrix(0, 8, 11)
partial <- c()
for(j in 1:500){
  data <- simdata_5pdense[[11]][[j]]
  partial[j] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3", "X4")]) $estimate
}

hist(partial)

count <- matrix(0, 8, 11)
indep <- matrix(NA, 8, 11)
for(j in 1:50){
  for(i in 1:11){
    data <- simdata_5pdense[[i]][[j]]
    #indep <- matrix(NA, 8, 11)
    indep[1,i] <- cor.test(data[,"X2"], data[,"X5"])$p.value
    indep[2,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X1")])$p.value
    indep[3,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X3")])$p.value
    indep[4,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,("X4")])$p.value
    indep[5,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X4")])$p.value
    indep[6,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3")])$p.value
    indep[7,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X3", "X4")])$p.value
    indep[8,i] <- ppcor::pcor.test(data[,"X2"], data[,"X5"], data[,c("X1", "X3", "X4")])$p.value
  }
  count <- ifelse(indep >= 0.05, count+1, count)
  #print(count)
}

rownames(count) <- c("marginal", "condX1", "condX3", "condX4", "condX1&X4", "condX1X3", "condX3X4", "condX1X3X4")
colnames(count) <- c("N=50" , "N=150", "N=500", "N=1000", "N=1500", "N=2000", "N=2500", "N=3000", "N=4000", "N=5000", "N=10000")

count;



#### Sample Bs
B5dense = matrix(c(0, 0, 0, 0, 0,
                   1.4, 0, 0.8, 0, 0,
                   0, 0, 0, 0.9, 0,
                   0, 0.7, 0, 0, 1,
                   1, 0, 0, 0, 0), 5, 5, byrow = T)


randomB <- function(B){
  samplevals <- runif(sum(B!=0), -0.8, 0.8)
  B[B!=0] <- samplevals  
  return(B)
}

# check the pattern with random sampling B

ranB <- randomB(B5dense)
# equilibrium check
equilibrium_check(ranB)

# generate data (sample size as specified above)
simdat <- N %>% future_map(function(z) {
  replicate(n = n,
            expr = gen_dat(ranB, N = z),  
            simplify = FALSE)
}, .options = furrr_options(seed=123))


## Run CCD algorithm
ccdres <- simdat %>%
  map_depth(2, ~ ccdKP(df = .x, dataType = "continuous", alpha = alpha)
  )
ccd_mat <- ccdres %>% 
  map_depth(2, ~CreateAdjMat(.x, length(.x$nodes)))


## Run FCI algorithm
fcires <- simdat %>%
  map_depth(2, ~fci(list(C = cor(.x), n = nrow(.x)), indepTest=gaussCItest,
                    alpha = alpha, doPdsep = TRUE, selectionBias= FALSE, labels = colnames(.x)) %>% .@amat # extract amat
  )


## Run CCI algorithm
ccires <- simdat %>%
  map_depth(2, ~cci(list(C = cor(.x), n = nrow(.x)), gaussCItest, alpha=alpha, labels = colnames(.x), p = ncol(.x)) %>% .$maag  # convert some logical matrix (0, 1 only) to a numeric matrix while keeping a matrix format (lost the row names but they are not needed)
  )


## evaluation
# CCD
# res_ccd <- ccd_mat %>% 
#   map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
#   do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame()

# UNCERTAINTY
uncer_ccd <- ccd_mat %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_ccd, na.rm=T)

# SHD
# SHD_ccd <- ccd_mat %>% 
#   map_depth(2, ~SHD(trueag_5pdense, .x)) %>% 
#   do.call("cbind", .) %>% apply(., 2, unlist) %>%  
#   as.data.frame %>% rename_with(~ paste0("N = ", N))
# 
# colMeans(SHD_ccd5pdense)

# FCI
# res_fci <- fcires %>% 
#   map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
#   do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY
uncer_fci <- fcires%>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_fci, na.rm=T)

# SHD
# SHD_fci5pdense <- fci_5pdense %>% 
#   map_depth(2, ~SHD(trueag_5pdense, .x)) %>% 
#   do.call("cbind", .) %>% apply(., 2, unlist) %>%  
#   as.data.frame %>% rename_with(~ paste0("N = ", N))
# 
# colMeans(SHD_fci5pdense)

# CCI
# res_cci5pdense <- cci_5pdense %>% 
#   map_depth(2, ~precision_recall(trueag_5pdense, .x)) %>% 
#   do.call("cbind", .) %>% t() %>%  apply(., 2, unlist) %>%  as.data.frame() 

# UNCERTAINTY 
uncer_cci <- ccires %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))

colMeans(uncer_cci)

# # SHD
# SHD_cci5pdense <- cci_5pdense %>% 
#   map_depth(2, ~SHD(trueag_5pdense, .x)) %>% 
#   do.call("cbind", .) %>% apply(., 2, unlist) %>%  
#   as.data.frame %>% rename_with(~ paste0("N = ", N))
# 
# colMeans(SHD_cci5pdense)

## CCI 5p dense case
# high frequency
par(mfrow=c(2,6))

ccd_mat %>% 
  map(~high_freq(.x, p = 5) %>% 
        plotAG)

## FCI 5p dense case
# high frequency 
par(mfrow=c(2,6))

fcires %>% 
  map(~high_freq(.x, p = 5) %>% 
        plotAG)

## CCI 5p dense case
# high frequency 
par(mfrow=c(2,6))

ccires %>% 
  map(~high_freq(.x, p = 5) %>% 
        plotAG)

## ======================================================================


generate_DCG_LE <- function(p,en){
  
  DCG=list();
  graph_p=matrix(0,p,p);
  N = p*p - p;
  
  while(!isCyclic(graph_p>0)){
    
    samplesB = rbinom(N,1, en/(2*(p-1)) );
    
    graph_p=matrix(0,p,p);
    
    graph_p[upper.tri(graph_p, diag=FALSE)] <- samplesB[1:(N/2)];
    graph_p[lower.tri(graph_p, diag=FALSE)] <- samplesB[(N/2+1):N];
    
  }
  
  DCG$graph_p = graph_p;
  
  DCG$weights = matrix((0.75*runif(p^2)+0.25)*sample(c(-1,1),p^2,replace=TRUE),p,p)
  
  nS = sample(0:3,1);
  nL = sample(0:3,1);
  # nS = 0
  # nL = 0
  
  pL = which(rowSums(DCG$graph_p)>=2); #variables with >=2 children
  DCG$L = sample(pL, min(length(pL),nL));
  
  pS = which(colSums(DCG$graph_p)>=2); #variables with >=2 parents
  DCG$S = sample(pS, min(length(pS),nS));
  
  graph = matrix(0,p+length(DCG$S), p+length(DCG$S));
  graph[1:p,1:p] = graph_p;
  
  for (s in seq_len(length(DCG$S))){
    graph[DCG$S[s],p+s]=1;
  }
  
  DCG$graph=graph;
  
  S_prob = runif(length(DCG$S),0.1,0.5);
  
  if (length(DCG$S)>0){
    DCG$elim_prop = runif(1,0.1,0.5);
  } else{
    DCG$elim_prop = 0;
  }
  
  DCG$S_prob = DCG$elim_prop*(S_prob / sum(S_prob));
  
  return(DCG)
}



a_DCG <- generate_DCG_LE(10,10)

sample_DCG = sample_DCG_LE(nsamps=1000, a_DCG) #generate Gaussian samples from the DCG

suffStat=list(); suffStat$C = cor(sample_DCG); suffStat$n = 1000; # get all of the parameters needed by Fisher's z test

G=cci(suffStat,gaussCItest,alpha=0.95,p=ncol(sample_DCG))$maag

# a <- fci(suffStat,gaussCItest,alpha=0.05,p=ncol(sample_DCG))
# a

uncertainty(G)

uncer <- c()
for(i in 1:30){
a_DCG <- generate_DCG_LE(10,10)
sample_DCG = sample_DCG_LE(nsamps=1000, a_DCG) #generate Gaussian samples from the DCG
suffStat=list(); suffStat$C = cor(sample_DCG); suffStat$n = 1000; # get all of the parameters needed by Fisher's z test
G=cci(suffStat,gaussCItest,alpha=0.05,p=ncol(sample_DCG))$maag
uncer[i] <- uncertainty(G)
}



B10dense = matrix(c(0, , 1, 0, 0, 0, 0, 1, 1, 0.7,
                    0.3, 0, 0.8, 0, 0, 1, 0, 1, 1, 1, 
                    0.4, 1, 0, 0, 0, 0, 0, 0, 0, 1, 
                    1, 0, 1.1, 0, 0, 0.9, 0, 0, 0, 1, 
                    1, 0.4, 0, 1, 0, 0, 0, 0, 0, 1, 
                    1, 1, 0, 0, 0.9, 0, 0.5, 0, 0, 1, 
                    1, 1, 0, 0, 1, 0, 0, 1, 0.8, 1, 
                    1, 0, 0, 1, 0, 0, 1, 0, 0, 0.4, 
                    1, 0, 1, 1, 1, 0, 0, 0, 0, 0.7, 
                    1, 0, 1, 0, 1, 0, 0, 1, 1, 0), 10, 10, byrow = T)


sum(B10dense!=0)
equilibrium_check(B10dense)
true10pdense <- qgraph(t(B10dense), layout = layout10, labels = colnames(B10dense), theme="colorblind")

simdata_10pdense <- gen_dat(B10dense, N = 100000)
colnames(simdata_10pdense) <-paste("X", 1:10, sep="")

# cci_10pdense  <- simdata_10pdense %>%
#   cci(list(C = cor(.), n = nrow(.)), gaussCItest, alpha=0.01, labels = colnames(.), p = ncol(.)) %>% .$maag

ccdKP(df = simdata_10pdense, dataType = "continuous", alpha = 0.1) %>% CreateAdjMat(., p = 10)
# %>% uncertainty

cci(list(C = cor(simdata_10pdense), n = nrow(simdata_10pdense)), gaussCItest, alpha=0.1, labels = colnames(simdata_10pdense), p = ncol(simdata_10pdense)) %>% .$maag %>% plotAG
# %>% uncertainty()
  
# UNCERTAINTY
uncer_cci10pdense <- cci_10pdense %>% 
  map_depth(2, ~uncertainty(.x)) %>% 
  do.call("cbind", .) %>% apply(., 2, unlist) %>%  
  as.data.frame %>% rename_with(~ paste0("N = ", N))
# average uncertainty
colMeans(uncer_cci10pdense, na.rm=T)




B5sparse = matrix(c(0, 0, 0, 0, 0,
                    1, 0, 0.8, 0, 0,
                    0, 0, 0, 0.9, 0,
                    0, 0.7, 0, 0, 1.5,
                    0, 0, 0, 0, 0), 5, 5, byrow = T)


colnames(B5sparse) <- c("X1", "X2", "X3", "X4", "X5")

# specify layout
layout5 = matrix(c(0,1,
                   0,0,
                   1,-1,
                   2,0,
                   2,1),5,2,byrow = T)
par(mfrow=c(1,2))
true5psparse <- qgraph(t(B5sparse), layout=layout5, labels = colnames(B5sparse), theme="colorblind")

## Data generating
# equilibrium check
equilibrium_check(B5sparse)


## =============================== Oisin

mcnally <- read.csv("data/McNally.csv")

# separate depression / OCD symptoms
# (original data contains both depression and OCD symptoms)
# (here we only use depression symptoms)
depression <- mcnally[,1:16]
unique(depression)
#depression <- depression[,-c(4,6)]
ocd <- mcnally[,17:26]

trans_dep <- huge::huge.npn(depression)

trans_ocd <- huge::huge.npn(ocd)

skimr::skim(depression)
## estimate GGM via graphical LASSO on depression symptoms
cordep <- cor(depression)
# found the optimal sparsity by gamma = 1
glassoFitdep <- EBICglasso(cordep, n = nrow(depression), gamma = 1)
qgraph(glassoFitdep, layout = "spring", theme="colorblind", vsize=6,
       nodeNames = colnames(depression), legend.cex = 0.4)


## estimate the PAG on depression symptoms by running CCD
# run CCD
# ccd_mcnally_dep <- ccdKP(df=depression, dataType = "discrete", depth = -1, alpha=0.05)

ccd_mcnally_transdep <- ccdKP(df=trans_dep, dataType = "continuous", depth = -1, alpha=0.01)

ccd_mcnally_ocd <- ccdKP(df=ocd, dataType = "discrete", depth = -1)
# 
# ccd_trans_ocd <- ccdKP(df=trans_ocd, dataType = "continuous", depth = -1)

# create an adjacency matrix for PAG
mat_mcnally_dep <- CreateAdjMat(ccd_mcnally_dep, p = ncol(depression))
mat_trans_dep <- CreateAdjMat(ccd_mcnally_transdep, p = ncol(depression))
# mat_ocd <- CreateAdjMat(ccd_mcnally_ocd, p = ncol(ocd))
# mat_trans_ocd <- CreateAdjMat(ccd_mcnally_ocd, p = ncol(ocd))

# plot the PAG
pag_mcnally_dep <- plotPAG(ccd_mcnally_dep, mat_mcnally_dep)
pag_trans_dep <- plotPAG(ccd_mcnally_transdep, mat_trans_dep)
# pag_ocd <- plotPAG(ccd_mcnally_ocd, mat_ocd)
# pag_trans_ocd <- plotPAG(ccd_trans_ocd, mat_trans_ocd)


## define sufficient statistics 
suffStat <- list(dm = depression, nlev = rep(4, ncol(depression)), adaptDF = FALSE)

## estimate the PAG on depression symptoms by running CCI
cci(list(C = cor(mcnally), n = nrow(mcnally)), gaussCItest, alpha=0.01, labels = colnames(mcnally), p = ncol(mcnally)) %>% .$maag %>% plotAG


cci(list(dm = depression, nlev = rep(4, ncol(depression)), adaptDF = FALSE), disCItest, alpha=0.01, labels = colnames(depression), p = ncol(depression)) %>% .$maag %>% plotAG

cci(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha=0.01, labels = colnames(depression), p = ncol(depression)) %>% .$maag %>% plotAG

cci(list(C = cor(trans_dep), n = nrow(trans_dep)), gaussCItest, alpha=0.01, labels = colnames(trans_dep), p = ncol(trans_dep), verbose=TRUE) %>% .$maag %>% plotAG

cci(list(dm = ocd, nlev = rep(4, ncol(ocd)), adaptDF = FALSE), disCItest, alpha=0.01, labels = colnames(ocd), p = ncol(ocd)) %>% .$maag %>% plotAG

cci(list(C = cor(ocd), n = nrow(ocd)), gaussCItest, alpha=0.01, labels = colnames(ocd), p = ncol(ocd)) %>% .$maag %>% plotAG

cci(list(C = cor(trans_ocd), n = nrow(trans_ocd)), gaussCItest, alpha=0.01, labels = colnames(trans_ocd), p = ncol(trans_ocd)) %>% .$maag %>% plotAG

a <- pc(list(dm = depression, nlev = rep(4, ncol(depression)), adaptDF = FALSE), disCItest, alpha=0.01, labels = colnames(depression), p = ncol(depression))

a <- pc(list(C = cor(depression), n = nrow(depression)), gaussCItest, alpha=0.01, labels = colnames(depression), p = ncol(depression)) 

b <- pc(list(dm = depression, nlev = rep(4, ncol(depression)), adaptDF = FALSE), disCItest, alpha=0.01, labels = colnames(depression), p = ncol(depression))

b <- fci(list(C = cor(trans_dep), n = nrow(trans_dep)), gaussCItest, alpha=0.01, labels = colnames(trans_dep), verbose=T) 
b@amat %>% plotAG

plot(a@graph)
dev.off()



## huge.npn

library(ggh4x)


MyTheme2 <-  theme(plot.title = element_text(, family = "Palatino", size = 14, hjust=0.5),
                  plot.subtitle = element_text(face = "italic", family = "Palatino", size = 15, hjust=0.5),
                  axis.text=element_text(face = "bold",family = "Palatino", size = 11),
                  #axis.text.x = element_text(angle = 45, hjust = 1.2, vjust =1.2),
                  axis.title = element_text(face = "bold",family = "Palatino", size = 12),
                  legend.text = element_text(face = "bold", family = "Palatino", size = 12),
                  legend.position="bottom",
                  strip.text = element_text(size=12, family = "Palatino"),
                  strip.background = element_rect(fill="#f0f0f0", linetype = "solid", color="gray"),
                  strip.placement = "outside",
                  panel.border = element_rect(color = "#DCDCDC", fill = NA)
)


# original data
p1 <- depression %>%
  tidyr::pivot_longer(where(is.numeric)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~name) +
  theme_minimal() +
  ggtitle("(a) Original data") +
  MyTheme2

## transform data
# some are not normal --> makes the data semiparametric Gaussian using huge package
# following Sacha's advie on non-normal variable
transformed_dat <- huge::huge.npn(depression) %>% as.data.frame()
## check the distributions of transformed data
p2 <- transformed_dat %>%
  tidyr::pivot_longer(where(is.numeric)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~name) +
  theme_minimal() +
  ggtitle("(b) Transformed data") +
  MyTheme2

ggpubr::ggarrange(p1,p2)
ggsave(filename = "results/transformdat.pdf", width = 25, height = 13, dpi = 300, units = "cm")
