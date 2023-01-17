library(pcalg)

set.seed(123)

## 4p sparse
## fci from pcalg
suffStat_4p = list()
suffStat_4p$C = cor(data4p)
suffStat_4p$n = 1000

res4p <- fci(suffStat_4p,indepTest=gaussCItest,
    alpha = 0.95, p=p, doPdsep = FALSE)

plotAG(res4p@amat)


# ## fci using tetrad
# fci_4p <- tetradrunner(algoId = 'fci', df = data4p, 
#                        dataType = 'continuous', verbose=TRUE)
# fci_4p$graph
# fci_4p$edges
# fci_4p$graph$getUnderLines()
# 
# matfci_4p <- CreateAdjMat(fci_4p, 4)
# plotPAG(fci_4p, matfci_4p)


## 4p dense
suffStat_4pdense = list()
suffStat_4pdense$C = cor(data4p_high)
suffStat_4pdense$n = 1000

res4p_dense <- fci(suffStat_4pdense,indepTest=gaussCItest,
             alpha = 0.95, p=p, doPdsep = FALSE)

plotAG(res4p_dense@amat)

## 5p sparse
suffStat_5p = list()
suffStat_5p$C = cor(data5p)
suffStat_5p$n = 1000

res5p <- fci(suffStat_5p,indepTest=gaussCItest,
                   alpha = 0.95, p=p, doPdsep = FALSE)

plotAG(res5p@amat)

## 5p dense
suffStat_5pdense = list()
suffStat_5pdense$C = cor(data5p_high)
suffStat_5pdense$n = 1000

res5p_high <- fci(suffStat_5pdense,indepTest=gaussCItest,
             alpha = 0.95, p=p, doPdsep = FALSE)

plotAG(res5p_high@amat)


## 6p sparse
suffStat_6p = list()
suffStat_6p$C = cor(data6p)
suffStat_6p$n = 1000

res6p <- fci(suffStat_6p,indepTest=gaussCItest,
                  alpha = 0.95, p=p, doPdsep = FALSE)

plotAG(res6p@amat)


## 6p dense
suffStat_6pdense = list()
suffStat_6pdense$C = cor(data6p_high)
suffStat_6pdense$n = 1000

res6p_high <- fci(suffStat_6pdense,indepTest=gaussCItest,
             alpha = 0.95, p=p, doPdsep = FALSE)

plotAG(res6p_high@amat)


# ## run ccd on depression symptoms
# suffStat_mcnallydep <- list(dm = depression, nlev = rep(4, ncol(depression)), adaptDF = FALSE)
# # doesn't work: warning: "n=408 is too small (n < n.min = 1440 )"
# res_mcnallydep <- fci(suffStat_mcnallydep,indepTest=disCItest,
#                   alpha = 0.95, p = ncol(depression))
# 
# plotAG(res_mcnallydep@amat)

######################
## fci using tetrad ##
######################
## run fci on depression symptoms
fci_mcdep <- tetradrunner(algoId = 'fci', df = depression,
                       dataType = 'discrete', verbose=TRUE)
fci_mcdep$graph
fci_mcdep$edges
fci_mcdep$graph$getUnderLines()

matfci_mcdep <- CreateAdjMat(fci_mcdep, ncol(depression))
plotPAG(fci_mcdep, matfci_mcdep)


## run fci on ocd symptoms
fci_mcocd <- tetradrunner(algoId = 'fci', df = ocd,
                          dataType = 'discrete', verbose=TRUE)
fci_mcocd$graph
fci_mcocd$edges
fci_mcocd$graph$getUnderLines()

matfci_mcocd <- CreateAdjMat(fci_mcocd, ncol(ocd))
plotPAG(fci_mcocd, matfci_mcocd)
