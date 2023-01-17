library(CCI)
source("code/R/plot_fnc.R")

###################
## basic 4 nodes ##
###################
suffStat_4p = list()
suffStat_4p$C = cor(data4p)
suffStat_4p$n = 1000

G4p = cci(suffStat_4p, gaussCItest, alpha=0.05, p=ncol(data4p)) 
plotAG(G4p$maag)

###################
## 4 nodes dense ##
###################
suffStat_4pdense = list()
suffStat_4pdense$C = cor(data4p_high)
suffStat_4pdense$n = 1000

G4p_high = cci(suffStat_4pdense, gaussCItest, alpha=0.05, p=ncol(data4p_high)) 
plotAG(G4p_high$maag)

####################
## 5 nodes sparse ##
####################
G5p = cci(suffStat_5p, gaussCItest, alpha=0.05, p=ncol(data5p)) 
plotAG(G5p$maag)


###################
## 5 nodes dense ##
###################

G5p_high = cci(suffStat_5pdense, gaussCItest, alpha=0.05, p=ncol(data5p_high)) 
plotAG(G5p_high$maag)

####################
## 6 nodes sparse ##
####################

G6p = cci(suffStat_6p, gaussCItest, alpha=0.05, p=ncol(data6p)) 
plotAG(G6p$maag)


###################
## 6 nodes dense ##
###################

G6p_high = cci(suffStat_6pdense, gaussCItest, alpha=0.05, p=ncol(data6p_high)) 
plotAG(G6p_high$maag)


########################
## mcnally depression ##
########################

G_mcdep = cci(suffStat_mcnallydep, disCItest, alpha=0.05, p=ncol(depression)) 
plotAG(G_mcdep$maag)


#################
## mcnally ocd ##
#################

G_mcocd = cci(suffStat_mcnallyocd, disCItest, alpha=0.05, p=ncol(ocd)) 
plotAG(G_mcocd$maag)
