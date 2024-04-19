
#############################################        
###          Pour les v?roniques          ###
#############################################
# Biomasse s?che a?rienne
Vp_A <- mybiom_Ind[mybiom_Ind$id_taxo == "Ver_persi" & mybiom_Ind$strata == 0, ]
Vp_A_DW <- myfun_Ind(mdata = Vp_A, vb = log1p(Vp_A$DW), seuil_prcRE = 30)
Vp_A_DW$select
myplot_Ind(mdata = Vp_A, vb = log1p(Vp_A$DW), mod = Vp_A_DW$mod10A, 
           ylabl = "V. persica aerial biomass (log(g))", plot_PAR = T)