#############################################        
###          Pour les g?raniums           ###
#############################################
# Biomasse s?che a?rienne
Gd_A <- mybiom_Ind[mybiom_Ind$id_taxo == "Ger_disse" & mybiom_Ind$strata == 0, ]
Gd_A_DW <- myfun_Ind(mdata = Gd_A, vb = asin(tan(Gd_A$DW)), seuil_prcRE = 30) # transf log et sqrt ne nrmalisent pas les r?sidus, asin(tan()) pour l'exemple
Gd_A_DW$select
pGd_A_DW <- myplot_Ind(mdata = Gd_A, vb = Gd_A$DW, mod = Gd_A_DW$mod10B, 
                       ylabl = "L. perenne\n Shoot:Root ratio (DW)", plot_PAR = T)
pGd_A_DW$pDIV
Gd_R <- mybiom_Ind[mybiom_Ind$id_taxo == "Ger_disse" & mybiom_Ind$strata > 0, ]

