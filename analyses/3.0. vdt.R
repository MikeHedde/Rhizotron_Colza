#####################"

## analyses biomasse plante / Rhizotron


# MZ
# 28/09/2016 : reprend SCRIPT_BIOM_PL_20160928.r au propre
# 6/12/2016 : complete stat effet traitement sur biomass pop
# 26/06/2017 : refait analyse apr?s elimination des individus trop petits (germintaion tardive)

# 11/08/2017 : reprends avec les stas de D. Makowsky
# 16/08/2017 : modif script + cr?ation de fonctions par MH

################################################################################

## library


# Source des fonctions stats et graphiques
source("analyses/1.0. fonctions stats.R")
source("analyses/1.1. fonctions graphiques.R")

# chargement des données
myvdt <- read.csv("data/derived-data/ew_com.csv", h=T, sep = ",")
myvdt_esp <- read.csv("data/derived-data/ew_sp.csv", h=T, sep = ",")


# repr?sentation de toutes les esp simultan?ment /!!\ (? am?liorer)
vdt_masssLoss <- ggplot(myvdt_esp, aes(x = id_taxo, y = -massLoss*100))+
  geom_boxplot()+
  facet_grid(div~leg)+
  labs(x = "Earthworm species", y = "Mass loss (%)")+
  theme_bw()

png("figures/ew/ew_massLoss_1boxplot.png", width = 1000, height = 480)
vdt_masssLoss
dev.off()

# A l'echelle communautaire
vdt_massLoss_T <- myfun_Com(mdata = myvdt, vb = myvdt$massLoss, seuil_prcRE = 30)
vdt_massLoss_T$select
pvdt_massLoss_T <- myplot_Com_vdt(mdata = myvdt, vb = -myvdt$massLoss, mod = vdt_massLoss_T$mod3A, 
                              ylabl = "Earthworm mass change (%FW mg)", plot_PAR = F)
pvdt_massLoss_T$moyFRic # pour voir les moyennes estim?es

png("figures/ew/ew_massLoss_2allEffects.png", width = 1000, height = 480)
grid.arrange(pvdt_massLoss_T$pDIV, pvdt_massLoss_T$pLEG, pvdt_massLoss_T$pINTER,
             ncol = 3)
dev.off()

# A l'?chelle sp?cifique
# A. icterica
myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_icter",]
vdt_massLoss_Ai <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
vdt_massLoss_Ai$select
pvdt_massLoss_Ai <- myplot_Com_vdt(mdata = myvdt_esp_temp, vb = -myvdt_esp_temp$massLoss, mod = vdt_massLoss_Ai$mod3A, 
                               ylabl = "A. icterica mass change (%FW mg)", plot_PAR = F)

png("figures/ew/ew_massLoss_Aict_Effects.png", width = 1000, height = 480)
grid.arrange(pvdt_massLoss_Ai$pDIV, pvdt_massLoss_Ai$pLEG, pvdt_massLoss_Ai$pINTER,
             ncol = 3)
dev.off()


# A. caliginosa
myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_calig",]
vdt_massLoss_Ac <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
vdt_massLoss_Ac$select
pvdt_massLoss_Ac <- myplot_Com_vdt(mdata = myvdt_esp_temp, vb = -myvdt_esp_temp$massLoss, mod = vdt_massLoss_Ac$mod3A, 
                               ylabl = "A. caliginosa mass change (%FW mg)", plot_PAR = T)

png("figures/ew/ew_massLoss_Acal_Effects.png", width = 1000, height = 480)
grid.arrange(pvdt_massLoss_Ac$pDIV, pvdt_massLoss_Ac$pLEG, pvdt_massLoss_Ac$pINTER,
             ncol = 3)
dev.off()

# L. terrestris
myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Lum_terre",]
vdt_massLoss_Lt <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 24)
vdt_massLoss_Lt$select
pvdt_massLoss_Lt <- myplot_Com_vdt(mdata = myvdt_esp_temp, vb = -myvdt_esp_temp$massLoss, mod = vdt_massLoss_Lt$mod3B, 
                                   ylabl = "L. terrestris mass change (%FW mg)", plot_PAR = F)

png("figures/ew/ew_massLoss_Lt_Effects.png", width = 1000, height = 480)
grid.arrange(pvdt_massLoss_Lt$pDIV, pvdt_massLoss_Lt$pLEG, pvdt_massLoss_Lt$pINTER,
             ncol = 3)
dev.off()

# A. longa
myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_longa",]
vdt_massLoss_Al <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
# pas possible car pr?sent dans un seul traitement  --> KW ??

# A. giardi
myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_giard",]
vdt_massLoss_Ag <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
# pas possible car pr?sent dans un seul traitement --> KW ??


#############################################        
###    Topologie du reseau de galeries    ###
#############################################   

# longueur totale des galeries
vdt_gal_T <- myfun_Com(mdata = myvdt, vb = myvdt$long, seuil_prcRE = 24)
vdt_gal_T$select
pvdt_gal_T <- myplot_Com_vdt(mdata = myvdt, vb = myvdt$long, mod = vdt_gal_T$mod3A, 
                             ylabl = "Total burrow length (m)", plot_PAR = F)
pvdt_gal_T$moyFRic # pour voir les moyennes estim?es
pvdt_gal_T$pDIV

png("figures/ew/ew_galLength.png", width = 1000, height = 480)
grid.arrange(pvdt_gal_T$pDIV, pvdt_gal_T$pLEG, pvdt_gal_T$pINTER,
             ncol = 3)
dev.off()


# longueur total des galeries
vdt_connec_T <- myfun_Com(mdata = myvdt, vb = myvdt$connec , seuil_prcRE = 24)
vdt_connec_T$select
pvdt_connec_T <- myplot_Com_vdt(mdata = myvdt, vb = myvdt$connec, mod = vdt_connec_T$mod3A, 
                                ylabl = "Burrow connectance", plot_PAR = T)
pvdt_connec_T$moyFRic # pour voir les moyennes estim?es
pvdt_connec_T$pDIV

# Specific Burrow length
myvdt$SBL_f <- myvdt$long/myvdt$FW_f
myvdt$SBL_mean <- myvdt$long/((myvdt$FW_f+myvdt$FW_i)/2)
vdt_SBL <- myfun_Com(mdata = myvdt, vb = myvdt$SBL_f , seuil_prcRE = 24)
vdt_SBL$select
pvdt_SBL <- myplot_Com_vdt(mdata = myvdt, vb = myvdt$SBL_f, mod = vdt_SBL$mod3B, 
                           ylabl = "Specific burrow length (m g-1)", plot_PAR = F)
pvdt_SBL$moyFRic # pour voir les moyennes estim?es
pvdt_SBL$pDIV


png("figures/ew/ew_galSBL.png", width = 1000, height = 480)
grid.arrange(pvdt_gal_T$pDIV, pvdt_gal_T$pLEG, pvdt_gal_T$pINTER,
             ncol = 3)
dev.off()
# -----------------------------------------------------------------------------