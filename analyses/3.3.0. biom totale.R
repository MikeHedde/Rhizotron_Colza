#############################################        
###            Biomasses seches           ###
############################################# 

# Chargement des données
mybiom_Com <- read.csv("data/derived-data/pl_biom_Com.csv", h = T, sep = ",")

# Source des fonctions stats et graphiques
source("analyses/1.0. fonctions stats.R")
source("analyses/1.1. fonctions graphiques.R")

# Comparaison des modeles testes et representations des differences
## Total
MC_DW_T <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_tot, seuil_prcRE = 30)
MC_DW_T$select
MC_DW_T <- myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_tot, mod = MC_DW_T$mod4A, 
           ylabl = "Total biomass of plant community (mg)", plot_PAR = T)
png("figures/pl/tot_mass_tot.png", width = 1000, height = 480)
grid.arrange(MC_DW_T$pDIV, MC_DW_T$pLEG, MC_DW_T$pINTER,
             ncol = 3)
dev.off()

## Racines
MC_DW_R <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_R, seuil_prcRE = 30)
MC_DW_R$select
MC_DW_R <- myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_R, mod = MC_DW_R$mod4A, 
           ylabl = "Root biomass of plant community (mg)", plot_PAR = T)
png("figures/pl/tot_mass_root.png", width = 1000, height = 480)
grid.arrange(MC_DW_R$pDIV, MC_DW_R$pLEG, MC_DW_R$pINTER,
             ncol = 3)
dev.off()

## Parties aériennes
MC_DW_A <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_A, seuil_prcRE = 30)
MC_DW_A$select
MC_DW_A <- myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_A, mod = MC_DW_A$mod4A, 
           ylabl = "Shoot biomass of plant community (mg)", plot_PAR = T)
png("figures/pl/tot_mass_shoot.png", width = 1000, height = 480)
grid.arrange(MC_DW_A$pDIV, MC_DW_A$pLEG, MC_DW_A$pINTER,
             ncol = 3)
dev.off()
# ---------------------------------------------------------------------------------