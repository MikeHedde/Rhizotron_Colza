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
library(tidyverse)

#------------------------------------------------------------
# 1 Chargement des donnees 
#------------------------------------------------------------
par <- read.csv("data/raw-data/rayonnement_plante.csv", sep = ";", dec=".", h=T)

vdt <- read.csv("data/raw-data/ew_weight.csv", sep=";", dec=".", header=T, na.strings = NA)
gal <- read.csv("data/raw-data/galeries_170828.csv", header = TRUE, sep = ";", dec=".", na.strings=NA)

enz <- read.csv("data/raw-data/Enz_180817.csv", h = T, sep = ";", dec = ".", na.strings = "")
anion <- read.csv("data/raw-data/anion.csv", sep=";", dec=".")

pl <- read.csv("data/raw-data/plant_weight.csv", h = T, sep = ";", dec = ".", na.strings = "")
Ncont <- read.csv("data/raw-data/N_160920.csv", h = T, sep = ";", dec = ".", na.strings = "")
Pcont <- read.csv("data/raw-data/P_180517.csv", h = T, sep = ";", dec = ".", na.strings = "")
lgrac <- read.csv("data/raw-data/root_length_170828.csv", h = T, sep = ";", dec = ".", na.strings = "") 
nod <- read.csv("data/raw-data/Talex_nod.csv", sep=";", dec=".")
BnGerm <- read.csv("data/raw-data/Bnap_GermDay.csv", sep=";", dec=".")

#------------------------------------------------------------
# 2 Calcul des biomasses vegetales seches
#-----------------------------------------------------------
# Masses aerienne et racinaires individuelle par strate
biom_Ind <- pl %>%
  group_by(div, leg, id_loca, id_taxo, distance, strata, RT) %>%
  summarise(DW = sum(DW, na.rm = T), FW = sum(FW, na.rm = T))

# Masse totale Communaute (aerien + racines)
biom_Com <- biom_Ind %>%
  group_by(div, leg, RT) %>%
  summarise(DW_tot = sum(DW, na.rm = T),
            DW_A = sum(DW[strata == 0], na.rm = T),
            DW_R = sum(DW[strata != 0], na.rm = T))
#-----------------------------------------------------------


#-----------------------------------------------------------
# 3    Mise en forme des df pour les traitements statistiques 
#     (1) ? l'?chelle de la communaut? (ajout facteurs s?rie,
#     position, PAR) pour les plantes, les activit?s enzymatiques 
#     et les vdt
#     (2) ? l'?chelle de l'individu pour les valeurs de biomasse 
#     des esp?ces v?g?tales et de [N] des organes de colza, 
#     ray-grass et tr?fle
#------------------------------------------------------------
# fichier source pour les codes des fonctions
source("analyses/1.0. fonctions stats.R")
  
# mise en forme pour les activites enzymatiques
myenz <- myprep(enz) 
write.csv(myenz, "data/derived-data/enz.csv")

# mise en forme pour les vdt (total)
vdt2 <- vdt %>%
  group_by(RT, div, leg) %>%
  summarise(FW_i = sum(FW[jour == 0]),
            FW_f = sum(FW[jour != 0])) %>%
  mutate(massLoss = (FW_i - FW_f)/FW_i) %>%
  left_join(gal)
myvdt <- myprep(vdt2) 
write.csv(myvdt, "data/derived-data/ew_com.csv")

# mise en forme pour les vdt (espèces)
vdt3 <- vdt %>%
  group_by(RT, div, leg, id_taxo) %>%
  summarise(FW_i = sum(FW[jour == 0]),
            FW_f = sum(FW[jour != 0])) %>%
  mutate(massLoss = (FW_i - FW_f)/FW_i)
myvdt_esp <- myprep(vdt3) 
write.csv(myvdt_esp, "data/derived-data/ew_sp.csv")

# mise en forme pour les plantes (total)
mybiom_Com <- myprep(biom_Com) 
write.csv(mybiom_Com, "data/derived-data/pl_biom_Com.csv")

# Masses vegetales seches 
mybiom_Ind_temp <- myprep(biom_Ind)
mybiom_Ind <- myloca(mybiom_Ind_temp)
write.csv(mybiom_Ind, "data/derived-data/pl_biom_Ind.csv")

# Date de germination du Colza modelisee           
myBnGerm <- myprep(BnGerm)
write.csv(myBnGerm, "data/derived-data/pl_Bn_germ.csv")

# Concentration en N dans les organes des plantes   
n_temp <- myprep(Ncont)
n <- myloca(n_temp)
myNcont <- n %>%
  mutate(organes = org) %>%
  mutate(organes = case_when(grepl('Leaf', organes) ~ "Limb",
                             organes %in% c("Other leaves", "Aerial") ~ "Leaves",
                             organes %in% c("Fine roots") ~"Fine roots",
                             organes %in% c("Taproot") ~"Taproot",
                             organes %in% c("Stem") ~"Stem")) %>%
  pivot_wider(id_cols = c(RT, leg, div, id_taxo, distance, id_loca, serie, position, PAR, numPL), 
              names_from = organes, values_from = N, values_fn = mean)
write.csv(myNcont, "data/derived-data/pl_Ncont.csv")

# Concentration en P dans les organes des plantes   
p_temp <- myprep(Pcont)
p <- myloca(p_temp)
myPcont <- p %>%
  rename(organes = org) %>%
  pivot_wider(id_cols = c(RT, leg, div, id_taxo, distance, serie, position, PAR, numPL), 
              names_from = organes, values_from = Pcont, values_fn = mean)
write.csv(myPcont, "data/derived-data/pl_Pcont.csv")


# Nodulations des trefles
Ta_temp <- myprep(nod)
ta <- myloca(Ta_temp)
# on ajoute les longueurs de racines dans la meme table
myTa <- ta %>%
  left_join(filter(mybiom_Ind, 
                   id_taxo == "Tri_alexa" & mybiom_Ind$strata != 0))
write.csv(myTa, "data/derived-data/pl_Ta_nod.csv")
# ----------------------------------------------------------------------------------------------