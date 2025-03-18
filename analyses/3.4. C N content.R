#-----------------------------------------------------------
# Effet des traitements a l'echelle specifique:   
#         [N] dans les organes de colza, ray-grass et trefle
#------------------------------------------------------------

librarian::shelf(ggplot2)

source("analyses/1.0. fonctions stats.R")
source("analyses/1.1. fonctions graphiques.R")

# Chargement données
Ncont <- read.csv("data/derived-data/pl_Ncont.csv")


#############################################        
###           Pour les trefles            ###
############################################# 

Ta_N <- Ncont[Ncont$id_taxo == "Tri_alexa" , ]
# dans les racines fines
Ta_N_temp <- Ta_N[!is.na(Ta_N$Fine.roots),]
Ta_N_fRoots <- myfun_Leg(mdata = Ta_N_temp, vb = log(Ta_N_temp$Fine.roots), seuil_prcRE = 30) # pas de transf efficaces
Ta_N_fRoots$resum
Ta_N_fRoots$select
pTa_N_fRoots <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Fine.roots, mod = Ta_N_fRoots$mod10A, 
                           ylabl = "T. alexandrium \n roots N content (mg g-1)", plot_PAR = T)
pTa_N_fRoots$pDIV
# dans les limbes matures
Ta_N_temp <- Ta_N[!is.na(Ta_N$Limb),]
Ta_N_limb <- myfun_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Limb, seuil_prcRE = 24) # pas de transf 
Ta_N_limb$resum
Ta_N_limb$select
pTa_N_limb <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Limb, mod = Ta_N_limb$mod10A, 
                         ylabl = "T. alexandrium LNC (mg g-1)", plot_PAR = F)
pTa_N_limb$pDIV
# dans les limbes d'autres feuilles
Ta_N_temp <- Ta_N[!is.na(Ta_N$Leaves),]
Ta_N_leaves <- myfun_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Leaves, seuil_prcRE = 30) # pas de transf 
Ta_N_leaves$resum
pTa_N_leaves <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Leaves, mod = Ta_N_leaves$mod9A, 
                           ylabl = "T. alexandrium \nleaves N content (mg g-1)", plot_PAR = T)
pTa_N_leaves$pDIV
# dans les tiges
Ta_N_temp <- Ta_N[!is.na(Ta_N$Stem),]
Ta_N_stem <- myfun_Leg(mdata = Ta_N_temp, vb = asin(tan(Ta_N_temp$Stem)), seuil_prcRE = 30) 
Ta_N_stem$resum
pTa_N_stem <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Stem, mod = Ta_N_stem$mod9A, 
                         ylabl = "T. alexandrium \nstem N content (mg g-1)", plot_PAR = T)             
pTa_N_stem$pDIV
# comparaison de la distribution des concentrations
Ta_meaN <- rbind.data.frame(pTa_N_limb$moyFRic, pTa_N_leaves$moyFRic, pTa_N_stem$moyFRic, pTa_N_fRoots$moyFRic)
Ta_meaN$organes <- rep(c("a.limb", "b.leaves", "c.stem", "d.roots"), each = 4)
B <- ggplot(data = Ta_meaN, aes(x = organes, y = x))+
  geom_point(aes(shape = div), size = 4)+
  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
  scale_shape_discrete(guide = FALSE)+
  ylab("N content (mg g-1)")   +
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  facet_grid(div~.)
B
#############################################        
###           Pour les colzas             ###
#############################################
Bn_N <- Ncont[Ncont$id_taxo == "Bra_napus" , ]
# dans les racines fines
Bn_N_temp <- Bn_N[!is.na(Bn_N$Fine.roots),]
Bn_N_fRoots <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Fine.roots, seuil_prcRE = 24) 
Bn_N_fRoots$resum # pas de transf efficaces
Bn_N_fRoots$select
pBn_N_fRoots <- myplot_Ind(mdata = Bn_N_temp, vb = log(Bn_N_temp$Fine.roots), mod = Bn_N_fRoots$mod10A, 
                           ylabl = "B. napus \nfine roots N content (mg g-1)", plot_PAR = T)
# v?rifier NA's

# dans les racines pivot
Bn_N_temp <- Bn_N[!is.na(Bn_N$Taproot),]
Bn_N_tRoots <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Taproot, seuil_prcRE = 30)
Bn_N_tRoots$resum
Bn_N_tRoots$select
pBn_N_tRoots <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Taproot, mod = Bn_N_tRoots$mod10A, 
                           ylabl = "B. napus \ntaproots N content (mg g-1)", plot_PAR = T)
pBn_N_tRoots$pDIV_INTER

# dans les limbes matures
Bn_N_temp <- Bn_N[!is.na(Bn_N$Limb),]
Bn_N_limb <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Limb, seuil_prcRE = 24) # pas de transf 
Bn_N_limb$resum
Bn_N_limb$select
pBn_N_limb <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Limb, mod = Bn_N_limb$mod10A, 
                         ylabl = "B. napus mature LNC (mg g-1)", plot_PAR = F)
pBn_N_limb$pDIV_INTER

# dans les limbes d'autres feuilles
Bn_N_temp <- Bn_N[!is.na(Bn_N$Leaves),]
Bn_N_leaves <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Leaves, seuil_prcRE = 30) # pas de transf 
Bn_N_leaves$resum
Bn_N_leaves$select
pBn_N_leaves <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Leaves, mod = Bn_N_leaves$mod10A, 
                           ylabl = "B. napus \nleaves N content (mg g-1)", plot_PAR = T)

# dans les tiges
Bn_N_temp <- Bn_N[!is.na(Bn_N$Stem),]
Bn_N_stem <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Stem, seuil_prcRE = 30) 
Bn_N_stem$resum
Bn_N_stem$select
pBn_N_stem <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Stem, mod = Bn_N_stem$mod10A, 
                         ylabl = "B. napus \nstem N content (mg g-1)", plot_PAR = T)             

# comparaison de la distribution des concentrations
Bn_meaN <- rbind.data.frame(pBn_N_limb$moyFRic, pBn_N_leaves$moyFRic, pBn_N_stem$moyFRic, pBn_N_tRoots$moyFRic, pBn_N_fRoots$moyFRic)
Bn_meaN$organes <- as.factor(rep(c("a.limb", "b.leaves", "c.stem", "d.taproots", "e.fine roots"), each = 4))
D <- ggplot(data = Bn_meaN, aes(x = organes, y = x))+
  geom_point(aes(shape = div), size = 4)+
  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
  scale_shape_discrete(guide = FALSE)+
  ylab("N content (mg g-1)")   +
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  facet_grid(div~.)
D


#############################################        
###           Pour les ray-grass          ###
#############################################
Lp_N <- Ncont[Ncont$id_taxo == "Lol_peren" , ]
# dans les racines fines
Lp_N_temp <- Lp_N[!is.na(Lp_N$Fine.roots),]
Lp_N_fRoots <- myfun_Ind(mdata = Lp_N_temp, vb = sqrt(Lp_N_temp$Fine.roots), seuil_prcRE = 30) 
Lp_N_fRoots$resum # pas de transf efficaces
Lp_N_fRoots$select
# test avec les options de filtrage
Lp_N_fRoots_opt2 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
Lp_N_fRoots_opt4 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
Lp_N_fRoots <- myfun_Ind(mdata = Lp_N_fRoots_opt2, vb = sqrt(Lp_N_fRoots_opt2$Fine.roots), seuil_prcRE = 30)
# aucune n'est satisfaisante => Kruskal-Wallis ?
pLp_N_fRoots <- myplot_Ind(mdata = Lp_N_temp, vb = Lp_N_temp$Fine.roots, mod = Lp_N_fRoots$mod10A, 
                           ylabl = "L. perenne RNC (mg g-1)", plot_PAR = F)
pLp_N_fRoots$pDIV

# dans les limbes matures
Lp_N_temp <- Lp_N[!is.na(Lp_N$Limb),]
Lp_N_limb <- myfun_Ind(mdata = Lp_N_temp, vb = log(Lp_N_temp$Limb), seuil_prcRE = 24) # pas de transf efficaces
Lp_N_limb$resum
Lp_N_limb$select
# test avec les options de filtrage
#Lp_N_limb_opt2 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
#Lp_N_limb_opt4 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
#Lp_N_limb <- myfun_Ind(mdata = Lp_N_limb_opt4, vb = sqrt(Lp_N_limb_opt4$Fine.roots), seuil_prcRE = 30)
# aucune n'est satisfaisante => Kruskal-Wallis ?
pLp_N_limb <- myplot_Ind(mdata = Lp_N_temp, vb = Lp_N_temp$Limb, mod = Lp_N_limb$mod10A, 
                         ylabl = "L. perenne LNC (mg g-1)", plot_PAR =F)
pLp_N_limb$pDIV

# dans les limbes d'autres feuilles
Lp_N_temp <- Lp_N[!is.na(Lp_N$Leaves),]
Lp_N_leaves <- myfun_Ind(mdata = Lp_N_temp, vb = Lp_N_temp$Leaves, seuil_prcRE = 24) # pas de transf 
Lp_N_leaves$resum
Lp_N_leaves$select
# test avec les options de filtrage
#Lp_N_leaves_opt2 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
#Lp_N_leaves_opt4 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
#Lp_N_leaves <- myfun_Ind(mdata = Lp_N_leaves_opt2, vb = sqrt(Lp_N_leaves_opt2$Fine.roots), seuil_prcRE = 30)
# aucune n'est satisfaisante => Kruskal-Wallis ?
pLp_N_leaves <- myplot_Ind(mdata = Lp_N_temp[!is.na(Lp_N_temp$Leaves),], vb = Lp_N_temp$Leaves[!is.na(Lp_N_temp$Leaves)], mod = Lp_N_leaves$mod10A, 
                           ylabl = "L. perenne leaves N content (mg g-1)", plot_PAR =F)





#-----------------------------------------------------------
# 3.3 Effet des traitements ? l'?chelle sp?cifique:   
#         [P] dans les organes de colza, ray-grass et tr?fle
#------------------------------------------------------------
#############################################        
###           Pour les tr?fles            ###
############################################# 
Ta_P <- myPcont[myPcont$id_taxo == "Tri_alexa" , ]
# dans les racines fines
Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Fine.roots) & Ta_P$`Fine roots` > 0.1),]
Ta_P_fRoots <- myfun_Leg(mdata = Ta_P_temp, vb = log(Ta_P_temp$Fine.roots), seuil_prcRE = 30) # pas de transf efficaces
Ta_P_fRoots$resum
Ta_P_fRoots$select
pTa_P_fRoots <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Fine.roots, mod = Ta_P_fRoots$mod10A, 
                           ylabl = "T. alexandrium \n roots P content (mg g-1)", plot_PAR = T)
pTa_P_fRoots$pDIV

# dans les limbes matures
Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Limb) & Ta_P$Limb > 0.01),]
Ta_P_limb <- myfun_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Limb, seuil_prcRE = 24) # pas de transf 
Ta_P_limb$resum
Ta_P_limb$select
pTa_P_limb <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Limb, mod = Ta_P_limb$mod10A, 
                         ylabl = "T. alexandrium LPC (mg g-1)", plot_PAR = F)

# dans les limbes d'autres feuilles
Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Leaves) & Ta_P$Leaves > 0.1),]
Ta_P_leaves <- myfun_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Leaves, seuil_prcRE = 30) # pas de transf 
Ta_P_leaves$resum
pTa_P_leaves <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Leaves, mod = Ta_P_leaves$mod9A, 
                           ylabl = "T. alexandrium \nleaves P content (mg g-1)", plot_PAR = T)
pTa_P_leaves$pDIV

# dans les tiges
Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Stem) & Ta_P$Stem > 0.1),]
Ta_P_stem <- myfun_Leg(mdata = Ta_P_temp, vb = asin(tan(Ta_P_temp$Stem)), seuil_prcRE = 30) 
Ta_P_stem$resum
pTa_P_stem <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Stem, mod = Ta_P_stem$mod9A, 
                         ylabl = "T. alexandrium \nstem N content (mg g-1)", plot_PAR = T)             

# comparaison de la distribution des concentrations
Ta_meaN <- rbind.data.frame(pTa_P_limb$moyFRic, pTa_P_leaves$moyFRic, pTa_P_stem$moyFRic, pTa_P_fRoots$moyFRic)
Ta_meaN$organes <- rep(c("a.limb", "b.leaves", "c.stem", "d.roots"), each = 4)
B <- ggplot(data = Ta_meaN, aes(x = div, y = x))+
  geom_point(aes(shape = div), size = 4)+
  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
  scale_shape_discrete(guide = FALSE)+
  ylab("P content (mg g-1)")   +
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  facet_grid(organes~.)
B
#############################################        
###           Pour les colzas             ###
#############################################
Bn_P <- myPcont[myPcont$id_taxo == "Bra_napus" , ]
# dans les racines fines
Bn_P_temp <- Bn_P[!is.na(Bn_P$Fine.roots),]
Bn_P_fRoots <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Fine.roots, seuil_prcRE = 24) 
Bn_P_fRoots$resum # pas de transf efficaces
Bn_P_fRoots$select
# test avec les options de filtrage
#Bn_P_fRoots_opt8 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt8==1]), ]
#Bn_P_fRoots_opt7 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt7==1]), ]
#Bn_P_fRoots_opt6 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt6==1]), ]
#Bn_P_fRoots_opt4 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt4==1]), ]
#Bn_P_fRoots <- myfun_Ind(mdata = Bn_P_fRoots_opt4, vb = log(Bn_P_fRoots_opt4$Fine.roots), seuil_prcRE = 30)  # il faut l'option 4 pour r?ussir ? obtenir des r?sidus normaux
#Bn_P_fRoots$resum
#Bn_P_fRoots$select
pBn_P_fRoots <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Fine.roots, mod = Bn_P_fRoots$mod10A, 
                           ylabl = "B. napus \nfine roots P content (mg g-1)", plot_PAR = T)
pBn_P_fRoots$pDIV_INTER

# dans les racines pivot
Bn_P_temp <- Bn_P[!is.na(Bn_P$Taproot),]
Bn_P_tRoots <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Taproot, seuil_prcRE = 30)
Bn_P_tRoots$resum
Bn_P_tRoots$select
pBn_P_tRoots <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Taproot, mod = Bn_P_tRoots$mod10A, 
                           ylabl = "B. napus \ntaproots P content (mg g-1)", plot_PAR = T)
pBn_P_tRoots$pDIV

# dans les limbes matures
Bn_P_temp <- Bn_P[!is.na(Bn_P$Limb),]
Bn_P_limb <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Limb, seuil_prcRE = 24) # pas de transf 
Bn_P_limb$resum
Bn_P_limb$select
pBn_P_limb <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Limb, mod = Bn_P_limb$mod10A, 
                         ylabl = "B. napus mature LPC (mg g-1)", plot_PAR = F)
pBn_P_limb$pDIV_INTER

# dans les limbes d'autres feuilles
Bn_P_temp <- Bn_P[!is.na(Bn_P$Leaves),]
Bn_P_leaves <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Leaves, seuil_prcRE = 30) # pas de transf 
Bn_P_leaves$resum
Bn_P_leaves$select
pBn_P_leaves <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Leaves, mod = Bn_P_leaves$mod10A, 
                           ylabl = "B. napus \nleaves P content (mg g-1)", plot_PAR = T)

# dans les tiges
Bn_P_temp <- Bn_P[!is.na(Bn_P$Stem),]
Bn_P_stem <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Stem, seuil_prcRE = 30) 
Bn_P_stem$resum
Bn_P_stem$select
pBn_P_stem <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Stem, mod = Bn_P_stem$mod10A, 
                         ylabl = "B. napus \nstem P content (mg g-1)", plot_PAR = T)             

# comparaison de la distribution des concentrations
Bn_meaN <- rbind.data.frame(pBn_P_limb$moyFRic, pBn_P_leaves$moyFRic, pBn_P_stem$moyFRic, pBn_P_tRoots$moyFRic, pBn_P_fRoots$moyFRic)
Bn_meaN$organes <- as.factor(rep(c("a.limb", "b.leaves", "c.stem", "d.taproots", "e.fine roots"), each = 4))
D <- ggplot(data = Bn_meaN, aes(x = div, y = x))+
  geom_point(aes(shape = div), size = 4)+
  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
  scale_shape_discrete(guide = FALSE)+
  ylab("P content (mg g-1)")   +
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  facet_grid(.~organes)
D


#############################################        
###           Pour les ray-grass          ###
#############################################
Lp_P <- myPcont[myPcont$id_taxo == "Lol_peren" , ]
# dans les racines fines
Lp_P_temp <- Lp_P[!is.na(Lp_P$Fine.roots),]
Lp_P_fRoots <- myfun_Ind(mdata = Lp_P_temp, vb = sqrt(Lp_P_temp$Fine.roots), seuil_prcRE = 30) 
Lp_P_fRoots$resum # pas de transf efficaces
Lp_P_fRoots$select
pLp_P_fRoots <- myplot_Ind(mdata = Lp_P_temp, vb = Lp_P_temp$Fine.roots, mod = Lp_P_fRoots$mod10A, 
                           ylabl = "L. perenne RPC (mg g-1)", plot_PAR = F)
pLp_P_fRoots$pDIV

# dans les limbes matures
Lp_P_temp <- Lp_P[!is.na(Lp_P$Limb),]
Lp_P_limb <- myfun_Ind(mdata = Lp_P_temp, vb = log(Lp_P_temp$Limb), seuil_prcRE = 24) # pas de transf efficaces
Lp_P_limb$resum
Lp_P_limb$select
pLp_P_limb <- myplot_Ind(mdata = Lp_P_temp, vb = Lp_P_temp$Limb, mod = Lp_P_limb$mod10A, 
                         ylabl = "L. perenne LPC (mg g-1)", plot_PAR =F)

# dans les limbes d'autres feuilles
Lp_P_temp <- Lp_P[!is.na(Lp_P$Leaves),]
Lp_P_leaves <- myfun_Ind(mdata = Lp_P_temp, vb = Lp_P_temp$Leaves, seuil_prcRE = 24) # pas de transf 
Lp_P_leaves$resum
Lp_P_leaves$select
# test avec les options de filtrage
#Lp_P_leaves_opt2 <- Lp_P_temp[which(Lp_P_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
#Lp_P_leaves_opt4 <- Lp_P_temp[which(Lp_P_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
#Lp_P_leaves <- myfun_Ind(mdata = Lp_P_leaves_opt2, vb = sqrt(Lp_P_leaves_opt2$Fine.roots), seuil_prcRE = 30)
# aucune n'est satisfaisante => Kruskal-Wallis ?
pLp_P_leaves <- myplot_Ind(mdata = Lp_P_temp[!is.na(Lp_P_temp$Leaves),], vb = Lp_P_temp$Leaves[!is.na(Lp_P_temp$Leaves)], mod = Lp_P_leaves$mod10A, 
                           ylabl = "L. perenne leaves P content (mg g-1)", plot_PAR =F)

pLp_P_leaves$pDIV_INTER

# comparaison de la distribution des concentrations
Lp_meaN <- rbind.data.frame(pLp_P_limb$moyFRic, pLp_P_leaves$moyFRic, pLp_P_fRoots$moyFRic)
Lp_meaN$organes <- as.factor(rep(c("a.limb", "b.leaves", "e.fine roots"), each = 4))
E <- ggplot(data = Lp_meaN, aes(x = div, y = x))+
  geom_point(aes(shape = div), size = 4)+
  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
  scale_shape_discrete(guide = FALSE)+
  ylab("P content (mg g-1)")   +
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  facet_grid(.~organes)
E

#-----------------------------------------------------------
# 3.4 Effet des traitements ? l'?chelle sp?cifique:   
#        quantit? de N dans les organes de colza, ray-grass et tr?fle
#------------------------------------------------------------                             
Pcont <- read.csv("data/derived-data/pl_Pcont.csv")
Pcont$qN <- Pcont$DM*Pcont$Ncont/100
qttN <- aggregate(qN~distance+RT+leg+div+id_taxo, sum, data = Pcont)
myqttN <- myprep(qttN)

#############################################        
###           Pour les colzas             ###
#############################################
Bn_qN <- myqttN[myqttN$id_taxo == "Bra_napus" , ]
moyFi <- aggregate(Bn_qN$qN, by = list(div = Bn_qN$div, leg = Bn_qN$leg), mean, na.rm=T)
seFi <- aggregate(Bn_qN$qN, by = list(div = Bn_qN$div, leg = Bn_qN$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("qN") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER

#############################################        
###           Pour les ray-grass          ###
#############################################

Lp_qN <- myqttN[myqttN$id_taxo == "Lol_peren" , ]
moyFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div, leg = Ta_qN$leg), mean, na.rm=T)
seFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div, leg = Ta_qN$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("qN") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER

#############################################        
###           Pour les tr?fles            ###
#############################################
Ta_qN <- myqttN[myqttN$id_taxo == "Tri_alexa" , ]
moyFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div), mean, na.rm=T)
seFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = div, y=x)) + 
  geom_bar(aes(fill = div), position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("qN") +    
  xlab("") +
  scale_fill_manual (values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                     name = "",
                     labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER



#############################################        
###           Bilan par RT                ###
#############################################
moyFi_0 <- aggregate(qttN$qN, by = list(div = qttN$div, leg = qttN$leg, RT = qttN$RT), sum, na.rm=T)
moyFi <- aggregate(moyFi_0$x, by = list(div = moyFi_0$div, leg = moyFi_0$leg), mean, na.rm=T)
seFi <- aggregate(moyFi_0$x, by = list(div = moyFi_0$div, leg = moyFi_0$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("Production de N totale par RT (mg)") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER

#-----------------------------------------------------------
# 3.5 Effet des traitements ? l'?chelle sp?cifique:   
#        quantit? de P dans les organes de colza, ray-grass et tr?fle
#------------------------------------------------------------  
Pcont$qP <- Pcont$DM*Pcont$Pcont/100
qttP <- aggregate(qP~distance+RT+leg+div+id_taxo, sum, data = Pcont)
myqttP <- myprep(qttP)


#############################################        
###           Pour les colzas             ###
#############################################
Bn_qP <- myqttP[myqttP$id_taxo == "Bra_napus" , ]
moyFi <- aggregate(Bn_qP$qP, by = list(div = Bn_qP$div, leg = Bn_qP$leg), mean, na.rm=T)
seFi <- aggregate(Bn_qP$qP, by = list(div = Bn_qP$div, leg = Bn_qP$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("qP") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER


#############################################        
###           Pour les ray-grass          ###
#############################################
Lp_qP <- myqttP[myqttP$id_taxo == "Lol_peren" , ]
moyFi <- aggregate(Lp_qP$qP, by = list(div = Lp_qP$div, leg = Lp_qP$leg), mean, na.rm=T)
seFi <- aggregate(Lp_qP$qP, by = list(div = Lp_qP$div, leg = Lp_qP$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("qP") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER

#############################################        
###           Pour les tr?fles            ###
#############################################
Ta_qP <- myqttP[myqttP$id_taxo == "Tri_alexa" , ]
moyFi <- aggregate(Ta_qP$qP, by = list(div = Ta_qP$div), mean, na.rm=T)
seFi <- aggregate(Ta_qP$qP, by = list(div = Ta_qP$div), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = div, y=x)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("qP") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER


#-----------------------------------------------------------
# 3.6 Effet des traitements ? l'?chelle sp?cifique:   
#        rapport N/P dans les organes de colza, ray-grass et tr?fle
#------------------------------------------------------------  
myqttN$code <- paste(myqttN$RT, myqttN$distance, myqttN$leg, myqttN$div, myqttN$id_taxo, myqttN$serie, myqttN$position, myqttN$PAR)
myqttP$code <- paste(myqttP$RT, myqttP$distance, myqttP$leg, myqttP$div, myqttP$id_taxo, myqttP$serie, myqttP$position, myqttP$PAR)

myqttNP <- merge(myqttN, myqttP, by = "code")
myqttNP <- myqttNP[, c(2:7, 16)]
colnames(myqttNP) <- c("RT",  "distance", "leg","div", "id_taxo", "qN", "qP" )       
myqttNP$rNP <- myqttNP$qN/myqttNP$qP

#############################################        
###           Pour les colzas             ###
#############################################
Bn_rNP <- myqttNP[myqttNP$id_taxo == "Bra_napus" , ]
moyFi <- aggregate(Bn_rNP$rNP, by = list(div = Bn_rNP$div, leg = Bn_rNP$leg), mean, na.rm=T)
seFi <- aggregate(Bn_rNP$rNP, by = list(div = Bn_rNP$div, leg = Bn_rNP$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("rNP") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER


#############################################        
###           Pour les ray-grass          ###
#############################################
Lp_rNP <- myqttNP[myqttNP$id_taxo == "Lol_peren" , ]
moyFi <- aggregate(Lp_rNP$rNP, by = list(div = Lp_rNP$div, leg = Lp_rNP$leg), mean, na.rm=T)
seFi <- aggregate(Lp_rNP$rNP, by = list(div = Lp_rNP$div, leg = Lp_rNP$leg), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("rNP") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER

#############################################        
###           Pour les tr?fles            ###
#############################################
Ta_rNP <- myqttNP[myqttNP$id_taxo == "Tri_alexa" , ]
moyFi <- aggregate(Ta_rNP$rNP, by = list(div = Ta_rNP$div), mean, na.rm=T)
seFi <- aggregate(Ta_rNP$rNP, by = list(div = Ta_rNP$div), std.error, na.rm=T)
moyFi$se <- seFi$x
DIV_INTER <- ggplot(data = moyFi, aes(x = div, y=x)) + 
  geom_bar(position = position_dodge(width = 0.8), stat="identity")+
  geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
  ylab("rNP") +    
  xlab("") +
  scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                    name = "",
                    labels = c("Ew0", "FR1", "FR2", "FR3"))+
  theme_classic() +
  theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  theme(axis.text = element_text(color="black", size=14)) +
  theme(axis.title = element_text(color="black", size=14))+
  scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))

DIV_INTER

