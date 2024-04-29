#############################################        
###           Pour les colzas             ###
#############################################

librarian::shelf(tidyverse, ade4, factoextra)

# Chargement des données
## Masses individuelles
mybiom_Ind <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",")
## Date de germination
myBnGerm <- read.csv("data/derived-data/pl_Bn_germ.csv", h = T, sep = ",")
## C et N content
Ncont <- read.csv("data/derived-data/pl_Ncont.csv", h = T, sep = ",")  
Pcont <- read.csv("data/derived-data/pl_Pcont.csv", h = T, sep = ",")    

# Prepartion des matrices de données
Pcont_pca <- Pcont %>%
  filter(id_taxo == "Bra_napus")%>%
  pivot_longer(cols = (c("Leaves", "Limb", "Stem", "Fine.roots", "Taproot")), names_to = "organes",
                   values_to = "values") %>%
  select(-c(Aerial, X, id_taxo, numPL)) %>%
  mutate(organes = paste("Pcont", organes, sep = "_"))
Ncont_pca <- Ncont  %>%
  filter(id_taxo == "Bra_napus")%>%
  pivot_longer(cols = (c("Leaves", "Limb", "Stem", "Fine.roots", "Taproot")), names_to = "organes",
               values_to = "values") %>%
  mutate(organes = paste("Ncont", organes, sep = "_"))%>%
  select(-c(id_loca, X, "NA.", id_taxo, numPL))
biom_pca <- mybiom_Ind %>%
  filter(id_taxo == "Bra_napus")  %>%
  mutate(organes = paste("DW_", strata, sep = "")) %>%
  rename(values = DW)%>%
  select(-c(id_loca, X,FW, strata, id_taxo, numPL))
germ_pca <- myBnGerm %>%
  mutate(organes = "germ")%>%
  rename(values = SGDm, distance = dist)%>%
  select(-c(X))

# Prépartion matrice pour facteurs suppl dans l'ACP
biom_Lol <- mybiom_Ind %>%
  filter(id_taxo %in% c("Lol_peren", "Tri_alexa"))  %>%
  mutate(organes = case_when(strata == 0 & id_taxo == "Lol_peren" ~"Lp_aerial",
                             strata != 0 & id_taxo == "Lol_peren" ~"Lp_root",
                             strata == 0 & id_taxo == "Tri_alexa" ~"Ta_aerial",
                             strata != 0 & id_taxo == "Tri_alexa" ~"Ta_root")) %>%
  select(-c(id_loca, X,FW, strata, id_taxo, numPL)) %>%
  rename(values = DW) %>%
  group_by(RT,  div,  leg, serie, position, PAR, organes) %>%
  summarise(values = sum(values)) %>%
  pivot_wider(names_from = organes, values_from = values, values_fill = 0)


# Concatenation des matrices
prep_pca <- bind_rows(Ncont_pca, Pcont_pca, biom_pca, germ_pca) %>%
  pivot_wider(names_from = "organes", values_from = "values") %>%
  mutate_all(~replace_na(., 0)) %>%
  mutate(DW_R = rowSums(across(c("DW_1", "DW_2", "DW_3", "DW_4"))),
         ShootRoot = DW_0/(DW_R+0.001)) %>%
  na.omit() %>%
  left_join(biom_Lol)

# Réalisation de l'ACP
sel <-c("ShootRoot", "DW_0", "germ", "Ncont_Leaves", "Pcont_Leaves", "Ncont_Taproot", "Pcont_Taproot")
pca <- dudi.pca(prep_pca[1:121,sel], scannf = FALSE, nf = 2)

fviz_eig(pca)
   ## Predict coordinates and compute cos2
    quanti.coord <- supcol(pca, scale(prep_pca[1:121, c("Lp_aerial", "Lp_root", "Ta_aerial", "Ta_root")])) %>%
      .$cosup
    quanti.cos2 <- quanti.coord^2
    ## Graph of variables including supplementary variables
    p <- fviz_pca_var(pca)
    fviz_add(p, quanti.cos2, color ="blue", geom="arrow")
    fviz_pca_var(pca, col.var = "contrib")

    ##Graph of objects
    fviz_pca_ind(pca, 
             col.ind = as.factor(paste(prep_pca$leg, prep_pca$div, sep = "_")[1:121]), 
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             palette = c("dodgerblue1",  "dodgerblue4", "slateblue3", "slateblue4", "khaki", "khaki4", "burlywood4", "chocolate4"))

    res <- cbind(pca$li[,1:2], prep_pca[1:121, c("leg", "div", "distance")])
    
    ggplot(res, aes(Axis1, Axis2, 
                    fill= leg))+
      geom_point(aes(shape = as.factor(distance)))+
      stat_ellipse(geom = "polygon",
                   alpha = 0.25)+
      facet_wrap(~div)+
      theme_bw()


 

