#############################################        
###           Pour les colzas             ###
#############################################

librarian::shelf(tidyverse, ade4, factoextra)

# Chargement des données
## Masses individuelles
mybiom_Ind <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",")
myLDMC <-read.csv("data/raw-data/LDMC_v3.csv", h = T, sep = ";")
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
ldmc_pca <- myLDMC %>% 
  filter(id_taxo == "Bra_napus")  %>%
  select(-c(FW,DW, strata, id_taxo)) %>%
  separate(id_loca, remove = T, into = c("loc", "distance"), sep="d") %>%
  rename(leg = id_leg, div = id_div, RT = id_RT) %>%
  mutate(distance = as.numeric(distance)) %>%
  select(-c(loc, division))
biom_pca <- mybiom_Ind %>%
  filter(id_taxo == "Bra_napus")  %>%
  mutate(organes = paste("DW_", strata, sep = "")) %>%
  rename(values = DW)%>%
  select(-c(id_loca, X,FW, strata, id_taxo, numPL))
germ_pca <- myBnGerm %>%
  mutate(organes = "germ")%>%
  rename(values = SGDm, distance = dist)%>%
  select(-c(X))


# Concatenation des matrices
prep_pca <- bind_rows(Ncont_pca, Pcont_pca, biom_pca, germ_pca) %>%
  pivot_wider(names_from = "organes", values_from = "values") %>%
  mutate_all(~replace_na(., 0)) %>%
  mutate(DW_R = rowSums(across(c("DW_1", "DW_2", "DW_3", "DW_4"))),
         ShootRoot = DW_0/(DW_R+0.001)) %>%
  left_join(ldmc_pca)%>%
  na.omit() %>%
  mutate(categ = as.factor(paste(div,leg, sep="_")))

# K-tables
    ## ACPs pour réduire la dimensionalité des tableaux
      acp_biom <- dudi.pca(prep_pca[, c(18, 24:26)], scannf = F, nf = 2, scale = T)
      acp_Ncont <- dudi.pca(prep_pca[, c(8:12)], scannf = F, nf = 2, scale = F)
      acp_Pcont <- dudi.pca(prep_pca[, c(13:17)], scannf = F, nf = 2, scale = F)
      
      ACP <- acp_biom
      fviz_pca_var(ACP)
      fviz_pca_ind(ACP, 
                   col.ind = prep_pca$categ, 
                   addEllipses = TRUE, # Concentration ellipses
                   ellipse.type = "confidence",
                   palette = c("dodgerblue1",  "dodgerblue4", "slateblue3", "slateblue4", "khaki", "khaki4", "burlywood4", "chocolate4"))
      
      synthComp <- as.data.frame(cbind(
                         biom1 = acp_biom$li[,1],
                         #biom2 = acp_biom$li[,2],
                         Ncont1 = acp_Ncont$li[,1], 
                         #Ncont2 = acp_Ncont$li[,2],
                         Pcont1 = acp_Pcont$li[,1], 
                         #Pcont2 = acp_Pcont$li[,2],
                         germ = prep_pca$germ)
      )
      
    kta1 <- ktab.within(withinpca(synthComp, prep_pca$categ, scann = F, nf = 2))
    statis1 <- statis(kta1, scann = FALSE)
    plot(statis1)
