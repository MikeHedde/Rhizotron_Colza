#############################################        
###           Pour les lolium             ###
#############################################

librarian::shelf(tidyverse, ade4)

# Chargement des données
## Masses individuelles
mybiom_Ind <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",")
## C et N content
Ncont <- read.csv("data/derived-data/pl_Ncont.csv", h = T, sep = ",")  
Pcont <- read.csv("data/derived-data/pl_Pcont.csv", h = T, sep = ",")    


Pcont_pca <- Pcont %>%
  filter(id_taxo == "Lol_peren")%>%
  pivot_longer(cols = (c("Leaves", "Limb", "Stem", "Fine.roots", "Taproot")), names_to = "organes",
               values_to = "values") %>%
  select(-c(Aerial, X, id_taxo, numPL)) %>%
  mutate(organes = paste("Pcont", organes, sep = "_"))
Ncont_pca <- Ncont  %>%
  filter(id_taxo == "Lol_peren")%>%
  pivot_longer(cols = (c("Leaves", "Limb", "Stem", "Fine.roots", "Taproot")), names_to = "organes",
               values_to = "values") %>%
  mutate(organes = paste("Ncont", organes, sep = "_"))%>%
  select(-c(id_loca, X, "NA.", id_taxo, numPL))
biom_pca <- mybiom_Ind %>%
  filter(id_taxo == "Lol_peren")  %>%
  mutate(organes = paste("DW_", strata, sep = "")) %>%
  rename(values = DW)%>%
  select(-c(id_loca, X,FW, strata, id_taxo, numPL))

sel <-c("leg", "div", "ShootRoot", "DW_0",  "Ncont_Leaves", "Pcont_Leaves")

prep_pca <- bind_rows(Ncont_pca, Pcont_pca, biom_pca) %>%
  pivot_wider(names_from = "organes", values_from = "values", values_fill = 0) %>%
  mutate_all(~replace_na(., 0)) %>%
  mutate(DW_R = rowSums(across(c("DW_1", "DW_2", "DW_3", "DW_4"))),
         ShootRoot = DW_0/(DW_R+0.001)) %>%
  na.omit() %>%
  select(all_of(sel))


pca <- dudi.pca(prep_pca[,-c(1:2)], scannf = FALSE, nf = 2)
s.corcircle(pca$co)
s.class(pca$li, as.factor(paste(prep_pca$leg, prep_pca$div, sep = "_")), 
        cellipse = F, col = 1:8)

pcares <- cbind(pca$li[, 1:2], prep_pca[,1:6]) 
ggplot(pcares, aes(x = div, y = Axis1))+
  geom_boxplot()+
  facet_grid(leg~.)
