librarian::shelf(dplyr, ggplot2, Rmisc, tidyr, stringr, ggpubr)


# masse
    ## téléchargement et préparation des données
    Bn <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
      filter(id_taxo == "Bra_napus") %>%
      select(div, leg, strata, DW) %>% 
      pivot_longer(cols = c(DW), 
                   names_to = "variable", values_to = "value") %>%
      summarySE(measurevar=c("value"), na.rm = T,
                groupvars=c("leg", "div", "strata", "variable")) %>%
      mutate(strata = case_when(strata == 0 ~ 10,
                                strata == 1 ~-3.5,
                                strata == 2 ~-18.5,
                                strata == 3 ~-40.5,
                                strata == 4 ~-70.5)) %>%
      mutate(variable = str_replace(variable, "DW", "Dry mass (g)"))
    
    # représentation graphique
    Bn_biom_p <- ggplot(Bn, aes(x = strata, y = value, colour = div))+
      geom_point(position=position_dodge(width= 2), size = 2)+
      geom_line(alpha = 0.2)+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                    width=.1, position=position_dodge(width= 2)) + 
      scale_y_continuous(trans='log10')+
      geom_vline(xintercept = 0)+
      labs(x= "Soil depth (cm)", y = "Mass (g; log scale)", colour = "Earthworm\ndiversity")+
      coord_flip()+
      facet_grid(~leg, scales = "free")+
      theme_bw()

# téléchargement et préparation des données sur les concentrations en N
Ncont <- read.csv("data/derived-data/pl_Ncont.csv") %>%
  filter(id_taxo == "Bra_napus") %>%
  select(div, leg, Fine.roots, Limb, Leaves, Stem, Taproot) %>%
  pivot_longer(cols = c(Fine.roots, Limb, Leaves, Stem, Taproot), 
               names_to = "organ", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("leg", "div", "organ")) %>%
  filter(!is.na(value)) %>%
  mutate(trait = "N content  (mg.g-1)")

# téléchargement et préparation des données sur les concentrations en P
Pcont <- read.csv("data/derived-data/pl_Pcont.csv") %>%
  filter(id_taxo == "Bra_napus") %>%
  select(div, leg, Fine.roots, Limb, Leaves, Stem, Taproot) %>%
  pivot_longer(cols = c(Fine.roots, Leaves), 
               names_to = "organ", values_to = "value") %>%
  mutate(value = ifelse(value<0,0.001,value)) %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("leg", "div", "organ")) %>%
  filter(!is.na(value)) %>%
  mutate(trait = "P content (mg.g-1)")

# téléchargement et préparation des données de germination du colza
Bn_germ <- read.csv("data/derived-data/pl_Bn_germ.csv", h = T, sep = ",") %>%
  select(div, leg, SGDm) %>%
  pivot_longer(cols = SGDm, 
               names_to = "variable", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("div", "leg")) %>%
  filter(!is.na(value)) %>%
  mutate(trait = "Germination \n(days after sowing)", 
         organ = "")
  
# téléchargement et préparation des données de Leaf dry matter content 
Bn_ldmc <- read.csv("data/raw-data/LDMC_v3.csv", h = T, sep = ";") %>%
  filter(id_taxo == "Bra_napus") %>%
  select(div, leg, LDMC) %>%
  pivot_longer(cols = LDMC, 
               names_to = "variable", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("div", "leg")) %>%
  filter(!is.na(value))%>%
  mutate(trait = "Leaf dry matter \ncontent (mg.g-1)",
         organ = "")

# profondeur max d'enracinement
##téléchargement et préparation des données
Bn_rootD <- read.csv("data/raw-data/root_length_250318.csv", h = T, sep = ";") %>%
  filter(id_taxo == "Bra_napus") %>%
  select(leg, div, RLMax) %>%
  pivot_longer(cols = RLMax, 
               names_to = "variable", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("div", "leg")) %>%
  filter(!is.na(value))%>%
  mutate(trait = "Maximum root \nlength (cm)",
         organ = "")

# fusion des dataset des 5 traits
Bn_trait <- bind_rows(Bn_rootD, Bn_germ, Bn_ldmc, Ncont, Pcont) %>% 
  mutate(organ = str_replace(organ, "Fine.roots", "Fine roots"),
         measure = paste(organ, trait, sep="\n")) %>%
  filter(!organ %in% c("Limb", "Stem", "Taproot"))

# Graph
Bn_trait_p <- ggplot(subset(Bn_trait, organ == ""), 
                     aes(x = div, y = value))+
  geom_bar(stat = "identity", position=position_dodge(width= 2),
           aes(fill = div), colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.1, position=position_dodge(width= 2),
                colour = "black")+
  labs(x= "", y = "", fill = "Earthworm diversity")+
  facet_grid(trait~leg, scales = "free")+
  theme_bw()+
  theme(axis.text.x=element_blank())

Bn_content_p <- ggplot(subset(Bn_trait, !organ == ""), 
                     aes(x = div, y = value))+
  geom_bar(stat = "identity", position=position_dodge(width= 2),
           aes(fill = div), colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.1, position=position_dodge(width= 2),
                colour = "black")+
  labs(x= "", y = "", fill = "Earthworm diversity")+
  facet_grid(measure~leg, scales = "free")+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.y=element_blank())

Bn_p <- ggarrange(Bn_biom_p, Bn_trait_p, Bn_content_p, 
                  common.legend = T, labels = c("A", "B"), ncol = 3)

png("figures/pl/Bn.png", width = 700, height = 500)
Bn_p
dev.off()


