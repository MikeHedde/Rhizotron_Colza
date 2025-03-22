librarian::shelf(dplyr, ggplot2, Rmisc, tidyr, stringr, ggpubr)

# Densité de masse racinaire
## téléchargement et préparation des données de masse 
Lp <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
  filter(id_taxo == "Lol_peren") %>%
  select(div, leg, strata, DW) %>% 
  mutate(vol = case_when(strata == 0 ~ 1,
                         strata == 1 ~ 600,
                         strata == 2 ~ 1200,
                         strata == 3 ~ 1200,
                         strata == 4 ~ 3000)) %>%
  mutate(strata = case_when(strata == 1 ~-5,
                            strata == 2 ~-20,
                            strata == 3 ~-40,
                            strata == 4 ~-75))

Lp_DWvol <- Lp %>%
  mutate(DWvol = DW*vol/1000) %>%
  filter(!is.na(strata)) %>%
  summarySE(measurevar=c("DWvol"), na.rm = T,
            groupvars=c("leg", "div", "strata"))

## Graphique du profil de densité de masse 
Lp_biom_p <- ggplot(Lp_DWvol,
                    aes(x = strata, y = DWvol, colour = div))+
  geom_point(position=position_dodge(width= 6), size = 2)+
  geom_errorbar(aes(ymin=DWvol-se, ymax=DWvol+se), 
                width=.1, position=position_dodge(width= 6))+
  scale_x_continuous(breaks=c(0, -10, -30, -50, -100), 
                     limits = c(-100,0))+
  geom_vline(xintercept = c(0, -10, -30, -50), color = "gray")+
  labs(x= "Soil depth (cm)", y = "Dry mass density \n (g cm-3)", 
       colour = "Earthworm\ndiversity")+
  facet_grid(.~leg)+
  coord_flip()+
  paletteer::scale_color_paletteer_d("beyonce::X115") +
  theme_bw()+
  theme(legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  guides(colour=guide_legend(nrow = 1,byrow=TRUE))
Lp_biom_p

# Masse aérienne
## préparation des données de masse 
Lp_aerial <- Lp %>%
  filter(vol == 1) %>%
  select(leg, div, DW) %>%
  pivot_longer(cols = DW, 
               names_to = "variable", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("div", "leg")) %>%
  filter(!is.na(value))%>%
  mutate(trait = "Aerial dry \nmass (g)")

# profondeur max d'enracinement
##téléchargement et préparation des données
Lp_rootD <- read.csv("data/raw-data/root_length_250318.csv", h = T, sep = ";") %>%
  filter(id_taxo == "Lol_peren") %>%
  select(leg, div, RLMax) %>%
  pivot_longer(cols = RLMax, 
               names_to = "variable", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("div", "leg")) %>%
  filter(!is.na(value))%>%
  mutate(trait = "Maximum root \nlength (cm)")

# téléchargement et préparation des données de Leaf dry matter content 
Lp_ldmc <- read.csv("data/raw-data/LDMC_v3.csv", h = T, sep = ";") %>%
  filter(id_taxo == "Lol_peren") %>%
  select(div, leg, LDMC) %>%
  pivot_longer(cols = LDMC, 
               names_to = "variable", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("div", "leg")) %>%
  filter(!is.na(value))%>%
  mutate(trait = "Leaf dry matter \ncontent (mg.g-1)")

# fusion des dataset des 4 traits
Lp_trait <- bind_rows(Lp_aerial, Lp_rootD, Lp_ldmc)

# Graph
Lp_trait_p <- ggplot(Lp_trait, aes(x = div, y = value))+
  geom_bar(stat = "identity", position=position_dodge(width= 2),
           aes(fill = div), colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.1, position=position_dodge(width= 2),
                colour = "black")+
  labs(x= "", y = "", fill = "Earthworm diversity")+
  paletteer::scale_fill_paletteer_d("beyonce::X115") +
  facet_grid(trait~leg, scales = "free")+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none",
        axis.title=element_text(size=9))
Lp_trait_p

# téléchargement et préparation des données sur les concentrations en N
Ncont <- read.csv("data/derived-data/pl_Ncont.csv") %>%
  filter(id_taxo == "Lol_peren") %>%
  select(div, leg, Fine.roots, Limb, Leaves, Stem, Taproot) %>%
  pivot_longer(cols = c(Fine.roots, Limb, Leaves, Stem, Taproot), 
               names_to = "organ", values_to = "value") %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("leg", "div", "organ")) %>%
  filter(!is.na(value)) %>%
  mutate(trait = "N content  (mg.g-1)",
         element = "N")

# téléchargement et préparation des données sur les concentrations en P
Pcont <- read.csv("data/derived-data/pl_Pcont.csv") %>%
  filter(id_taxo == "Lol_peren") %>%
  select(div, leg, Fine.roots, Limb, Leaves, Stem, Taproot) %>%
  pivot_longer(cols = c(Fine.roots, Leaves), 
               names_to = "organ", values_to = "value") %>%
  mutate(value = ifelse(value<0,0.001,value)) %>%
  summarySE(measurevar=c("value"), na.rm = T,
            groupvars=c("leg", "div", "organ")) %>%
  filter(!is.na(value)) %>%
  mutate(trait = "P content (mg.g-1)",
         element = "P")

## Fusion datasets N et P
cont <- bind_rows(Ncont, Pcont) %>% 
  mutate(organ = str_replace(organ, "Fine.roots", "Fine roots")) %>%
  filter(organ %in% c("Fine roots", "Leaves")) %>%
  mutate(measure = paste(organ, trait, sep = "\n"))

## Graph
Lp_Ncont_p <- ggplot(subset(cont, element == "N"),
                     aes(x = div, y = value))+
  geom_bar(stat = "identity", position=position_dodge(width= 2),
           aes(fill = div), colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.1, position=position_dodge(width= 2),
                colour = "black")+
  labs(x= "", y = "", fill = "Earthworm diversity")+
  facet_grid(measure~leg, scales = "free")+
  paletteer::scale_fill_paletteer_d("beyonce::X115") +
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())
Lp_Ncont_p

Lp_Pcont_p <- ggplot(subset(cont, element == "P"),
                     aes(x = div, y = value))+
  geom_bar(stat = "identity", position=position_dodge(width= 2),
           aes(fill = div), colour="black")+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.1, position=position_dodge(width= 2),
                colour = "black")+
  labs(x= "", y = "", fill = "Earthworm diversity")+
  facet_grid(measure~leg, scales = "free")+
  paletteer::scale_fill_paletteer_d("beyonce::X115") +
  theme_bw()+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())
Lp_Pcont_p

# Figure globale
Lp_content_p <- ggarrange(Lp_Ncont_p, Lp_Pcont_p, labels = c("C", "D"), ncol = 1)

Lp_p <- ggarrange(Lp_biom_p, Lp_trait_p, Lp_content_p,
                  common.legend = T, labels = c("A", "B", ""), ncol = 3,
                  widths = c(0.3, .3, 0.3))
Lp_p

png("figures/pl/Lp.png", width = 700, height = 500)
Lp_p
dev.off()
