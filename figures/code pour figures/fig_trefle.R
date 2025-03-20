librarian::shelf(dplyr, ggplot2, Rmisc, tidyr, stringr, ggpubr)

# profondeur max d'enracinement
  ##téléchargement et préparation des données
      Ta_rootD <- read.csv("data/raw-data/root_length_250318.csv", h = T, sep = ";") %>%
        #filter(id_taxo == "Tri_alex") %>%
        select(div, RLMax) %>%
        summarySE(measurevar=c("RLMax"), na.rm = T,
                  groupvars=c("div"))
  ## graph
      Ta_rootD_p <- ggplot(Ta_rootD, aes(x = div, y = RLMax))+
        geom_bar(stat = "identity", aes(fill = div), colour="black")+
        geom_errorbar(aes(ymin=RLMax-se, ymax=RLMax+se), width=0.2)+
        labs(x= "", y = "Maximum root length (cm)", colour = "Earthworm\ndiversity")+
        lims(y=c(0,50))+
        paletteer::scale_fill_paletteer_d("beyonce::X115") +
        theme_bw()+
        theme(axis.title.x=element_blank())+ 
        theme(legend.position = "none",
              axis.title=element_text(size=9))
      
# Nodulation
    ## téléchargement et préparation des données de nodulation
      Ta_nod <- read.csv("data/derived-data/pl_Ta_nod.csv", h = T, sep = ",") %>%
        select(div, leg, strata, NOD, DW) %>%
        mutate(NOD_dm = NOD/DW) %>%
        select(-DW) %>%
        pivot_longer(cols = c(NOD, NOD_dm), 
                     names_to = "variable", values_to = "value") %>%
        summarySE(measurevar=c("value"), na.rm = T,
                  groupvars=c("leg", "div", "strata", "variable")) %>%
        mutate(strata = case_when(strata == 1 ~-3.5,
                                  strata == 2 ~-18.5,
                                  strata == 3 ~-40.5,
                                  strata == 4 ~-70.5)) %>%
        filter(!is.na(value)) 
      
      ## téléchargement et préparation des données de masse + fusion
      Ta <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
        filter(id_taxo == "Tri_alexa") %>%
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
        bind_rows(Ta_nod)
      
      ## Graphique du profil de nodulation
      Ta_biom_nod_labs <- as_labeller(c(DW = "Dry mass profile (g)",
                                        NOD = "No nodules per \nroot section",
                                        NOD_dm = "No nodules per \nroot section dry mass"))
      
      Ta_biom_nod_p <- ggplot(subset(Ta, variable %in% c("DW")),
                              aes(x = strata, y = value, colour = div))+
        geom_point(position=position_dodge(width= 2), size = 2)+
        geom_line(alpha = 0.2)+
        geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                      width=.1, position=position_dodge(width= 2))+
        geom_vline(xintercept = 0)+
        labs(x= "Soil depth (cm)", y = "", colour = "Earthworm\ndiversity")+
        coord_flip()+
        paletteer::scale_color_paletteer_d("beyonce::X115") +
        facet_grid(.~variable, scales = "free", labeller = Ta_biom_nod_labs)+
        theme_bw()+
        theme(axis.title.x=element_blank())
      
      
      
      
# Nodulation
    ## téléchargement et préparation des données de nodulation
        Ta_nod <- read.csv("data/derived-data/pl_Ta_nod.csv", h = T, sep = ",") %>%
          select(div, leg, strata, NOD, DW) %>%
          mutate(NOD_dm = NOD/DW) %>%
          select(-DW) %>%
          pivot_longer(cols = c(NOD, NOD_dm), 
                       names_to = "variable", values_to = "value") %>%
          summarySE(measurevar=c("value"), na.rm = T,
                    groupvars=c("leg", "div", "strata", "variable")) %>%
          mutate(strata = case_when(strata == 1 ~-3.5,
                                    strata == 2 ~-18.5,
                                    strata == 3 ~-40.5,
                                    strata == 4 ~-70.5)) %>%
          filter(!is.na(value)) 

    ## téléchargement et préparation des données de masse + fusion
        Ta <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
          filter(id_taxo == "Tri_alexa") %>%
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
          bind_rows(Ta_nod)

    ## Graphique du profil de nodulation
    Ta_biom_nod_labs <- as_labeller(c(DW = "Dry mass profile (g)",
                                      NOD = "No nodules per \nroot section",
                                      NOD_dm = "No nodules per \nroot section dry mass"))
    
      Ta_biom_nod_p <- ggplot(subset(Ta, variable %in% c("DW")),
                              aes(x = strata, y = value, colour = div))+
        geom_point(position=position_dodge(width= 2), size = 2)+
        geom_line(alpha = 0.2)+
        geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                      width=.1, position=position_dodge(width= 2))+
        geom_vline(xintercept = 0)+
        labs(x= "Soil depth (cm)", y = "", colour = "Earthworm\ndiversity")+
        coord_flip()+
        paletteer::scale_color_paletteer_d("beyonce::X115") +
        facet_grid(.~variable, scales = "free", labeller = Ta_biom_nod_labs)+
        theme_bw()+
        theme(axis.title.x=element_blank())
    
    ## Graphique nodule ~ mass
      Ta_rel_nod <- read.csv("data/derived-data/pl_Ta_nod.csv", h = T, sep = ",") %>%
        tibble() %>%
        filter(!is.na(DW)) %>%
        #select(RT, div, leg, strata, distance, NOD, DW) %>%
        group_by(RT, div, leg) %>%
        dplyr::summarise(NODs = sum(NOD, na.rm = T), 
                  DWs = sum(DW, na.rm = T))  
      
      Ta_rel_nod_p <- ggplot(Ta_rel_nod, aes(x = DWs, y = NODs, colour = div))+
            geom_point(size = 2)+
            geom_smooth(method = lm, se = FALSE)+
            labs(x= "Total root dry mass\nper rhizotron (g)", 
                 y = "Total nodule number\nper rhizotron", colour = "Earthworm\ndiversity")+
            theme_bw()+
        paletteer::scale_color_paletteer_d("beyonce::X115") +
            theme(legend.position = "none",
                  axis.title=element_text(size=9))

# Concentrations en nutriments
    ## téléchargement et préparation des données N
      Ncont <- read.csv("data/derived-data/pl_Ncont.csv") %>%
        filter(id_taxo == "Tri_alexa") %>%
        select(div, Fine.roots, Limb, Leaves, Stem) %>%
        pivot_longer(cols = c(Fine.roots, Limb, Leaves, Stem), 
                     names_to = "variable", values_to = "value") %>%
        summarySE(measurevar=c("value"), na.rm = T,
                  groupvars=c("div", "variable")) %>%
        filter(!is.na(value)) %>%
        mutate(element = "N content  (mg.g-1)")

    ## téléchargement et préparation des données P 
      Pcont <- read.csv("data/derived-data/pl_Pcont.csv") %>%
        filter(id_taxo == "Tri_alexa") %>%
        select(div, Fine.roots, Leaves) %>%
        pivot_longer(cols = c(Fine.roots, Leaves), 
                     names_to = "variable", values_to = "value") %>%
        mutate(value = ifelse(value<0,0.001,value)) %>%
        summarySE(measurevar=c("value"), na.rm = T,
                  groupvars=c("div", "variable")) %>%
        filter(!is.na(value)) %>%
        mutate(element = "P content (mg.g-1)")

    ## Fusion datasets N et P
      cont <- bind_rows(Ncont, Pcont) %>% 
          mutate(variable = str_replace(variable, "Fine.roots", "Fine roots"))

    ## Graph
      Ta_nutrient_content <- ggplot(subset(cont, variable %in% c("Fine roots", "Leaves")), 
             aes(x = div, y = value))+
        geom_bar(stat = "identity", position=position_dodge(width= 2),
                 aes(fill = div), colour="black")+
        geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                      width=.1, position=position_dodge(width= 2),
                      colour = "black")+
        paletteer::scale_fill_paletteer_d("beyonce::X115") +
        labs(x= "", y = "", fill = "Earthworm diversity")+
        facet_grid(element~variable, scales = "free")+
        theme_bw()+
        theme(#axis.text.x = element_blank(),
              axis.title = element_blank())

# Figure globale
  gp1 <- ggarrange(Ta_rel_nod_p, Ta_rootD_p, ncol = 1)
  ta_p <- ggarrange(Ta_biom_nod_p, gp1, Ta_nutrient_content,
            common.legend = T, labels = c("A", "B"), ncol = 3,
            widths = c(0.33, 0.25, .41))

  png("figures/pl/Ta.png", width = 700, height = 500)
  ta_p
  dev.off()
