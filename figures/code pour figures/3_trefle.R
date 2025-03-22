librarian::shelf(dplyr, ggplot2, Rmisc, tidyr, stringr, ggpubr)

# Densité de masse racinaire
      ## téléchargement et préparation des données de masse 
      Ta <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
        filter(id_taxo == "Tri_alexa") %>%
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
      
        Ta_DWvol <- Ta %>%
        mutate(DWvol = DW*vol/1000) %>%
          filter(!is.na(strata)) %>%
        summarySE(measurevar=c("DWvol"), na.rm = T,
                  groupvars=c("leg", "div", "strata"))
      
      ## Graphique du profil de densité de masse 
      Ta_biom_p <- ggplot(Ta_DWvol,
                          aes(x = strata, y = DWvol, colour = div))+
        geom_point(position=position_dodge(width= 6), size = 2)+
        geom_errorbar(aes(ymin=DWvol-se, ymax=DWvol+se), 
                      width=.1, position=position_dodge(width= 6))+
        scale_x_continuous(breaks=c(0, -10, -30, -50, -100), 
                           limits = c(-100,0))+
        geom_vline(xintercept = c(0, -10, -30, -50), color = "gray")+
        labs(x= "Soil depth (cm)", y = "Dry mass density \n (g cm-3)", 
             colour = "Earthworm\ndiversity")+
        coord_flip()+
        paletteer::scale_color_paletteer_d("beyonce::X115") +
        theme_bw()+
        theme(legend.position="none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())+
        guides(colour=guide_legend(nrow = 1,byrow=TRUE))
      Ta_biom_p
      
# Masse aérienne
      ## préparation des données de masse 
      Ta_aerial <- Ta %>%
        filter(vol == 1) %>%
        select(div, DW) %>%
        summarySE(measurevar=c("DW"), na.rm = T,
                  groupvars=c("div"))
        
      Ta_aerial_p <- ggplot(Ta_aerial, aes(x = div, y = DW))+
        geom_bar(stat = "identity", aes(fill = div), colour="black")+
        geom_errorbar(aes(ymin=DW-se, ymax=DW+se), width=0.2)+
        labs(x= "", y = "Aerial dry \nmass (g)", colour = "Earthworm\ndiversity")+
        lims(y=c(0,50))+
        paletteer::scale_fill_paletteer_d("beyonce::X115") +
        theme_bw()+
        theme(axis.title.x=element_blank(),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title=element_text(size=9))
      Ta_aerial_p
      
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
        labs(x= "", y = "Maximum root \nlength (cm)", colour = "Earthworm\ndiversity")+
        lims(y=c(0,50))+
        paletteer::scale_fill_paletteer_d("beyonce::X115") +
        theme_bw()+
        theme(axis.title.x=element_blank())+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.title=element_text(size=9))
      Ta_rootD_p
      
# Nodulation
    ## téléchargement et préparation des données de nodulation
        Ta_nod <- read.csv("data/derived-data/pl_Ta_nod.csv", h = T, sep = ",") %>%
          select(div, leg, strata, NOD, DW) %>% 
          mutate(vol = case_when(strata == 0 ~ 1,
                                 strata == 1 ~ 600,
                                 strata == 2 ~ 1200,
                                 strata == 3 ~ 1200,
                                 strata == 4 ~ 3000)) %>%
          mutate(strata = case_when(strata == 1 ~-5,
                                    strata == 2 ~-20,
                                    strata == 3 ~-40,
                                    strata == 4 ~-75)) %>%
          mutate(NOD_dm_vol = NOD/DW/vol*1000) %>%
          select(-DW) %>%
          summarySE(measurevar=c("NOD_dm_vol"), na.rm = T,
                    groupvars=c("leg", "div", "strata"))

    ## Graphique du profil de nodulation
       Ta_nod_p <- ggplot(Ta_nod, aes(x = strata, y = NOD_dm_vol, colour = div))+
         geom_errorbar(aes(ymin=NOD_dm_vol-se, ymax=NOD_dm_vol+se), 
                       width=.1, position=position_dodge(width= 6))+
         scale_x_continuous(breaks=c(0, -10, -30, -50, -100), 
                            limits = c(-100,0))+
         geom_vline(xintercept = c(0, -10, -30, -50), color = "gray")+
         geom_point(position=position_dodge(width= 6), size = 2)+
        labs(x= "", y = "Nodules \n(N mg-1 root cm-3)", colour = "Earthworm\ndiversity")+
        coord_flip()+
        paletteer::scale_color_paletteer_d("beyonce::X115") +
         theme_bw()+
         theme(legend.position="none",
               axis.text.y = element_blank(),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank())+
         guides(colour=guide_legend(nrow = 1,byrow=TRUE))
       Ta_nod_p
       
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
                 y = "Total nodule \nnumber per rhizotron", colour = "Earthworm\ndiversity")+
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
        theme(axis.text.x = element_blank(),
              axis.title = element_blank())

# Figure globale
  gp1 <- ggarrange(Ta_aerial_p, Ta_rel_nod_p, Ta_rootD_p, ncol = 1,
                   labels = c("C", "D", "E"))
  ta_p <- ggarrange(Ta_biom_p, Ta_nod_p, gp1, Ta_nutrient_content,
            common.legend = T, labels = c("A", "B", "", "F"), ncol = 4,
            widths = c(0.2, .165, 0.25, .41))
  ta_p
  
  png("figures/pl/Ta.png", width = 700, height = 500)
  ta_p
  dev.off()
