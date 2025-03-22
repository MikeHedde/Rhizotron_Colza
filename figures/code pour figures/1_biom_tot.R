librarian::shelf(dplyr, ggplot2, Rmisc, tidyr, 
                 stringr, ggpubr, paletteer, scatterpie)


# masse
    ## téléchargement et préparation des données
    biom <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
      select(RT, div, leg, strata, DW) %>% 
      group_by(RT, div, leg, strata) %>%
      dplyr::summarise(DWs= sum(DW)) %>%
      mutate(volume = case_when(strata == "0" ~ 1,
                                strata == "1" ~ 600,
                                strata == "2" ~ 1200,
                                strata == "3" ~ 1200,
                                strata == "4" ~ 3000)) %>%
      mutate(DWvol = DWs/volume)%>%
      pivot_longer(cols = c(DWvol), 
                   names_to = "variable", values_to = "value") %>%
      summarySE(measurevar=c("value"), na.rm = T,
                groupvars=c("leg", "div", "strata", "variable")) %>%
      mutate(strata = case_when(strata == 0 ~ 10,
                                strata == 1 ~ -5,
                                strata == 2 ~ -20,
                                strata == 3 ~ -40,
                                strata == 4 ~ -75)) %>%
      mutate(variable = str_replace(variable, "DWvol", "Dry mass density (g.cm3)"))
    
    ## représentation graphique
    biom_p <- ggplot(subset(biom, strata < 0), aes(x = strata, y = value, colour = div))+
      scale_x_continuous(breaks=c(0, -10, -30, -50, -100), 
                         limits = c(-100,0))+
      geom_vline(xintercept = c(0, -10, -30, -50), color = "gray")+
      geom_point(position=position_dodge(width= 6), size = 2)+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                    width=.1, position=position_dodge(width= 6)) + 
      scale_y_continuous(trans='log10')+
      paletteer::scale_color_paletteer_d("beyonce::X115") +
      labs(x= "Soil depth (cm)", y = "Dry mass density \n(g cm-3; log scale)", colour = "Earthworm\ndiversity")+
      coord_flip()+
      facet_grid(~leg, scales = "free")+
      theme_bw()+
      theme(legend.position="none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank())+
      guides(colour=guide_legend(ncol = 1,bycol=TRUE))
    biom_p

# masse par esp
    # téléchargement et préparation des données 
      biom_tax <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
        select(RT, leg, strata, id_taxo, DW) %>% 
        group_by(RT, leg, strata, id_taxo) %>%
        dplyr::summarise(DWs= sum(DW)) %>%
        group_by(leg, strata, id_taxo) %>%
        dplyr::summarise(DWm = mean(DWs)) %>%
        mutate(strata = case_when(strata == 0 ~ 10,
                                  strata == 1 ~-3.5,
                                  strata == 2 ~-18.5,
                                  strata == 3 ~-40.5,
                                  strata == 4 ~-70.5),
               leg_po = case_when(leg == "Leg-" ~ 1,
                                  leg == "Leg+" ~ 20)) %>%
        pivot_wider(id_cols = c(leg, leg_po, strata), names_from = id_taxo, 
                    values_from = DWm, values_fn = sum, values_fill = 0) %>%
        mutate(all=Bra_napus+Ger_disse+Lol_peren+Tri_alexa+Ver_persi,
               pos = 0)
      
      
      ## Graph
      sp_lab <- c(Bra_napus = "B. napus",
                           Ger_disse = "G. dissectum",
                           Lol_peren = "L. perenne",
                           Tri_alexa = "T. alexandrinum",
                           Ver_persi = "V. persica")
      
      
      biom_tax_p <- ggplot() +
        scale_y_continuous(breaks=c(0, -10, -30, -50, -100), 
                           limits = c(-100,0))+
        geom_hline(yintercept = c(0, -10, -30, -50), color = "gray")+
        geom_scatterpie(data = subset(biom_tax, strata < 0), 
                        aes(x = pos, y = strata, fill = id_taxo, r = log(all)*.5), 
                        cols = c("Bra_napus", "Ger_disse", "Lol_peren", "Tri_alexa", "Ver_persi"), 
                        colour = "white")+
        labs(fill="Plant species", x = "Dry mass \nproportion")+ 
        facet_grid(.~leg)+
        paletteer::scale_fill_paletteer_d("calecopal::lupinus") +
        theme_bw()+
        theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)),
              axis.text.y = element_blank(),
              legend.position="right",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())+
        guides(fill=guide_legend(ncol=1,bycol=TRUE))
      biom_tax_p

      
# Masse totale et R/S ratio
      content <- read.csv("data/raw-data/P_180517.csv", h = T, sep = ";") %>%
        mutate(Nmass = Ncont*DM/1000,
               ABG = case_when(org %in% c("Leaves", "Limb", "Stem", "Aerial") ~ "Aerial",
                               org %in% c("Fine roots", "Taproot") ~ "Root")) %>%
        group_by(RT, leg, div, ABG) %>%
        dplyr::summarise(Nqt = sum(Nmass, na.rm = T), 
                         DMtot = sum(DM, na.rm = T))
      Ncont <- content %>%
        summarySE(measurevar=c("Nqt"), na.rm = T,
                  groupvars=c("leg", "div", "ABG")) %>%
        mutate(variable = "Total plant N mass") %>%
        dplyr::rename("mean" = "Nqt")
      
      DM <- content %>%
        summarySE(measurevar=c("DMtot"), na.rm = T,
                  groupvars=c("leg", "div", "ABG")) %>%
        mutate(variable = "Total plant mass") %>%
        dplyr::rename("mean" = "DMtot")
      
      qt <- bind_rows(Ncont, DM) %>%
        mutate_at(c("mean", "sd", "se", "ci"), ~ ifelse(ABG == "Aerial", .x , .x*-1))
        
      qt_p <- ggplot(qt, aes(x = div, y = mean, fill = div))+
        geom_bar(stat = "identity", position = "stack")+
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)+
        paletteer::scale_fill_paletteer_d("beyonce::X115") +
        facet_grid(variable~leg, scales = "free")+
        labs(x="Total quantity\n(g)", y="")+
        geom_hline(yintercept = 0, colour = "white")+
        theme_bw()+
        #theme(axis.title.x=element_blank())+ 
        theme(legend.position = "none",
              axis.title.x=element_text(size=12))
        
      

biom_tot_p <- ggarrange(qt_p, labels = "A", widths = c(1.5,2), common.legend = T,
                        ggarrange(biom_p, biom_tax_p, widths = c(1.5, 1.1), align = "h",
                                  ncol = 2, labels = c("B", "C"))
                        )
biom_tot_p

png("figures/pl/biom_tot.png", width = 700, height = 500)
biom_tot_p
dev.off()


