librarian::shelf(dplyr, ggplot2, Rmisc, tidyr, 
                 stringr, ggpubr, paletteer, scatterpie)


# masse
    ## téléchargement et préparation des données
    biom <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
      select(RT, div, leg, strata, DW) %>% 
      group_by(RT, div, leg, strata) %>%
      dplyr::summarise(DWs= sum(DW)) %>%
      pivot_longer(cols = c(DWs), 
                   names_to = "variable", values_to = "value") %>%
      summarySE(measurevar=c("value"), na.rm = T,
                groupvars=c("leg", "div", "strata", "variable")) %>%
      mutate(strata = case_when(strata == 0 ~ 10,
                                strata == 1 ~-3.5,
                                strata == 2 ~-18.5,
                                strata == 3 ~-40.5,
                                strata == 4 ~-70.5)) %>%
      mutate(variable = str_replace(variable, "DW", "Dry mass (g)"))
    
    ## représentation graphique
    biom_p <- ggplot(biom, aes(x = strata, y = value, colour = div))+
      geom_point(position=position_dodge(width= 2), size = 2)+
      geom_line(alpha = 0.2)+
      geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                    width=.1, position=position_dodge(width= 2)) + 
      scale_y_continuous(trans='log10')+
      paletteer::scale_color_paletteer_d("beyonce::X115") +
      geom_vline(xintercept = 0)+
      xlim(-80,20)+
      labs(x= "Soil depth (cm)", y = "Mass (g; log scale)", colour = "Earthworm\ndiversity")+
      coord_flip()+
      facet_grid(~leg, scales = "free")+
      theme_bw()+
      theme(legend.position="left")+
      guides(colour=guide_legend(ncol = 1,bycol=TRUE))

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
        geom_hline(yintercept = 0)+
        geom_scatterpie(data = biom_tax, 
                        aes(x = pos, y = strata, fill = id_taxo, r = log(all)/1.3), 
                        cols = c("Bra_napus", "Ger_disse", "Lol_peren", "Tri_alexa", "Ver_persi"), 
                        colour = "white")+
        labs(fill="Plant species")+ 
        ylim(-80,20)+
        facet_grid(.~leg)+
        paletteer::scale_fill_paletteer_d("calecopal::lupinus") +
        theme_bw()+
        theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)),
              axis.text = element_blank(),
              legend.position="right",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())+
        guides(fill=guide_legend(ncol=1,bycol=TRUE))
      
      
biom_tot_p <- ggarrange(biom_p, biom_tax_p, widths = c(1, .8), align = "h",
                        ncol = 2, labels = c("A", "B"))
biom_tot_p

png("figures/pl/biom_tot.png", width = 700, height = 500)
biom_tot_p
dev.off()


