# Source des fonctions stats et graphiques
source("analyses/1.0. fonctions stats.R")
source("analyses/1.1. fonctions graphiques.R")


# Chargement des données
myTa <- read.csv("data/derived-data/pl_Ta_nod.csv", h = T, sep = ",")


# Nodulation
  ## Nodule systeme racinaire total
  vb <- c("NOD", "DW")
  myTa_nod <- aggregate(myTa[,vb], list(RT = myTa$RT, numPL = myTa$numPL,  div = myTa$div,  leg = myTa$leg, 
                                        id_loca = myTa$id_loca, id_taxo = myTa$id_taxo, distance = myTa$distance, serie = myTa$serie, position = myTa$position, PAR = myTa$PAR), sum)
  myTa_nod_temp <- myTa_nod[!is.na(myTa_nod$NOD),]
  Ta_nod <- myfun_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NOD), seuil_prcRE = 24) ## log.transformed
  Ta_nod$select
  pTa_nod <- myplot_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NOD), mod = Ta_nod$mod10A, 
                        ylabl = "T. alexandrium nodules (total, ind-1, log transf.)", plot_PAR = F)


      ### nodules strate 1
      myTa_nod1 <- myTa[myTa$strata == 1, ]
      myTa_nod1_temp <- myTa_nod1[!is.na(myTa_nod1$NOD),]
      Ta_nod1 <- myfun_Leg(mdata = myTa_nod1_temp, vb = sqrt(myTa_nod1_temp$NOD), seuil_prcRE = 24) # 
      Ta_nod1$resum
      Ta_nod1$select
      pTa_nod1 <- myplot_Leg(mdata = myTa_nod1_temp, vb = myTa_nod1_temp$NOD, mod = Ta_nod1$mod10A, 
                             ylabl = "T. alexandrium nodules (0-10 cm, ind-1)", plot_PAR = F)
      
      ### nodules strate 2
      myTa_nod2 <- myTa[myTa$strata == 2, ]
      myTa_nod2_temp <- myTa_nod2[!is.na(myTa_nod2$NOD),]
      Ta_nod2 <- myfun_Leg(mdata = myTa_nod2_temp, vb = sqrt(myTa_nod2_temp$NOD), seuil_prcRE = 24) # 
      Ta_nod2$resum
      Ta_nod2$select
      pTa_nod2 <- myplot_Leg(mdata = myTa_nod2_temp, vb = myTa_nod2_temp$NOD, mod = Ta_nod2$mod10A, 
                             ylabl = "T. alexandrium nodules (20-30 cm, ind-1)", plot_PAR = F)
      
      ### nodules strate 3
      myTa_nod3 <- myTa[myTa$strata == 3, ]
      myTa_nod3_temp <- myTa_nod3[!is.na(myTa_nod3$NOD),]
      Ta_nod3 <- myfun_Leg(mdata = myTa_nod3_temp, vb = sqrt(myTa_nod3_temp$NOD), seuil_prcRE = 24) # 
      Ta_nod3$resum
      Ta_nod3$select
      pTa_nod3 <- myplot_Leg(mdata = myTa_nod3_temp, vb = myTa_nod3_temp$NOD, mod = Ta_nod3$mod10A, 
                             ylabl = "T. alexandrium nodules (30-50 cm, ind-1)", plot_PAR = F)
      
      ### nodules strate 4
      myTa_nod4 <- myTa[myTa$strata == 4, ]
      myTa_nod4_temp <- myTa_nod4[!is.na(myTa_nod4$NOD),]
      Ta_nod4 <- myfun_Leg(mdata = myTa_nod4_temp, vb = sqrt(myTa_nod4_temp$NOD), seuil_prcRE = 24) # 
      Ta_nod4$resum
      Ta_nod4$select
      pTa_nod4 <- myplot_Leg(mdata = myTa_nod4_temp, vb = myTa_nod4_temp$NOD, mod = Ta_nod4$mod10A, 
                             ylabl = "T. alexandrium nodules (50-100 cm, ind-1)", plot_PAR = F)


  ## Profil racinaire nodules
  decalage <- as.numeric(rep(c(1.5, 0.5, -0.5, -1.5), 4)) 
  Ta_DistriNOD <- rbind.data.frame(pTa_nod1$moyFRic, pTa_nod2$moyFRic, pTa_nod3$moyFRic, pTa_nod4$moyFRic)
  Ta_DistriNOD$prof <- as.numeric(rep(c("-5", "-20", "-40", "-70"), each = 4))+decalage
  Ta_DistriNOD$strate <- rep(c("R1", "R2", "R3", "R4"), each = 4)
  Ta_DistriNOD$semin <- Ta_DistriNOD$x - 2*Ta_DistriNOD$se
  Ta_DistriNOD$semax <- Ta_DistriNOD$x + 2*Ta_DistriNOD$se

  Ta_nod_PR <- ggplot(data = Ta_DistriNOD, aes(x = prof, y = x^2, group = div, shape = div))+
    geom_line(aes(linetype = div))+
    geom_point(size = 4)+
    geom_errorbar(aes(ymin = semin^2, ymax = semax^2), width = 0.4)+
    coord_flip()+
    scale_shape_discrete(name="Earthworm FRic")+
    xlim(c(-90, 0))+
    ylab ("T. alexandrium nodules per individual")+
    xlab("Soil depth (cm)")+
    scale_linetype_discrete(guide = FALSE)+
    theme(legend.justification = c(0,0), legend.position = c(0.5,0.05),
          panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")
  png("figures/pl/Ta_nodule_profile.png", width = 500, height = 1000)
  Ta_nod_PR
  dev.off()

      ### nodules par unit? de masse racinaire
      myTa_nod$NODms <- myTa_nod$NOD/myTa_nod$DW
      myTa_nod_temp <- myTa_nod[!is.na(myTa_nod$NODms),]
      Ta_nodms <- myfun_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NODms), seuil_prcRE = 24) ## log.transformed
      Ta_nodms$select
      pTa_nodms <- myplot_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NODms), mod = Ta_nod$mod10A, 
                              ylabl = "T. alexandrium nodules (total, mg-1, log transf.)", plot_PAR = F)
      png("figures/pl/Ta_nodule_perMass.png", width = 500, height = 500)
      pTa_nodms$pDIV
      dev.off()

      ### relation nodules - masse racinaire
      vb <- c("NOD","DW")
      Ta_RT_nod <- aggregate(myTa_nod[,vb], by=list(div= myTa_nod$div, RT=myTa_nod$RT), mean, na.rm=T)
      #Ta_se_nod <- aggregate(myTa_nod[,vb], by=list(div= myTa_nod$div), mean, na.rm=T)
      
      rel_ta <- ggplot(Ta_RT_nod, aes(x=DW, y=NOD)) + 
        geom_point(aes(color=div),size=8)+ 
        scale_color_manual(values=c("black", "#ffeda0", "#feb24c","#f03b20"),
                           name="",labels=c("Ew0","FR1","FR2","FR3"))+
        xlab("T. alexandrium root mass (mg ind-1)")+
        ylab("T. alexandrium nodules (ind-1)")+
        xlim(0,100)+
        ylim(0,100)+
        theme_classic() +
        theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
              text =  element_text(face = "plain",
                                   color = "black", 
                                   hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                   margin = margin(), debug = FALSE))+
        theme(axis.text = element_text(color="black", size=14)) +
        theme(axis.title = element_text(color="black", size=14))+
        guides(shape = FALSE)#, alpha=FALSE +
      
            # on ajoute les regressions
            reg_ta_div0 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div0",])
            summary(reg_ta_div0)
            reg_ta_div1 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div1",])
            summary(reg_ta_div1)
            reg_ta_div2 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div2",])
            summary(reg_ta_div2)
            reg_ta_div3 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div3",])
            summary(reg_ta_div3)
            
        final <- rel_ta + geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div0",],aes(x =DW, y =NOD ), 
                           method = "lm", alpha = 0.1, color="black", se=F)+
        annotate(geom="text", x=77, y=45, label="r2=0.84", size=5) +
        geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div1",],aes(x =DW, y =NOD ), 
                    method = "lm", alpha = 0.1, color="#ffeda0", se=F)+
        annotate(geom="text", x=75, y=75, label="r2=0.79", size=5, color="#ffeda0") +
        geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div2",],aes(x =DW, y =NOD ), 
                    method = "lm", alpha = 0.1, color="#feb24c", se=F)+
        annotate(geom="text", x=90, y=90, label="r2=0.89", size=5,color="#feb24c") +
        geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div3",],aes(x =DW, y =NOD ), 
                    method = "lm", alpha = 0.1, color="#f03b20", se=F)+
        annotate(geom="text", x=75, y=60, label="r2=0.89", size=5,color="#f03b20")
        
        png("figures/pl/Ta_nodule_perMass_lm.png", width = 500, height = 500)
        final
        dev.off()


# Biomasse 
    ## Biomasse seche aerienne
    SR_ratio <- mybiom_Ind %>%
          filter(id_taxo == "Tri_alexa") %>%
          mutate(shootroot = case_when(strata == 0 ~ "shoot",
                                       strata != 0 ~ "root")) %>%
          pivot_wider(id_cols = c(RT, div, leg, id_loca, id_taxo, distance, serie, position, PAR, numPL),
                      values_from = DW, names_from = shootroot, values_fn = sum) %>%
          mutate(shootroot = shoot/root)
    Ta_A_DW <- myfun_Leg(mdata = SR_ratio, vb = log(SR_ratio$shoot), seuil_prcRE = 30) # transf log
    Ta_A_DW$resum
    Ta_A_DW$select
    pTa_A_DW <- myplot_Leg(mdata = SR_ratio, vb = SR_ratio$shoot, mod = Ta_A_DW$mod10A, 
                           ylabl = "T. alexandrium \nshoot mass (DW, mg)", plot_PAR = T)
    png("figures/pl/mass_Ta_shoot.png", width = 1000, height = 480)
    pTa_A_DW$pDIV
    dev.off()
    
    ## Biomasse seche racine totale
    Ta_R_DW <- myfun_Leg(mdata = SR_ratio, vb = log1p(SR_ratio$root), seuil_prcRE = 30) ## log.transformed
    pTa_R_DW <- myplot_Leg(mdata = SR_ratio, vb = log1p(SR_ratio$root), mod = Ta_R_DW$mod10A, 
                           ylabl = "T. alexandrium \nroot mass (total, DW mg, log transf.)", plot_PAR = T)
    png("figures/pl/mass_Ta_root.png", width = 1000, height = 480)
    Ta_R_DW$pDIV
    dev.off()    
    
        ### Biomasse s?che strate 1            
        Ta_R1 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 1, ]
        Ta_R1_DW <- myfun_Leg(mdata = Ta_R1, vb = sqrt(Ta_R1$DW), seuil_prcRE = 30) # idem
        Ta_R1_DW$resum
        Ta_R1_DW$select
        pTa_R1_DW <- myplot_Leg(mdata = Ta_R1, vb = sqrt(Ta_R1$DW), mod = Ta_R1_DW$mod10A, 
                                ylabl = "T. alexandrium \nroot mass (0-10 cm, DW mg, sqrt transf.)", plot_PAR = T)
        
        ### Biomasse s?che strate 2            
        Ta_R2 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 2, ]
        Ta_R2_DW <- myfun_Leg(mdata = Ta_R2, vb = sqrt(Ta_R2$DW), seuil_prcRE = 30) # exemple avec transf. sqrt()
        Ta_R2_DW$resum
        pTa_R2_DW <- myplot_Leg(mdata = Ta_R2, vb = sqrt(Ta_R2$DW), mod = Ta_R2_DW$mod10A, 
                                ylabl = "T. alexandrium \nroot mass (10-30 cm, DW mg, sqrt transf.)", plot_PAR = T)
        
        ### Biomasse s?che strate 3            
        Ta_R3 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 3, ]
        Ta_R3_DW <- myfun_Leg(mdata = Ta_R3, vb = sqrt(Ta_R3$DW), seuil_prcRE = 30) # exemple avec transf. sqrt()
        Ta_R3_DW$resum
        pTa_R3_DW <- myplot_Leg(mdata = Ta_R3, vb = sqrt(Ta_R3$DW), mod = Ta_R3_DW$mod10A, 
                                ylabl = "T. alexandrium \nroot mass (30-50 cm, DW mg, sqrt transf.)", plot_PAR = T)
        
        ### Biomasse s?che strate 4            
        Ta_R4 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 4, ]
        Ta_R4_DW <- myfun_Leg(mdata = Ta_R4, vb = Ta_R4$DW, seuil_prcRE = 30) 
        # n'arrive pas ? cr?er les mod?les car design hyper d?s?quilibr? (cf summary(Ta_R4$div))
        pTa_R4_DW <- myplot_Leg(mdata = Ta_R4, vb = sqrt(Ta_R4$DW), mod = Ta_R4_DW$mod10A, 
                                ylabl = "T. alexandrium \nroot mass (50-90 cm, DW mg, sqrt transf.)", plot_PAR = T)
        
    
    ## Profil racinaire      
    Ta_DistriDW <- rbind.data.frame(pTa_A_DW$moyFRic, pTa_R1_DW$moyFRic, pTa_R2_DW$moyFRic, pTa_R3_DW$moyFRic, pTa_R4_DW$moyFRic)
    Ta_DistriDW$prof <- as.numeric(rep(c("10", "-5", "-20", "-40", "-70"), each = 4))
    Ta_DistriDW$strate <- rep(c("A", "R1", "R2", "R3", "R4"), each = 4)
    Ta_DistriDW$semin <- Ta_DistriDW$x - 2*Ta_DistriDW$se
    Ta_DistriDW$semax <- Ta_DistriDW$x + 2*Ta_DistriDW$se
    Ta_PR <- ggplot(data = Ta_DistriDW[Ta_DistriDW$strate != "A", ], aes(x = prof, y = x^2, group = div, shape = div))+
      geom_line(aes(linetype = div))+
      geom_point(size = 4)+
      geom_errorbar(aes(ymin = semin^2, ymax = semax^2), width = 0.4)+
      coord_flip()+
      scale_shape_discrete(name="Earthworm FRic")+
      xlim(c(-90, 0))+
      ylim(c(0, 20))+
      ylab ("T. alexandrium \nroot profile (DW, mg)")+
      xlab("Soil depth\n(cm)")+
      scale_linetype_discrete(guide = FALSE)+
      # scale_shape_discrete(guide = FALSE)+
      theme(legend.justification = c(0,0), legend.position = c(0.5,0.05),
            panel.background = element_rect(fill="white", colour = "black"),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")
    
    png("figures/pl/mass_Ta_rootProfile.png", width = 500, height = 1000)
    Ta_PR
    dev.off()
    
    # Shoot/root ratio
    Ta_SR_DW <- myfun_Leg(mdata = SR_ratio, vb = sqrt(SR_ratio$shootroot), seuil_prcRE = 30)  # aucune transf. ne permet de normaliser les r?sidus
    Ta_SR_DW$resum
    pTa_SR_DW <- myplot_Leg(mdata = SR_ratio, vb = sqrt(SR_ratio$shootroot), mod = Ta_SR_DW$mod9A, 
                            ylabl = "T. alexandrium\nShoot:Root ratio (DW)", plot_PAR = T)
    png("figures/pl/mass_Ta_shootroot.png", width = 1000, height = 480)
    pTa_SR_DW$pDIV
    dev.off()

# Exemple de tableau r?capitulatif
Ta_recap <- matrix(NA, ncol = 16, nrow = 10, 
                   dimnames = list(c("Shoot", "Root S1", "Root S2", "Root S3", "Root S4", "Ta_nod", "Ta_nod1", "Ta_nod2", "Ta_nod3", "Ta_nod4"), 
                                   c(rownames(Ta_A_DW$select)[6:17], "Div0 (mean)", "Div1 (mean +/- se)", "Div2 (mean +/- se)", "Div3 (mean +/- se)")))
Ta_recap[1, 1:12] <- Ta_A_DW$select[6:17,1]
Ta_recap[2, 1:12] <- Ta_R1_DW$select[6:17,1]
Ta_recap[3, 1:12] <- Ta_R2_DW$select[6:17,1]
Ta_recap[4, 1:12] <- Ta_R3_DW$select[6:17,1]
Ta_recap[5, 1:12] <- Ta_R4_DW$select[6:17,1]
Ta_recap[6, 1:12] <- Ta_nod$select[6:17,1]
Ta_recap[7, 1:12] <- Ta_nod1$select[6:17,1]
Ta_recap[8, 1:12] <- Ta_nod2$select[6:17,1]
Ta_recap[9, 1:12] <- Ta_nod3$select[6:17,1]
Ta_recap[10, 1:12] <- Ta_nod4$select[6:17,1]
Ta_recap[1, 13:16] <- paste(round(pTa_A_DW$moyFRic$x, 2), "+/-", 2*round(pTa_A_DW$moyFRic$se, 2)) 
Ta_recap[2, 13:16] <- paste(round(pTa_R1_DW$moyFRic$x, 2), "+/-", 2*round(pTa_R1_DW$moyFRic$se, 2))              
Ta_recap[3, 13:16] <- paste(round(pTa_R2_DW$moyFRic$x, 2), "+/-", 2*round(pTa_R2_DW$moyFRic$se, 2))              
Ta_recap[4, 13:16] <- paste(round(pTa_R3_DW$moyFRic$x, 2), "+/-", 2*round(pTa_R3_DW$moyFRic$se, 2))
Ta_recap[5, 13:16] <- paste(round(pTa_R4_DW$moyFRic$x, 2), "+/-", 2*round(pTa_R4_DW$moyFRic$se, 2))              
Ta_recap[6, 13:16] <- paste(round(pTa_nod$moyFRic$x, 2), "+/-", 2*round(pTa_nod$moyFRic$se, 2))              
Ta_recap[7, 13:16] <- paste(round(pTa_nod1$moyFRic$x, 2), "+/-", 2*round(pTa_nod1$moyFRic$se, 2))              
Ta_recap[8, 13:16] <- paste(round(pTa_nod2$moyFRic$x, 2), "+/-", 2*round(pTa_nod2$moyFRic$se, 2))              
Ta_recap[9, 13:16] <- paste(round(pTa_nod3$moyFRic$x, 2), "+/-", 2*round(pTa_nod3$moyFRic$se, 2))              
Ta_recap[10, 13:16] <- paste(round(pTa_nod4$moyFRic$x, 2), "+/-", 2*round(pTa_nod4$moyFRic$se, 2))
Ta_recap <- as.data.frame(t(Ta_recap))
write.csv(Ta_recap, "Ta_recap.csv")
