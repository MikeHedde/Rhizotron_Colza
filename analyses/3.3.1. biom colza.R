#############################################        
###           Pour les colzas             ###
#############################################

# Chargement des données
  ## Masses individuelles
  mybiom_Ind <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",")
  ## Date de germination
  myBnGerm <- read.csv("data/derived-data/pl_Bn_germ.csv", h = T, sep = ",")

# Source des fonctions stats et graphiques
source("analyses/1.0. fonctions stats.R")
source("analyses/1.1. fonctions graphiques.R")

# Date de germination
myBnGerm_m <- myfun_Ind(mdata = myBnGerm, vb = myBnGerm$SGDm, seuil_prcRE = 30)


# Biomasses
SR_ratio <- mybiom_Ind %>%
  filter(id_taxo == "Bra_napus") %>%
  mutate(shootroot = case_when(strata == 0 ~ "shoot",
                               strata != 0 ~ "root")) %>%
  pivot_wider(id_cols = c(RT, div, leg, id_loca, id_taxo, distance, serie, position, PAR, numPL),
              values_from = DW, names_from = shootroot, values_fn = sum) %>%
  mutate(shootroot = shoot/root)

  
  ## Biomasse seche aerienne
  Bn_A_all_DW <- myfun_Ind(mdata = SR_ratio, vb = sqrt(SR_ratio$shoot), seuil_prcRE = 30) # transf. sqrt
  Bn_A_all_DW$select    # ni la diversit? des VdT (ou leur pr?sence) ni la pr?sence de l?gumineuse ne change pas la biomasse a?rienne du colza
  pBn_A_all_DW <- myplot_Ind(mdata = Bn_A_all, vb = sqrt(Bn_A_all$DW), mod = Bn_A_all_DW$mod10A, 
                             ylabl = "B. napus\nshoot mass (DW, mg)", plot_PAR = F)
  png("figures/pl/mass_Bn_shoot.png", width = 1000, height = 480)
  grid.arrange(pBn_A_all_DW$pDIV, pBn_A_all_DW$pLEG, pBn_A_all_DW$pDIV_INTER,
               ncol = 3)
  dev.off()

  ## Biomasse seche racinaire totale
  Bn_R_DW <- myfun_Ind(mdata = SR_ratio, vb = sqrt(SR_ratio$root), seuil_prcRE = 30) ## sqrt transformed
  pBn_R_DW <- myplot_Ind(mdata = Bn_R, vb = sqrt(Bn_R$DW), mod = Bn_R_DW$mod10A, 
                         ylabl = "B. napus \nroot mass (total, mg, sqrt transf.)", plot_PAR = T)
  png("figures/pl/mass_Bn_root.png", width = 1000, height = 480)
  grid.arrange(pBn_R_DW$pDIV, pBn_R_DW$pLEG, pBn_R_DW$pDIV_INTER,
               ncol = 3)
  dev.off()
  
      ### Biomasse seche strate 1   
      Bn_R1_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 1, ]
      Bn_R1_all_DW <- myfun_Ind(mdata = Bn_R1_all, vb = sqrt(Bn_R1_all$DW), seuil_prcRE = 30)
      Bn_R1_all_DW$select
      pBn_R1_all_DW <- myplot_Ind(mdata = Bn_R1_all, vb = sqrt(Bn_R1_all$DW), mod = Bn_R1_all_DW$mod10A, 
                                  ylabl = "B. napus\nroot mass (0-10 cm, DW mg)", plot_PAR = F)
      
      ### Biomasse seche strate 2  
      Bn_R2_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 2, ]
      Bn_R2_all_DW <- myfun_Ind(mdata = Bn_R2_all, vb = sqrt(Bn_R2_all$DW), seuil_prcRE = 30)
      Bn_R2_all_DW$select
      pBn_R2_all_DW <- myplot_Ind(mdata = Bn_R2_all, vb = sqrt(Bn_R2_all$DW), mod = Bn_R2_all_DW$mod10A, 
                                  ylabl = "B. napus\nroot mass (10-30 cm, DW mg)", plot_PAR = F)
      
      
      ### Biomasse seche strate 3  
      Bn_R3_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 3, ]
      Bn_R3_all_DW <- myfun_Ind(mdata = Bn_R3_all, vb = sqrt(Bn_R3_all$DW), seuil_prcRE = 30)
      Bn_R3_all_DW$select
      pBn_R3_all_DW <- myplot_Ind(mdata = Bn_R3_all, vb = sqrt(Bn_R3_all$DW), mod = Bn_R3_all_DW$mod10A, 
                                  ylabl = "B. napus\nroot mass (30-50 cm, DW mg)", plot_PAR = T)
      
      ### Biomasse seche strate 4      
      Bn_R4_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 4, ]
      Bn_R4_all_DW <- myfun_Ind(mdata = Bn_R4_all, vb = sqrt(Bn_R4_all$DW), seuil_prcRE = 30)
      Bn_R4_all_DW$select
      pBn_R4_all_DW <- myplot_Ind(mdata = Bn_R4_all, vb = sqrt(Bn_R4_all$DW), mod = Bn_R4_all_DW$mod10A, 
                                  ylabl = "B. napus\nroot mass (50-90 cm, DW mg)", plot_PAR = T)
                                    
      ### Profil racinaire      
      Bn_DistriDW <- rbind.data.frame(pBn_A_all_DW$moyFRic, pBn_R1_all_DW$moyFRic, pBn_R2_all_DW$moyFRic, pBn_R3_all_DW$moyFRic, pBn_R4_all_DW$moyFRic)
      Bn_DistriDW$prof <- as.numeric(rep(c("10", "-5", "-20", "-40", "-70"), each = 4))
      Bn_DistriDW$strate <- rep(c("A", "R1", "R2", "R3", "R4"), each = 4)
      Bn_DistriDW$semin <- Bn_DistriDW$x - 2*Bn_DistriDW$se
      Bn_DistriDW$semax <- Bn_DistriDW$x + 2*Bn_DistriDW$se
      Bn_PR <- ggplot(data = Bn_DistriDW[Bn_DistriDW$strate != "A", ], aes(x = prof, y = x^2, group = div, shape = div))+
        geom_line(aes(linetype = div))+
        geom_point(size = 4)+
        geom_errorbar(aes(ymin = semin^2, ymax = semax^2), width = 0.4)+
        coord_flip()+
        scale_shape_discrete(name="Earthworm FRic")+
        xlim(c(-90, 0))+
        ylab ("B. napus \nroot mass profile (DW mg)")+
        xlab("Soil depth\n(cm)")+
        scale_linetype_discrete(guide=FALSE)+
        #scale_shape_discrete(guide = FALSE)+
        theme(legend.justification=c(0,0), legend.position=c(0.1,0.05), legend.box = "horizontal")+
        theme(panel.background = element_rect(fill="white", colour = "black"))+
        geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")
      png("figures/pl/mass_Bn_rootProfile.png", width = 500, height = 1000)
      Bn_PR
      dev.off()

  ## Shoot/root ratio
  Bn_SR_DW <- myfun_Ind(mdata = SR_ratio, vb = asin(tan(SR_ratio$shootroot)), seuil_prcRE = 30)  # asin(tan()), transf sqrt et log n'am?liorent pas la distribution des r?sidus
  Bn_SR_DW$select
  pBn_SR_DW <- myplot_Ind(mdata = SR_ratio, vb = SR_ratio$shootroot, mod = Bn_SR_DW$mod10A, 
                          ylabl = "B. napus \nShoot:Root ratio (DW)", plot_PAR = T)
  png("figures/pl/mass_Bn_shootroot.png", width = 1000, height = 480)
  grid.arrange(pBn_SR_DW$pDIV, pBn_SR_DW$pLEG, pBn_SR_DW$pDIV_INTER,
               ncol = 3)
  dev.off()
  
  
