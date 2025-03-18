#####################"

## analyses biomasse plante / Rhizotron


# MZ
# 28/09/2016 : reprend SCRIPT_BIOM_PL_20160928.r au propre
# 6/12/2016 : complete stat effet traitement sur biomass pop
# 26/06/2017 : refait analyse apr?s elimination des individus trop petits (germintaion tardive)

# 11/08/2017 : reprends avec les stas de D. Makowsky
# 16/08/2017 : modif script + cr?ation de fonctions par MH

################################################################################


#-----------------------------------------------------------
# 2.2 Effet des traitements sur les anions du sol
#-----------------------------------------------------------
myanion <- myprep(anion)

#############################################        
###                anions                ###
#############################################
# strate1
myanion_s1 <- myanion[myanion$strate == "1",]

# NO3
myanion_s1_NO3 <- myfun_Com(mdata = myanion_s1, vb = myanion_s1$NO3, seuil_prcRE = 30)  
myanion_s1_NO3$select
pmyanion_s1_NO3 <- myplot_Com(mdata = myanion_s1, vb = myanion_s1$NO3, mod = myanion_s1_NO3$mod4A, 
                              ylabl = "NO3", plot_PAR = F) 

# HPO42
myanion_s1_HPO42 <- myfun_Com(mdata = myanion_s1, vb = myanion_s1$HPO42, seuil_prcRE = 30)  
myanion_s1_HPO42$select
pmyanion_s1_HPO42 <- myplot_Com(mdata = myanion_s1, vb = myanion_s1$HPO42, mod = myanion_s1_HPO42$mod3A, 
                                ylabl = "HPO42", plot_PAR = F) 

# SO42
myanion_s1_SO42 <- myfun_Com(mdata = myanion_s1, vb = myanion_s1$SO42, seuil_prcRE = 30)  
myanion_s1_SO42$select
pmyanion_s1_SO42 <- myplot_Com(mdata = myanion_s1, vb = myanion_s1$SO42, mod = myanion_s1_SO42$mod4A, 
                               ylabl = "SO42", plot_PAR = F) 

# Cl
myanion_s1_Cl <- myfun_Com(mdata = myanion_s1, vb = myanion_s1$Cl, seuil_prcRE = 30)  
myanion_s1_Cl$select
pmyanion_s1_Cl <- myplot_Com(mdata = myanion_s1, vb = myanion_s1$Cl, mod = myanion_s1_Cl$mod4A, 
                             ylabl = "Cl", plot_PAR = F)     

# strate2
myanion_s2 <- myanion[myanion$strate == "2",]

# NO3
myanion_s2_NO3 <- myfun_Com(mdata = myanion_s2, vb = myanion_s2$NO3, seuil_prcRE = 30)  
myanion_s2_NO3$select
pmyanion_s2_NO3 <- myplot_Com(mdata = myanion_s2, vb = myanion_s2$NO3, mod = myanion_s2_NO3$mod1A, 
                              ylabl = "NO3cosidase activity (AU)", plot_PAR = F) 

# HPO42
myanion_s2_HPO42 <- myfun_Com(mdata = myanion_s2, vb = myanion_s2$HPO42, seuil_prcRE = 30)  
myanion_s2_HPO42$select
pmyanion_s2_HPO42 <- myplot_Com(mdata = myanion_s2, vb = myanion_s2$HPO42, mod = myanion_s2_HPO42$mod3A, 
                                ylabl = "HPO42 activity (AU)", plot_PAR = F) 

# SO42
myanion_s2_SO42 <- myfun_Com(mdata = myanion_s2, vb = myanion_s2$SO42, seuil_prcRE = 30)  
myanion_s2_SO42$select
pmyanion_s2_SO42 <- myplot_Com(mdata = myanion_s2, vb = myanion_s2$SO42, mod = myanion_s2_SO42$mod1A, 
                               ylabl = "SO42 activity (AU)", plot_PAR = F) 

# Cl
myanion_s2_Cl <- myfun_Com(mdata = myanion_s2, vb = myanion_s2$Cl, seuil_prcRE = 30)  
myanion_s2_Cl$select
pmyanion_s2_Cl <- myplot_Com(mdata = myanion_s2, vb = myanion_s2$Cl, mod = myanion_s2_Cl$mod3A, 
                             ylabl = "Cl activity (AU)", plot_PAR = F)  

# strate3
myanion_s3 <- myanion[myanion$strate == "3",]

# NO3
myanion_s3_NO3 <- myfun_Com(mdata = myanion_s3, vb = myanion_s3$NO3, seuil_prcRE = 30)  
myanion_s3_NO3$select
pmyanion_s3_NO3 <- myplot_Com(mdata = myanion_s3, vb = myanion_s3$NO3, mod = myanion_s3_NO3$mod1A, 
                              ylabl = "NO3cosidase activity (AU)", plot_PAR = F) 

# HPO42
myanion_s3_HPO42 <- myfun_Com(mdata = myanion_s3, vb = myanion_s3$HPO42, seuil_prcRE = 30)  
myanion_s3_HPO42$select
pmyanion_s3_HPO42 <- myplot_Com(mdata = myanion_s3, vb = myanion_s3$HPO42, mod = myanion_s3_HPO42$mod3A, 
                                ylabl = "HPO42 activity (AU)", plot_PAR = F) 

# SO42
myanion_s3_SO42 <- myfun_Com(mdata = myanion_s3, vb = myanion_s3$SO42, seuil_prcRE = 30)  
myanion_s3_SO42$select
pmyanion_s3_SO42 <- myplot_Com(mdata = myanion_s3, vb = myanion_s3$SO42, mod = myanion_s3_SO42$mod2A, 
                               ylabl = "SO42 activity (AU)", plot_PAR = T) 

# Cl
myanion_s3_Cl <- myfun_Com(mdata = myanion_s3, vb = myanion_s3$Cl, seuil_prcRE = 30)  
myanion_s3_Cl$select
pmyanion_s3_Cl <- myplot_Com(mdata = myanion_s3, vb = myanion_s3$Cl, mod = myanion_s3_Cl$mod1A, 
                             ylabl = "Cl activity (AU)", plot_PAR = F)              

# strate4
myanion_s4 <- myanion[myanion$strate == "4",]

# NO3
myanion_s4_NO3 <- myfun_Com(mdata = myanion_s4, vb = myanion_s4$NO3, seuil_prcRE = 30)  
myanion_s4_NO3$select # pas de transf. ad?quate
pmyanion_s4_NO3 <-myplot_Com(mdata = myanion_s4, vb = myanion_s4$NO3, mod = myanion_s4_NO3$mod4A, 
                             ylabl = "NO3cosidase activity (AU)", plot_PAR = F) 

# HPO42
myanion_s4_HPO42 <- myfun_Com(mdata = myanion_s4, vb = myanion_s4$HPO42, seuil_prcRE = 30)  
myanion_s4_HPO42$select
pmyanion_s4_HPO42 <-myplot_Com(mdata = myanion_s4, vb = myanion_s4$HPO42, mod = myanion_s4_HPO42$mod3A, 
                               ylabl = "HPO42 activity (AU)", plot_PAR = F) 

# SO42
myanion_s4_SO42 <- myfun_Com(mdata = myanion_s4, vb = myanion_s4$SO42, seuil_prcRE = 30)  
myanion_s4_SO42$select
pmyanion_s4_SO42 <-myplot_Com(mdata = myanion_s4, vb = myanion_s4$SO42, mod = myanion_s4_SO42$mod2A, 
                              ylabl = "SO42 activity (AU)", plot_PAR = T) 

# Cl
myanion_s4_Cl <- myfun_Com(mdata = myanion_s4, vb = myanion_s4$Cl, seuil_prcRE = 30)  
myanion_s4_Cl$select
pmyanion_s4_Cl <- myplot_Com(mdata = myanion_s4, vb = myanion_s4$Cl, mod = myanion_s4_Cl$mod3A, 
                             ylabl = "Cl activity (AU)", plot_PAR = F) 

# Profil anionique en fonction de la profondeur      
anion_Distri <- rbind.data.frame(pmyanion_s1_NO3$moyFRic, pmyanion_s2_NO3$moyFRic, pmyanion_s3_NO3$moyFRic, pmyanion_s4_NO3$moyFRic,
                                 pmyanion_s1_HPO42$moyFRic, pmyanion_s2_HPO42$moyFRic, pmyanion_s3_HPO42$moyFRic, pmyanion_s4_HPO42$moyFRic,
                                 pmyanion_s1_SO42$moyFRic, pmyanion_s2_SO42$moyFRic, pmyanion_s3_SO42$moyFRic, pmyanion_s4_SO42$moyFRic,
                                 pmyanion_s1_Cl$moyFRic, pmyanion_s2_Cl$moyFRic, pmyanion_s3_Cl$moyFRic, pmyanion_s4_Cl$moyFRic)
anion_Distri$prof <- as.numeric(rep(rep(c("-5", "-20", "-40", "-70"), each = 4), 4))               
decalage <- as.numeric(rep(c(1.5, 0.5, -0.5, -1.5), 16)) 
anion_Distri$prof <- anion_Distri$prof+decalage
anion_Distri$anion <- rep(c("NO3", "HPO42", "SO42", "Cl"), each = 16)
anion_Distri$semin <- anion_Distri$x - anion_Distri$se
anion_Distri$semax <- anion_Distri$x + anion_Distri$se
anion_PR <- ggplot(data = anion_Distri, aes(x = prof, y = x^2, group = div, shape = div))+
  geom_line(aes(linetype = div))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = semin^2, ymax = semax^2), width = 0.4)+
  coord_flip()+
  xlim(c(-90, 0))+
  ylab ("Anion concentration")+
  xlab("Soil depth\n(cm)")+
  scale_linetype_discrete(guide=FALSE)+
  theme(legend.justification=c(0,0), legend.position=c(1,1))+
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")+
  #scale_shape_discrete(guide = FALSE)+
  facet_grid(~anion, scales = "free_x")
anion_PR

