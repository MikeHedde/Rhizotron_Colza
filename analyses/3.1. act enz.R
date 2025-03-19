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
# 2.1 Effet des traitements sur les enzymes du sol
#-----------------------------------------------------------
#############################################        
###         Activit?s enzymatiques        ###
#############################################
myenz <- read.csv("data/derived-data/enz.csv", h = T, sep = ",")

# strate1
myenz_s1 <- myenz[myenz$strate == "1",]

# humidity
myenz_s1_hum <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$hum, seuil_prcRE = 30)  
myenz_s1_hum$select
pmyenz_s1_hum <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$hum, mod = myenz_s1_hum$mod4A, 
                            ylabl = "Humidity", plot_PAR = F) 

# glucosidase
myenz_s1_glu <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$glu, seuil_prcRE = 30)  
myenz_s1_glu$select
pmyenz_s1_glu <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$glu, mod = myenz_s1_glu$mod4A, 
                            ylabl = "Glucosidase activity (AU)", plot_PAR = F) 

# phosphatase
myenz_s1_phos <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$phos, seuil_prcRE = 30)  
myenz_s1_phos$select
pmyenz_s1_phos <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$phos, mod = myenz_s1_phos$mod3A, 
                             ylabl = "phosphatase activity (AU)", plot_PAR = F) 
pmyenz_s1_phos$pDIV
# urease
myenz_s1_ure <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$ure, seuil_prcRE = 30)  
myenz_s1_ure$select
pmyenz_s1_ure <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$ure, mod = myenz_s1_ure$mod4A, 
                            ylabl = "urease activity (AU)", plot_PAR = F) 

# arylsulfatase
myenz_s1_ary <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$ary, seuil_prcRE = 30)  
myenz_s1_ary$select
pmyenz_s1_ary <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$ary, mod = myenz_s1_ary$mod4A, 
                            ylabl = "arylsulfatase activity (AU)", plot_PAR = F)     

# strate2
myenz_s2 <- myenz[myenz$strate == "2",]

# glucosidase
myenz_s2_glu <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$glu, seuil_prcRE = 30)  
myenz_s2_glu$select
pmyenz_s2_glu <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$glu, mod = myenz_s2_glu$mod1A, 
                            ylabl = "Glucosidase activity (AU)", plot_PAR = F) 

# phosphatase
myenz_s2_phos <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$phos, seuil_prcRE = 30)  
myenz_s2_phos$select
pmyenz_s2_phos <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$phos, mod = myenz_s2_phos$mod3A, 
                             ylabl = "phosphatase activity (AU)", plot_PAR = F) 

# urease
myenz_s2_ure <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$ure, seuil_prcRE = 30)  
myenz_s2_ure$select
pmyenz_s2_ure <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$ure, mod = myenz_s2_ure$mod1A, 
                            ylabl = "urease activity (AU)", plot_PAR = F) 

# arylsulfatase
myenz_s2_ary <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$ary, seuil_prcRE = 30)  
myenz_s2_ary$select
pmyenz_s2_ary <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$ary, mod = myenz_s2_ary$mod3A, 
                            ylabl = "arylsulfatase activity (AU)", plot_PAR = F)  

# strate3
myenz_s3 <- myenz[myenz$strate == "3",]

# glucosidase
myenz_s3_glu <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$glu, seuil_prcRE = 30)  
myenz_s3_glu$select
pmyenz_s3_glu <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$glu, mod = myenz_s3_glu$mod1A, 
                            ylabl = "Glucosidase activity (AU)", plot_PAR = F) 

# phosphatase
myenz_s3_phos <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$phos, seuil_prcRE = 30)  
myenz_s3_phos$select
pmyenz_s3_phos <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$phos, mod = myenz_s3_phos$mod3A, 
                             ylabl = "phosphatase activity (AU)", plot_PAR = F) 

# urease
myenz_s3_ure <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$ure, seuil_prcRE = 30)  
myenz_s3_ure$select
pmyenz_s3_ure <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$ure, mod = myenz_s3_ure$mod2A, 
                            ylabl = "urease activity (AU)", plot_PAR = T) 

# arylsulfatase
myenz_s3_ary <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$ary, seuil_prcRE = 30)  
myenz_s3_ary$select
pmyenz_s3_ary <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$ary, mod = myenz_s3_ary$mod1A, 
                            ylabl = "arylsulfatase activity (AU)", plot_PAR = F)              

# strate4
myenz_s4 <- myenz[myenz$strate == "4",]

# glucosidase
myenz_s4_glu <- myfun_Com(mdata = myenz_s4, vb = myenz_s4$glu, seuil_prcRE = 30)  
myenz_s4_glu$select # pas de transf. ad?quate
pmyenz_s4_glu <-myplot_Com(mdata = myenz_s4, vb = myenz_s4$glu, mod = myenz_s4_glu$mod4A, 
                           ylabl = "Glucosidase activity (AU)", plot_PAR = F) 

# phosphatase
myenz_s4_phos <- myfun_Com(mdata = myenz_s4, vb = myenz_s4$phos, seuil_prcRE = 30)  
myenz_s4_phos$select
pmyenz_s4_phos <-myplot_Com(mdata = myenz_s4, vb = myenz_s4$phos, mod = myenz_s4_phos$mod3A, 
                            ylabl = "phosphatase activity (AU)", plot_PAR = F) 

# urease
myenz_s4_ure <- myfun_Com(mdata = myenz_s4, vb = myenz_s4$ure, seuil_prcRE = 30)  
myenz_s4_ure$select
pmyenz_s4_ure <-myplot_Com(mdata = myenz_s4, vb = myenz_s4$ure, mod = myenz_s4_ure$mod2A, 
                           ylabl = "urease activity (AU)", plot_PAR = T) 

# arylsulfatase
myenz_s4_ary <- myfun_Com(mdata = myenz_s4, vb = myenz_s4$ary, seuil_prcRE = 30)  
myenz_s4_ary$select
pmyenz_s4_ary <- myplot_Com(mdata = myenz_s4, vb = myenz_s4$ary, mod = myenz_s4_ary$mod3A, 
                            ylabl = "arylsulfatase activity (AU)", plot_PAR = F) 

# Profil enzymatique en fonction de la profondeur      
enz_Distri <- rbind.data.frame(pmyenz_s1_glu$moyFRic, pmyenz_s2_glu$moyFRic, pmyenz_s3_glu$moyFRic, pmyenz_s4_glu$moyFRic,
                               pmyenz_s1_phos$moyFRic, pmyenz_s2_phos$moyFRic, pmyenz_s3_phos$moyFRic, pmyenz_s4_phos$moyFRic,
                               pmyenz_s1_ure$moyFRic, pmyenz_s2_ure$moyFRic, pmyenz_s3_ure$moyFRic, pmyenz_s4_ure$moyFRic,
                               pmyenz_s1_ary$moyFRic, pmyenz_s2_ary$moyFRic, pmyenz_s3_ary$moyFRic, pmyenz_s4_ary$moyFRic)
enz_Distri$prof <- as.numeric(rep(rep(c("-5", "-20", "-40", "-70"), each = 4), 4))               
decalage <- as.numeric(rep(c(1.5, 0.5, -0.5, -1.5), 16)) 
enz_Distri$prof <- enz_Distri$prof+decalage
enz_Distri$enz <- rep(c("glu", "phos", "ure", "ary"), each = 16)
enz_Distri$semin <- enz_Distri$x - enz_Distri$se
enz_Distri$semax <- enz_Distri$x + enz_Distri$se
enz_PR <- ggplot(data = enz_Distri, aes(x = prof, y = x^2, group = div, shape = div))+
  geom_line(aes(linetype = div))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = semin^2, ymax = semax^2), width = 0.4)+
  coord_flip()+
  xlim(c(-75, 0))+
  ylab ("Enzymatic activity (AU)")+
  xlab("Soil depth\n(cm)")+
  scale_linetype_discrete(guide=FALSE)+
  #theme(legend.justification=c(0,0), legend.position=c(1,1))+
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")+
  #scale_shape_discrete(guide = FALSE)+
  facet_grid(.~enz, scales = "free_x")
enz_PR



#############################################        
###         Stoechiom enzymatique         ###
#############################################
myenz$CN <- log(myenz$glu)/log(myenz$ure)
myenz$CP <- log(myenz$glu)/log(myenz$phos)
myenz$CS <- log(myenz$glu)/log(myenz$ary)
myenz$NP <- log(myenz$ure)/log(myenz$phos)

# strate1
myenz_s1 <- myenz[myenz$strate == "1",]

# C:N
myenz_s1_CN <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$CN, seuil_prcRE = 30)  
myenz_s1_CN$select ## pas de transf. efficace
pmyenz_s1_CN <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$CN, mod = myenz_s1_CN$mod4A, 
                           ylabl = "Enzymatic C:N ratio", plot_PAR = F) 

# C:P
myenz_s1_CP <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$CP, seuil_prcRE = 30)  
myenz_s1_CP$select
pmyenz_s1_CP <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$CP, mod = myenz_s1_CP$mod1A, 
                           ylabl = "Enzymatic C:P ratio", plot_PAR = F) 

# C:S
myenz_s1_CS <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$CS, seuil_prcRE = 30)  
myenz_s1_CS$select
pmyenz_s1_CS <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$CS, mod = myenz_s1_CS$mod3A, 
                           ylabl = "Enzymatic C:S ratio", plot_PAR = F) 

# N:P
myenz_s1_NP <- myfun_Com(mdata = myenz_s1, vb = myenz_s1$NP, seuil_prcRE = 30)  
myenz_s1_NP$select
pmyenz_s1_NP <- myplot_Com(mdata = myenz_s1, vb = myenz_s1$NP, mod = myenz_s1_NP$mod4A, 
                           ylabl = "Enzymatic N:P ratio", plot_PAR = F)     

# strate2
myenz_s2 <- myenz[myenz$strate == "2",]

# C:N
myenz_s2_CN <- myfun_Com(mdata = myenz_s2, vb = sqrt(myenz_s2$CN), seuil_prcRE = 30)  
myenz_s2_CN$select# pas de transf. efficace => KW ??
pmyenz_s2_CN <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$CN, mod = myenz_s2_CN$mod4A, 
                           ylabl = "Enzymatic C:N ratio", plot_PAR = F) 

# C:P
myenz_s2_CP <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$CP, seuil_prcRE = 30)  
myenz_s2_CP$select
pmyenz_s2_CP <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$CP, mod = myenz_s2_CP$mod3A, 
                           ylabl = "Enzymatic C:P ratio", plot_PAR = F) 

# C:S
myenz_s2_CS <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$CS, seuil_prcRE = 30)  
myenz_s2_CS$select
pmyenz_s2_CS <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$CS, mod = myenz_s2_CS$mod3A, 
                           ylabl = "Enzymatic C:S ratio", plot_PAR = F) 

# N:P
myenz_s2_NP <- myfun_Com(mdata = myenz_s2, vb = myenz_s2$NP, seuil_prcRE = 30)  
myenz_s2_NP$select
pmyenz_s2_NP <- myplot_Com(mdata = myenz_s2, vb = myenz_s2$NP, mod = myenz_s2_NP$mod2A, 
                           ylabl = "Enzymatic N:P ratio", plot_PAR = T)  

# strate3
myenz_s3 <- myenz[myenz$strate == "3",]

# C:N
myenz_s3_CN <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$CN, seuil_prcRE = 30)  
myenz_s3_CN$select #pas de transf. efficaces => KW ??
pmyenz_s3_CN <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$CN, mod = myenz_s3_CN$mod4A, 
                           ylabl = "Enzymatic C:N ratio", plot_PAR = F) 

# C:P
myenz_s3_CP <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$CP, seuil_prcRE = 30)  
myenz_s3_CP$select
pmyenz_s3_CP <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$CP, mod = myenz_s3_CP$mod3A, 
                           ylabl = "Enzymatic C:P ratio", plot_PAR = F) 

# C:S
myenz_s3_CS <- myfun_Com(mdata = myenz_s3, vb = myenz_s3$CS, seuil_prcRE = 30)  
myenz_s3_CS$select
pmyenz_s3_CS <- myplot_Com(mdata = myenz_s3, vb = myenz_s3$CS, mod = myenz_s3_CS$mod3A, 
                           ylabl = "Enzymatic C:S ratio", plot_PAR = T) 

# N:P
myenz_s3_NP <- myfun_Com(mdata = myenz_s3, vb = log1p(myenz_s3$NP), seuil_prcRE = 30)  
myenz_s3_NP$select
pmyenz_s3_NP <- myplot_Com(mdata = myenz_s3, vb = log1p(myenz_s3$NP), mod = myenz_s3_NP$mod1A, 
                           ylabl = "Enzymatic N:P ratio", plot_PAR = F)              

# strate4
myenz_s4 <- myenz[myenz$strate == "4",]

# C:N
myenz_s4_CN <- myfun_Com(mdata = myenz_s4[-39,], vb = myenz_s4$CN[-39], seuil_prcRE = 30)  
myenz_s4_CN$select # pas de transf. ad?quate
pmyenz_s4_CN <-myplot_Com(mdata = myenz_s4[-39,], vb = myenz_s4$CN[-39], mod = myenz_s4_CN$mod4A, 
                          ylabl = "Enzymatic C:N ratio", plot_PAR = F) 

# C:P
myenz_s4_CP <- myfun_Com(mdata = myenz_s4[-39,], vb = myenz_s4$CP[-39], seuil_prcRE = 30)  
myenz_s4_CP$select
pmyenz_s4_CP <-myplot_Com(mdata = myenz_s4[-39,], vb = myenz_s4$CP[-39], mod = myenz_s4_CP$mod1A, 
                          ylabl = "Enzymatic C:P ratio", plot_PAR = F) 

# C:S
myenz_s4_CS <- myfun_Com(mdata = myenz_s4[-39,], vb = myenz_s4$CS[-39], seuil_prcRE = 30)  
myenz_s4_CS$select
pmyenz_s4_CS <-myplot_Com(mdata = myenz_s4[-39,], vb = myenz_s4$CS[-39], mod = myenz_s4_CS$mod3A, 
                          ylabl = "Enzymatic C:S ratio", plot_PAR = T) 

# N:P
myenz_s4_NP <- myfun_Com(mdata = myenz_s4, vb = myenz_s4$NP, seuil_prcRE = 30)  
myenz_s4_NP$select
pmyenz_s4_NP <- myplot_Com(mdata = myenz_s4, vb = myenz_s4$NP, mod = myenz_s4_NP$mod2A, 
                           ylabl = "Enzymatic N:P ratio", plot_PAR = F) 

# Profil enzymatique en fonction de la profondeur      
enzSto_Distri <- rbind.data.frame(pmyenz_s1_CN$moyFRic, pmyenz_s2_CN$moyFRic, pmyenz_s3_CN$moyFRic, pmyenz_s4_CN$moyFRic,
                                  pmyenz_s1_CP$moyFRic, pmyenz_s2_CP$moyFRic, pmyenz_s3_CP$moyFRic, pmyenz_s4_CP$moyFRic,
                                  pmyenz_s1_CS$moyFRic, pmyenz_s2_CS$moyFRic, pmyenz_s3_CS$moyFRic, pmyenz_s4_CS$moyFRic,
                                  pmyenz_s1_NP$moyFRic, pmyenz_s2_NP$moyFRic, pmyenz_s3_NP$moyFRic, pmyenz_s4_NP$moyFRic)
enzSto_Distri$prof <- as.numeric(rep(rep(c("-5", "-20", "-40", "-70"), each = 4), 4))               
decalage <- as.numeric(rep(c(1.5, 0.5, -0.5, -1.5), 16)) 
enzSto_Distri$prof <- enz_Distri$prof+decalage
enzSto_Distri$enz <- rep(c("C:N", "C:P", "C:S", "N:P"), each = 16)
enzSto_Distri$semin <- enzSto_Distri$x - enzSto_Distri$se
enzSto_Distri$semax <- enzSto_Distri$x + enzSto_Distri$se
enzSto_PR <- ggplot(data = enzSto_Distri, aes(x = prof, y = x^2, group = div, shape = div))+
  geom_line(aes(linetype = div))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = semin^2, ymax = semax^2), width = 0.4)+
  coord_flip()+
  xlim(c(-90, 0))+
  ylab ("Enzymatic stoechiometry")+
  xlab("Soil depth\n(cm)")+
  scale_linetype_discrete(guide=FALSE)+
  theme(legend.justification=c(0,0), legend.position=c(0.1, 0))+
  theme(panel.background = element_rect(fill="white", colour = "black"))+
  geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")+
  #scale_shape_discrete(guide = FALSE)+
  facet_grid(~enz, scales = "free_x")
enzSto_PR
