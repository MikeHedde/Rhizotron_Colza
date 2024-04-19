




#-----------------------------------------------------------
# 4.1 Relations entre traits (moyenne par rhizotron)
#------------------------------------------------------------                
# Rappel des df cr??es plus haut
# myvdt -> biom vdt
# myburrow <- r?seau de galeries des vdt
# mybiom_Com -> biom DW communaut? v?g?tale          
# mybiom_Com -> biom DW our chaque esp par strate    
mybiom_Ind2 <- dcast(mybiom_Ind, strata + RT + div + leg + position + serie ~ id_taxo, value.var = "DW", fun.aggregate = mean, na.rm = T)
# Ta_N -> [N] dans les organes de T. alexandrium
# Bn_N -> [N] dans les organes de B. napus
# Lp_N -> [N] dans les organes de L. perenne
############################################        
### Regressions entre [N] dans Ta et Bn  ###
############################################ 
colnames(Bn_N)[11] <- "fRoots"
colnames(Ta_N)[11] <- "fRoots"
k <- aggregate(cbind(Leaves, fRoots, Limb, Stem) ~ leg+div+PAR+ RT, data = Bn_N, mean, na.rm = TRUE)
k <- k[k$leg == "Leg+", ]
k[22:42, 1:4] <- k[1:21, 1:4]
k$esp <- rep(c("Bn", "Ta"), each = 21)
k$Leaves[22:42] <- aggregate(Leaves ~ RT, data = Ta_N, mean)$Leaves
#k$Limb[22:42] <- aggregate(Limb ~ RT, data = Ta_N, mean)$Limb   ## il n'a pas de limbes matures dans tous les RT
k$Stem[22:42] <- aggregate(Stem ~ RT, data = Ta_N, mean)$Stem
k$fRoots[22:42] <- aggregate(fRoots ~ RT, data = Ta_N, mean)$fRoots
k2 <- melt(data = k, id = c("div", "PAR", "RT", "esp"), measure = c("Leaves", "Stem", "fRoots", "Limb"), variable.name = "organes", na.rm = T)
k2 <- dcast(k2, div + PAR + RT + organes ~ esp)

# pour les organes a?riens
summary(lm(Bn~Ta, data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div0"),])) #
summary(lm(Bn~Ta, data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div1"),])) # R^2 = 0.44*
summary(lm(sqrt(Bn)~sqrt(Ta), data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div2"),])) # R^2 = 0.18$
summary(lm(Bn~Ta, data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div3"),])) # 

# pour les racines fines
summary(lm(Bn~Ta, data = k2[which(k2$organes == "fRoots" & k2$div == "Div0"),])) #
summary(lm(sqrt(Bn)~sqrt(Ta), data = k2[which(k2$organes == "fRoots" & k2$div == "Div1") ,]) # 
        summary(lm(sqrt(Bn)~sqrt(Ta), data = k2[which(k2$organes == "fRoots" & k2$div == "Div2"),])) # R^2 = 0.65$ si l'outlier est supprim? (RT 42)
        summary(lm(Bn~Ta, data = k2[which(k2$organes == "fRoots" & k2$div == "Div3"),])) # 
        
        # Repr?sentations graphiques
        Reg1 <- ggplot(data = k2[k2$organes == c("Leaves", "Stem"),], aes(y = Bn, x = Ta))+
          geom_point(aes(shape = div))+
          geom_smooth(method = "lm")+
          facet_grid(div~.)+
          ylab("B. napus\n Shoot organs N content (mg g-1)")+
          xlab("T. alexandrium\n Shoot organs N content (mg g-1)")+
          theme(panel.background = element_rect(fill="white", colour = "black"))+
          scale_shape_discrete(guide = FALSE)
        
        Reg2 <- ggplot(data = k2[which(k2$organes == "fRoots"),], aes(y = Bn, x = Ta))+ # voir pourquoi il y a un Ta ? 4 mg g-1 (RT20)??
          geom_point(aes(shape = div))+
          geom_smooth(method = "lm")+
          facet_grid(div~.)+
          ylab("B. napus\n fine roots N content (mg g-1)")+
          xlab("T. alexandrium\n fine roots N content (mg g-1)")+
          theme(panel.background = element_rect(fill="white", colour = "black"))+
          scale_shape_discrete(guide = FALSE)
        
        
        
        #-----------------------------------------------------------
        # 5.1 Liste des figures ? proposer
        #------------------------------------------------------------    
        
        # Figure Racine
        grid.arrange(pTa_A_DW$pDIV, pBn_A_all_DW$pDIV, pLp_A_all_DW$pDIV,
                     Ta_PR, Bn_PR, Lp_PR, 
                     ncol = 3, nrow = 2)
        # Figure activit? enzymatiques
        enz_PR        
        enzSto_PR
        
        # Figure N content
        grid.arrange(pBn_N_limb$pDIV, pBn_N_leaves$pDIV, pBn_N_stem$pDIV, pBn_N_fRoots$pDIV, pBn_N_tRoots$pDIV,
                     pTa_N_limb$pDIV, pTa_N_leaves$pDIV, pTa_N_stem$pDIV, pTa_N_fRoots$pDIV,
                     ncol = 5, nrow = 2)
        
        grid.arrange(Reg1, Reg2, ncol = 2)
        
        
        
        
        ####################################################
        
        dev.off()
        
        
        ###########################
        #          ACP            #
        ###########################
        # Calcul Shanon sur biomasses racinaires
        require(vegan)
        biom_R_ind_strata <- aggregate(DW_R ~ strata + RT + id_taxo, sum, data = biom_R_Ind, na.rm=T)
        biom_R_shan <- aggregate(DW_R ~ RT + strata, diversity, data = biom_R_ind_strata)
        biom_R_shan$S <- aggregate(DW_R ~ RT + strata, specnumber, data = biom_R_ind_strata)$DW_R
        biom_R_shan$J  <- biom_R_shan$DW_R/log(biom_R_shan$S)
        biom_R_J_mean <- aggregate(J ~ RT, mean, data = biom_R_shan)
        
        mydf_PCA <- mybiom_Com[,1:3]
        mydf_PCA$div2 <- as.numeric(substr(mydf_PCA$div, 4,4))
        mydf_PCA$div2 <- ifelse(mydf_PCA$div2 > 1, "high FD", "low FD")
        
        mydf_PCA$DW_A <- mybiom_Com$DW_A
        mydf_PCA$DW_R <- mybiom_Com$DW_R
        mydf_PCA$DW_evenness <- biom_R_J_mean$J
        mydf_PCA$tot_Root_length <- lgrac$longR
        
        pivot_vs_fines <- aggregate(DW ~ div + leg + organe + RT, mean, data = pl[pl$id_taxo == "Bra_napus",], na.rm = T)
        mydf_PCA$Bn_tR_DW <- pivot_vs_fines$DW[pivot_vs_fines$organe == "Taproot"]
        mydf_PCA$Bn_fR_DW <- pivot_vs_fines$DW[pivot_vs_fines$organe == "Fine roots"]
        mydf_PCA$Bn_R_DW <- aggregate(Bn_R$DW, list(Bn_R$RT), mean, na.omit= T)$x
        mydf_PCA$Lp_R_DW <- aggregate(Lp_R$DW, list(Lp_R$RT), mean, na.omit= T)$x
        mydf_PCA$Max_Root_Bn <- lgrac$Max_root_depth_Bn
        mydf_PCA$Max_Root_Lp <- lgrac$Max_root_depth_Lp
        
        mydf_PCA$Bn_qN <- aggregate(Bn_qN$qN, list(Bn_qN$RT), mean, na.omit= T)$x
        mydf_PCA$Lp_qN <- aggregate(Lp_qN$qN, list(Lp_qN$RT), mean, na.omit= T)$x
        mydf_PCA$Bn_N_limb <- aggregate(Bn_N$Limb, list(Bn_N$RT), mean, na.rm= T)$x
        mydf_PCA$Lp_N_limb <- aggregate(Lp_N$Limb, list(Lp_N$RT), mean, na.rm= T)$x
        
        mydf_PCA$Bn_qP <- aggregate(Bn_qP$qP, list(Bn_qP$RT), mean, na.omit= T)$x
        mydf_PCA$Lp_qP <- aggregate(Lp_qP$qP, list(Lp_qP$RT), mean, na.omit= T)$x
        mydf_PCA$Bn_N_limb <- aggregate(Bn_N$Limb, list(Bn_N$RT), mean, na.rm= T)$x
        mydf_PCA$Lp_N_limb <- aggregate(Lp_N$Limb, list(Lp_N$RT), mean, na.rm= T)$x
        
        mydf_PCA$Bn_NP <- mydf_PCA$Bn_qN/mydf_PCA$Bn_qP
        mydf_PCA$Lp_NP <- mydf_PCA$Lp_qN/mydf_PCA$Lp_qP
        
        mydf_PCA <- merge(mydf_PCA, myvdt[,c(1, 12, 14, 18)], "RT", all = T)
        mydf_PCA[is.na(mydf_PCA)] <- 0
        
        mydf_PCA$NO3_s4  <- myanion_s4$NO3
        mydf_PCA$NO3_s3  <- myanion_s3$NO3
        mydf_PCA$NO3_s2  <- myanion_s2$NO3
        mydf_PCA$NO3_s1  <- myanion_s1$NO3
        
        mydf_PCA$HPO42_s4  <- myanion_s4$HPO42
        mydf_PCA$HPO42_s3  <- myanion_s3$HPO42
        mydf_PCA$HPO42_s2  <- myanion_s2$HPO42
        mydf_PCA$HPO42_s1  <- myanion_s1$HPO42
        
        mydf_PCA$URE_s4 <- myenz_s4$ure
        mydf_PCA$URE_s3 <- myenz_s3$ure
        mydf_PCA$URE_s2 <- myenz_s2$ure
        mydf_PCA$URE_s1 <- myenz_s1$ure
        
        
        mydf_PCA <- mydf_PCA[-c(5, 2, 25),]
        
        pdf()
        mydf_PCA_soil <- mydf_PCA[, c("long", "NO3_s4", "NO3_s3", "NO3_s2", "NO3_s1")]
        soilPCA <- dudi.pca(mydf_PCA_soil, scannf = F, nf = 4)
        par(mfrow = c(nr = 2, nc = 2))
        s.corcircle(soilPCA$co)
        s.class(soilPCA$li, mydf_PCA$leg, cellipse= F, col = c("red", "blue"))
        s.class(soilPCA$li, mydf_PCA$div, cellipse= F)
        s.class(soilPCA$li, as.factor(paste(mydf_PCA$leg, mydf_PCA$div2, sep = "/")), cellipse= F, col = c("firebrick1", "firebrick4", "dodgerblue", "dodgerblue4"))
        
        mydf_PCA_plant <- mydf_PCA[, c("tot_Root_length", "Lp_R_DW", "DW_evenness",
                                       "Bn_R_DW", "Bn_qN","Lp_qN")]
        plantPCA <- dudi.pca(mydf_PCA_plant, scannf = F, nf = 4)
        par(mfrow = c(nr = 2, nc = 2))
        s.corcircle(plantPCA$co)
        s.class(plantPCA$li, as.factor(mydf_PCA$leg), cellipse= F, col = c("red", "blue"))
        s.class(plantPCA$li, as.factor(mydf_PCA$div), cellipse= F)
        s.class(plantPCA$li, as.factor(paste(mydf_PCA$leg, mydf_PCA$div2, sep = "/")), cellipse= F, col = c("firebrick1", "firebrick4", "dodgerblue", "dodgerblue4"))
        
        
        coi1 <- coinertia(plantPCA, soilPCA, scannf = F, nf = 2)
        coi1$RV
        coi1$eig/sum(coi1$eig)
        rd <- randtest(coi1, nrepet = 999, fixed = 2)
        rd
        #plot(rd)
        plot(coi1)
        
        s.arrow(coi1$l1, clab = 1.2)
        s.arrow(coi1$c1, clab = 1.2)
        
        colnames(coi1$lX) <- colnames(coi1$lY) 
        s.match(coi1$lX, coi1$lY)
        dev.off()
        
        