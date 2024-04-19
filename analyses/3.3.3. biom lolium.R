#############################################        
###          Pour les ray-grass           ###
#############################################
# Biomasse s?che a?rienne
Lp_A_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 0, ]
Lp_A_all_DW <- myfun_Ind(mdata = Lp_A_all, vb = log1p(Lp_A_all$DW), seuil_prcRE = 24) 
Lp_A_all_DW$resum
Lp_A_all_DW$select
pLp_A_all_DW <- myplot_Ind(mdata = Lp_A_all, vb = sqrt(Lp_A_all$DW), mod = Lp_A_all_DW$mod9A, 
                           ylabl = "L. perenne\n shoot mass (DW mg)", plot_PAR = T)

# Biomasse s?che strate 1
Lp_R1_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 1, ]
# test si les diff?rentes options de "nettoyage" du dataset apporte quelquechose
#Lp_R1_opt2 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 1) , ]
#Lp_R1_opt4 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 1) , ]
Lp_R1_all_DW <- myfun_Ind(mdata = Lp_R1_all, vb = sqrt(Lp_R1_all$DW), seuil_prcRE = 24) 
#Lp_R1_opt2_DW <- myfun_Ind(mdata = Lp_R1_opt2, vb = log1p(Lp_R1_opt2$DW), seuil_prcRE = 30)
#Lp_R1_opt4_DW <- myfun_Ind(mdata = Lp_R1_opt4, vb = log1p(Lp_R1_opt4$DW), seuil_prcRE = 30)
#Lp_R1_comb <- cbind.data.frame(Lp_R1_all_DW$select[,1], Lp_R1_opt2_DW$select[,1], Lp_R1_opt4_DW$select[,1])
#colnames(Lp_R1_comb) <- c("all", "opt2", "opt4")
#   Lp_R1_comb
Lp_R1_all_DW$select
pLp_R1_all_DW <- myplot_Ind(mdata = Lp_R1_all, vb = sqrt(Lp_R1_all$DW), mod = Lp_R1_all_DW$mod9A, 
                            ylabl = "L. perenne\nroot mass (0-10 cm, DW mg)", plot_PAR = F)        # , plot_PAR = T       
# Biomasse s?che strate 2
Lp_R2_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 2, ]
# test si les diff?rentes options de "nettoyage" du dataset apporte quelquechose
#Lp_R2_opt2 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 2) , ]
#Lp_R2_opt4 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 2) , ]
Lp_R2_all_DW <- myfun_Ind(mdata = Lp_R2_all, vb = sqrt(Lp_R2_all$DW), seuil_prcRE = 24) 
#    Lp_R2_opt2_DW <- myfun_Ind(mdata = Lp_R2_opt2, vb = log1p(Lp_R2_opt2$DW), seuil_prcRE = 30)
#   Lp_R2_opt4_DW <- myfun_Ind(mdata = Lp_R2_opt4, vb = log1p(Lp_R2_opt4$DW), seuil_prcRE = 30)
#      Lp_R2_comb <- cbind.data.frame(Lp_R2_all_DW$select[,1], Lp_R2_opt2_DW$select[,1], Lp_R2_opt4_DW$select[,1])
#     colnames(Lp_R2_comb) <- c("all", "opt2", "opt4")
#        Lp_R2_comb
Lp_R2_all_DW$select
pLp_R2_all_DW <- myplot_Ind(mdata = Lp_R2_all, vb = sqrt(Lp_R2_all$DW), mod = Lp_R2_all_DW$mod9A, 
                            ylabl = "L. perenne\nroot mass (10-30 cm, DW mg)", plot_PAR = F)

# Biomasse s?che strate 3
Lp_R3_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 3, ]
# test si les diff?rentes options de "nettoyage" du dataset apporte quelquechose
#Lp_R3_opt2 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 3) , ]
#Lp_R3_opt4 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 3) , ]
Lp_R3_all_DW <- myfun_Ind(mdata = Lp_R3_all, vb = log1p(Lp_R3_all$DW), seuil_prcRE = 24) 
#      Lp_R3_opt2_DW <- myfun_Ind(mdata = Lp_R3_opt2, vb = log1p(Lp_R3_opt2$DW), seuil_prcRE = 30)
#     Lp_R3_opt4_DW <- myfun_Ind(mdata = Lp_R3_opt4, vb = log1p(Lp_R3_opt4$DW), seuil_prcRE = 30)
#        Lp_R3_comb <- cbind.data.frame(Lp_R3_all_DW$select[,1], Lp_R3_opt2_DW$select[,1], Lp_R3_opt4_DW$select[,1])
#       colnames(Lp_R3_comb) <- c("all", "opt2", "opt4")
#          Lp_R3_comb
Lp_R3_all_DW$select
pLp_R3_all_DW <- myplot_Ind(mdata = Lp_R3_all, vb = sqrt(Lp_R3_all$DW), mod = Lp_R3_all_DW$mod9A, 
                            ylabl = "L. perenne\nroot mass (30-50 cm, DW mg)", plot_PAR = F)          
# Biomasse s?che strate 4
Lp_R4_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 4, ]
# test si les diff?rentes options de "nettoyage" du dataset apporte quelquechose
#Lp_R4_opt2 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 4) , ]
#Lp_R4_opt4 <- mybiom_Ind[which(mybiom_Ind$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1] & mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 4) , ]
Lp_R4_all_DW <- myfun_Ind(mdata = Lp_R4_all, vb = log1p(Lp_R4_all$DW), seuil_prcRE = 24) 
#     Lp_R4_opt2_DW <- myfun_Ind(mdata = Lp_R4_opt2, vb = log1p(Lp_R4_opt2$DW), seuil_prcRE = 30)
#    Lp_R4_opt4_DW <- myfun_Ind(mdata = Lp_R4_opt4, vb = log1p(Lp_R4_opt4$DW), seuil_prcRE = 30)
#       Lp_R4_comb <- cbind.data.frame(Lp_R4_all_DW$select[,1], Lp_R4_opt2_DW$select[,1], Lp_R4_opt4_DW$select[,1])
#      colnames(Lp_R4_comb) <- c("all", "opt2", "opt4")
#          Lp_R4_comb
Lp_R4_all_DW$select
pLp_R4_all_DW <- myplot_Ind(mdata = Lp_R4_all, vb = sqrt(Lp_R4_all$DW), mod = Lp_R4_all_DW$mod9A, 
                            ylabl = "L. perenne\nroot mass (50-90 cm, DW mg)", plot_PAR = F)

# Biomasse s?che racine totale
Lp_R <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata != 0, ]
Lp_R <- aggregate(Lp_R$DW, list(RT = Lp_R$RT, numPL = Lp_R$numPL,  div = Lp_R$div,  leg = Lp_R$leg, 
                                id_loca = Lp_R$id_loca, id_taxo = Lp_R$id_taxo, distance = Lp_R$distance, 
                                serie = Lp_R$serie, position = Lp_R$position, PAR = Lp_R$PAR), sum)
colnames(Lp_R)[11] <- "DW"
Lp_R_DW <- myfun_Ind(mdata = Lp_R, vb = log1p(Lp_R$DW), seuil_prcRE = 24) ## log transformed
Lp_R_DW$select
pLp_R_DW <- myplot_Ind(mdata = Lp_R, vb = log(Lp_R$DW), mod = Lp_R_DW$mod9A, 
                       ylabl = "L. perenne root\n dry mass (total, mg, log transf.)", plot_PAR = T)
pLp_R_DW$moyLeg


# Profil racinaire   
Lp_DistriDW <- rbind.data.frame(pLp_A_all_DW$moyFRic, pLp_R1_all_DW$moyFRic, pLp_R2_all_DW$moyFRic, pLp_R3_all_DW$moyFRic, pLp_R4_all_DW$moyFRic)
Lp_DistriDW$prof <- as.numeric(rep(c("10", "-5", "-20", "-40", "-70"), each = 4))
Lp_DistriDW$strate <- rep(c("A", "R1", "R2", "R3", "R4"), each = 4)
Lp_DistriDW$semin <- Lp_DistriDW$x - 2*Lp_DistriDW$se
Lp_DistriDW$semax <- Lp_DistriDW$x + 2*Lp_DistriDW$se
Lp_PR <- ggplot(data = Lp_DistriDW[Lp_DistriDW$strate != "A", ], aes(x = prof, y = x^2, group = div, shape = div))+
  geom_line(aes(linetype = div))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin = semin^2, ymax = semax^2, linetype = div), width = 0.4)+
  coord_flip()+
  scale_shape_discrete(name="Earthworm FRic")+
  xlim(c(-90, 0))+
  ylab ("L. perenne\n root mass profile (DW mg)")+
  xlab("Soil depth\n(cm)")+
  scale_linetype_discrete(guide=FALSE)+
  theme(legend.justification=c(0,0), legend.position=c(0.5,0.05),
        panel.background = element_rect(fill="white", colour = "black"),
        text =  element_text(face = "plain",
                             color = "black", 
                             hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                             margin = margin(), debug = FALSE))+
  geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")
Lp_PR    

# Shoot/root ratio
Lp_SR <- merge(Lp_R, Lp_A_all, by = "id_loca")
Lp_SR <- Lp_SR[, -c(12:17, 19:22)]
colnames(Lp_SR) <- c("id_loca",  "RT", "numPL", "div", "leg", "id_taxo", "distance",
                     "serie", "position", "PAR", "DW_R", "DW_A") 
Lp_SR$SR <- Lp_SR$DW_A/Lp_SR$DW_R
Lp_SR_DW <- myfun_Ind(mdata = Lp_SR, vb = Lp_SR$SR, seuil_prcRE = 30)  # asin fonctionne, mais pb ensuite dans myplot_Ind; transf sqrt et log n'am?liorent pas la distribution des r?sidus
Lp_SR_DW$select
pLp_SR_DW <- myplot_Ind(mdata = Lp_SR, vb = Lp_SR$SR, mod = Lp_SR_DW$mod10A, 
                        ylabl = "L. perenne\n Shoot:Root ratio (DW)", plot_PAR = T)
pLp_SR_DW
