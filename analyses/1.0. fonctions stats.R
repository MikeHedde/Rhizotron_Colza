#------------------------------------------------------------
#     Fonctions creees pour automatiser le travail de mise en
#     forme des donnees, de selection du modele le plus adapte 
#     et de representation graphique de sortie des modeles
#------------------------------------------------------------

# fonction pour ajouter le facteur id_loca
myloca <- function(X){
  # on introduit une colonne loc
  X$numPL <- NULL
  
  # pour remplacer le tr?s long code (avec risque d'erreur de saisie) que tu avais ?crit  
  # pour attribuer une valeur de localisation, je te propose:
  pos <- seq(-24, 24, by = 3)    
  for (j in pos){    
    X$numPL[X$serie=="s1" & X$leg=="Leg-" & X$distance == j] <- which(pos == j)
    X$numPL[X$serie=="s2" & X$leg=="Leg+" & X$distance == j] <- which(pos == j)
    X$numPL[X$serie=="s1" & X$leg=="Leg+" & X$distance == j] <- which(pos == -j)
    X$numPL[X$serie=="s2" & X$leg=="Leg-" & X$distance == j] <- which(pos == -j)
  }
  X$numPL <- as.factor(X$numPL)
  X
}

# Mise en forme des df avec ajout des facteurs serie, id_loca et PAR
myprep <- function(mydata){
  
  # on ajoute une colonne s?rie
  mydata$serie = NULL
  mydata$serie[mydata$RT< 21] = "s1"
  mydata$serie[mydata$RT> 20] = "s2"
  mydata$serie <- as.factor(mydata$serie)
  
  # on ajoute une colonne position du rhizotron dans la chambre
  mydata$position=NULL
  mydata$position[mydata$RT==1|mydata$RT==20|mydata$RT==21|mydata$RT==40]="1"
  mydata$position[mydata$RT==2|mydata$RT==19|mydata$RT== 22|mydata$RT== 39]="2"
  mydata$position[mydata$RT==3|mydata$RT==18|mydata$RT== 23|mydata$RT== 38]="3"
  mydata$position[mydata$RT==4|mydata$RT==17|mydata$RT== 24|mydata$RT== 37]="4"
  mydata$position[mydata$RT==5|mydata$RT==16|mydata$RT==25|mydata$RT==36]="5"
  mydata$position[mydata$RT==6|mydata$RT==15|mydata$RT==26|mydata$RT==35]="6"
  mydata$position[mydata$RT==7|mydata$RT==14|mydata$RT==27|mydata$RT==34]="7"
  mydata$position[mydata$RT==8|mydata$RT==13|mydata$RT==28|mydata$RT==33]="8"
  mydata$position[mydata$RT==9|mydata$RT==12|mydata$RT==29|mydata$RT==32]="9"
  mydata$position[mydata$RT==10|mydata$RT==11|mydata$RT==30|mydata$RT==31]="10"
  mydata$position[mydata$RT==41|mydata$RT==42]="11"
  mydata <- mydata[order(mydata$position, decreasing=F),]
  mydata$position <- as.factor(mydata$position)
  
  # on ajoute le rayonnement en covariable (? modifier quand les valeurs de PAR de la 2e s?rie seront dispo)
  parmoy <- aggregate(par$PAR, by=list(RT=par$RT), mean, na.rm=T)
  parmoy2 <- parmoy
  
  # manipulations pour attribuer une valeur mod?lis?e ? chaque position de plante
  RTtemp = seq(39,1,-2)
  parmoy2$RT <- parmoy$RT + RTtemp
  parmoy3 <- rbind(parmoy, parmoy2, parmoy2[which(parmoy2$RT == c(30:31)),])
  parmoy3$RT[41:42] <- c(42, 41)
  
  mydata <- merge(mydata, parmoy3, by=c("RT"))
  colnames(mydata)[ncol(mydata)] <- "PAR"
  
  mydata
}

# Creation d'une fonction pour retenir le modele le plus adapt? 
# parmi les mod?les d?finis (et en intervertissant l'ordre leg / div)  
# Les arguments sont :
# mdata => df contenant l'ensemble des valeurs pour le calcul des mod?les
# vb => variable d?pendante
# seuil_prcRE => valeur minimale souhait?e pour le % de 
#                variance expliquee par le random effect

# rq : contraste par d?faut = options(contrasts = c("contr.treatment", "contr.treatment"))
options(contrasts =c ("contr.sum","contr.poly"))  # A l'?chelle de la communaut?

# A l'echelle de la communaute        
myfun_Com <- function(mdata, vb, seuil_prcRE){
  
  # Diff?rents modeles complets avec et sans covariable pour s?lectionner le mod?le le plus adapt?
  mod1A <- lmer(vb ~ -1 + div*leg + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod2A <- lmer(vb ~ -1 + div*leg + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod3A <- lm(vb ~ -1 + div*leg + serie, data = mdata, na.action = na.omit) 
  mod4A <- lm(vb ~ -1 + div*leg + serie + PAR, data = mdata, na.action = na.omit)
  mod1B <- lmer(vb ~ -1 + leg*div + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod2B <- lmer(vb ~ -1 + leg*div + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod3B <- lm(vb ~ -1 + leg*div + serie, data = mdata, na.action = na.omit) 
  mod4B <- lm(vb ~ -1 + leg*div + serie + PAR, data = mdata, na.action = na.omit) 
  
  # Les param?tres ? comparer sont collect?s dans une table
  resum <- matrix(NA, ncol = 8, nrow = 17, 
                  dimnames = list(c("position (%)", "residual (%)", "AIC", "BIC", "p-val", "R2", "Shapiro", 
                                    "Type II: fac1", "Type II: fac2", "Type II: serie", "Type II: PAR", "Type II: inter",
                                    "Type III: fac1", "Type III: fac2", "Type III: serie", "Type III: PAR", "Type III: inter"), 
                                  c("mod1A", "mod1B", "mod2A", "mod2B", "mod3A", "mod3B", "mod4A", "mod4B")))
  
  resum[3,] <- round(t(AIC(mod1A, mod1B, mod2A, mod2B, mod3A, mod3B, mod4A, mod4B)[,2]), 0)
  resum[4,] <- round(t(BIC(mod1A, mod1B, mod2A, mod2B, mod3A, mod3B, mod4A, mod4B)[,2]), 0)
  
  # pour les lmer
  lmerlist <- list(mod1A, mod1B, mod2A, mod2B)
  allRE_prc <- lapply(lmerlist, function(x) 
    as.data.frame(VarCorr(x))$vcov/sum(as.data.frame(VarCorr(x))$vcov))
  for (i in c(1:4)){
    resum[1:2,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  
  # pour les lm
  lmlist <- list(mod3A, mod3B, mod4A, mod4B)
  for (i in c(1:4)){
    resum[5, i+4] <- round(lapply(lmlist, function(x) anova(x)$'Pr(>F)'[1])[[i]], 3)
    resum[6, i+4] <- round(lapply(lmlist, function(x) round(summary(x)[[9]], 2))[[i]], 2)
  }
  
  # les lmer sont recalcul?s ont posant REML = T
  mod1A <- lmer(vb ~ -1 + div*leg + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod2A <- lmer(vb ~ -1 + div*leg + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod1B <- lmer(vb ~ -1 + leg*div + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
  mod2B <- lmer(vb ~ -1 + leg*div + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
  
  y <- list(mod1A, mod1B, mod2A, mod2B, mod3A, mod3B, mod4A, mod4B)
  
  # v?rification de la distribution des r?sidus
  resum[7, ] <- round(unlist(lapply(y, function(x) shapiro.test(residuals(x))[[2]])), 3)         
  
  # extraction des p-values pour chaque facteur, anova de type II
  for(i in c(1:4)){
    resum[8:10, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[3]][1:3])), 3)
    resum[8:10, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=2)[[4]][1:3])), 3)
  }
  for(i in c(3:4)){
    resum[11, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[3]][4])), 3)
    resum[11, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=2)[[4]][4])), 3)
    resum[12, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[3]][5])), 3)
    resum[12, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=2)[[4]][5])), 3)
  }
  for(i in c(1:2)){
    resum[12, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[3]][4])), 3)
    resum[12, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=2)[[4]][4])), 3)
  }
  
  # extraction des p-values pour chaque facteur, anova de type III
  for(i in c(1:4)){
    resum[13:15, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[3]][1:3])), 3)
    resum[13:15, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=3)[[4]][1:3])), 3)
  }
  for(i in c(3:4)){
    resum[16, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[3]][4])), 3)
    resum[16, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=3)[[4]][4])), 3)
    resum[17, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[3]][5])), 3)
    resum[17, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=3)[[4]][5])), 3)
  }
  for(i in c(1:2)){
    resum[17, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[3]][4])), 3)
    resum[17, i+4] <- round(unlist(lapply(y[i+4], function(x) Anova(x, type=3)[[4]][4])), 3)
  }
  
  # R?gles de d?cision du mod?le le plus adapt?
  if((resum[11, 3] < 0.1) && (resum[1, 3] > seuil_prcRE)){sel <- resum[, 3:4]
  }else{
    if((resum[11, 7] < 0.1) && (resum[1, 3] < seuil_prcRE)){sel <- resum[, 7:8]
    }else{
      if((resum[1, 1] > seuil_prcRE)){sel <- resum[, 1:2]
      }else{sel <- resum[, 5:6]
      }}}
  
  z <- list(mod1A = summary(mod1A)$call, mod2A = summary(mod2A)$call, mod3A = summary(mod3A)$call,
            mod4A = summary(mod3A)$call, mod1B = summary(mod1B)$call, mod2B = summary(mod2B)$call, 
            mod3B = summary(mod3B)$call, mod4B = summary(mod4B)$call, 
            resum = resum, select = sel)
  z
}

# A l'?chelle de l'esp?ce
myfun_Ind <- function(mdata, vb, seuil_prcRE){
  # on ?crit les modeles 
  mod5A <- lmer(vb ~ -1+ div*leg + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod5B <- lmer(vb ~ -1+ leg*div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6A <- lmer(vb ~ -1+div*leg + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6B <- lmer(vb ~ -1+ leg*div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7A <- lmer(vb ~ -1+div*leg + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7B <- lmer(vb ~ -1+ leg*div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod8A <- lmer(vb ~ -1+div*leg + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = F) 
  mod8B <- lmer(vb ~ -1+ leg*div + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = F) 
  mod9A <- lm(vb ~ -1+div*leg + serie + PAR, data= mdata, na.action = na.omit)
  mod9B <- lm(vb ~ -1+ leg*div + serie + PAR, data= mdata, na.action = na.omit)
  mod10A <- lm(vb ~ -1+div*leg + serie, data= mdata, na.action = na.omit)
  mod10B <- lm(vb ~ -1+ leg*div + serie, data= mdata, na.action = na.omit)
  
  # cr?ation d'une df pour collecter tous les param?tres des mod?les 
  # et des tests stats
  resum <- matrix(NA, ncol = 12, nrow = 23, 
                  dimnames = list(c("position (%)", "numPL(%)", "residual (%)", "AIC", "BIC", "p-val", "R2", "Shapiro", 
                                    "Type I: fac1", "Type I: fac2", "Type I: serie", "Type I: PAR", "Type I: inter",
                                    "Type II: fac1", "Type II: fac2", "Type II: serie", "Type II: PAR","Type II: inter", 
                                    "Type III: fac1", "Type III: fac2", "Type III: serie", "Type III: PAR","Type III: inter"), 
                                  c("mod5A", "mod5B", "mod6A", "mod6B", "mod7A", "mod7B", "mod8A", "mod8B", "mod9A", "mod9B", "mod10A", "mod10B")))
  
  # insertion des AIC et BIC
  resum[4,] <- round(t(AIC(mod5A, mod5B, mod6A, mod6B, mod7A, mod7B, mod8A, mod8B, mod9A, mod9B, mod10A, mod10B)[,2]), 0)
  resum[5,] <- round(t(BIC(mod5A, mod5B, mod6A, mod6B, mod7A, mod7B, mod8A, mod8B, mod9A, mod9B, mod10A, mod10B)[,2]), 0)
  
  # pour les lmer
  # liste contenant tous les lmer
  lmerlist <- list(mod5A, mod5B, mod6A, mod6B, mod7A, mod7B, mod8A, mod8B) 
  # calcul du % total de variance expliqu? par les effets al?atoires
  allRE_prc <- lapply(lmerlist, function(x) 
    as.data.frame(VarCorr(x))$vcov/sum(as.data.frame(VarCorr(x))$vcov))
  # insertion allRE_prc des mod?les ? deux RE (numPL, pos + residuelle)
  for (i in c(1, 2, 5,6)){
    resum[1:3,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  # insertion allRE_prc des mod?les ? 1 RE (numPl + residuelle)
  for (i in c(3, 4, 7, 8)){
    resum[2:3,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  
  # pour les lm
  # liste contenant les lm
  lmlist <- list(mod9A, mod9B, mod10A, mod10B)
  # insertion des p-val et des R? des lm
  for (i in c(1:4)){
    resum[6, i+8] <- round(lapply(lmlist, function(x) anova(x)$'Pr(>F)'[1])[[i]], 3)
    resum[7, i+8] <- round(lapply(lmlist, function(x) round(summary(x)[[9]], 2))[[i]], 2)
  }
  
  # on r??crit les mod?les avec REML = TRUE
  mod5A <- lmer(vb ~ -1+ div*leg + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod5B <- lmer(vb ~ -1+ leg*div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6A <- lmer(vb ~ -1+div*leg + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6B <- lmer(vb ~ -1+ leg*div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7A <- lmer(vb ~ -1+div*leg + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7B <- lmer(vb ~ -1+ leg*div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod8A <- lmer(vb ~ -1+div*leg + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = F) 
  mod8B <- lmer(vb ~ -1+ leg*div + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = F) 
  
  # cr?ation d'une liste avec tous les mod?les (lmer + lm)
  y <- list(mod5A, mod5B, mod6A, mod6B, mod7A, mod7B, mod8A, mod8B, mod9A, mod9B, mod10A, mod10B)
  
  # v?rification de la distribution des r?sidus
  resum[8, ] <- round(unlist(lapply(y, function(x) shapiro.test(residuals(x))[[2]])), 3)         
  
  # extraction des p-values pour chaque facteur, anova de type I
  for(i in c(1:4)){
    resum[c(9:11, 13), i] <- round(unlist(lapply(y[i], function(x) anova(x)[[4]])), 3)
  }
  for(i in c(5:8)){
    resum[9:13, i] <- round(unlist(lapply(y[i], function(x) anova(x)[[4]])), 3)
  }
  for(i in c(9:10)){
    resum[9:13, i] <- round(unlist(lapply(y[i], function(x) anova(x)[[5]][1:5])), 3)
  }
  for(i in c(11:12)){
    resum[c(9:11,13), i] <- round(unlist(lapply(y[i], function(x) anova(x)[[5]][1:4])), 3)
  }
  
  # extraction des p-values pour chaque facteur, anova de type II
  for(i in c(1:4)){
    resum[c(14:16, 18), i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[3]])), 3)
  }
  for(i in c(5:8)){
    resum[14:18, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[3]])), 3)
  }
  for(i in c(9:10)){
    resum[14:18, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[4]][1:5])), 3)
  }
  for(i in c(11:12)){
    resum[c(14:16,18), i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=2)[[4]][1:4])), 3)
  }
  
  # extraction des p-values pour chaque facteur, anova de type III
  for(i in c(1:4)){
    resum[c(19:21, 23), i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[3]])), 3)
  }
  for(i in c(5:8)){
    resum[19:23, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[3]])), 3)
  }
  for(i in c(9:10)){
    resum[19:23, i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[4]][1:5])), 3)
  }
  for(i in c(11:12)){
    resum[c(19:21, 23), i] <- round(unlist(lapply(y[i], function(x) Anova(x, type=3)[[4]][1:4])), 3)
  }
  
  # R?gles de d?cision du mod?le le plus adapt?
  if((resum[17, 5] < 0.1) && (resum[3, 5] < (100-seuil_prcRE))){sel <- resum[, 5:6]
  }else{
    if((resum[17, 7] < 0.1) && (resum[3, 7] < (100-seuil_prcRE))){sel <- resum[, 7:8]
    }else{
      if((resum[17, 9] < 0.1)){sel <- resum[, 9:10]
      }else{
        if((resum[1, 1] > 5) && (resum[3, 1] < (100-seuil_prcRE))){sel <- resum[, 1:2]
        }else{
          if(resum[3, 3] < (100-seuil_prcRE)){sel <- resum[, 3:4]
          }else{sel <- resum[, 11:12]
          }}}}}
  
  z <- list(mod5A = summary(mod5A)$call, mod6A = summary(mod6A)$call, mod7A = summary(mod7A)$call,
            mod8A = summary(mod8A)$call, mod9A = summary(mod9A)$call, mod10A = summary(mod10A)$call, 
            mod5B = summary(mod5B)$call, mod6B = summary(mod6B)$call, mod7B = summary(mod7B)$call,
            mod8A = summary(mod8B)$call, mod9B = summary(mod9B)$call, mod10B = summary(mod10B)$call, 
            resum = resum, select = sel)
  z
}

# A l'?chelle de l'esp?ce pour la l?gumineuse
myfun_Leg <- function(mdata, vb, seuil_prcRE){
  # on ?crit les modeles 
  mod5A <- lmer(vb ~ -1+ div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6A <- lmer(vb ~ -1+div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7A <- lmer(vb ~ -1+div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod8A <- lmer(vb ~ -1+div + serie + PAR  + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod9A <- lm(vb ~ -1+div + serie + PAR, data = mdata, na.action = na.omit)
  mod10A <- lm(vb ~ -1+div + serie, data = mdata, na.action = na.omit)
  
  # cr?ation d'une df pour collecter tous les param?tres des mod?les 
  # et des tests stats
  resum <- matrix(NA, ncol = 6, nrow = 17, 
                  dimnames = list(c("position (%)", "numPL(%)", "residual (%)", "AIC", "BIC", "p-val", "R2", "Shapiro", 
                                    "Type I: div", "Type I: serie", "Type I: PAR",
                                    "Type II: div", "Type II: serie", "Type II: PAR", 
                                    "Type III: div", "Type III: serie", "Type III: PAR"), 
                                  c("mod5A", "mod6A", "mod7A", "mod8A", "mod9A", "mod10A")))
  
  # insertion des AIC et BIC
  resum[4,] <- round(t(AIC(mod5A, mod6A, mod7A, mod8A, mod9A, mod10A)[,2]), 0)
  resum[5,] <- round(t(BIC(mod5A, mod6A, mod7A, mod8A, mod9A, mod10A)[,2]), 0)
  
  # pour les lmer
  # liste contenant tous les lmer
  lmerlist <- list(mod5A, mod6A, mod7A, mod8A) 
  # calcul du % total de variance expliqu? par les effets al?atoires
  allRE_prc <- lapply(lmerlist, function(x) 
    as.data.frame(VarCorr(x))$vcov/sum(as.data.frame(VarCorr(x))$vcov))
  # insertion allRE_prc des mod?les ? deux RE (numPL, pos + residuelle)
  for (i in c(1, 3)){
    resum[1:3,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  # insertion allRE_prc des mod?les ? 1 RE (numPl + residuelle)
  for (i in c(2, 4)){
    resum[2:3,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  
  # pour les lm
  # liste contenant les lm
  lmlist <- list(mod9A, mod10A)
  # insertion des p-val et des R? des lm
  for (i in c(1:2)){
    resum[6, i+4] <- round(lapply(lmlist, function(x) anova(x)$'Pr(>F)'[1])[[i]], 3)
    resum[7, i+4] <- round(lapply(lmlist, function(x) round(summary(x)[[9]], 2))[[i]], 2)
  }
  
  # on r??crit les mod?les avec REML = TRUE
  mod5A <- lmer(vb ~ -1+ div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = T) 
  mod6A <- lmer(vb ~ -1+div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = T) 
  mod7A <- lmer(vb ~ -1+div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = T) 
  mod8A <- lmer(vb ~ -1+div + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = T) 
  
  
  # cr?ation d'une liste avec tous les mod?les (lmer + lm)
  y <- list(mod5A, mod6A, mod7A, mod8A, mod9A, mod10A)
  
  # v?rification de la distribution des r?sidus
  resum[8, ] <- round(unlist(lapply(y, function(x) shapiro.test(residuals(x))[[2]])), 3)         
  
  # extraction des p-values pour chaque facteur, anova de type I
  #resum[9:10, 1] <- round(anova(mod5A)[[5]], 3)
  #resum[9:10, 2] <- round(anova(mod6A)[[5]], 3)
  #resum[9:11, 3] <- round(anova(mod7A)[[5]], 3)
  #resum[9:11, 4] <- round(anova(mod8A)[[5]], 3)
  resum[9:11, 5] <- round(anova(mod9A)[[5]][1:3], 3)
  resum[9:10, 6] <- round(anova(mod10A)[[5]][1:2], 3)
  
  # extraction des p-values pour chaque facteur, anova de type II
  resum[12:13, 1] <- round(Anova(mod5A, type=2)[[3]], 3)
  resum[12:13, 2] <- round(Anova(mod6A, type=2)[[3]], 3)
  resum[12:14, 3] <- round(Anova(mod7A, type=2)[[3]], 3)
  resum[12:14, 4] <- round(Anova(mod8A, type=2)[[3]], 3)
  resum[12:14, 5] <- round(Anova(mod9A, type=2)[[4]][1:3], 3)
  resum[12:13, 6] <- round(Anova(mod10A, type=2)[[4]][1:2], 3)
  
  # extraction des p-values pour chaque facteur, anova de type III
  resum[15:16, 1] <- round(Anova(mod5A, type=3)[[3]], 3)
  resum[15:16, 2] <- round(Anova(mod6A, type=3)[[3]], 3)
  resum[15:17, 3] <- round(Anova(mod7A, type=3)[[3]], 3)
  resum[15:17, 4] <- round(Anova(mod8A, type=3)[[3]], 3)
  resum[15:17, 5] <- round(Anova(mod9A, type=3)[[4]][1:3], 3)
  resum[15:16, 6] <- round(Anova(mod10A, type=3)[[4]][1:2], 3)
  
  # R?gles de d?cision du mod?le le plus adapt?
  if((resum[17, 3] < 0.1) && (resum[3, 3] < (100-seuil_prcRE))){sel <- as.data.frame(resum[, 3])
  }else{
    if((resum[17, 4] < 0.1) && (resum[3, 4] < (100-seuil_prcRE))){sel <- as.data.frame(resum[, 4])
    }else{
      if((resum[17, 5] < 0.1)){sel <- as.data.frame(resum[, 5])
      }else{
        if((resum[1, 1] > 5) && (resum[3, 1] < (100-seuil_prcRE))){sel <- as.data.frame(resum[, 1])
        }else{
          if(resum[3, 2] < (100-seuil_prcRE)){sel <- as.data.frame(resum[, 2])
          }else{sel <- as.data.frame(resum[, 6])
          }}}}}
  
  z <- list(mod5A = summary(mod5A)$call, mod6A = summary(mod6A)$call, mod7A = summary(mod7A)$call,
            mod8A = summary(mod8A)$call, mod9A = summary(mod10A)$call, mod10A = summary(mod10A)$call, 
            resum = resum, select = sel)
  suppressWarnings(z)
  z
}        
#------------------------------------------------------------