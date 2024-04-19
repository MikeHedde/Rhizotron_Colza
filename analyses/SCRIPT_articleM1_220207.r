#####################"

## analyses biomasse plante / Rhizotron


# MZ
# 28/09/2016 : reprend SCRIPT_BIOM_PL_20160928.r au propre
# 6/12/2016 : complete stat effet traitement sur biomass pop
# 26/06/2017 : refait analyse après elimination des individus trop petits (germintaion tardive)

# 11/08/2017 : reprends avec les stas de D. Makowsky
# 16/08/2017 : modif script + création de fonctions par MH

################################################################################

## library
require(ggplot2)
library(lme4)
require(car)
library(plotrix)
library(multcomp)
library(PMCMR)
library(gridExtra)
library(reshape2)
library(Hmisc)
library(SensoMineR)
library(reshape2)

## répertoire de travail
setwd("D:/2015 JENNY M ZWICKE/2.Manipes/1.Rhizotron/1. Analyses Mike")
#setwd("C:/1. Mike/5. RH/4. Post-doctorant/2015 JENNY M ZWICKE/2.Manipes/1.Rhizotron/1. Analyses Mike")
# setwd("F:/5. RH/4. Post-doctorant/2015 JENNY M ZWICKE/Manip/170111_depuis sauvegarde réseau/1. Interactions Vdt - colza")

pdf()

#------------------------------------------------------------
# 0   Fonctions créées pour automatiser le travail de mise en
#     forme des données, de sélection du modèle le plus adapté 
#     et de représentation graphique de sortie des modèles
#------------------------------------------------------------
# Mise en forme des df avec ajout des facteurs serie, id_loca et PAR
      myprep <- function(mydata){
        
          # on ajoute une colonne série
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
          
          # on ajoute le rayonnement en covariable (à modifier quand les valeurs de PAR de la 2e série seront dispo)
          parmoy <- aggregate(par$PAR, by=list(RT=par$RT), mean, na.rm=T)
          parmoy2 <- parmoy
          
          # manipulations pour attribuer une valeur modélisée à chaque position de plante
          RTtemp = seq(39,1,-2)
          parmoy2$RT <- parmoy$RT + RTtemp
          parmoy3 <- rbind(parmoy, parmoy2, parmoy2[which(parmoy2$RT == c(30:31)),])
          parmoy3$RT[41:42] <- c(42, 41)
          
          mydata <- merge(mydata, parmoy3, by=c("RT"))
          colnames(mydata)[ncol(mydata)] <- "PAR"
          
          mydata
}

# Création d'une fonction pour retenir le modèle le plus adapté 
# parmi les modèles définis (et en intervertissant l'ordre leg / div)  
    # Les arguments sont :
          # mdata => df contenant l'ensemble des valeurs pour le calcul des modèles
          # vb => variable dépendante
          # seuil_prcRE => valeur minimale souhaitée pour le % de 
          #                variance expliquée par le random effect

    # rq : contraste par défaut = options(contrasts = c("contr.treatment", "contr.treatment"))
        options(contrasts =c ("contr.sum","contr.poly"))  # A l'échelle de la communauté

    # A l'échelle de la communauté        
        myfun_Com <- function(mdata, vb, seuil_prcRE){
          
          # Différents modeles complets avec et sans covariable pour sélectionner le modèle le plus adapté
          mod1A <- lmer(vb ~ -1 + div*leg + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod2A <- lmer(vb ~ -1 + div*leg + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod3A <- lm(vb ~ -1 + div*leg + serie, data = mdata, na.action = na.omit) 
          mod4A <- lm(vb ~ -1 + div*leg + serie + PAR, data = mdata, na.action = na.omit)
          mod1B <- lmer(vb ~ -1 + leg*div + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod2B <- lmer(vb ~ -1 + leg*div + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod3B <- lm(vb ~ -1 + leg*div + serie, data = mdata, na.action = na.omit) 
          mod4B <- lm(vb ~ -1 + leg*div + serie + PAR, data = mdata, na.action = na.omit) 
          
          # Les paramètres à comparer sont collectés dans une table
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
          
          # les lmer sont recalculés ont posant REML = T
          mod1A <- lmer(vb ~ -1 + div*leg + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod2A <- lmer(vb ~ -1 + div*leg + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod1B <- lmer(vb ~ -1 + leg*div + serie + (1|position), data = mdata, na.action = na.omit, REML = F) 
          mod2B <- lmer(vb ~ -1 + leg*div + serie + PAR + (1|position), data = mdata, na.action = na.omit, REML = F) 
          
          y <- list(mod1A, mod1B, mod2A, mod2B, mod3A, mod3B, mod4A, mod4B)
          
          # vérification de la distribution des résidus
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
          
          # Règles de décision du modèle le plus adapté
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

    # A l'échelle de l'espèce
        myfun_Ind <- function(mdata, vb, seuil_prcRE){
  # on écrit les modeles 
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
  
  # création d'une df pour collecter tous les paramètres des modèles 
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
  # calcul du % total de variance expliqué par les effets aléatoires
  allRE_prc <- lapply(lmerlist, function(x) 
    as.data.frame(VarCorr(x))$vcov/sum(as.data.frame(VarCorr(x))$vcov))
  # insertion allRE_prc des modèles à deux RE (numPL, pos + residuelle)
  for (i in c(1, 2, 5,6)){
    resum[1:3,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  # insertion allRE_prc des modèles à 1 RE (numPl + residuelle)
  for (i in c(3, 4, 7, 8)){
    resum[2:3,i] <- round(allRE_prc[[i]] * 100, 0)
  }
  
  # pour les lm
  # liste contenant les lm
  lmlist <- list(mod9A, mod9B, mod10A, mod10B)
  # insertion des p-val et des R² des lm
  for (i in c(1:4)){
    resum[6, i+8] <- round(lapply(lmlist, function(x) anova(x)$'Pr(>F)'[1])[[i]], 3)
    resum[7, i+8] <- round(lapply(lmlist, function(x) round(summary(x)[[9]], 2))[[i]], 2)
  }
  
  # on réécrit les modèles avec REML = TRUE
  mod5A <- lmer(vb ~ -1+ div*leg + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod5B <- lmer(vb ~ -1+ leg*div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6A <- lmer(vb ~ -1+div*leg + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod6B <- lmer(vb ~ -1+ leg*div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7A <- lmer(vb ~ -1+div*leg + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod7B <- lmer(vb ~ -1+ leg*div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
  mod8A <- lmer(vb ~ -1+div*leg + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = F) 
  mod8B <- lmer(vb ~ -1+ leg*div + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = F) 
  
  # création d'une liste avec tous les modèles (lmer + lm)
  y <- list(mod5A, mod5B, mod6A, mod6B, mod7A, mod7B, mod8A, mod8B, mod9A, mod9B, mod10A, mod10B)
  
  # vérification de la distribution des résidus
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
  
  # Règles de décision du modèle le plus adapté
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

    # A l'échelle de l'espèce pour la légumineuse
        myfun_Leg <- function(mdata, vb, seuil_prcRE){
          # on écrit les modeles 
          mod5A <- lmer(vb ~ -1+ div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
          mod6A <- lmer(vb ~ -1+div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
          mod7A <- lmer(vb ~ -1+div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
          mod8A <- lmer(vb ~ -1+div + serie + PAR  + (1|numPL), data = mdata, na.action = na.omit, REML = F) 
          mod9A <- lm(vb ~ -1+div + serie + PAR, data = mdata, na.action = na.omit)
          mod10A <- lm(vb ~ -1+div + serie, data = mdata, na.action = na.omit)
          
          # création d'une df pour collecter tous les paramètres des modèles 
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
          # calcul du % total de variance expliqué par les effets aléatoires
          allRE_prc <- lapply(lmerlist, function(x) 
            as.data.frame(VarCorr(x))$vcov/sum(as.data.frame(VarCorr(x))$vcov))
          # insertion allRE_prc des modèles à deux RE (numPL, pos + residuelle)
          for (i in c(1, 3)){
            resum[1:3,i] <- round(allRE_prc[[i]] * 100, 0)
          }
          # insertion allRE_prc des modèles à 1 RE (numPl + residuelle)
          for (i in c(2, 4)){
            resum[2:3,i] <- round(allRE_prc[[i]] * 100, 0)
          }
          
          # pour les lm
          # liste contenant les lm
          lmlist <- list(mod9A, mod10A)
          # insertion des p-val et des R² des lm
          for (i in c(1:2)){
            resum[6, i+4] <- round(lapply(lmlist, function(x) anova(x)$'Pr(>F)'[1])[[i]], 3)
            resum[7, i+4] <- round(lapply(lmlist, function(x) round(summary(x)[[9]], 2))[[i]], 2)
          }
          
          # on réécrit les modèles avec REML = TRUE
          mod5A <- lmer(vb ~ -1+ div + serie + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = T) 
          mod6A <- lmer(vb ~ -1+div + serie + (1|numPL), data = mdata, na.action = na.omit, REML = T) 
          mod7A <- lmer(vb ~ -1+div + serie + PAR  + (1|position) + (1|numPL), data = mdata, na.action = na.omit, REML = T) 
          mod8A <- lmer(vb ~ -1+div + serie + PAR  + (1|numPL), data= mdata, na.action = na.omit, REML = T) 

          
          # création d'une liste avec tous les modèles (lmer + lm)
          y <- list(mod5A, mod6A, mod7A, mod8A, mod9A, mod10A)
          
          # vérification de la distribution des résidus
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
          
          # Règles de décision du modèle le plus adapté
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
        
# Création d'une fonction pour représenter les différences entre modalités
    # Les arguments sont :  
          # mdata, df contenant les données de calcul des modèles
          # mod => modèle selectionné (il faut utiliser le $modxx généré 
          #       par la fonction myfun_Com)
          # ylbl => ylabel des graphiques
          # plot_PAR => T/F pour représenter un effet du PAR    
    # Les outputs sont : 
          # moyFRic => moyenne et se de vb par modalités du facteur 'div'
          # moyLeg => moyenne et se de vb par modalités du facteur 'leg'
          # moyInteraction => moyenne et se de vb par modalités croisées des facteurs 'div et 'leg'
          # pDIV => graph des moy (+/- se) de vb en fonction des modalités du facteur 'div'
          # pDIV_PAR => droite de réponse de vb en fonction du PAR (modulo les niveaux de 'div')
          # pDIV_INTER => graph des moy (+/- se) de vb en fonction des modalités croisées des facteurs 'div et 'leg'
          # pLEG => graph des moy (+/- se) de vb en fonction des modalités du facteur 'leg' 
          # pLEG_PAR => droite de réponse de vb en fonction du PAR (modulo les niveaux de 'leg')
        
      # A l'échelle de la communauté vdt  
      myplot_Com_vdt <- function(mdata, vb, mod, ylabl, plot_PAR){
          mod <- eval(mod)
          newdata <- data.frame(mdata) 
          newdata$predict <- predict(mod)
          
          # Effet FRic 
          moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
          seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
          moyF$se <- seF$x
          
          DIV <- ggplot(data = moyF, aes(x = div,y=x, fill=div)) + 
            geom_bar(position = position_dodge(width = 0.8), stat="identity")+
            geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
            ylab(ylabl) +    
            xlab("") +
            scale_fill_manual(values = c("#ffeda0", "#feb24c","#f03b20"), 
                              name = "",
                              labels = c("FR1", "FR2", "FR3"))+
            theme_classic() +
            theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
                  text =  element_text(face = "plain",
                                       color = "black", 
                                       hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                       margin = margin(), debug = FALSE))+
            theme(axis.text = element_text(color="black", size=14)) +
            theme(axis.title = element_text(color="black", size=14))+
            scale_x_discrete(labels=c("Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
          
          DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
            geom_point()+
            geom_smooth(method = "lm")+
            ylab(ylabl)+
            theme(text =  element_text(face = "plain",
                                       color = "black", 
                                       hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                       margin = margin(), debug = FALSE))+
            xlab("PAR (µmol/m²/s)")
          
          # Effet leg 
          moyL <- aggregate(newdata$predict, by = list(leg = newdata$leg), mean, na.rm = T)
          seL <- aggregate(newdata$predict, by = list(leg = newdata$leg), std.error, na.rm = T)
          moyL$se <- seL$x
          LEG <- ggplot(moyL, aes(x = leg, y = x)) + 
            geom_pointrange(aes(ymin = x-1.96*se, ymax = x+1.96*se))+
            theme(text =  element_text(face = "plain",
                                       color = "black", 
                                       hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                       margin = margin(), debug = FALSE))+
            ylab("")
          
          LEG_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = leg, col = leg)) + 
            geom_point()+
            geom_smooth(method = "lm")+
            theme(text =  element_text(face = "plain",
                                       color = "black", 
                                       hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                       margin = margin(), debug = FALSE))+
            xlab("PAR (µmol/m²/s)")
          
          # Effet interaction
          moyFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), mean, na.rm=T)
          seFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), std.error, na.rm=T)
          moyFi$se <- seFi$x
          
          INTER <- ggplot(data = moyFi, aes(x = leg,y=x, fill=div)) + 
            geom_bar(position = position_dodge(width = 0.8), stat="identity")+
            geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
            ylab(ylabl) +    
            xlab("") +
            scale_fill_manual(values = c("#ffeda0", "#feb24c","#f03b20"), 
                              name = "",
                              labels = c("FR1", "FR2", "FR3"))+
            theme_classic() +
            theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
                  text =  element_text(face = "plain",
                                       color = "black", 
                                       hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                       margin = margin(), debug = FALSE))+
            theme(axis.text = element_text(color="black", size=14)) +
            theme(axis.title = element_text(color="black", size=14))+
            scale_x_discrete(labels=c("Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))        
          
          if (plot_PAR == FALSE) {A <- grid.arrange(DIV, LEG, INTER, ncol = 3)
          }else{A <- grid.arrange(DIV, LEG, INTER, DIV_PAR, LEG_PAR, ncol = 3, nrow = 2)
          }
          z <- list(moyFRic = moyF, moyLeg = moyL, moyInter = moyFi, pDIV = DIV, pINTER = INTER, 
                    pLEG = LEG, pDIV_PAR = DIV_PAR, pLEG_PAR = LEG_PAR)
          print(A)
          z
          
        }
        
       # A l'échelle de la communauté plante
      myplot_Com <- function(mdata, vb, mod, ylabl, plot_PAR){
        mod <- eval(mod)
        newdata <- data.frame(mdata) 
        newdata$predict <- predict(mod)
        
        # Effet FRic 
        moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
        seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
        moyF$se <- seF$x
        DIV <- ggplot(data = moyF, aes(x = div, y = x, shape = div)) + 
          geom_point(size = 4)+
          geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
          ylab(ylabl) +    
          xlab("Earthworm FRic") +
          theme(panel.background = element_rect(fill="white", colour = "black"),
                text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                    margin = margin(), debug = FALSE))+
          geom_vline(xintercept = 1.5)
        
            
        DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
          geom_point()+
          geom_smooth(method = "lm")+
          ylab(ylabl)+
          theme(text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                     margin = margin(), debug = FALSE))+
          xlab("PAR (µmol/m²/s)")
        
        # Effet leg 
        moyL <- aggregate(newdata$predict, by = list(leg = newdata$leg), mean, na.rm = T)
        seL <- aggregate(newdata$predict, by = list(leg = newdata$leg), std.error, na.rm = T)
        moyL$se <- seL$x
        LEG <- ggplot(moyL, aes(x = leg, y = x)) + 
          geom_pointrange(aes(ymin = x-1.96*se, ymax = x+1.96*se))+
          theme(text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                     margin = margin(), debug = FALSE))+
          ylab("")
        
        LEG_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = leg, col = leg)) + 
          geom_point()+
          geom_smooth(method = "lm")+
          theme(text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                     margin = margin(), debug = FALSE))+
          xlab("PAR (µmol/m²/s)")
        
        # Effet interaction
        moyFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), mean, na.rm=T)
        seFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), std.error, na.rm=T)
        moyFi$se <- seFi$x
        INTER <- ggplot(data = moyFi, aes(x = div, y = x)) + 
          geom_point()+
          geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
          ylab(ylabl)+
          facet_grid(leg~.)+
          theme(panel.background = element_rect(fill="white", colour = "black"),
                text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                     margin = margin(), debug = FALSE))
        
        
        if (plot_PAR == FALSE) {A <- grid.arrange(DIV, LEG, INTER, ncol = 3)
        }else{A <- grid.arrange(DIV, LEG, INTER, DIV_PAR, LEG_PAR, ncol = 3, nrow = 2)
        }
        z <- list(moyFRic = moyF, moyLeg = moyL, moyInter = moyFi, pDIV = DIV, pINTER = INTER, 
                  pLEG = LEG, pDIV_PAR = DIV_PAR, pLEG_PAR = LEG_PAR)
        print(A)
        z
        
}
    
      # A l'échelle de l'espèce
      myplot_Ind <- function(mdata, vb, mod, ylabl, plot_PAR){
  mod <- eval(mod)
  newdata <- data.frame(mdata) 
  newdata$predict <- predict(mod)
  
  # Effet FRic 
  moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
  seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
  moyF$se <- seF$x
  DIV <- ggplot(data = moyF, aes(x = div, y=x, fill=div)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("Ew0", "FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
  
    
  
  DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    ylab(ylabl)+
    xlab("PAR (µmol/m²/s)")+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  moyFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), mean, na.rm=T)
  seFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), std.error, na.rm=T)
  moyFi$se <- seFi$x
  
  DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("Ew0", "FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
  
  
  # Effet leg 
  moyL <- aggregate(newdata$predict, by = list(leg = newdata$leg), mean, na.rm = T)
  seL <- aggregate(newdata$predict, by = list(leg = newdata$leg), std.error, na.rm = T)
  moyL$se <- seL$x
  LEG <- ggplot(moyL, aes(x = leg, y = x)) + 
    geom_pointrange(aes(ymin = x-1.96*se, ymax = x+1.96*se))+
    ylab(ylabl)+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  LEG_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = leg, col = leg)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    xlab("PAR (µmol/m²/s)")+
    ylab(ylabl)+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  if (plot_PAR == FALSE) {A <- grid.arrange(DIV, LEG, DIV_INTER, ncol = 2, nrow = 2)
  }else{A <- grid.arrange(DIV, DIV_PAR, DIV_INTER, LEG, LEG_PAR, ncol = 3, nrow = 2)
  }
  print(A)
  z <- list(moyFRic = moyF, moyInteraction = moyFi, moyLeg = moyL, 
            pDIV = DIV, pDIV_PAR = DIV_PAR, pDIV_INTER = DIV_INTER, pLEG = LEG, pLEG_PAR = LEG_PAR)
  z
}

      # A l'échelle de l'espèce pour le trèfle      
      myplot_Leg <- function(mdata, vb, mod, ylabl, plot_PAR){
        mod <- eval(mod)
        newdata <- data.frame(mdata) 
        newdata$predict <- predict(mod)
        
        # Effet FRic 
        moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
        seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
        moyF$se <- seF$x
        
        DIV <- ggplot(data = moyF, aes(x = div, y = x, fill = div)) + 
          geom_bar(position = position_dodge(width = 0.8), stat="identity")+
          geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
          ylab(ylabl) +    
          xlab("") +
          scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                            name = "",
                            labels = c("Ew0","FR1", "FR2", "FR3"))+
          theme_classic() +
          theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
                text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                     margin = margin(), debug = FALSE))+
          theme(axis.text = element_text(color="black", size=14)) +
          theme(axis.title = element_text(color="black", size=14))+
          scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
        
        
        DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
          geom_point()+
          geom_smooth(method = "lm")+
          ylab(ylabl)+
          xlab("PAR (µmol/m²/s)")+
          theme(panel.background = element_rect(fill="white", colour = "black"),
                text =  element_text(face = "plain",
                                     color = "black", 
                                     hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                     margin = margin(), debug = FALSE))
        
        if (plot_PAR == FALSE) {A <- DIV
        }else{A <- grid.arrange(DIV, DIV_PAR, ncol = 2, nrow = 1)
        }
        print(A)
        z <- list(moyFRic = moyF, pDIV = DIV, pDIV_PAR = DIV_PAR)
        z
      }
    
      
    #------------------------------------------------------------
    # 1.0 Chargement des données 
    #------------------------------------------------------------
    pl <- read.csv("DATA_RT_FW_DW_all_data.csv", h = T, sep = ";", dec = ".", na.strings = "")
    par <- read.csv("DOC_RT_rayonnement_plante.csv", sep = ";", dec=".", h=T)
    # filter_bn <- read.csv("DOC_filtre_BN.csv", sep=";", dec=".", header=T, na.strings = NA)
    # filter_lp <- read.csv("DOC_filtre_LP.csv", sep=";", dec=".", header=T, na.strings = NA)
    vdt <- read.csv("DOC_RT_EW_FW.csv", sep=";", dec=".", header=T, na.strings = NA)
    gal <- read.csv("DATA_RT_GAL_170828.csv", header = TRUE, sep = ";", dec=".", na.strings=NA)
    Ncont <- read.csv("DATA_RT_N_160920.csv", h = T, sep = ";", dec = ".", na.strings = "")
    Pcont <- read.csv("DATA_RT_P_180517.csv", h = T, sep = ";", dec = ".", na.strings = "")
    enz <- read.csv("DATA_RT_ENZ_180817.csv", h = T, sep = ";", dec = ".", na.strings = "")
    lgrac <- read.csv("DATA_RT_LG_RAC_170828.csv", h = T, sep = ";", dec = ".", na.strings = "") 
    nod <- read.csv("DATA_RT_Tri_alexa_NOD.csv", sep=";", dec=".")
    anion <- read.csv("DATA_RT_anion.csv", sep=";", dec=".")
    BnGerm <- read.csv("DATA_RT_Bra_napus_GermDay.csv", sep=";", dec=".")
    
    #------------------------------------------------------------
    # 1.1 Calcul des biomasses végétales sèches
    #-----------------------------------------------------------
      # Masses aérienne et racinaires individuelle par strate
          biom_Ind <- aggregate(DW ~ div + leg + id_loca + id_taxo + distance + strata + RT, sum, data = pl, na.rm = T)
          
      # Masse racinaire individuelle totale 
          biom_R_Ind <- aggregate(DW ~ div + leg + id_loca + id_taxo + distance + strata + RT, sum, data = pl[pl$strata > 0,], na.rm=T)
          colnames(biom_R_Ind)[8] <- "DW_R"
          
      # Masse totale Communauté (aérien + racines)
          biom_T_Com <- aggregate(DW ~ div + leg + RT, sum, data = biom_Ind, na.rm = T)
          colnames(biom_T_Com)[4] <- "DW_tot"
          
      # Masse aérienne Communauté
          biom_A_Com <- aggregate(biom_Ind$DW, by=list(RT=biom_Ind$RT, div=biom_Ind$div, leg=biom_Ind$leg, strata=biom_Ind$strata), sum, na.rm=T)
          biom_A_Com <- biom_A_Com[biom_A_Com$strata == 0, -4]
          colnames(biom_A_Com)[4] <- "DW_A"
          
      # Masse racinaire Communauté
          biom_R_Com <- aggregate(biom_R_Ind$DW, by=list(RT=biom_R_Ind$RT, div=biom_R_Ind$div, leg=biom_R_Ind$leg), sum, na.rm=T)
          colnames(biom_R_Com)[4] <- "DW_R"
    
      # Création d'une df avec les valeurs tot, aeriennes et racinaires
          biom_Com <- merge(biom_T_Com, biom_A_Com, by=c("RT", "div", "leg"), sort=T)
          biom_Com <- merge(biom_Com, biom_R_Com, by=c("RT", "div", "leg"), sort=T)
    
    #-----------------------------------------------------------
    # 1.2 Mise en forme des df pour les traitements statistiques 
    #     (1) à l'échelle de la communauté (ajout facteurs série,
    #     position, PAR) pour les plantes, les activités enzymatiques 
    #     et les vdt
    #     (2) à l'échelle de l'individu pour les valeurs de biomasse 
    #     des espèces végétales et de [N] des organes de colza, 
    #     ray-grass et trèfle
    #------------------------------------------------------------
      # Communautés   
          # mise en forme pour les plantes
          #tempoCom <- merge (biom_Com, lgrac, by=c("RT", "div", "leg"), sort=T)
          mybiom_Com <- myprep(biom_Com) 
          
          # mise en forme pour les activités enzymatiques
          myenz <- myprep(enz) 
          
          # mise en forme pour les vdt
          FW_i <- aggregate(vdt$FW[vdt$jour == "0"], list(div = vdt$div[vdt$jour == "0"], leg = vdt$leg[vdt$jour == "0"], RT = vdt$RT[vdt$jour == "0"]), mean, na.omit = T)
          FW_f <- aggregate(vdt$FW[vdt$jour != "0"], list(div = vdt$div[vdt$jour != "0"], leg = vdt$leg[vdt$jour != "0"], RT = vdt$RT[vdt$jour != "0"]), mean, na.omit = T)
          FW_f$FWi <- FW_i$x
          colnames(FW_f)[4] <- "FWf"
          FW_f$massLoss <- (FW_f$FWi - FW_f$FWf)/FW_f$FWi
          # on ajoute les longueurs de galeries dans la même table
          tempovdt <- merge(FW_f, gal,by=c("RT", "div", "leg"), sort=T )
          myvdt <- myprep(tempovdt) 
      
       # Individus 
          # mise en forme pour les lombrics
              myvdt_esp <- myprep(vdt)
              myvdt_esp_i <- aggregate(FW ~ PAR + div+leg+id_taxo+serie+position+RT, sum, data= myvdt_esp[myvdt_esp$jour == "0",])
              myvdt_esp_i$FWf <- aggregate(FW ~ PAR+div+leg+id_taxo+serie+position+RT, sum, data= myvdt_esp[myvdt_esp$jour != "0",])$FW
              myvdt_esp_i$massLoss <- (myvdt_esp_i$FW - myvdt_esp_i$FWf)/myvdt_esp_i$FW
              myvdt_esp <- myvdt_esp_i
          
          # Petite fonction pour ajouter le facteur id_loca
          myloca <- function(X){
            # on introduit une colonne loc
            X$numPL <- NULL
            
            # pour remplacer le très long code (avec risque d'erreur de saisie) que tu avais écrit  
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
          
          # Masses végétales sèches 
          mybiom_Ind_temp <- myprep(biom_Ind)
          mybiom_Ind <- myloca(mybiom_Ind_temp)
          
          # Date de germination du Colza modélisée           
          myBnGerm <- myprep(BnGerm)
          
          # Concentration en N dans les organes des plantes   
          n_temp <- myprep(Ncont)
          n <- myloca(n_temp)
          n$organes <- NULL
          n$organes[n$org=="Leaf3"|n$org=="Leaf2"|n$org=="Leaf4"|n$org=="Leaf5"|n$org=="Leaf7"|n$org=="Leaf1"|n$org=="Limb"]="Limb"
          n$organes[n$org=="Fine roots"]="Fine roots"
          n$organes[n$org=="Taproot"]="Taproot"
          n$organes[n$org=="Other Leaves"|n$org=="Other leaves"|n$org=="Aerial"]="Leaves"
          n$organes[n$org=="Stem"]="Stem"
          n$organes <- as.factor(n$organes)
          myNcont <- dcast(n, RT+leg+div+id_taxo + id_loca + numPL + distance + serie + position + PAR~organes, value.var = "N")
          
          # Concentration en P dans les organes des plantes   
          p_temp <- myprep(Pcont)
          p <- myloca(p_temp)
          myPcont <- dcast(p, RT+leg+div+id_taxo + numPL + distance + serie + position + PAR~org, value.var = "Pcont", fun.aggregate = mean)
          
          # Nodulations des trefles
          Ta_temp <- myprep(nod)
          ta <- myloca(Ta_temp)
          # on ajoute les longueurs de galeries dans la même table
          Ta_Root <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata != 0, ]
          myTa <- merge(ta, Ta_Root ,by=c("RT", "div", "leg", "id_taxo", "distance","strata","serie","position","PAR","numPL"), sort=T )
          
        #-----------------------------------------------------------
        # 2.0 Effet des traitements sur les lombrics  
        #-----------------------------------------------------------
                #############################################        
                ###        Modifications de masse         ###
                #############################################
                # A l'échelle communautaire
                      vdt_massLoss_T <- myfun_Com(mdata = myvdt, vb = myvdt$massLoss, seuil_prcRE = 30)
                      vdt_massLoss_T$select
                      pvdt_massLoss_T <- myplot_Com(mdata = myvdt, vb = -myvdt$massLoss, mod = vdt_massLoss_T$mod3A, 
                                 ylabl = "Earthworm mass change (%FW mg)", plot_PAR = T)
                      pvdt_massLoss_T$moyFRic # pour voir les moyennes estimées
                      
                  # A l'échelle spécifique
                      # A. icterica
                      myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_icter",]
                      vdt_massLoss_Ai <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
                      vdt_massLoss_Ai$select
                      pvdt_massLoss_Ai <- myplot_Com(mdata = myvdt_esp_temp, vb = -myvdt_esp_temp$massLoss, mod = vdt_massLoss_Ai$mod3A, 
                                                    ylabl = "A. icterica mass change (%FW mg)", plot_PAR = F)
                      
                      # A. caliginosa
                      myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_calig",]
                      vdt_massLoss_Ac <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
                      vdt_massLoss_Ac$select
                      pvdt_massLoss_Ac <- myplot_Com(mdata = myvdt_esp_temp, vb = -myvdt_esp_temp$massLoss, mod = vdt_massLoss_Ac$mod3A, 
                                                     ylabl = "A. caliginosa mass change (%FW mg)", plot_PAR = T)
                      
                      # L. terrestris
                      myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Lum_terre",]
                      vdt_massLoss_Lt <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 24)
                      vdt_massLoss_Lt$select
                      pvdt_massLoss_Lt <- myplot_Com_vdt(mdata = myvdt_esp_temp, vb = -myvdt_esp_temp$massLoss, mod = vdt_massLoss_Lt$mod3B, 
                                                     ylabl = "L. terrestris mass change (%FW mg)", plot_PAR = F)
                      
                      # A. longa
                      myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_longa",]
                      vdt_massLoss_Al <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
                         # pas possible car présent dans un seul traitement  --> KW ??
                      
                      # A. giardi
                      myvdt_esp_temp <- myvdt_esp[myvdt_esp$id_taxo == "Apo_giard",]
                      vdt_massLoss_Ag <- myfun_Com(mdata = myvdt_esp_temp, vb = myvdt_esp_temp$massLoss, seuil_prcRE = 30)
                      # pas possible car présent dans un seul traitement --> KW ??
                      
                      
                      # représentation de toutes les esp simultanément /!!\ (à améliorer)
                      ggplot(myvdt_esp, aes(x = id_taxo, y = -massLoss))+
                        geom_boxplot()+
                        facet_grid(div~leg)
                      
                      #############################################        
                      ###    Topologie du réseau de galeries    ###
                      #############################################   
                      # 
                      # longueur total des galeries
                          vdt_gal_T <- myfun_Com(mdata = myvdt, vb = myvdt$long, seuil_prcRE = 24)
                          vdt_gal_T$select
                          pvdt_gal_T <- myplot_Com_vdt(mdata = myvdt, vb = myvdt$long, mod = vdt_gal_T$mod3A, 
                                                    ylabl = "Total burrow length (m)", plot_PAR = F)
                          pvdt_gal_T$moyFRic # pour voir les moyennes estimées
                          pvdt_gal_T$pDIV
                          
                    # longueur total des galeries
                          vdt_connec_T <- myfun_Com(mdata = myvdt, vb = myvdt$connec , seuil_prcRE = 24)
                          vdt_connec_T$select
                          pvdt_connec_T <- myplot_Com_vdt(mdata = myvdt, vb = myvdt$connec, mod = vdt_connec_T$mod3A, 
                                                       ylabl = "Burrow connectance", plot_PAR = T)
                          pvdt_connec_T$moyFRic # pour voir les moyennes estimées
                          pvdt_connec_T$pDIV
                          
                   # Specific Burrow length
                          myvdt$SBL_f <- myvdt$long/myvdt$FWf
                          myvdt$SBL_mean <- myvdt$long/((myvdt$FWf+myvdt$FWi)/2)
                          vdt_SBL <- myfun_Com(mdata = myvdt, vb = myvdt$SBL_f , seuil_prcRE = 24)
                          vdt_SBL$select
                          pvdt_SBL <- myplot_Com_vdt(mdata = myvdt, vb = myvdt$SBL_f, mod = vdt_connec_T$mod3B, 
                                                          ylabl = "Specific burrow length (m g-1)", plot_PAR = F)
                          pvdt_SBL$moyFRic # pour voir les moyennes estimées
                          pvdt_SBL$pDIV
        #-----------------------------------------------------------
        # 2.1 Effet des traitements sur les enzymes du sol
        #-----------------------------------------------------------
                #############################################        
                ###         Activités enzymatiques        ###
                #############################################
            # strate1
                myenz_s1 <- myenz[myenz$strate == "1",]
                
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
                  myenz_s4_glu$select # pas de transf. adéquate
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
                    xlim(c(-90, 0))+
                    ylab ("Enzymatic activity (AU)")+
                    xlab("Soil depth\n(cm)")+
                    scale_linetype_discrete(guide=FALSE)+
                    theme(legend.justification=c(0,0), legend.position=c(1,1))+
                    theme(panel.background = element_rect(fill="white", colour = "black"))+
                    geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")+
                    #scale_shape_discrete(guide = FALSE)+
                    facet_grid(~enz, scales = "free_x")
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
                  myenz_s4_CN$select # pas de transf. adéquate
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
                  myanion_s4_NO3$select # pas de transf. adéquate
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
                  
                  
        #-----------------------------------------------------------
        # 2.3 Effet des traitements sur la communauté végétale 
        #-----------------------------------------------------------
                  #############################################        
                  ###            Biomasses sèches           ###
                  #############################################   
            # Comparaison des modèles testés et représentations des différences
            MC_DW_T <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_tot, seuil_prcRE = 30)
                MC_DW_T$select
                myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_tot, mod = MC_DW_T$mod4A, 
                       ylabl = "Total biomass of plant community (mg)", plot_PAR = T)
            
            MC_DW_R <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_R, seuil_prcRE = 30)
                MC_DW_R$select
                myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_R, mod = MC_DW_R$mod4A, 
                       ylabl = "Root biomass of plant community (mg)", plot_PAR = T)
        
            MC_DW_A <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_A, seuil_prcRE = 30)
                MC_DW_A$select
                myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_A, mod = MC_DW_A$mod4A, 
                           ylabl = "Shoot biomass of plant community (mg)", plot_PAR = T)
               
                      #############################################        
                      ###               Biomasse des racines    ###
                      ############################################# 
                      
                COM_LG_R <- myfun_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_tot, seuil_prcRE = 30)
                MC_DW_T$select
                pCOM_LG_R <- myplot_Com(mdata = mybiom_Com, vb = mybiom_Com$DW_tot, mod = MC_DW_T$mod4A, 
                           ylabl = "Total biomass of plant community (mg)", plot_PAR = T)
                
                          #############################################        
                          ###             traits rac communauté     ###
                          ############################################# 
                          # (à faire)  
                
            #-----------------------------------------------------------
            # 3.0 Effet des traitements à l'échelle spécifique:   
            #            traits particuliers à une espèce végétale
            #------------------------------------------------------------
                          #############################################        
                          ###    Pour les trèfles : nodulation      ###
                          #############################################
                # nodules strate 1
                myTa_nod1 <- myTa[myTa$strata == 1, ]
                myTa_nod1_temp <- myTa_nod1[!is.na(myTa_nod1$NOD),]
                Ta_nod1 <- myfun_Leg(mdata = myTa_nod1_temp, vb = sqrt(myTa_nod1_temp$NOD), seuil_prcRE = 24) # 
                Ta_nod1$resum
                Ta_nod1$select
                pTa_nod1 <- myplot_Leg(mdata = myTa_nod1_temp, vb = myTa_nod1_temp$NOD, mod = Ta_nod1$mod10A, 
                                           ylabl = "T. alexandrium nodules (0-10 cm, ind-1)", plot_PAR = F)
                
                # nodules strate 2
                myTa_nod2 <- myTa[myTa$strata == 2, ]
                myTa_nod2_temp <- myTa_nod2[!is.na(myTa_nod2$NOD),]
                Ta_nod2 <- myfun_Leg(mdata = myTa_nod2_temp, vb = sqrt(myTa_nod2_temp$NOD), seuil_prcRE = 24) # 
                Ta_nod2$resum
                Ta_nod2$select
                pTa_nod2 <- myplot_Leg(mdata = myTa_nod2_temp, vb = myTa_nod2_temp$NOD, mod = Ta_nod2$mod10A, 
                                       ylabl = "T. alexandrium nodules (20-30 cm, ind-1)", plot_PAR = F)
                
                # nodules strate 3
                myTa_nod3 <- myTa[myTa$strata == 3, ]
                myTa_nod3_temp <- myTa_nod3[!is.na(myTa_nod3$NOD),]
                Ta_nod3 <- myfun_Leg(mdata = myTa_nod3_temp, vb = sqrt(myTa_nod3_temp$NOD), seuil_prcRE = 24) # 
                Ta_nod3$resum
                Ta_nod3$select
                pTa_nod3 <- myplot_Leg(mdata = myTa_nod3_temp, vb = myTa_nod3_temp$NOD, mod = Ta_nod3$mod10A, 
                                       ylabl = "T. alexandrium nodules (30-50 cm, ind-1)", plot_PAR = F)
                
                # nodules strate 4
                myTa_nod4 <- myTa[myTa$strata == 4, ]
                myTa_nod4_temp <- myTa_nod4[!is.na(myTa_nod4$NOD),]
                Ta_nod4 <- myfun_Leg(mdata = myTa_nod4_temp, vb = sqrt(myTa_nod4_temp$NOD), seuil_prcRE = 24) # 
                Ta_nod4$resum
                Ta_nod4$select
                pTa_nod4 <- myplot_Leg(mdata = myTa_nod4_temp, vb = myTa_nod4_temp$NOD, mod = Ta_nod4$mod10A, 
                                       ylabl = "T. alexandrium nodules (50-100 cm, ind-1)", plot_PAR = F)
                
                # Nodule systèle racinaire total
                vb <- c("NOD", "DW")
                myTa_nod <- aggregate(myTa[,vb], list(RT = myTa$RT, numPL = myTa$numPL,  div = myTa$div,  leg = myTa$leg, 
                                                id_loca = myTa$id_loca, id_taxo = myTa$id_taxo, distance = myTa$distance, serie = myTa$serie, position = myTa$position, PAR = myTa$PAR), sum)
                myTa_nod_temp <- myTa_nod[!is.na(myTa_nod$NOD),]
                Ta_nod <- myfun_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NOD), seuil_prcRE = 24) ## log.transformed
                Ta_nod$select
                pTa_nod <- myplot_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NOD), mod = Ta_nod$mod10A, 
                                       ylabl = "T. alexandrium nodules (total, ind-1, log transf.)", plot_PAR = F)
                
                
                # Profil racinaire      
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
                Ta_nod_PR
                
                # nodules par unité de masse racinaire
                myTa_nod$NODms <- myTa_nod$NOD/myTa_nod$DW
                myTa_nod_temp <- myTa_nod[!is.na(myTa_nod$NODms),]
                Ta_nodms <- myfun_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NODms), seuil_prcRE = 24) ## log.transformed
                Ta_nodms$select
                pTa_nodms <- myplot_Leg(mdata = myTa_nod_temp, vb = log1p(myTa_nod_temp$NODms), mod = Ta_nod$mod10A, 
                                      ylabl = "T. alexandrium nodules (total, mg-1, log transf.)", plot_PAR = F)
                
                # relation nodules - masse racinaire
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
                
                
                rel_ta + geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div0",],aes(x =DW, y =NOD ), 
                                  method = "lm", alpha = 0.1, color="black", se=F)+
                  annotate(geom="text", x=77, y=45, label="r²=0.84", size=5) +
                  geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div1",],aes(x =DW, y =NOD ), 
                                method = "lm", alpha = 0.1, color="#ffeda0", se=F)+
                  annotate(geom="text", x=75, y=75, label="r²=0.79", size=5, color="#ffeda0") +
                  geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div2",],aes(x =DW, y =NOD ), 
                              method = "lm", alpha = 0.1, color="#feb24c", se=F)+
                  annotate(geom="text", x=90, y=90, label="r²=0.89", size=5,color="#feb24c") +
                  geom_smooth(data = Ta_RT_nod[Ta_RT_nod$div=="Div3",],aes(x =DW, y =NOD ), 
                              method = "lm", alpha = 0.1, color="#f03b20", se=F)+
                  annotate(geom="text", x=75, y=60, label="r²=0.89", size=5,color="#f03b20")
                # + geom_abline(intercept=0, slope=0.76425, size=0.8, color="black")+
                
                # on ajoute les regressions
                reg_ta_div0 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div0",])
                summary(reg_ta_div0)
                reg_ta_div1 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div1",])
                summary(reg_ta_div1)
                reg_ta_div2 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div2",])
                summary(reg_ta_div2)
                reg_ta_div3 <- lm(NOD~DW, data=Ta_RT_nod[Ta_RT_nod$div=="Div3",])
                summary(reg_ta_div3)
                
                              #############################################        
                              ###     Pour les ray-grass : clonalité    ###
                              #############################################
                              # (à faire)
            #-----------------------------------------------------------
            # 3.1 Effet des traitements à l'échelle spécifique:   
            #            biomasses végétales sèches
            #------------------------------------------------------------
                          #############################################        
                          ###           Pour les trèfles            ###
                          #############################################
                                # Biomasse sèche aérienne
                                    Ta_A <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 0, ]
                                    Ta_A_DW <- myfun_Leg(mdata = Ta_A, vb = log(Ta_A$DW), seuil_prcRE = 30) # transf log
                                    Ta_A_DW$resum
                                    Ta_A_DW$select
                                    pTa_A_DW <- myplot_Leg(mdata = Ta_A, vb = Ta_A$DW, mod = Ta_A_DW$mod10A, 
                                               ylabl = "T. alexandrium \nshoot mass (DW, mg)", plot_PAR = T)
                        
                                # Biomasse sèche strate 1            
                                    Ta_R1 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 1, ]
                                    Ta_R1_DW <- myfun_Leg(mdata = Ta_R1, vb = sqrt(Ta_R1$DW), seuil_prcRE = 30) # idem
                                    Ta_R1_DW$resum
                                    Ta_R1_DW$select
                                    pTa_R1_DW <- myplot_Leg(mdata = Ta_R1, vb = sqrt(Ta_R1$DW), mod = Ta_R1_DW$mod10A, 
                                               ylabl = "T. alexandrium \nroot mass (0-10 cm, DW mg, sqrt transf.)", plot_PAR = T)
                      
                                # Biomasse sèche strate 2            
                                    Ta_R2 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 2, ]
                                    Ta_R2_DW <- myfun_Leg(mdata = Ta_R2, vb = sqrt(Ta_R2$DW), seuil_prcRE = 30) # exemple avec transf. sqrt()
                                    Ta_R2_DW$resum
                                    pTa_R2_DW <- myplot_Leg(mdata = Ta_R2, vb = sqrt(Ta_R2$DW), mod = Ta_R2_DW$mod10A, 
                                               ylabl = "T. alexandrium \nroot mass (10-30 cm, DW mg, sqrt transf.)", plot_PAR = T)
                                    
                                # Biomasse sèche strate 3            
                                    Ta_R3 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 3, ]
                                    Ta_R3_DW <- myfun_Leg(mdata = Ta_R3, vb = sqrt(Ta_R3$DW), seuil_prcRE = 30) # exemple avec transf. sqrt()
                                    Ta_R3_DW$resum
                                    pTa_R3_DW <- myplot_Leg(mdata = Ta_R3, vb = sqrt(Ta_R3$DW), mod = Ta_R3_DW$mod10A, 
                                               ylabl = "T. alexandrium \nroot mass (30-50 cm, DW mg, sqrt transf.)", plot_PAR = T)
                                    
                                # Biomasse sèche strate 4            
                                    Ta_R4 <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata == 4, ]
                                    Ta_R4_DW <- myfun_Leg(mdata = Ta_R4, vb = Ta_R4$DW, seuil_prcRE = 30) 
                                    # n'arrive pas à créer les modèles car design hyper déséquilibré (cf summary(Ta_R4$div))
                                    pTa_R4_DW <- myplot_Leg(mdata = Ta_R4, vb = sqrt(Ta_R4$DW), mod = Ta_R4_DW$mod10A, 
                                                            ylabl = "T. alexandrium \nroot mass (50-90 cm, DW mg, sqrt transf.)", plot_PAR = T)
                                    
                                        # Biomasse sèche racine totale
                                            Ta_R <- mybiom_Ind[mybiom_Ind$id_taxo == "Tri_alexa" & mybiom_Ind$strata != 0, ]
                                            Ta_R <- aggregate(Ta_R$DW, list(RT = Ta_R$RT, numPL = Ta_R$numPL,  div = Ta_R$div,  leg = Ta_R$leg, 
                                                                            id_loca = Ta_R$id_loca, id_taxo = Ta_R$id_taxo, distance = Ta_R$distance, serie = Ta_R$serie, position = Ta_R$position, PAR = Ta_R$PAR), sum)
                                            colnames(Ta_R)[11] <- "DW"
                                            Ta_R_DW <- myfun_Leg(mdata = Ta_R, vb = log1p(Ta_R$DW), seuil_prcRE = 30) ## log.transformed
                                            pTa_R_DW <- myplot_Leg(mdata = Ta_R, vb = log1p(Ta_R$DW), mod = Ta_R_DW$mod10A, 
                                                                    ylabl = "T. alexandrium \nroot mass (total, DW mg, log transf.)", plot_PAR = T)
                                            
                                            
                                       # Profil racinaire      
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
                                              
                                            Ta_PR
                                        
                                        # Shoot/root ratio
                                            Ta_R <- Ta_R[order(Ta_R$id_loca),]
                                            Ta_A <- Ta_A[order(Ta_A$id_loca),]
                                            Ta_R$id_loca == Ta_A$id_loca ## pour vérifier que les plantes sont dans le même ordre dans chacune des df
                                            Ta_A$SR <- Ta_A$DW/Ta_R$DW
                                            Ta_SR_DW <- myfun_Leg(mdata = Ta_A, vb = sqrt(Ta_A$SR), seuil_prcRE = 30)  # aucune transf. ne permet de normaliser les résidus
                                            Ta_SR_DW$resum
                                            pTa_SR_DW <- myplot_Leg(mdata = Ta_A, vb = sqrt(Ta_A$SR), mod = Ta_SR_DW$mod9A, 
                                                                   ylabl = "T. alexandrium\nShoot:Root ratio (DW)", plot_PAR = T)
                                            
                                            pTa_SR_DW$pDIV
                                            
                                                # Exemple de tableau récapitulatif
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
                                                    
                              #############################################        
                              ###           Pour les colzas             ###
                              #############################################
                                # Date de germination
                                  myBnGerm_m <- myfun_Ind(mdata = myBnGerm, vb = myBnGerm$SGDm, seuil_prcRE = 30)
                                                    
                                # Biomasse sèche aérienne
                                  Bn_A_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 0, ]
                                          Bn_A_all_DW <- myfun_Ind(mdata = Bn_A_all, vb = sqrt(Bn_A_all$DW), seuil_prcRE = 30) # transf. sqrt
                                          Bn_A_all_DW$select    # ni la diversité des VdT (ou leur présence) ni la présence de légumineuse ne change pas la biomasse aérienne du colza
                                  pBn_A_all_DW <- myplot_Ind(mdata = Bn_A_all, vb = sqrt(Bn_A_all$DW), mod = Bn_A_all_DW$mod10A, 
                                               ylabl = "B. napus\nshoot mass (DW, mg)", plot_PAR = F)
                                  pBn_A_all_DW$pDIV_INTER
                                  
                                # Biomasse sèche strate 1   
                                  Bn_R1_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 1, ]
                                  Bn_R1_all_DW <- myfun_Ind(mdata = Bn_R1_all, vb = sqrt(Bn_R1_all$DW), seuil_prcRE = 30)
                                  Bn_R1_all_DW$select
                                  pBn_R1_all_DW <- myplot_Ind(mdata = Bn_R1_all, vb = sqrt(Bn_R1_all$DW), mod = Bn_R1_all_DW$mod10A, 
                                                           ylabl = "B. napus\nroot mass (0-10 cm, DW mg)", plot_PAR = F)
                         
                                # Biomasse sèche strate 2  
                                    Bn_R2_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 2, ]
                                    Bn_R2_all_DW$select
                                    pBn_R2_all_DW <- myplot_Ind(mdata = Bn_R2_all, vb = sqrt(Bn_R2_all$DW), mod = Bn_R2_all_DW$mod10A, 
                                               ylabl = "B. napus\nroot mass (10-30 cm, DW mg)", plot_PAR = F)
                                    
                                
                                 # Biomasse sèche strate 3  
                                    Bn_R3_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 3, ]
                                    Bn_R3_all_DW$select
                                    pBn_R3_all_DW <- myplot_Ind(mdata = Bn_R3_all, vb = sqrt(Bn_R3_all$DW), mod = Bn_R3_all_DW$mod10A, 
                                               ylabl = "B. napus\nroot mass (30-50 cm, DW mg)", plot_PAR = T)
                                
                                # Biomasse sèche strate 4      
                                    Bn_R4_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata == 4, ]
                                    Bn_R4_all_DW$select
                                    pBn_R4_all_DW <- myplot_Ind(mdata = Bn_R4_all, vb = sqrt(Bn_R4_all$DW), mod = Bn_R4_all_DW$mod10A, 
                                               ylabl = "B. napus\nroot mass (50-90 cm, DW mg)", plot_PAR = T)
                                    
                                # Biomasse sèche racine totale
                                      Bn_R <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata != 0, ]
                                      Bn_R <- aggregate(Bn_R$DW, list(RT = Bn_R$RT, numPL = Bn_R$numPL,  div = Bn_R$div,  leg = Bn_R$leg, 
                                                                      id_loca = Bn_R$id_loca, id_taxo = Bn_R$id_taxo, distance = Bn_R$distance, 
                                                                      serie = Bn_R$serie, position = Bn_R$position, PAR = Bn_R$PAR), sum)
                                      colnames(Bn_R)[11] <- "DW"
                                      Bn_R_DW <- myfun_Ind(mdata = Bn_R, vb = sqrt(Bn_R$DW), seuil_prcRE = 30) ## sqrt transformed
                                      pBn_R_DW <- myplot_Ind(mdata = Bn_R, vb = sqrt(Bn_R$DW), mod = Bn_R_DW$mod10A, 
                                                             ylabl = "B. napus \nroot mass (total, mg, sqrt transf.)", plot_PAR = T)
                                      pBn_R_DW

                                # Biomasse sèche racines fines totale
                                      Bn_R <- mybiom_Ind[mybiom_Ind$id_taxo == "Bra_napus" & mybiom_Ind$strata != 0 & mybiom_Ind$, ]
                                      Bn_R <- aggregate(Bn_R$DW, list(RT = Bn_R$RT, numPL = Bn_R$numPL,  div = Bn_R$div,  leg = Bn_R$leg, 
                                                                      id_loca = Bn_R$id_loca, id_taxo = Bn_R$id_taxo, distance = Bn_R$distance, 
                                                                      serie = Bn_R$serie, position = Bn_R$position, PAR = Bn_R$PAR), sum)
                                      colnames(Bn_R)[11] <- "DW"
                                      Bn_R_DW <- myfun_Ind(mdata = Bn_R, vb = sqrt(Bn_R$DW), seuil_prcRE = 30) ## sqrt transformed
                                      pBn_R_DW <- myplot_Ind(mdata = Bn_R, vb = sqrt(Bn_R$DW), mod = Bn_R_DW$mod10A, 
                                                             ylabl = "B. napus \nroot mass (total, mg, sqrt transf.)", plot_PAR = T)
                                      pBn_R_DW                                      
                                # Profil racinaire      
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
                                          #scale_shape_discrete(name="Earthworm FRic")+
                                          xlim(c(-90, 0))+
                                          ylab ("B. napus \nroot mass profile (DW mg)")+
                                          xlab("Soil depth\n(cm)")+
                                          scale_linetype_discrete(guide=FALSE)+
                                          scale_shape_discrete(guide = FALSE)+
                                          theme(legend.justification=c(0,0), legend.position=c(1,1))+
                                          theme(panel.background = element_rect(fill="white", colour = "black"))+
                                          geom_vline(xintercept = c(0, -10, -30, -50, -90), colour = "gray")
                                        Bn_PR
                                        
                                    # Shoot/root ratio
                                        Bn_A_all$R_DW <- NULL
                                        for (i in levels(Bn_A_all$id_loca))
                                          Bn_A_all$R_DW[Bn_A_all$id_loca == i] <- Bn_R$DW[Bn_R$id_loca == i]
                                        Bn_A_all$SR <- Bn_A_all$DW/Bn_A_all$R_DW
                                        Bn_SR_DW <- myfun_Ind(mdata = Bn_A_all, vb = asin(tan(Bn_A_all$SR)), seuil_prcRE = 30)  # asin(tan()), transf sqrt et log n'améliorent pas la distribution des résidus
                                        Bn_SR_DW$select
                                        pBn_SR_DW <- myplot_Ind(mdata = Bn_A_all, vb = Bn_A_all$SR, mod = Bn_SR_DW$mod10A, 
                                                               ylabl = "B. napus \nShoot:Root ratio (DW)", plot_PAR = T)
                                        pBn_SR_DW
                                        
                                 # Tableau récapitulatif
                                    
                                    
                                    #############################################        
                                    ###          Pour les ray-grass           ###
                                    #############################################
                                  # Biomasse sèche aérienne
                                      Lp_A_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 0, ]
                                      Lp_A_all_DW <- myfun_Ind(mdata = Lp_A_all, vb = log1p(Lp_A_all$DW), seuil_prcRE = 24) 
                                      Lp_A_all_DW$resum
                                      Lp_A_all_DW$select
                                      pLp_A_all_DW <- myplot_Ind(mdata = Lp_A_all, vb = sqrt(Lp_A_all$DW), mod = Lp_A_all_DW$mod9A, 
                                                 ylabl = "L. perenne\n shoot mass (DW mg)", plot_PAR = T)
                       
                                  # Biomasse sèche strate 1
                                      Lp_R1_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 1, ]
                                      # test si les différentes options de "nettoyage" du dataset apporte quelquechose
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
                                # Biomasse sèche strate 2
                                      Lp_R2_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 2, ]
                                      # test si les différentes options de "nettoyage" du dataset apporte quelquechose
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
                                      
                                # Biomasse sèche strate 3
                                      Lp_R3_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 3, ]
                                      # test si les différentes options de "nettoyage" du dataset apporte quelquechose
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
                                # Biomasse sèche strate 4
                                     Lp_R4_all <- mybiom_Ind[mybiom_Ind$id_taxo == "Lol_peren" & mybiom_Ind$strata == 4, ]
                                     # test si les différentes options de "nettoyage" du dataset apporte quelquechose
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
                                
                                 # Biomasse sèche racine totale
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
                                     Lp_SR_DW <- myfun_Ind(mdata = Lp_SR, vb = Lp_SR$SR, seuil_prcRE = 30)  # asin fonctionne, mais pb ensuite dans myplot_Ind; transf sqrt et log n'améliorent pas la distribution des résidus
                                     Lp_SR_DW$select
                                     pLp_SR_DW <- myplot_Ind(mdata = Lp_SR, vb = Lp_SR$SR, mod = Lp_SR_DW$mod10A, 
                                                             ylabl = "L. perenne\n Shoot:Root ratio (DW)", plot_PAR = T)
                                     pLp_SR_DW
                                     
                                           #############################################        
                                           ###          Pour les géraniums           ###
                                           #############################################
                                # Biomasse sèche aérienne
                                    Gd_A <- mybiom_Ind[mybiom_Ind$id_taxo == "Ger_disse" & mybiom_Ind$strata == 0, ]
                                    Gd_A_DW <- myfun_Ind(mdata = Gd_A, vb = asin(tan(Gd_A$DW)), seuil_prcRE = 30) # transf log et sqrt ne nrmalisent pas les résidus, asin(tan()) pour l'exemple
                                    Gd_A_DW$select
                                    pGd_A_DW <- myplot_Ind(mdata = Gd_A, vb = Gd_A$DW, mod = Gd_A_DW$mod10B, 
                                                            ylabl = "L. perenne\n Shoot:Root ratio (DW)", plot_PAR = T)
                                    pGd_A_DW$pDIV
                                    Gd_R <- mybiom_Ind[mybiom_Ind$id_taxo == "Ger_disse" & mybiom_Ind$strata > 0, ]
                          
                                                #############################################        
                                                ###          Pour les véroniques          ###
                                                #############################################
                                # Biomasse sèche aérienne
                                    Vp_A <- mybiom_Ind[mybiom_Ind$id_taxo == "Ver_persi" & mybiom_Ind$strata == 0, ]
                                    Vp_A_DW <- myfun_Ind(mdata = Vp_A, vb = log1p(Vp_A$DW), seuil_prcRE = 30)
                                    Vp_A_DW$select
                                    myplot_Ind(mdata = Vp_A, vb = log1p(Vp_A$DW), mod = Vp_A_DW$mod10A, 
                                               ylabl = "V. persica aerial biomass (log(g))", plot_PAR = T)
                      
                                    
            #-----------------------------------------------------------
            # 3.2 Effet des traitements à l'échelle spécifique:   
            #         [N] dans les organes de colza, ray-grass et trèfle
            #------------------------------------------------------------
                        #############################################        
                        ###           Pour les trèfles            ###
                        ############################################# 
                                  Ta_N <- myNcont[myNcont$id_taxo == "Tri_alexa" , ]
                              # dans les racines fines
                                  Ta_N_temp <- Ta_N[!is.na(Ta_N$'Fine roots'),]
                                  Ta_N_fRoots <- myfun_Leg(mdata = Ta_N_temp, vb = log(Ta_N_temp$'Fine roots'), seuil_prcRE = 30) # pas de transf efficaces
                                  Ta_N_fRoots$resum
                                  Ta_N_fRoots$select
                                  pTa_N_fRoots <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$'Fine roots', mod = Ta_N_fRoots$mod10A, 
                                                         ylabl = "T. alexandrium \n roots N content (mg g-1)", plot_PAR = T)
                                  pTa_N_fRoots$pDIV
                              # dans les limbes matures
                                  Ta_N_temp <- Ta_N[!is.na(Ta_N$Limb),]
                                  Ta_N_limb <- myfun_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Limb, seuil_prcRE = 24) # pas de transf 
                                  Ta_N_limb$resum
                                  Ta_N_limb$select
                                  pTa_N_limb <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Limb, mod = Ta_N_limb$mod10A, 
                                                             ylabl = "T. alexandrium LNC (mg g-1)", plot_PAR = F)
                                  pTa_N_limb$pDIV
                              # dans les limbes d'autres feuilles
                                  Ta_N_temp <- Ta_N[!is.na(Ta_N$Leaves),]
                                  Ta_N_leaves <- myfun_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Leaves, seuil_prcRE = 30) # pas de transf 
                                  Ta_N_leaves$resum
                                  pTa_N_leaves <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Leaves, mod = Ta_N_leaves$mod9A, 
                                                           ylabl = "T. alexandrium \nleaves N content (mg g-1)", plot_PAR = T)
                                  pTa_N_leaves$pDIV
                              # dans les tiges
                                  Ta_N_temp <- Ta_N[!is.na(Ta_N$Stem),]
                                  Ta_N_stem <- myfun_Leg(mdata = Ta_N_temp, vb = asin(tan(Ta_N_temp$Stem)), seuil_prcRE = 30) 
                                  Ta_N_stem$resum
                                  pTa_N_stem <- myplot_Leg(mdata = Ta_N_temp, vb = Ta_N_temp$Stem, mod = Ta_N_stem$mod9A, 
                                                             ylabl = "T. alexandrium \nstem N content (mg g-1)", plot_PAR = T)             
                                  pTa_N_stem$pDIV
                              # comparaison de la distribution des concentrations
                                  Ta_meaN <- rbind.data.frame(pTa_N_limb$moyFRic, pTa_N_leaves$moyFRic, pTa_N_stem$moyFRic, pTa_N_fRoots$moyFRic)
                                  Ta_meaN$organes <- rep(c("a.limb", "b.leaves", "c.stem", "d.roots"), each = 4)
                                  B <- ggplot(data = Ta_meaN, aes(x = organes, y = x))+
                                    geom_point(aes(shape = div), size = 4)+
                                    geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
                                    scale_shape_discrete(guide = FALSE)+
                                    ylab("N content (mg g-1)")   +
                                    theme(panel.background = element_rect(fill="white", colour = "black"))+
                                    facet_grid(div~.)
                                  B
                            #############################################        
                            ###           Pour les colzas             ###
                            #############################################
                                  Bn_N <- myNcont[myNcont$id_taxo == "Bra_napus" , ]
                                  # dans les racines fines
                                  Bn_N_temp <- Bn_N[!is.na(Bn_N$'Fine roots'),]
                                  Bn_N_fRoots <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$'Fine roots', seuil_prcRE = 24) 
                                  Bn_N_fRoots$resum # pas de transf efficaces
                                  Bn_N_fRoots$select
                                  pBn_N_fRoots <- myplot_Ind(mdata = Bn_N_fRoots, vb = log(Bn_N_fRoots$'Fine roots'), mod = Bn_N_fRoots$mod10A, 
                                                            ylabl = "B. napus \nfine roots N content (mg g-1)", plot_PAR = T)
                                  # vérifier NA's

                                  # dans les racines pivot
                                  Bn_N_temp <- Bn_N[!is.na(Bn_N$Taproot),]
                                  Bn_N_tRoots <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Taproot, seuil_prcRE = 30)
                                  Bn_N_tRoots$resum
                                  Bn_N_tRoots$select
                                  pBn_N_tRoots <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Taproot, mod = Bn_N_tRoots$mod10A, 
                                                             ylabl = "B. napus \ntaproots N content (mg g-1)", plot_PAR = T)
                                  pBn_N_tRoots$pDIV_INTER
                                  
                                  # dans les limbes matures
                                  Bn_N_temp <- Bn_N[!is.na(Bn_N$Limb),]
                                  Bn_N_limb <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Limb, seuil_prcRE = 24) # pas de transf 
                                  Bn_N_limb$resum
                                  Bn_N_limb$select
                                  pBn_N_limb <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Limb, mod = Bn_N_limb$mod10A, 
                                                           ylabl = "B. napus mature LNC (mg g-1)", plot_PAR = F)
                                  pBn_N_limb$pDIV_INTER
                                  
                                  # dans les limbes d'autres feuilles
                                  Bn_N_temp <- Bn_N[!is.na(Bn_N$Leaves),]
                                  Bn_N_leaves <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Leaves, seuil_prcRE = 30) # pas de transf 
                                  Bn_N_leaves$resum
                                  Bn_N_leaves$select
                                  pBn_N_leaves <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Leaves, mod = Bn_N_leaves$mod10A, 
                                                             ylabl = "B. napus \nleaves N content (mg g-1)", plot_PAR = T)
                                  
                                  # dans les tiges
                                  Bn_N_temp <- Bn_N[!is.na(Bn_N$Stem),]
                                  Bn_N_stem <- myfun_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Stem, seuil_prcRE = 30) 
                                  Bn_N_stem$resum
                                  Bn_N_stem$select
                                  pBn_N_stem <- myplot_Ind(mdata = Bn_N_temp, vb = Bn_N_temp$Stem, mod = Bn_N_stem$mod10A, 
                                                           ylabl = "B. napus \nstem N content (mg g-1)", plot_PAR = T)             
                                  
                                  # comparaison de la distribution des concentrations
                                  Bn_meaN <- rbind.data.frame(pBn_N_limb$moyFRic, pBn_N_leaves$moyFRic, pBn_N_stem$moyFRic, pBn_N_tRoots$moyFRic, pBn_N_fRoots$moyFRic)
                                  Bn_meaN$organes <- as.factor(rep(c("a.limb", "b.leaves", "c.stem", "d.taproots", "e.fine roots"), each = 4))
                                  D <- ggplot(data = Bn_meaN, aes(x = organes, y = x))+
                                    geom_point(aes(shape = div), size = 4)+
                                    geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
                                    scale_shape_discrete(guide = FALSE)+
                                    ylab("N content (mg g-1)")   +
                                    theme(panel.background = element_rect(fill="white", colour = "black"))+
                                    facet_grid(div~.)
                                  D
                                  
                              
                                  #############################################        
                                  ###           Pour les ray-grass          ###
                                  #############################################
                                  Lp_N <- myNcont[myNcont$id_taxo == "Lol_peren" , ]
                                # dans les racines fines
                                  Lp_N_temp <- Lp_N[!is.na(Lp_N$'Fine roots'),]
                                  Lp_N_fRoots <- myfun_Ind(mdata = Lp_N_temp, vb = sqrt(Lp_N_temp$'Fine roots'), seuil_prcRE = 30) 
                                  Lp_N_fRoots$resum # pas de transf efficaces
                                  Lp_N_fRoots$select
                                    # test avec les options de filtrage
                                    Lp_N_fRoots_opt2 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
                                    Lp_N_fRoots_opt4 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
                                    Lp_N_fRoots <- myfun_Ind(mdata = Lp_N_fRoots_opt2, vb = sqrt(Lp_N_fRoots_opt2$'Fine roots'), seuil_prcRE = 30)
                                    # aucune n'est satisfaisante => Kruskal-Wallis ?
                                  pLp_N_fRoots <- myplot_Ind(mdata = Lp_N_temp, vb = Lp_N_temp$'Fine roots', mod = Lp_N_fRoots$mod10A, 
                                             ylabl = "L. perenne RNC (mg g-1)", plot_PAR = F)
                                  pLp_N_fRoots$pDIV
                                
                                # dans les limbes matures
                                  Lp_N_temp <- Lp_N[!is.na(Lp_N$Limb),]
                                  Lp_N_limb <- myfun_Ind(mdata = Lp_N_temp, vb = log(Lp_N_temp$Limb), seuil_prcRE = 24) # pas de transf efficaces
                                  Lp_N_limb$resum
                                  Lp_N_limb$select
                                    # test avec les options de filtrage
                                    #Lp_N_limb_opt2 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
                                    #Lp_N_limb_opt4 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
                                    #Lp_N_limb <- myfun_Ind(mdata = Lp_N_limb_opt4, vb = sqrt(Lp_N_limb_opt4$'Fine roots'), seuil_prcRE = 30)
                                    # aucune n'est satisfaisante => Kruskal-Wallis ?
                                  pLp_N_limb <- myplot_Ind(mdata = Lp_N_temp, vb = Lp_N_temp$Limb, mod = Lp_N_limb$mod10A, 
                                                           ylabl = "L. perenne LNC (mg g-1)", plot_PAR =F)
                                  pLp_N_limb$pDIV
                                
                                # dans les limbes d'autres feuilles
                                  Lp_N_temp <- Lp_N[!is.na(Lp_N$Leaves),]
                                  Lp_N_leaves <- myfun_Ind(mdata = Lp_N_temp, vb = Lp_N_temp$Leaves, seuil_prcRE = 24) # pas de transf 
                                  Lp_N_leaves$resum
                                  Lp_N_leaves$select
                                    # test avec les options de filtrage
                                    #Lp_N_leaves_opt2 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
                                    #Lp_N_leaves_opt4 <- Lp_N_temp[which(Lp_N_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
                                    #Lp_N_leaves <- myfun_Ind(mdata = Lp_N_leaves_opt2, vb = sqrt(Lp_N_leaves_opt2$'Fine roots'), seuil_prcRE = 30)
                                    # aucune n'est satisfaisante => Kruskal-Wallis ?
                                pLp_N_leaves <- myplot_Ind(mdata = Lp_N_temp[!is.na(Lp_N_temp$Leaves),], vb = Lp_N_temp$Leaves[!is.na(Lp_N_temp$Leaves)], mod = Lp_N_leaves$mod10A, 
                                                             ylabl = "L. perenne leaves N content (mg g-1)", plot_PAR =F)
                                    
                        

                                
                                
                                #-----------------------------------------------------------
                                # 3.3 Effet des traitements à l'échelle spécifique:   
                                #         [P] dans les organes de colza, ray-grass et trèfle
                                #------------------------------------------------------------
                                #############################################        
                                ###           Pour les trèfles            ###
                                ############################################# 
                                Ta_P <- myPcont[myPcont$id_taxo == "Tri_alexa" , ]
                                # dans les racines fines
                                Ta_P_temp <- Ta_P[which(!is.na(Ta_P$'Fine roots') & Ta_P$`Fine roots` > 0.1),]
                                Ta_P_fRoots <- myfun_Leg(mdata = Ta_P_temp, vb = log(Ta_P_temp$'Fine roots'), seuil_prcRE = 30) # pas de transf efficaces
                                Ta_P_fRoots$resum
                                Ta_P_fRoots$select
                                pTa_P_fRoots <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$'Fine roots', mod = Ta_P_fRoots$mod10A, 
                                                           ylabl = "T. alexandrium \n roots P content (mg g-1)", plot_PAR = T)
                                pTa_P_fRoots$pDIV
                                
                                # dans les limbes matures
                                Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Limb) & Ta_P$Limb > 0.01),]
                                Ta_P_limb <- myfun_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Limb, seuil_prcRE = 24) # pas de transf 
                                Ta_P_limb$resum
                                Ta_P_limb$select
                                pTa_P_limb <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Limb, mod = Ta_P_limb$mod10A, 
                                                         ylabl = "T. alexandrium LPC (mg g-1)", plot_PAR = F)
                                
                                # dans les limbes d'autres feuilles
                                Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Leaves) & Ta_P$Leaves > 0.1),]
                                Ta_P_leaves <- myfun_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Leaves, seuil_prcRE = 30) # pas de transf 
                                Ta_P_leaves$resum
                                pTa_P_leaves <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Leaves, mod = Ta_P_leaves$mod9A, 
                                                           ylabl = "T. alexandrium \nleaves P content (mg g-1)", plot_PAR = T)
                                pTa_P_leaves$pDIV
                                
                                # dans les tiges
                                Ta_P_temp <- Ta_P[which(!is.na(Ta_P$Stem) & Ta_P$Stem > 0.1),]
                                Ta_P_stem <- myfun_Leg(mdata = Ta_P_temp, vb = asin(tan(Ta_P_temp$Stem)), seuil_prcRE = 30) 
                                Ta_P_stem$resum
                                pTa_P_stem <- myplot_Leg(mdata = Ta_P_temp, vb = Ta_P_temp$Stem, mod = Ta_P_stem$mod9A, 
                                                         ylabl = "T. alexandrium \nstem N content (mg g-1)", plot_PAR = T)             
                                
                                # comparaison de la distribution des concentrations
                                Ta_meaN <- rbind.data.frame(pTa_P_limb$moyFRic, pTa_P_leaves$moyFRic, pTa_P_stem$moyFRic, pTa_P_fRoots$moyFRic)
                                Ta_meaN$organes <- rep(c("a.limb", "b.leaves", "c.stem", "d.roots"), each = 4)
                                B <- ggplot(data = Ta_meaN, aes(x = div, y = x))+
                                  geom_point(aes(shape = div), size = 4)+
                                  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
                                  scale_shape_discrete(guide = FALSE)+
                                  ylab("P content (mg g-1)")   +
                                  theme(panel.background = element_rect(fill="white", colour = "black"))+
                                  facet_grid(organes~.)
                                B
                                #############################################        
                                ###           Pour les colzas             ###
                                #############################################
                                Bn_P <- myPcont[myPcont$id_taxo == "Bra_napus" , ]
                                # dans les racines fines
                                Bn_P_temp <- Bn_P[!is.na(Bn_P$'Fine roots'),]
                                Bn_P_fRoots <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$'Fine roots', seuil_prcRE = 24) 
                                Bn_P_fRoots$resum # pas de transf efficaces
                                Bn_P_fRoots$select
                                # test avec les options de filtrage
                                #Bn_P_fRoots_opt8 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt8==1]), ]
                                #Bn_P_fRoots_opt7 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt7==1]), ]
                                #Bn_P_fRoots_opt6 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt6==1]), ]
                                #Bn_P_fRoots_opt4 <- Bn_P_temp[which(Bn_P_temp$id_loca %in% filter_bn$id_loca[filter_bn$opt4==1]), ]
                                #Bn_P_fRoots <- myfun_Ind(mdata = Bn_P_fRoots_opt4, vb = log(Bn_P_fRoots_opt4$'Fine roots'), seuil_prcRE = 30)  # il faut l'option 4 pour réussir à obtenir des résidus normaux
                                #Bn_P_fRoots$resum
                                #Bn_P_fRoots$select
                                pBn_P_fRoots <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$'Fine roots', mod = Bn_P_fRoots$mod10A, 
                                                           ylabl = "B. napus \nfine roots P content (mg g-1)", plot_PAR = T)
                                pBn_P_fRoots$pDIV_INTER
                                
                                # dans les racines pivot
                                Bn_P_temp <- Bn_P[!is.na(Bn_P$Taproot),]
                                Bn_P_tRoots <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Taproot, seuil_prcRE = 30)
                                Bn_P_tRoots$resum
                                Bn_P_tRoots$select
                                pBn_P_tRoots <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Taproot, mod = Bn_P_tRoots$mod10A, 
                                                           ylabl = "B. napus \ntaproots P content (mg g-1)", plot_PAR = T)
                                pBn_P_tRoots$pDIV
                                
                                # dans les limbes matures
                                Bn_P_temp <- Bn_P[!is.na(Bn_P$Limb),]
                                Bn_P_limb <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Limb, seuil_prcRE = 24) # pas de transf 
                                Bn_P_limb$resum
                                Bn_P_limb$select
                                pBn_P_limb <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Limb, mod = Bn_P_limb$mod10A, 
                                                         ylabl = "B. napus mature LPC (mg g-1)", plot_PAR = F)
                                pBn_P_limb$pDIV_INTER
                                
                                # dans les limbes d'autres feuilles
                                Bn_P_temp <- Bn_P[!is.na(Bn_P$Leaves),]
                                Bn_P_leaves <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Leaves, seuil_prcRE = 30) # pas de transf 
                                Bn_P_leaves$resum
                                Bn_P_leaves$select
                                pBn_P_leaves <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Leaves, mod = Bn_P_leaves$mod10A, 
                                                           ylabl = "B. napus \nleaves P content (mg g-1)", plot_PAR = T)
                                
                                # dans les tiges
                                Bn_P_temp <- Bn_P[!is.na(Bn_P$Stem),]
                                Bn_P_stem <- myfun_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Stem, seuil_prcRE = 30) 
                                Bn_P_stem$resum
                                Bn_P_stem$select
                                pBn_P_stem <- myplot_Ind(mdata = Bn_P_temp, vb = Bn_P_temp$Stem, mod = Bn_P_stem$mod10A, 
                                                         ylabl = "B. napus \nstem P content (mg g-1)", plot_PAR = T)             
                                
                                # comparaison de la distribution des concentrations
                                Bn_meaN <- rbind.data.frame(pBn_P_limb$moyFRic, pBn_P_leaves$moyFRic, pBn_P_stem$moyFRic, pBn_P_tRoots$moyFRic, pBn_P_fRoots$moyFRic)
                                Bn_meaN$organes <- as.factor(rep(c("a.limb", "b.leaves", "c.stem", "d.taproots", "e.fine roots"), each = 4))
                                D <- ggplot(data = Bn_meaN, aes(x = div, y = x))+
                                  geom_point(aes(shape = div), size = 4)+
                                  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
                                  scale_shape_discrete(guide = FALSE)+
                                  ylab("P content (mg g-1)")   +
                                  theme(panel.background = element_rect(fill="white", colour = "black"))+
                                  facet_grid(.~organes)
                                D
                                
                                
                                #############################################        
                                ###           Pour les ray-grass          ###
                                #############################################
                                Lp_P <- myPcont[myPcont$id_taxo == "Lol_peren" , ]
                                # dans les racines fines
                                Lp_P_temp <- Lp_P[!is.na(Lp_P$'Fine roots'),]
                                Lp_P_fRoots <- myfun_Ind(mdata = Lp_P_temp, vb = sqrt(Lp_P_temp$'Fine roots'), seuil_prcRE = 30) 
                                Lp_P_fRoots$resum # pas de transf efficaces
                                Lp_P_fRoots$select
                                pLp_P_fRoots <- myplot_Ind(mdata = Lp_P_temp, vb = Lp_P_temp$'Fine roots', mod = Lp_P_fRoots$mod10A, 
                                                           ylabl = "L. perenne RPC (mg g-1)", plot_PAR = F)
                                pLp_P_fRoots$pDIV
                                
                                # dans les limbes matures
                                Lp_P_temp <- Lp_P[!is.na(Lp_P$Limb),]
                                Lp_P_limb <- myfun_Ind(mdata = Lp_P_temp, vb = log(Lp_P_temp$Limb), seuil_prcRE = 24) # pas de transf efficaces
                                Lp_P_limb$resum
                                Lp_P_limb$select
                                pLp_P_limb <- myplot_Ind(mdata = Lp_P_temp, vb = Lp_P_temp$Limb, mod = Lp_P_limb$mod10A, 
                                                         ylabl = "L. perenne LPC (mg g-1)", plot_PAR =F)
                                
                                # dans les limbes d'autres feuilles
                                Lp_P_temp <- Lp_P[!is.na(Lp_P$Leaves),]
                                Lp_P_leaves <- myfun_Ind(mdata = Lp_P_temp, vb = Lp_P_temp$Leaves, seuil_prcRE = 24) # pas de transf 
                                Lp_P_leaves$resum
                                Lp_P_leaves$select
                                # test avec les options de filtrage
                                #Lp_P_leaves_opt2 <- Lp_P_temp[which(Lp_P_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt2==1]), ]
                                #Lp_P_leaves_opt4 <- Lp_P_temp[which(Lp_P_temp$id_loca %in% filter_lp$id_loca[filter_lp$opt4==1]), ]
                                #Lp_P_leaves <- myfun_Ind(mdata = Lp_P_leaves_opt2, vb = sqrt(Lp_P_leaves_opt2$'Fine roots'), seuil_prcRE = 30)
                                # aucune n'est satisfaisante => Kruskal-Wallis ?
                                pLp_P_leaves <- myplot_Ind(mdata = Lp_P_temp[!is.na(Lp_P_temp$Leaves),], vb = Lp_P_temp$Leaves[!is.na(Lp_P_temp$Leaves)], mod = Lp_P_leaves$mod10A, 
                                                           ylabl = "L. perenne leaves P content (mg g-1)", plot_PAR =F)
                                
                                pLp_P_leaves$pDIV_INTER
                                
                                # comparaison de la distribution des concentrations
                                Lp_meaN <- rbind.data.frame(pLp_P_limb$moyFRic, pLp_P_leaves$moyFRic, pLp_P_fRoots$moyFRic)
                                Lp_meaN$organes <- as.factor(rep(c("a.limb", "b.leaves", "e.fine roots"), each = 4))
                                E <- ggplot(data = Lp_meaN, aes(x = div, y = x))+
                                  geom_point(aes(shape = div), size = 4)+
                                  geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
                                  scale_shape_discrete(guide = FALSE)+
                                  ylab("P content (mg g-1)")   +
                                  theme(panel.background = element_rect(fill="white", colour = "black"))+
                                  facet_grid(.~organes)
                                E
   
                                #-----------------------------------------------------------
                                # 3.4 Effet des traitements à l'échelle spécifique:   
                                #        quantité de N dans les organes de colza, ray-grass et trèfle
                                #------------------------------------------------------------                             
                                
                                Pcont$qN <- Pcont$DM*Pcont$Ncont/100
                                qttN <- aggregate(qN~distance+RT+leg+div+id_taxo, sum, data = Pcont)
                                myqttN <- myprep(qttN)
                                
                                #############################################        
                                ###           Pour les colzas             ###
                                #############################################
   Bn_qN <- myqttN[myqttN$id_taxo == "Bra_napus" , ]
    moyFi <- aggregate(Bn_qN$qN, by = list(div = Bn_qN$div, leg = Bn_qN$leg), mean, na.rm=T)
    seFi <- aggregate(Bn_qN$qN, by = list(div = Bn_qN$div, leg = Bn_qN$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("qN") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    #############################################        
    ###           Pour les ray-grass          ###
    #############################################
    
    Lp_qN <- myqttN[myqttN$id_taxo == "Lol_peren" , ]
    moyFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div, leg = Ta_qN$leg), mean, na.rm=T)
    seFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div, leg = Ta_qN$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("qN") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    #############################################        
    ###           Pour les trèfles            ###
    #############################################
    Ta_qN <- myqttN[myqttN$id_taxo == "Tri_alexa" , ]
    moyFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div), mean, na.rm=T)
    seFi <- aggregate(Ta_qN$qN, by = list(div = Ta_qN$div), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = div, y=x)) + 
      geom_bar(aes(fill = div), position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("qN") +    
      xlab("") +
      scale_fill_manual (values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    
    
    #############################################        
    ###           Bilan par RT                ###
    #############################################
    moyFi_0 <- aggregate(qttN$qN, by = list(div = qttN$div, leg = qttN$leg, RT = qttN$RT), sum, na.rm=T)
    moyFi <- aggregate(moyFi_0$x, by = list(div = moyFi_0$div, leg = moyFi_0$leg), mean, na.rm=T)
    seFi <- aggregate(moyFi_0$x, by = list(div = moyFi_0$div, leg = moyFi_0$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("Production de N totale par RT (mg)") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    #-----------------------------------------------------------
    # 3.5 Effet des traitements à l'échelle spécifique:   
    #        quantité de P dans les organes de colza, ray-grass et trèfle
    #------------------------------------------------------------  
    Pcont$qP <- Pcont$DM*Pcont$Pcont/100
    qttP <- aggregate(qP~distance+RT+leg+div+id_taxo, sum, data = Pcont)
    myqttP <- myprep(qttP)
    
   
    #############################################        
    ###           Pour les colzas             ###
    #############################################
    Bn_qP <- myqttP[myqttP$id_taxo == "Bra_napus" , ]
    moyFi <- aggregate(Bn_qP$qP, by = list(div = Bn_qP$div, leg = Bn_qP$leg), mean, na.rm=T)
    seFi <- aggregate(Bn_qP$qP, by = list(div = Bn_qP$div, leg = Bn_qP$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("qP") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    
    #############################################        
    ###           Pour les ray-grass          ###
    #############################################
    Lp_qP <- myqttP[myqttP$id_taxo == "Lol_peren" , ]
    moyFi <- aggregate(Lp_qP$qP, by = list(div = Lp_qP$div, leg = Lp_qP$leg), mean, na.rm=T)
    seFi <- aggregate(Lp_qP$qP, by = list(div = Lp_qP$div, leg = Lp_qP$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("qP") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    #############################################        
    ###           Pour les trèfles            ###
    #############################################
    Ta_qP <- myqttP[myqttP$id_taxo == "Tri_alexa" , ]
    moyFi <- aggregate(Ta_qP$qP, by = list(div = Ta_qP$div), mean, na.rm=T)
    seFi <- aggregate(Ta_qP$qP, by = list(div = Ta_qP$div), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = div, y=x)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("qP") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER

    
    #-----------------------------------------------------------
    # 3.6 Effet des traitements à l'échelle spécifique:   
    #        rapport N/P dans les organes de colza, ray-grass et trèfle
    #------------------------------------------------------------  
    myqttN$code <- paste(myqttN$RT, myqttN$distance, myqttN$leg, myqttN$div, myqttN$id_taxo, myqttN$serie, myqttN$position, myqttN$PAR)
    myqttP$code <- paste(myqttP$RT, myqttP$distance, myqttP$leg, myqttP$div, myqttP$id_taxo, myqttP$serie, myqttP$position, myqttP$PAR)
    
    myqttNP <- merge(myqttN, myqttP, by = "code")
    myqttNP <- myqttNP[, c(2:7, 16)]
    colnames(myqttNP) <- c("RT",  "distance", "leg","div", "id_taxo", "qN", "qP" )       
    myqttNP$rNP <- myqttNP$qN/myqttNP$qP
    
    #############################################        
    ###           Pour les colzas             ###
    #############################################
    Bn_rNP <- myqttNP[myqttNP$id_taxo == "Bra_napus" , ]
    moyFi <- aggregate(Bn_rNP$rNP, by = list(div = Bn_rNP$div, leg = Bn_rNP$leg), mean, na.rm=T)
    seFi <- aggregate(Bn_rNP$rNP, by = list(div = Bn_rNP$div, leg = Bn_rNP$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("rNP") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    
    #############################################        
    ###           Pour les ray-grass          ###
    #############################################
    Lp_rNP <- myqttNP[myqttNP$id_taxo == "Lol_peren" , ]
    moyFi <- aggregate(Lp_rNP$rNP, by = list(div = Lp_rNP$div, leg = Lp_rNP$leg), mean, na.rm=T)
    seFi <- aggregate(Lp_rNP$rNP, by = list(div = Lp_rNP$div, leg = Lp_rNP$leg), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("rNP") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    #############################################        
    ###           Pour les trèfles            ###
    #############################################
    Ta_rNP <- myqttNP[myqttNP$id_taxo == "Tri_alexa" , ]
    moyFi <- aggregate(Ta_rNP$rNP, by = list(div = Ta_rNP$div), mean, na.rm=T)
    seFi <- aggregate(Ta_rNP$rNP, by = list(div = Ta_rNP$div), std.error, na.rm=T)
    moyFi$se <- seFi$x
    DIV_INTER <- ggplot(data = moyFi, aes(x = div, y=x)) + 
      geom_bar(position = position_dodge(width = 0.8), stat="identity")+
      geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
      ylab("rNP") +    
      xlab("") +
      scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                        name = "",
                        labels = c("Ew0", "FR1", "FR2", "FR3"))+
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
            text =  element_text(face = "plain",
                                 color = "black", 
                                 hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                                 margin = margin(), debug = FALSE))+
      theme(axis.text = element_text(color="black", size=14)) +
      theme(axis.title = element_text(color="black", size=14))+
      scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
    
    DIV_INTER
    
    
    
    
    
    
    #-----------------------------------------------------------
    # 4.1 Relations entre traits (moyenne par rhizotron)
    #------------------------------------------------------------                
    # Rappel des df créées plus haut
    # myvdt -> biom vdt
    # myburrow <- réseau de galeries des vdt
    # mybiom_Com -> biom DW communauté végétale          
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
    
    # pour les organes aériens
    summary(lm(Bn~Ta, data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div0"),])) #
    summary(lm(Bn~Ta, data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div1"),])) # R^2 = 0.44*
    summary(lm(sqrt(Bn)~sqrt(Ta), data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div2"),])) # R^2 = 0.18$
    summary(lm(Bn~Ta, data = k2[which(k2$organes == c("Leaves", "Stem") & k2$div == "Div3"),])) # 
    
    # pour les racines fines
    summary(lm(Bn~Ta, data = k2[which(k2$organes == "fRoots" & k2$div == "Div0"),])) #
    summary(lm(sqrt(Bn)~sqrt(Ta), data = k2[which(k2$organes == "fRoots" & k2$div == "Div1") ,]) # 
            summary(lm(sqrt(Bn)~sqrt(Ta), data = k2[which(k2$organes == "fRoots" & k2$div == "Div2"),])) # R^2 = 0.65$ si l'outlier est supprimé (RT 42)
            summary(lm(Bn~Ta, data = k2[which(k2$organes == "fRoots" & k2$div == "Div3"),])) # 
            
            # Représentations graphiques
            Reg1 <- ggplot(data = k2[k2$organes == c("Leaves", "Stem"),], aes(y = Bn, x = Ta))+
              geom_point(aes(shape = div))+
              geom_smooth(method = "lm")+
              facet_grid(div~.)+
              ylab("B. napus\n Shoot organs N content (mg g-1)")+
              xlab("T. alexandrium\n Shoot organs N content (mg g-1)")+
              theme(panel.background = element_rect(fill="white", colour = "black"))+
              scale_shape_discrete(guide = FALSE)
            
            Reg2 <- ggplot(data = k2[which(k2$organes == "fRoots"),], aes(y = Bn, x = Ta))+ # voir pourquoi il y a un Ta à 4 mg g-1 (RT20)??
              geom_point(aes(shape = div))+
              geom_smooth(method = "lm")+
              facet_grid(div~.)+
              ylab("B. napus\n fine roots N content (mg g-1)")+
              xlab("T. alexandrium\n fine roots N content (mg g-1)")+
              theme(panel.background = element_rect(fill="white", colour = "black"))+
              scale_shape_discrete(guide = FALSE)
            
            
            
            #-----------------------------------------------------------
            # 5.1 Liste des figures à proposer
            #------------------------------------------------------------    
            
            # Figure Racine
            grid.arrange(pTa_A_DW$pDIV, pBn_A_all_DW$pDIV, pLp_A_all_DW$pDIV,
                        Ta_PR, Bn_PR, Lp_PR, 
                         ncol = 3, nrow = 2)
            # Figure activité enzymatiques
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

