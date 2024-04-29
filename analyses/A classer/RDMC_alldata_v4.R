################################################################


### analyse RDMC RHIZOTRON

# MZ 
# 18 avril 2016
# 21 avril 2016 : reprends les analyses après vérifcation des données DW et FW + calcul des FW et DW manquants

############################################################"#

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/root_trait/")

tab <- read.csv("RDMC_alldata_v4.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
head(tab)


### vérification des poids echantillons
#---------------------------------------------------
tab$serie <- NULL
tab$serie[tab$id_RT<21]="s1"
tab$serie[tab$id_RT>20]="s2"

par(mfrow=c(1,2))
hist(tab$RDMC[tab$serie=="s1"])
hist(tab$RDMC[tab$serie=="s2"])

### vérification des données manquantes
summary(tab)
tab[is.na(tab$RDMC)==TRUE,]
# d'abord on vérifie les points NA
# puis on recalcule les FW et DW manquants à partir des relations FW/DW
bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]
gd <- tab[tab$id_taxo=="Ger_disse",]
vp <- tab[tab$id_taxo=="Ver_persi",]

par(mfrow=c(1,1))
plot(DW~FW, data=bn)
reg_bn<- lm(DW~-1+FW, data=bn)
summary(reg_bn)
abline(reg_bn)
bn[bn$FW>2000 & bn$DW<100,]
# on recalclue le FW pour RT16_d+15_P2   et  RT20_d0_R4
bn[bn$FW>5000 & bn$DW<400,]
# on recalclue le FW pour RT14_d0_R4  (et  RT27_d+15_R4)
# on considère que DW est correcte, mais incertitude sur mesure de FW (racines mal séchées avant pesées et donc surestimation du poids)
# DW = 0.102169*FW, donc FW = DW/0.102169

#
plot(DW~FW, data=lp)
reg_lp <- lm(DW~-1+FW, data=lp)
summary(reg_lp)
abline(reg_lp)
lp[lp$FW>1000 & lp$DW<100,]
# on recalclue le FW pour RT17_d+9_R1 
lp[lp$FW>1500 & lp$DW>150,]
# on recalclue le FW pour RT17_d+9_R2
lp[lp$FW<250 & lp$DW>50,]
# on recalclue le FW pour RT10_d+9_R1 et RT10_d+9_R2 
# DW = 0.12357*FW, donc FW = DW/0.12357

# 
plot(DW~FW, data=ta)
reg_ta<- lm(DW~-1+FW, data=ta)
summary(reg_ta)
abline(reg_ta)
ta[ta$FW>1000 & ta$DW>150,]
# vérifier poids TARE pour RT42_d+21_R2
# tare OK --> outlier mais poids corrects
ta[ta$FW<100 & ta$DW>20,]
# on recalclue le FW pour RT18_d+6_R3 et RT18_d-3_R2 
# DW = 0.105761*FW, donc FW = DW/0.105761

#
plot(DW~FW, data=gd)
reg_gd <- lm(DW~-1+FW, data=gd)
summary(reg_gd)
abline(reg_gd)
gd[gd$FW<20 & gd$DW>5,]
gd[gd$FW>45 & gd$DW<2,]
gd[gd$FW>100 & gd$DW<5,]
# DW = 0.109599*FW 



# on corrige la base de données, puis on relance à la requête pour calcul RDMC



###############################################################################################

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/root_trait/")

tab <- read.csv("RDMC_alldata_v5.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
head(tab)



# vérification des données entre séries 1 et 2
#-----------------------------------------------
tab$serie <- NULL
tab$serie[tab$id_RT<21]="s1"
tab$serie[tab$id_RT>20]="s2"
par(mfrow=c(1,1))
boxplot(tab$RDMC~tab$serie)
tab[tab$RDMC >800,]
tab[tab$RDMC >400 & tab$RDMC < 800,]

tab$serie <- as.factor(tab$serie)
kruskal.test(RDMC~serie, data=tab, na.action=na.omit)

# P < 0.05
mean(tab$RDMC[tab$serie=="s1"], na.rm=TRUE)
mean(tab$RDMC[tab$serie=="s2"], na.rm=TRUE)

tab[tab$RDMC>8000,]

bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]
gd <- tab[tab$id_taxo=="Ger_disse",]
vp <- tab[tab$id_taxo=="Ver_persi",]

par (mfrow=c(2,3))
data = vp
main = "Ver_persi"
boxplot(RDMC~serie, data=data, main=main)


kruskal.test(RDMC~serie, data=vp, na.action=na.omit)
# Bn */lp ***/ ta*** / gd *** / vp ***
mean(ta$RLMax[ta$serie=="s1"]) 
mean(ta$RLMax[ta$serie=="s2"]) 

# on vérifie les points suspects :
bn[is.na(bn$RDMC)==FALSE & bn$RDMC>500,]

# on corriges les points suspects : tous les poids secs ont été vérifiés... 
# on recalcule les FW pour corriger RDMC bizarres
bn[bn$id_loca=="RT04_d-15"&bn$division=="Taproot"&bn$strata==2,]
bn[bn$id_loca=="RT13_d+15"&bn$division=="Fine roots"&bn$strata==1,]
bn[bn$id_loca=="RT13_d-15"&bn$division=="Fine roots"&bn$strata==1,]

ta[ta$id_loca=="RT18_d-3"&ta$division=="Fine roots"&ta$strata==3,]

gd[gd$id_loca=="RT24_d+3"&gd$division=="Fine roots"&gd$strata==1,]
gd[gd$id_loca=="RT24_d+3"&gd$division=="Fine roots"&gd$strata==2,]


#################################################################################################

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/root_trait/")

tab <- read.csv("RDMC_alldata_v6.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
head(tab)

bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]
gd <- tab[tab$id_taxo=="Ger_disse",]
vp <- tab[tab$id_taxo=="Ver_persi",]

hist(bn$RDMC)
bn[bn$RDMC<50,]

id_loca id_div id_leg division strata   id_taxo    DW    FW      RDMC id_RT
1041 RT25_d+15   Div2   Leg-  Taproot      1 Bra_napus  78.4 143.8  545.2017    25
1119 RT26_d-15   Div3   Leg-  Taproot      1 Bra_napus  92.2 154.4  597.1503    26
1750 RT39_d+15   Div2   Leg+  Taproot      1 Bra_napus 189.2 171.1 1105.7861    39


hist(lp$RDMC)
lp[lp$RDMC<50,]
664   RT17_d+9   Div2   Leg+ Fine roots      4 Lol_peren 3.00000e-01  9.2 3.260870e+01    17
866  RT21_d-12   Div2   Leg- Fine roots      2 Lol_peren 2.22045e-13  0.8 2.775560e-10    21
906   RT22_d+9   Div3   Leg- Fine roots      4 Lol_peren 7.00000e-01 14.1 4.964539e+01    22
929  RT22_d-24   Div3   Leg- Fine roots      4 Lol_peren 5.00000e-01 12.3 4.065041e+01    22
1260  RT29_d+9   Div2   Leg- Fine roots      3 Lol_peren 4.00000e-01 10.8 3.703704e+01    29
1723 RT38_d-12   Div1   Leg+ Fine roots      4 Lol_peren 2.00000e-01 13.3 1.503759e+01    38
1900 RT41_d-24   Div1   Leg- Fine roots      4 Lol_peren 6.00000e-01 14.6 4.109589e+01    41


hist(ta$RDMC)
ta[ta$RDMC<50,]
ta[ta$RDMC>300,]

# pour le moment, on élimine les points > RDMC 600 mg/g
tab2 <- tab[tab$RDMC < 600,]
r <- tab2[tab2$division=="Fine roots",]

r$strata <- as.factor(r$strata)
mod <- lm(RDMC~strata, data=r, na.action=na.omit)
anova(mod)
shapiro.test(residuals(mod))
kruskal.test(RDMC~strata, data=r, na.action=na.omit)
par(mfrow=c(1,1))
plot(RDMC~strata, data=r)


bn <- r[r$id_taxo=="Bra_napus",]
lp <- r[r$id_taxo=="Lol_peren",]
ta <- r[r$id_taxo=="Tri_alexa",]
gd <- r[r$id_taxo=="Ger_disse",]
vp <- r[r$id_taxo=="Ver_persi",]

par(mfrow=c(2,3))
plot(RDMC~strata, data=bn, main= "Bra_napus")
plot(RDMC~strata, data=lp, main= "Lol_peren")
plot(RDMC~strata, data=ta, main= "Tri_alexa")
plot(RDMC~strata, data=gd, main= "Ger_disse")
plot(RDMC~strata, data=vp, main= "Ver_persi")

library(nlme)

# effet strate sur RDMC BN
bn$strata <- as.factor (bn$strata)
mod_bn <- lm(RDMC~strata,data=bn)
shapiro.test(residuals(mod_bn))
anova (mod_bn)
kruskal.test(RDMC~strata, data=bn, na.action=na.omit)
mod_bn <- lme(RDMC~strata, random=~1|serie,data=bn)
shapiro.test(residuals(mod_bn))
anova (mod_bn)

# effet strate sur RDMC LP
lp$strata <- as.factor (lp$strata)
mod_lp <- lm(RDMC~strata,data=lp)
shapiro.test(residuals(mod_lp))
anova (mod_bn)
kruskal.test(RDMC~strata, data=lp, na.action=na.omit)

# effet strate sur RDMC TA
ta$strata <- as.factor (ta$strata)
mod_ta <- lm(RDMC~strata,data=ta)
shapiro.test(residuals(mod_ta))
anova (mod_ta)
kruskal.test(RDMC~strata, data=ta, na.action=na.omit)

# effet strate sur RDMC GD
gd <- gd[gd$strata!=4,]
gd$strata <- as.factor (gd$strata)
mod_gd <- lm(RDMC~strata,data=gd)
shapiro.test(residuals(mod_gd))
anova (mod_gd)
kruskal.test(RDMC~strata, data=gd, na.action=na.omit)

# effet strate sur RDMC VP
# vp <- vp[vp$strata!=4,]
gd$strata <- as.factor (gd$strata)
mod_vp <- lm(RDMC~strata,data=vp)
shapiro.test(residuals(mod_vp))
anova (mod_vp)
kruskal.test(RDMC~strata, data=vp, na.action=na.omit)

### effet des traitements
#-----------------------------
par(mfrow=c(1,1))
plot(RDMC~id_div, data=r)

tab2$div <- as.factor(tab2$div)
mod <- lm(RDMC~div*leg, data= tab2, na.action=na.omit)
shapiro.test(residuals(mod))
kruskal.test(RDMC~div, data=tab2, na.action=na.omit)

par(mfrow=c(2,3))
plot(RDMC~id_div, data=bn, main= "Bra_napus")
plot(RDMC~id_div, data=lp, main= "Lol_peren")
plot(RDMC~id_div, data=ta, main= "Tri_alexa")
plot(RDMC~id_div, data=gd, main= "Ger_disse")
plot(RDMC~id_div, data=vp, main= "Ver_persi")

# effet diversité sur BN
bn$id_div <- as.factor(bn$id_div)
mod_bn <- lm(RDMC~id_div,data=bn)
shapiro.test(residuals(mod_bn))
anova (mod_bn)
kruskal.test(RDMC~id_div, data=bn, na.action=na.omit)
mod_bn <- lme(RDMC~id_div, random=~1|serie,data=bn)
shapiro.test(residuals(mod_bn))
anova (mod_bn)
bn_r4 <- bn[bn$strata==4,]
mod_bn_r4 <- lm(log(RDMC)~id_div,data=bn_r4)
shapiro.test(residuals(mod_bn_r4))
anova (mod_bn_r4)
kruskal.test(RDMC~id_div, data=bn_r4, na.action=na.omit)
bn_r3 <- bn[bn$strata==3,]
mod_bn_r3 <- lm(log(RDMC)~id_div,data=bn_r3)
shapiro.test(residuals(mod_bn_r3))
anova (mod_bn_r3)
kruskal.test(RDMC~id_div, data=bn_r3, na.action=na.omit)
bn_r2 <- bn[bn$strata==2,]
mod_bn_r2 <- lm(log(RDMC)~id_div,data=bn_r2)
shapiro.test(residuals(mod_bn_r2))
anova (mod_bn_r2)
kruskal.test(RDMC~id_div, data=bn_r2, na.action=na.omit)
bn_r1 <- bn[bn$strata==1,]
mod_bn_r1 <- lm(log(RDMC)~id_div,data=bn_r1)
shapiro.test(residuals(mod_bn_r1))
anova (mod_bn_r1)
kruskal.test(RDMC~id_div, data=bn_r1, na.action=na.omit)
par(mfrow=c(2,2))
plot(RDMC~id_div, data=bn_r1, main=" BN strate 1")
plot(RDMC~id_div, data=bn_r2, main=" BN strate 2")
plot(RDMC~id_div, data=bn_r3, main=" BN strate 3")
plot(RDMC~id_div, data=bn_r4, main=" BN strate 4")
# sur BN 0 uniquement
bn$id_loca <- as.character(bn$id_loca)
dist <- strsplit(x=bn$id_loca, split="_")
dist <- unlist(dist)
z <- seq(from=2, to =length(dist), by=2)
bn$distance <- dist[z]

library(nlme)
bn0 <- bn[bn$distance=="d0",]
hist(bn0$RDMC)
mod_bn0 <- lm(log(RDMC)~id_leg*id_div, data=bn0, na.action=na.omit) # , random=~1|serie
shapiro.test(residuals(mod_bn0))
fligner.test(log(RDMC)~id_leg, data=bn0, na.action=na.omit)
fligner.test(log(RDMC)~id_div, data=bn0, na.action=na.omit)
anova(mod_bn0)

# pour faire une anova de type III
library(car)
#  using Anova() {car} in combination with options(contrasts=c("contr.sum", "contr.poly")) in R gives the same Type-III ANOVA tables as calculated in SPSS.
Anova(mod_bn0, type="III")

# on calcule un RDMC moye par plante
bn_ind <- aggregate (bn$RDMC, by=list(id_loca=bn$id_loca, id_RT=bn$id_RT, id_leg=bn$id_leg, id_div=bn$id_div, id_taxo=bn$id_taxo, serie=bn$serie, distance=bn$distance, division=bn$division), mean, na.rm=TRUE)
colnames(bn_ind) <- c("id_loca",)