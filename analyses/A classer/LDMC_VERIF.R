###############################################################""

### analyse RDMC RHIZOTRON

# 21 avril 2016 : vérification des données 


###############################################################""

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/leaf_trait/")

tab <- read.csv("LDMC_all_data_v3.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
head(tab)
summary(tab)

bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]

par(mfrow=c(1,1))
plot(DW~FW, data=bn)
reg_bn<- lm(DW~-1+FW, data=bn)
summary(reg_bn)
abline(reg_bn)
bn[bn$FW>1500 & bn$DW<100,]
# RT10_d0
# on retire RT 10 pour calculer la regression 
bn[bn$id_loca=="RT10_d0",c("DW","FW","DMC")]=NA
plot(DW~FW, data=bn)
reg_bn<- lm(DW~-1+FW, data=bn)
summary(reg_bn)
abline(reg_bn)

#
plot(DW~FW, data=lp)
reg_lp <- lm(DW~-1+FW, data=lp)
summary(reg_lp)
abline(reg_lp)
lp[lp$FW<50 & lp$DW>19,]
lp[lp$FW<40 & lp$DW>10,]
lp[lp$DW>40,]


# 
plot(DW~FW, data=ta)
reg_ta<- lm(DW~-1+FW, data=ta)
summary(reg_ta)
abline(reg_ta)


## on vérifie les points suspects --> repesée des échantillons et calculs des FW
bn[bn$id_loca=="RT10_d0",]



##############################################################################################

# analyses des données
#########################"
###############################################################""

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/leaf_trait/")

tab <- read.csv("LDMC_all_data_v3.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
head(tab)
summary(tab)
tab$serie <- NULL
tab$serie[tab$id_RT<21]="s1"
tab$serie[tab$id_RT>20]="s2"
# on crée une colonne distance
tab$id_loca <- as.character(tab$id_loca)
dist <- strsplit(x=tab$id_loca, split="_")
dist <- unlist(dist)
z <- seq(from=2, to =length(dist), by=2)
tab$distance <- dist[z]

bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]
summary(bn)
summary(lp)
lp[is.na(lp$LDMC)==TRUE,] 
# RT41 _LP+21 manquant !!!!
summary(ta)

hist(bn$LDMC)
hist(lp$LDMC)
lp[lp$LDMC>350,] 

hist(ta$LDMC)
lp[lp$LDMC<150,] 


##############################################################################################
# vérification des petites éh < 15 mg

pt_ech <- tab[is.na(tab$DW)==FALSE & tab$DW < 15,]
# on crée une colonne distance
pt_ech$id_loca <- as.character(pt_ech$id_loca)
dist <- strsplit(x=pt_ech$id_loca, split="_")
dist <- unlist(dist)
z <- seq(from=2, to =length(dist), by=2)
pt_ech$distance <- dist[z]
pt_ech

nrow(bn)
pt_bn <- pt_ech[pt_ech$id_taxo=="Bra_napus",]
nrow(pt_bn)
nrow(pt_bn)/nrow(bn)
#
nrow(lp)
pt_lp <- pt_ech[pt_ech$id_taxo=="Lol_peren",]
nrow(pt_lp)
nrow(pt_lp)/nrow(lp)
#
nrow(ta)
pt_ta <- pt_ech[pt_ech$id_taxo=="Tri_alexa",]
nrow(pt_ta)
nrow(pt_ta)/nrow(lp)


################################################################################################

library(nlme)
library(car)

# analyse de données pour BN
#####################################
boxplot (bn$LDMC~bn$serie)
bn[bn$LDMC>300,]
# on exclue les bordures
bn0 <- bn[bn$distance=="d0",]
boxplot (bn0$LDMC~bn0$serie)
hist(bn0$LDMC)
# effet des traitements sur LDMC
mod1 <- lme(LDMC~id_div*id_leg,random=~1|serie, data=bn0, na.action=na.omit)
shapiro.test(residuals(mod1))
bartlett.test(LDMC~id_div, data=bn0, na.action=na.omit)
bartlett.test(LDMC~id_leg, data=bn0, na.action=na.omit)
anova(mod1)
# pour faire une anova de type III
Anova(mod1, type="III")

boxplot(LDMC~paste(id_leg, id_div, sep="_"), data=bn0, ylab="LDMC (mg/g DM", main="Bn0")


## choix du modele
######################
# modele linaire mixte
mod1 <- lme(LDMC~id_div+id_leg+id_div:id_leg,random=~1|serie, data=bn0, na.action=na.omit)
summary(mod1)
Anova(mod1, type="III")
# modele polynomial
mod2 <- lme(LDMC~id_div+id_leg+id_div^2+id_leg^2+id_div:id_leg,random=~1|serie, data=bn0, na.action=na.omit)
summary(mod2)
Anova(mod2, type="III")


# analyse de données pour LP
#####################################
boxplot (lp$LDMC~lp$serie)
lp[lp$LDMC>350,]
# id_loca id_div id_leg division strata   id_taxo  DW   FW     LDMC id_RT serie
# RT10_d-12   Div1   Leg-    Leaf3      0 Lol_peren 9.9 28.2 351.0638    10    s1
# on exclue les bordures
lp0 <- lp[lp$distance=="d+9"|lp$distance=="d-12"|lp$distance=="d-3",]
boxplot (lp0$LDMC~lp0$serie)
hist(lp0$LDMC)
# effet des traitements sur LDMC
mod1 <- lme(LDMC~id_leg*id_div,random=~1|serie, data=lp0, na.action=na.omit)
shapiro.test(residuals(mod1))
bartlett.test(LDMC~id_div, data=lp0, na.action=na.omit)
bartlett.test(LDMC~id_leg, data=lp0, na.action=na.omit)
anova(mod1)
# pour faire une anova de type III
Anova(mod1, type="III")

boxplot(LDMC~paste(id_leg, id_div, sep="_"), data=lp0, ylab="LDMC (mg/g DM", main="Lp (-12/-3/+9)")
boxplot(LDMC~id_leg, data=lp0, ylab="LDMC (mg/g DM)", main="Lp (-12/-3/+9)")
boxplot(LDMC~id_leg, data=lp, ylab="LDMC (mg/g DM)", main="Lp")

# analyse de données pour TA
#####################################
boxplot (ta$LDMC~ta$serie)
ta[ta$LDMC<150,]
# id_loca id_div id_leg division strata   id_taxo  DW   FW     LDMC id_RT serie
# RT19_d-18   Div0   Leg+    Leaf3      0 Tri_alexa 4.7 34.4 136.6279    19    s1
# on exclue les bordures
ta0 <- ta[ta$distance=="d+12"|ta$distance=="d+6"|ta$distance=="d-3"|ta$distance=="d-9",]
boxplot (ta0$LDMC~ta0$serie)
hist(ta0$LDMC)
# effet des traitements sur LDMC
mod1 <- lme(LDMC~id_div,random=~1|serie, data=ta0, na.action=na.omit)
shapiro.test(residuals(mod1))
bartlett.test(LDMC~id_div, data=ta0, na.action=na.omit)
anova(mod1)
# pour faire une anova de type III
Anova(mod1, type="III")

boxplot(LDMC~id_div, data=ta0, ylab="LDMC (mg/g DM", main="Ta (-9/-3/+6/+12)")
boxplot(LDMC~id_div, data=lp, ylab="LDMC (mg/g DM)", main="Ta")


####################################################################### 
# on analyse à partir des moyennes par RT (bordures exclues)
moylp <- aggregate (lp0$LDMC, by=list(serie=lp0$serie, id_RT=lp0$id_RT, id_div=lp0$id_div, id_leg=lp0$id_leg), mean, na.rm=TRUE)
mod2 <- lme(x~id_leg*id_div,random=~1|serie, data=moylp, na.action = na.omit)
Anova(mod2, type="III")

moyta <- aggregate (ta0$LDMC, by=list(serie=ta0$serie, id_RT=ta0$id_RT, id_div=ta0$id_div, id_leg=ta0$id_leg), mean, na.rm=TRUE)
mod2 <- lme(x~id_div,random=~1|serie, data=moyta, na.action = na.omit)
Anova(mod2, type="III")
