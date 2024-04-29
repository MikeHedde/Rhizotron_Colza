################################################################


### analyse maximum rooting depth RHIZOTRON

# MZ 
# 8 mars 2016

############################################################"#

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/root_trait/")

tab <- read.csv("RLMax.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
head(tab)

## on vérifie les profondeurs d'enracinement dans la serie 1 et 2
tab$serie <- NULL
tab$serie[tab$id_RT<21]="s1"
tab$serie[tab$id_RT>20]="s2"

par(mfrow=c(1,2))
hist(tab$RLMax[tab$serie=="s1"])
hist(tab$RLMax[tab$serie=="s2"])

# vérification des données entre séries 1 et 2
#-----------------------------------------------
tab$serie <- as.factor(tab$serie)
kruskal.test(RLMax~serie, data=tab, na.action=na.omit)
par(mfrow=c(1,1))
boxplot(tab$RLMax~tab$serie)
# P < 0.01
mean(tab$RLMax[tab$serie=="s1"])
mean(tab$RLMax[tab$serie=="s2"])

bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]
gd <- tab[tab$id_taxo=="Ger_disse",]
vp <- tab[tab$id_taxo=="Ver_persi",]

par (mfrow=c(2,3))
data = vp
main = "Ver_persi"
boxplot(RLMax~serie, data=data, main=main)


kruskal.test(RLMax~serie, data=vp, na.action=na.omit)
# NS sur Bn (P=0.5121), lp (P = 0.08697), 
# significatif sur Ta (0=005593), gd (P< 0.0001), vp (P< 0.0001)
mean(ta$RLMax[ta$serie=="s1"]) 
mean(ta$RLMax[ta$serie=="s2"]) 

#  effet des traitements sur la profondeur max totale ?
#-----------------------------------------------------------
# on fait un premier essai avec tous les individus
par(mfrow=c(1,2))
plot(tab$RLMax~tab$id_leg, xlab="", ylab="Maximum rooting depth (cm)")
plot(tab$RLMax~tab$id_div, xlab="", ylab="Maximum rooting depth (cm)")

tab[tab$RLMax>150,]
# on supprime le point pour le moment
tab$RLMax[tab$RLMax>150]=NA

tab$id_div <- as.factor(tab$id_div)
tab$id_leg <- as.factor(tab$id_leg)
kruskal.test(RLMax~id_div, data=tab, na.action=na.omit)
kruskal.test(RLMax~id_leg, data=tab, na.action=na.omit)

mean(tab$RLMax[tab$id_leg=="Leg+"])
mean(tab$RLMax[tab$id_leg=="Leg-"], na.rm=TRUE)

# modele mixte
hist(tab$RLMax)
library(lme4)

tab$tmt <- paste(div, leg, sep="_")
# essai 1
# mm1 <- glmer(log(RLMax) ~ -1 + leg * div + (1|serie), data=tab, family= gaussian)
# summary(mm1)
# shapiro.test(residuals(mm1))
# anova(mm1)
# essai 2
# library(nlme)
# mm2 <- lme (log(RLMax)~div*leg,  random=~1|serie/div, data=tab, na.action=na.omit)
# shapiro.test(residuals(mm2))
# bartlett.test(log(RLMax)~div,  data=tab, na.action=na.omit)
# bartlett.test(log(RLMax)~leg,  data=tab, na.action=na.omit)
# anova(mm2)



##comparaisons multiples après GLMM (pas conseillée d'après David M)
library(multcomp)
mm1_glht <- glht(mm1, linfct = mcp(id_div = "Tukey"))
summary(mm1_glht)

##calculs des IC du GLMM (conseillée d'après David M)
mm1W <- confint(mm1, parm = "beta_", method="Wald")
mm1tab <- cbind(est=fixef(mm1), mm1W)
rtab <- as.data.frame(exp(mm1tab))
print(rtab,digits=3)

#  effet des traitements sur la profondeur max totale de chaque sp?
#-----------------------------------------------------------
bn <- tab[tab$id_taxo=="Bra_napus",]
lp <- tab[tab$id_taxo=="Lol_peren",]
ta <- tab[tab$id_taxo=="Tri_alexa",]
gd <- tab[tab$id_taxo=="Ger_disse",]
vp <- tab[tab$id_taxo=="Ver_persi",]

data = gd

mm1 <- glmer(log(RLMax) ~ -1 + div*leg + (1|serie), data=data, family= poisson)
shapiro.test(residuals(mm1))
bartlett.test(RLMax ~div, data=data, na.action=na.omit)
bartlett.test(RLMax ~leg, data=data, na.action=na.omit)
anova(mm1)
mm1_glht <- glht(mm1, linfct = mcp(div = "Tukey"))
summary(mm1_glht)

kruskal.test(RLMax ~div, data=data, na.action=na.omit)
kruskal.test(RLMax ~leg, data=data, na.action=na.omit)

plot(RLMax~div, data=data, xlab="", ylab="Maximum rooting depth (cm)")
plot(RLMax~leg, data=data, xlab="", ylab="Maximum rooting depth (cm)")

mm2 <- lme (log(RLMax)~div*leg,  random=~1|serie/div, data=data, na.action=na.omit)
shapiro.test(residuals(mm2))
bartlett.test(log(RLMax)~div,  data=data, na.action=na.omit)
bartlett.test(log(RLMax)~leg,  data=data, na.action=na.omit)
anova(mm2)

anova(mm1,mm2)

# essai ACP pour voir la dsitribution des données
#---------------------------------------------------
setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/DATA/DW/essai_ACP/")

######librairies#######
require(gplots)
require(vegan)
require(ade4)
require(FD)
require(entropart)
require(ggplot2)
require(kruskalmc)
require(reshape2)


######données brutes####
brac <-read.csv("brac.csv", h=T, sep = ";")
brac2 <-read.csv("brac2.csv", h=T, sep = ";")

# répartition des masses racinaires par strate et pour chaque loc 

tab_acp <- read.csv("Table_ACP_RAC.csv", header=TRUE, sep=";", dec=".")
colnames(tab_acp) <- c("RT","strata","div","leg","-24","-21","-18","-15","-12","-9","-6","-3","0","3","6","9","12","15","18","21","24")
head(tab_acp)
s1 <- tab_acp[tab_acp$RT<21,]
s2 <- tab_acp[tab_acp$RT>20,]

par(mfrow = c(2,3))

# série 1
pca1 <- dudi.pca(s1[,-c(1:4)], scannf=F, nf=3)
# on regarde les eigenvalmues:
pca1$eig
sumI <- sum(pca1$eig)
pve <- 100 * pca1$eig/sum(pca1$eig)
# on trace les projections
s.corcircle (pca1$co)
s.class(pca1$li, as.factor(s1$div))
s.class(pca1$li, as.factor(s1$strata))

# série 2
pca2 <- dudi.pca(s2[,-c(1:4)], scannf=F, nf=3)
# on regarde les eigenvalmues:
pca2$eig
sumI <- sum(pca2$eig)
pve2 <- 100 * pca2$eig/sum(pca2$eig)
# on trace les projections
s.corcircle (pca2$co)
s.class(pca2$li, as.factor(s2$div))
s.class(pca2$li, as.factor(s2$strata))

# séries 1 et 2
par(mfrow = c(1,3))
pca3 <- dudi.pca(tab_acp[,-c(1:4)], scannf=F, nf=3)
# on regarde les eigenvalmues:
pca3$eig
sumI <- sum(pca3$eig)
pve3 <- 100 * pca3$eig/sum(pca2$eig)
# on trace les projections
s.corcircle (pca3$co)
s.class(pca3$li, as.factor(tab_acp$div))
s.class(pca3$li, as.factor(tab_acp$strata))



