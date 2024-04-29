#####################"

## analyses densité plante / Rhizotron


# MZ
# 28/09/2016 : reprend SCRIPT_BIOM_PL_20160928.r au propre
# 12/12:2016 : reprend stat après élimintaion des 4 RT douteux


#######################


#### Librairies ####
require(ggplot2)
require(vegan)
require(ade4)
require(reshape2)
require(multcomp)
require(car)
require(nlme)
require(plotrix)
require(PMCMR)

##Données
##################
setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/Données/Plantes/densité")

# nombre de plantes
###################
nb <- read.csv("DOC_RT_Plant_Nb.csv", sep=";", dec=".", header=T, na.strings=NA)
summary(nb)
# on ajoute la colonne série
nb$serie <- NULL
nb$serie[nb$RT<21]="s1"
nb$serie[nb$RT>20]="s2"
# on supprime les 4 RT douteux
nb$Pl_Nb[nb$RT==8|nb$RT==12|nb$RT==18|nb$RT==40]=NA



p <- ggplot(nb, aes(div, Pl_Nb)) +
  ylab("Nombre d'individus par rhizotron")+
  xlab("Earthworm functional diversity") +
  geom_boxplot() +
  facet_grid(id_taxo ~ leg)
p  

# Listing des RT avec outliers
nb[nb$id_taxo=="Bra_napus" & nb$Pl_Nb<3,]
# RT id_div id_leg   id_taxo Pl_Nb
#   5   Div0   Leg- Bra_napus     2
#   9   Div0   Leg- Bra_napus     2
#  16   Div1   Leg+ Bra_napus     2
#  17   Div2   Leg+ Bra_napus     2
#  36   Div3   Leg+ Bra_napus     2
nb[nb$id_taxo=="Lol_peren" & nb$leg=="Leg-" & nb$Pl_Nb<5,]
# RT id_div id_leg   id_taxo Pl_Nb
# 8   Div3   Leg- Lol_peren     4
nb[nb$id_taxo=="Tri_alexa" & nb$leg=="Leg+" & nb$Pl_Nb<5,]
# RT id_div id_leg   id_taxo Pl_Nb
#  14   Div3   Leg+ Tri_alexa     4
#  40   Div3   Leg+ Tri_alexa     4
#  42   Div2   Leg+ Tri_alexa     4


# on suave la table
write.table(nb, "DATA_RT_NB_plantes_SP.csv", sep=";", dec=".", col.names=T, row.names=F)


# stat : effet des traitements sur nb de plantes/esp
##-----------------------------------------------------
bn<- nb[nb$id_taxo=="Bra_napus",]
lp<- nb[nb$id_taxo=="Lol_peren",]
gd<- nb[nb$id_taxo=="Ger_disse",]
vp<- nb[nb$id_taxo=="Ver_persi",]

data = gd[gd$leg=="Leg-",]
data$serie <- as.factor(data$serie)
mod <- lm(log(Pl_Nb)~serie+div, data=data, na.action=na.omit)
shapiro.test(residuals(mod))
fligner.test(Pl_Nb~leg, data=data, na.action=na.omit)
fligner.test(Pl_Nb~div, data=data, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
# comparaison multiple
comp <- glht(mod, linfct=mcp(div = "Tukey"))
summary(comp)
# tests non paramétriques si conditions ANOVA non respectée (ou P-valie model NS)
kruskal.test(Pl_Nb~leg,data=data, na.action=na.omit)
kruskal.test(Pl_Nb~div,data=data, na.action=na.omit)
kruskal.test(Pl_Nb~serie,data=data, na.action=na.omit)
require(PMCMR)
posthoc.kruskal.nemenyi.test(Pl_Nb~div, data=data, dist="Tukey")


# nombre total de plantes par rhizotron
#--------------------------------------
nbtot <- aggregate(nb$Pl_Nb, by=list(RT=nb$RT, div=nb$div, leg=nb$leg, serie=nb$serie), sum, na.rm=T)
colnames(nbtot)=c("RT", "div", "leg", "serie","Nb")
summary(nbtot)
nbtot[nbtot$Nb==10,]

# on élimine les RT douteux
den <- nbtot
den$Nb[den$RT==12 | den$RT==8 | den$RT==18 |den$RT==40]=NA

# on sauve la table 
write.table(den, "DATA_RT_Plant_density.csv", sep=";", dec=".", col.names=T, row.names=F)

# stat : effet des traitement sur le nombre de plantes/RT
den$serie <- as.factor(den$serie)
mod <- lm(Nb~serie+leg*div, data=den, na.action=na.omit )
shapiro.test(residuals(mod))
fligner.test(Nb~leg,data=den, na.action=na.omit)
fligner.test(Nb~div,data=den, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
# test non paramétriques
kruskal.test(Nb~leg,data=den, na.action=na.omit)
kruskal.test(Nb~div,data=den, na.action=na.omit)
kruskal.test(Nb~serie,data=den, na.action=na.omit)

# visualisation des données
boxplot(den$Nb~den$serie)

# moyenne par traitement sans les RT douteux
nbmoy <- aggregate(den$Nb, by=list(div=den$div, leg=den$leg), mean, na.rm=T)
se <- aggregate(den$Nb, by=list(div=den$div, leg=den$leg), std.error, na.rm=T)
colnames(nbmoy)=c("div", "leg", "Nb")
nbmoy$se <- se$x


p <- ggplot(nbmoy, aes(x=leg, y=Nb, fill=div)) + 
  geom_bar(position=position_dodge(), stat="identity", 
           colour="black", # use black outlines, 
           size=0.3)+   # for thinner lines
  geom_errorbar(aes(ymin=Nb-se, ymax=Nb+se),
                size=0.3, # for thinner lines
                width=0.2,  # Width of the error bars
                position=position_dodge(.9))+
  ylab("Plant density")+
  xlab("Earthworm functional diversity gradient")+
  scale_y_continuous(breaks=0:27*1) +
  scale_fill_manual(values=c("white", "azure3", "azure4", "black")) + # pour changer les couleurs des barres
  theme_bw() 
p



