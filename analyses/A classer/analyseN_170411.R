# 160601 : analyses colza
# 160606 : analyses LP et TA 
# 160609 : modi analyses stat




# chermin d'acc�s pour retrouver les donn�es
setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/Donn�es/Plantes/analyses_N") # ordi Marine


#### Librairies ####
require(ggplot2)
require(multcomp)
require(car)
require(nlme)
require(plotrix)

require(conover.test) # pour comparaison multiple apr�s KW
require(PMCMR)


##### 1.lire la table :
n <- read.csv("N_160920.csv", header=TRUE, sep=";", dec=".", na.strings=NA)
summary(n)

# on cr�e une colonne serie
n$serie<- NULL
n$serie[n$RT<21]="s1"
n$serie[n$RT>20]="s2"
summary(n)

# on recr�e une colonne organes pour supprimer les leaf1, leaf2...
unique(n$org)
n$organes <- NULL
n$organes[n$org=="Leaf3"|n$org=="Leaf2"|n$org=="Leaf4"|n$org=="Leaf5"|n$org=="Leaf7"|n$org=="Leaf1"|n$org=="Limb"]="Limb"
n$organes[n$org=="Fine roots"]="Fine roots"
n$organes[n$org=="Taproot"]="Taproot"
n$organes[n$org=="Other Leaves"|n$org=="Other leaves"|n$org=="Aerial"]="Leaves"
n$organes[n$org=="Stem"]="Stem"
n$organes <- as.factor(n$organes)

## effet organes
#-------------------------
bn <- n[n$id_taxo=="Bra_napus",]
lp <- n[n$id_taxo=="Lol_peren",]
ta <- n[n$id_taxo=="Tri_alexa",]

# graphique
p <- ggplot(bn, aes(organes, N)) +
  ylab("Concentration en azote (%)")+
  xlab("Organes") +
  geom_boxplot()
p  
# verfi des outliers:
bn[bn$N>10,]
# RT 23 BN-15 Fine roots
# on remplace par la moyenne des points pour le m�me traitement : Leg-Div0
mean(bn$N[n$leg=="Leg-" & n$div=="Div0"], na.rm=T)
n$N[n$RT==23 & n$distance==-15] = mean(bn$N[n$leg=="Leg-" & n$div=="Div0"], na.rm=T)
bn <- n[n$id_taxo=="Bra_napus",]

# stat
data= bn
mod <- lm(asin(sqrt(N/100))~organes, data=data, na.action=na.omit)
shapiro.test(residuals(mod))
fligner.test(asin(sqrt(N/100))~organes,data=data, na.action=na.omit)
# tests non param�triques
kruskal.test(asin(sqrt(N/100))~organes,data=data, na.action=na.omit )
# comparaison multiple
posthoc.kruskal.nemenyi.test(asin(sqrt(N/100))~organes, data=data, dist="Tukey")
conover.test(asin(sqrt(data$N/100)), data$organes)

## stat effet traitement
tempo <- lp#[bn$organe=="Leaves",]
mod <- lm(asin(sqrt(N/100))~serie+leg*div, data=tempo, na.action=na.omit)
summary(mod)
shapiro.test(residuals(mod))
fligner.test(asin(sqrt(N/100))~leg, data=tempo, na.action=na.omit)
fligner.test(asin(sqrt(N/100))~div, data=tempo, na.action=na.omit)
Anova(mod, type="III")

kruskal.test(asin(sqrt(N/100))~leg,data=data, na.action=na.omit )
kruskal.test(asin(sqrt(N/100))~div,data=data, na.action=na.omit )


## repr�sentations graphiques

# graphique : effet de la leg
p <- ggplot(bn, aes(leg, N)) +
  ylab("Concentration en azote (%)")+
  xlab("Plante de service") +
  geom_boxplot()+
  facet_grid(. ~leg)
p  

#graphique :  effet de la diversit� dans les organes de colza
p <- ggplot(bn[bn$organes=="Fine roots",], aes(div, N)) +
  ylab("Concentration en azote dans les racines fines (%)")+
  xlab("Diversit� fonctionnelle des vers de terre") +
  geom_boxplot()
p  

#graphique : effet de tous les traitements
p <- ggplot(bn[bn$organes=="Fine roots",], aes(leg, N)) +
  ylab("Concentration en azote dans les racines fines (%)")+
  xlab("Plante de service") +
  geom_boxplot()+
  facet_grid( . ~ div)
p  

# v�rif des outliers
bn[bn$organes=="Fine roots" & bn$div=="Div0" & bn$Ncontent > 2.25,] # RT23
bn[bn$organes=="Fine roots" & bn$div=="Div2" & bn$Ncontent > 2.75,] # RT42
bn[bn$organes=="Fine roots" & bn$div=="Div3" & bn$Ncontent > 2.25,] # RT

#########################################################################

# faire une boucle pour ajouter DW dans le jeu de donn�es




#########################################################################

##### bilan d'azote pour le colza
n$qN = n$Ncontent * n$DW / 100
bn <- n[n$id_taxo=="Bra_napus",]
vb <- c("qN", "DW")
bilan <- aggregate(bn[,vb], by=list(leg=bn$leg, div=bn$div, RT=bn$RT, id_taxo=bn$id_taxo, serie=bn$serie), sum, na.rm=FALSE)

# pour le moment, on �limine les colza incomplets
bil<- bilan[bilan$RT < 20|bilan$RT!=33|bilan$RT!=38|bilan$RT!=40,]

hist(bil$qN)

## effet des traitements sur le bilan d'azote 
p <- ggplot(bil, aes(leg, qN)) +
  ylab("Quantit� totale d'azote export�e par le colza central (mg)")+
  xlab("Plante de service") +
  geom_boxplot()
p  



##################################################################################

 
## calucls � faire :
# bilan N / plante (mg) = somme de toutes les qt� N ds organes x masse de chq organes
# efficience (NUE) = (qt� N plante enti�re (mg)/ nb jours de croissance)/ MS racines fines (mg)
# efficience C = (Ctot (mg)/nb jours)/surface foliaire totale
# --> r�cup�rer date de germination pour calucler age des plantes

# rq : efficience renseigne sur la capacit� de la plante � pomper N par unit� de surface  ou masse racinaire produite
# efficience C renseigne sur la capacit� de la platne � produire la biomasse avec une qt� N et surface foliaire donn�es

# les plantes avec + de biomasse produites pour une m�me qt� N entr�e => plantes les + efficientes (NUE forte)

# si les traitements modifient lg racinaires, morpho racinaires, distrib des raicnes dans le sols 
# + preturbation de la qt� N dispo dans le sol par les vdt
# --> quelles cons�quences obsev�es sur la NUE ?