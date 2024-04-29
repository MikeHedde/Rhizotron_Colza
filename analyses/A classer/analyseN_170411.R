# 160601 : analyses colza
# 160606 : analyses LP et TA 
# 160609 : modi analyses stat




# chermin d'accés pour retrouver les données
setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/Données/Plantes/analyses_N") # ordi Marine


#### Librairies ####
require(ggplot2)
require(multcomp)
require(car)
require(nlme)
require(plotrix)

require(conover.test) # pour comparaison multiple après KW
require(PMCMR)


##### 1.lire la table :
n <- read.csv("N_160920.csv", header=TRUE, sep=";", dec=".", na.strings=NA)
summary(n)

# on crée une colonne serie
n$serie<- NULL
n$serie[n$RT<21]="s1"
n$serie[n$RT>20]="s2"
summary(n)

# on recrée une colonne organes pour supprimer les leaf1, leaf2...
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
# on remplace par la moyenne des points pour le même traitement : Leg-Div0
mean(bn$N[n$leg=="Leg-" & n$div=="Div0"], na.rm=T)
n$N[n$RT==23 & n$distance==-15] = mean(bn$N[n$leg=="Leg-" & n$div=="Div0"], na.rm=T)
bn <- n[n$id_taxo=="Bra_napus",]

# stat
data= bn
mod <- lm(asin(sqrt(N/100))~organes, data=data, na.action=na.omit)
shapiro.test(residuals(mod))
fligner.test(asin(sqrt(N/100))~organes,data=data, na.action=na.omit)
# tests non paramétriques
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


## représentations graphiques

# graphique : effet de la leg
p <- ggplot(bn, aes(leg, N)) +
  ylab("Concentration en azote (%)")+
  xlab("Plante de service") +
  geom_boxplot()+
  facet_grid(. ~leg)
p  

#graphique :  effet de la diversité dans les organes de colza
p <- ggplot(bn[bn$organes=="Fine roots",], aes(div, N)) +
  ylab("Concentration en azote dans les racines fines (%)")+
  xlab("Diversité fonctionnelle des vers de terre") +
  geom_boxplot()
p  

#graphique : effet de tous les traitements
p <- ggplot(bn[bn$organes=="Fine roots",], aes(leg, N)) +
  ylab("Concentration en azote dans les racines fines (%)")+
  xlab("Plante de service") +
  geom_boxplot()+
  facet_grid( . ~ div)
p  

# vérif des outliers
bn[bn$organes=="Fine roots" & bn$div=="Div0" & bn$Ncontent > 2.25,] # RT23
bn[bn$organes=="Fine roots" & bn$div=="Div2" & bn$Ncontent > 2.75,] # RT42
bn[bn$organes=="Fine roots" & bn$div=="Div3" & bn$Ncontent > 2.25,] # RT

#########################################################################

# faire une boucle pour ajouter DW dans le jeu de données




#########################################################################

##### bilan d'azote pour le colza
n$qN = n$Ncontent * n$DW / 100
bn <- n[n$id_taxo=="Bra_napus",]
vb <- c("qN", "DW")
bilan <- aggregate(bn[,vb], by=list(leg=bn$leg, div=bn$div, RT=bn$RT, id_taxo=bn$id_taxo, serie=bn$serie), sum, na.rm=FALSE)

# pour le moment, on élimine les colza incomplets
bil<- bilan[bilan$RT < 20|bilan$RT!=33|bilan$RT!=38|bilan$RT!=40,]

hist(bil$qN)

## effet des traitements sur le bilan d'azote 
p <- ggplot(bil, aes(leg, qN)) +
  ylab("Quantité totale d'azote exportée par le colza central (mg)")+
  xlab("Plante de service") +
  geom_boxplot()
p  



##################################################################################

 
## calucls à faire :
# bilan N / plante (mg) = somme de toutes les qté N ds organes x masse de chq organes
# efficience (NUE) = (qté N plante entière (mg)/ nb jours de croissance)/ MS racines fines (mg)
# efficience C = (Ctot (mg)/nb jours)/surface foliaire totale
# --> récupérer date de germination pour calucler age des plantes

# rq : efficience renseigne sur la capacité de la plante à pomper N par unité de surface  ou masse racinaire produite
# efficience C renseigne sur la capacité de la platne à produire la biomasse avec une qté N et surface foliaire données

# les plantes avec + de biomasse produites pour une même qté N entrée => plantes les + efficientes (NUE forte)

# si les traitements modifient lg racinaires, morpho racinaires, distrib des raicnes dans le sols 
# + preturbation de la qté N dispo dans le sol par les vdt
# --> quelles conséquences obsevées sur la NUE ?