###########################""

# leaf trait
# manip Rhizotron
# post-doc JENNY

# 28/08/2015
# MZ
###############################"


# open file :
setwd("C:/Users/MZ/Documents/JENNY/DATA/rhizotron/RESULTS/leaf_trait/")
tab <- read.csv ("Leaf_trait.csv", sep=";", dec=".", header=TRUE, na.strings=NA)
summary(tab)

colnames(tab) <- c("id_loca", "id_div", "id_leg", "id_taxo", "division", "FW", "DW", "LA", "LDMC", "SLA", "LT")
# on crÈe une colonne id_RT
tab$id_loca <- as.character(tab$id_loca)
id_loca <- strmrt(x=tab$id_loca, mrt="_")
id_loca <- unlist(id_loca)
z <- seq(from=1, to=length(id_loca), by=2)
tab$id_RT <- id_loca[z]
z <- seq(from=2, to=length(id_loca), by=2)
tab$loca <- id_loca[z]
head(tab)

# check data for each species (outliers)
#*****************************************
# Bra_napus
#-------------
i = "Bra_napus"
tempo <- tab[tab$id_taxo==i,]
plot (tempo$DW ~tempo$FW)
reg_BN <- lm(tempo$DW ~tempo$FW)
summary(reg_BN)
tempo[tempo$FW> 1500 & tempo$DW<100,]
# id_loca id_div id_assem.v√©g√©tal   id_taxo division   FW   DW       LA    LDMC      SLA        LT
# 70 RT10_d0            Div1             Leg- Bra_napus               Leaf4 1873 22.6 4355.525 12.0662 192.7223 0.4300285
#

# on supprime l'outlier
#-----------------------
tab[tab$id_loca=="RT10_d0", c("FW", "DW", "LDMC", "LA", "SLA", "LT")]= NA
#


# check leaf traits per treatment
par (mfrow = c(3,2))
i = "Bra_napus"
v="LDMC"
tempo <- tab[tab$id_taxo==i,]
boxplot (tempo[,v]~tempo$id_leg, ylab=v, main=i)
boxplot (tempo[,v]~tempo$id_div, ylab=v, main=i)
# check the 2  points > 300 mg g-1
tab[tab$id_taxo==i & tab$LDMC>300,]
# id_loca id_div id_assem.v√©g√©tal   id_taxo division    FW    DW       LA     LDMC      SLA        LT
# 46 RT07_d+15            Div2             Leg- Bra_napus               Leaf4 596.9 189.5 1266.059 317.4736 6.681050 0.4714630
# 72 RT10_d-15            Div1             Leg- Bra_napus               Leaf3 691.4 209.2 1538.939 302.5745 7.356305 0.4492706
#
#
v="SLA"
boxplot (tempo[,v]~tempo$id_leg, ylab=v, main=i)
boxplot (tempo[,v]~tempo$id_div, ylab=v, main=i)
# check the 3  points > 300 mg g-1
tab[tab$id_taxo==i & tab$SLA>14,]
# id_loca id_div id_assem.v√©g√©tal   id_taxo division     FW    DW       LA     LDMC      SLA        LT
# 11     RT02_d0            Div1             Leg- Bra_napus               Leaf3 1588.2 262.4 3854.121 165.2185 14.68796 0.4120784
# 13   RT02_d-15            Div1             Leg- Bra_napus               Leaf3  198.6  34.9  594.770 175.7301 17.04212 0.3339106
# 24   RT04_d+15            Div3             Leg- Bra_napus               Leaf2  295.4  66.5  978.259 225.1185 14.71066 0.3019650
# 76   RT11_d+15            Div0             Leg+ Bra_napus               Leaf4  871.1 130.7 2153.454 150.0402 16.47631 0.4045129
#
#
v="LT"
boxplot (tempo[,v]~tempo$id_leg, ylab=v, main=i)
boxplot (tempo[,v]~tempo$id_div, ylab=v, main=i)


# Lol_peren
#-------------
i = "Lol_peren"
tempo <- tab[tab$id_taxo==i,]
par (mfrow = c(1,1))
plot (tempo$DW ~tempo$FW)
reg_LP <- lm(tempo$DW ~tempo$FW)
summary(reg_LP)
tempo[tempo$FW<50 & tempo$DW>20,]
tempo[tempo$FW>100,]
# check leaf traits per treatment
par (mfrow = c(3,2))
#
v="LDMC"
tempo <- tab[tab$id_taxo==i,]
boxplot (tempo[,v]~tempo$id_leg, ylab=v, main=i)
boxplot (tempo[,v]~tempo$id_div, ylab=v, main=i)
#
v="SLA"
boxplot (tempo[,v]~tempo$id_leg, ylab=v, main=i)
boxplot (tempo[,v]~tempo$id_div, ylab=v, main=i)
#
v="LT"
boxplot (tempo[,v]~tempo$id_leg, ylab=v, main=i)
boxplot (tempo[,v]~tempo$id_div, ylab=v, main=i)

# check the 3  points > 400 mg g-1
tab[tab$id_taxo==i & tab$LDMC> 400,]
# id_loca id_div id_assem.v√©g√©tal   id_taxo division   FW   DW      LA     LDMC      SLA        LT
# 36  RT05_d-12            Div0             Leg- Lol_peren               Leaf5 44.2 22.5 272.774 509.0498 12.12329 0.1620389
# 38   RT05_d-3            Div0             Leg- Lol_peren               Leaf4 43.4 19.5 346.175 449.3088 17.75256 0.1253701
# 69   RT10_d+9            Div1             Leg- Lol_peren               Leaf7 36.2 22.2 306.997 613.2597 13.82869 0.1179165
# 159  RT19_d+9            Div0             Leg+ Lol_peren               Leaf3 30.8 12.8 210.462 415.5844 16.44234 0.1463447
#

# Tri_alexa
#-------------
i = "Tri_alexa"
tempo <- tab[tab$id_taxo==i,]
plot (tempo$DW ~tempo$FW)
reg_TA <- lm(tempo$DW ~tempo$FW)
summary(reg_TA)

par (mfrow=c(1,3))
v="LDMC"
boxplot (tempo[,v]~tempo$id_div, ylab="LDMC (mg/g)", main=i)

v="SLA"
boxplot (tempo[,v]~tempo$id_div, ylab="SLA (mg/g)", main=i)

v="LT"
boxplot (tempo[,v]~tempo$id_div, ylab="Leaf Thickness (mm)", main=i)
tab[tab$id_taxo==i & tab$LT>0.25,]
# id_loca id_div id_leg   id_taxo division  FW  DW      LA     LDMC  SLA      LT        id_RT oca
# 163 RT19_d-18            Div0             Leg+ Tri_alexa        Leaf3 34.4 4.7 121.838 136.6279  25.92298 0.2823421  RT19 d-18


# on fait les moyenne de traits par esp et par rhizotron (=moyenne par populations)
#-----------------------------------------------------------------------------------
vb <- c("LDMC", "SLA", "LT")
mrt <- aggregate(tab[,vb], by=list(RT= tab$id_RT, id_taxo=tab$id_taxo, id_div=tab$id_div, id_leg=tab$id_leg), mean, na.rm=TRUE)
summary(mrt)
mrt <- mrt[order(mrt$RT, decreasing=FALSE),]
sdrt <- aggregate(tab[,vb], by=list(RT= tab$id_RT, id_taxo=tab$id_taxo, id_div=tab$id_div, id_leg=tab$id_leg), sd, na.rm=TRUE)
sdrt <- sdrt[order(sdrt$RT, decreasing=FALSE),]
par (mfrow = c(1,3))
plot(sdrt$LDMC~mrt$LDMC)
sdrt[sdrt$LDMC>60,]
# div  leg   id_taxo      LDMC      SLA         LT
# 9  Div0 Leg- Lol_peren  93.60817 6.837579 0.07016581
# 10 Div1 Leg- Lol_peren 108.34371 7.296353 0.05691758
# 13 Div0 Leg+ Lol_peren  81.93123 4.103495 0.04872875

plot(sdrt$SLA~mrt$SLA)
plot(sdrt$LT~mrt$LT)


# analyses des traits foliaires
#*******************************

# 1.caractÈrisation du traitement trËfle
#-------------------------------------
ta <- tab[tab$id_taxo=="Tri_alexa",]
TA <- mrt[mrt$id_taxo=="Tri_alexa",]
#
# nombre de trËfles par rhizotrons
ta$nb = 1
nbTA_rt <- aggregate(ta$nb, by=list(id_RT=ta$id_RT, id_div=ta$id_div), sum, na.rm=TRUE)
# nombre de TA moyen par RT
mTA <- aggregate(nbTA_rt$x, by=list(id_div= nbTA_rt$id_div), mean, na.rm =TRUE)


# on calcule les moyennes de traits pour chaque traitement de diversitÈ
meanTA <- aggregate(TA[,vb], by=list(id_taxo=TA$id_taxo, id_div=TA$id_div), mean, na.rm=TRUE)
sdTA <-  aggregate(TA[,vb], by=list(id_taxo=TA$id_taxo, id_div=TA$id_div), sd, na.rm=TRUE)
meanTA <- meanTA[order(meanTA$id_div, decreasing=FALSE),]
sdTA <- sdTA[order(sdTA$id_div, decreasing=FALSE),]
library(plotrix)
seTA <-  aggregate(TA[,vb], by=list(id_taxo=TA$id_taxo, div=TA$id_div), std.error, na.rm=TRUE)


# stat sur moyennes... ‡ partir des donnÈes sur individus
mod1 <- lm(LDMC~id_div, data=ta, na.action=na.omit)
anova(mod1)
library(multcomp)
a<- glht(mod1, lincft=mcp (id_div = "Tukey"))
summary(a)

mod2 <- lm(SLA~id_div, data=ta, na.action=na.omit)
anova(mod2)

mod3 <- lm(LT~id_div, data=ta, na.action=na.omit)
anova(mod3)

## relations competition colza et lolium selon les traitements Leg et diversit√©
#*******************************************************************************

# relation LDMC BN- LDMC LP
#----------------------------
par(mfrow=c(1,2))

i = c("Bra_napus", "Lol_peren")
xlim=c(min(mrt$LDMC), max(mrt$LDMC))
ylim= c(min(mrt$LDMC), max(mrt$LDMC))
plot(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-"]
     ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"], 
     xlab="LDMC Bra_napus (mg g-1)", ylab="LDMC Lol_peren (mg g-1)", 
     xlim=xlim, ylim=ylim,pch=16)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"], 
       xlab="", ylab="", 
       xlim=xlim, ylim=ylim, pch=1)
lm1 <- lm(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-"]
          ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"])
summary(lm1)
lm2 <- lm(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+"]
          ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"])
summary(lm2)
abline(lm2, lty="dotted")
text(310,350,"R≤=0.36*")
legend(270, 210, legend=c("Leg-", "Leg+"), pch=c(16, 1), lty=c("solid", "dotted"))

# relation SLA BN-  LP
#----------------------------
xlim=c(min(mrt$SLA[mrt$id_taxo!="Tri_alexa"]), max(mrt$SLA[mrt$id_taxo!="Tri_alexa"]))
ylim= c(min(mrt$SLA[mrt$id_taxo!="Tri_alexa"]), max(mrt$SLA[mrt$id_taxo!="Tri_alexa"]))
# par(mfrow=c(1,1))
plot(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-"]
     ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"], 
     xlab="SLA Bra_napus (mg g-1)", ylab="SLA Lol_peren (mg g-1)", 
     xlim=xlim, ylim=ylim,pch=16)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"], 
       xlab="", ylab="", 
       xlim=xlim, ylim=ylim, pch=1)
lm1 <- lm(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-"]
          ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"])
summary(lm1)
lm2 <- lm(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+"]
          ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"])
summary(lm2)
abline(lm2, lty="dotted")
text(17,25,"R≤=0.53**")
legend(17,11, legend=c("Leg-", "Leg+"), pch=c(16, 1), lty=c("solid", "dotted"))



# pour visualiser les traitement div:
#
xlim=c(min(mrt$LDMC), max(mrt$LDMC))
ylim= c(min(mrt$LDMC), max(mrt$LDMC))
plot(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div0"]
     ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div0"], 
     xlab="LDMC Bra_napus (mg g-1)", ylab="LDMC Lol_peren (mg g-1)", xlim=xlim, ylim=ylim,
     pch=15)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div1"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=16)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div2"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=17)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div3"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=23, bg="black")
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div0"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div0"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=0)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div1"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=1)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div2"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=2)
points(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div3"]
       ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=5)
abline(1,1, col="grey")
lm2 <- lm(mrt$LDMC[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+"]
          ~mrt$LDMC[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"])
summary(lm2)
abline(lm2, lty="dotted")
text(310,350,"R≤=0.36*")

lm0<- lm(mrt$LDMC[mrt$id_taxo=="Lol_peren"]~mrt$LDMC[mrt$id_taxo=="Bra_napus"])
summary(lm0)
abline(lm0, col="red")
text(320,340,"R≤=0.26*", col="red")


# SLA
xlim=c(min(mrt$SLA[mrt$id_taxo!="Tri_alexa"]), max(mrt$SLA[mrt$id_taxo!="Tri_alexa"]))
ylim= c(min(mrt$SLA[mrt$id_taxo!="Tri_alexa"]), max(mrt$SLA[mrt$id_taxo!="Tri_alexa"]))
plot(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div0"]
     ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div0"], 
     xlab="SLA Bra_napus (mg g-1)", ylab="SLA Lol_peren (mg g-1)", xlim=xlim, ylim=ylim,
     pch=15)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div1"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=16)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div2"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=17)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg-" & mrt$id_div=="Div3"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg-"& mrt$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=23, bg="black")
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div0"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div0"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=0)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div1"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=1)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div2"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=2)
points(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+" & mrt$id_div=="Div3"]
       ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"& mrt$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=5)
abline(1,1, col="grey")
lm2 <- lm(mrt$SLA[mrt$id_taxo=="Lol_peren"& mrt$id_leg=="Leg+"]
          ~mrt$SLA[mrt$id_taxo=="Bra_napus"&mrt$id_leg=="Leg+"])
summary(lm2)
abline(lm2, lty="dotted")
text(17,25,"R≤=0.53**")
legend(17, 17, legend=c("Div0", "Div1", "Div2", "Div3", "Leg-", "Leg+"), 
       pch=c(15, 16,17, 23, 16, 1), lty=c("blank","blank","blank","blank","solid", "dotted"), pt.bg="black")


lm0<- lm(mrt$SLA[mrt$id_taxo=="Lol_peren"]~mrt$SLA[mrt$id_taxo=="Bra_napus"])
summary(lm0)

# mÍme chose mais ‡ partir des moyenne par tmt
vb <- c("LDMC", "SLA", "LT")
m <- aggregate(tab[,vb], by=list(id_div=tab$id_div, id_leg=tab$id_leg, id_taxo=tab$id_taxo), mean, na.rm=TRUE)
summary(m)
# graphiques 
# LDMC BN-LP
xlim=c(min(m$LDMC), max(m$LDMC))
ylim= c(min(m$LDMC), max(m$LDMC))
plot(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div0"]
     ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div0"], 
     xlab="LDMC Bra_napus (mg g-1)", ylab="LDMC Lol_peren (mg g-1)", xlim=xlim, ylim=ylim,
     pch=15)
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div1"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=16)
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div2"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=17)
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div3"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=23, bg="black")
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div0"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div0"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=0)
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div1"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=1)
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div2"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=2)
points(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div3"]
       ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=5)
abline(1,1, col="grey")
lm2 <- lm(m$LDMC[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+"]
          ~m$LDMC[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"])
summary(lm2)

lm0<- lm(m$LDMC[m$id_taxo=="Lol_peren"]~m$LDMC[m$id_taxo=="Bra_napus"])
summary(lm0)

# SLA
xlim=c(min(m$SLA[m$id_taxo!="Tri_alexa"]), max(m$SLA[m$id_taxo!="Tri_alexa"]))
ylim= c(min(m$SLA[m$id_taxo!="Tri_alexa"]), max(m$SLA[m$id_taxo!="Tri_alexa"]))
plot(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div0"]
     ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div0"], 
     xlab="SLA Bra_napus (mg g-1)", ylab="SLA Lol_peren (mg g-1)", xlim=xlim, ylim=ylim,
     pch=15)
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div1"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=16)
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div2"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=17)
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-" & m$id_div=="Div3"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"& m$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=23, bg="black")
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div0"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div0"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=0)
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div1"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div1"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=1)
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div2"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div2"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=2)
points(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+" & m$id_div=="Div3"]
       ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"& m$id_div=="Div3"], 
       xlab="", ylab="", xlim=xlim, ylim=ylim,
       pch=5)
abline(1,1, col="grey")
lm1 <- lm(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg-"]
          ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg-"])
summary(lm1)

lm2 <- lm(m$SLA[m$id_taxo=="Lol_peren"& m$id_leg=="Leg+"]
          ~m$SLA[m$id_taxo=="Bra_napus"&m$id_leg=="Leg+"])
summary(lm2)
abline(lm2, lty="dotted")
text(14,21,"R≤=0.85+")
legend(16, 15, legend=c("Div0", "Div1", "Div2", "Div3", "Leg-", "Leg+"), 
       pch=c(15, 16,17, 23, 16, 1), lty=c("blank","blank","blank","blank","solid", "dotted"), pt.bg="black")


lm0<- lm(m$SLA[m$id_taxo=="Lol_peren"]~m$SLA[m$id_taxo=="Bra_napus"])
summary(lm0)
abline(lm0, col="red")
text(14,19,"R≤=0.67**", col="red")
