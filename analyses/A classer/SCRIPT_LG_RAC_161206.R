####
 # modil le 25 aout 2016 pour elimier RT douteux
# modif 1/10/2016 : mise à jour script ) partir 
# modif 8/11/2016 : trace la relation mass root- burrow length
# 29/11/2016 : modif de la colonne biomDelta 

#####################

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/Données/Plantes/longueur_racines")
# setwd("G:/PRO_jenny/JENNY/RHIZOTRON/Données/Plantes/longueur_racines/")



#### Librairies ####
require(ggplot2)
require(vegan)
require(ade4)
require(reshape2)
require(multcomp)
require(nlme)
require(car)
require(plotrix)
require(PMCMR)

##Données
##### 
rac <- read.csv("RAC_LgTOT.csv", header = TRUE, sep = ";", dec=".", na.strings=NA)
summary(rac)

rac$serie <- NULL
rac$serie[rac$RT<21]="s1"
rac$serie[rac$RT>20]="s2"
rac$serie <- as.factor(rac$serie)

# on elim les RT douteux
rac[rac$RT==12 | rac$RT==8 | rac$RT==18 | rac$RT==40,c(4:18)]=NA
# on crée une colonne traitment
rac$tmt <- paste(rac$leg, rac$div, sep="")
rac$tmt <- as.factor(rac$tmt)

# on ajoute les données vdt
vdt <- read.csv("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/Données/vdt/masse/DATA_RT_VDT_biomasse_COMT.csv", sep=";", dec=".", header=T, na.strings=NA)
summary(vdt)
result <- NULL
for(i in unique (rac$RT))
{
  tempo <- rac[rac$RT==i,]
  if(tempo$div=="Div0")
  {
    tempo$biomDelta = 0
  } else {
    tempo2 <- vdt[vdt$RT==i,]
    tempo$biomDelta <- tempo2$deltaBiom
  }
  result <- rbind(result, tempo)
}
summary(result)
rac <- result
# on elim les RT douteux
rac[rac$RT==12 | rac$RT==8 | rac$RT==18 | rac$RT==40,c(4:18)]=NA


#####
# on calcule les moyennes de longueurs
moy <- aggregate(rac$long, by=list(leg=rac$leg, div=rac$div), mean, na.rm=T)
colnames(moy) <- c("leg", "div", "Root_length")
se <- aggregate(rac$long, by=list(leg=rac$leg, div=rac$div), std.error, na.rm=T)
moy$se <- se$x
write.table(moy, "DATA_RT_LONG_RAC_MOY.csv", sep=";", dec=".", col.names=T, row.names = F)


#############################################################"
##Effet des traitements sur les longueur de racines
#-----------------------------------------------------
p <- ggplot(rac, aes(leg, long)) +
  ylab("Longueur de racines (m)")+
  xlab("Présence de légumineuse") +
  geom_boxplot() +
  facet_grid( . ~ div)
p  

rac[rac$div=="Div3" & rac$long<40,]


# test stat
# mod <- lme(long~leg*div, random=~1|serie, data=rac, na.action=na.omit)
mod <- lm(long~serie+leg*div, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
fligner.test(long~div, data=rac, na.action=na.omit)
fligner.test(long~leg, data=rac, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
# summary(aov(mod)) # Effet "diversité"
# comparaison multiple
# comp <- glht(mod, linfct=mcp(div = "Tukey"))
# summary(comp)
#tests non parmétriques
kruskal.test(long~serie, data=rac, na.action=na.omit)
kruskal.test(long~leg, data=rac, na.action=na.omit)
kruskal.test(long~div, data=rac, na.action=na.omit)

# graphique
moy <- aggregate(rac$long, by=list(leg=rac$leg, div=rac$div), mean, na.rm=T)
colnames(moy) <- c("leg", "div", "Root_length")
se <- aggregate(rac$long, by=list(leg=rac$leg, div=rac$div), std.error, na.rm=T)
moy$se <- se$x

png("Total root length.png")
p <- ggplot(moy, aes(x=leg, y=Root_length, fill=div)) + 
  geom_bar(position=position_dodge(), stat="identity", 
           colour="black", # use black outlines, 
           size=0.3)+   # for thinner lines
  geom_errorbar(aes(ymin=Root_length-se, ymax=Root_length+se),
                size=0.3, # for thinner lines
                width=0.2, # Width of the error bars
                position=position_dodge(.9))+
  ylab("Total root lenght (m)")+
  xlab("")+
  scale_y_continuous(breaks=0:60*5) +
  scale_fill_manual(values=c("white", "azure3", "azure4", "black")) + # pour changer les couleurs des barres
  theme_bw() 
p
dev.off()

png("Total root length_color.png")
p <- ggplot(moy, aes(x=leg, y=Root_length, fill=div)) + 
  geom_bar(position=position_dodge(), stat="identity", 
           colour="black", # use black outlines, 
           size=0.3)+   # for thinner lines
  geom_errorbar(aes(ymin=Root_length-se, ymax=Root_length+se),
                size=0.3, # for thinner lines
                width=0.2, # Width of the error bars
                position=position_dodge(.9))+
  ylab("Total root lenght (m)")+
  xlab("")+
  scale_y_continuous(breaks=0:60*5) +
  scale_fill_manual(values=c("black", "orange", "purple", "green")) + # pour changer les couleurs des barres
  theme_bw() 
p
dev.off()

# effet des traitements sur la connectivité
p <- ggplot(rac, aes(leg, connec)) +
  ylab("Connectivité)")+
  xlab("Présence de légumineuse") +
  geom_boxplot() +
  facet_grid( . ~ div)
p  

# summary(aov(connec~leg*div, data = mydata)) #

##Effet des traitements sur les masses racinaires totales
#--------------------------------------------------------
rac$biomR = rac$DW/1000 # en g
summary(rac)

p <- ggplot(rac, aes(leg, biomR)) +
  ylab("Biomasse racinaire (g)")+
  xlab("Présence de légumineuse") +
  geom_boxplot() +
  facet_grid( . ~ div)
p  

rac[rac$div=="Div3" & rac$leg=="Leg-" & rac$biomRac<6,]

#summary(aov(biomMoy~leg*div, data = mydata)) # Effet "diversité"
# mod <- lme(biomR~nb_pl+leg*div, random=~1|serie, data=rac, na.action=na.omit)
mod <- lm(biomR~serie+leg*div, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
bartlett.test(biomR~div, data=rac, na.action=na.omit)
bartlett.test(biomR~leg, data=rac, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
# summary(aov(mod)) # Effet "diversité"
# comparaison multiple
comp <- glht(mod, linfct=mcp(div = "Tukey"))
summary(comp)
#
mod2 <- lme(biomR~tmt, random=~1|serie, data=rac, na.action=na.omit)
Anova(mod2)
comp <- glht(mod2, linfct=mcp(tmt = "Tukey"))
summary(comp)

##################################################################################################

# effet des traitements sur le specific root length
#-------------------------------------------------------

# on calcule le SRL = root lenth/ Root mass (m/g)
rac$SRL = rac$long/rac$biomR

# effet des traitments sur le SRL
mod <- lm(SRL~serie+leg*div, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
fligner.test(SRL~div, data=rac, na.action=na.omit)
fligner.test(SRL~leg, data=rac, na.action=na.omit)
fligner.test(SRL~serie, data=rac, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
# comparaison multiple
comp <- glht(mod, linfct=mcp(div = "Tukey"))
summary(comp)

mean(rac$SRL[rac$div=="Div0"],na.rm=T)
sd(rac$SRL[rac$div=="Div0"],na.rm=T)
mean(rac$SRL[rac$div!="Div0"],na.rm=T)
sd(rac$SRL[rac$div!="Div0"],na.rm=T)

mod2 <- lm(SRL~serie+tmt, data=rac, na.action=na.omit)
shapiro.test(residuals(mod2))
bartlett.test(SRL~tmt, data=rac, na.action=na.omit)
summary(mod2)
Anova(mod2, type="III")
comp <- glht(mod2, linfct=mcp(tmt = "Tukey"))
summary(comp)


# graphique
p <- ggplot(rac, aes(leg, SRL)) +
  ylab("SRL(m/g MS)")+
  xlab("") +
  geom_boxplot() +
  facet_grid( . ~ div)
p  

moy <- aggregate(rac$SRL, by=list(div=rac$div, leg=rac$leg), mean, na.rm=T)
se <- aggregate(rac$SRL, by=list(div=rac$div, leg=rac$leg), std.error, na.rm=T)
moy$se <- se$x

write.table(moy,"DATA_RT_SRL_MOY.csv", sep=";", dec=".", col.names=T, row.names=F)


png("Specific Root lenght.png")
p <- ggplot(moy, aes(x=leg, y=Root_length, fill=div)) + 
  geom_bar(position=position_dodge(), stat="identity", 
           colour="black", # use black outlines, 
           size=0.3)+   # for thinner lines
  geom_errorbar(aes(ymin=Root_length-se, ymax=Root_length+se),
                size=0.3, # for thinner lines
                width=0.2, # Width of the error bars
                position=position_dodge(.9))+
  ylab("Total root lenght (m)")+
  xlab("")+
  scale_y_continuous(breaks=0:60*5) +
  scale_fill_manual(values=c("white", "azure3", "azure4", "black")) + # pour changer les couleurs des barres
  theme_bw() 
p
dev.off()


# on calcule le specific burrow lenght
#------------------------------------
rac$SBL <- rac$gal/rac$biomFin
moy <- aggregate(rac$SBL, by=list(div=rac$div, leg=rac$leg), mean, na.rm=T)
se <- aggregate(rac$SBL, by=list(div=rac$div, leg=rac$leg), std.error, na.rm=T)
moy$se <- se$x

write.table(moy,"DATA_RT_SBL_MOY.csv", sep=";", dec=".", col.names=T, row.names=F)

# effet des traitements sur SBL
mod <- lm(SBL~serie+leg*div, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
bartlett.test(SBL~div, data=rac, na.action=na.omit)
bartlett.test(SBL~leg, data=rac, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
# comparaison multiple
comp <- glht(mod, linfct=mcp(div = "Tukey"))
summary(comp)

mod <- lm(SBL~tmt, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
bartlett.test(SBL~tmt, data=rac, na.action=na.omit)
summary(mod)
Anova(mod, type="III")
comp <- glht(mod, linfct=mcp(tmt = "Tukey"))
summary(comp)

mod1 <- lm(SBL~serie+div, data=rac[rac$leg=="Leg-",], na.action=na.omit)
summary(mod1)
Anova(mod1, type="III")
comp <- glht(mod1, linfct=mcp(div = "Tukey"))
summary(comp)

mod2 <- lm(SBL~serie+div, data=rac[rac$leg=="Leg+",], na.action=na.omit)
summary(mod2)
Anova(mod2, type="III")
comp <- glht(mod, linfct=mcp(div = "Tukey"))
summary(comp)

# graphique
p <- ggplot(rac, aes(div, SBL)) +
  ylab("Specific Burrow Length (m/g MF)")+
  xlab("") +
  geom_boxplot() +
  facet_grid( . ~ leg)
p  

# ensuite on met un 0 dans la colonne SBL pour les Div0
rac$SBL[rac$div=="Div0"]=0


## on sauve la table avec les traits rac :
write.table(rac, "DATA_RT_RAC.csv",sep=";", dec=".", col.names=T, row.names = F)

###########################################################################################
###########################################################################################
rac <- read.csv("DATA_RT_RAC.csv",sep=";", dec=".", header=T, na.strings=NA)

##  Relation entre la longueur de galeries et la longueur des racines
#------------------------------------------------------------------
## tous traitements confondus
reg <- lm(long~gal, data = rac)
coeff <- coefficients(reg)
summary(reg)  # relation légèrement significative
# Equation de la droite de regression :  
eq = paste0("y = ", round(coeff[2],1), "*x + ", round(coeff[1],1))

p <- ggplot(rac, aes(gal, long, col = div, pch = leg)) +
  xlab("Longueur de galeries (m)")+
  ylab("Longueur de racines (m)") +
  geom_point() +
  geom_abline(intercept = 51.3, slope = -0.5, linetype="dashed", size=1.5)+
  ggtitle(eq)
p

par(mfrow=c(1,1))
xlim=c(0,45)
ylim=c(25,65)
plot(rac$long[rac$div=="Div0" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
      xlim=xlim, ylim=ylim, xlab="Total burrow length (m)", ylab="Total root lenght (m)", las=1)
points(rac$long[rac$div=="Div0" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(long~gal, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
abline (reg)
text(40,40, "R²=0.38***")
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))

# relation sans les point Div0
reg <- lm(long~gal, data=rac[rac$div!="Div0",], na.action=na.omit)
shapiro.test(residuals(mod))
plot(rac$long[rac$div!="Div0"]~rac$gal[rac$div!="Div0"])
summary(reg)
abline (reg)
text(40,40, "R²=0.20*")

# on trace la relation avec les moyennes seulement
vb <- c("long", "gal", "biomR")
moy <- aggregate(rac[,vb], by=list(leg=rac$leg, div=rac$div), mean, na.rm=T)
se <-  aggregate(rac[,vb], by=list(leg=rac$leg, div=rac$div), std.error, na.rm=T)

plot(long~gal, data=moy)
reg <- lm(long~gal, data=rac, na.action=na.omit)
summary(reg)
abline(reg)
text(20, 40, "R²=0.38***")


## relation entre la logueur des racines et la connectance des galeries
#-----------------------------------------------------------
plot(rac$long ~rac$connec)
par(mfrow=c(1,1))
xlim=c(0,0.5)
ylim=c(25,65)
plot(rac$long[rac$div=="Div0" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Burrow connectance", ylab="Total root lenght (m)", las=1)
points(rac$long[rac$div=="Div0" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(long~connec, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))

# on trace la relation avec les moyennes seulement
vb <- c("long", "gal", "biomR","connec")
moy <- aggregate(rac[,vb], by=list(leg=rac$leg, div=rac$div), mean, na.rm=T)
se <-  aggregate(rac[,vb], by=list(leg=rac$leg, div=rac$div), std.error, na.rm=T)
plot(long~connec, data=moy)
reg <- lm(long~connec, data=rac, na.action=na.omit)
summary(reg)


## relation longueur racines - masse vdt
#-------------------------------------------
plot(rac$long ~rac$biomFin)
par(mfrow=c(1,1))

xlim=c(0,7)
ylim=c(25,65)
plot(rac$long[rac$div=="Div0" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Earthworm final biomass (g)", ylab="Total root lenght (m)", las=1)
points(rac$long[rac$div=="Div0" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(long~biomFin, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
abline(reg)
text(6,30,"R²=0.14*")
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))


plot(long~biomDelta, data=rac, xlab="Delta vdt", ylab="longeur des racines (m)")

xlim=c(-100,0)
ylim=c(25,65)
plot(rac$long[rac$div=="Div0" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Delta masse des vers de terre", ylab="Longueur totale des racines (m)", las=1)
points(rac$long[rac$div=="Div0" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(long~biomDelta, data=rac, na.action=na.omit)
summary(reg)
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))
dev.off()

# relation entre la longueur des racine et SBL
#------------------------------------------------
plot(long~SBL, data=rac, xlab="Specific burrow length (m/g MF)", ylab="Longueur racines (m)")
xlim=c(0,25)
ylim=c(25,65)
plot(rac$long[rac$div=="Div0" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Specific burrow length (m/g FM)", ylab="Total root length (m)", las=1)
points(rac$long[rac$div=="Div0" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$long[rac$div=="Div1" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$long[rac$div=="Div2" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$long[rac$div=="Div3" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(long~SBL, data=rac, na.action=na.omit)
summary(reg)
abline(reg)
text(20,36,"R²=0.25***")

############################################################################################################################
############################################################################################################################

# relation entre la masse des racines et la longueur des galeries
#------------------------------------------------------------------
plot(rac$biomR ~rac$gal)
par(mfrow=c(1,1))
xlim=c(0,45)
ylim=c(4,9.5)
plot(rac$biomR[rac$div=="Div0" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Total burrow length (m)", ylab="Total root mass (g DM)", las=1)
points(rac$biomR[rac$div=="Div0" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(biomR~gal, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))
reg0 <- lm(biomR~gal, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg0)
reg1 <- lm(biomR~gal, data=rac[rac$div=="Div1",], na.action=na.omit)
summary(reg1)
reg2 <- lm(biomR~gal, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
reg3 <- lm(biomR~gal, data=rac[rac$div=="Div3",], na.action=na.omit)
summary(reg3)


# relation entre la masse des racines et la connectance des galeries
#------------------------------------------------------------------
plot(rac$biomR ~rac$connec)
par(mfrow=c(1,1))
xlim=c(0,0.5)
ylim=c(4,9.5)
plot(rac$biomR[rac$div=="Div0" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Burrow connectance", ylab="Total root mass (g DM)", las=1)
points(rac$biomR[rac$div=="Div0" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(biomR~connec, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))


# relation entre la masse des racines et la masse vdt
#------------------------------------------------------------------
plot(rac$biomR ~rac$biomFin)
par(mfrow=c(1,1))
xlim=c(0,7)
ylim=c(4,9.5)
plot(rac$biomR[rac$div=="Div0" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Final earthworm biomass (g DM)", ylab="Total root mass (g DM)", las=1)
points(rac$biomR[rac$div=="Div0" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(biomR~biomFin, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
#legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))
# on regarde la regression seulement pour div 3
reg2 <- lm(biomR~biomFin, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
abline(reg2, col="purple")
text(6,7, "R²=0.28.", col="purple")
reg3 <- lm(biomR~biomFin, data=rac[rac$div=="Div3",], na.action=na.omit)
shapiro.test(residuals(reg3))
summary(reg3)


xlim=c(-100,0)
plot(rac$biomR[rac$div=="Div0" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Earthworm delta biomass", ylab="Total root mass (g DM)", las=1)
points(rac$biomR[rac$div=="Div0" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(biomR~biomDelta, data=rac, na.action=na.omit)
shapiro.test(residuals(mod))
summary(reg)
# on ritre Div0
reg0 <- lm(biomR~biomDelta, data=rac[rac$div!="Div0",], na.action=na.omit)
shapiro.test(residuals(reg0))
summary(reg0)



# relation entre la masse des racine et SBL
#------------------------------------------------
plot(biomR~SBL, data=rac, xlab="Specific burrow length (m/g MF)", ylab="biomRueur racines (m)")
xlim=c(0,25)
ylim=c(4,9.5)
plot(rac$biomR[rac$div=="Div0" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Specific burrow length (m/g FM)", ylab="Total root mass (g DM)", las=1)
points(rac$biomR[rac$div=="Div0" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(biomR~SBL, data=rac, na.action=na.omit)
summary(reg)
# on retire Div0
reg <- lm(biomR~SBL, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg)
reg1 <- lm(biomR~SBL, data=rac[rac$div=="Div1",], na.action=na.omit)
summary(reg1)
reg2 <- lm(biomR~SBL, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
abline(reg2, col="purple")
text(22,6.5,"R²=0.32*", col="purple")
reg3 <- lm(biomR~SBL, data=rac[rac$div=="Div3",], na.action=na.omit)
summary(reg3)

############################################################################################################################
############################################################################################################################

##  Relation entre le SRL et la longueur des galeries
#------------------------------------------------------------------
## tous traitements confondus
plot(SRL~gal, data=rac, xlab="Longueur totale des galeries (m)", ylab="SRL (m/g MS)")
text(40, 6.5, "R²=0.61***")


xlim=c(0,50)
ylim=c(4,10)
plot(rac$SRL[rac$div=="Div0" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Total burrow length (m)", ylab="Specific root length (m/g DM)", las=1)
points(rac$SRL[rac$div=="Div0" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg-"]~rac$gal[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg+"]~rac$gal[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(SRL~gal, data=rac, na.action=na.omit)
summary(reg)
abline(reg)
text(40, 6.5, "R²=0.61***")
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))
reg0 <- lm(SRL~gal, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg0)
# abline(reg0, lty="dashed")
reg1 <- lm(SRL~gal, data=rac[rac$div=="Div1",], na.action=na.omit)
summary(reg1)
abline(reg1, col="orange")
text(10, 8.5, "R²=0.29.", col="orange")
reg2 <- lm(SRL~gal, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
abline(reg2, col="purple")
text(10, 6, "R²=0.39*", col="purple")
reg3 <- lm(SRL~gal, data=rac[rac$div=="Div3",], na.action=na.omit)
summary(reg3)

##  Relation entre le SRL et la connectance
#------------------------------------------------------------------
## tous traitements confondus
plot(SRL~connec, data=rac)
xlim=c(0,0.5)
ylim=c(4,10)
plot(rac$SRL[rac$div=="Div0" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Burrow connectance", ylab="Specific root length (m/g DM)", las=1)
points(rac$SRL[rac$div=="Div0" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg-"]~rac$connec[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg+"]~rac$connec[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(SRL~connec, data=rac, na.action=na.omit)
summary(reg)
# abline(reg)
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))
# on fait la reg sans les div0
reg0 <- lm(SRL~connec, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg0)
reg1 <- lm(SRL~connec, data=rac[rac$div=="Div1",], na.action=na.omit)
summary(reg1)
abline(reg1, col="orange")
text(0.45,7.5, "R²=0.34*", col="orange")
reg2 <- lm(SRL~connec, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
reg3 <- lm(SRL~connec, data=rac[rac$div=="Div3",], na.action=na.omit)
summary(reg3)
reg <- lm(SRL~connec, data=rac[rac$leg=="Leg+",], na.action=na.omit)
summary(reg)
# abline(reg, col="grey")
reg <- lm(SRL~connec, data=rac[rac$leg=="Leg-",], na.action=na.omit)
summary(reg)
# abline(reg, col="red")

##  Relation entre le SRL et la SBL (Specific Burrow Lenght)
#------------------------------------------------------------------
## tous traitements confondus
plot(SRL~SBL, data=rac, xlab="Specific burrow length (m/g MF)", ylab="SRL (m/g MS)")
xlim=c(0,25)
ylim=c(4,10)
plot(rac$SRL[rac$div=="Div0" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Specific burrow length (m/g FM)", ylab="Specific root length (m/g DM)", las=1)
points(rac$SRL[rac$div=="Div0" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg-"]~rac$SBL[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg+"]~rac$SBL[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(SRL~SBL, data=rac, na.action=na.omit)
shapiro.test(residuals(reg))
summary(reg)
abline(reg)
text(22,6,"R²=0.34***")
reg0 <- lm(SRL~SBL, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg0)
reg1 <- lm(SRL~SBL, data=rac[rac$div=="Div1",], na.action=na.omit)
summary(reg1)
reg2 <- lm(SRL~SBL, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
reg3 <- lm(SRL~SBL, data=rac[rac$div=="Div3",], na.action=na.omit)
summary(reg3)


##  Relation entre le SRL et la masse finale des vdt
#------------------------------------------------------------------
## tous traitements confondus
plot(SRL~biomFin, data=rac)
xlim=c(0,7)
ylim=c(4,10)
plot(rac$SRL[rac$div=="Div0" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Final earthworm biomass (g FM)", ylab="Specific root length (m/g DM)", las=1)
points(rac$SRL[rac$div=="Div0" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg-"]~rac$biomFin[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg+"]~rac$biomFin[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(SRL~biomFin, data=rac, na.action=na.omit)
summary(reg)
abline(reg)
text(6,6, "R²=0.43***")
legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))
# on fait la reg sans les div0
reg <- lm(SRL~connec, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg)
reg1 <- lm(SRL~connec, data=rac[rac$div=="Div1",], na.action=na.omit)
summary(reg1)
reg2 <- lm(SRL~connec, data=rac[rac$div=="Div2",], na.action=na.omit)
summary(reg2)
reg3 <- lm(SRL~connec, data=rac[rac$div=="Div3",], na.action=na.omit)
summary(reg3)


xlim=c(-100,0)
ylim=c(4,10)
plot(rac$SRL[rac$div=="Div0" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Earthworm biomass change (%)", ylab="Specific root length (m/g DM)", las=1)
points(rac$SRL[rac$div=="Div0" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div1" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div2" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg-"]~rac$biomDelta[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$SRL[rac$div=="Div3" & rac$leg=="Leg+"]~rac$biomDelta[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
reg <- lm(SRL~biomDelta, data=rac, na.action=na.omit)
summary(reg)
# abline(reg)
# text(0.2,8,"R²=0.35***")
# on teste la reg sans div0
reg <- lm(SRL~biomDelta, data=rac[rac$div!="Div0",], na.action=na.omit)
summary(reg)
abline(reg, lty="dashed")

############################################################################################################################


############################################################################################################################

### relations entre traits racinaires
#-------------------------------------

# Relation entre la longeur des racines sur la masse racinaire
#-----------------------------------------------------------
rac$biomR <- rac$DW/1000 # masse racinaire en g
summary(rac)

reg <-lm(biomR~long, data = rac) 
summary(reg)
ylim=c(4,10)
xlim=c(25,65)
plot(rac$biomR[rac$div=="Div0" & rac$leg=="Leg-"]~rac$long[rac$div=="Div0" & rac$leg=="Leg-"], pch=16, cex=1.5, col="black", 
     xlim=xlim, ylim=ylim, xlab="Total root length (m)", ylab="Total root mass (g DM)", las=1)
points(rac$biomR[rac$div=="Div0" & rac$leg=="Leg+"]~rac$long[rac$div=="Div0" & rac$leg=="Leg+"], pch=17, cex=1.5, col="black")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg-"]~rac$long[rac$div=="Div1" & rac$leg=="Leg-"], pch=16, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div1" & rac$leg=="Leg+"]~rac$long[rac$div=="Div1" & rac$leg=="Leg+"], pch=17, cex=1.5, col="orange")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg-"]~rac$long[rac$div=="Div2" & rac$leg=="Leg-"], pch=16, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div2" & rac$leg=="Leg+"]~rac$long[rac$div=="Div2" & rac$leg=="Leg+"], pch=17, cex=1.5, col="purple")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg-"]~rac$long[rac$div=="Div3" & rac$leg=="Leg-"], pch=16, cex=1.5, col="green")
points(rac$biomR[rac$div=="Div3" & rac$leg=="Leg+"]~rac$long[rac$div=="Div3" & rac$leg=="Leg+"], pch=17, cex=1.5, col="green")
# abline (reg, col="grey")
# text(60,7, "R²=0.27***")

reg0 <-lm(rac$biomR[rac$div=="Div0"]~rac$long[rac$div=="Div0"]-1, data = rac) 
summary(reg0)
abline(reg0, col='black')
text(60,7,"R²=0.98***")
reg1 <-lm(rac$biomR[rac$div=="Div1"]~rac$long[rac$div=="Div1"]-1, data = rac) 
summary(reg1)
abline(reg1, col='orange')
text(60,9.2,"R²=0.98***", col='orange')
reg2 <-lm(rac$biomR[rac$div=="Div2"]~rac$long[rac$div=="Div2"]-1, data = rac) 
summary(reg2)
abline(reg2, col='purple')
text(50,9.5,"R²=0.98***", col='purple')
reg3 <-lm(rac$biomR[rac$div=="Div3"]~rac$long[rac$div=="Div3"]-1, data = rac) 
summary(reg3)
abline(reg3, col='green')
text(57,9.7,"R²=0.99***", col='green')

legend("topright", legend=c("Leg-","Leg+","Div0","Div1","Div2","Div3"), pch=c(1,2,16,16,16,16), col=c("black","black","black","orange","purple","green"))


mod <- lm(rac$biomR~rac$long-1+rac$div)
summary(mod)
anova(mod)




