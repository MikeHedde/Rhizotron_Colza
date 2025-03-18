#######################################"

# humidité du sol - rhizotron

# MZ 
# 21/06/2017
###########################################

setwd("C:/Users/MZ/Documents/PRO/JENNY/RHIZOTRON/Données/sol/enzym")
tab <- read.table ("BCE-BDD-JENNY-GLOBAL-20160204.csv", sep=";", dec=".", header=TRUE, na.strings="")
summary(tab)

# on retire les RT douteux
tab$HUM[tab$CODE=="RT_8" | tab$CODE=="RT_12" | tab$CODE=="RT18" | tab$CODE=="RT40"] = NA


## effet profondeur et effeet traitement sur les humidité du sol
tab$PLAN <- as.factor(tab$PLAN)
tab$STRATE <- as.factor(tab$STRATE)
tab$VDT <- as.factor(tab$VDT)
mod <- lm(asin(sqrt(HUM/100))~STRATE*PLAN*VDT, data=tab, na.action=na.omit)
shapiro.test(residuals(mod))
fligner.test(asin(sqrt(HUM/100))~STRATE, data=tab, na.action=na.omit)
fligner.test(asin(sqrt(HUM/100))~PLAN, data=tab, na.action=na.omit)
fligner.test(asin(sqrt(HUM/100))~VDT, data=tab, na.action=na.omit)
anova(mod)
library(multcomp)
mod <- lm(asin(sqrt(HUM/100))~STRATE+VDT, data=tab, na.action=na.omit)
anova(mod)
comp <- glht(mod,linfct=mcp(STRATE="Tukey") )
summary(comp)
mod <- lm(asin(sqrt(HUM/100))~VDT, data=tab, na.action=na.omit)
anova(mod)
comp <- glht(mod,linfct=mcp(VDT="Tukey") )
summary(comp)

a <- glht(mod, linfct=mcp(STRATE="Tukey"))
summary(a)


kruskal.test(asin(sqrt(HUM/100))~STRATE, data=tab, na.action=na.omit)
library(PMCMR)
posthoc.kruskal.nemenyi.test(x=asin(sqrt(tab$HUM/100)), g=tab$STRATE, dist="Tukey")
# strate 4 différentes de autres
kruskal.test(asin(sqrt(HUM/100))~PLAN, data=tab, na.action=na.omit)
kruskal.test(asin(sqrt(HUM/100))~VDT, data=tab, na.action=na.omit)


##### humidité moyenne dans les RT pour chaque strate en fonction des traitements
moy <- aggregate(tab$HUM, by=list(strate=tab$STRATE, leg=tab$PLAN, div=tab$VDT), mean, na.rm=T)
library(plotrix)
se <- aggregate(tab$HUM, by=list(strate=tab$STRATE, leg=tab$PLAN, div=tab$VDT), std.error, na.rm=T)
moy$se <- se$x

## représentation graphique avec tous les traitements
library(ggthemes) 
p <- ggplot(moy, aes(x=strate, y=x, fill="white")) + 
  geom_bar(position=position_dodge(), stat="identity", 
           colour="black", # use black outlines, 
           size=0.3)+   # for thinner lines
  geom_errorbar(aes(ymin=x-se, ymax=x+se),
                size=0.3, # for thinner lines
                width=0.2, # Width of the error bars
                position=position_dodge(.9))+
  ylab("Soil moisture (%)")+
  ylim(0,15) +
  scale_fill_manual(values=c("white")) + # pour changer les couleurs des barres
  facet_grid( div ~ leg) +
  theme_tufte(base_size = 15, base_family ="Arial") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
p

#theme_bw() ; theme_tufte() ;   theme_classic()+

## représentation graphique pour effet strates
moy2 <- aggregate(tab$HUM, by=list(strate=tab$STRATE), mean, na.rm=T)
se2 <- aggregate(tab$HUM, by=list(strate=tab$STRATE), std.error, na.rm=T)
moy2$se <- se2$x

p <- ggplot(moy2, aes(x=strate, y=x)) + # fill = "white"
  geom_bar(position=position_dodge(), stat="identity", 
           colour="black", # use black outlines, 
           size=0.3)+   # for thinner lines
  geom_errorbar(aes(ymin=x-se, ymax=x+se),
                size=0.3, # for thinner lines
                width=0.2, # Width of the error bars
                position=position_dodge(.9))+
  ylab("Soil moisture (%)")+
  ylim(0,15) +
  theme_tufte(base_size = 15, base_family ="sans") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
p
