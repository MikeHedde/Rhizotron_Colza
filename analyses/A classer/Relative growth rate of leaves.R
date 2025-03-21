setwd("G:/INRA/5. RH/4. Post-doctorant/JENNY/RT")

require(ggplot2)
require(vegan)
library(PerformanceAnalytics)


#Chargement jeu de donn�es => il faut les valeurs de longueur de feuille en croissance
#(i.e. sans la s�nescence)
bn0<-read.csv("lg.csv", h=T, sep=";")
lg2<-read.csv("lg2.csv", h=T, sep=";")


#############################################################

#creation d'une fonction pour tester si les donn�es de longeur d'une feuille 
#sont croissantes dans le temps. Sort une liste des feuilles pour lesquelles 
#il y a des donn�es abh�rentes (seuil = tol [% de longueur inf�rieure entre j et j+1]) 
#et cr��e une nouvelle data.frame sans les plantes pour lesquelles les donn�es 
#ne sont pas consolid�es
VERIF<-function(mdat, tol){
  fac<-aggregate(mdat, list(mdat$loca, mdat$RT, mdat$div, mdat$leg, 
                            mdat$leaf, mdat$sp), mean)[,1:6]
  colnames(fac)<-c("loca", "RT", "div", "leg", "leaf", "sp")
  
  verif<-matrix(NA, nrow=nrow(fac), ncol=1)
  
  
  for(i in 1:nrow(fac)){
    sub<-mdat[which(mdat$loca==fac$loca[i] & mdat$leaf==fac$leaf[i]),]
    delta<-matrix(NA, nrow = nrow(sub)-1, ncol=1)
    for(j in 1:nrow(delta))
      delta[j]<-ifelse((sub$LL[j+1]-sub$LL[j])/sub$LL[j]*100 < tol, -1, 0)
    verif[i]<-sum(delta)
  }
  
  fac$verif<-verif[,1]
  tbc<-table(fac[fac$verif!=0,c(1,5)])
  to_be_check<-tbc[rowSums(tbc)>0,]
  print(to_be_check)
  mdat
}


#############################################################

#creation d'une fonction pour ajuster les donn�es de longeur d'une feuille 
#� une courbe sigmo�de
SIGM<-function(mdat){
  fit <- nls(LL ~ theta1/(1 + exp(-(theta2 + theta3*day))),
             start=list(theta1 = 7, theta2 = -13, theta3 = 0.5),
             data=mdat, algo="port")
  
  
  #extraction des param�tres du mod�le  
  LLm<-max(mdat$LL)
  T1 <- coef(fit)[[1]]
  T2 <- coef(fit)[[2]]
  T3 <- coef(fit)[[3]]
  ddl <-  df.residual(fit)
  Res.dev <- summary(fit)[[4]] 
  pval <- 1-pchisq(Res.dev,ddl) 
  ymoy <- mean(mdat$LL)
  yobs <- mdat$LL
  SStot <- sum((yobs-ymoy)^2)
  SSerr <- sum(residuals(fit)^2)
  Rsquare <- 1 - (SSerr/SStot)
  
  #cr�ation d'une data.frame contenant les valeurs mod�lis�es de la longeur de la feuille
  #en fonction du temps
  sigm1 <- data.frame(day = seq(0, 41, by=1))
  sigm1$LL <- predict(fit, newdata = sigm1)
  
  #cr�ation d'une data.frame contenant les valeurs mod�lis�s de la vitesse de croissance 
  # de la feuille pour chacun des intervelle de jours de la p�riode d'observation
  sigm2 <- data.frame(dayinter = seq(1, 40, by=1))
  for(i in 1:nrow(sigm2))
    sigm2$LLd[i]<-sigm1$LL[sigm2$dayinter[i]+1]-sigm1$LL[sigm2$dayinter[i]]
  
  #cr�ation d'un objet avec la liste des param�tres
  resultfit<-list()
  resultfit$Lmax <- LLm
  resultfit$T1 <- T1
  resultfit$T2 <- T2
  resultfit$T3 <- T3
  resultfit$pval <- pval
  resultfit$R2 <- Rsquare
  resultfit$LLd<-sigm2$LLd
  
  resultfit
}

#############################################################

#cr�ation d'une fonction pour calculer simultan�ment les param�tres de r�gressions  
# de la longueur de plusieurs feuilles en fonction du temps 
FUN <- function(DF, loc, leafn){
  
  mdat <- DF[which(DF$loca==loc & DF$leaf==leafn ),]
  
  #creation d'une liste de param�tres avec valeur NA pour les cas ou il n'est pas possible
  #d'ajuster les donn�es � une courbe sigmo�de (nb observation inf � 4)
  resultnotfit <- list()
  resultnotfit$Lmax <- NA
  resultnotfit$T1 <- NA
  resultnotfit$T2 <- NA
  resultnotfit$T3 <- NA
  resultnotfit$pval <- NA
  resultnotfit$R2 <- NA
  resultnotfit$LLd <- rep(NA, 40)
  
  
  #code pour tester le nombre d'observation et attribuer �ventuellement des valeurs NA par d�faut
  options(show.error.messages = FALSE)
  if(length(unique(mdat$LL))>3){
    tt<-try(SIGM(mdat))
    ifelse(class(tt) != "try-error",
           resultfit<-SIGM(mdat),
           resultfit<-resultnotfit)
  }
  else 
    resultfit<-resultnotfit
  resultfit  
}


#############################################################

#produire les matrices de param�tres des mod�les sigmo�des et des vitesses  
#de croissance journali�res pour toutes les feuilles d'un jeu de donn�es en exemple

MODLEAF<-function(mydata){
  
  library(reshape2, quietly =  T)
  
  #Cr�ation d'une matrice pour accueillir les param�tres des mod�les
  fac<-aggregate(mydata, list(mydata$loca, mydata$RT, mydata$div, mydata$leg, 
                              mydata$leaf, mydata$sp), mean)[,1:6]
  parm<-matrix(NA, ncol=6, nrow=nrow(fac))
  fitParam<-cbind(fac, parm)
  colnames(fitParam)<-c("loca", "RT", "div", "leg", "leaf", "sp", "Lmax", "T1", "T2", "T3", "pval", "R2")
  
  #Cr�ation d'une matrice pour accueillir les donn�es mod�lis�es de vitesse de croissance
  fitGrwSpeed <- matrix(NA, ncol = 40, nrow = nrow(fac)) 
  colnames (fitGrwSpeed) <- 1:40
  
  #Collecte et insertion des param�tres des mod�les ("fitParam") et des vitesses de croissance 
  #journali�res ("fitGrwSpeed")pour toutes les feuilles d'un jeu de donn�es
  for(i in 1 : nrow(fitParam))  {
    k <- FUN(DF = mydata, loc = fitParam$loca[i], leafn = fitParam$leaf[i] )
    fitParam[i,]$Lmax<-k$Lmax
    fitParam[i,]$T1<-k$T1
    fitParam[i,]$T2<-k$T2
    fitParam[i,]$T3<-k$T3
    fitParam[i,]$pval<-k$pval[1]
    fitParam[i,]$R2<-k$R2
    fitGrwSpeed[i,]<-k$LLd
    fitGrwSpeed<-as.data.frame(fitGrwSpeed)
  }
  
  fitGrwSpeed <- cbind(fac, fitGrwSpeed)
  colnames(fitGrwSpeed)[1:6] <- c("loca", "RT", "div", "leg", "leaf", "sp")
  fitGS <- melt(fitGrwSpeed, id = c("loca", "RT", "div", "leg", "leaf", "sp"))
  colnames(fitGS)[7:8] <- c("day", "GrwSpd")
  class(fitGS$day)<-"numeric"
  class(fitGS$div)<-"numeric"
  
  moddata<-list()
  moddata$fitGS<-fitGS
  moddata$fitParam<-fitParam
  moddata
}


############################################################################
############################################################################
lg2_verif<-VERIF(mdat = lg2, tol = "-5")
mod<-MODLEAF(bn0)

RGR <- mod$fitGS
PAR <- mod$fitParam
RGR <- RGR[is.na(RGR$GrwSpd) == FALSE, ]
PAR <- PAR[is.na(PAR$Lmax) == FALSE,]

############################################################################
#calcul des indices (somme, moyenne, max, CV) des RGR individuels des feuilles d'un individu pour un intervalle donn�
RGR_indice<-aggregate(RGR$GrwSpd, list(RGR$loca, RGR$RT, RGR$div, RGR$leg, RGR$day,
                                       RGR$sp), sum)
colnames(RGR_indice)<-c("loca", "RT", "div", "leg", "day", "sp", "sum")

RGR_indice$mean <- aggregate(RGR$GrwSpd, list(RGR$loca, RGR$RT, RGR$div, RGR$leg, RGR$day,
                                              RGR$sp), mean)$x
RGR_indice$max <- aggregate(RGR$GrwSpd, list(RGR$loca, RGR$RT, RGR$div, RGR$leg, RGR$day,
                                             RGR$sp), max)$x
RGR_indice$sd <- aggregate(RGR$GrwSpd, list(RGR$loca, RGR$RT, RGR$div, RGR$leg, RGR$day,
                                             RGR$sp), sd)$x
RGR_indice$cv <- RGR_indice$sd/RGR_indice$mean
RGR_indice$ngrowing <- aggregate(RGR$GrwSpd[RGR$GrwSpd>.1], 
                          list(RGR$loca[RGR$GrwSpd>.1], RGR$RT[RGR$GrwSpd>.1], 
                               RGR$div[RGR$GrwSpd>.1], RGR$leg[RGR$GrwSpd>.1], 
                               RGR$day[RGR$GrwSpd>.1], RGR$sp[RGR$GrwSpd>.1]), 
                          specnumber)$x
RGR_indice$Hgrowing <- aggregate(RGR$GrwSpd, list(RGR$loca[RGR$GrwSpd>.1], RGR$RT[RGR$GrwSpd>.1], 
                               RGR$div[RGR$GrwSpd>.1], RGR$leg[RGR$GrwSpd>.1], 
                               RGR$day[RGR$GrwSpd>.1], RGR$sp[RGR$GrwSpd>.1]), diversity)$x

RGR_indagg<-aggregate(RGR_indice[,-c(1:6)], list(RGR_indice$div, RGR_indice$leg, RGR_indice$day,
                                                 RGR_indice$sp), mean)
colnames(RGR_indagg)[1:4]<-c("div", "leg", "day", "sp")

chart.Correlation(RGR_indices[,-c(1:6)], histogram=TRUE, pch=19)
chart.Correlation(RGR_indagg[,c(5:8)], histogram=TRUE, pch=19)

############################################################################
#repr�sentation des RGR individuels des feuilles d'un individu
p <- ggplot(RGR[which(RGR$sp=="Bra_napus"),], aes(x=day, y=10*GrwSpd, col=leg)) +
  geom_point() +
  facet_grid(leaf~div, scales="fixed") +
  labs(y = "Relative Growth Rate (RGR, mm d-1)") 
p



############################################################################
#repr�sentation des indices des RGR des feuilles d'un individu
p <- ggplot(RGR_indice[which(RGR_indice$sp=="Bra_napus"),], aes(x=day, y=max)) + 
  geom_point() +
  facet_grid(leg~div, scales="fixed") +
  ylab("max Relative Growth Rate (RGR, mm d-1)")
p


p <- ggplot(PAR, aes(x=div, y=10*Lmax)) +  
     geom_point() +
     geom_smooth(aes(col = leg), method = "lm", se = T, 
                 legendTitle = "Legume \nTreatment") +
     facet_grid(leaf~leg, scales="fixed", margin = T) + 
     labs(y = "max Leaf Length (LL, mm)", 
          x = "Function Richness of earthworm assemblages (FRic)") +
     scale_color_discrete (name = "Legume \nTreatment", 
                          labels = c(italic("-T. alexandrium"), italic("+T. alexandrium"), "all"))
p


############################################################################
#repr�sentation des indices par traitement exp�rimental
p <- ggplot(RGR_indagg[which(RGR_indagg$sp=="Bra_napus"),], aes(x=day, y=mean, col=leg)) + 
  geom_point() +
  facet_grid(~div, scales="fixed") +
  geom_ribbon(aes(ymin=mean-.99*(mean-sd), ymax=mean+.99*(mean-sd)), alpha=0.3) +
  ylab("mean Relative Growth Rate (mm d-1)")
p

p <- ggplot(RGR_indagg[which(RGR_indagg$sp=="Bra_napus"),], aes(x=day, y=sum, col=leg)) + 
  geom_line() +
  facet_grid(~div, scales="fixed") +
  ylab("cumulated Relative Growth Rate  (mm d-1)")
p


