#------------------------------------------------------------
# Cr?ation d'une fonction pour repr?senter les diff?rences entre modalit?s
# Les arguments sont :  
# mdata, df contenant les donn?es de calcul des mod?les
# mod => mod?le selectionn? (il faut utiliser le $modxx g?n?r? 
#       par la fonction myfun_Com)
# ylbl => ylabel des graphiques
# plot_PAR => T/F pour repr?senter un effet du PAR    
# Les outputs sont : 
# moyFRic => moyenne et se de vb par modalit?s du facteur 'div'
# moyLeg => moyenne et se de vb par modalit?s du facteur 'leg'
# moyInteraction => moyenne et se de vb par modalit?s crois?es des facteurs 'div et 'leg'
# pDIV => graph des moy (+/- se) de vb en fonction des modalit?s du facteur 'div'
# pDIV_PAR => droite de r?ponse de vb en fonction du PAR (modulo les niveaux de 'div')
# pDIV_INTER => graph des moy (+/- se) de vb en fonction des modalit?s crois?es des facteurs 'div et 'leg'
# pLEG => graph des moy (+/- se) de vb en fonction des modalit?s du facteur 'leg' 
# pLEG_PAR => droite de r?ponse de vb en fonction du PAR (modulo les niveaux de 'leg')
#------------------------------------------------------------

# Librairies
librarian::shelf(grid, plotrix, gridExtra)

# A l'?chelle de la communaut? vdt  
myplot_Com_vdt <- function(mdata, vb, mod, ylabl, plot_PAR){
  mod <- eval(mod)
  newdata <- data.frame(mdata) 
  newdata$predict <- predict(mod)
  
  # Effet FRic 
  moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
  seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
  moyF$se <- seF$x
  
  DIV <- ggplot(data = moyF, aes(x = div,y=x, fill=div)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
  
  DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    ylab(ylabl)+
    theme(text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    xlab("PAR (?mol/m?/s)")
  
  # Effet leg 
  moyL <- aggregate(newdata$predict, by = list(leg = newdata$leg), mean, na.rm = T)
  seL <- aggregate(newdata$predict, by = list(leg = newdata$leg), std.error, na.rm = T)
  moyL$se <- seL$x
  LEG <- ggplot(moyL, aes(x = leg, y = x)) + 
    geom_pointrange(aes(ymin = x-1.96*se, ymax = x+1.96*se))+
    theme(text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    ylab("")
  
  LEG_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = leg, col = leg)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    theme(text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    xlab("PAR (?mol/m?/s)")
  
  # Effet interaction
  moyFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), mean, na.rm=T)
  seFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), std.error, na.rm=T)
  moyFi$se <- seFi$x
  
  INTER <- ggplot(data = moyFi, aes(x = leg,y=x, fill=div, alpha = leg)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))  +
    facet_grid(div~.)
  
  #if (plot_PAR == FALSE) {A <- grid.arrange(DIV, LEG, INTER, ncol = 3)
  #}else{A <- grid.arrange(DIV, LEG, INTER, DIV_PAR, LEG_PAR, ncol = 3, nrow = 2)
  #}
  # A <- arrangeGrob(DIV, LEG, INTER, DIV_PAR, LEG_PAR, ncol = 3, nrow = 2)
  # grid.draw(A)
  z <- list(moyFRic = moyF, moyLeg = moyL, moyInter = moyFi, pDIV = DIV, pINTER = INTER, 
            pLEG = LEG, pDIV_PAR = DIV_PAR, pLEG_PAR = LEG_PAR)
  return(z)
  
  }

# A l'?chelle de la communaut? plante
myplot_Com <- function(mdata, vb, mod, ylabl, plot_PAR){
  mod <- eval(mod)
  newdata <- data.frame(mdata) 
  newdata$predict <- predict(mod)
  
  # Effet FRic 
  moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
  seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
  moyF$se <- seF$x
  DIV <- ggplot(data = moyF, aes(x = div, y = x, shape = div)) + 
    geom_point(size = 4)+
    geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
    ylab(ylabl) +    
    xlab("Earthworm FRic") +
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  
  DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    ylab(ylabl)+
    theme(text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    xlab("PAR (?mol/m?/s)")
  
  # Effet leg 
  moyL <- aggregate(newdata$predict, by = list(leg = newdata$leg), mean, na.rm = T)
  seL <- aggregate(newdata$predict, by = list(leg = newdata$leg), std.error, na.rm = T)
  moyL$se <- seL$x
  LEG <- ggplot(moyL, aes(x = leg, y = x)) + 
    geom_pointrange(aes(ymin = x-1.96*se, ymax = x+1.96*se))+
    theme(text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    ylab("")
  
  LEG_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = leg, col = leg)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    theme(text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    xlab("PAR (?mol/m?/s)")
  
  # Effet interaction
  moyFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), mean, na.rm=T)
  seFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), std.error, na.rm=T)
  moyFi$se <- seFi$x
  INTER <- ggplot(data = moyFi, aes(x = div, y = x)) + 
    geom_point()+
    geom_pointrange(aes(ymin = x - 1.96*se, ymax = x + 1.96*se))+
    ylab(ylabl)+
    facet_grid(leg~.)+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  z <- list(moyFRic = moyF, moyLeg = moyL, moyInter = moyFi, pDIV = DIV, pINTER = INTER, 
            pLEG = LEG, pDIV_PAR = DIV_PAR, pLEG_PAR = LEG_PAR)
  z
  
}

# A l'?chelle de l'esp?ce
myplot_Ind <- function(mdata, vb, mod, ylabl, plot_PAR){
  mod <- eval(mod)
  newdata <- data.frame(mdata) 
  newdata$predict <- predict(mod)
  
  # Effet FRic 
  moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
  seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
  moyF$se <- seF$x
  DIV <- ggplot(data = moyF, aes(x = div, y=x, fill=div)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("Ew0", "FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
  
  
  
  DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    ylab(ylabl)+
    xlab("PAR (?mol/m?/s)")+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  moyFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), mean, na.rm=T)
  seFi <- aggregate(newdata$predict, by = list(div = newdata$div, leg = newdata$leg), std.error, na.rm=T)
  moyFi$se <- seFi$x
  
  DIV_INTER <- ggplot(data = moyFi, aes(x = leg, y=x, fill=div)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("Ew0", "FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
  
  
  # Effet leg 
  moyL <- aggregate(newdata$predict, by = list(leg = newdata$leg), mean, na.rm = T)
  seL <- aggregate(newdata$predict, by = list(leg = newdata$leg), std.error, na.rm = T)
  moyL$se <- seL$x
  LEG <- ggplot(moyL, aes(x = leg, y = x)) + 
    geom_pointrange(aes(ymin = x-1.96*se, ymax = x+1.96*se))+
    ylab(ylabl)+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  LEG_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = leg, col = leg)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    xlab("PAR (?mol/m?/s)")+
    ylab(ylabl)+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  if (plot_PAR == FALSE) {A <- grid.arrange(DIV, LEG, DIV_INTER, ncol = 2, nrow = 2)
  }else{A <- grid.arrange(DIV, DIV_PAR, DIV_INTER, LEG, LEG_PAR, ncol = 3, nrow = 2)
  }
  print(A)
  z <- list(moyFRic = moyF, moyInteraction = moyFi, moyLeg = moyL, 
            pDIV = DIV, pDIV_PAR = DIV_PAR, pDIV_INTER = DIV_INTER, pLEG = LEG, pLEG_PAR = LEG_PAR)
  z
}

# A l'?chelle de l'esp?ce pour le tr?fle      
myplot_Leg <- function(mdata, vb, mod, ylabl, plot_PAR){
  mod <- eval(mod)
  newdata <- data.frame(mdata) 
  newdata$predict <- predict(mod)
  
  # Effet FRic 
  moyF <- aggregate(newdata$predict, by = list(div = newdata$div), mean, na.rm=T)
  seF <- aggregate(newdata$predict, by = list(div = newdata$div), std.error, na.rm=T)
  moyF$se <- seF$x
  
  DIV <- ggplot(data = moyF, aes(x = div, y = x, fill = div)) + 
    geom_bar(position = position_dodge(width = 0.8), stat="identity")+
    geom_errorbar(aes(ymin = x - 1.96*se, ymax = x + 1.96*se), width = 0, position=position_dodge(width = 0.8))+
    ylab(ylabl) +    
    xlab("") +
    scale_fill_manual(values = c("black", "#ffeda0", "#feb24c","#f03b20"), 
                      name = "",
                      labels = c("Ew0","FR1", "FR2", "FR3"))+
    theme_classic() +
    theme(panel.background = element_rect(fill = "white", colour = "black", size=0.5),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))+
    theme(axis.text = element_text(color="black", size=14)) +
    theme(axis.title = element_text(color="black", size=14))+
    scale_x_discrete(labels=c("Div0"="Ew0", "Div1" = "FR1", "Div2" = "FR2","Div3" = "FR3"))
  
  
  DIV_PAR <- ggplot(data = newdata, aes(x = PAR, y = predict, group = div, col = div)) + 
    geom_point()+
    geom_smooth(method = "lm")+
    ylab(ylabl)+
    xlab("PAR (?mol/m?/s)")+
    theme(panel.background = element_rect(fill="white", colour = "black"),
          text =  element_text(face = "plain",
                               color = "black", 
                               hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9,
                               margin = margin(), debug = FALSE))
  
  if (plot_PAR == FALSE) {A <- DIV
  }else{A <- grid.arrange(DIV, DIV_PAR, ncol = 2, nrow = 1)
  }
  print(A)
  z <- list(moyFRic = moyF, pDIV = DIV, pDIV_PAR = DIV_PAR)
  z
}
#------------------------------------------------------------