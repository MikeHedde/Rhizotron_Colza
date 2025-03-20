librarian::shelf(tidyverse, dplyr, ggpubr, Rmisc, ade4)


##################    17/03/2025  ######################
myenz <- read.csv("data/derived-data/enz.csv", h = T, sep = ",") 
anion <- read.csv("data/raw-data/anion.csv", h = T, sep = ";") 

nutrient <- left_join(myenz, anion)

myenz_long <- myenz %>%
  mutate(ratioCN = log(glu)/log(ure),
         ratioCP = log(glu)/log(phos),
         ratioCS = log(glu)/log(ary)) %>%
  pivot_longer(cols = c(phos, ary, glu, ure, ratioCN, ratioCP, ratioCS), 
               names_to = "variable", values_to = "value") %>%
  select(-(c(X, serie, position, PAR))) %>%
  summarySE(measurevar="value", groupvars=c("leg", "div", "strate" ,"variable")) %>%
  mutate(strate = case_when(strate == 1 ~-3.5,
                            strate == 2 ~-18.5,
                            strate == 3 ~-40.5,
                            strate == 4 ~-70.5)) 
ratio <- myenz_long %>%
  filter(variable %in% c("ratioCN", "ratioCP", "ratioCS"),
         !is.na(sd))

ggplot(ratio, aes(x = strate, y = value, colour = div))+
  geom_point(position=position_dodge(width= 4))+
  geom_line(aes(linetype = div))+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), 
                width=.1, position=position_dodge(width= 4))+
  coord_flip()+
  facet_grid(leg~variable, scales = "free")+
  theme_bw()

# BCA 
nutrient_pca_df <- nutrient %>%
  filter(glu != 0) %>%
  mutate(ratioCN = log(glu)/log(ure),
         ratioCP = log(glu)/log(phos),
         ratioCS = log(glu)/log(ary),
         code_legstrate = as.factor(paste(leg, "(", strate, ")", sep = "")),
         code_divstrate = as.factor(paste(div, "(", strate, ")", sep = ""))) %>%
  select(code_legstrate, code_divstrate, glu, ratioCP, ratioCN, ratioCS, SO42, NO3, HPO42) 

nutrient_pca <- dudi.pca(nutrient_pca_df[, -c(1:2)], scannf = FALSE, nf = 4)
nutrient_bca_leg <- bca(nutrient_pca, as.factor(nutrient_pca_df$code_legstrate), scan = FALSE, nf = 2)
nutrient_bca_div <- bca(nutrient_pca, as.factor(nutrient_pca_df$code_divstrate), scan = FALSE, nf = 2)

nutrient_bca_leg$ratio
nutrient_bca_div$ratio

leg_eig <- nutrient_bca_leg$eig/sum(nutrient_bca_leg$eig)*100
div_eig <- nutrient_bca_div$eig/sum(nutrient_bca_div$eig)*100

nutrient_bca_leg_baryc <- nutrient_bca_leg$li[,1:2] %>%
  as_tibble(rownames = "code")
nutrient_bca_leg_fac <- nutrient_bca_leg$co[,1:2] %>%
  as_tibble(rownames = "fac")

nutrient_bca_div_baryc <- nutrient_bca_div$li[,1:2]%>%
  as_tibble(rownames = "code")
nutrient_bca_div_fac <- nutrient_bca_div$co[,1:2] %>%
  as_tibble(rownames = "fac")

nutrient_pca_points <- nutrient_pca$li[,1:2]
nutrient_pca_fac <- nutrient_pca$co[,1:2]%>%
  as_tibble(rownames = "nutrient")

# Factorial plan BCA leg
BCA_leg  <- ggplot(nutrient_pca_points, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+  geom_point(alpha=0.2)+
  geom_point(data = nutrient_bca_leg_baryc, aes(x = Axis1, y = Axis2), col = "red") +
  geom_text(data = nutrient_bca_leg_baryc, aes(x = Axis1, y = Axis2, label = code), 
            col = "red", nudge_y = .1, size = 10/.pt, hjust=0,vjust=0) +
  annotate("text", x = -2, y = 4, , size = 18/.pt, fontface = "bold", fill = "white",
           label = paste("Between Legume \nclass inertia % = ", round(nutrient_bca_leg$ratio, 3)*100))+
  geom_segment(data = nutrient_bca_leg_fac, aes(x = 0, y = 0, xend = Comp1*10, yend = Comp2*10),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data = nutrient_bca_leg_fac, aes(x = Comp1*10, y = Comp2*10, label = fac), 
            col = "purple", nudge_y = .1, size = 14/.pt, hjust=0,vjust=0) +
  labs(x = paste("Axis 1 (", round(leg_eig[1], 2), "%)", sep =""),
       y = paste("Axis 2 (", round(leg_eig[2], 2), "%)", sep =""))+
  xlim(-8,8)+
  theme_bw()

# Factorial plan BCA div
BCA_div <- ggplot(nutrient_pca_points, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)+  geom_point(alpha=0.2)+
  geom_point(data = nutrient_bca_div_baryc, aes(x = Axis1, y = Axis2), col = "red") +
  geom_text(data = nutrient_bca_div_baryc, aes(x = Axis1, y = Axis2, label = code), 
            col = "red", nudge_y = .1, size = 10/.pt, hjust=0,vjust=0) +
  annotate("text", x = -2, y = 4, , size = 18/.pt, fontface = "bold", fill = "white",
           label = paste("Between earthworm diversity \nclass inertia % = ", round(nutrient_bca_div$ratio, 3)*100))+
  geom_segment(data = nutrient_bca_div_fac, aes(x = 0, y = 0, xend = Comp1*10, yend = Comp2*10),
               arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(data = nutrient_bca_div_fac, aes(x = Comp1*10, y = Comp2*10, label = fac), 
            col = "purple", nudge_y = .1, size = 14/.pt, hjust=0,vjust=0) +
  labs(x = paste("Axis 1 (", round(div_eig[1], 2), "%)", sep =""),
       y = paste("Axis 2 (", round(div_eig[2], 2), "%)", sep =""))+
  xlim(-8,8)+
  theme_bw()

nutrient_p <- ggarrange(BCA_leg, BCA_div, ncol = 2, labels = c("A", "B"))
nutrient_p
