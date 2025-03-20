librarian::shelf(tidyverse, dplyr, ggpubr, Rmisc, ade4)


##################    17/03/2025  ######################
myenz <- read.csv("data/derived-data/enz.csv", h = T, sep = ",") 
anion <- read.csv("data/raw-data/anion.csv", h = T, sep = ";") 

nutrient <- left_join(myenz, anion)


# BCA 
nutrient_pca_df <- nutrient %>%
  filter(glu != 0) %>%
  mutate(ratioCN = log(glu)/log(ure),
         ratioCP = log(glu)/log(phos),
         ratioCS = log(glu)/log(ary)) %>%
  select(leg, strate, div, glu, ratioCP, ratioCN, ratioCS, SO42, NO3, HPO42, hum, NH4) 

nutrient_pca <- dudi.pca(nutrient_pca_df[, -c(1:3)], scannf = FALSE, nf = 4)
nutrient_bca_leg <- bca(nutrient_pca, as.factor(nutrient_pca_df$leg), scan = FALSE, nf = 2)
nutrient_bca_div <- bca(nutrient_pca, as.factor(nutrient_pca_df$div), scan = FALSE, nf = 2)
nutrient_bca_strate <- bca(nutrient_pca, as.factor(nutrient_pca_df$strate), scan = FALSE, nf = 2)

nutrient_bca_strate$ratio
nutrient_bca_leg$ratio
nutrient_bca_div$ratio

strate_eig <- nutrient_bca_strate$eig/sum(nutrient_bca_strate$eig)*100
nutrient_bca_strate_baryc <- nutrient_bca_strate$li[,1:2]%>%
  as_tibble(rownames = "strate")
nutrient_bca_strate_fac <- nutrient_bca_strate$co[,1:2] %>%
  as_tibble(rownames = "fac") %>%
  mutate(expr = case_when(fac == "glu" ~ "Glucosidase",
                          fac == "ratioCP" ~ "Glu:Phos ratio",
                          fac == "ratioCN" ~ "Glu:Ure ratio",
                          fac == "ratioCS" ~ "Glu:AryS ratio",
                          fac == "SO42" ~ "Sulfate",
                          fac == "NO3" ~ "Nitrate",
                          fac == "HPO42" ~ "Orthophosphate",
                          fac == "hum" ~ "Humidity",
                          fac == "NH4" ~ "Ammonium"))

nutrient_bca_points <- cbind(nutrient_bca_strate$ls, strate = as.factor(nutrient_pca_df$strate))
nutrient_bca_fac <- nutrient_bca_strate$co[,1:2]%>%
  as_tibble(rownames = "nutrient")

# Factorial plan BCA strate
BCA_strate  <- ggplot(nutrient_pca_points, aes(x = -Axis1, y = -Axis2))+
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)+  
  geom_point(alpha=0.2, aes(colour = strate))+
  stat_ellipse(aes(linetype = strate, colour = strate))+
  #geom_point(data = nutrient_bca_strate_baryc, size = 3, aes(x = Axis1, y = Axis2, colour = strate)) +
  #geom_text(data = nutrient_bca_strate_baryc, aes(x = Axis1, y = Axis2, label = strate), 
  #          nudge_y = .1, size = 13/.pt, hjust=0,vjust=0) +
  annotate("text", x = -5.5, y = -4, , size = 14/.pt, fontface = "bold", hjust = 0,
           label = paste("Between class inertia \nSoil strata = ", round(nutrient_bca_strate$ratio, 3)*100, "%\n",
                         "Legume presence = ", round(nutrient_bca_leg$ratio, 3)*100, "%\n",
                         "Earthworm diversity = ", round(nutrient_bca_div$ratio, 3)*100, "%\n",
                   sep =""))+
  geom_segment(data = nutrient_bca_strate_fac, aes(x = 0, y = 0, xend = Comp1*10, yend = Comp2*10),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(data = nutrient_bca_strate_fac, aes(x = Comp1*10, y = Comp2*10, label = expr), 
            col = "gray30", nudge_y = .1, size = 12/.pt, hjust=0,vjust=0) +
  labs(x = paste("Axis 1 (", round(strate_eig[1], 2), "%)", sep =""),
       y = paste("Axis 2 (", round(strate_eig[2], 2), "%)", sep =""))+
  xlim(-6,9)+
  theme_bw()
BCA_strate

