librarian::shelf(tidyverse, dplyr, ggpubr, Rmisc, ade4)


##################    17/03/2025  ######################
myenz <- read.csv("data/derived-data/enz.csv", h = T, sep = ",") 
anion <- read.csv("data/raw-data/anion.csv", h = T, sep = ";") 

nutrient0 <- left_join(myenz, anion)

biom <- read.csv("data/derived-data/pl_biom_Ind.csv", h = T, sep = ",") %>%
  select(RT, div, leg, strata, DW) %>% 
  group_by(RT, div, leg, strata) %>%
  dplyr::summarise(DWs= sum(DW)) %>%
  filter(!strata == 0) %>%
  mutate(volume = case_when(strata == "1" ~ 600,
                             strata == "2" ~ 1200,
                             strata == "3" ~ 1800,
                             strata == "4" ~ 2400)) %>%
  mutate(DWvol = DWs/volume) %>%
  select(-c(DWs, volume)) %>%
  dplyr::rename("strate" = "strata")

nutrient <- left_join(nutrient0, biom, by = join_by(RT, leg, div, strate))

# BCA 
nutrient_pca_df <- nutrient %>%
  filter(glu != 0) %>%
  mutate(ratioCN = log(glu)/log(ure),
         ratioCP = log(glu)/log(phos),
         ratioCS = log(glu)/log(ary)) %>%
  select(RT,leg, strate, div, DWvol, glu, ratioCP, ratioCN, ratioCS, SO42, NO3, HPO42, hum, NH4) 

nutrient_pca <- dudi.pca(nutrient_pca_df[, c(6:14)], scannf = FALSE, nf = 4)
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

nutrient_bca_points <- cbind(nutrient_bca_strate$ls, 
                             strate = as.factor(nutrient_pca_df$strate)) %>%
                       mutate(strcode = case_when(strate == "1" ~ "0-10 cm",
                                         strate == "2" ~ "10-30 cm",
                                         strate == "3" ~ "30-50 cm",
                                         strate == "4" ~ "50-100 cm"))
          
nutrient_bca_fac <- nutrient_bca_strate$co[,1:2]%>%
  as_tibble(rownames = "nutrient")

quanti.coord <- supcol(nutrient_pca, nutrient_pca_df$DWvol) %>%
  .$cosup


# Factorial plan BCA strate
BCA_strate  <- ggplot(nutrient_bca_points, aes(x = CS1, y = CS2))+
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)+  
  geom_point(alpha=0.2, aes(colour = strcode))+
  stat_ellipse(linetype = 2, aes(colour = strcode))+
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
  geom_segment(data = quanti.coord, aes(x = 0, y = 0, xend = -Comp1*10, yend = -Comp2*10),
               arrow = arrow(length = unit(0.2, "cm")), colour = "seagreen")+
  geom_text(data = quanti.coord, aes(x = -Comp1*10, y = -(Comp2+.01)*10, label = "Root mass density"), 
            col = "seagreen", nudge_y = .1, size = 12/.pt, hjust=1,vjust=1) +
  labs(x = paste("Axis 1 (", round(strate_eig[1], 2), "%)", sep =""),
       y = paste("Axis 2 (", round(strate_eig[2], 2), "%)", sep =""),
       colour = "Soil layer")+
  lims(x=c(-6,9), y=c(-5.5, 5.5))+
  theme_bw()
BCA_strate

