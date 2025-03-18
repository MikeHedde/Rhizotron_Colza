
require(ggplot2)


#Chargement jeu de donn?es 
stoe<-read.csv("data/raw-data/stoe.csv", h=T, sep=";")




############################################################################
#repr?sentation des rapports stoechiom?triques
p <- ggplot(stoe, aes(y=as.numeric(strata), x=ratio, col=div)) +
  geom_point( size = 5) +
  geom_path() +
  facet_grid(leg~stoe, scales = "free") +
  labs(y = "Rhizotron depth (cm)", x ="ecoenzymatic ratio") +
  scale_color_discrete (name = "Earthworm functional \ndiversity level", 
                      labels = c("div0", "div1", "div2", "div3"))
p


p <- ggplot(stoe, aes(y=as.numeric(strata), x=ratio, col=leg)) +
  geom_point( size = 5) +
  geom_path() +
  facet_grid(div~stoe, scales = "free") +
  labs(y = "Rhizotron depth (cm)", x ="ecoenzymatic ratio") +
  theme_bw()
p
