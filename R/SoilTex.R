#Soil texture triangle at surface and depth of brood cells (based on https://saryace.github.io/flipbook_soiltexture_en/#25)

library(tidyverse)
library(ggtern)

data(USDA)

USDA_text <- USDA  %>% group_by(Label) %>%
  summarise_if(is.numeric, mean, na.rm = TRUE) %>% 
  cbind(abbv = c("Cl", "SaCl", "SaClLo", "SaLo", "LoSa", "Sa", "ClLo", "Lo", "SiLo", "SiCl", "SiClLo", "Si"))

soil <- read_csv("./data/SoilTex.csv")

tri <- ggplot(USDA, aes(x = Sand, y = Clay, z = Silt)) +
  coord_tern(L = "x", T = "y", R = "z") +
  geom_polygon(aes(fill = Label),alpha = 0.0,size = 0.5,color = "black") +
  geom_text(data = USDA_text, aes(label = abbv), color = "black", size = 3) +
  geom_point(data = soil, aes(x = sand, y = clay, z = silt, color = depth), size = 1.75, alpha = 0.8) +
  labs(yarrow = "Clay (%)",
       zarrow = "Silt (%)",
       xarrow = "Sand (%)") +
  theme_bw() +
  theme_showarrows()+
  theme_clockwise() +
  theme(legend.position = "none")

ggsave("./SoilTexFig.jpg", tri, dpi = 600, height = 10, units = "cm")

