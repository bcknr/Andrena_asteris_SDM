library(tidyverse)

nest <- read_csv("./data/NestDepth.csv") %>%
  mutate(spp = str_replace(spp, "Andrena", "A.")) %>% 
  mutate(group = ifelse(group == "A", "Larkin Gr. A (=Callandrena s.s.)", group)) %>% 
  mutate(group = ifelse(group == "B", "Larkin Gr. BCD", group))

nest_norange <- nest %>% 
  filter(min_depth == max_depth) 

head(nest)

focal <- nest %>% 
  filter(spp == "A. asteris")
#536F72 #A3B18A

plot <- ggplot(data = nest, aes(color = group)) +
  geom_hline(yintercept = 0) +
  geom_linerange(aes(x = fct_reorder(spp, index), ymin = min_depth, ymax = max_depth), lwd = 5) +
  geom_linerange(data = focal, aes(x = fct_reorder(spp, index), ymin = min_depth, ymax = max_depth), lwd = 5, color = "#EC4E20") +
  geom_point(data = nest_norange, aes(fct_reorder(spp, index), median_depth), size = 3) +
  scale_color_manual(values = c("#101C42","#E6874C", "#8F8F8F"), name = "") +
  labs(x = "Species", y = "Nest Depth (cm)") +
  scale_x_discrete(expand = c(0.06,0.06)) +
  theme_bw() +
  theme(text = element_text(size = 16), legend.position = c(0.85,0.1), legend.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")) +
  scale_y_reverse(breaks = seq(0,300,25),  minor_breaks = seq(0,300,25))
  
plot

ggsave("./NestDepthFig.jpg", plot, dpi = 600, height = 20, width = 30, units = "cm")


  



