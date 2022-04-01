library(tidyverse)

nest <- read_csv("./data/NestDepth.csv")
nest_norange <- nest %>% 
  filter(min_depth == max_depth)

head(nest)



plot <- ggplot(data = nest, aes(index, median_depth, color = group)) +
#  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = "#A37B73", alpha = 0.15, color = NA) +
#  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#EBEBFF", alpha = 0.3, color = NA) +
  geom_hline(yintercept = 0) +
  geom_linerange(aes(ymin = min_depth, ymax = max_depth), lwd = 5) +
  geom_point(data = nest_norange, size = 3)+
  scale_color_manual(values = c("#87B38D", "#FFB01F", "#336699"), name = "Larkin Group") +
  labs(x = "Species", y = "Nest Depth (cm)") +
  scale_x_continuous(n.breaks = nrow(nest), limits = c(1,41), expand = c(0.014,0.014), minor_breaks = seq(1,41,1)) +
  theme_bw()+
  theme(text = element_text(size = 16), legend.position = c(0.9,0.1), legend.background = element_rect(fill = NA)) +
  scale_y_reverse(breaks = seq(0,300,25),  minor_breaks = seq(0,300,25))
  
plot

ggsave("./NestDepthFig.jpg", plot, dpi = 600)
