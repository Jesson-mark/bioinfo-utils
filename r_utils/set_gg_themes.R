# 用来修改ggplot2的默认主题

library(ggplot2)
classic_theme <- theme(
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), # 图片与四周边缘的距离
    plot.title = element_text(hjust=0.5, 
                                size = 20, 
                                face = "bold",
                                color = "#22292F", 
                                margin = margin(b = 2, t=2)),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text = element_text(size=15, color = "#22292F"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=15)
) 
bw_theme <- theme(
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), # 图片与四周边缘的距离
    plot.title = element_text(hjust=0.5, 
                                size = 20, 
                                face = "bold",
                                color = "#22292F", 
                                margin = margin(b = 2, t=2)),
    axis.title.x = element_text(size=15),
    axis.title.y = element_text(size=15),
    axis.text = element_text(size=15, color = "#22292F"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size=1.5), 
    legend.text = element_text(size=15)
) 

set_bw_theme <- function() {
  theme_set(theme_bw(base_size = 12) +
              bw_theme)
}

set_classic_theme <- function() {
  theme_set(theme_classic(base_size = 12) +
              classic_theme)
}

# margin: see https://r-charts.com/ggplot2/margins/
