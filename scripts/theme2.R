ggplot2::theme_set(theme_bw() + 
                     theme(plot.title = element_text(face = "bold", color = "#000000", size = 35), 
                           panel.border = element_blank(),
                           panel.grid.major.y = element_blank(),
                           axis.text = element_text(face = 'bold', size = 20, color = "#000000"),
                           axis.title = element_text(face = "bold", size = 25), 
                           axis.text.x = element_text(angle = .5, hjust = .5),
                           axis.text.y = element_text(angle = .1, hjust = .5, vjust = 0),
                           axis.ticks = element_blank(),
                           legend.text = element_text(face = "bold", size = 20),
                           legend.title = element_text(face = "bold", size = 25),
                           legend.title.align = 0.5))
pal <- 
  c("#6A0213", # 1
    "#359B73", # 2
    "#FF7F52", # 3
    "#FFB2FD", # 4
    "#064273", # 5
    "#D55E00", # 6
    "#9F0162", # 7
    "#FFD3CD", # 8
    "#F0E442", # 9
    "#00AF8E", # 10
    "#004002", # 11
    "#DC3220", # 12
    "#994F00", # 13
    "#000000", # 14
    "Grey50", # 15
    "Grey100" #16
  )

pal.abundance<- 
  c("#359B73", # 1
    "#6A0213", # 2
    "#FF7F52", # 3
    "#FFB2FD", # 4
    "#064273", # 5
    "#D55E00", # 6
    "#9F0162", # 7
    "#FFD3CD", # 8
    "#F0E442", # 9
    "#00AF8E", # 10
    "#004002", # 11
    "#DC3220", # 12
    "#000000", # 13
    "grey50",  # 14
    "#994F00", # 15
    "Yellow"   #16
  )
