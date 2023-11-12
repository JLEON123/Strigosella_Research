library(tidyverse)

pal.discrete2 <- 
  c("#359B73", # 1
    "#F0E442", # 2
    "#FFB2FD", # 3
    "#FF7F52", # 4
    "#064273", # 5
    "#000000", # 6
    "#D55E00", # 7
    "#9F0162", # 8
    "#FFD3CD", # 9
    "#6A0213", # 10
    "#00AF8E", # 11
    "#004002"
  )


ggplot2::theme_set(theme_minimal() +
                     theme(axis.text = element_text(face='bold'),
                           legend.title = element_text(face='bold')))





