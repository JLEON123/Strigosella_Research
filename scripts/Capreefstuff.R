# Download gg plot and ggmosiac in R and then run this code
library(ggplot2)
library(tidyverse)
library(ggmosaic)


df <- read_csv("./data/Trial1_data")
attach(df)
Drought_survive <- chisq.test(Alive, Drought)
Temp_survive <- chisq.test(Alive, Heat)


heat_mosiac <- 
  df %>%
  ggplot() +
        geom_mosaic(aes(x = product(Alive, Heat)),
                    fill=c("#000000", "#275d38", "#000000", "#275d38")) +
        labs(x = "Drought", y = "Alive") +
        theme_minimal() +
        theme(
                axis.title = element_text(face = "bold", size = 20),
                axis.text = element_text(face = "bold.italic", size = 15),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank() 
        )

Drought_mosiac <- 
  df %>%
  ggplot() +
  geom_mosaic(aes(x = product(Alive, Heat)),
              fill=c("#000000", "#275d38", "#000000", "#275d38")) +
  labs(x = "Drought", y = "Alive") +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold.italic", size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank() 
  )

ggsave(filename = "./figures/Drought_mosiac.png", plot = Drought_mosiac, dpi = 300)
ggsave(filename = "./figures/Heat_mosiac.png", plot = heat_mosiac, dpi = 300)
citation()
