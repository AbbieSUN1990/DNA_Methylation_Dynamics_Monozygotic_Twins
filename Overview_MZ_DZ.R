
library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)

##################Sample information
metadata <- read.csv("sampleInfo2.csv")
metadata <- metadata[order(metadata$Tissue, metadata$Age), ]
metadata <- metadata %>% 
  mutate(ID = row_number())

layer_gap <- 0.05
radius_params <- list(
  Tissue = list(inner = 0.1, outer = 0.3),    
  Twin_Type = list(inner = 0.35, outer = 0.55),    
  Sex = list(inner = 0.6, outer = 0.8), 
  Ethnicity = list(inner = 0.85, outer = 1.05),  
  Age = list(inner = 1.1, outer = 1.4),
  label = list(inner = 1.45, outer = 1.65)
)

max_age_radius <- radius_params$Age$inner + 
  (radius_params$Age$outer - radius_params$Age$inner) 

box_data <- metadata %>% 
  mutate(
    xmin = ID - 0.4,
    xmax = ID + 0.4,
    ymin = radius_params$Age$inner - 0.02,
    ymax = max_age_radius + 0.05
  )

blank_rings <- data.frame(
  ymin = c(radius_params$Tissue$outer, 
           radius_params$Twin_Type$outer, 
           radius_params$Sex$outer, 
           radius_params$Ethnicity$outer),
  ymax = c(radius_params$Twin_Type$inner, 
           radius_params$Sex$inner, 
           radius_params$Ethnicity$inner, 
           radius_params$Age$inner)
)

ggplot(metadata) +
  geom_rect(
    aes(xmin = ID - 0.5, xmax = ID + 0.5,
        ymin = radius_params$Tissue$inner, 
        ymax = radius_params$Tissue$outer, fill = Tissue),
    color = "white", linewidth = 0.3) +
  geom_rect(
    aes(xmin = ID - 0.5, xmax = ID + 0.5,
        ymin = radius_params$Twin_Type$inner, 
        ymax = radius_params$Twin_Type$outer, fill = Twin_Type),
    color = "white", linewidth = 0.3) +
  geom_rect(
    aes(xmin = ID - 0.5, xmax = ID + 0.5,
        ymin = radius_params$Sex$inner, 
        ymax = radius_params$Sex$outer, fill = Sex),
    color = "white", linewidth = 0.3) 
  geom_rect(
    aes(xmin = ID - 0.5, xmax = ID + 0.5,
        ymin = radius_params$Ethnicity$inner, 
        ymax = radius_params$Ethnicity$outer, fill = Ethnicity),
    color = "white", linewidth = 0.3) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    data = box_data, fill = NA, color = "black", linewidth = 0.3) + 
  geom_segment(
    aes(x = ID, y = radius_params$Age$inner,
        xend = ID, yend = radius_params$Age$inner + (radius_params$Age$outer - radius_params$Age$inner)*(Age/max(Age))),
    color = "#008585", linewidth = 0.6) +
  geom_point(
    aes(x = ID, y = radius_params$Age$inner + (radius_params$Age$outer - radius_params$Age$inner)*(Age/max(Age))),
    color = "#008585", size = 1.5) +
  geom_text(
    aes(x = ID, y = radius_params$label$inner,
        label = Sample_ID), 
    vjust = 0.5, hjust = 1, size = 3) +
geom_rect(
  aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
  data = blank_rings, fill = "white", inherit.aes = FALSE) +
coord_polar(start = 210 * pi/180, direction = -1, clip = "off") +
  scale_y_continuous(limits = c(0, radius_params$Age$outer * 1.1), 
                     expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none",
        legend.box = "vertical",
        legend.spacing.y = unit(0.2, "cm"),
        plot.margin = unit(rep(0.5, 4), "lines")
  ) +
  scale_fill_manual(
    values = c( 
      "blood" = "#CD5C5C", "saliva" = "#8973b8",
      "DZ" = "#f5e502", "MZ" = "#cdeba0",
      "M" = "#9ac7a8", "F" = "#ffc2bd",
      "Han" = "#f0c1ad", "HaNi" = "#fdffb6", "Yi" = "#b8e0d4"
      ))

