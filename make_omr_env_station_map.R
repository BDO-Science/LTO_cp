library(tidyverse)
library(lubridate) #today()
library(readr)
library(here)
library(viridis)
library(deltamapr)
library(sf)
library(ggspatial)
library(readxl)
library(cowplot)


stations_env <- read_excel(here("data_raw/omr_station_locations.xlsx"))
stations_sf <- st_as_sf(stations_env, coords = c("Longitude", "Latitude"), crs = 4326, remove = FALSE)
delta_4326 <- st_transform(WW_Delta, crs = st_crs(stations_sf))

(omr_map <- ggplot() +
    geom_sf(data = delta_4326, fill = "gray90", color = "gray80") +
    geom_sf(data = stations_sf, aes(color = `Data type`, shape = `Data type`), size =  3, inherit.aes = FALSE)+
    geom_text(data = stations_sf, aes(x = Longitude, y = Latitude, label = Station), size = 3, vjust = 1, hjust = 1.2)+
    scale_color_viridis(discrete=TRUE) +
    scale_shape_manual(values = c(18, 8, 15, 19))+
    scale_x_continuous(limits = c(-121.9, -121.3)) +
    scale_y_continuous( limits = c(37.8, 38.5))+
    # labs(fill = "Region") +
    annotation_north_arrow(location = "tl", which_north = "true",
                           pad_x = unit(.1, "in"), pad_y = unit(0.2, "in"),
                           style = north_arrow_fancy_orienteering) +
    annotation_scale(location = "bl", bar_cols = c("black", "white", "black", "white")) +
    guides(color = guide_legend(nrow = 2))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
          legend.position = "top",
          legend.box = "vertical",
          legend.title = element_blank()))


png(filename = here("figures", "map_omr.png"), width = 5, height = 6, units = "in", pointsize = 12, family = "sans", res = 300)
omr_map
dev.off()
