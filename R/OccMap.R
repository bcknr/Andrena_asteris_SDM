if(!require(pacman)) install.packages("pacman")
library(pacman)

pacman::p_load(tidyverse, lubridate, raster, sf, ggmap, 
               tmap, rnaturalearth, grid)

pacman::p_load_gh("ropensci/rnaturalearthhires", "ropensci/rnaturalearthdata")

setwd("./R")

gbif <- read_tsv("../occ/Aa_GBIF.txt")
scan <- read_csv("../occ/Aa_SCAN.csv")
amnh <- read_tsv("../occ/Aa_AMNH.txt")

t_scan <- scan %>% 
  dplyr::select(id, date = eventDate, lat = decimalLatitude, lon = decimalLongitude) %>% 
  drop_na() %>% 
  mutate(lat = round(lat,5), lon = round(lon,5)) %>%
  mutate(source = "SCAN")


t_gbif <- gbif %>% 
  dplyr::select(gbifID, eventDate, verbatimEventDate, lat = decimalLatitude, lon = decimalLongitude) %>% 
  separate(eventDate, c("date", NA), sep = " ") %>% 
  mutate(verbatimEventDate = ifelse(is.na(date), verbatimEventDate, NA)) %>% 
  unite(date, date, verbatimEventDate, na.rm = TRUE) %>% 
  mutate(date = parse_date_time(date, orders = c("mdy", "ymd", "dmy"))) %>% 
  drop_na() %>% 
  mutate(lat = round(lat,5), lon = round(lon,5)) %>%
  mutate(source = "GBIF")

t_amnh <- amnh %>% 
  mutate(date = parse_date_time(amnh$Start_Date, orders = c("mdy", "dmy"))) %>%
  filter(year(date) <= Det_Date) %>%    
  dplyr::select(PBIUSI, date, lat = Lat, lon = Lon) %>% 
  drop_na() %>% 
  mutate(lat = round(lat,5), lon = round(lon,5)) %>%
  mutate(source = "AMNH")

#Parsing errors for year only (old data) or "no date provided"

occ.joined <- full_join(t_scan, t_gbif, by = c("date", "lat", "lon")) %>% 
  full_join(t_amnh, by = c("date", "lat", "lon")) %>%
  filter(!duplicated(.[c("date", "lat", "lon")])) %>% 
  unite(data_source, source.x, source.y, source, sep = "/", na.rm = TRUE) %>% 
  dplyr::select(date, lat, lon, data_source, scanID = id, gbifID, PBIUSI)

occ <- occ.joined %>% 
  mutate(lat = floor(lat*10000)/10000, lon = floor(lon*10000)/10000) %>% 
  distinct(lat, lon, .keep_all = TRUE) %>% 
  mutate("spp" = "Aasteris")

qmplot(x = lon, y = lat, data = occ, maptype = "toner-lite", mapcolor = "bw", source = "stamen", force = T)





pal <- c("#ffffff","#dad7cd","#a3b18a","#588157","#3a5a40","#344e41")

ref <- ne_countries(continent = c("North America", "South America"), scale = 50, returnclass = "sf") %>% 
  st_crop(xmin = -170, ymin = 0, xmax = -11, ymax = 83)

adm <- ne_states(country = c("United States of America", "Canada", "Mexico"), returnclass = "sf") %>% 
  st_crop(xmin = -102, ymin = 24, xmax = -56, ymax = 52) 


occs <- st_as_sf(occ, coords = c("lon", "lat"), crs = st_crs(4326))

bb <- extent(-100.7498, -57.85033, 26.00973, 50.74983) %>% 
  st_bbox()

cortland <- matrix(c(42.5563, -76.2282), nrow =1) %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c(2, 1), crs = st_crs(4326))

refmap <- tm_shape(ref) +
  tm_polygons(alpha = 0.3, col = "#405F63", border.alpha = 0) +
  tm_shape(st_as_sfc(bb)) +
  tm_polygons(alpha = 0, border.col = "#EC4E20", lwd = 3) +
  tm_layout(inner.margins = 0, bg.color = NA, frame = FALSE)

occmap <- tm_shape(adm, bbox = bb) +
  tm_polygons(border.col = "#405F63", border.alpha = 0.2, col = "#FFFFFF") +
  tm_shape(occs) +
  tm_dots(shape = 1, size = 0.08, alpha = 0.85) +
  tm_shape(cortland) +
  tm_dots(shape = 19, size = 0.1, col = "#EC4E20") +
  tm_compass(type = "arrow", show.labels = 0, color.light = "#000000", position = c(0.88, 0.09)) +
  tm_scale_bar(breaks = c(0, 200), color.light = "#000000", text.size = 0.7,
               position = c(0.872, 0.01)) +
  tm_layout(inner.margins = 0, bg.color = "#E9EBEC", frame = "#D3D7D9",
            legend.outside = T, legend.outside.position = "right")


occmap
print(refmap, vp = grid::viewport(0.54, 0.243, width = 0.25, height = 0.3))

ext <- extent(bb)
asp <- (ext[4] - ext[3])/(abs(ext[1]) - abs(ext[2]))
vp <- grid::viewport(0.719, 0.165, width = 0.25, height = 0.3)

tmap_save(occmap, filename = "../pred/dist/Aa_occ_figure.png",
          dpi = 600, insets_tm = refmap, insets_vp = vp,
          height = asp*200, width = 200, units = "mm")
