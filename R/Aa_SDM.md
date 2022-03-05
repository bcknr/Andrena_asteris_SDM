Andrena asteris
================
Mark Buckner
2022-03-05

## Load Pkgs and functions

Project directory structure:

Andrena_asteris_SDM -

 occ - occurrence datasets  
 pred - .txt/s containing URLs to environmental data  
 R - code

### Pkgs

``` r
if(!require(pacman)) install.packages("pacman")
library(pacman)

pacman::p_load(tidyverse, lubridate, ggmap, spThin, 
               CoordinateCleaner, curl, raster, sf, 
               RStoolbox, rgdal, ENMeval, rJava, rasterVis, tmap, 
               rnaturalearth, grid)

pacman::p_load_gh("SEEG-Oxford/seegSDM", "ropensci/rnaturalearthhires", "ropensci/rnaturalearthdata")
```

### Functions

``` r
#envLoad(): Downloads and processes environmental raster data from source. 
#--paths = urls for data sources read in from .txt file;
#--filedir = path to dir where the processed raw files are or will be stored.
#--ref = anticipated path of a reference .tif with the correct resolution and projection (eg. worldclim layer) in the paths .txt ["../pred/env/NAME.tif"]
#--bb = bounding box for study area (any object an extent can be extracted from)

envLoad <- function(paths, filedir, ref, bb, type = "cont") {
  
  savepaths <- sapply(paths, function (x) gsub("\\?.*", "", x))
  dest <- sapply(savepaths, function(x) paste(filedir,"/", gsub(".*/(\\w+)", "\\1", x), sep = "")) %>% 
    sort()
  names <- sapply(savepaths, function(x) paste( gsub(".*/(\\w+)", "\\1", x), sep = "")) %>% 
    sort() %>% 
    sapply(function(x) paste("../pred/env/",x, sep = ""))
  
  names <- str_replace(names, ".zip$", ".tif")
  
  for (i in 1:length(paths)) {
    if(file.exists(dest[i]) | file.exists(names[i])){
      if (file.exists(names[i])) {
        print(paste("File exists:",names[i]))
      } else {
        unproc_cont <- TRUE
        print(paste("File exists but is unprocessed:",dest[i]))
      }
      
    } else {
      
      print(paste("Downloading: ", paths[i]))
      
      unproc_cont <- TRUE
      
      file <- try(curl::curl_download(url = paths[i], destfile = dest[i]))
      
      if(class(file) == "try-error") {
        download.file(url = paths[i], destfile = dest[i], method = "curl")
      }
      
      if(grepl(".*zip", dest[i])) {
        unzip(dest[i], exdir = filedir)
        
      }
    }
  }
  
  if(!all(file.exists(names))) {
    print("Processing files")
    env.files <- list.files(filedir, recursive = TRUE, pattern = ".*tif$") %>% 
      sapply(function(x) paste(filedir, "/", x, sep = ""))
    
    r <- raster(ref) %>% 
      raster::crop(extent(bb)) %>% 
      projectRaster(crs = "+init=epsg:4326")
    
    method <- ifelse(type == "cont", "bilinear", "ngb")
    print(paste("Method:", method))
    
    env.stack <- lapply(env.files, raster) %>% 
      lapply(projectRaster, to = r, method = "ngb") %>% 
      stack()
    
    print("Saving stack")
    
    bylayer <- ifelse(nlayers(env.stack) > 1, TRUE, FALSE)
    writeRaster(env.stack, names, format = "GTiff", bylayer = bylayer, overwrite = TRUE)
    
  } else {
    print("Loading processed files")
    env.stack <- stack(names)
  }
  return(env.stack)
  
}
```

## Load and process occ Data

### Load and filter by date

``` r
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

date_start <- as.Date("1971-01-01")
date_end <- as.Date("2010-12-31")

occ <- occ.joined %>% 
  filter(date >= date_start & date <= date_end) %>%
  mutate(lat = floor(lat*10000)/10000, lon = floor(lon*10000)/10000) %>% 
  distinct(lat, lon, .keep_all = TRUE) %>% 
  mutate("spp" = "Aasteris")

qmplot(x = lon, y = lat, data = occ, maptype = "toner-lite", mapcolor = "bw", source = "stamen", force = T)
```

![](Aa_SDM_files/figure-gfm/Occ-1.png)<!-- -->

### Check for georeferencing errors

``` r
flags <- clean_coordinates(x = occ, lon = "lon", 
                           lat = "lat", species = "spp",
                           tests = c("capitals", "centroids", 
                                     "equal", "gbif", "institutions", 
                                     "outliers", "seas", "zeros"))
```

    ## Testing coordinate validity

    ## Flagged 0 records.

    ## Testing equal lat/lon

    ## Flagged 0 records.

    ## Testing zero coordinates

    ## Flagged 0 records.

    ## Testing country capitals

    ## Flagged 0 records.

    ## Testing country centroids

    ## Flagged 0 records.

    ## Testing sea coordinates

    ## OGR data source with driver: ESRI Shapefile 
    ## Source: "C:\Users\mabuc\AppData\Local\Temp\RtmpwtqNJT", layer: "ne_50m_land"
    ## with 1420 features
    ## It has 3 fields
    ## Integer64 fields read as strings:  scalerank

    ## Flagged 22 records.

    ## Testing geographic outliers

    ## Flagged 0 records.

    ## Testing GBIF headquarters, flagging records around Copenhagen

    ## Flagged 0 records.

    ## Testing biodiversity institutions

    ## Flagged 2 records.

    ## Flagged 24 of 86 records, EQ = 0.28.

``` r
summary(flags)
```

    ##     .val     .equ     .zer     .cap     .cen     .sea     .otl     .gbf 
    ##        0        0        0        0        0       22        0        0 
    ##    .inst .summary 
    ##        2       24

``` r
plot(flags, lon = "lon", lat = "lat")
```

![](Aa_SDM_files/figure-gfm/Clean%20coordinates-1.png)<!-- -->

``` r
occ.flagged <- occ[!flags$.summary,]

write_csv(occ.flagged, file = "../occ/Aa_flagged.csv")
```

22 occurrences were flagged for potentially being in the sea. These
points are located on various Island which matches the metadata. The
georeferencing is sufficiently accurate for this analysis.

Two points match research institutions. The [occurrence
remarks](https://bugguide.net/node/view/370287) for 34630476 - Cornell
University provide no additional details this observation was excluded.
The metadata for 56951087 - Connecticut Agricultural Exp. list the
address of the station.

``` r
occs <- occ %>% 
  filter(is.na(scanID) | scanID != 56951087 & scanID != 34630476) %>% 
  write_csv(file = "../occ/Aa_combined.csv") %>% 
  dplyr::select(spp, lat, lon)
```

### Spatial Thinning w/ spThin

``` r
thinned <-
  thin( loc.data = occs, 
        lat.col = "lat", long.col = "lon", 
        spec.col = "spp", 
        thin.par = 10, reps = 100, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 5, 
        out.dir = "../occ/Aa_thinned_full/", out.base = "Aa_thinned", 
        write.log.file = TRUE,
        log.file = "../occ/Aa_thinned_full_log_file.txt" )
```

    ## ********************************************** 
    ##  Beginning Spatial Thinning.
    ##  Script Started at: Sat Mar 05 07:24:00 2022
    ## lat.long.thin.count
    ##  49 
    ## 100 
    ## [1] "Maximum number of records after thinning: 49"
    ## [1] "Number of data.frames with max records: 100"
    ## [1] "Writing new *.csv files"

    ## Warning in dir.create(out.dir, recursive = TRUE): '..\occ\Aa_thinned_full'
    ## already exists

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col =
    ## "spp", : Created new output directory: ../occ/Aa_thinned_full/

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col = "spp", : ../occ/Aa_thinned_full/Aa_thinned_thin1_new.csv ' already exists. Renaming file 
    ##                                    to avoid overwriting.

    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin1_new.csv"

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col = "spp", : ../occ/Aa_thinned_full/Aa_thinned_thin2_new.csv ' already exists. Renaming file 
    ##                                    to avoid overwriting.

    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin2_new.csv"

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col = "spp", : ../occ/Aa_thinned_full/Aa_thinned_thin3_new.csv ' already exists. Renaming file 
    ##                                    to avoid overwriting.

    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin3_new.csv"

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col = "spp", : ../occ/Aa_thinned_full/Aa_thinned_thin4_new.csv ' already exists. Renaming file 
    ##                                    to avoid overwriting.

    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin4_new.csv"

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col = "spp", : ../occ/Aa_thinned_full/Aa_thinned_thin5_new.csv ' already exists. Renaming file 
    ##                                    to avoid overwriting.

    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin5_new.csv"

``` r
plotThin(thinned)
```

![](Aa_SDM_files/figure-gfm/Thin%20occ-1.png)<!-- -->![](Aa_SDM_files/figure-gfm/Thin%20occ-2.png)<!-- -->![](Aa_SDM_files/figure-gfm/Thin%20occ-3.png)<!-- -->

#### Load thinned occ data

``` r
occt <- read_csv("../occ/Aa_thinned_full/Aa_thinned_thin1.csv") %>% 
  dplyr::select(lon, lat)


qmplot(x = lon, y = lat, data = occt, maptype = "toner-lite", mapcolor = "bw", source = "stamen", force = T)
```

![](Aa_SDM_files/figure-gfm/load%20thin-1.png)<!-- -->

Retained 49 after spatial thinning.

## Environmental Covariates

Bioclimatic Variables - [CHELSA Climate
v2.1](https://chelsa-climate.org/)

-   Bio 1 : Mean Annual Temperature
-   Bio 2 : Annual Mean Diurnal Range
-   Bio 3 : Isothermality
-   Bio 4 : Temperature Seasonality
-   Bio 5 : Max Temperature of Warmest Month
-   Bio 7 : Annual Temperature Range
-   Bio 8 : Mean Temperature of Wettest Quarter
-   Bio 9 : Mean Temperature of Driest Quarter
-   Bio 10 : Mean Temperature of Warmest Quarter
-   Bio 11 : Mean Temperature of Coldest Quarter
-   Bio 12 : Annual Precipitation
-   Bio 13 : Precipitation of Wettest Month
-   Bio 14 : Precipitation of Driest Month
-   Bio 15 : Precipitation Seasonality
-   Bio 16 : Precipitation of Wettest Quarter
-   Bio 17 : Precipitation of Driest Quarter
-   Bio 18 : Precipitation of Warmest Quarter
-   Bio 19 : Precipitation of Coldest Quarter

Other Climatologies - [CHELSA Climate v2.1](https://chelsa-climate.org/)

-   gddlgd0 : Last growing degree day above 0C

Topography - [WorldClim SRTM 5
arc-min](https://www.worldclim.org/data/worldclim21.html)

-   DEM

Other

-   [Soil sand at 0cm](https://zenodo.org/record/2525662#.YiISa3rMJOQ)
-   [Soil clay at 0cm](https://zenodo.org/record/2525663#.YiISZnrMJOQ)
-   [Land
    use](http://www.cec.org/north-american-environmental-atlas/land-cover-30m-2015-landsat-and-rapideye/)
    (excluded)

### Define study area

``` r
proj.res <- c(0.08333333,0.08333333)
proj.crs <- st_crs(4326)

occs.sf <- sf::st_as_sf(occt, coords = c("lon", "lat"), crs = proj.crs) %>% 
  st_cast("MULTIPOINT") %>% 
  st_union()

sa.bb <- st_bbox(occs.sf)
sa.bb[1:2] <- sa.bb[1:2] - 7
sa.bb[3:4] <- sa.bb[3:4] + 7
```

### Download and process env

``` r
dir.create("../pred/cat")
dir.create("../pred/cont")
dir.create("../pred/env")
```

    ## Warning in dir.create("../pred/env"): '..\pred\env' already exists

``` r
# Download and process continuous data
paths.cont <- read_delim("../pred/cov_paths_cont.txt", delim = "\\n", col_names = FALSE) %>% 
  pull()
```

    ## Rows: 23 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\\n"
    ## chr (1): X1
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
env.cont <- envLoad(paths.cont, filedir = "../pred/cont", ref = "../pred/cont/wc2.1_5m_elev.tif", bb = sa.bb)
```

    ## [1] "File exists: ../pred/env/CHELSA_bio1_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio10_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio11_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio12_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio13_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio14_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio15_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio16_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio17_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio18_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio19_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio2_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio3_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio4_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio5_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio6_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio7_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio8_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_bio9_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/CHELSA_gddlgd0_1981-2010_V.2.1.tif "
    ## [1] "File exists: ../pred/env/sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif"
    ## [1] "File exists: ../pred/env/sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2.tif"
    ## [1] "File exists: ../pred/env/wc2.1_5m_elev.tif"
    ## [1] "Loading processed files"

``` r
# Download and process categorical data
paths.cat <- read_delim("../pred/cov_paths_cat.txt", delim = "\\n", col_names = FALSE) %>% 
  pull()
```

    ## Rows: 1 Columns: 1
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: "\\n"
    ## chr (1): X1
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
env.cat <- envLoad(paths.cat, filedir = "../pred/cat", ref = "../pred/cont/wc2.1_5m_elev.tif", bb = sa.bb, type = "cat")%>% 
  `names<-`("LU")
```

    ## [1] "File exists: ../pred/env/north_america_2015_v2.tif"
    ## [1] "Loading processed files"

``` r
unlink("../pred/cat", recursive = TRUE)
unlink("../pred/cont", recursive = TRUE)
```

### PCA

``` r
#GDD data modification to allow proper masking
env.cont$CHELSA_gddlgd0_1981.2010_V.2.1[is.na(env.cont$CHELSA_gddlgd0_1981.2010_V.2.1)] <- 0

mmask <- masterMask(env.cont)

env.m <- mask(env.cont, mmask)
env.cat.m <- mask(env.cat, mmask)

dir.create("../pred/PCA")
pca <- rasterPCA(env.m, spca = TRUE, maskCheck = FALSE, filename = "../pred/PCA/envPCA.tif", overwrite = TRUE)

summary(pca$model)
```

    ## Importance of components:
    ##                          Comp.1    Comp.2     Comp.3     Comp.4     Comp.5
    ## Standard deviation     3.273579 2.4283101 1.40026590 1.15590054 0.98093012
    ## Proportion of Variance 0.465927 0.2563778 0.08524977 0.05809157 0.04183582
    ## Cumulative Proportion  0.465927 0.7223048 0.80755458 0.86564614 0.90748197
    ##                            Comp.6     Comp.7    Comp.8     Comp.9     Comp.10
    ## Standard deviation     0.75240800 0.71645352 0.6015950 0.55287585 0.335507138
    ## Proportion of Variance 0.02461382 0.02231764 0.0157355 0.01329007 0.004894132
    ## Cumulative Proportion  0.93209578 0.95441342 0.9701489 0.98343899 0.988333127
    ##                            Comp.11     Comp.12     Comp.13     Comp.14
    ## Standard deviation     0.302936735 0.258285028 0.212672731 0.154748337
    ## Proportion of Variance 0.003990029 0.002900485 0.001966508 0.001041176
    ## Cumulative Proportion  0.992323156 0.995223641 0.997190149 0.998231325
    ##                             Comp.15      Comp.16      Comp.17      Comp.18
    ## Standard deviation     0.1174087452 0.0954031748 0.0855662141 0.0713953354
    ## Proportion of Variance 0.0005993397 0.0003957289 0.0003183294 0.0002216215
    ## Cumulative Proportion  0.9988306646 0.9992263936 0.9995447230 0.9997663445
    ##                             Comp.19      Comp.20      Comp.21      Comp.22
    ## Standard deviation     0.0572040760 4.254521e-02 1.466703e-02 8.749498e-03
    ## Proportion of Variance 0.0001422742 7.869977e-05 9.353122e-06 3.328422e-06
    ## Cumulative Proportion  0.9999086187 9.999873e-01 9.999967e-01 1.000000e+00
    ##                             Comp.23
    ## Standard deviation     1.046154e-07
    ## Proportion of Variance 4.758424e-16
    ## Cumulative Proportion  1.000000e+00

``` r
loadings(pca$model)
```

    ## 
    ## Loadings:
    ##                                                               Comp.1 Comp.2
    ## CHELSA_bio1_1981.2010_V.2.1                                    0.271  0.184
    ## CHELSA_bio10_1981.2010_V.2.1                                   0.237  0.246
    ## CHELSA_bio11_1981.2010_V.2.1                                   0.282  0.145
    ## CHELSA_bio12_1981.2010_V.2.1                                   0.210 -0.282
    ## CHELSA_bio13_1981.2010_V.2.1                                   0.193 -0.185
    ## CHELSA_bio14_1981.2010_V.2.1                                   0.193 -0.282
    ## CHELSA_bio15_1981.2010_V.2.1                                  -0.140  0.248
    ## CHELSA_bio16_1981.2010_V.2.1                                   0.186 -0.219
    ## CHELSA_bio17_1981.2010_V.2.1                                   0.203 -0.276
    ## CHELSA_bio18_1981.2010_V.2.1                                   0.125 -0.214
    ## CHELSA_bio19_1981.2010_V.2.1                                   0.211 -0.258
    ## CHELSA_bio2_1981.2010_V.2.1                                           0.292
    ## CHELSA_bio3_1981.2010_V.2.1                                    0.245  0.194
    ## CHELSA_bio4_1981.2010_V.2.1                                   -0.290       
    ## CHELSA_bio5_1981.2010_V.2.1                                    0.208  0.285
    ## CHELSA_bio6_1981.2010_V.2.1                                    0.280  0.148
    ## CHELSA_bio7_1981.2010_V.2.1                                   -0.286       
    ## CHELSA_bio8_1981.2010_V.2.1                                           0.256
    ## CHELSA_bio9_1981.2010_V.2.1                                    0.280       
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                -0.246       
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2         0.251
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2        -0.136
    ## wc2.1_5m_elev                                                 -0.156       
    ##                                                               Comp.3 Comp.4
    ## CHELSA_bio1_1981.2010_V.2.1                                                
    ## CHELSA_bio10_1981.2010_V.2.1                                               
    ## CHELSA_bio11_1981.2010_V.2.1                                               
    ## CHELSA_bio12_1981.2010_V.2.1                                         -0.175
    ## CHELSA_bio13_1981.2010_V.2.1                                   0.349 -0.281
    ## CHELSA_bio14_1981.2010_V.2.1                                  -0.201       
    ## CHELSA_bio15_1981.2010_V.2.1                                   0.388       
    ## CHELSA_bio16_1981.2010_V.2.1                                   0.340 -0.266
    ## CHELSA_bio17_1981.2010_V.2.1                                  -0.186       
    ## CHELSA_bio18_1981.2010_V.2.1                                   0.480 -0.200
    ## CHELSA_bio19_1981.2010_V.2.1                                  -0.200       
    ## CHELSA_bio2_1981.2010_V.2.1                                          -0.170
    ## CHELSA_bio3_1981.2010_V.2.1                                                
    ## CHELSA_bio4_1981.2010_V.2.1                                          -0.137
    ## CHELSA_bio5_1981.2010_V.2.1                                                
    ## CHELSA_bio6_1981.2010_V.2.1                                                
    ## CHELSA_bio7_1981.2010_V.2.1                                          -0.139
    ## CHELSA_bio8_1981.2010_V.2.1                                    0.328       
    ## CHELSA_bio9_1981.2010_V.2.1                                   -0.138       
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                       -0.105
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2 -0.185 -0.458
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.294  0.665
    ## wc2.1_5m_elev                                                        -0.145
    ##                                                               Comp.5 Comp.6
    ## CHELSA_bio1_1981.2010_V.2.1                                                
    ## CHELSA_bio10_1981.2010_V.2.1                                               
    ## CHELSA_bio11_1981.2010_V.2.1                                               
    ## CHELSA_bio12_1981.2010_V.2.1                                               
    ## CHELSA_bio13_1981.2010_V.2.1                                          0.181
    ## CHELSA_bio14_1981.2010_V.2.1                                         -0.202
    ## CHELSA_bio15_1981.2010_V.2.1                                          0.414
    ## CHELSA_bio16_1981.2010_V.2.1                                          0.139
    ## CHELSA_bio17_1981.2010_V.2.1                                         -0.180
    ## CHELSA_bio18_1981.2010_V.2.1                                         -0.156
    ## CHELSA_bio19_1981.2010_V.2.1                                  -0.108       
    ## CHELSA_bio2_1981.2010_V.2.1                                   -0.534 -0.133
    ## CHELSA_bio3_1981.2010_V.2.1                                   -0.261       
    ## CHELSA_bio4_1981.2010_V.2.1                                                
    ## CHELSA_bio5_1981.2010_V.2.1                                                
    ## CHELSA_bio6_1981.2010_V.2.1                                                
    ## CHELSA_bio7_1981.2010_V.2.1                                                
    ## CHELSA_bio8_1981.2010_V.2.1                                    0.153 -0.747
    ## CHELSA_bio9_1981.2010_V.2.1                                           0.262
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                       -0.112
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.250       
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2 -0.143       
    ## wc2.1_5m_elev                                                 -0.694       
    ##                                                               Comp.7 Comp.8
    ## CHELSA_bio1_1981.2010_V.2.1                                           0.101
    ## CHELSA_bio10_1981.2010_V.2.1                                          0.198
    ## CHELSA_bio11_1981.2010_V.2.1                                   0.120       
    ## CHELSA_bio12_1981.2010_V.2.1                                               
    ## CHELSA_bio13_1981.2010_V.2.1                                               
    ## CHELSA_bio14_1981.2010_V.2.1                                               
    ## CHELSA_bio15_1981.2010_V.2.1                                               
    ## CHELSA_bio16_1981.2010_V.2.1                                               
    ## CHELSA_bio17_1981.2010_V.2.1                                               
    ## CHELSA_bio18_1981.2010_V.2.1                                               
    ## CHELSA_bio19_1981.2010_V.2.1                                         -0.110
    ## CHELSA_bio2_1981.2010_V.2.1                                   -0.507       
    ## CHELSA_bio3_1981.2010_V.2.1                                          -0.184
    ## CHELSA_bio4_1981.2010_V.2.1                                   -0.266       
    ## CHELSA_bio5_1981.2010_V.2.1                                   -0.179  0.212
    ## CHELSA_bio6_1981.2010_V.2.1                                    0.155       
    ## CHELSA_bio7_1981.2010_V.2.1                                   -0.380       
    ## CHELSA_bio8_1981.2010_V.2.1                                                
    ## CHELSA_bio9_1981.2010_V.2.1                                                
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                 0.140  0.305
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.178 -0.718
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2 -0.112 -0.454
    ## wc2.1_5m_elev                                                  0.599       
    ##                                                               Comp.9 Comp.10
    ## CHELSA_bio1_1981.2010_V.2.1                                                 
    ## CHELSA_bio10_1981.2010_V.2.1                                          0.168 
    ## CHELSA_bio11_1981.2010_V.2.1                                                
    ## CHELSA_bio12_1981.2010_V.2.1                                          0.103 
    ## CHELSA_bio13_1981.2010_V.2.1                                  -0.120        
    ## CHELSA_bio14_1981.2010_V.2.1                                                
    ## CHELSA_bio15_1981.2010_V.2.1                                          0.149 
    ## CHELSA_bio16_1981.2010_V.2.1                                                
    ## CHELSA_bio17_1981.2010_V.2.1                                          0.108 
    ## CHELSA_bio18_1981.2010_V.2.1                                         -0.186 
    ## CHELSA_bio19_1981.2010_V.2.1                                          0.252 
    ## CHELSA_bio2_1981.2010_V.2.1                                    0.215        
    ## CHELSA_bio3_1981.2010_V.2.1                                    0.239 -0.314 
    ## CHELSA_bio4_1981.2010_V.2.1                                   -0.186  0.234 
    ## CHELSA_bio5_1981.2010_V.2.1                                           0.202 
    ## CHELSA_bio6_1981.2010_V.2.1                                                 
    ## CHELSA_bio7_1981.2010_V.2.1                                   -0.117  0.185 
    ## CHELSA_bio8_1981.2010_V.2.1                                   -0.213  0.312 
    ## CHELSA_bio9_1981.2010_V.2.1                                           0.519 
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                 0.823  0.300 
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.166  0.116 
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.151  0.295 
    ## wc2.1_5m_elev                                                 -0.182  0.174 
    ##                                                               Comp.11 Comp.12
    ## CHELSA_bio1_1981.2010_V.2.1                                    0.153         
    ## CHELSA_bio10_1981.2010_V.2.1                                   0.336   0.153 
    ## CHELSA_bio11_1981.2010_V.2.1                                                 
    ## CHELSA_bio12_1981.2010_V.2.1                                   0.122         
    ## CHELSA_bio13_1981.2010_V.2.1                                   0.125  -0.595 
    ## CHELSA_bio14_1981.2010_V.2.1                                                 
    ## CHELSA_bio15_1981.2010_V.2.1                                  -0.136         
    ## CHELSA_bio16_1981.2010_V.2.1                                                 
    ## CHELSA_bio17_1981.2010_V.2.1                                   0.103         
    ## CHELSA_bio18_1981.2010_V.2.1                                           0.669 
    ## CHELSA_bio19_1981.2010_V.2.1                                  -0.116   0.227 
    ## CHELSA_bio2_1981.2010_V.2.1                                                  
    ## CHELSA_bio3_1981.2010_V.2.1                                   -0.259  -0.138 
    ## CHELSA_bio4_1981.2010_V.2.1                                    0.122   0.114 
    ## CHELSA_bio5_1981.2010_V.2.1                                    0.359   0.130 
    ## CHELSA_bio6_1981.2010_V.2.1                                    0.116         
    ## CHELSA_bio7_1981.2010_V.2.1                                                  
    ## CHELSA_bio8_1981.2010_V.2.1                                   -0.201  -0.157 
    ## CHELSA_bio9_1981.2010_V.2.1                                   -0.609         
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                               
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.166         
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.282         
    ## wc2.1_5m_elev                                                  0.132         
    ##                                                               Comp.13 Comp.14
    ## CHELSA_bio1_1981.2010_V.2.1                                                  
    ## CHELSA_bio10_1981.2010_V.2.1                                                 
    ## CHELSA_bio11_1981.2010_V.2.1                                                 
    ## CHELSA_bio12_1981.2010_V.2.1                                   0.124   0.210 
    ## CHELSA_bio13_1981.2010_V.2.1                                                 
    ## CHELSA_bio14_1981.2010_V.2.1                                  -0.455  -0.476 
    ## CHELSA_bio15_1981.2010_V.2.1                                  -0.689         
    ## CHELSA_bio16_1981.2010_V.2.1                                   0.193         
    ## CHELSA_bio17_1981.2010_V.2.1                                  -0.278  -0.189 
    ## CHELSA_bio18_1981.2010_V.2.1                                          -0.209 
    ## CHELSA_bio19_1981.2010_V.2.1                                  -0.251   0.705 
    ## CHELSA_bio2_1981.2010_V.2.1                                                  
    ## CHELSA_bio3_1981.2010_V.2.1                                            0.150 
    ## CHELSA_bio4_1981.2010_V.2.1                                                  
    ## CHELSA_bio5_1981.2010_V.2.1                                                  
    ## CHELSA_bio6_1981.2010_V.2.1                                                  
    ## CHELSA_bio7_1981.2010_V.2.1                                                  
    ## CHELSA_bio8_1981.2010_V.2.1                                            0.125 
    ## CHELSA_bio9_1981.2010_V.2.1                                    0.271  -0.318 
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                               
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2  0.111         
    ## wc2.1_5m_elev                                                                
    ##                                                               Comp.15 Comp.16
    ## CHELSA_bio1_1981.2010_V.2.1                                                  
    ## CHELSA_bio10_1981.2010_V.2.1                                   0.145   0.229 
    ## CHELSA_bio11_1981.2010_V.2.1                                                 
    ## CHELSA_bio12_1981.2010_V.2.1                                  -0.286  -0.248 
    ## CHELSA_bio13_1981.2010_V.2.1                                   0.502  -0.152 
    ## CHELSA_bio14_1981.2010_V.2.1                                           0.290 
    ## CHELSA_bio15_1981.2010_V.2.1                                  -0.146  -0.108 
    ## CHELSA_bio16_1981.2010_V.2.1                                  -0.626   0.390 
    ## CHELSA_bio17_1981.2010_V.2.1                                          -0.180 
    ## CHELSA_bio18_1981.2010_V.2.1                                   0.280  -0.114 
    ## CHELSA_bio19_1981.2010_V.2.1                                   0.119         
    ## CHELSA_bio2_1981.2010_V.2.1                                   -0.151  -0.338 
    ## CHELSA_bio3_1981.2010_V.2.1                                    0.213   0.536 
    ## CHELSA_bio4_1981.2010_V.2.1                                    0.172   0.345 
    ## CHELSA_bio5_1981.2010_V.2.1                                                  
    ## CHELSA_bio6_1981.2010_V.2.1                                                  
    ## CHELSA_bio7_1981.2010_V.2.1                                            0.127 
    ## CHELSA_bio8_1981.2010_V.2.1                                                  
    ## CHELSA_bio9_1981.2010_V.2.1                                                  
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                               
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## wc2.1_5m_elev                                                                
    ##                                                               Comp.17 Comp.18
    ## CHELSA_bio1_1981.2010_V.2.1                                                  
    ## CHELSA_bio10_1981.2010_V.2.1                                                 
    ## CHELSA_bio11_1981.2010_V.2.1                                                 
    ## CHELSA_bio12_1981.2010_V.2.1                                   0.753         
    ## CHELSA_bio13_1981.2010_V.2.1                                           0.114 
    ## CHELSA_bio14_1981.2010_V.2.1                                           0.492 
    ## CHELSA_bio15_1981.2010_V.2.1                                   0.173         
    ## CHELSA_bio16_1981.2010_V.2.1                                  -0.310  -0.130 
    ## CHELSA_bio17_1981.2010_V.2.1                                          -0.778 
    ## CHELSA_bio18_1981.2010_V.2.1                                                 
    ## CHELSA_bio19_1981.2010_V.2.1                                  -0.283   0.160 
    ## CHELSA_bio2_1981.2010_V.2.1                                   -0.215         
    ## CHELSA_bio3_1981.2010_V.2.1                                    0.337  -0.192 
    ## CHELSA_bio4_1981.2010_V.2.1                                           -0.166 
    ## CHELSA_bio5_1981.2010_V.2.1                                                  
    ## CHELSA_bio6_1981.2010_V.2.1                                   -0.112         
    ## CHELSA_bio7_1981.2010_V.2.1                                    0.188         
    ## CHELSA_bio8_1981.2010_V.2.1                                                  
    ## CHELSA_bio9_1981.2010_V.2.1                                                  
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                               
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## wc2.1_5m_elev                                                                
    ##                                                               Comp.19 Comp.20
    ## CHELSA_bio1_1981.2010_V.2.1                                    0.214   0.451 
    ## CHELSA_bio10_1981.2010_V.2.1                                   0.429   0.128 
    ## CHELSA_bio11_1981.2010_V.2.1                                           0.278 
    ## CHELSA_bio12_1981.2010_V.2.1                                   0.147  -0.125 
    ## CHELSA_bio13_1981.2010_V.2.1                                                 
    ## CHELSA_bio14_1981.2010_V.2.1                                                 
    ## CHELSA_bio15_1981.2010_V.2.1                                                 
    ## CHELSA_bio16_1981.2010_V.2.1                                                 
    ## CHELSA_bio17_1981.2010_V.2.1                                           0.120 
    ## CHELSA_bio18_1981.2010_V.2.1                                                 
    ## CHELSA_bio19_1981.2010_V.2.1                                                 
    ## CHELSA_bio2_1981.2010_V.2.1                                    0.256         
    ## CHELSA_bio3_1981.2010_V.2.1                                   -0.155         
    ## CHELSA_bio4_1981.2010_V.2.1                                    0.340  -0.409 
    ## CHELSA_bio5_1981.2010_V.2.1                                   -0.616  -0.275 
    ## CHELSA_bio6_1981.2010_V.2.1                                           -0.429 
    ## CHELSA_bio7_1981.2010_V.2.1                                   -0.381   0.471 
    ## CHELSA_bio8_1981.2010_V.2.1                                                  
    ## CHELSA_bio9_1981.2010_V.2.1                                                  
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                               
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## wc2.1_5m_elev                                                                
    ##                                                               Comp.21 Comp.22
    ## CHELSA_bio1_1981.2010_V.2.1                                    0.758         
    ## CHELSA_bio10_1981.2010_V.2.1                                  -0.510   0.353 
    ## CHELSA_bio11_1981.2010_V.2.1                                  -0.297  -0.821 
    ## CHELSA_bio12_1981.2010_V.2.1                                                 
    ## CHELSA_bio13_1981.2010_V.2.1                                                 
    ## CHELSA_bio14_1981.2010_V.2.1                                                 
    ## CHELSA_bio15_1981.2010_V.2.1                                                 
    ## CHELSA_bio16_1981.2010_V.2.1                                                 
    ## CHELSA_bio17_1981.2010_V.2.1                                                 
    ## CHELSA_bio18_1981.2010_V.2.1                                                 
    ## CHELSA_bio19_1981.2010_V.2.1                                                 
    ## CHELSA_bio2_1981.2010_V.2.1                                                  
    ## CHELSA_bio3_1981.2010_V.2.1                                                  
    ## CHELSA_bio4_1981.2010_V.2.1                                    0.212  -0.447 
    ## CHELSA_bio5_1981.2010_V.2.1                                                  
    ## CHELSA_bio6_1981.2010_V.2.1                                    0.102         
    ## CHELSA_bio7_1981.2010_V.2.1                                                  
    ## CHELSA_bio8_1981.2010_V.2.1                                                  
    ## CHELSA_bio9_1981.2010_V.2.1                                                  
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                               
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2                
    ## wc2.1_5m_elev                                                                
    ##                                                               Comp.23
    ## CHELSA_bio1_1981.2010_V.2.1                                          
    ## CHELSA_bio10_1981.2010_V.2.1                                         
    ## CHELSA_bio11_1981.2010_V.2.1                                         
    ## CHELSA_bio12_1981.2010_V.2.1                                         
    ## CHELSA_bio13_1981.2010_V.2.1                                         
    ## CHELSA_bio14_1981.2010_V.2.1                                         
    ## CHELSA_bio15_1981.2010_V.2.1                                         
    ## CHELSA_bio16_1981.2010_V.2.1                                         
    ## CHELSA_bio17_1981.2010_V.2.1                                         
    ## CHELSA_bio18_1981.2010_V.2.1                                         
    ## CHELSA_bio19_1981.2010_V.2.1                                         
    ## CHELSA_bio2_1981.2010_V.2.1                                          
    ## CHELSA_bio3_1981.2010_V.2.1                                          
    ## CHELSA_bio4_1981.2010_V.2.1                                          
    ## CHELSA_bio5_1981.2010_V.2.1                                   -0.373 
    ## CHELSA_bio6_1981.2010_V.2.1                                    0.784 
    ## CHELSA_bio7_1981.2010_V.2.1                                    0.496 
    ## CHELSA_bio8_1981.2010_V.2.1                                          
    ## CHELSA_bio9_1981.2010_V.2.1                                          
    ## CHELSA_gddlgd0_1981.2010_V.2.1                                       
    ## sol_clay.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2        
    ## sol_sand.wfraction_usda.3a1a1a_m_250m_b0..0cm_1950..2017_v0.2        
    ## wc2.1_5m_elev                                                        
    ## 
    ##                Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7 Comp.8 Comp.9
    ## SS loadings     1.000  1.000  1.000  1.000  1.000  1.000  1.000  1.000  1.000
    ## Proportion Var  0.043  0.043  0.043  0.043  0.043  0.043  0.043  0.043  0.043
    ## Cumulative Var  0.043  0.087  0.130  0.174  0.217  0.261  0.304  0.348  0.391
    ##                Comp.10 Comp.11 Comp.12 Comp.13 Comp.14 Comp.15 Comp.16 Comp.17
    ## SS loadings      1.000   1.000   1.000   1.000   1.000   1.000   1.000   1.000
    ## Proportion Var   0.043   0.043   0.043   0.043   0.043   0.043   0.043   0.043
    ## Cumulative Var   0.435   0.478   0.522   0.565   0.609   0.652   0.696   0.739
    ##                Comp.18 Comp.19 Comp.20 Comp.21 Comp.22 Comp.23
    ## SS loadings      1.000   1.000   1.000   1.000   1.000   1.000
    ## Proportion Var   0.043   0.043   0.043   0.043   0.043   0.043
    ## Cumulative Var   0.783   0.826   0.870   0.913   0.957   1.000

``` r
propvar <- summary(pca$model)$sdev^2/sum(summary(pca$model)$sdev^2)
cumprop <- cumsum(propvar)
plot(1:23,cumprop, type = "o", xlab = "PCs", ylim = c(0,1))
points(1:23, propvar, type = "o", col = "red")
```

![](Aa_SDM_files/figure-gfm/PCA-1.png)<!-- -->

``` r
features <- cumprop < 0.96
env.pca <- subset(pca$map, which(features))

plot(env.pca)
```

![](Aa_SDM_files/figure-gfm/PCA-2.png)<!-- -->

``` r
env <- env.pca #addLayer(env.pca, env.cat.m)

plot(env[[1]], main = "PC 1 | A. asteris occ.")
points(occt, pch = 16, cex = 0.5)
```

![](Aa_SDM_files/figure-gfm/PCA-3.png)<!-- -->

### Extract env. and MESS

``` r
occp <- st_as_sf(occt, coords = c("lon", "lat"), crs = proj.crs)
occe <- raster::extract(env, occp)


occ.sim <- similarity(env, occe)
occ.mess <- occ.sim$similarity_min
occ.sp <- as_Spatial(occp)

myScale <- seq(cellStats(occ.mess, min), cellStats(occ.mess, max), length.out = 100)
rasterVis::levelplot(occ.mess, main = "Environmental similarity", at = myScale, margin = FALSE) + 
  latticeExtra::layer(sp.points(occ.sp, col="black"))
```

![](Aa_SDM_files/figure-gfm/MESS-1.png)<!-- -->

``` r
### Remove occ which correspond to NA env
occ.nna <- raster::extract(env, occp, sp = TRUE) %>% 
  as.data.frame() %>% 
  filter(!is.na(PC1))

occe <- occ.nna[,1:(ncol(occ.nna)-2)]
occp <- occ.nna[,(ncol(occ.nna)-1):ncol(occ.nna)]
```

In total 0 observations fall occur where the environmental data is NA.
These correspond to islands and shoreline sites which are lost during
masking due to the resolution used.

## SDM Model

### Assign random background points

Model performance improved using background points sampled from full
study area over buffered convex hull.

``` r
 # crs(env) <- crs(proj.crs)
 # sa.bg <- st_convex_hull(occs.sf) %>% 
 #   st_buffer(occs.sf, dist = 100000) %>% 
 #   st_sf()
 # 
 # env.bg <- raster::mask(env, sa.bg)

(points <- sum(!is.na(getValues(env[[1]])))*0.3)
```

    ## [1] 24664.5

``` r
bg <- dismo::randomPoints(env, n = points) %>% 
  as.data.frame()
colnames(bg) <- colnames(occp)

plot(env[[1]], main="PC1 | Background Points")
points(bg, pch = 20, cex = 0.05)
```

![](Aa_SDM_files/figure-gfm/Bkg.%20points-1.png)<!-- -->

### Partition for Model Evaluation

Use block partitioning since n = 38.

``` r
block <- get.block(occp, bg)
evalplot.grps(pts = occp, pts.grp = block$occs.grp, envs = env)
```

    ## Plotting first raster in stack...

![](Aa_SDM_files/figure-gfm/Partition-1.png)<!-- -->

### Tune Model

``` r
e.mx <- ENMevaluate(occs = occp, envs = env, bg = bg, #categoricals = names(env.cat), 
                    algorithm = 'maxent.jar', partitions = 'jackknife', parallel = TRUE, 
                    tune.args = list(fc = c("L","Q","LQ","LQH","H"), rm = 1:5))
```

    ## Package ecospat is not installed, so Continuous Boyce Index (CBI) cannot be calculated.

    ## *** Running initial checks... ***

    ## * Clamping predictor variable rasters...

    ## * Model evaluations with k-1 jackknife (leave-one-out) cross validation...

    ## 
    ## *** Running ENMeval v2.0.3 with maxent.jar v3.4.3 from dismo package v1.3.5 ***

    ##   |                                                                              |                                                                      |   0%

    ## 
    ## Of 16 total cores using 16...

    ## Running in parallel using doSNOW...

    ##   |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

    ## ENMevaluate completed in 107 minutes 14.2 seconds.

``` r
e.mx
```

    ## An object of class:  ENMevaluation 
    ##  occurrence/background points:  38 / 24664 
    ##  partition method:  jackknife 
    ##  partition settings:  none 
    ##  clamp:  left: PC1, PC2, PC3, PC4, PC5, PC6, PC7
    ##          right: PC1, PC2, PC3, PC4, PC5, PC6, PC7 
    ##  categoricals:   
    ##  algorithm:  maxent.jar 
    ##  tune settings:  fc: L,Q,LQ,LQH,H
    ##                  rm: 1,2,3,4,5 
    ##  overlap:  TRUE 
    ## Refer to ?ENMevaluation for information on slots.

### Model Selection

AUC: Radosavljevic and Anderson 2013

``` r
evalplot.stats(e = e.mx, stats = c("auc.val"), color = "fc", x.var = "rm", 
               error.bars = FALSE)
```

![](Aa_SDM_files/figure-gfm/Model%20selection-1.png)<!-- -->

``` r
res <- eval.results(e.mx)
kable(res)
```

| fc  | rm  | tune.args   | auc.train | cbi.train | auc.diff.avg | auc.diff.sd | auc.val.avg | auc.val.sd | cbi.val.avg | cbi.val.sd | or.10p.avg | or.10p.sd | or.mtp.avg | or.mtp.sd |      AICc |  delta.AICc |     w.AIC | ncoef |
|:----|:----|:------------|----------:|:----------|-------------:|------------:|------------:|-----------:|------------:|-----------:|-----------:|----------:|-----------:|----------:|----------:|------------:|----------:|------:|
| L   | 1   | fc.L_rm.1   | 0.8175569 | NA        |    0.1909747 |   1.0985148 |   0.7960414 |  1.5540718 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  820.8047 |  36.5554860 | 0.0000000 |     7 |
| Q   | 1   | fc.Q_rm.1   | 0.8236931 | NA        |    0.1187292 |   0.5266742 |   0.8013123 |  0.8611596 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |  816.6736 |  32.4244165 | 0.0000000 |     6 |
| LQ  | 1   | fc.LQ_rm.1  | 0.9028368 | NA        |    0.1222380 |   0.7520367 |   0.8765540 |  1.0204791 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |  791.2117 |   6.9625192 | 0.0057452 |    14 |
| LQH | 1   | fc.LQH_rm.1 | 0.9065904 | NA        |    0.1211160 |   0.7402393 |   0.8738343 |  1.0005964 |          NA |         NA |  0.1052632 |  1.866752 |  0.0789474 | 1.6402571 |  975.3283 | 191.0791315 | 0.0000000 |    28 |
| H   | 1   | fc.H_rm.1   | 0.9061268 | NA        |    0.1223988 |   0.7739388 |   0.8740851 |  1.0365462 |          NA |         NA |  0.1315789 |  2.056171 |  0.0789474 | 1.6402571 |        NA |          NA |        NA |    63 |
| L   | 2   | fc.L_rm.2   | 0.8178802 | NA        |    0.1910074 |   1.1021421 |   0.7946944 |  1.5547957 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |  821.6490 |  37.3997423 | 0.0000000 |     7 |
| Q   | 2   | fc.Q_rm.2   | 0.8194519 | NA        |    0.1202551 |   0.5168669 |   0.7998633 |  0.8632881 |          NA |         NA |  0.1315789 |  2.056171 |  0.0263158 | 0.9736842 |  815.1566 |  30.9073569 | 0.0000000 |     5 |
| LQ  | 2   | fc.LQ_rm.2  | 0.9029904 | NA        |    0.1179281 |   0.7510209 |   0.8792823 |  1.0037580 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  785.1648 |   0.9156036 | 0.1181344 |    12 |
| LQH | 2   | fc.LQH_rm.2 | 0.9030523 | NA        |    0.1179854 |   0.7474564 |   0.8781396 |  1.0004895 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  784.8075 |   0.5583368 | 0.1412394 |    12 |
| H   | 2   | fc.H_rm.2   | 0.9000594 | NA        |    0.1282928 |   0.7539119 |   0.8717207 |  1.0513227 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |        NA |          NA |        NA |    61 |
| L   | 3   | fc.L_rm.3   | 0.8174641 | NA        |    0.1918554 |   1.1040352 |   0.7922953 |  1.5578492 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |  823.1486 |  38.8993793 | 0.0000000 |     7 |
| Q   | 3   | fc.Q_rm.3   | 0.8175356 | NA        |    0.1203622 |   0.5060811 |   0.7975491 |  0.8575384 |          NA |         NA |  0.1578947 |  2.218032 |  0.0263158 | 0.9736842 |  813.2862 |  29.0369880 | 0.0000001 |     4 |
| LQ  | 3   | fc.LQ_rm.3  | 0.9022275 | NA        |    0.1153430 |   0.7456404 |   0.8791756 |  0.9888773 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  786.1512 |   1.9019577 | 0.0721427 |    11 |
| LQH | 3   | fc.LQH_rm.3 | 0.9022062 | NA        |    0.1152448 |   0.7460926 |   0.8789121 |  0.9886031 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  786.1524 |   1.9031661 | 0.0720991 |    11 |
| H   | 3   | fc.H_rm.3   | 0.8959889 | NA        |    0.1325025 |   0.7620095 |   0.8641729 |  1.0724649 |          NA |         NA |  0.1578947 |  2.218032 |  0.0263158 | 0.9736842 |        NA |          NA |        NA |    48 |
| L   | 4   | fc.L_rm.4   | 0.8159223 | NA        |    0.1920012 |   1.1128587 |   0.7885492 |  1.5617733 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |  825.2254 |  40.9762100 | 0.0000000 |     7 |
| Q   | 4   | fc.Q_rm.4   | 0.8154368 | NA        |    0.1211667 |   0.4936931 |   0.7924073 |  0.8539213 |          NA |         NA |  0.1578947 |  2.218032 |  0.0263158 | 0.9736842 |  814.4309 |  30.1816423 | 0.0000001 |     4 |
| LQ  | 4   | fc.LQ_rm.4  | 0.9011627 | NA        |    0.1143027 |   0.7329058 |   0.8794039 |  0.9760434 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  784.2495 |   0.0002983 | 0.1866948 |     9 |
| LQH | 4   | fc.LQH_rm.4 | 0.9011595 | NA        |    0.1142941 |   0.7327920 |   0.8794114 |  0.9759228 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  784.2492 |   0.0000000 | 0.1867227 |     9 |
| H   | 4   | fc.H_rm.4   | 0.8873769 | NA        |    0.1348794 |   0.7699321 |   0.8576388 |  1.0894834 |          NA |         NA |  0.1315789 |  2.056171 |  0.0263158 | 0.9736842 | 1178.6737 | 394.4244829 | 0.0000000 |    31 |
| L   | 5   | fc.L_rm.5   | 0.8138695 | NA        |    0.1927027 |   1.1054704 |   0.7877431 |  1.5589134 |          NA |         NA |  0.1315789 |  2.056171 |  0.0526316 | 1.3582619 |  824.6020 |  40.3527634 | 0.0000000 |     6 |
| Q   | 5   | fc.Q_rm.5   | 0.8101415 | NA        |    0.1229337 |   0.4764793 |   0.7921256 |  0.8545870 |          NA |         NA |  0.1315789 |  2.056171 |  0.0263158 | 0.9736842 |  818.6785 |  34.4293138 | 0.0000000 |     5 |
| LQ  | 5   | fc.LQ_rm.5  | 0.9007850 | NA        |    0.1125002 |   0.7197214 |   0.8792524 |  0.9591015 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  785.3330 |   1.0838097 | 0.1086052 |     8 |
| LQH | 5   | fc.LQH_rm.5 | 0.9007860 | NA        |    0.1124973 |   0.7197280 |   0.8792513 |  0.9590934 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 | 1.3582619 |  785.3328 |   1.0836047 | 0.1086164 |     8 |
| H   | 5   | fc.H_rm.5   | 0.8817555 | NA        |    0.1401378 |   0.7767999 |   0.8441608 |  1.1102636 |          NA |         NA |  0.1578947 |  2.218032 |  0.0263158 | 0.9736842 |        NA |          NA |        NA |    40 |

``` r
opt.auc <- res %>% 
  slice_max(auc.val.avg)
kable(opt.auc)
```

| fc  | rm  | tune.args   | auc.train | cbi.train | auc.diff.avg | auc.diff.sd | auc.val.avg | auc.val.sd | cbi.val.avg | cbi.val.sd | or.10p.avg | or.10p.sd | or.mtp.avg | or.mtp.sd |     AICc | delta.AICc |     w.AIC | ncoef |
|:----|:----|:------------|----------:|:----------|-------------:|------------:|------------:|-----------:|------------:|-----------:|-----------:|----------:|-----------:|----------:|---------:|-----------:|----------:|------:|
| LQH | 4   | fc.LQH_rm.4 | 0.9011595 | NA        |    0.1142941 |    0.732792 |   0.8794114 |  0.9759228 |          NA |         NA |  0.1052632 |  1.866752 |  0.0526316 |  1.358262 | 784.2492 |          0 | 0.1867227 |     9 |

``` r
dismo::response(eval.models(e.mx)[[opt.auc$tune.args]])
```

![](Aa_SDM_files/figure-gfm/Model%20selection-2.png)<!-- -->

``` r
plot(eval.models(e.mx)[[opt.auc$tune.args]])
```

![](Aa_SDM_files/figure-gfm/Model%20selection-3.png)<!-- -->

## Predictions

``` r
dir.create("../pred/dist")
pred <- eval.predictions(e.mx)[[opt.auc$tune.args]]
plot(pred)
```

![](Aa_SDM_files/figure-gfm/Prediction-1.png)<!-- -->

``` r
writeRaster(pred, "../pred/dist/Aa_dist.tif", overwrite = TRUE)
```

### Threshold

``` r
pred.vals <- raster::extract(pred, occp)
or.min.threshold <- pred > min(pred.vals)
writeRaster(or.min.threshold, "../pred/dist/Aa_dist_min.tif", overwrite = TRUE)
plot(or.min.threshold)
```

![](Aa_SDM_files/figure-gfm/Threshold-1.png)<!-- -->

``` r
n10 <- ceiling(length(pred.vals) * 0.1)
or.10.threshold <- pred > sort(pred.vals)[n10]
writeRaster(or.10.threshold, "../pred/dist/Aa_dist_p10.tif", overwrite = TRUE)
plot(or.10.threshold)
```

![](Aa_SDM_files/figure-gfm/Threshold-2.png)<!-- -->

### ENM Null

``` r
mod.null <- ENMnulls(e.mx, mod.settings = list(fc = as.character(opt.auc$fc), rm = as.numeric(opt.auc$rm)), no.iter = 100, parallel = TRUE, quiet = TRUE)
null.results(mod.null)
```

    ##      fc rm   tune.args auc.train cbi.train auc.diff.avg auc.diff.sd auc.val.avg
    ## 1   LQH  4 fc.LQH_rm.4 0.5462671        NA    0.2144757   0.1542071   0.5015540
    ## 2   LQH  4 fc.LQH_rm.4 0.6308534        NA    0.2494490   0.1114616   0.8106387
    ## 3   LQH  4 fc.LQH_rm.4 0.6714213        NA    0.2853100   0.1907195   0.4845615
    ## 4   LQH  4 fc.LQH_rm.4 0.5870836        NA    0.2772284   0.1545894   0.4172505
    ## 5   LQH  4 fc.LQH_rm.4 0.6309569        NA    0.2656044   0.1626522   0.4062255
    ## 6   LQH  4 fc.LQH_rm.4 0.6157531        NA    0.2751955   0.1664584   0.4752271
    ## 7   LQH  4 fc.LQH_rm.4 0.6171418        NA    0.2166481   0.1201522   0.6298083
    ## 8   LQH  4 fc.LQH_rm.4 0.6356014        NA    0.3454793   0.1860429   0.3412341
    ## 9   LQH  4 fc.LQH_rm.4 0.6301673        NA    0.2098605   0.1271759   0.7055612
    ## 10  LQH  4 fc.LQH_rm.4 0.6080768        NA    0.2234829   0.1056521   0.7765137
    ## 11  LQH  4 fc.LQH_rm.4 0.6386167        NA    0.3035891   0.1867330   0.3466863
    ## 12  LQH  4 fc.LQH_rm.4 0.5619841        NA    0.2331627   0.1542314   0.6449300
    ## 13  LQH  4 fc.LQH_rm.4 0.6011809        NA    0.2831078   0.1790843   0.3850242
    ## 14  LQH  4 fc.LQH_rm.4 0.6172623        NA    0.3223758   0.1751218   0.3525104
    ## 15  LQH  4 fc.LQH_rm.4 0.6399515        NA    0.2835470   0.1678589   0.4381845
    ## 16  LQH  4 fc.LQH_rm.4 0.6174613        NA    0.2238371   0.1429847   0.7244700
    ## 17  LQH  4 fc.LQH_rm.4 0.6674932        NA    0.3559148   0.2081398   0.3225215
    ## 18  LQH  4 fc.LQH_rm.4 0.6288011        NA    0.2544396   0.1806469   0.5685151
    ## 19  LQH  4 fc.LQH_rm.4 0.5826007        NA    0.2521973   0.1832372   0.4624383
    ## 20  LQH  4 fc.LQH_rm.4 0.6855987        NA    0.2592031   0.1787824   0.5452167
    ## 21  LQH  4 fc.LQH_rm.4 0.5742346        NA    0.2415869   0.1747990   0.4697305
    ## 22  LQH  4 fc.LQH_rm.4 0.6682102        NA    0.2283478   0.1747793   0.6198855
    ## 23  LQH  4 fc.LQH_rm.4 0.6579006        NA    0.3306057   0.1755167   0.3727823
    ## 24  LQH  4 fc.LQH_rm.4 0.6788186        NA    0.1775624   0.1265215   0.6299844
    ## 25  LQH  4 fc.LQH_rm.4 0.6415407        NA    0.3379878   0.1903372   0.3409620
    ## 26  LQH  4 fc.LQH_rm.4 0.6677397        NA    0.2648937   0.1836357   0.4908966
    ## 27  LQH  4 fc.LQH_rm.4 0.6104401        NA    0.2399924   0.1516970   0.4730339
    ## 28  LQH  4 fc.LQH_rm.4 0.6459121        NA    0.2035798   0.1406667   0.5894938
    ## 29  LQH  4 fc.LQH_rm.4 0.6971155        NA    0.3592055   0.2314476   0.3986382
    ## 30  LQH  4 fc.LQH_rm.4 0.6434693        NA    0.2033198   0.1344350   0.5655883
    ## 31  LQH  4 fc.LQH_rm.4 0.5634624        NA    0.2157822   0.1242452   0.6173850
    ## 32  LQH  4 fc.LQH_rm.4 0.6474672        NA    0.2257689   0.1576051   0.5744986
    ## 33  LQH  4 fc.LQH_rm.4 0.6216278        NA    0.3416114   0.1518974   0.3792972
    ## 34  LQH  4 fc.LQH_rm.4 0.6217521        NA    0.3276567   0.1584757   0.5311540
    ## 35  LQH  4 fc.LQH_rm.4 0.6338425        NA    0.2296226   0.1621395   0.6234033
    ## 36  LQH  4 fc.LQH_rm.4 0.6008832        NA    0.2868941   0.1540383   0.5392667
    ## 37  LQH  4 fc.LQH_rm.4 0.6036648        NA    0.3991908   0.1629376   0.2250462
    ## 38  LQH  4 fc.LQH_rm.4 0.6277917        NA    0.2451167   0.1977900   0.5375334
    ## 39  LQH  4 fc.LQH_rm.4 0.5419283        NA    0.2012380   0.1439398   0.4295948
    ## 40  LQH  4 fc.LQH_rm.4 0.6429731        NA    0.4128451   0.1786572   0.2344729
    ## 41  LQH  4 fc.LQH_rm.4 0.6234433        NA    0.2672568   0.1725473   0.5653920
    ## 42  LQH  4 fc.LQH_rm.4 0.6103537        NA    0.2438538   0.1493484   0.4914685
    ## 43  LQH  4 fc.LQH_rm.4 0.6769082        NA    0.2849277   0.1817902   0.6050823
    ## 44  LQH  4 fc.LQH_rm.4 0.7442800        NA    0.4964236   0.1262752   0.2472008
    ## 45  LQH  4 fc.LQH_rm.4 0.5982772        NA    0.3680045   0.1714360   0.3054671
    ## 46  LQH  4 fc.LQH_rm.4 0.6559721        NA    0.2470879   0.1755631   0.4917219
    ## 47  LQH  4 fc.LQH_rm.4 0.6480530        NA    0.2394558   0.1901142   0.4563070
    ## 48  LQH  4 fc.LQH_rm.4 0.6299646        NA    0.2838458   0.1608958   0.6051090
    ## 49  LQH  4 fc.LQH_rm.4 0.6287003        NA    0.2558563   0.1752806   0.4850800
    ## 50  LQH  4 fc.LQH_rm.4 0.6121926        NA    0.2306153   0.1537243   0.5852046
    ## 51  LQH  4 fc.LQH_rm.4 0.6191551        NA    0.2773315   0.1450815   0.4480636
    ## 52  LQH  4 fc.LQH_rm.4 0.6248880        NA    0.3341182   0.1699793   0.3375039
    ## 53  LQH  4 fc.LQH_rm.4 0.6386903        NA    0.3178343   0.1697786   0.3994929
    ## 54  LQH  4 fc.LQH_rm.4 0.6675338        NA    0.3787473   0.1873598   0.3380449
    ## 55  LQH  4 fc.LQH_rm.4 0.6599369        NA    0.2306443   0.1219292   0.4465474
    ## 56  LQH  4 fc.LQH_rm.4 0.6391150        NA    0.2317974   0.1485739   0.6785412
    ## 57  LQH  4 fc.LQH_rm.4 0.6098266        NA    0.2688363   0.1748330   0.5589870
    ## 58  LQH  4 fc.LQH_rm.4 0.6233883        NA    0.2203361   0.1249042   0.6376543
    ## 59  LQH  4 fc.LQH_rm.4 0.6458902        NA    0.1299404   0.1309496   0.5881889
    ## 60  LQH  4 fc.LQH_rm.4 0.6366641        NA    0.2376799   0.1193850   0.7776970
    ## 61  LQH  4 fc.LQH_rm.4 0.6609185        NA    0.2967341   0.2008619   0.4079796
    ## 62  LQH  4 fc.LQH_rm.4 0.5795433        NA    0.3336565   0.1659929   0.3047842
    ## 63  LQH  4 fc.LQH_rm.4 0.6435536        NA    0.2909134   0.1689784   0.4603364
    ## 64  LQH  4 fc.LQH_rm.4 0.6151118        NA    0.2843849   0.1702068   0.5418152
    ## 65  LQH  4 fc.LQH_rm.4 0.6142257        NA    0.2989556   0.1647740   0.3655152
    ## 66  LQH  4 fc.LQH_rm.4 0.6220946        NA    0.2980050   0.1856213   0.5512888
    ## 67  LQH  4 fc.LQH_rm.4 0.6262126        NA    0.3205086   0.1418746   0.4540738
    ## 68  LQH  4 fc.LQH_rm.4 0.6334461        NA    0.4065039   0.1685285   0.2589087
    ## 69  LQH  4 fc.LQH_rm.4 0.6794796        NA    0.3415186   0.1992294   0.4106822
    ## 70  LQH  4 fc.LQH_rm.4 0.6635241        NA    0.3471024   0.2038951   0.3478792
    ## 71  LQH  4 fc.LQH_rm.4 0.6271628        NA    0.3474969   0.1752899   0.3090569
    ## 72  LQH  4 fc.LQH_rm.4 0.6729081        NA    0.3296904   0.2167880   0.4146823
    ## 73  LQH  4 fc.LQH_rm.4 0.6686818        NA    0.2118527   0.1570110   0.6057886
    ## 74  LQH  4 fc.LQH_rm.4 0.6068540        NA    0.2606381   0.1818421   0.4927878
    ## 75  LQH  4 fc.LQH_rm.4 0.6236770        NA    0.3249709   0.1486435   0.5453010
    ## 76  LQH  4 fc.LQH_rm.4 0.6373395        NA    0.1818190   0.1249690   0.5969370
    ## 77  LQH  4 fc.LQH_rm.4 0.5726053        NA    0.2767058   0.1031968   0.8112596
    ## 78  LQH  4 fc.LQH_rm.4 0.6195857        NA    0.2044135   0.1460151   0.5770988
    ## 79  LQH  4 fc.LQH_rm.4 0.5471521        NA    0.2122342   0.1521231   0.4168706
    ## 80  LQH  4 fc.LQH_rm.4 0.6518322        NA    0.2948063   0.2049531   0.4061919
    ## 81  LQH  4 fc.LQH_rm.4 0.5891722        NA    0.2345766   0.1633207   0.4797921
    ## 82  LQH  4 fc.LQH_rm.4 0.6019294        NA    0.3815800   0.1347718   0.2239312
    ## 83  LQH  4 fc.LQH_rm.4 0.6728057        NA    0.2533772   0.1674133   0.4806883
    ## 84  LQH  4 fc.LQH_rm.4 0.6916655        NA    0.3284178   0.1816126   0.4029819
    ## 85  LQH  4 fc.LQH_rm.4 0.6563690        NA    0.2694389   0.1723127   0.4910977
    ## 86  LQH  4 fc.LQH_rm.4 0.6273687        NA    0.2872507   0.1637317   0.5372955
    ## 87  LQH  4 fc.LQH_rm.4 0.5845618        NA    0.1573195   0.1404826   0.5526919
    ## 88  LQH  4 fc.LQH_rm.4 0.6191359        NA    0.2436953   0.1192353   0.6974362
    ## 89  LQH  4 fc.LQH_rm.4 0.6204451        NA    0.3537186   0.1683200   0.3418081
    ## 90  LQH  4 fc.LQH_rm.4 0.6817384        NA    0.3828538   0.1894168   0.3546545
    ## 91  LQH  4 fc.LQH_rm.4 0.6749385        NA    0.3137788   0.2092393   0.3916357
    ## 92  LQH  4 fc.LQH_rm.4 0.6838590        NA    0.1666399   0.1296979   0.6681307
    ## 93  LQH  4 fc.LQH_rm.4 0.6405458        NA    0.3466691   0.1704988   0.3734118
    ## 94  LQH  4 fc.LQH_rm.4 0.6319150        NA    0.3268515   0.1687659   0.3614489
    ## 95  LQH  4 fc.LQH_rm.4 0.5782784        NA    0.2561382   0.1492963   0.4798236
    ## 96  LQH  4 fc.LQH_rm.4 0.6375326        NA    0.1947821   0.1184724   0.6839411
    ## 97  LQH  4 fc.LQH_rm.4 0.6107591        NA    0.3539775   0.1315880   0.2553130
    ## 98  LQH  4 fc.LQH_rm.4 0.5931669        NA    0.2311998   0.1227764   0.8163768
    ## 99  LQH  4 fc.LQH_rm.4 0.6693167        NA    0.2950264   0.1619143   0.4156020
    ## 100 LQH  4 fc.LQH_rm.4 0.6803145        NA    0.4581535   0.1969992   0.2625503
    ##     auc.val.sd cbi.val.avg cbi.val.sd or.10p.avg or.10p.sd or.mtp.avg or.mtp.sd
    ## 1    0.2648603          NA         NA 0.10526316 0.3110117 0.00000000 0.0000000
    ## 2    0.2045491          NA         NA 0.05263158 0.2262943 0.02631579 0.1622214
    ## 3    0.2886992          NA         NA 0.34210526 0.4807829 0.18421053 0.3928595
    ## 4    0.2696038          NA         NA 0.23684211 0.4308515 0.02631579 0.1622214
    ## 5    0.2161575          NA         NA 0.15789474 0.3695370 0.05263158 0.2262943
    ## 6    0.2887435          NA         NA 0.26315789 0.4462583 0.13157895 0.3425700
    ## 7    0.2524785          NA         NA 0.05263158 0.2262943 0.00000000 0.0000000
    ## 8    0.2599695          NA         NA 0.36842105 0.4888515 0.18421053 0.3928595
    ## 9    0.2345275          NA         NA 0.05263158 0.2262943 0.00000000 0.0000000
    ## 10   0.1830783          NA         NA 0.00000000 0.0000000 0.00000000 0.0000000
    ## 11   0.1959164          NA         NA 0.21052632 0.4131550 0.02631579 0.1622214
    ## 12   0.2652986          NA         NA 0.07894737 0.2732763 0.02631579 0.1622214
    ## 13   0.2489931          NA         NA 0.39473684 0.4953554 0.23684211 0.4308515
    ## 14   0.2559858          NA         NA 0.44736842 0.5038966 0.05263158 0.2262943
    ## 15   0.2637390          NA         NA 0.26315789 0.4462583 0.05263158 0.2262943
    ## 16   0.2424931          NA         NA 0.07894737 0.2732763 0.02631579 0.1622214
    ## 17   0.2256563          NA         NA 0.42105263 0.5003555 0.07894737 0.2732763
    ## 18   0.3094120          NA         NA 0.21052632 0.4131550 0.21052632 0.4131550
    ## 19   0.2844698          NA         NA 0.34210526 0.4807829 0.13157895 0.3425700
    ## 20   0.2851064          NA         NA 0.21052632 0.4131550 0.10526316 0.3110117
    ## 21   0.2802910          NA         NA 0.34210526 0.4807829 0.28947368 0.4596059
    ## 22   0.2857825          NA         NA 0.15789474 0.3695370 0.07894737 0.2732763
    ## 23   0.2440247          NA         NA 0.50000000 0.5067117 0.07894737 0.2732763
    ## 24   0.2137603          NA         NA 0.02631579 0.1622214 0.02631579 0.1622214
    ## 25   0.2447687          NA         NA 0.39473684 0.4953554 0.21052632 0.4131550
    ## 26   0.2680105          NA         NA 0.28947368 0.4596059 0.02631579 0.1622214
    ## 27   0.2470694          NA         NA 0.23684211 0.4308515 0.02631579 0.1622214
    ## 28   0.2432035          NA         NA 0.10526316 0.3110117 0.07894737 0.2732763
    ## 29   0.3040488          NA         NA 0.42105263 0.5003555 0.28947368 0.4596059
    ## 30   0.2328278          NA         NA 0.10526316 0.3110117 0.02631579 0.1622214
    ## 31   0.2502246          NA         NA 0.02631579 0.1622214 0.00000000 0.0000000
    ## 32   0.2681329          NA         NA 0.15789474 0.3695370 0.07894737 0.2732763
    ## 33   0.2880199          NA         NA 0.31578947 0.4710691 0.00000000 0.0000000
    ## 34   0.3556858          NA         NA 0.28947368 0.4596059 0.13157895 0.3425700
    ## 35   0.2830162          NA         NA 0.10526316 0.3110117 0.05263158 0.2262943
    ## 36   0.3203267          NA         NA 0.21052632 0.4131550 0.02631579 0.1622214
    ## 37   0.2062490          NA         NA 0.65789474 0.4807829 0.28947368 0.4596059
    ## 38   0.3061020          NA         NA 0.21052632 0.4131550 0.13157895 0.3425700
    ## 39   0.2237899          NA         NA 0.13157895 0.3425700 0.00000000 0.0000000
    ## 40   0.1820450          NA         NA 0.47368421 0.5060094 0.26315789 0.4462583
    ## 41   0.3173199          NA         NA 0.21052632 0.4131550 0.13157895 0.3425700
    ## 42   0.2624175          NA         NA 0.15789474 0.3695370 0.10526316 0.3110117
    ## 43   0.3329562          NA         NA 0.26315789 0.4462583 0.07894737 0.2732763
    ## 44   0.1278976          NA         NA 0.97368421 0.1622214 0.15789474 0.3695370
    ## 45   0.2781857          NA         NA 0.44736842 0.5038966 0.23684211 0.4308515
    ## 46   0.2587161          NA         NA 0.18421053 0.3928595 0.05263158 0.2262943
    ## 47   0.2392329          NA         NA 0.21052632 0.4131550 0.07894737 0.2732763
    ## 48   0.3266102          NA         NA 0.21052632 0.4131550 0.07894737 0.2732763
    ## 49   0.2754673          NA         NA 0.13157895 0.3425700 0.02631579 0.1622214
    ## 50   0.2792400          NA         NA 0.18421053 0.3928595 0.05263158 0.2262943
    ## 51   0.2617944          NA         NA 0.15789474 0.3695370 0.00000000 0.0000000
    ## 52   0.2446068          NA         NA 0.52631579 0.5060094 0.02631579 0.1622214
    ## 53   0.2692204          NA         NA 0.39473684 0.4953554 0.10526316 0.3110117
    ## 54   0.2669727          NA         NA 0.47368421 0.5060094 0.02631579 0.1622214
    ## 55   0.1497798          NA         NA 0.05263158 0.2262943 0.00000000 0.0000000
    ## 56   0.2703447          NA         NA 0.10526316 0.3110117 0.02631579 0.1622214
    ## 57   0.3158278          NA         NA 0.23684211 0.4308515 0.10526316 0.3110117
    ## 58   0.2559728          NA         NA 0.07894737 0.2732763 0.00000000 0.0000000
    ## 59   0.1770347          NA         NA 0.07894737 0.2732763 0.00000000 0.0000000
    ## 60   0.2286472          NA         NA 0.02631579 0.1622214 0.00000000 0.0000000
    ## 61   0.2531241          NA         NA 0.26315789 0.4462583 0.13157895 0.3425700
    ## 62   0.2531159          NA         NA 0.39473684 0.4953554 0.07894737 0.2732763
    ## 63   0.2820277          NA         NA 0.31578947 0.4710691 0.10526316 0.3110117
    ## 64   0.3249528          NA         NA 0.18421053 0.3928595 0.02631579 0.1622214
    ## 65   0.2316564          NA         NA 0.36842105 0.4888515 0.10526316 0.3110117
    ## 66   0.3458367          NA         NA 0.26315789 0.4462583 0.13157895 0.3425700
    ## 67   0.3094202          NA         NA 0.18421053 0.3928595 0.00000000 0.0000000
    ## 68   0.2297468          NA         NA 0.42105263 0.5003555 0.02631579 0.1622214
    ## 69   0.2936227          NA         NA 0.42105263 0.5003555 0.05263158 0.2262943
    ## 70   0.2530347          NA         NA 0.39473684 0.4953554 0.05263158 0.2262943
    ## 71   0.2142291          NA         NA 0.42105263 0.5003555 0.28947368 0.4596059
    ## 72   0.3034526          NA         NA 0.47368421 0.5060094 0.26315789 0.4462583
    ## 73   0.2594357          NA         NA 0.18421053 0.3928595 0.07894737 0.2732763
    ## 74   0.2986708          NA         NA 0.28947368 0.4596059 0.15789474 0.3695370
    ## 75   0.3503049          NA         NA 0.18421053 0.3928595 0.07894737 0.2732763
    ## 76   0.2209699          NA         NA 0.10526316 0.3110117 0.00000000 0.0000000
    ## 77   0.1815251          NA         NA 0.02631579 0.1622214 0.02631579 0.1622214
    ## 78   0.2507886          NA         NA 0.07894737 0.2732763 0.00000000 0.0000000
    ## 79   0.2245125          NA         NA 0.02631579 0.1622214 0.00000000 0.0000000
    ## 80   0.2637357          NA         NA 0.23684211 0.4308515 0.23684211 0.4308515
    ## 81   0.2688824          NA         NA 0.18421053 0.3928595 0.02631579 0.1622214
    ## 82   0.1353729          NA         NA 0.39473684 0.4953554 0.02631579 0.1622214
    ## 83   0.2385607          NA         NA 0.18421053 0.3928595 0.05263158 0.2262943
    ## 84   0.2412040          NA         NA 0.47368421 0.5060094 0.07894737 0.2732763
    ## 85   0.2765025          NA         NA 0.44736842 0.5038966 0.15789474 0.3695370
    ## 86   0.3224685          NA         NA 0.21052632 0.4131550 0.07894737 0.2732763
    ## 87   0.2093100          NA         NA 0.10526316 0.3110117 0.00000000 0.0000000
    ## 88   0.2639206          NA         NA 0.10526316 0.3110117 0.00000000 0.0000000
    ## 89   0.2763923          NA         NA 0.39473684 0.4953554 0.07894737 0.2732763
    ## 90   0.2765194          NA         NA 0.52631579 0.5060094 0.39473684 0.4953554
    ## 91   0.2497750          NA         NA 0.42105263 0.5003555 0.15789474 0.3695370
    ## 92   0.2109525          NA         NA 0.02631579 0.1622214 0.02631579 0.1622214
    ## 93   0.2805285          NA         NA 0.47368421 0.5060094 0.07894737 0.2732763
    ## 94   0.2447820          NA         NA 0.42105263 0.5003555 0.13157895 0.3425700
    ## 95   0.2840833          NA         NA 0.13157895 0.3425700 0.00000000 0.0000000
    ## 96   0.2248110          NA         NA 0.02631579 0.1622214 0.02631579 0.1622214
    ## 97   0.1322783          NA         NA 0.34210526 0.4807829 0.02631579 0.1622214
    ## 98   0.1513527          NA         NA 0.00000000 0.0000000 0.00000000 0.0000000
    ## 99   0.2199497          NA         NA 0.26315789 0.4462583 0.10526316 0.3110117
    ## 100  0.2758956          NA         NA 0.63157895 0.4888515 0.36842105 0.4888515
    ##     ncoef
    ## 1       2
    ## 2       4
    ## 3       6
    ## 4       4
    ## 5       4
    ## 6       4
    ## 7       5
    ## 8       5
    ## 9       6
    ## 10      3
    ## 11      4
    ## 12      4
    ## 13      7
    ## 14      6
    ## 15      6
    ## 16      6
    ## 17      4
    ## 18      5
    ## 19      4
    ## 20      7
    ## 21      5
    ## 22      6
    ## 23      6
    ## 24      6
    ## 25      7
    ## 26      4
    ## 27      7
    ## 28      6
    ## 29      6
    ## 30      4
    ## 31      3
    ## 32      5
    ## 33      5
    ## 34      4
    ## 35      2
    ## 36      1
    ## 37      5
    ## 38      8
    ## 39      2
    ## 40      5
    ## 41      6
    ## 42      5
    ## 43      3
    ## 44      7
    ## 45      5
    ## 46      6
    ## 47      5
    ## 48      5
    ## 49      3
    ## 50      5
    ## 51      5
    ## 52      5
    ## 53      5
    ## 54      6
    ## 55      5
    ## 56      4
    ## 57      4
    ## 58      3
    ## 59      6
    ## 60      6
    ## 61      7
    ## 62      4
    ## 63      6
    ## 64      3
    ## 65      4
    ## 66      5
    ## 67      4
    ## 68      5
    ## 69      7
    ## 70      5
    ## 71      7
    ## 72      7
    ## 73      5
    ## 74      3
    ## 75      4
    ## 76      6
    ## 77      3
    ## 78      4
    ## 79      1
    ## 80      7
    ## 81      3
    ## 82      3
    ## 83      7
    ## 84      4
    ## 85      7
    ## 86      5
    ## 87      4
    ## 88      5
    ## 89      5
    ## 90     10
    ## 91      5
    ## 92      9
    ## 93      2
    ## 94      2
    ## 95      2
    ## 96      6
    ## 97      2
    ## 98      6
    ## 99      4
    ## 100     6

``` r
evalplot.nulls(mod.null, stats = "auc.val", plot.type = "histogram")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](Aa_SDM_files/figure-gfm/Null%20model-1.png)<!-- -->

## Visualize

``` r
pal <- c("#ffffff","#dad7cd","#a3b18a","#588157","#3a5a40","#344e41")

ref <- ne_countries(continent = c("North America", "South America"), scale = 50, returnclass = "sf") %>% 
  st_crop(xmin = -170, ymin = 0, xmax = -11, ymax = 83)
```

    ## Warning: attribute variables are assumed to be spatially constant throughout all
    ## geometries

``` r
adm <- ne_states(country = c("United States of America", "Canada"), returnclass = "sf") %>% 
  st_crop(xmin = -102, ymin = 24, xmax = -56, ymax = 52) 
```

    ## Warning: attribute variables are assumed to be spatially constant throughout all
    ## geometries

``` r
bb <- sa.bb %>% 
  st_as_sfc()

refmap <- tm_shape(ref) +
  tm_polygons(alpha = 0.3, col = "#405F63", border.alpha = 0) +
  tm_shape(bb) +
  tm_polygons(alpha = 0, border.col = "#EC4E20", lwd = 3) +
  tm_layout(inner.margins = 0, bg.color = NA, frame = FALSE)



names(pred) <- "Suitability"

predmap <- tm_shape(pred) +
    tm_raster(palette = pal, style = "cont", colorNA = "#E9EBEC", showNA = FALSE, 
            alpha = 1, legend.reverse = TRUE, breaks = c(0, 0.25, 0.50, 0.75, 1)) +
  tm_shape(adm) +
    tm_polygons(alpha = 0.1, border.col = "#405F63", border.alpha = 0.2, col = "#EAEAEC") +
  tm_compass(type = "arrow", show.labels = 0, color.light = "#000000", position = c(0.875, 0.1)) +
  tm_scale_bar(breaks = c(0, 200), color.light = "#000000", text.size = 0.7,
               position = c(0.872, 0.01)) +
  tm_layout(inner.margins = 0, bg.color = "#ffffff", frame = "#D3D7D9",
            legend.outside = T, legend.outside.position = "right")


predmap
print(refmap, vp = grid::viewport(0.54, 0.243, width = 0.25, height = 0.3))
```

![](Aa_SDM_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ext <- extent(pred)
asp <- (ext[4] - ext[3])/(abs(ext[1]) - abs(ext[2]))
vp <- grid::viewport(0.54, 0.225, width = 0.25, height = 0.3)

tmap_save(predmap, filename = "../pred/dist/Aa_figure.png",
          dpi = 600, insets_tm = refmap, insets_vp = vp,
          height = asp*200, width = 200, units = "mm")
```

    ## Map saved to Z:\Github\Andrena_asteris_SDM\pred\dist\Aa_figure.png

    ## Resolution: 4724.409 by 2724.562 pixels

    ## Size: 7.874016 by 4.540937 inches (600 dpi)

## Session Info

``` r
sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 22000)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.1252 
    ## [2] LC_CTYPE=English_United States.1252   
    ## [3] LC_MONETARY=English_United States.1252
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.1252    
    ## 
    ## attached base packages:
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] rnaturalearthdata_0.2.0  rnaturalearthhires_0.2.0 seegSDM_0.1-9           
    ##  [4] snowfall_1.84-6.1        snow_0.4-4               rgeos_0.5-9             
    ##  [7] rnaturalearth_0.1.0      tmap_3.3-3               rasterVis_0.51.2        
    ## [10] lattice_0.20-45          rJava_1.0-6              ENMeval_2.0.3           
    ## [13] magrittr_2.0.2           rgdal_1.5-28             RStoolbox_0.2.6         
    ## [16] sf_1.0-6                 raster_3.5-15            sp_1.4-6                
    ## [19] curl_4.3.2               CoordinateCleaner_2.0-20 spThin_0.2.0            
    ## [22] knitr_1.37               fields_13.3              viridis_0.6.2           
    ## [25] viridisLite_0.4.0        spam_2.8-0               ggmap_3.0.0             
    ## [28] lubridate_1.8.0          forcats_0.5.1            stringr_1.4.0           
    ## [31] dplyr_1.0.8              purrr_0.3.4              readr_2.1.2             
    ## [34] tidyr_1.2.0              tibble_3.1.6             ggplot2_3.3.5           
    ## [37] tidyverse_1.3.1          pacman_0.5.1            
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] utf8_1.2.2               tidyselect_1.1.2         rgbif_3.7.0             
    ##   [4] htmlwidgets_1.5.4        pROC_1.18.0              munsell_0.5.0           
    ##   [7] codetools_0.2-18         units_0.8-0              future_1.24.0           
    ##  [10] withr_2.5.0              colorspace_2.0-3         highr_0.9               
    ##  [13] uuid_1.0-3               rstudioapi_0.13          stats4_4.1.2            
    ##  [16] wk_0.6.0                 listenv_0.8.0            labeling_0.4.2          
    ##  [19] conditionz_0.1.0         RgoogleMaps_1.4.5.3      oai_0.3.2               
    ##  [22] bit64_4.0.5              farver_2.1.0             parallelly_1.30.0       
    ##  [25] vctrs_0.3.8              generics_0.1.2           ipred_0.9-12            
    ##  [28] xfun_0.29                R6_2.5.1                 doParallel_1.0.17       
    ##  [31] bitops_1.0-7             assertthat_0.2.1         scales_1.1.1            
    ##  [34] vroom_1.5.7              nnet_7.3-16              gtable_0.3.0            
    ##  [37] wellknown_0.7.4          globals_0.14.0           lwgeom_0.2-8            
    ##  [40] timeDate_3043.102        rlang_1.0.1              splines_4.1.2           
    ##  [43] lazyeval_0.2.2           ModelMetrics_1.2.2.2     dichromat_2.0-0         
    ##  [46] hexbin_1.28.2            broom_0.7.12             s2_1.0.7                
    ##  [49] yaml_2.3.5               reshape2_1.4.4           abind_1.4-5             
    ##  [52] modelr_0.1.8             crosstalk_1.2.0          backports_1.4.1         
    ##  [55] caret_6.0-90             tools_4.1.2              lava_1.6.10             
    ##  [58] ellipsis_0.3.2           RColorBrewer_1.1-2       proxy_0.4-26            
    ##  [61] Rcpp_1.0.8               plyr_1.8.6               base64enc_0.1-3         
    ##  [64] classInt_0.4-3           rpart_4.1-15             tmaptools_3.1-1         
    ##  [67] zoo_1.8-9                haven_2.4.3              fs_1.5.2                
    ##  [70] leafem_0.1.6             data.table_1.14.2        reprex_2.0.1            
    ##  [73] whisker_0.4              hms_1.1.1                evaluate_0.15           
    ##  [76] XML_3.99-0.9             leaflet_2.1.0            jpeg_0.1-9              
    ##  [79] readxl_1.3.1             gridExtra_2.3            compiler_4.1.2          
    ##  [82] maps_3.4.0               KernSmooth_2.23-20       crayon_1.5.0            
    ##  [85] rangeModelMetadata_0.1.4 htmltools_0.5.2          tzdb_0.2.0              
    ##  [88] DBI_1.1.2                dbplyr_2.1.1             MASS_7.3-54             
    ##  [91] Matrix_1.3-4             cli_3.2.0                parallel_4.1.2          
    ##  [94] dotCall64_1.0-1          gower_1.0.0              pkgconfig_2.0.3         
    ##  [97] geosphere_1.5-14         terra_1.5-21             recipes_0.2.0           
    ## [100] xml2_1.3.3               foreach_1.5.2            hardhat_0.2.0           
    ## [103] prodlim_2019.11.13       rvest_1.0.2              digest_0.6.29           
    ## [106] rmarkdown_2.12           cellranger_1.1.0         leafsync_0.1.0          
    ## [109] PresenceAbsence_1.1.10   gtools_3.9.2             rjson_0.2.21            
    ## [112] lifecycle_1.0.1          nlme_3.1-153             dismo_1.3-5             
    ## [115] jsonlite_1.8.0           fansi_1.0.2              pillar_1.7.0            
    ## [118] fastmap_1.1.0            httr_1.4.2               survival_3.2-13         
    ## [121] glue_1.6.2               gbm_2.1.8                png_0.1-7               
    ## [124] iterators_1.0.14         bit_4.0.4                class_7.3-19            
    ## [127] stringi_1.7.6            doSNOW_1.0.20            latticeExtra_0.6-29     
    ## [130] stars_0.5-5              e1071_1.7-9              future.apply_1.8.1
