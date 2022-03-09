Andrena asteris
================
Mark Buckner
2022-03-09

-   [Load Pkgs and functions](#load-pkgs-and-functions)
    -   [Pkgs](#pkgs)
    -   [Functions](#functions)
-   [Load and process occ Data](#load-and-process-occ-data)
    -   [Load and filter by date](#load-and-filter-by-date)
    -   [Check for georeferencing
        errors](#check-for-georeferencing-errors)
    -   [Spatial Thinning w/ spThin](#spatial-thinning-w-spthin)
-   [Environmental Covariates](#environmental-covariates)
    -   [Define study area](#define-study-area)
    -   [Download and process env](#download-and-process-env)
    -   [PCA](#pca)
    -   [Extract env. and MESS](#extract-env-and-mess)
-   [SDM Model](#sdm-model)
    -   [Assign random background
        points](#assign-random-background-points)
    -   [Partition for Model
        Evaluation](#partition-for-model-evaluation)
    -   [Tune Model](#tune-model)
    -   [Model Selection](#model-selection)
-   [Predictions](#predictions)
    -   [Threshold](#threshold)
    -   [ENM Null](#enm-null)
-   [Visualize](#visualize)
-   [Session Info](#session-info)

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

date_start <- as.Date("1981-01-01")
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
    ## Source: "C:\Users\mabuc\AppData\Local\Temp\RtmpGikklm", layer: "ne_50m_land"
    ## with 1420 features
    ## It has 3 fields
    ## Integer64 fields read as strings:  scalerank

    ## Flagged 17 records.

    ## Testing geographic outliers

    ## Flagged 0 records.

    ## Testing GBIF headquarters, flagging records around Copenhagen

    ## Flagged 0 records.

    ## Testing biodiversity institutions

    ## Flagged 2 records.

    ## Flagged 19 of 69 records, EQ = 0.28.

``` r
summary(flags)
```

    ##     .val     .equ     .zer     .cap     .cen     .sea     .otl     .gbf 
    ##        0        0        0        0        0       17        0        0 
    ##    .inst .summary 
    ##        2       19

``` r
plot(flags, lon = "lon", lat = "lat")
```

![](Aa_SDM_files/figure-gfm/Clean%20coordinates-1.png)<!-- -->

``` r
occ.flagged <- occ[!flags$.summary,]

write_csv(occ.flagged, file = "../occ/Aa_flagged.csv")
```

17 occurrences were flagged for potentially being in the sea. These
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
    ##  Script Started at: Wed Mar 09 12:51:54 2022
    ## lat.long.thin.count
    ##  42 
    ## 100 
    ## [1] "Maximum number of records after thinning: 42"
    ## [1] "Number of data.frames with max records: 100"
    ## [1] "Writing new *.csv files"

    ## Warning in thin(loc.data = occs, lat.col = "lat", long.col = "lon", spec.col =
    ## "spp", : Created new output directory: ../occ/Aa_thinned_full/

    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin1.csv"
    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin2.csv"
    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin3.csv"
    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin4.csv"
    ## [1] "Writing file: ../occ/Aa_thinned_full/Aa_thinned_thin5.csv"

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

Retained 42 after spatial thinning.

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

n = 31.

``` r
part <- get.jackknife(occp, bg)
evalplot.grps(pts = occp, pts.grp = part$occs.grp, envs = env)
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

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

    ## ENMevaluate completed in 12 minutes 9.4 seconds.

``` r
e.mx
```

    ## An object of class:  ENMevaluation 
    ##  occurrence/background points:  31 / 24664 
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
| L   | 1   | fc.L_rm.1   | 0.8305275 | NA        |    0.1884625 |   1.0532924 |   0.8069473 |  1.4227606 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  666.2063 |  28.6753851 | 0.0000002 |     7 |
| Q   | 1   | fc.Q_rm.1   | 0.8263775 | NA        |    0.1249702 |   0.5140650 |   0.7965927 |  0.8156526 |          NA |         NA |  0.1612903 |  2.014515 |  0.0967742 | 1.6193420 |  669.7335 |  32.2026015 | 0.0000000 |     6 |
| LQ  | 1   | fc.LQ_rm.1  | 0.9108129 | NA        |    0.1255044 |   0.7379657 |   0.8727844 |  0.9602475 |          NA |         NA |  0.1290323 |  1.836161 |  0.0967742 | 1.6193420 |  647.6438 |  10.1129153 | 0.0018648 |    13 |
| LQH | 1   | fc.LQH_rm.1 | 0.9142462 | NA        |    0.1215326 |   0.7380632 |   0.8729361 |  0.9417823 |          NA |         NA |  0.1290323 |  1.836161 |  0.0967742 | 1.6193420 |        NA |          NA |        NA |    35 |
| H   | 1   | fc.H_rm.1   | 0.9106566 | NA        |    0.1288751 |   0.7875312 |   0.8753544 |  1.0157591 |          NA |         NA |  0.1290323 |  1.836161 |  0.0967742 | 1.6193420 |        NA |          NA |        NA |    34 |
| L   | 2   | fc.L_rm.2   | 0.8312285 | NA        |    0.1865278 |   1.0612227 |   0.8067812 |  1.4193852 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  667.6615 |  30.1305628 | 0.0000001 |     7 |
| Q   | 2   | fc.Q_rm.2   | 0.8207273 | NA        |    0.1248953 |   0.5013026 |   0.7898583 |  0.8047486 |          NA |         NA |  0.1612903 |  2.014515 |  0.0645161 | 1.3455906 |  667.9824 |  30.4514879 | 0.0000001 |     5 |
| LQ  | 2   | fc.LQ_rm.2  | 0.9075000 | NA        |    0.1197332 |   0.7178189 |   0.8780696 |  0.9299334 |          NA |         NA |  0.1290323 |  1.836161 |  0.0967742 | 1.6193420 |  637.5329 |   0.0019411 | 0.2925537 |    10 |
| LQH | 2   | fc.LQH_rm.2 | 0.9075170 | NA        |    0.1194056 |   0.7190189 |   0.8774078 |  0.9291183 |          NA |         NA |  0.1290323 |  1.836161 |  0.0967742 | 1.6193420 |  637.5309 |   0.0000000 | 0.2928378 |    10 |
| H   | 2   | fc.H_rm.2   | 0.9042761 | NA        |    0.1338717 |   0.7662581 |   0.8759731 |  1.0249937 |          NA |         NA |  0.1612903 |  2.014515 |  0.0967742 | 1.6193420 | 1477.6630 | 840.1320348 | 0.0000000 |    28 |
| L   | 3   | fc.L_rm.3   | 0.8307550 | NA        |    0.1850816 |   1.0502066 |   0.8079060 |  1.4064172 |          NA |         NA |  0.1290323 |  1.836161 |  0.0322581 | 0.9677419 |  665.7786 |  28.2476956 | 0.0000002 |     6 |
| Q   | 3   | fc.Q_rm.3   | 0.8151322 | NA        |    0.1241723 |   0.4855546 |   0.7878245 |  0.7943491 |          NA |         NA |  0.1935484 |  2.163937 |  0.0322581 | 0.9677419 |  666.7142 |  29.1832210 | 0.0000001 |     4 |
| LQ  | 3   | fc.LQ_rm.3  | 0.9060666 | NA        |    0.1185052 |   0.7036650 |   0.8788269 |  0.9168804 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  638.6781 |   1.1472161 | 0.1650108 |     9 |
| LQH | 3   | fc.LQH_rm.3 | 0.9060640 | NA        |    0.1184925 |   0.7035506 |   0.8788413 |  0.9167628 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  638.6816 |   1.1506930 | 0.1647242 |     9 |
| H   | 3   | fc.H_rm.3   | 0.8954333 | NA        |    0.1384746 |   0.7852651 |   0.8685468 |  1.0550269 |          NA |         NA |  0.1612903 |  2.014515 |  0.0322581 | 0.9677419 |        NA |          NA |        NA |    63 |
| L   | 4   | fc.L_rm.4   | 0.8300540 | NA        |    0.1841075 |   1.0309542 |   0.8074835 |  1.3889395 |          NA |         NA |  0.1290323 |  1.836161 |  0.0322581 | 0.9677419 |  667.5083 |  29.9773848 | 0.0000001 |     6 |
| Q   | 4   | fc.Q_rm.4   | 0.8107167 | NA        |    0.1254770 |   0.4749307 |   0.7807082 |  0.7924748 |          NA |         NA |  0.2258065 |  2.290095 |  0.0322581 | 0.9677419 |  668.1473 |  30.6163615 | 0.0000001 |     4 |
| LQ  | 4   | fc.LQ_rm.4  | 0.9043323 | NA        |    0.1179340 |   0.6877551 |   0.8786268 |  0.9041840 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  644.4165 |   6.8856120 | 0.0093634 |     9 |
| LQH | 4   | fc.LQH_rm.4 | 0.9043231 | NA        |    0.1179392 |   0.6877757 |   0.8786242 |  0.9042183 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  644.4211 |   6.8901229 | 0.0093423 |     9 |
| H   | 4   | fc.H_rm.4   | 0.8844509 | NA        |    0.1442853 |   0.8177661 |   0.8563101 |  1.1000911 |          NA |         NA |  0.1290323 |  1.836161 |  0.0322581 | 0.9677419 |        NA |          NA |        NA |    46 |
| L   | 5   | fc.L_rm.5   | 0.8292863 | NA        |    0.1829738 |   1.0114504 |   0.8066949 |  1.3710579 |          NA |         NA |  0.1290323 |  1.836161 |  0.0322581 | 0.9677419 |  669.5577 |  32.0267637 | 0.0000000 |     6 |
| Q   | 5   | fc.Q_rm.5   | 0.8032224 | NA        |    0.1274915 |   0.4529281 |   0.7789438 |  0.7927073 |          NA |         NA |  0.2580645 |  2.396668 |  0.0322581 | 0.9677419 |  672.9879 |  35.4569296 | 0.0000000 |     5 |
| LQ  | 5   | fc.LQ_rm.5  | 0.9013895 | NA        |    0.1170262 |   0.6690236 |   0.8784332 |  0.8888014 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  641.9493 |   4.4183780 | 0.0321506 |     7 |
| LQH | 5   | fc.LQH_rm.5 | 0.9013895 | NA        |    0.1170274 |   0.6690352 |   0.8784319 |  0.8888130 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 | 1.3455906 |  641.9493 |   4.4183241 | 0.0321514 |     7 |
| H   | 5   | fc.H_rm.5   | 0.8729858 | NA        |    0.1660713 |   0.8725472 |   0.8312756 |  1.2171696 |          NA |         NA |  0.1290323 |  1.836161 |  0.0322581 | 0.9677419 |        NA |          NA |        NA |    37 |

``` r
opt.auc <- res %>% 
  slice_max(auc.val.avg)
kable(opt.auc)
```

| fc  | rm  | tune.args   | auc.train | cbi.train | auc.diff.avg | auc.diff.sd | auc.val.avg | auc.val.sd | cbi.val.avg | cbi.val.sd | or.10p.avg | or.10p.sd | or.mtp.avg | or.mtp.sd |     AICc | delta.AICc |     w.AIC | ncoef |
|:----|:----|:------------|----------:|:----------|-------------:|------------:|------------:|-----------:|------------:|-----------:|-----------:|----------:|-----------:|----------:|---------:|-----------:|----------:|------:|
| LQH | 3   | fc.LQH_rm.3 |  0.906064 | NA        |    0.1184925 |   0.7035506 |   0.8788413 |  0.9167628 |          NA |         NA |  0.1290323 |  1.836161 |  0.0645161 |  1.345591 | 638.6816 |   1.150693 | 0.1647242 |     9 |

``` r
dismo::response(eval.models(e.mx)[[opt.auc[1,]$tune.args]])
```

![](Aa_SDM_files/figure-gfm/Model%20selection-2.png)<!-- -->

``` r
plot(eval.models(e.mx)[[opt.auc[1,]$tune.args]])
```

![](Aa_SDM_files/figure-gfm/Model%20selection-3.png)<!-- -->

## Predictions

``` r
dir.create("../pred/dist")
pred <- eval.predictions(e.mx)[[opt.auc[1,]$tune.args]]
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
mod.null <- ENMnulls(e.mx, mod.settings = list(fc = as.character(opt.auc[1,]$fc), rm = as.numeric(opt.auc[1,]$rm)), no.iter = 100, parallel = TRUE, quiet = TRUE)
null.results(mod.null)
```

    ##      fc rm   tune.args auc.train cbi.train auc.diff.avg auc.diff.sd auc.val.avg
    ## 1   LQH  3 fc.LQH_rm.3 0.6866754        NA    0.2553649  0.12368041   0.4385757
    ## 2   LQH  3 fc.LQH_rm.3 0.6535665        NA    0.3532272  0.19828341   0.3592490
    ## 3   LQH  3 fc.LQH_rm.3 0.7850289        NA    0.5189130  0.23977874   0.2666038
    ## 4   LQH  3 fc.LQH_rm.3 0.6560398        NA    0.2655061  0.17260451   0.4258584
    ## 5   LQH  3 fc.LQH_rm.3 0.5532492        NA    0.2701460  0.17108906   0.4937705
    ## 6   LQH  3 fc.LQH_rm.3 0.6416528        NA    0.1432912  0.12217917   0.6199345
    ## 7   LQH  3 fc.LQH_rm.3 0.6287478        NA    0.2450802  0.18065630   0.4004452
    ## 8   LQH  3 fc.LQH_rm.3 0.6863739        NA    0.3160379  0.21709667   0.4526069
    ## 9   LQH  3 fc.LQH_rm.3 0.6565433        NA    0.2097405  0.14372447   0.6309523
    ## 10  LQH  3 fc.LQH_rm.3 0.6479962        NA    0.2764569  0.14082206   0.4220086
    ## 11  LQH  3 fc.LQH_rm.3 0.7308968        NA    0.2312305  0.18719745   0.6257050
    ## 12  LQH  3 fc.LQH_rm.3 0.6844068        NA    0.1634452  0.09974911   0.7617796
    ## 13  LQH  3 fc.LQH_rm.3 0.6480432        NA    0.3419601  0.18590174   0.3667328
    ## 14  LQH  3 fc.LQH_rm.3 0.5847121        NA    0.2138807  0.16241027   0.6073930
    ## 15  LQH  3 fc.LQH_rm.3 0.6187136        NA    0.1160818  0.11268603   0.5501881
    ## 16  LQH  3 fc.LQH_rm.3 0.6350414        NA    0.2849097  0.14360044   0.6489280
    ## 17  LQH  3 fc.LQH_rm.3 0.6644326        NA    0.1756319  0.13374301   0.6719151
    ## 18  LQH  3 fc.LQH_rm.3 0.6676918        NA    0.3613072  0.18812794   0.4065348
    ## 19  LQH  3 fc.LQH_rm.3 0.7018281        NA    0.3437777  0.19890082   0.4317851
    ## 20  LQH  3 fc.LQH_rm.3 0.6395301        NA    0.2105711  0.15480217   0.5828777
    ## 21  LQH  3 fc.LQH_rm.3 0.6454314        NA    0.2854699  0.19933439   0.4235048
    ## 22  LQH  3 fc.LQH_rm.3 0.6295188        NA    0.2863686  0.16380782   0.4724772
    ## 23  LQH  3 fc.LQH_rm.3 0.7045341        NA    0.3146656  0.17751348   0.4087608
    ## 24  LQH  3 fc.LQH_rm.3 0.6723951        NA    0.3088576  0.18978850   0.4603235
    ## 25  LQH  3 fc.LQH_rm.3 0.6807290        NA    0.3147269  0.22485639   0.4236016
    ## 26  LQH  3 fc.LQH_rm.3 0.6703822        NA    0.3884767  0.20215794   0.3267555
    ## 27  LQH  3 fc.LQH_rm.3 0.6568807        NA    0.2481228  0.14055639   0.4788159
    ## 28  LQH  3 fc.LQH_rm.3 0.7064659        NA    0.2940260  0.16476704   0.5479646
    ## 29  LQH  3 fc.LQH_rm.3 0.6354017        NA    0.3181959  0.18686991   0.5358221
    ## 30  LQH  3 fc.LQH_rm.3 0.6411179        NA    0.3015158  0.19669062   0.3706558
    ## 31  LQH  3 fc.LQH_rm.3 0.6991155        NA    0.2051148  0.09913567   0.8044139
    ## 32  LQH  3 fc.LQH_rm.3 0.6585058        NA    0.3953545  0.18450991   0.2883333
    ## 33  LQH  3 fc.LQH_rm.3 0.6038597        NA    0.2903304  0.10449706   0.7772671
    ## 34  LQH  3 fc.LQH_rm.3 0.6441902        NA    0.3249425  0.15737626   0.4526618
    ## 35  LQH  3 fc.LQH_rm.3 0.6263491        NA    0.1977765  0.16407647   0.5319520
    ## 36  LQH  3 fc.LQH_rm.3 0.6302336        NA    0.3378611  0.18238776   0.5247376
    ## 37  LQH  3 fc.LQH_rm.3 0.6098251        NA    0.2420630  0.16613868   0.4463296
    ## 38  LQH  3 fc.LQH_rm.3 0.5012183        NA    0.2618720  0.15737480   0.5616976
    ## 39  LQH  3 fc.LQH_rm.3 0.6530290        NA    0.2415691  0.17681761   0.6352932
    ## 40  LQH  3 fc.LQH_rm.3 0.6060119        NA    0.3080292  0.14372739   0.4031624
    ## 41  LQH  3 fc.LQH_rm.3 0.6706470        NA    0.2503845  0.12159849   0.5755731
    ## 42  LQH  3 fc.LQH_rm.3 0.6790261        NA    0.2358597  0.14427227   0.5819426
    ## 43  LQH  3 fc.LQH_rm.3 0.6262902        NA    0.1390239  0.12925665   0.6163012
    ## 44  LQH  3 fc.LQH_rm.3 0.6591499        NA    0.2854728  0.17074250   0.5366565
    ## 45  LQH  3 fc.LQH_rm.3 0.6806950        NA    0.2950692  0.19357938   0.4858760
    ## 46  LQH  3 fc.LQH_rm.3 0.7147226        NA    0.2541091  0.20892898   0.5678853
    ## 47  LQH  3 fc.LQH_rm.3 0.6509049        NA    0.4351221  0.19859885   0.2565742
    ## 48  LQH  3 fc.LQH_rm.3 0.6764182        NA    0.2517119  0.16434459   0.5370084
    ## 49  LQH  3 fc.LQH_rm.3 0.6599399        NA    0.3953013  0.16900112   0.4050399
    ## 50  LQH  3 fc.LQH_rm.3 0.6672079        NA    0.3997717  0.18108409   0.2928991
    ## 51  LQH  3 fc.LQH_rm.3 0.6495552        NA    0.2429678  0.16885453   0.5839947
    ## 52  LQH  3 fc.LQH_rm.3 0.6952925        NA    0.2067056  0.11746796   0.6275413
    ## 53  LQH  3 fc.LQH_rm.3 0.6799377        NA    0.2636721  0.16568150   0.6735342
    ## 54  LQH  3 fc.LQH_rm.3 0.6887046        NA    0.2187524  0.17333895   0.5325537
    ## 55  LQH  3 fc.LQH_rm.3 0.6862313        NA    0.3927720  0.19166228   0.3156880
    ## 56  LQH  3 fc.LQH_rm.3 0.6273784        NA    0.2970127  0.18270939   0.5049962
    ## 57  LQH  3 fc.LQH_rm.3 0.7106603        NA    0.3584956  0.22886234   0.3750536
    ## 58  LQH  3 fc.LQH_rm.3 0.6314381        NA    0.2913297  0.16424403   0.4013882
    ## 59  LQH  3 fc.LQH_rm.3 0.6118746        NA    0.2915611  0.17982749   0.3804252
    ## 60  LQH  3 fc.LQH_rm.3 0.6330187        NA    0.3309709  0.18003182   0.5143084
    ## 61  LQH  3 fc.LQH_rm.3 0.7013899        NA    0.3932059  0.23740634   0.3305235
    ## 62  LQH  3 fc.LQH_rm.3 0.6595574        NA    0.2694291  0.13950442   0.5283514
    ## 63  LQH  3 fc.LQH_rm.3 0.7407826        NA    0.5418149  0.20720639   0.2072251
    ## 64  LQH  3 fc.LQH_rm.3 0.6219624        NA    0.2509791  0.13132370   0.7575754
    ## 65  LQH  3 fc.LQH_rm.3 0.6398414        NA    0.3395146  0.19039614   0.3312612
    ## 66  LQH  3 fc.LQH_rm.3 0.6060871        NA    0.1898853  0.11635949   0.4492247
    ## 67  LQH  3 fc.LQH_rm.3 0.6487090        NA    0.4067766  0.19185687   0.2598132
    ## 68  LQH  3 fc.LQH_rm.3 0.6148867        NA    0.3116381  0.15461270   0.5180118
    ## 69  LQH  3 fc.LQH_rm.3 0.6372897        NA    0.3299778  0.15475445   0.4002608
    ## 70  LQH  3 fc.LQH_rm.3 0.5992030        NA    0.1966010  0.14080830   0.4918452
    ## 71  LQH  3 fc.LQH_rm.3 0.6524221        NA    0.4212669  0.18554055   0.2504185
    ## 72  LQH  3 fc.LQH_rm.3 0.6205970        NA    0.2059023  0.13444392   0.4652504
    ## 73  LQH  3 fc.LQH_rm.3 0.6384550        NA    0.2524364  0.19642737   0.4880222
    ## 74  LQH  3 fc.LQH_rm.3 0.5806347        NA    0.2798934  0.16302174   0.5356821
    ## 75  LQH  3 fc.LQH_rm.3 0.7300630        NA    0.3503514  0.18725188   0.4315850
    ## 76  LQH  3 fc.LQH_rm.3 0.6643730        NA    0.3099728  0.18397346   0.4194242
    ## 77  LQH  3 fc.LQH_rm.3 0.6366501        NA    0.3331431  0.18112961   0.3202042
    ## 78  LQH  3 fc.LQH_rm.3 0.6653625        NA    0.2089020  0.18281146   0.6218944
    ## 79  LQH  3 fc.LQH_rm.3 0.6962852        NA    0.4032466  0.16955218   0.3129650
    ## 80  LQH  3 fc.LQH_rm.3 0.6296018        NA    0.3137377  0.15927736   0.4674725
    ## 81  LQH  3 fc.LQH_rm.3 0.6966782        NA    0.1967630  0.14395885   0.6993018
    ## 82  LQH  3 fc.LQH_rm.3 0.6057797        NA    0.1892193  0.12453075   0.5741881
    ## 83  LQH  3 fc.LQH_rm.3 0.6694693        NA    0.3434241  0.18843818   0.3711378
    ## 84  LQH  3 fc.LQH_rm.3 0.5826037        NA    0.2868463  0.16356795   0.4277385
    ## 85  LQH  3 fc.LQH_rm.3 0.6819767        NA    0.3290030  0.18734476   0.4331702
    ## 86  LQH  3 fc.LQH_rm.3 0.6818119        NA    0.3177579  0.19222162   0.5417410
    ## 87  LQH  3 fc.LQH_rm.3 0.6196814        NA    0.2737937  0.17668449   0.4297513
    ## 88  LQH  3 fc.LQH_rm.3 0.6950806        NA    0.3060388  0.22751798   0.5346241
    ## 89  LQH  3 fc.LQH_rm.3 0.6967364        NA    0.3044614  0.20835337   0.4685934
    ## 90  LQH  3 fc.LQH_rm.3 0.6023995        NA    0.2139143  0.12880580   0.6799626
    ## 91  LQH  3 fc.LQH_rm.3 0.7005450        NA    0.4501312  0.18912866   0.2630528
    ## 92  LQH  3 fc.LQH_rm.3 0.6629285        NA    0.2752466  0.20198445   0.4453690
    ## 93  LQH  3 fc.LQH_rm.3 0.6265047        NA    0.3021282  0.19397377   0.4141455
    ## 94  LQH  3 fc.LQH_rm.3 0.6916500        NA    0.3314668  0.21594724   0.4333481
    ## 95  LQH  3 fc.LQH_rm.3 0.6011393        NA    0.2632110  0.17112516   0.6115587
    ## 96  LQH  3 fc.LQH_rm.3 0.6034661        NA    0.2722118  0.12995755   0.6277244
    ## 97  LQH  3 fc.LQH_rm.3 0.6604382        NA    0.2487374  0.19035927   0.4446503
    ## 98  LQH  3 fc.LQH_rm.3 0.6709158        NA    0.3831328  0.19142998   0.3280372
    ## 99  LQH  3 fc.LQH_rm.3 0.6598497        NA    0.3175997  0.18079343   0.4214259
    ## 100 LQH  3 fc.LQH_rm.3 0.6570298        NA    0.2995081  0.18145911   0.5334443
    ##     auc.val.sd cbi.val.avg cbi.val.sd or.10p.avg or.10p.sd or.mtp.avg or.mtp.sd
    ## 1    0.1318823          NA         NA 0.41935484 0.5016103 0.00000000 0.0000000
    ## 2    0.2768927          NA         NA 0.45161290 0.5058794 0.35483871 0.4863735
    ## 3    0.2413947          NA         NA 0.80645161 0.4016097 0.54838710 0.5058794
    ## 4    0.2159477          NA         NA 0.22580645 0.4250237 0.00000000 0.0000000
    ## 5    0.3078102          NA         NA 0.22580645 0.4250237 0.12903226 0.3407771
    ## 6    0.1897287          NA         NA 0.03225806 0.1796053 0.00000000 0.0000000
    ## 7    0.1941266          NA         NA 0.25806452 0.4448027 0.03225806 0.1796053
    ## 8    0.3047543          NA         NA 0.29032258 0.4614144 0.09677419 0.3005372
    ## 9    0.2576995          NA         NA 0.12903226 0.3407771 0.00000000 0.0000000
    ## 10   0.2128321          NA         NA 0.25806452 0.4448027 0.03225806 0.1796053
    ## 11   0.2806836          NA         NA 0.25806452 0.4448027 0.19354839 0.4016097
    ## 12   0.1759217          NA         NA 0.03225806 0.1796053 0.00000000 0.0000000
    ## 13   0.2672514          NA         NA 0.38709677 0.4951376 0.22580645 0.4250237
    ## 14   0.2695894          NA         NA 0.19354839 0.4016097 0.12903226 0.3407771
    ## 15   0.1463521          NA         NA 0.03225806 0.1796053 0.00000000 0.0000000
    ## 16   0.3227464          NA         NA 0.25806452 0.4448027 0.12903226 0.3407771
    ## 17   0.2185822          NA         NA 0.12903226 0.3407771 0.03225806 0.1796053
    ## 18   0.3162766          NA         NA 0.45161290 0.5058794 0.29032258 0.4614144
    ## 19   0.2871921          NA         NA 0.38709677 0.4951376 0.06451613 0.2497310
    ## 20   0.2570579          NA         NA 0.16129032 0.3738783 0.06451613 0.2497310
    ## 21   0.2742724          NA         NA 0.35483871 0.4863735 0.19354839 0.4016097
    ## 22   0.2910378          NA         NA 0.38709677 0.4951376 0.09677419 0.3005372
    ## 23   0.2068296          NA         NA 0.48387097 0.5080005 0.22580645 0.4250237
    ## 24   0.2953143          NA         NA 0.41935484 0.5016103 0.12903226 0.3407771
    ## 25   0.2906432          NA         NA 0.48387097 0.5080005 0.32258065 0.4751910
    ## 26   0.2749793          NA         NA 0.61290323 0.4951376 0.35483871 0.4863735
    ## 27   0.2248800          NA         NA 0.12903226 0.3407771 0.03225806 0.1796053
    ## 28   0.2973525          NA         NA 0.22580645 0.4250237 0.03225806 0.1796053
    ## 29   0.3557105          NA         NA 0.29032258 0.4614144 0.25806452 0.4448027
    ## 30   0.2317068          NA         NA 0.38709677 0.4951376 0.12903226 0.3407771
    ## 31   0.2062701          NA         NA 0.06451613 0.2497310 0.03225806 0.1796053
    ## 32   0.2295153          NA         NA 0.54838710 0.5058794 0.38709677 0.4951376
    ## 33   0.2592164          NA         NA 0.09677419 0.3005372 0.06451613 0.2497310
    ## 34   0.3099324          NA         NA 0.54838710 0.5058794 0.09677419 0.3005372
    ## 35   0.2455718          NA         NA 0.22580645 0.4250237 0.06451613 0.2497310
    ## 36   0.3751302          NA         NA 0.41935484 0.5016103 0.25806452 0.4448027
    ## 37   0.2413487          NA         NA 0.12903226 0.3407771 0.03225806 0.1796053
    ## 38   0.3051967          NA         NA 0.09677419 0.3005372 0.03225806 0.1796053
    ## 39   0.3026261          NA         NA 0.16129032 0.3738783 0.12903226 0.3407771
    ## 40   0.2798751          NA         NA 0.35483871 0.4863735 0.06451613 0.2497310
    ## 41   0.2612135          NA         NA 0.06451613 0.2497310 0.00000000 0.0000000
    ## 42   0.2628626          NA         NA 0.16129032 0.3738783 0.03225806 0.1796053
    ## 43   0.1879859          NA         NA 0.06451613 0.2497310 0.00000000 0.0000000
    ## 44   0.3109050          NA         NA 0.41935484 0.5016103 0.19354839 0.4016097
    ## 45   0.2938277          NA         NA 0.38709677 0.4951376 0.16129032 0.3738783
    ## 46   0.3034080          NA         NA 0.29032258 0.4614144 0.22580645 0.4250237
    ## 47   0.2707722          NA         NA 0.67741935 0.4751910 0.38709677 0.4951376
    ## 48   0.2664258          NA         NA 0.25806452 0.4448027 0.03225806 0.1796053
    ## 49   0.3507740          NA         NA 0.54838710 0.5058794 0.06451613 0.2497310
    ## 50   0.2319926          NA         NA 0.64516129 0.4863735 0.32258065 0.4751910
    ## 51   0.2898546          NA         NA 0.16129032 0.3738783 0.06451613 0.2497310
    ## 52   0.2286428          NA         NA 0.06451613 0.2497310 0.00000000 0.0000000
    ## 53   0.3107756          NA         NA 0.12903226 0.3407771 0.03225806 0.1796053
    ## 54   0.2310052          NA         NA 0.29032258 0.4614144 0.03225806 0.1796053
    ## 55   0.2282406          NA         NA 0.67741935 0.4751910 0.09677419 0.3005372
    ## 56   0.3247739          NA         NA 0.38709677 0.4951376 0.16129032 0.3738783
    ## 57   0.2570585          NA         NA 0.67741935 0.4751910 0.35483871 0.4863735
    ## 58   0.2373048          NA         NA 0.41935484 0.5016103 0.00000000 0.0000000
    ## 59   0.2542266          NA         NA 0.32258065 0.4751910 0.06451613 0.2497310
    ## 60   0.3599388          NA         NA 0.38709677 0.4951376 0.25806452 0.4448027
    ## 61   0.2736257          NA         NA 0.67741935 0.4751910 0.48387097 0.5080005
    ## 62   0.2757072          NA         NA 0.12903226 0.3407771 0.00000000 0.0000000
    ## 63   0.2268158          NA         NA 0.80645161 0.4016097 0.74193548 0.4448027
    ## 64   0.2511243          NA         NA 0.03225806 0.1796053 0.03225806 0.1796053
    ## 65   0.2406761          NA         NA 0.48387097 0.5080005 0.16129032 0.3738783
    ## 66   0.1535382          NA         NA 0.03225806 0.1796053 0.03225806 0.1796053
    ## 67   0.2220700          NA         NA 0.41935484 0.5016103 0.19354839 0.4016097
    ## 68   0.3332383          NA         NA 0.25806452 0.4448027 0.19354839 0.4016097
    ## 69   0.2797403          NA         NA 0.45161290 0.5058794 0.16129032 0.3738783
    ## 70   0.2063225          NA         NA 0.16129032 0.3738783 0.03225806 0.1796053
    ## 71   0.2196752          NA         NA 0.70967742 0.4614144 0.32258065 0.4751910
    ## 72   0.1910446          NA         NA 0.12903226 0.3407771 0.03225806 0.1796053
    ## 73   0.2866838          NA         NA 0.25806452 0.4448027 0.06451613 0.2497310
    ## 74   0.3223361          NA         NA 0.22580645 0.4250237 0.09677419 0.3005372
    ## 75   0.2629303          NA         NA 0.61290323 0.4951376 0.12903226 0.3407771
    ## 76   0.2636532          NA         NA 0.38709677 0.4951376 0.09677419 0.3005372
    ## 77   0.2057992          NA         NA 0.51612903 0.5080005 0.12903226 0.3407771
    ## 78   0.2779464          NA         NA 0.12903226 0.3407771 0.09677419 0.3005372
    ## 79   0.2107297          NA         NA 0.41935484 0.5016103 0.00000000 0.0000000
    ## 80   0.3157784          NA         NA 0.48387097 0.5080005 0.19354839 0.4016097
    ## 81   0.2444305          NA         NA 0.06451613 0.2497310 0.03225806 0.1796053
    ## 82   0.2291756          NA         NA 0.03225806 0.1796053 0.03225806 0.1796053
    ## 83   0.2545177          NA         NA 0.45161290 0.5058794 0.22580645 0.4250237
    ## 84   0.2820007          NA         NA 0.32258065 0.4751910 0.12903226 0.3407771
    ## 85   0.2886234          NA         NA 0.32258065 0.4751910 0.09677419 0.3005372
    ## 86   0.3455864          NA         NA 0.35483871 0.4863735 0.25806452 0.4448027
    ## 87   0.2550801          NA         NA 0.22580645 0.4250237 0.12903226 0.3407771
    ## 88   0.3523978          NA         NA 0.29032258 0.4614144 0.29032258 0.4614144
    ## 89   0.2891177          NA         NA 0.41935484 0.5016103 0.09677419 0.3005372
    ## 90   0.2408582          NA         NA 0.06451613 0.2497310 0.00000000 0.0000000
    ## 91   0.2131166          NA         NA 0.64516129 0.4863735 0.29032258 0.4614144
    ## 92   0.2575544          NA         NA 0.41935484 0.5016103 0.09677419 0.3005372
    ## 93   0.2889192          NA         NA 0.35483871 0.4863735 0.19354839 0.4016097
    ## 94   0.2973505          NA         NA 0.48387097 0.5080005 0.29032258 0.4614144
    ## 95   0.3219906          NA         NA 0.19354839 0.4016097 0.03225806 0.1796053
    ## 96   0.3100007          NA         NA 0.16129032 0.3738783 0.06451613 0.2497310
    ## 97   0.2270166          NA         NA 0.29032258 0.4614144 0.03225806 0.1796053
    ## 98   0.2570247          NA         NA 0.51612903 0.5080005 0.25806452 0.4448027
    ## 99   0.2851812          NA         NA 0.51612903 0.5080005 0.06451613 0.2497310
    ## 100  0.3285413          NA         NA 0.25806452 0.4448027 0.06451613 0.2497310
    ##     ncoef
    ## 1       8
    ## 2       6
    ## 3       9
    ## 4       5
    ## 5       3
    ## 6       7
    ## 7       6
    ## 8      18
    ## 9       9
    ## 10      5
    ## 11      6
    ## 12      6
    ## 13      5
    ## 14      7
    ## 15      4
    ## 16      6
    ## 17      8
    ## 18      8
    ## 19     10
    ## 20      6
    ## 21      8
    ## 22      6
    ## 23      7
    ## 24      7
    ## 25      4
    ## 26      8
    ## 27      7
    ## 28      7
    ## 29      7
    ## 30      8
    ## 31      7
    ## 32      5
    ## 33      4
    ## 34      7
    ## 35      5
    ## 36      7
    ## 37      3
    ## 38      3
    ## 39      6
    ## 40      6
    ## 41      7
    ## 42      9
    ## 43      8
    ## 44      5
    ## 45      9
    ## 46      7
    ## 47      8
    ## 48      6
    ## 49      7
    ## 50      7
    ## 51      4
    ## 52      8
    ## 53      9
    ## 54      5
    ## 55      7
    ## 56      4
    ## 57     10
    ## 58      4
    ## 59      5
    ## 60     10
    ## 61      7
    ## 62      4
    ## 63      5
    ## 64      8
    ## 65      5
    ## 66      5
    ## 67      5
    ## 68      6
    ## 69      8
    ## 70      5
    ## 71      7
    ## 72      4
    ## 73      6
    ## 74      5
    ## 75      9
    ## 76      8
    ## 77      4
    ## 78      4
    ## 79      7
    ## 80      5
    ## 81      6
    ## 82      5
    ## 83      8
    ## 84      6
    ## 85      5
    ## 86      7
    ## 87      9
    ## 88      8
    ## 89      8
    ## 90      5
    ## 91     10
    ## 92      8
    ## 93      7
    ## 94      7
    ## 95      4
    ## 96      6
    ## 97      5
    ## 98      5
    ## 99      7
    ## 100     5

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

occv <- st_as_sf(occp, coords = c("coords.x1", "coords.x2"), crs = proj.crs)
cortland <- matrix(c(42.5563, -76.2282), nrow =1) %>% 
  as.data.frame() %>% 
  st_as_sf(coords = c(2, 1), crs = proj.crs)

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
  tm_shape(occv) +
    tm_dots(shape = 1, size = 0.08, alpha = 0.85)+
  tm_shape(cortland) +
    tm_dots(shape = 19, size = 0.1,alpha = 0.85) +
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
