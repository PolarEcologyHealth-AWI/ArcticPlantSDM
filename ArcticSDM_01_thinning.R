library(sf)
sf_use_s2(FALSE)
library(stars)
library(terra)
library(tidyverse)

{
  wd        <- "Arctic_SDM/"
  data      <- "ArcticSDM_data/"
  gbif_data <- "ArcticSDM_data/"
  env_data  <- "CHELSA data 5 km/Modern/"
  wd_out    <- "ArcticSDM/"
}

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"


### GBIF Grid
gbif_fls <- tibble(fls = list.files(glue::glue("{gbif_data}GBIF/Plantea"))) %>%
  mutate(ID = as.numeric(sapply(strsplit(fls, "_"), function(x) x[[1]])), exists = TRUE)

gbif_grid <- st_read(glue::glue("{data}GBIF/occurance_grid.shp")) %>% st_as_sf() %>%
  left_join(gbif_fls %>% dplyr::select(ID, exists), by = join_by(grid_id==ID))

plot(gbif_grid %>% dplyr::select(exists))

for(i in 1:nrow(gbif_fls)){
  
  cat("\b\b\b\b\b\b")
  cat(sprintf("%6d", i))
  flush.console()
  
  
  gbif <- tryCatch(read_csv(glue::glue("{gbif_data}GBIF/Plantea/{gbif_fls$fls[i]}"), progress = F, show_col_types = FALSE) %>%
                     dplyr::select(phylum, order, family, genus, species, scientificName, year, countryCode, 
                                   decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, basisOfRecord) %>%
                     filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species), 
                            !is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 5000, year > 1970, 
                            basisOfRecord != "FOSSIL_SPECIMEN"), error = function(e) NULL) %>% suppressWarnings()
  
  
  if(!is.null(gbif)) {
    if(nrow(gbif)>0) {
      
      gbifTab_sf <- gbif %>% filter(!is.na(as.numeric(decimalLongitude))) %>%
        st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
        st_transform(st_crs(thinnin_grid)) %>% dplyr::select(species, year, countryCode, coordinateUncertaintyInMeters)
     
      if(exists('gbif_all_sf')) {
        gbif_all_sf <- rbind(gbif_all_sf, gbifTab_sf)
      } else gbif_all_sf <- gbifTab_sf
    }
  }
  
}

# save(gbif_all_sf, file = 'gbif_all_sf_29052025.rda')


### Thinning ###
env_list   <- tibble(fls = list.files(env_data, pattern = ".tif")) %>%
  mutate(var = sapply(strsplit(fls, "_"), function(x) x[2]))

thinRast   <- read_stars(glue::glue("{env_data}/{env_list$fls[1]}")) %>%
  setNames('cells') %>% mutate(cells = NA) %>% rast()


for(i in unique(gbif_all_sf$species)) {
  
  cat("\b\b\b\b\b\b")
  cat(sprintf("%6d", which(unique(gbif_all_sf$species)==i)))
  flush.console()
  
  sp_sf_tab <- gbif_all_sf %>% filter(species == i) %>% vect()
  
  occ <- rasterize(
    sp_sf_tab,
    thinRast,
    field = 1,      
    background = NA
  ) %>% st_as_stars() %>% st_as_sf(as_point = TRUE) %>% st_set_crs(proj) %>%
    st_centroid() %>% mutate(species = i) %>% dplyr::select(species) %>% suppressWarnings()
  
  if(nrow(presence_raster)>0) {
    if(exists("OccTab")) {
      OccTab <- OccTab %>% bind_rows(occ)
    } else OccTab <- occ
  }
  
}


### add environmental variables




