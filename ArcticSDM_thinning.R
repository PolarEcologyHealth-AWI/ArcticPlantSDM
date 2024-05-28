library(sf)
sf_use_s2(FALSE)
library(stars)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(mgcv)
library(earth)

data     <- "/Volumes/projects/bioing/data/ArcticSDM_data/"
env_data <- "/Volumes/projects/bioing/data/ArcticSDM/environment/CHELSA data 5 km/Modern/"

### Projection
proj <- "+proj=laea +lon_0=-170 +lat_0=90"

# ### Maps
# ecoreg   <- st_read(glue::glue("{data}/Ecoregions/tnc_terr_ecoregions.shp")) %>%
#    filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0) %>%
#    st_union() %>%
#    suppressWarnings() 
#  
# map_wrld <- st_read(glue::glue("{data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
#    dplyr::select(c('name', 'admin', 'geonunit', 'region')) %>% st_intersection(ecoreg) %>% 
#    suppressWarnings()
# 
#  
# ### Subset
# map  <- map_wrld %>%
#    st_transform(proj) %>%
#    mutate(region = if_else(is.na(region), admin, region))
# plot(map %>% dplyr::select(region))
# save(map, file = glue::glue("{data}GBIF/gbif_regions.rda")) 
load(glue::glue("{data}GBIF/gbif_regions.rda"))

### Resolution
res <- 5000

## Thinning
res_thin <- 50000


#########################
## EnvData calibration ##
#########################

env_delete <- c("bio3", "swe", "scd")
env_list   <- tibble(fls = list.files(env_data, pattern = ".tif")) %>%
  mutate(var = sapply(strsplit(fls, "_"), function(x) x[2])) %>%
  filter(!(var%in%env_delete))

envData    <- lapply(env_list$fls, function(x) {
  st <- read_stars(glue::glue("{env_data}{x}"))
}) %>% Reduce("c", .) %>% setNames(env_list$var) %>% merge()
  

##############
## GBIF Occ ##
##############

gbif_grid <- st_read(glue::glue("{data}GBIF/occurance_grid.shp"))

ids <- which(c(gbif_grid %>% st_transform(st_crs(map)) %>%
                 st_intersects(map %>% st_geometry() %>% st_union(), sparse = FALSE)))


### Read in data files
gbif_list <- tibble(fls = list.files(glue::glue("{data}/GBIF/Plantea"))) %>%
  mutate(ids = as.numeric(gsub("_gbif_plantea.csv", "", fls)))
# plot(gbif_grid[gbif_list$ids,])

for(i in 297:nrow(gbif_list)) {
  
  gbif <- tryCatch(read_csv(glue::glue("{data}/GBIF/Plantea/{gbif_list$fls[i]}"), progress = F, show_col_types = FALSE) %>%
                     dplyr::select(phylum, order, family, genus, species, scientificName, year, countryCode, 
                                   decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, basisOfRecord) %>%
                     filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species), 
                            !is.na(coordinateUncertaintyInMeters) | coordinateUncertaintyInMeters < 5000, year > 1970, 
                            basisOfRecord != "FOSSIL_SPECIMEN"), error = function(e) NULL)
  
  if(!is.null(gbif)) {
    
    gbifTab_sf <- gbif %>% filter(!is.na(as.numeric(decimalLongitude))) %>%
      st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>%
      st_transform(st_crs(envData))
    
    
    gbif_thin <- gbifTab_sf %>% group_split(species) %>%
         parallel::mclapply(function(x) { 
              
            if(nrow(x) > 1) {
              distM   <- st_distance(x)
              hc      <- hclust(as.dist(distM), method="complete")
              x$clust <- cutree(hc, h = res_thin)
            } else x$clust <- 1
              
            x[rep(1, length(unique(x$clust))),1:5] %>% st_drop_geometry() %>% bind_cols(
                x %>% st_extract(envData, .) %>% st_as_sf() %>% 
                mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], clust = x$clust) %>%
                st_drop_geometry() %>% group_by(clust) %>% summarise(across(everything(), ~median(.x, na.rm = T)))) %>%
                dplyr::select(-clust)
    
           }, mc.cores = 4) %>% Reduce("rbind", .)
    
    
    # ggplot() +
    #   geom_sf(data = gbif_thin %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(envData)) %>% dplyr::select(species), 
    #           mapping = aes(geometry = geometry, color = as.factor(species)), show.legend = F)
    
    save(gbif_thin, file = glue::glue("{data}GBIF/GBIF_dataComp/ClusterThin/{gsub('.csv', '', gbif_list$fls[i])}.rda"))
       
  } else {
    gbif_thin <- NULL
    save(gbif_thin, file = glue::glue("{data}GBIF/GBIF_dataComp/ClusterThin/{gsub('.csv', '', gbif_list$fls[i])}_error.rda"))
  }
}
