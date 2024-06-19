library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(terra)


{
  sdm_wd <- "/Volumes/GeoDatasets/ArcticSDM/SDM_Results/"
  out_wd <- "/Volumes/GeoDatasets/ArcticSDM/Dispersal_Results/"
}

speciesDisp <- read_csv("~/Downloads/species_meters.csv")

spResults <- (tibble(species = list.files(sdm_wd)) %>%
  left_join(tibble(fls = list.files(sdm_wd, pattern = "restricted_50.tif", recursive = T)) %>%
                   mutate(species = sapply(strsplit(fls, "/"), function(x) x[[1]]), predict = TRUE) %>%
                   dplyr::select(-fls), by = "species") %>% filter(!is.na(predict)) %>%
  left_join(speciesDisp, by = "species")) %>% filter(!is.na(meters_year))


#### Array Map
load("~/Documents/sdm_arrays/grid_25km.rda")
##############


for(sp in spResults$species) {

  if(!file.exists(glue::glue("{out_wd}/{sp}"))) {
    dir.create(glue::glue("{out_wd}/{sp}"))
  }
  
  if(!file.exists(glue::glue("{out_wd}/{sp}/disp_stars.rda"))) {
    load(glue::glue("{sdm_wd}/{sp}/Predictions/pred_stars.rda"))
    
    ### Maxent threshold
    maxentTrheshold <- as.numeric(strsplit(
      strsplit(readLines(glue::glue("{sdm_wd}/{sp}/MaxentModelOutput/species.html"))[12], "</th><th>")[[1]][5], "</td><td>")[[1]][20])
    
    current <- read_stars(glue::glue("{sdm_wd}/{sp}/{sp}_MaxEnt_calibration_restricted_50.tif")) %>% setNames("present") %>%
      mutate(present = ifelse(present>=maxentTrheshold, 1, 0))
    
    distance <- current %>%
      mutate(present = ifelse(present==1, 1, NA)) %>% rast() %>% terra::distance() %>% st_as_stars() %>% setNames("distance")
    
    sppList <- lapply(1:length(pred_stars), function(spp) {
      tmp <- pred_stars[spp] %>% split(drop = TRUE)
      lapply(1:length(tmp), function(x) {
        p_dist <- distance %>% 
          mutate(distance = 1 - pgamma(distance, shape = 1, rate = 1 / ((spResults %>% 
                                                                           filter(species==sp) %>% 
                                                                           pull(meters_year))*30*x)))
        (rast(tmp[x]) * rast(p_dist)) %>% st_as_stars() %>% st_set_crs(st_crs(pred_stars)) %>%
          setNames("tmp") %>% mutate(tmp = ifelse(tmp>=maxentTrheshold, 1, 0)) %>%
          setNames(names(tmp)[x]) 
      }) %>% do.call("c", .) %>% merge(name = 'years') %>% setNames(names(pred_stars)[spp])
    })
    
    disp_stars <- do.call('c', sppList)
    write_stars(current, glue::glue("{out_wd}/{sp}/current.tif"))
    save(disp_stars, file = glue::glue("{out_wd}/{sp}/disp_stars.rda"))
    
    rm(pred_stars); gc()
    
    # spp_list <- lapply(tibble(spp = rep(1:3, each = 3), years = rep(1:3, 3)) %>% group_split(spp, years), function(x) {
    #   out <- ggplot() + geom_stars(data = disp_stars[x$spp][,,,x$years], downsample = 5, show.legend = F) +
    #     scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
    #     coord_equal() +
    #     theme_void()
    # })
    # 
    # library(gridExtra)
    # spp_out <- do.call('grid.arrange', c(spp_list, nrow = 3, ncol = 3, left = "spp: [5.85, 3.70, 1.26]", top = "years: [2026, 2056, 2086]"))
    # ggsave(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_dispersal.png"), spp_out, units = "cm", width = 20, height = 20, bg = "grey90") 
    # dev.off()
    
    curr_extr <- st_extract(current, grid$geometry %>% st_transform(st_crs(current))) %>% st_as_sf() %>%
      st_drop_geometry()
    
    extr_pred <- st_extract(disp_stars, grid$geometry %>% st_transform(st_crs(disp_stars))) %>% st_as_sf() %>%
      st_drop_geometry()
    
    listOut   <- lapply(1:3, function(x) {
      out <- as.numeric(unlist(cbind(curr_extr[,1], extr_pred[,rbind(c(1:3), c(4:6), c(7:9))[x,]])))
      ifelse(is.na(out), 0, out)
    })
    
    binarArray <- array(dim = c(1, nrow(grid), 4, 3))
    binarArray[1, , , 1] <- listOut[[1]]
    binarArray[1, , , 2] <- listOut[[2]]
    binarArray[1, , , 3] <- listOut[[3]]
    
    save(binarArray, file = glue::glue("{out_wd}/{sp}/binArray.rda"))
  }
  
}