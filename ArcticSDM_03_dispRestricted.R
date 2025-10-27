suppressMessages({
  library(tidyverse)
  library(sf)
  sf_use_s2(FALSE)
  library(stars)
  library(terra)
})

run <- as.numeric(commandArgs(trailingOnly = T))

root <- ""

{
  sdm_wd <- glue::glue("{root}Results/SDM_biasCorrection/species/")
  out_wd <- glue::glue("{root}Results/SDM_biasCorrection/dispConstraint2")
}

speciesDisp <- suppressMessages(read_delim(glue::glue("{root}data/species_meters_all.csv"), delim = ";"))

spResults <- (tibble(species = list.files(sdm_wd)) %>%
                left_join(tibble(fls = list.files(sdm_wd, pattern = "distRestriction.tif", recursive = T)) %>%
                            mutate(species = sapply(strsplit(fls, "/"), function(x) x[[1]]), predict = TRUE) %>%
                            dplyr::select(-fls), by = "species") %>% filter(!is.na(predict)) %>%
                left_join(speciesDisp, by = "species")) %>% filter(!is.na(meters_year_50), !is.na(meters_year_99))

sp <- spResults$species[run]


if(!file.exists(glue::glue("{out_wd}/{sp}_binArray_dispRestricted.rda"))) {
    
  #### Array Map
  load("data/grid_25km.rda")
  ##############
    
  load(glue::glue("{sdm_wd}/{sp}/Predictions/pred_stars.rda"))
    
  ### Maxent threshold
  load(glue::glue("{sdm_wd}/{sp}/{sp}_maxent_mod.rda"))
  maxentTrheshold <- as.numeric(maxent_mod@results['X10.percentile.training.presence.Cloglog.threshold',1])
    
  current <- read_stars(glue::glue("{sdm_wd}/{sp}/{sp}_MaxEnt_calibration_distRestriction.tif")) %>% setNames("present") %>%
      mutate(present = ifelse(present>=maxentTrheshold, 1, 0))
    
  distance <- current %>%
      mutate(present = ifelse(present==1, 1, NA)) %>% rast() %>% terra::distance() %>% st_as_stars() %>% setNames("distance") %>%
      suppressWarnings()
    
  ### Sigmoidal dispersal
  sigm_params <-lapply(1:3, function(x) {
      p50 <- as.numeric(spResults %>% filter(species==sp) %>% dplyr::select(meters_year_50)) * 30 * x
      p99 <- as.numeric(spResults %>% filter(species==sp) %>% dplyr::select(meters_year_99)) * 30 * x
      s <- seq(0.0001, 0.1, length = 1000)
      c2  <- s[which.min(abs(0.99 - sapply(s, function(x) (1 / (1 + exp(-x*(p99 - p50)))))))]
      tibble(p50, c2)
    }) %>% Reduce("rbind",.)
    
  
  sppList <- lapply(1:length(pred_stars), function(spp) {
    
    tmp <- pred_stars[spp] %>% split(drop = TRUE)
    
    lapply(1:length(tmp), function(x) {
      
      p_dist <- distance %>% 
        mutate(runif = runif(prod(dim(.))),
               p     = ifelse(distance < 25, 1, 
                              runif >= (1 / (1 + exp(-as.numeric(sigm_params[x,2])*(distance - as.numeric(sigm_params[x,1]))))))) %>%
        dplyr::select(p) %>% st_set_crs(st_crs(pred_stars))
      
      (rast(tmp[x]) * rast(p_dist)) %>% st_as_stars() %>% st_set_crs(st_crs(pred_stars)) %>%
        setNames("tmp") %>% mutate(tmp = ifelse(tmp>=maxentTrheshold, 1, 0)) %>%
        setNames(names(tmp)[x]) 
      }) %>% do.call("c", .) %>% merge(name = 'years') %>% setNames(names(pred_stars)[spp])
    
    })
        
  disp_stars <- do.call('c', sppList)
        
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
        
  save(binarArray, file = glue::glue("{out_wd}/{sp}_binArray_dispRestricted.rda"))

}