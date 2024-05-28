library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(terra)


{
  sdm_wd <- "/Volumes/projects/p_ecohealth/Arctic_SDM/Results/"
  out_wd <- "/Volumes/projects/p_ecohealth/Arctic_SDM/migResults/"
}


spResults <- tibble(species = list.files(sdm_wd)) [c(3,15),]%>%
  mutate(DispRate = c(500, 1500))

for(sp in spResults$species) {
  
  if(!file.exists(glue::glue("{out_wd}/{sp}"))) {
    dir.create(glue::glue("{out_wd}/{sp}"))
  }
  
  load(glue::glue("{sdm_wd}/{sp}/Predictions/pred_stars.rda"))
  distance <- read_stars(glue::glue("{sdm_wd}/{sp}/{sp}_MaxEnt_calibration_restricted_50.tif")) %>% setNames("present") %>%
    mutate(present = ifelse(present>0.5, 1, NA)) %>% rast() %>% terra::distance() %>% st_as_stars() %>% setNames("distance")

  sppList <- lapply(1:length(pred_stars), function(spp) {
    tmp <- pred_stars[spp] %>% split(drop = TRUE)
    lapply(1:length(tmp), function(x) {
      p_dist <- distance %>% 
        mutate(distance = 1 - pgamma(distance, shape = 1, rate = 1 / ((spResults %>% 
               filter(species==sp) %>% 
               pull(DispRate))*30*x)))
      (rast(tmp[x]) * rast(p_dist)) %>% st_as_stars() %>% st_set_crs(st_crs(pred_stars)) %>%
          setNames(names(tmp)[x]) 
    }) %>% do.call("c", .) %>% merge(name = 'years') %>% setNames(names(pred_stars)[spp])
  })
  
  disp_stars <- do.call('c', sppList)
  save(disp_stars, file = glue::glue("{out_wd}/{sp}/disp_stars.rda"))
  
  spp_list <- lapply(tibble(spp = rep(1:3, each = 3), years = rep(1:3, 3)) %>% group_split(spp, years), function(x) {
    out <- ggplot() + geom_stars(data = disp_stars[x$spp][,,,x$years], downsample = 5, show.legend = F) +
      scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
      coord_equal() +
      theme_void()
  })
  
  library(gridExtra)
  spp_out <- do.call('grid.arrange', c(spp_list, nrow = 3, ncol = 3, left = "spp: [5.85, 3.70, 1.26]", top = "years: [2026, 2056, 2086]"))
  ggsave(glue::glue("{out_wd}/{sp}/{sp}_MaxEnt_dispersal.png"), spp_out, units = "cm", width = 20, height = 20, bg = "grey90") 
  dev.off()
  
}
