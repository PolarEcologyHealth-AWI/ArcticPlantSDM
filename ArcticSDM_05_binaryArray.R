### Distribution extraction
library(sf)
sf_use_s2(FALSE)
library(stars)
library(abind)
library(tidyverse)
library(terra)

root <- ""

{
  unres_wd <- glue::glue("{root}Results/SDM_biasCorrection/unconstraint")
  restr_wd <- glue::glue("{root}Results/SDM_biasCorrection/dispConstraint2")
  out_wd   <- glue::glue("{root}Results/SDM_biasCorrection/")
}

spTable <- tibble(species  = sapply(strsplit(list.files(unres_wd), '_'), function(x) x[1])) %>%
  bind_rows(tibble(species = sapply(strsplit(list.files(restr_wd), '_'), function(x) x[1]))) %>%
  filter(!duplicated(species))

save(spTable, file = glue::glue("{out_wd}/spTable.rda"))

### Grid
load("/data/grid_25km.rda")
########

predArray <- array(dim = c(nrow(spTable), nrow(grid), 4, 3, 2))

for(sp in spTable$species) {
  
  ### unrestricted
  tryCatch(load(glue::glue("{unres_wd}/{sp}_binArray_Unconstraint.rda")), error = function(e) NULL)
  if(!is.null(binarArray)) unc <- binarArray else unc <- array(dim = c(1, nrow(grid), 4, 3))
  
  ### restricted
  tryCatch(load(glue::glue("{restr_wd}/{sp}_binArray_dispRestricted.rda")), error = function(e) NULL)
  if(!is.null(binarArray)) con <- binarArray else con <- array(dim = c(1, nrow(grid), 4, 3))
  
  out <- abind::abind(unc, con, along = 5)
  
  predArray[which(spTable$species==sp),,,,] <- out
  
}

save(predArray, file = glue::glue("{out_wd}/predArray_Revision2_2.rda"))
