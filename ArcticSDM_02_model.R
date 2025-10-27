suppressMessages({suppressWarnings({
library(sf)
sf_use_s2(FALSE)
library(stars)
library(terra)
library(tidyverse)
library(flexsdm)
library(e1071)
library(pROC)
library(patchwork)
library(modEvA)
library(dismo)

root <- ""

{
  wd        <- glue::glue("{root}ArcticSDM/")
  gbif_data <- glue::glue("{root}ArcticSDM_data/")
  map_data  <- glue::glue("{root}ArcticSDM_data/")
  env_data  <- glue::glue("{root}CHELSA data 5 km/Modern/")
  out_wd    <- glue::glue("{root}Results/SDM_biasCorrection")
  
  source(glue::glue("{root}/functions.R"))
}

##############
## GBIF Occ ##
##############

minSpecies <- 150

{
  ## Read in thinned occurrences
  load(glue::glue("{wd}Results/OccTab_all.rda"))
  
  spList <- OccTab_all %>% st_drop_geometry() %>% na.omit() %>% group_by(species) %>% summarise(count = n()) %>%
    arrange(desc(count)) %>% filter(count >= (100*minSpecies)/75)
}

### Species 
sp <- spList$species[as.numeric(commandArgs(trailingOnly = T))] 

if(!dir.exists(glue::glue("{out_wd}/species/{sp}"))) {

  ### Projection
  proj <- "+proj=laea +lon_0=-170 +lat_0=90"
  
  ### Maps
  ecoreg   <- st_read(glue::glue("{map_data}Ecoregions/tnc_terr_ecoregions.shp")) %>%
    filter(WWF_MHTNAM%in%c("Boreal Forests/Taiga", "Tundra"), st_coordinates(st_centroid(.))[,2]>0) %>% suppressWarnings()
  
  map_wrld <- st_read(glue::glue("{map_data}/ne_10m_admin_1_states_provinces/ne_10m_admin_1_states_provinces.shp")) %>%
    st_intersection(ecoreg %>% st_union()) %>% dplyr::select(c('name', 'admin'))  %>% suppressWarnings()
  
  map   <- map_wrld %>% st_transform(proj) %>% mutate(id = as.numeric(as.factor(name)))
  admin <- map %>% dplyr::select(id, name, admin) %>% st_drop_geometry() %>% filter(!duplicated(id))
  # plot(map %>% dplyr::select(name))
  
  
  
  #######################
  ## SDM PreProcessing ##
  #######################
  
  ### Current
  {
    env_list <- list.files(env_data, pattern = ".tif")
    # env_curr <- lapply(env_list, function(x) {
    #   read_stars(glue::glue("{env_data}{x}")) %>% setNames(strsplit(x, "_")[[1]][2])
    # }) %>% Reduce("c", .) %>% merge()
    # 
    # predTab_curr <- env_curr %>% st_as_sf() %>% filter(apply(., 1, function(x) all(!is.na(x))))
    # 
    # save(predTab_curr, file = glue::glue("{wd}/Results/env/predTab_curr.rds"))
    load(glue::glue("{wd}/Results/env/predTab_curr.rds"))
  }
  
  # ### Future
  {
    fut_wd   <- glue::glue("{gsub('Modern/', '', env_data)}/Future_CMIP6")
    fls      <- list.files(fut_wd, recursive = T)[grep(".tif", list.files(fut_wd, recursive = T))]
    
    futureList   <- tibble(fls  = fls,
                           year = sapply(strsplit(fls, "/"), function(x) round(mean(as.numeric(unlist(strsplit(x[[1]], "-")))),0)),
                           spp  = sapply(strsplit(fls, "/"), function(x) x[[2]]),
                           var  = sapply(strsplit(fls, "_"), function(x) x[2])) %>% arrange(year, spp, var) %>%
      filter(!is.na(var)) %>% group_split(spp)
    
    # env_future <- lapply(futureList, function(x) {
    #   year_stars <- x %>% group_split(year) %>% lapply(function(y) {
    #     read_stars(glue::glue("{fut_wd}/{y$fls}")) %>% setNames(y$var) %>% st_as_sf() %>% st_drop_geometry() %>%
    #       as.matrix()
    #   }) %>% abind::abind(., along = 3)
    # }) %>% abind::abind(., along = 4)
    #
    # save(env_future, file = glue::glue("{wd}/Results/env/env_future.rds"))
    load(glue::glue("{wd}/Results/env/env_future.rds"))
    
    predTab_future0 <- read_stars(glue::glue("{fut_wd}/{fls[1]}")) %>% 
      setNames("pred_future")# %>% st_as_sf()
  }
  
  
  #######################
  ## SDM               ##
  #######################
  
  ### Occurance
  occ_sp_all <- OccTab_all %>% na.omit() %>% filter(species == sp) %>%
    st_intersection(map) %>% dplyr::select(species, name, admin, names(predTab_curr %>% st_drop_geometry())) %>%
    mutate(id = 1:nrow(.), .before = "species") %>% suppressWarnings()
  
  ### Spatial thinning with distance
  {
    distance <- 150
    
    ssPts <- sapply(3:4, function(x) (st_bbox(map)[x]*2)/1000/distance)
    
    grd   <- st_sample(map %>% st_union(), ssPts[1]*ssPts[2], type = "random") %>%
      st_intersection(map %>% st_buffer(distance*1000)) %>% suppressWarnings()
    
    occ_sp <- lapply(1:length(grd), function(x) {
      occ_sp_all %>% mutate(dist = as.numeric(st_distance(., grd[x]))) %>%
        arrange(dist) %>% filter(dist<(distance*1000)*2) %>%
        slice(1) %>% dplyr::select(-dist)
    }) %>% do.call("rbind", .) %>% 
    mutate(id = 1:nrow(.), .before = "species") %>%
    mutate(training = (1:nrow(.))%in%sample(1:nrow(.), round(nrow(.)*0.75), 0), .after = "species")
  }
  
  
  if(sum(occ_sp$training)>=minSpecies) {
  
    if(!file.exists(glue::glue("{out_wd}/species/{sp}"))) {
      dir.create(glue::glue("{out_wd}/species/{sp}"))
      dir.create(glue::glue("{out_wd}/species/{sp}/Predictions"))
    }
    
    
    ### Co-linearity
    {
      ### Remove specific environmental variables
      env_delete <- c("bio3", "swe", "scd")
      
      env  <- occ_sp %>% st_drop_geometry() %>% 
        dplyr::select(all_of(names(predTab_curr %>% st_drop_geometry()))) %>%
        dplyr::select(-all_of(env_delete))
      
      vif  <- custom_vif(env, c("bio1", "bio12"), th = "10")
      
      collSpecies <- tibble(species = sp) %>%
        bind_cols(
          tibble(var = names(predTab_curr %>% st_drop_geometry())) %>%
            mutate(incl = ifelse(var%in%vif$removed_variables, FALSE, TRUE)) %>%
            pull(incl) %>% as.matrix(nrow = 1) %>% t() %>%
            as_tibble() %>% setNames(c(names(predTab_curr %>% st_drop_geometry())))
        )
      
      bioClim_filter <- vif$kept_variables
      
    }
    
    ### Environmental data
    {
      env_curr <- lapply(env_list[sapply(strsplit(env_list, "_"), function(x) x[2]) %in% bioClim_filter], function(x) {
        read_stars(glue::glue("{env_data}{x}")) %>% setNames(strsplit(x, "_")[[1]][2])
      }) %>% Reduce("c", .) %>% merge()
      
      
      bfr <- occ_sp %>% st_buffer(150000) %>% st_union() %>%
        st_sym_difference(occ_sp %>% st_buffer(50000) %>% st_union() %>% 
                            st_simplify(dTolerance = 10000)) %>%
        st_intersection(map %>% st_union())
      
      ## svm
      svm_tab <- occ_sp %>% mutate(p=1) %>%
        dplyr::select(p, bioClim_filter) %>%
        bind_rows(
          st_sample(bfr, 50000) %>% st_sfc() %>% 
            st_extract(env_curr %>% st_set_crs(st_crs(occ_sp)), .) %>% st_as_sf() %>%
            mutate(p = 0, .before = names(.)[1])) %>% na.omit() %>% suppressWarnings()
      
      svm_mod <- svm(svm_tab %>% st_drop_geometry() %>% dplyr::select(bioClim_filter), 
                     svm_tab %>% st_drop_geometry() %>% pull(p), type = "one-classification")
      
      abs_svm <- svm_tab %>% filter(p == 0 & !predict(svm_mod)) %>%
        mutate(training = (1:nrow(.))%in%sample(1:nrow(.), round(nrow(.)*0.75), 0))
      
      modelTab <- occ_sp %>% mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2], p = 1, .before = 'training') %>%
        st_drop_geometry() %>% dplyr::select(species, name, admin, lon, lat, training, p, bioClim_filter) %>%
        bind_rows(
          abs_svm %>%  mutate(species = sp) %>% st_intersection(map) %>%
            mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
            st_drop_geometry() %>% dplyr::select(species, name, admin, lon, lat, training, p, bioClim_filter)
        ) %>% suppressWarnings()
      
      save(modelTab, file = glue::glue("{out_wd}/species/{sp}/{sp}_modelTab.rds"))
      
    }
    
    ### MaxEnt Model
    {
      args_list <- c('responsecurves=FALSE',
                     'jackknife=FALSE',
                     'pictures=FALSE',
                     'autofeature=FALSE',
                     'linear=TRUE',
                     'quadratic=TRUE',
                     'product=TRUE',
                     'threshold=TRUE',
                     'hinge=TRUE',
                     'betamultiplier=1')
      
      
      maxentTab   <- modelTab %>% filter(training) %>%
        dplyr::select(p, bioClim_filter)
      
      maxent_mod  <- dismo::maxent(maxentTab[,-1], p = maxentTab$p, nbg = 0, args = args_list)
      save(maxent_mod, file = glue::glue("{out_wd}/species/{sp}/{sp}_maxent_mod.rda"))
      
      pred_maxent <- predTab_curr %>% mutate(pred = dismo::predict(maxent_mod, predTab_curr %>% st_drop_geometry() %>%
                                                                     dplyr::select(bioClim_filter), type = "response"))
      
      rastOut     <- st_rasterize(pred_maxent %>% dplyr::select(pred), st_as_stars(st_bbox(env_curr),
                                                                                   nx = st_dimensions(env_curr)$x$delta, ny = st_dimensions(env_curr)$x$delta, values = NA_real_)) %>%
        setNames(glue::glue("{sp}_MaxEnt_current"))
      
      write_stars(rastOut, glue::glue("{out_wd}/species/{sp}/{sp}_MaxEnt_calibration.tif"))
      
      
      ### Evaluation
      evalRows <- tibble(species = sp, admin = c("all", unique(map$admin)))
      evalTab  <- evalRows %>% bind_cols(lapply(1:nrow(evalRows), function(x) {
        if(evalRows$admin[x]!="all") {
          tmp <- modelTab %>% filter(admin==evalRows$admin[x])
        } else {
          tmp <- modelTab
        }
        
        if(nrow(tmp)>0) {
          tibble(occ_training    = sum(tmp %>% filter(training) %>% pull(p)), 
                 env_training    = sum(!(tmp %>% filter(training) %>% pull(p))), 
                 occ_eval        = sum(tmp %>% filter(!training) %>% pull(p)), 
                 env_eval        = sum(!(tmp %>% filter(!training) %>% pull(p))),
                 auc_training    = tryCatch(with(tmp %>% filter(training) %>% 
                                                   st_as_sf(coords = c("lon", "lat"), crs = st_crs(rastOut)) %>% 
                                                   mutate(pred = st_extract(rastOut,.) %>% st_drop_geometry() %>% setNames("pred") %>% pull(pred)), as.numeric(auc(p, pred))), error = function(e) NA), 
                 auc_eval        = tryCatch(with(tmp %>% filter(!training) %>% 
                                                   st_as_sf(coords = c("lon", "lat"), crs = st_crs(rastOut)) %>% 
                                                   mutate(pred = st_extract(rastOut,.) %>% st_drop_geometry() %>% setNames("pred") %>% pull(pred)), as.numeric(auc(p, pred))), error = function(e) NA), 
                 boyes_training  = tryCatch(with(tmp %>% filter(training, p==1), 
                                                 modEvA::Boyce(obs = cbind(lon,lat), pred = rast(rastOut), verbosity = 0, plot = FALSE)$Boyce), error = function(e) NA), 
                 boyes_eval      = tryCatch(with(tmp %>% filter(!training, p==1), 
                                                 modEvA::Boyce(obs = cbind(lon,lat), pred = rast(rastOut), verbosity = 0, plot = FALSE)$Boyce), error = function(e) NA)) %>%
            suppressMessages() %>% suppressWarnings() %>% invisible()
        } else {
          tibble(occ_training = NA, env_training = NA, 
                 occ_eval = NA, env_eval = NA,
                 auc_training = NA, auc_eval = NA, boyes_training = NA, boyes_eval = NA)
        }
        
      }) %>% Reduce("rbind", .))
      
      write.table(evalTab, file = glue::glue("{out_wd}/species/{sp}/{sp}_evalTab.txt"))
      
    }
    
    ### Distance restricted predictions
    {
      if(!file.exists(glue::glue("{out_wd}/species/{sp}/{sp}_MaxEnt_calibration_restricted.tif"))) {
        
        model_dist <- modelTab %>% filter(p==1) %>% st_as_sf(coords = c("lon", "lat"), crs = proj) %>%
          filter(apply(st_distance(.), 1, function(x) min(x[x>0])) < 1500*1000)
        
        
        ### geographical distance
        distRast <- model_dist %>% mutate(p = 1) %>% dplyr::select(p) %>% st_rasterize(
          rastOut %>% setNames('dist') %>% mutate(dist = NA)) %>% rast() %>% distance() %>% suppressWarnings()
        
        mean <- 1000*1000
        sd   <- 250*1000
        
        rastOut <-  rastOut %>% setNames("vals") %>% mutate(vals = vals * (1 - pnorm(distRast[], mean, sd))) %>%
          setNames(names(rastOut))
        
        write_stars(rastOut, glue::glue("{out_wd}/species/{sp}/{sp}_MaxEnt_calibration_distRestriction.tif"))
        
        outMap <- ggplot() +
          geom_stars(data = rastOut, downsample = 5, show.legend = F) +
          scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
          labs(title = glue::glue("{sp}: Current projection (MaxEnt)"), x = "", y = "") +
          geom_sf(data = model_dist, pch = 16, cex = 0.25) +
          theme_void() +
          theme(plot.title    = element_text(size = 9, face = "bold.italic", hjust = 0.5),
                plot.subtitle = element_text(size = 8, face = "bold.italic", hjust = 0.5))
        
        ggsave(glue::glue("{out_wd}/species/{sp}/{sp}_MaxEnt_calibration_distRestriction.png"), outMap, units = "cm", width = 18, height = 15, bg = "white")
        
      }
    }
    
    ### Predictions
    {
      ## future
      if(!file.exists(glue::glue("{out_wd}/species/{sp}/Predictions/pred_stars.rda"))) {
        envs    <- names(maxent_mod@presence)
        
        predList <- list()
        for(spp in 1:length(futureList %>% Reduce("rbind",.) %>% pull(spp) %>% unique())) {
          predList[[spp]] <- matrix(nrow = dim(env_future)[1], ncol = 3)
          predList[[spp]][,1] <- dismo::predict(maxent_mod, env_future[,envs,1,spp], type = "response"); gc()
          predList[[spp]][,2] <- dismo::predict(maxent_mod, env_future[,envs,2,spp], type = "response"); gc()
          predList[[spp]][,3] <- dismo::predict(maxent_mod, env_future[,envs,3,spp], type = "response"); gc()
        }
        predTab_future  <- predTab_future0 %>% st_as_sf() %>% bind_cols(do.call("cbind", predList) %>% as_tibble()) %>% dplyr::select(-pred_future)
        invisible(gc())
        
        starsList <- parallel::mclapply(1:ncol(predTab_future %>% st_drop_geometry()), function(x) {
          st_rasterize(predTab_future[,x], st_as_stars(st_bbox(env_curr), nx = st_dimensions(env_curr)$x$delta, ny = st_dimensions(env_curr)$x$delta, values = NA_real_))
        }, mc.cores = ncol(predTab_future %>% st_drop_geometry()))
        
        
        pred_stars <- lapply(tibble(ind   = 1:length(starsList),
                                    spp   = rep(futureList %>% Reduce('rbind',.) %>% pull(spp) %>% unique(), each = 3),
                                    years = rep(futureList %>% Reduce('rbind',.) %>% pull(year) %>% unique(), 3)) %>% group_split(spp),
                             function(spp) {
                               do.call("c", starsList[spp$ind]) %>% setNames(as.factor(spp$years)) %>% merge(name = "years") %>% setNames(spp$spp[1])
                             }) %>% do.call("c", .)
        
        save(pred_stars, file = glue::glue("{out_wd}/species/{sp}/Predictions/pred_stars.rda"))
      }
    }
    
    ### Maps Future
    {
      spp_list <- lapply(tibble(spp = rep(1:3, each = 3), years = rep(1:3, 3)) %>% group_split(spp, years), function(x) {
        out <- ggplot() + geom_stars(data = pred_stars[x$spp][,,,x$years], downsample = 5, show.legend = F) +
          scale_fill_gradient2(low = 'white', mid = "orange", high = 'darkred', breaks = seq(0, 1, length = 100), na.value = "grey70") +
          coord_equal() +
          theme_void()
      })
      
      library(gridExtra)
      spp_out <- do.call('grid.arrange', c(spp_list, nrow = 3, ncol = 3, left = "spp: [5.85, 3.70, 1.26]", top = "years: [2026, 2056, 2086]"))
      ggsave(glue::glue("{out_wd}/species/{sp}/Predictions/{sp}_MaxEnt_predictions.png"), spp_out, units = "cm", width = 20, height = 20, bg = "grey90")
      dev.off()
    }
    
  }
}

})})

print(glue::glue("{as.numeric(commandArgs(trailingOnly = T))}-{sp}"))
