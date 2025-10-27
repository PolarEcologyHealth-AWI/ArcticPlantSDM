correct_colinvar <- function(env_layer,
                             method,
                             proj = NULL) {
  
  . <- NULL
  
  if (!any(c("pearson", "vif", "pca", "fa") %in% method)) {
    stop(
      "argument 'method' was misused, select one of the available methods: pearson, vif, pca, fa"
    )
  }
  
  
  if (any(method %in% "pearson")) {
    if (is.na(method["th"])) {
      th <- 0.7
    } else {
      th <- as.numeric(method["th"])
    }
    
    h <- env_layer %>% stats::na.omit()
    h <- abs(stats::cor(h, method = "pearson"))
    diag(h) <- 0
    
    cor_var <- h>th
    cor_var <- apply(cor_var,2, function(x) colnames(h)[x])
    if(length(cor_var)==0){
      cor_var <- 'No pair of variables reached the specified correlation threshold.'
    }
    
    result <- list(
      cor_table = h,
      cor_variables = cor_var
    )
  }
  
  
  if (any(method %in% "vif")) {
    
    if (is.null(method["th"])) {
      th <- 10
    } else {
      th <- as.numeric(method["th"])
    }
    
    x <- as.data.frame(env_layer) %>% stats::na.omit()
    
    LOOP <- TRUE
    if (nrow(x) > 200000) {
      x <- x[sample(1:nrow(x), 200000), ]
    }
    n <- list()
    n$variables <- colnames(x)
    exc <- c()
    
    while (LOOP) {
      v <- rep(NA, ncol(x))
      names(v) <- colnames(x)
      for (i in 1:ncol(x)) {
        suppressWarnings(v[i] <- 1 / (1 - summary(lm(x[, i] ~ ., data = x[-i]))$r.squared))
      }
      if (v[which.max(v)] >= th) {
        ex <- names(v[which.max(v)])
        exc <- c(exc, ex)
        x <- x[, -which(colnames(x) == ex)]
      } else {
        LOOP <- FALSE
      }
    }
    if (length(exc) > 0) {
      n$excluded <- exc
    }
    
    v <- rep(NA, ncol(x))
    names(v) <- colnames(x)
    for (i in 1:ncol(x)) {
      v[i] <- 1 / (1 - summary(stats::lm(x[, i] ~ ., data = x[-i]))$r.squared)
    }
    
    # n$corMatrix <- stats::cor(x, method = "pearson")
    n$results <- data.frame(Variables = names(v), VIF = as.vector(v))
    
    result <- list(
      removed_variables = n$excluded,
      vif_table = dplyr::tibble(n$results)
    )
  }
  
  # if (any(method %in% "pca")) {
  #   
  #   # mean
  #   means <- t(terra::global(env_layer, 'mean', na.rm=T)) %>% c()
  #   names(means) <- names(env_layer)
  #   # SD
  #   stds <- t(terra::global(env_layer, 'sd', na.rm=T)) %>% c()
  #   names(stds) <- names(env_layer)
  #   
  #   # Standardize raster values
  #   env_layer <- terra::scale(env_layer, center = means, scale = stds)
  #   vnmes <- names(means)
  #   
  #   
  #   if(is.null(maxcell)){
  #     p0 <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)
  #   } else {
  #     # Raster random sample
  #     set.seed(10)
  #     p0 <- terra::as.data.frame(env_layer[[1]], cells=TRUE)[,1] %>%
  #       sample(., size = maxcell, replace = FALSE) %>%
  #       sort()
  #     p0 <- env_layer[p0] %>%
  #       stats::na.omit()
  #   }
  #   
  #   p <- stats::prcomp(p0,
  #                      retx = TRUE,
  #                      scale. = FALSE,
  #                      center = FALSE
  #   )
  #   cof <- p$rotation
  #   
  #   cvar <- summary(p)$importance["Cumulative Proportion", ]
  #   naxis <- Position(function(x) {
  #     x >= 0.95
  #   }, cvar)
  #   cvar <- data.frame(cvar)
  #   
  #   
  #   # p <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)
  #   p <- stats::prcomp(p0, retx = TRUE, scale. = FALSE, center = FALSE, rank. = naxis)
  #   env_layer <- terra::predict(env_layer, p)
  #   
  #   rm(p0)
  #   
  #   result <- list(
  #     env_layer = env_layer,
  #     coefficients = data.frame(cof) %>% dplyr::tibble(variable = rownames(.), .),
  #     cumulative_variance = dplyr::tibble(PC = 1:nrow(cvar), cvar)
  #   )
  #   
  #   if (!is.null(proj)) {
  #     dpca <- file.path(dirname(proj), "Projection_PCA")
  #     dir.create(dpca)
  #     subfold <- list.files(proj)
  #     subfold <- as.list(file.path(dpca, subfold))
  #     sapply(subfold, function(x) {
  #       dir.create(x)
  #     })
  #     
  #     proj <- list.files(proj, full.names = TRUE)
  #     for (i in 1:length(proj)) {
  #       scen <- terra::rast(list.files(proj[i], full.names = TRUE))
  #       scen <- scen[[names(means)]]
  #       scen <- terra::scale(scen, center = means, scale = stds)
  #       scen <- terra::predict(scen, p)
  #       terra::writeRaster(
  #         scen,
  #         file.path(subfold[[i]], "pcs.tif"),
  #         overwrite=TRUE
  #       )
  #     }
  #     
  #     result$proj <- dpca
  #   }
  # }
  # 
  # if (any(method %in% "fa")) {
  #   p <- terra::scale(env_layer, center = TRUE, scale = TRUE)
  #   
  #   if(is.null(maxcell)){
  #     p <- terra::as.data.frame(p, xy = FALSE, na.rm = TRUE)
  #   } else {
  #     # Raster random sample
  #     set.seed(10)
  #     p <- terra::as.data.frame(env_layer[[1]], cells=TRUE)[,1] %>%
  #       sample(., size = maxcell, replace = FALSE) %>%
  #       sort()
  #     p <- env_layer[p] %>%
  #       stats::na.omit()
  #   }
  #   
  #   if (nrow(p) > 200000) {
  #     p <- p[sample(1:nrow(p), 200000), ]
  #   }
  #   
  #   e <- eigen(stats::cor(p))
  #   len <- length(e$values)
  #   a <- NULL
  #   r <- NULL
  #   
  #   for (j in 1:len) {
  #     a[j] <- 1 / len * sum(1 / (j:len))
  #     r[j] <- e$values[j] / (sum(e$values))
  #   }
  #   
  #   ns <- length(which(r > a))
  #   
  #   fit <-
  #     tryCatch(
  #       stats::factanal(
  #         x = p,
  #         factors = ns,
  #         rotation = "varimax",
  #         lower = 0.001
  #       ),
  #       error = function(e) {
  #         stats::factanal(
  #           x = p,
  #           factors = ns - 1,
  #           rotation = "varimax",
  #           lower = 0.001
  #         )
  #       }
  #     )
  #   
  #   sel <-
  #     row.names(fit$loadings)[apply(fit$loadings, 2, which.max)]
  #   rem <-
  #     row.names(fit$loadings)[!row.names(fit$loadings) %in% sel]
  #   
  #   env_layer <- terra::subset(env_layer, sel)
  #   
  #   h <- fit$loadings %>%
  #     matrix() %>%
  #     data.frame()
  #   colnames(h) <- paste("Factor", 1:ncol(h), sep = "_")
  #   
  #   result <- list(
  #     env_layer = env_layer,
  #     number_factors = fit$factors,
  #     removed_variables = rem,
  #     uniqueness = fit$uniquenesses,
  #     loadings = dplyr::tibble(Variable = rownames(fit$loadings), h)
  #   )
  # }
  
  return(result)
}

occThinning <- function(occ_sp, threshold = 0.3, adjust = 0.25, cores = 3) {
  
  occ_ppp     <- spatstat.geom::as.ppp(occ_sp$geometry)
  occ_density <- st_as_stars(density(occ_ppp, adjust = adjust, dimyx = 500)) %>% st_set_crs(st_crs(occ_sp))
  
  # ggplot() +
  #   geom_stars(data = occ_density) +
  #   scale_fill_viridis_b() +
  #   geom_sf(data = occ_sp$geometry, size = 0.1)
  
  occ_sp_dens <- occ_sp %>% mutate(dens = st_extract(occ_density, occ_sp) %>% pull(v)) %>%
    mutate(dens = ifelse(is.na(dens), 0, dens))
  
  seq <- seq(1, 0.01, by = -0.025)
  test <- parallel::mclapply(seq(1, 0.01, by = -0.025), function(y) {
    tmp <- occ_sp_dens[sample(1:nrow(occ_sp_dens), nrow(occ_sp_dens)*y, prob = scales::rescale(occ_sp_dens %>% pull(dens), c(1, 0.1))),]
    
    tmp_ppp     <- spatstat.geom::as.ppp(tmp$geometry)
    tmp_density <- st_as_stars(density(tmp_ppp, adjust = adjust, dimyx = 500)) %>% st_set_crs(st_crs(occ_sp))
    
    tmp_dens <- tmp %>% mutate(dens = st_extract(tmp_density, tmp) %>% pull(v)) %>%
      mutate(dens = ifelse(is.na(dens), 0, dens))
    
    median(tmp_dens %>% pull(dens))
  }, mc.cores = cores)
  
  thresh_cp = beast(unlist(test), season='none', quiet = T, print.progress = F)
  tab <- tibble(percent = seq, data = unlist(test), p = thresh_cp$trend$cpOccPr)
  
  # plot(tab$percent, tab$data)
  # abline(h = tab$percent[which.max(tab$p)])
  # plot(tab$percent, tab$p, type = "o")
  
  if(any(tab$p>threshold)) {
    thresh <- tab %>% filter(p>threshold) %>% slice(nrow(.)) %>% pull(percent)
  } else thresh <- 1
  
  occ_sp_dens[sample(1:nrow(occ_sp_dens), nrow(occ_sp_dens)*thresh, prob = scales::rescale(occ_sp_dens %>% pull(dens), c(1, 0.1))),]
  
}

custom_vif <- function(env, key_vars, th = 10) {
  x <- as.data.frame(env) %>% stats::na.omit()
  exc <- c()
  LOOP <- TRUE
  
  while (LOOP) {
    v <- rep(NA, ncol(x))
    names(v) <- colnames(x)
    for (i in 1:ncol(x)) {
      suppressWarnings(v[i] <- 1 / (1 - summary(lm(x[, i] ~ ., data = x[-i]))$r.squared))
    }
    max_vif <- max(v)
    max_var <- names(which.max(v))
    
    if (max_vif >= th) {
      # If max_var is a key variable, find the highest VIF non-key variable that is collinear with it
      if (max_var %in% key_vars) {
        # Find all non-key variables with VIF above threshold
        non_key_vif <- v[!(names(v) %in% key_vars) & v >= th]
        if (length(non_key_vif) == 0) {
          # No non-key variable left to remove, break loop
          LOOP <- FALSE
        } else {
          # Remove the non-key variable with highest VIF
          ex <- names(which.max(non_key_vif))
          exc <- c(exc, ex)
          x <- x[, -which(colnames(x) == ex)]
        }
      } else {
        # Remove the variable with highest VIF (not a key variable)
        exc <- c(exc, max_var)
        x <- x[, -which(colnames(x) == max_var)]
      }
    } else {
      LOOP <- FALSE
    }
  }
  list(
    removed_variables = exc,
    kept_variables = colnames(x)
  )
}