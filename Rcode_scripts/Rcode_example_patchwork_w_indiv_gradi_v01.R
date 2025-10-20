# -loop over species
for (spec.lat in unq.spc_n) {
  print(spec.lat)
  # get the color for the species
  clf.spc <- unique(df_e13$col.f_spc[match(spec.lat, df_e13$Latinsk_navn)])
  # get an index number for the species
  idxNosp <- which(spec.lat == unq.spc_n)
  # use the index number for the species appendix number
  no.spc.app.plot <- as.character(idxNosp)
  # pad with zeros to 2 characters for 'no.spc.app.plot'
  no.spc.app.plot <- ifelse(nchar(no.spc.app.plot) < 2, stringr::str_pad(no.spc.app.plot, 2, pad = "0"), idxNosp)
  no.spc.app.plot <- as.character(no.spc.app.plot)
  # get the latin species name without underscore
  spec.lat.w_undersc <- paste(sub(' ', '_', spec.lat))
  spec.lat.no_undersc <- spec.lat
  # subset the dataframe based on variable value in column
  df_e14 <- df_e13[which(df_e13$LatNm_wu == spec.lat.w_undersc), ]
  # use two functions together
  genusl <- substr(spec.lat.w_undersc, 1, 1)
  spcl <- gsub(".*_", "", spec.lat.w_undersc)
  # and paste them together
  short.spec.lat <- paste(genusl, "_", spcl, sep = "")
  # pad with zeros to 2 characters for
  TwN <- ifelse(nchar(idxNosp) < 2, stringr::str_pad(idxNosp, 2, pad = "0"), idxNosp)
  TwN <- as.character(TwN)
  sbs.AssIDNo <- TwN
  
  # Initialize a list to hold all plots for this species
  all_plots <- list()
  plot_counter <- 1
  
  # --- Calculate global max value for color scale ---
  zvm_global <- ceiling(max(df_e14$l10cp_L, na.rm = TRUE))
  
  # loop over years sampled
  for (yr_smpl in yrs) {
    print(yr_smpl)
    # subset the original dataframe per year sampled
    df_e15 <- df_e14[which(df_e14$Dato_inds.yy == yr_smpl), ]
    # iterate over the season in the vector
    for (season in categories.of.seasons) {
      print(season)
      # use match to match the season with a data frame and get the name for the season
      spcfc_seaon_name <- seaons_nms_df$names.of.seasons[match(season, seaons_nms_df$categories.of.seasons)]
      spcfc_seaon_name <- as.character(spcfc_seaon_name)
      
      # iterate over the machines in the vector
      for (mch_tp in mchns) {
        print(mch_tp)
        # subset in the data frame to get the data for the season and the machine type
        df_e15_ssn <- df_e15[which((df_e15$ssn.smpl == season & df_e15$mch == mch_tp)), ]
        
        # Check if there are enough rows for interpolation
        if (nrow(df_e15_ssn) >= 3) {
          # --- ipdw interpolation ---
          # get minimum and maximum to define range for jitter of points
          M27_jmin_lon <- min(df_e15_ssn$lon.m)
          M27_jmax_lon <- max(df_e15_ssn$lon.m)
          jit_lon <- (M27_jmax_lon - M27_jmin_lon) / 80000
          # get minimum and maximum to define range for jitter of points
          M27_jmin_lat <- min(df_e15_ssn$lat.m)
          M27_jmax_lat <- max(df_e15_ssn$lat.m)
          jit_lat <- (M27_jmax_lat - M27_jmin_lat) / 80000
          # jitter points to work around overlapping points
          df_e15_ssn$jit.lok_pos_lon <- jitter(df_e15_ssn$lon.m, jit_lon)
          df_e15_ssn$jit.lok_pos_lat <- jitter(df_e15_ssn$lat.m, jit_lat)
          # make SpatialPointsDataFrame with decimal degree coordinates
          pnts2 <- SpatialPointsDataFrame(df_e15_ssn[, c("jit.lok_pos_lon", "jit.lok_pos_lat")], df_e15_ssn, proj4string = crs(epsg4326nCRS2))
          pnts3 <- SpatialPointsDataFrame(df_e15_ssn[, c("lon.m", "lat.m")], df_e15_ssn, proj4string = crs(epsg4326nCRS2))
          df_e15_ssn$lok_pos_lon.f.pnts4 <- df_e15_ssn$jit.lok_pos_lon + 0
          df_e15_ssn$lok_pos_lat.f.pnts4 <- df_e15_ssn$jit.lok_pos_lat + 0.16
          pnts4 <- SpatialPointsDataFrame(df_e15_ssn[, c("lok_pos_lon.f.pnts4", "lok_pos_lat.f.pnts4")], df_e15_ssn, proj4string = crs(epsg4326nCRS2))
          
          # --- Use coastline as barrier ---
          bbox_k2 <- raster::buffer(as(extent(sp::spTransform(pnts2, projection(coastline10))), "SpatialPolygons"), width = 7)
          projection(bbox_k2) <- projection(coastline10)
          library(sf)
          st_is_valid(coastline10, reason = TRUE)[!st_is_valid(coastline10)]
          sf_use_s2(use_s2 = FALSE)
          csl10_crop <- st_crop(coastline10, bbox_k2)
          csl10_crop <- st_crop(coastline10, xmin = 6, xmax = 17, ymin = 54, ymax = 59)
          pols2 <- csl10_crop
          sf_use_s2(use_s2 = TRUE)
          if (is.data.frame(pols2) == T) {
            pols2 <- st_transform(pols2, projection(pnts2))
          } else {
            pols2 <- sp::spTransform(pols2, projection(pnts2))
          }
          costras2 <- costrasterGen(pnts2, pols2, extent = "polys", projstr = projection(pols2), resolution = res_fac)
          pnts <- pnts2
          costras <- costras2
          library(spatstat)
          W <- owin(range(coordinates(pnts)[, 1]), range(coordinates(pnts)[, 2]))
          kat.pp <- ppp(coordinates(pnts)[, 1], coordinates(pnts)[, 2], window = W)
          mean.neighdist <- mean(nndist(kat.pp))
          gridsize <- mean.neighdist * 1 * r_mnd
          r_mnd1.lonlat <- (mean.neighdist * 1 * r_mnd)
          grainscale.fac <- gridsize / res(costras)[1]
          gridras <- aggregate(costras, fact = grainscale.fac)
          gridpol <- rasterToPolygons(gridras)
          gridpol$value <- row.names(gridpol)
          fulldataset.over <- over(pnts, gridpol)
          fulldataset.over <- cbind(data.frame(fulldataset.over), setNames(data.frame(pnts), c(colnames(data.frame(pnts)))))
          set.seed(2)
          gridlev <- unique(fulldataset.over$value)
          for (i in seq_along(gridlev)) {
            activesub <- subset(fulldataset.over, fulldataset.over$value == gridlev[i])
            # Skip if no rows in this grid cell
            if (nrow(activesub) == 0) next
            selectnum <- gdata::resample(seq_len(nrow(activesub)), 1)
            if (i == 1) {
              training <- activesub[selectnum, ]
            } else {
              training <- rbind(training, activesub[selectnum, ])
            }
          }
          validate <- fulldataset.over[!(row.names(fulldataset.over) %in% row.names(training)), ]
          xy <- cbind(training$jit.lok_pos_lon, training$jit.lok_pos_lat)
          training <- SpatialPointsDataFrame(xy, training)
          xy <- cbind(validate$jit.lok_pos_lon, validate$jit.lok_pos_lat)
          validate <- SpatialPointsDataFrame(xy, validate)
          projection(training) <- projection(pnts)
          projection(validate) <- projection(pnts)
          training_sf <- st_as_sf(training)
          pl <- c("l10cp_L")
          
          # Try to perform interpolation
          res.ipdw3 <- tryCatch({
            ipdw::ipdw(training_sf, costras, range = mean.neighdist * 8 * r_mnd, pl, overlapped = TRUE, dist_power = 1.0)
          }, error = function(e) {
            message("Interpolation failed for ", spec.lat, " in ", season, " with ", mch_tp, ": ", e$message)
            NULL
          })
          
          # --- Plot with ipdw raster ---
          if (!is.null(res.ipdw3)) {
            # Convert raster to data frame for ggplot
            df_raster <- as.data.frame(res.ipdw3, xy = TRUE)
            names(df_raster) <- c("x", "y", "value")
            
            # Create ggplot
            p <- ggplot() +
              # Add local basemap
              geom_sf(
                data = coastline10,
                fill = "white",
                color = "black"
              ) +
              # Add raster layer
              geom_tile(
                data = df_raster,
                aes(x = x, y = y, fill = value),
                alpha = 0.7
              ) +
              scale_fill_gradientn(
                colors = c("white", clf.spc, "black"),
                values = scales::rescale(c(0, zvm_global / 2, zvm_global)),
                name = "log10(eDNA\ncopies/L)",
                limits = c(0, zvm_global)
              ) +
              # Add sampling points
              geom_point(
                data = df_e15_ssn,
                aes(x = lon.m, y = lat.m, color = cl.f.evl, shape = mch),
                size = 3,
                stroke = 1.5
              ) +
              scale_color_identity(
                name = "Point\nColor",
                guide = "legend"
              ) +
              scale_shape_manual(
                values = c("qPCR" = 17, "ddPCR" = 15),
                name = "Method"
              ) +
              # Crop to your area of interest
              coord_sf(
                xlim = c(6, 17),
                ylim = c(54, 59),
                expand = FALSE
              ) +
              # Add labels and BW theme
              labs(
                title = paste(
                  mch_tp, "in", spcfc_seaon_name, yr_smpl
                ),
                x = "Longitude",
                y = "Latitude"
              ) +
              theme_bw() +
              theme(
                plot.title = element_text(size = 10, hjust = 0.5),
                legend.position = "right",
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "white")
              )
            all_plots[[plot_counter]] <- p
            plot_counter <- plot_counter + 1
          } else {
            # --- Plot without ipdw raster ---
            p <- ggplot() +
              # Add local basemap
              geom_sf(
                data = coastline10,
                fill = "white",
                color = "black"
              ) +
              # Add sampling points
              geom_point(
                data = df_e15_ssn,
                aes(x = lon.m, y = lat.m, color = cl.f.evl, shape = mch),
                size = 3,
                stroke = 1.5
              ) +
              scale_color_identity(
                name = "Point\nColor",
                guide = "legend"
              ) +
              scale_shape_manual(
                values = c("qPCR" = 17, "ddPCR" = 15),
                name = "Method"
              ) +
              # Crop to your area of interest
              coord_sf(
                xlim = c(6, 17),
                ylim = c(54, 59),
                expand = FALSE
              ) +
              # Add labels and BW theme
              labs(
                title = paste(
                  mch_tp, "in", spcfc_seaon_name, yr_smpl
                ),
                x = "Longitude",
                y = "Latitude"
              ) +
              theme_bw() +
              theme(
                plot.title = element_text(size = 10, hjust = 0.5),
                legend.position = "right",
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "white")
              )
            all_plots[[plot_counter]] <- p
            plot_counter <- plot_counter + 1
          }
        } else {
          # --- Plot without ipdw raster (not enough data) ---
          p <- ggplot() +
            # Add local basemap
            geom_sf(
              data = coastline10,
              fill = "white",
              color = "black"
            ) +
            # Add sampling points
            geom_point(
              data = df_e15_ssn,
              aes(x = lon.m, y = lat.m, color = cl.f.evl, shape = mch),
              size = 3,
              stroke = 1.5
            ) +
            scale_color_identity(
              name = "Point\nColor",
              guide = "legend"
            ) +
            scale_shape_manual(
              values = c("qPCR" = 17, "ddPCR" = 15),
              name = "Method"
            ) +
            # Crop to your area of interest
            coord_sf(
              xlim = c(6, 17),
              ylim = c(54, 59),
              expand = FALSE
            ) +
            # Add labels and BW theme
            labs(
              title = paste(
                mch_tp, "in", spcfc_seaon_name, yr_smpl
              ),
              x = "Longitude",
              y = "Latitude"
            ) +
            theme_bw() +
            theme(
              plot.title = element_text(size = 10, hjust = 0.5),
              legend.position = "right",
              panel.grid = element_blank(),
              panel.background = element_rect(fill = "white")
            )
          all_plots[[plot_counter]] <- p
          plot_counter <- plot_counter + 1
        }
      }
    }
  }
  
  # Combine all plots for this species into a single layout
  if (length(all_plots) > 0) {
    # Use patchwork to arrange plots, with only one legend for each aesthetic
    combined_plot <- patchwork::wrap_plots(all_plots, ncol = 2,
                                           guides = "collect", 
                                           axis_titles="collect",
                                           axes="collect") +
      plot_layout(guides = "collect", 
                  axis_titles="collect",
                  axes="collect") &
      theme(legend.position = "right")
    
    # Save the combined plot
    flnm <- paste(
      wd00_wd13, "/Fig13_v", sbs.AssIDNo,
      "_", short.spec.lat,
      "_res", res_fac, "_mnd", r_mnd,
      ".png", sep = ""
    )
    ggsave(flnm, plot = combined_plot,
           width = 20, height = 15, units = "cm", dpi = 300)
  }
}
