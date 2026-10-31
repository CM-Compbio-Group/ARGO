library(ggplot2)
library(ggnewscale)
library(sf)
library(scales)

plot_roi_split <- function(object, gene, ROI, metadata) {
  # --- 1. PREPARE DATA ---
  coords <- Embeddings(object, reduction = 'X_spatial_coords')
  
  plot_data <- data.frame(
    coord_x    = metadata[Cells(object), 'coord_x'],
    Expression = FetchData(object, vars = gene, layer = "data")[[1]],
    x          = coords[, 1],
    y          = coords[, 2],
    # IMPORTANT: Change 'tissue' below if your metadata column has a different name
    domain     = object$tissue 
  )
  
  plot_data$in_roi <- plot_data$coord_x %in% ROI
  df_outside <- plot_data[!plot_data$in_roi, ]
  df_inside  <- plot_data[plot_data$in_roi, ]
  
  
  # --- 2. GENERATE EXACT TOPOLOGICAL DOMAIN BORDERS ---
  generate_domain_outlines <- function(df, spot_radius = 65, smooth_iters = 3) {
    
    # Convert spots to spatial points
    pts <- st_as_sf(df, coords = c("x", "y"))
    pts_union <- st_union(pts)
    
    # Create a tight mask around the tissue so edges don't expand to infinity
    # Adjust spot_radius up slightly if you see tiny holes inside the tissue
    tissue_mask <- st_buffer(pts_union, dist = spot_radius)
    tissue_mask <- st_buffer(tissue_mask, dist = -spot_radius * 0.3)
    
    # Create Voronoi territories around every spot so adjacent domains share 1 exact border
    voro <- st_collection_extract(st_voronoi(pts_union), "POLYGON")
    voro_sf <- st_sf(geometry = voro)
    
    # Match each Voronoi cell to its spot's tissue domain
    join_iobject <- st_intersects(voro_sf, pts)
    voro_sf$domain <- sapply(join_iobject, function(i) if(length(i) > 0) as.character(pts$domain[i[1]]) else NA)
    voro_sf <- voro_sf[!is.na(voro_sf$domain), ]
    
    # Clip territories to the outer tissue boundary
    voro_clipped <- st_intersection(voro_sf, tissue_mask)
    
    # Dissolve individual spot territories into unified anatomical domain layers
    domain_polys <- aggregate(voro_clipped["geometry"], by = list(domain = voro_clipped$domain), FUN = st_union)
    
    # Extract boundary lines and apply Chaikin's corner-cutting algorithm to make them smooth
    lines_sf <- st_cast(domain_polys, "MULTILINESTRING")
    
    chaikin_smooth <- function(mat, iters = 3) {
      for (k in seq_len(iters)) {
        n <- nrow(mat)
        if (n < 3) return(mat)
        q <- 0.75 * mat[-n, , drop = FALSE] + 0.25 * mat[-1, , drop = FALSE]
        r <- 0.25 * mat[-n, , drop = FALSE] + 0.75 * mat[-1, , drop = FALSE]
        res <- matrix(0, nrow = 2 * (n - 1), ncol = 2)
        res[seq(1, nrow(res), by = 2), ] <- q
        res[seq(2, nrow(res), by = 2), ] <- r
        mat <- res
      }
      return(mat)
    }
    
    # Smooth all line segments
    smoothed_geoms <- lapply(st_geometry(lines_sf), function(multils) {
      coords_list <- lapply(multils, function(ls) chaikin_smooth(as.matrix(ls), iters = smooth_iters))
      st_multilinestring(coords_list)
    })
    
    st_geometry(lines_sf) <- st_sfc(smoothed_geoms)
    return(lines_sf)
  }
  
  # Compute the outlines
  spatial_outlines <- generate_domain_outlines(plot_data, spot_radius = 65, smooth_iters = 1)
  
  
  # --- 3. BUILD THE PLOT ---
  p <- ggplot() +
    
    # Layer 1: Outside ROI (Background grayscale)
    geom_point(
      data = df_outside, 
      aes(x = x, y = y, color = Expression,
          size = Expression > 0,
          alpha = Expression > 0), 
      shape = 16,  # Use a borderless circle
      stroke = 0   # Force border thickness to zero
    ) +
    scale_color_gradient(
      low = "#f0f0f0", high = "gray25",
      name = "Expression, background", 
      limits = c(0, 3),
      oob = scales::squish,
      guide = guide_colorbar(
        title.position = "right",
        barheight = unit(1, "in"),
        barwidth = unit(0.08, "in"),
        frame.colour = "black",
        ticks.colour = "black",
        order = 2
      )
    ) + 
    scale_size_manual(
      values = c('FALSE' = 0.25, 'TRUE' = 1),
      guide = "none"
    ) + 
    scale_alpha_manual(
      values = c("FALSE" = 0.15, "TRUE" = 0.5), 
      guide = "none"
    ) +
    
    # Reset color scale for the ROI
    ggnewscale::new_scale_color() +
    
    # Layer 2: Inside ROI (RdPu palette + low alpha for zeros)
    geom_point(
      data = df_inside, 
      aes(x = x, y = y, color = Expression, 
          alpha = Expression > 0, 
          size=Expression > 0), 
      shape = 16,  # Use a borderless circle
      stroke = 0   # Force border thickness to zero
    ) +
    scale_color_gradientn(
      colors = c("#fde0dd", "#fa9fb5", "#c51b8a", "#49006a"), 
      limits = c(0, 3),
      name = "Expression, ROI",
      oob = scales::squish,
      guide = guide_colorbar(
        title.position = "right",
        barheight = unit(1, "in"),
        barwidth = unit(0.08, "in"),
        frame.colour = "black",
        ticks.colour = "black",
        order = 1
      )
    ) +
    scale_size_manual(
      values = c('FALSE' = 0.25, 'TRUE' = 1.5),
      guide = "none"
    ) + 
    scale_alpha_manual(
      values = c("FALSE" = 0.15, "TRUE" = 0.5), 
      guide = "none"
    ) +
    
    # Layer 3: Smooth Anatomical Domain Outlines
    geom_sf(
      data = spatial_outlines, 
      color = "grey20", 
      linewidth = 0.15,
      linetype = "dashed", # Change to "solid" if you prefer solid lines over dashed
      inherit.aes = FALSE
    ) +
    coord_sf(expand = FALSE) +
    # Layer 4: Reference-matched styling & side-by-side vertical colorbars
    theme_void(base_family = 'Arial') + 
    theme(
      legend.position = "right",
      legend.box = "horizontal", # Places the two colorbars side-by-side like your reference image
      legend.title = element_text(angle = 90, hjust = 0.5, size = 10, margin = margin(l = 5)),
      legend.text = element_text(size = 9),
      legend.spacing.y = unit(0.2, "cm")
    )
  
  return(p)
}
