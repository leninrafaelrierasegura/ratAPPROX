save_plotly_figure_fixed <- function(fig,
                                     dpi = 600,
                                     scale = 2,
                                     viewer_change = 1) {
  
  folder <- here::here("data_files/plotlypdf")
  
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  fig_name <- deparse(substitute(fig))
  
  file_name_pdf <- file.path(folder, paste0(fig_name, ".pdf"))
  file_name_png <- file.path(folder, paste0(fig_name, ".png"))
  
  # ---- CAMERA (correct Plotly location) ----
  cam <- fig$x$layoutAttrs[[1]]$scene$camera
  
  if (!is.null(cam)) {
    
    eye <- cam$eye
    center <- cam$center
    
    if (is.null(eye)) eye <- list(x = 1, y = 1, z = 1)
    if (is.null(center)) center <- list(x = 0, y = 0, z = 0)
    
    # vector center -> eye
    vx <- eye$x - center$x
    vy <- eye$y - center$y
    vz <- eye$z - center$z
    
    # norm (distance from center)
    norm_v <- sqrt(vx^2 + vy^2 + vz^2)
    
    if (norm_v > 0) {
      
      # normalize direction
      vx <- vx / norm_v
      vy <- vy / norm_v
      vz <- vz / norm_v
      
      # apply zoom scaling on radius
      norm_v <- norm_v * viewer_change
      
      fig$x$layoutAttrs[[1]]$scene$camera$eye <- list(
        x = center$x + vx * norm_v,
        y = center$y + vy * norm_v,
        z = center$z + vz * norm_v
      )
    }
  }
  
  # Save PDF
  plotly::save_image(
    fig,
    file_name_pdf,
    width = NULL,
    height = NULL,
    scale = scale
  )
  
  # Convert PDF to PNG
  cmd <- paste(
    "convert -density", dpi,
    "-background white -trim +repage",
    "-bordercolor white -border 10",
    shQuote(file_name_pdf),
    shQuote(file_name_png)
  )
  
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  # Copy to presentation folder
  pres_dir <- path.expand("~/Desktop/leninPresentations/data_files")
  if (!dir.exists(pres_dir)) dir.create(pres_dir, recursive = TRUE)
  
  pres_png <- file.path(pres_dir, basename(file_name_png))
  pres_pdf <- file.path(pres_dir, basename(file_name_pdf))
  
  file.copy(file_name_png, pres_png, overwrite = TRUE)
  file.copy(file_name_pdf, pres_pdf, overwrite = TRUE)
  
  cat(
    "Files saved:\n",
    "  PDF:", file_name_pdf, "\n",
    "  PNG:", file_name_png, "\n",
    "  PDF copied to:", pres_pdf, "\n",
    "  PNG copied to:", pres_png, "\n",
    "  viewer change:", viewer_change, "\n"
  )
}