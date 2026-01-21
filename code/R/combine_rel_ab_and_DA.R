library(magick)

fname_out <- "./figures/rel_ab_DA_combined.png"

fname_rel <- "./figures/rel_ab_heatmap.png"
fname_high  <- "./figures/DA_ancom_high.png"

# Read images
img_left  <- image_read(fname_rel)
img_right <- image_read(fname_high)

# Get image heights
height_left  <- image_info(img_left)$height
height_right <- image_info(img_right)$height

# Get image heights
height_left  <- image_info(img_left)$height
height_right <- image_info(img_right)$height

# If the right image is shorter, add vertical padding to center it
if (height_right < height_left) {
  # New geometry string: WIDTHxHEIGHT
  new_geometry <- paste0(image_info(img_right)$width, "x", height_left)
  
  # Add transparent padding, center the original image
  img_right <- image_extent(img_right,
                            geometry = new_geometry,
                            gravity = "center") 
}

# Combine horizontally
combined <- image_append(c(img_left, img_right))

# Save the combined image
image_write(combined, path = fname_out)
