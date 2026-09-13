#
# Auxiliary functions to be used for the application to
# the Weather Images data set:
# https://www.kaggle.com/datasets/jehanbhathena/weather-dataset
#

extract_rgb_channels <- function(img_path)
{
  if (!file.exists(img_path))
    stop("File does not exist: ", img_path)

  header <- readBin(img_path, what = "raw", n = 8)

  is_jpeg <- length(header) >= 3L &&
    identical(header[1:3], as.raw(c(0xff, 0xd8, 0xff)))

  is_png <- length(header) >= 8L &&
    identical(header[1:8], as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a)))

  is_bmp <- length(header) >= 2L && identical(header[1:2], as.raw(c(0x42, 0x4d)))

  if (!is_jpeg && !is_png && !is_bmp)
    stop("unrecognized image format: ", img_path)

  img_magick <- magick::image_read(img_path)

  # Use the first frame (if there are multiple frames in the image).
  img_magick <- img_magick[1]

  img_magick <- magick::image_convert(img_magick, colorspace = "sRGB")

  # Remove alpha/transparency by compositing over white.
  # Transparent pixels become white and semitransparent pixels are blended with white.
  # Flatten to combine the image and the background into one layer and to remove the alpha/transparency channel.
  # Blending with white is necessary because otherwise the color of transparent or semitransparent pixels could not be determined.

  img_magick <- magick::image_background(img_magick, color = "white", flatten = TRUE)

  # Extract Red-Green-Blue channels.

  img <- magick::image_data(img_magick, channels = "rgb")
  d <- dim(img)

  if (length(d) != 3L)
    stop("unexpected image array dimension: ", img_path)

  if (d[1L] == 3L) {
    # [channels, width, height] -> [height, width, channels]
    img <- aperm(img, c(3L, 2L, 1L))
  } else if (d[3L] == 3L) {
    img <- img # already [height, width, channels].
  } else {
    stop("failed to extract RGB channels: ", img_path)
  }

  if (length(dim(img)) != 3L || dim(img)[3L] != 3L)
    stop("failed to extract RGB channels: ", img_path)

  img_dim <- dim(img)

  img <- as.numeric(img) / 255
  dim(img) <- img_dim

  # output
  # img[,,1] # red channel
  # img[,,2] # green channel
  # img[,,3] # blue channel
  img
}

compute_rgb_histogram_features <- function(img_path)
{
  # Normalized histograms for the red, green, and blue channels.
  # Each histogram contains (for a given color) the proportion of pixels
  # falling into each of 32 intensity bins in the intensity range 0-255.
  # The output feature vector concatenates the red, green, and blue histograms.

  bins_per_channel <- 32L
  n_of_features <- 3L * bins_per_channel
  breaks <- seq(0, 256, length.out = bins_per_channel + 1L)

  img <- try(extract_rgb_channels(img_path), silent = TRUE)

  if (inherits(img, "try-error")) {
    return(structure(list(err = TRUE, msg = "read error"), class = "procerr"))
  }

  r <- as.integer(round(as.vector(img[,,1L]) * 255))
  g <- as.integer(round(as.vector(img[,,2L]) * 255))
  b <- as.integer(round(as.vector(img[,,3L]) * 255))

  # Clipping (just for safety).
  r <- pmin(255L, pmax(0L, r))
  g <- pmin(255L, pmax(0L, g))
  b <- pmin(255L, pmax(0L, b))

  hist_r <- hist(r, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts
  hist_g <- hist(g, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts
  hist_b <- hist(b, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts

  if (sum(hist_r) == 0 || sum(hist_g) == 0 || sum(hist_b) == 0)
    return(structure(list(err = TRUE, msg = "empty histogram"), class = "procerr"))

  # Histogram normalization (convert bin counts into proportions).
  hist_r <- hist_r / sum(hist_r)
  hist_g <- hist_g / sum(hist_g)
  hist_b <- hist_b / sum(hist_b)

  features <- c(hist_r, hist_g, hist_b)

  if (any(!is.finite(features)))
    return(structure(list(err = TRUE, msg = "non-finite features"), class = "procerr"))

  list(err = FALSE, features = features)
}

resize_and_fill <- function(img_path, input_dir, output_dir,
  target_width = 150L, target_height = 100L,
  exclude_level = "dataset")
{
  # Resize an image keeping the aspect ratio and filling the whole target rectangle.
  #
  # Resizing is not needed to extract features, ... resized images are used to create
  # the collage containing all the images assigned to a cluster. If all images have
  # the same size the collage looks nicer.

  rel_path <- sub(
    paste0("^", normalizePath(input_dir, mustWork = TRUE), .Platform$file.sep),
    "", normalizePath(img_path, mustWork = TRUE))

  out_path <- file.path(output_dir,
    paste0(tools::file_path_sans_ext(rel_path), ".png"))

  if (nzchar(exclude_level))
  {
    exclude_level <- paste0(.Platform$file.sep, exclude_level, .Platform$file.sep)
    out_path <- sub(exclude_level, "/", out_path, fixed = TRUE)
  }

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

  img <- try(magick::image_read(img_path), silent = TRUE)

  if (inherits(img, "try-error"))
    return(list(ok = FALSE, input = img_path, output = out_path, msg = as.character(img)))

  # For animated GIFs, keep the first frame.
  # img <- img[1]

  # Convert to RGB.

  img <- magick::image_convert(img, colorspace = "sRGB")

  # Resize
  # '^' preserves the aspect ratio, covering also the entire rectangle

  img <- magick::image_resize(img, geometry = sprintf("%dx%d^", target_width, target_height))

  # Center crop to exactly target_width x target_height.

  img <- magick::image_extent(img,
    geometry = sprintf("%dx%d", target_width, target_height),
    gravity = "center")

  ok <- try(magick::image_write(img, path = out_path, format = "png"), silent = TRUE)

  if (inherits(ok, "try-error"))
    return(list(ok = FALSE, input = img_path, output = out_path, msg = as.character(ok)))

  list(ok = TRUE, input = img_path, output = out_path, msg = NA)
}

montage_golden_layout <- function(n,
  tile_width = 150, tile_height = 100,
  gap_x = 2, gap_y = 2, phi = (1 + sqrt(5)) / 2)
{
  # The collage is created using 'montage'.
  # Choose a layout that is close to the golden ratio to
  # make the collage look esthetically nicer.

  # Dimensions of possible layouts.

  cand <- lapply(seq_len(n), function(rows) {
    cols <- ceiling(n / rows)

    # Approximate final dimensions, including tile gaps.
    width <- cols * tile_width + max(0L, cols - 1L) * gap_x
    height <- rows * tile_height + max(0L, rows - 1L) * gap_y

    ratio <- width / height
    empty_tiles <- cols * rows - n

    data.frame(rows = rows, cols = cols,
      width = width, height = height,
      ratio = ratio,
      empty_tiles = empty_tiles,
      ratio_error = abs(log(ratio / phi)))
  })

  cand <- do.call(rbind, cand)

  # Favour layouts close to the golden ratio, but avoid too many empty tiles.
  # The second term penalizes unused cells.

  cand$score <- cand$ratio_error + 0.02 * cand$empty_tiles / n

  best <- cand[which.min(cand$score),]

  list(cols = best$cols, rows = best$rows,
    width = best$width, height = best$height,
    ratio = best$ratio,
    empty_tiles = best$empty_tiles)
}
