load_channel_map <- function(path) {
  channel_map <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("raw_channel", "standard_channel")

  if (!all(required_cols %in% colnames(channel_map))) {
    stop("Channel map must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (anyDuplicated(channel_map$raw_channel)) {
    stop("Channel map contains duplicate raw_channel values.")
  }
  if (anyDuplicated(channel_map$standard_channel)) {
    stop("Channel map contains duplicate standard_channel values.")
  }

  channel_map
}

standardize_channels <- function(fs, channel_map, allow_unmapped = "Original_ID") {
  raw_channels <- colnames(fs)
  missing_channels <- setdiff(raw_channels, c(channel_map$raw_channel, allow_unmapped))

  if (length(missing_channels) > 0L) {
    stop("Channel map is missing raw channels: ", paste(missing_channels, collapse = ", "))
  }

  mapped_channels <- raw_channels
  mapped_idx <- match(raw_channels, channel_map$raw_channel)
  mapped_channels[!is.na(mapped_idx)] <- channel_map$standard_channel[mapped_idx[!is.na(mapped_idx)]]

  if (anyDuplicated(mapped_channels)) {
    stop("Channel standardization produced duplicate channels: ", paste(unique(mapped_channels[duplicated(mapped_channels)]), collapse = ", "))
  }

  colnames(fs) <- mapped_channels
  fs
}
