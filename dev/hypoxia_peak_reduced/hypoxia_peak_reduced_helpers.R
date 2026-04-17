library(ggplot2)
library(dplyr)
library(tibble)

interpolate_half_height_crossing <- function(x_left, y_left, x_right, y_right, target_y) {
  if (!all(is.finite(c(x_left, y_left, x_right, y_right, target_y)))) {
    return(NA_real_)
  }
  if (y_left == y_right) {
    return(mean(c(x_left, x_right)))
  }
  x_left + (target_y - y_left) * (x_right - x_left) / (y_right - y_left)
}

compute_peak_prominence <- function(x, y, peak_idx) {
  peak_y <- y[peak_idx]
  left_higher <- which(y[seq_len(peak_idx)] > peak_y)
  left_boundary <- if (length(left_higher) == 0) 1 else max(left_higher)
  right_higher_rel <- which(y[peak_idx:length(y)] > peak_y)
  right_boundary <- if (length(right_higher_rel) == 0) length(y) else peak_idx + min(right_higher_rel) - 1

  left_base <- min(y[left_boundary:peak_idx], na.rm = TRUE)
  right_base <- min(y[peak_idx:right_boundary], na.rm = TRUE)
  peak_y - max(left_base, right_base)
}

summarise_threshold_peak <- function(values, side = c("below", "above"), threshold = 12000, density_n = 1024) {
  side <- match.arg(side)
  values <- values[is.finite(values)]
  values <- if (identical(side, "below")) values[values < threshold] else values[values > threshold]

  if (length(values) < 10 || length(unique(values)) < 2) {
    return(tibble(
      peak_location = NA_real_,
      peak_height = NA_real_,
      half_height_left = NA_real_,
      half_height_right = NA_real_,
      peak_width_fwhm = NA_real_,
      peak_n = length(values)
    ))
  }

  density_obj <- density(values, n = density_n, from = min(values), to = max(values), na.rm = TRUE)
  x <- density_obj$x
  y <- density_obj$y
  if (length(x) < 3 || all(!is.finite(y))) {
    return(tibble(
      peak_location = NA_real_,
      peak_height = NA_real_,
      half_height_left = NA_real_,
      half_height_right = NA_real_,
      peak_width_fwhm = NA_real_,
      peak_n = length(values)
    ))
  }

  local_max_idx <- which(y >= dplyr::lag(y, default = -Inf) & y >= dplyr::lead(y, default = -Inf))
  if (length(local_max_idx) == 0) {
    local_max_idx <- which.max(y)
  }
  peak_idx <- local_max_idx[which.max(y[local_max_idx])]
  peak_x <- x[peak_idx]
  peak_y <- y[peak_idx]
  half_height <- peak_y / 2

  left_idx <- which(y[seq_len(peak_idx)] < half_height)
  if (length(left_idx) == 0) {
    left_cross <- x[1]
  } else {
    left_low_idx <- max(left_idx)
    left_high_idx <- min(left_low_idx + 1, peak_idx)
    left_cross <- interpolate_half_height_crossing(
      x[left_low_idx], y[left_low_idx],
      x[left_high_idx], y[left_high_idx],
      half_height
    )
  }

  right_range <- peak_idx:length(y)
  right_idx_rel <- which(y[right_range] < half_height)
  if (length(right_idx_rel) == 0) {
    right_cross <- x[length(x)]
  } else {
    right_low_idx <- right_range[min(right_idx_rel)] - 1
    right_low_idx <- max(right_low_idx, peak_idx)
    right_high_idx <- min(right_low_idx + 1, length(x))
    right_cross <- interpolate_half_height_crossing(
      x[right_low_idx], y[right_low_idx],
      x[right_high_idx], y[right_high_idx],
      half_height
    )
  }

  tibble(
    peak_location = peak_x,
    peak_height = peak_y,
    half_height_left = left_cross,
    half_height_right = right_cross,
    peak_width_fwhm = right_cross - left_cross,
    peak_n = length(values)
  )
}

find_threshold_peak_candidates <- function(values, side = c("below", "above"), threshold = 12000, density_n = 2048) {
  side <- match.arg(side)
  values <- values[is.finite(values)]
  values <- if (identical(side, "below")) values[values < threshold] else values[values > threshold]

  if (length(values) < 10 || length(unique(values)) < 2) {
    return(tibble(
      candidate_rank = integer(),
      peak_location = numeric(),
      peak_height = numeric(),
      prominence = numeric(),
      half_height_left = numeric(),
      half_height_right = numeric(),
      peak_width_fwhm = numeric(),
      peak_n = integer()
    ))
  }

  density_obj <- density(values, n = density_n, from = min(values), to = max(values), na.rm = TRUE)
  x <- density_obj$x
  y <- density_obj$y
  if (length(x) < 3 || all(!is.finite(y))) {
    return(tibble(
      candidate_rank = integer(),
      peak_location = numeric(),
      peak_height = numeric(),
      prominence = numeric(),
      half_height_left = numeric(),
      half_height_right = numeric(),
      peak_width_fwhm = numeric(),
      peak_n = integer()
    ))
  }

  local_max_idx <- which(y >= dplyr::lag(y, default = -Inf) & y >= dplyr::lead(y, default = -Inf))
  if (length(local_max_idx) == 0) {
    local_max_idx <- which.max(y)
  }

  cand_tbl <- bind_rows(lapply(local_max_idx, function(peak_idx) {
    peak_x <- x[peak_idx]
    peak_y <- y[peak_idx]
    half_height <- peak_y / 2
    prominence <- compute_peak_prominence(x, y, peak_idx)

    left_idx <- which(y[seq_len(peak_idx)] < half_height)
    if (length(left_idx) == 0) {
      left_cross <- x[1]
    } else {
      left_low_idx <- max(left_idx)
      left_high_idx <- min(left_low_idx + 1, peak_idx)
      left_cross <- interpolate_half_height_crossing(
        x[left_low_idx], y[left_low_idx],
        x[left_high_idx], y[left_high_idx],
        half_height
      )
    }

    right_range <- peak_idx:length(y)
    right_idx_rel <- which(y[right_range] < half_height)
    if (length(right_idx_rel) == 0) {
      right_cross <- x[length(x)]
    } else {
      right_low_idx <- right_range[min(right_idx_rel)] - 1
      right_low_idx <- max(right_low_idx, peak_idx)
      right_high_idx <- min(right_low_idx + 1, length(x))
      right_cross <- interpolate_half_height_crossing(
        x[right_low_idx], y[right_low_idx],
        x[right_high_idx], y[right_high_idx],
        half_height
      )
    }

    tibble(
      peak_location = peak_x,
      peak_height = peak_y,
      prominence = prominence,
      half_height_left = left_cross,
      half_height_right = right_cross,
      peak_width_fwhm = right_cross - left_cross,
      peak_n = length(values)
    )
  }))

  cand_tbl %>%
    arrange(desc(prominence), desc(peak_height), peak_location) %>%
    mutate(candidate_rank = row_number()) %>%
    select(candidate_rank, everything())
}

select_separated_peaks <- function(candidate_tbl, n_keep = 2, min_separation = 35000, min_prominence = 0) {
  if (nrow(candidate_tbl) == 0) {
    return(candidate_tbl)
  }

  filtered_tbl <- candidate_tbl %>%
    filter(is.finite(prominence), prominence >= min_prominence) %>%
    arrange(desc(prominence), desc(peak_height), peak_location)

  if (nrow(filtered_tbl) == 0) {
    return(filtered_tbl)
  }

  keep_idx <- integer()
  for (i in seq_len(nrow(filtered_tbl))) {
    loc_i <- filtered_tbl$peak_location[[i]]
    if (!length(keep_idx) || all(abs(loc_i - filtered_tbl$peak_location[keep_idx]) >= min_separation)) {
      keep_idx <- c(keep_idx, i)
    }
    if (length(keep_idx) >= n_keep) {
      break
    }
  }

  filtered_tbl[keep_idx, , drop = FALSE] %>%
    arrange(peak_location) %>%
    mutate(selected_rank = row_number())
}

extract_indexed_vector <- function(row_df, prefix, n) {
  vapply(seq_len(n), function(i) row_df[[sprintf("%s.%d", prefix, i)]], numeric(1))
}

prepare_peak_input <- function(project_root, threshold = 12000, tumor_peak_min_separation = 35000, tumor_peak_min_prominence = 0) {
  metadata_path <- file.path(project_root, "processed_data", "hypoxia-sum159", "sample_metadata.csv")
  dna_path <- file.path(project_root, "processed_data", "hypoxia-sum159", "filtered_dna_area_vectors.rds")

  metadata_tbl <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  dna_vectors <- readRDS(dna_path)$raw
  sample_tbl <- tibble(sample_id = names(dna_vectors), raw_dna = unname(dna_vectors))

  peak_parts <- lapply(seq_len(nrow(sample_tbl)), function(i) {
    x <- sample_tbl$raw_dna[[i]]
    below_peak <- summarise_threshold_peak(x, side = "below", threshold = threshold)
    above_peak <- summarise_threshold_peak(x, side = "above", threshold = threshold)
    tumor_candidates <- find_threshold_peak_candidates(x, side = "above", threshold = threshold)
    tumor_selected <- select_separated_peaks(
      tumor_candidates,
      n_keep = 2,
      min_separation = tumor_peak_min_separation,
      min_prominence = tumor_peak_min_prominence
    )
    sum_below <- sum(x[x < threshold], na.rm = TRUE)
    sum_above <- sum(x[x > threshold], na.rm = TRUE)
    list(
      peak_tbl = tibble(
        sample_id = sample_tbl$sample_id[[i]],
        cen_peak = below_peak$peak_location,
        g1_peak = above_peak$peak_location,
        tumor_peak_1 = if (nrow(tumor_selected) >= 1) tumor_selected$peak_location[[1]] else NA_real_,
        tumor_peak_2 = if (nrow(tumor_selected) >= 2) tumor_selected$peak_location[[2]] else NA_real_,
        sum_below = sum_below,
        sum_above = sum_above,
        ratio_above_below = sum_above / pmax(sum_below, 1),
        n_below = sum(x < threshold, na.rm = TRUE),
        n_above = sum(x > threshold, na.rm = TRUE)
      )
    )
  })
  peak_tbl <- bind_rows(lapply(peak_parts, `[[`, "peak_tbl"))

  out <- metadata_tbl %>%
    inner_join(peak_tbl, by = "sample_id") %>%
    transmute(
      sample_name = sample_id,
      condition,
      latest_match_date,
      relative_day,
      cen_peak,
      g1_peak,
      tumor_peak_1,
      tumor_peak_2,
      sum_below,
      sum_above,
      ratio_above_below,
      n_below,
      n_above
    ) %>%
    filter(
      is.finite(cen_peak),
      is.finite(g1_peak),
      is.finite(ratio_above_below),
      ratio_above_below > 0,
      !is.na(latest_match_date)
    ) %>%
    arrange(as.Date(latest_match_date), sample_name)

  ratio_scale_obj <- scale(out$ratio_above_below)
  out$date_label <- as.character(out$latest_match_date)
  out$date_id <- match(out$date_label, unique(out$date_label))
  out$ratio_scaled <- as.numeric(ratio_scale_obj)
  attr(out, "date_levels") <- unique(out$date_label)
  attr(out, "ratio_center") <- attr(ratio_scale_obj, "scaled:center")
  attr(out, "ratio_scale") <- attr(ratio_scale_obj, "scaled:scale")
  out
}

compute_fit_table <- function(input_tbl, opt_row, date_levels) {
  log_phi_raw <- extract_indexed_vector(opt_row, "log_phi_raw", length(date_levels))
  log_phi <- log_phi_raw - mean(log_phi_raw)
  phi_date <- exp(log_phi)
  names(phi_date) <- date_levels

  alpha <- opt_row$alpha[[1]]
  beta <- opt_row$beta[[1]]
  M_cen <- exp(opt_row$log_M_cen[[1]])
  R_2n <- exp(opt_row$log_R_2n[[1]])
  R_4n <- R_2n * exp(opt_row$log_R_4n_over_2n[[1]])
  rho <- exp(opt_row$log_rho[[1]])
  sigma_g1 <- opt_row$sigma_g1[[1]]
  p_4n <- opt_row$p_4n[[1]]

  fit_tbl <- input_tbl %>%
    mutate(
      phi_date = unname(phi_date[date_label]),
      R_2n = R_2n,
      R_4n = R_4n,
      ell = alpha + beta * ratio_scaled,
      u = exp(-ell),
      mu_cen = M_cen * u / (1 + u),
      distortion = (1 + u) / (1 + rho * u),
      mu_g1_2n = mu_cen * phi_date * R_2n * distortion,
      mu_g1_4n = mu_cen * phi_date * R_4n * distortion,
      log_w4 = log(p_4n) + dlnorm(g1_peak, meanlog = log(mu_g1_4n), sdlog = sigma_g1, log = TRUE),
      log_w2 = log1p(-p_4n) + dlnorm(g1_peak, meanlog = log(mu_g1_2n), sdlog = sigma_g1, log = TRUE),
      prob_4n = plogis(log_w4 - log_w2),
      assigned_state = ifelse(prob_4n >= 0.5, "4N", "2N"),
      fitted_g1 = ifelse(prob_4n >= 0.5, mu_g1_4n, mu_g1_2n),
      phi1_counterfactual_peak = ifelse(assigned_state == "4N", mu_cen * R_4n * distortion, mu_cen * R_2n * distortion)
    )

  list(
    fit_tbl = fit_tbl,
    date_tbl = tibble(
      date_label = date_levels,
      log_phi = unname(log_phi),
      phi_date = unname(phi_date)
    ),
    param_tbl = tibble(
      alpha = alpha,
      beta = beta,
      M_cen = M_cen,
      R_2n = R_2n,
      R_4n = R_4n,
      rho = rho,
      sigma_cen = opt_row$sigma_cen[[1]],
      sigma_g1 = sigma_g1,
      p_4n = p_4n
    )
  )
}

build_sample_peak_summary <- function(fit_tbl) {
  summary_tbl <- fit_tbl %>%
    transmute(
      sample_name,
      condition,
      latest_match_date,
      relative_day,
      assigned_state,
      cen_peak,
      observed_g1_peak = g1_peak,
      observed_second_peak = tumor_peak_2,
      fitted_g1_peak = fitted_g1,
      phi1_counterfactual_peak,
      phi_date,
      ratio_above_below,
      ratio_scaled,
      prob_4n,
      abs_log2_peak1_peak2_ratio = abs(log2(tumor_peak_1 / tumor_peak_2))
    )

  max_abs_log2_ratio <- suppressWarnings(max(summary_tbl$abs_log2_peak1_peak2_ratio, na.rm = TRUE))
  summary_tbl %>%
    mutate(
      is_primary_outlier = is.finite(abs_log2_peak1_peak2_ratio) &
        is.finite(max_abs_log2_ratio) &
        abs_log2_peak1_peak2_ratio == max_abs_log2_ratio
    )
}

plot_phi_counterfactual <- function(fit_tbl) {
  plot_df <- fit_tbl %>%
    transmute(
      sample_name,
      observed = g1_peak,
      fitted_phi1 = phi1_counterfactual_peak
    ) %>%
    mutate(sample_name = factor(sample_name, levels = rev(unique(sample_name))))

  long_df <- bind_rows(
    plot_df %>% transmute(sample_name, series = "Measured fluorescence", peak = observed),
    plot_df %>% transmute(sample_name, series = "Predicted if phiF = 1", peak = fitted_phi1)
  )

  ggplot(long_df, aes(x = peak, y = sample_name, color = series)) +
    geom_point(size = 2.2, alpha = 0.9) +
    labs(
      title = "Measured G0/G1 peak vs counterfactual phiF = 1 prediction",
      subtitle = "This is the checkpoint figure to revisit",
      x = "G0/G1 peak location",
      y = NULL,
      color = NULL
    ) +
    scale_color_manual(
      values = c(
        "Measured fluorescence" = "black",
        "Predicted if phiF = 1" = "#d95f02"
      )
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}
