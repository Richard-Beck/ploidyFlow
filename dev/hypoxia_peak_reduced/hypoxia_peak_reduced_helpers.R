library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)

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

prepare_peak_input <- function(
  project_root,
  threshold = 12000,
  tumor_peak_min_separation = 35000,
  tumor_peak_min_prominence = 0,
  ratio_source = c("fluorescence_sum", "cell_count")
) {
  ratio_source <- match.arg(ratio_source)
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
        n_above = sum(x > threshold, na.rm = TRUE),
        count_ratio_above_below = n_above / pmax(n_below, 1)
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
      n_above,
      count_ratio_above_below
    ) %>%
    add_two_peak_model_columns() %>%
    mutate(
      modeling_ratio = if (identical(ratio_source, "cell_count")) count_ratio_above_below else ratio_above_below
    ) %>%
    filter(
      is.finite(cen_peak),
      is.finite(g1_peak),
      is.finite(modeling_ratio),
      modeling_ratio > 0,
      !is.na(latest_match_date)
    ) %>%
    arrange(as.Date(latest_match_date), sample_name)

  ratio_scale_obj <- scale(out$modeling_ratio)
  out$date_label <- as.character(out$latest_match_date)
  out$date_id <- match(out$date_label, unique(out$date_label))
  out$ratio_scaled <- as.numeric(ratio_scale_obj)
  out$ratio_source <- ratio_source
  attr(out, "date_levels") <- unique(out$date_label)
  attr(out, "ratio_center") <- attr(ratio_scale_obj, "scaled:center")
  attr(out, "ratio_scale") <- attr(ratio_scale_obj, "scaled:scale")
  attr(out, "ratio_source") <- ratio_source
  out
}

add_two_peak_model_columns <- function(input_tbl) {
  ratio_tbl <- input_tbl %>%
    mutate(
      abs_log2_peak1_peak2_ratio = ifelse(
        is.finite(tumor_peak_1) & is.finite(tumor_peak_2) & tumor_peak_1 > 0 & tumor_peak_2 > 0,
        abs(log2(tumor_peak_1 / tumor_peak_2)),
        NA_real_
      )
    )

  max_abs_log2_ratio <- suppressWarnings(max(ratio_tbl$abs_log2_peak1_peak2_ratio, na.rm = TRUE))
  ratio_tbl %>%
    rowwise() %>%
    mutate(
      is_primary_outlier = is.finite(abs_log2_peak1_peak2_ratio) &&
        is.finite(max_abs_log2_ratio) &&
        abs_log2_peak1_peak2_ratio == max_abs_log2_ratio,
      dropped_outlier_peak = {
        candidate_peaks <- c(tumor_peak_1, tumor_peak_2)
        if (is_primary_outlier && sum(is.finite(candidate_peaks)) >= 2) {
          finite_idx <- which(is.finite(candidate_peaks))
          candidate_peaks[[finite_idx[[which.max(abs(candidate_peaks[finite_idx] - g1_peak))]]]]
        } else {
          NA_real_
        }
      },
      modeled_peak_1 = {
        candidate_peaks <- c(tumor_peak_1, tumor_peak_2)
        if (is_primary_outlier && is.finite(dropped_outlier_peak)) {
          drop_idx <- which(is.finite(candidate_peaks) & candidate_peaks == dropped_outlier_peak)[[1]]
          candidate_peaks[[drop_idx]] <- NA_real_
        }
        modeled_peaks <- sort(candidate_peaks[is.finite(candidate_peaks)])
        if (length(modeled_peaks) >= 1) modeled_peaks[[1]] else NA_real_
      },
      modeled_peak_2 = {
        candidate_peaks <- c(tumor_peak_1, tumor_peak_2)
        if (is_primary_outlier && is.finite(dropped_outlier_peak)) {
          drop_idx <- which(is.finite(candidate_peaks) & candidate_peaks == dropped_outlier_peak)[[1]]
          candidate_peaks[[drop_idx]] <- NA_real_
        }
        modeled_peaks <- sort(candidate_peaks[is.finite(candidate_peaks)])
        if (length(modeled_peaks) >= 2) modeled_peaks[[2]] else NA_real_
      },
      modeled_peak_count = sum(is.finite(c(modeled_peak_1, modeled_peak_2)))
    ) %>%
    ungroup()
}

build_stan_data_for_spec <- function(prepared, spec) {
  input_tbl <- prepared$input_tbl
  sigma_delta_dna_prior_scale <- if (isTRUE(spec$sample_delta_dna) && !is.null(spec$delta_dna_prior_scale)) {
    spec$delta_dna_prior_scale
  } else {
    0.1
  }

  if (isTRUE(spec$two_peak_adjacent_8n)) {
    peak_upper <- input_tbl$modeled_peak_2
    list(
      N = nrow(input_tbl),
      N_date = length(prepared$date_levels),
      date_id = input_tbl$date_id,
      x_ratio = input_tbl$ratio_scaled,
      y_cen = input_tbl$cen_peak,
      y_peak_lower = input_tbl$modeled_peak_1,
      has_upper_peak = as.integer(is.finite(peak_upper)),
      y_peak_upper = ifelse(is.finite(peak_upper), peak_upper, input_tbl$modeled_peak_1),
      sigma_delta_dna_prior_scale = sigma_delta_dna_prior_scale
    )
  } else {
    stan_data <- prepared$stan_data
    if (isTRUE(spec$sample_delta_dna)) {
      stan_data$sigma_delta_dna_prior_scale <- sigma_delta_dna_prior_scale
    }
    stan_data
  }
}

compute_fit_table <- function(input_tbl, opt_row, date_levels) {
  compute_fit_table_from_spec(input_tbl, opt_row, date_levels, get_ablation_registry()[["baseline"]])
}

compute_fit_table_from_spec <- function(input_tbl, opt_row, date_levels, spec) {
  log_phi <- if (isTRUE(spec$include_date_effect)) {
    log_phi_raw <- extract_indexed_vector(opt_row, "log_phi_raw", length(date_levels))
    log_phi_raw - mean(log_phi_raw)
  } else {
    rep(0, length(date_levels))
  }
  phi_date <- exp(log_phi)
  names(phi_date) <- date_levels

  alpha <- opt_row$alpha[[1]]
  beta <- if (isTRUE(spec$include_x_ratio_effect)) opt_row$beta[[1]] else 0
  M_cen <- exp(opt_row$log_M_cen[[1]])
  R_2n <- exp(opt_row$log_R_2n[[1]])
  R_4n <- R_2n * exp(opt_row$log_R_4n_over_2n[[1]])
  R_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) R_4n * exp(opt_row$log_R_8n_over_4n[[1]]) else NA_real_
  rho <- if (isTRUE(spec$rho_fixed_one)) 1 else exp(opt_row$log_rho[[1]])
  sigma_u <- if (isTRUE(spec$latent_u) && "sigma_u" %in% names(opt_row)) opt_row$sigma_u[[1]] else NA_real_
  sigma_delta_dna <- if (isTRUE(spec$sample_delta_dna) && "sigma_delta_dna" %in% names(opt_row)) opt_row$sigma_delta_dna[[1]] else NA_real_
  sigma_g1 <- opt_row$sigma_g1[[1]]
  p_4n <- opt_row$p_4n[[1]]
  eta <- alpha + beta * input_tbl$ratio_scaled
  log_u <- if (isTRUE(spec$latent_u) && "log_u.1" %in% names(opt_row)) {
    extract_indexed_vector(opt_row, "log_u", nrow(input_tbl))
  } else {
    -eta
  }
  delta_dna <- if (isTRUE(spec$sample_delta_dna) && "delta_dna.1" %in% names(opt_row)) {
    extract_indexed_vector(opt_row, "delta_dna", nrow(input_tbl))
  } else {
    rep(0, nrow(input_tbl))
  }
  delta_dna_multiplier <- exp(delta_dna)
  u <- exp(log_u)
  counterfactual_mode <- if (isTRUE(spec$sample_delta_dna)) "delta_dna_only" else "phi1_only"
  reference_ratio_scaled <- 0
  reference_eta <- alpha + beta * reference_ratio_scaled
  reference_u <- exp(-reference_eta)
  reference_mu_cen <- M_cen * reference_u / (1 + reference_u)
  reference_distortion <- (1 + reference_u) / (1 + rho * reference_u)
  phi1_only_counterfactual_2n <- (M_cen * u / (1 + u)) * R_2n * ((1 + u) / (1 + rho * u))
  phi1_only_counterfactual_4n <- (M_cen * u / (1 + u)) * R_4n * ((1 + u) / (1 + rho * u))
  delta_only_counterfactual_2n <- reference_mu_cen * R_2n * reference_distortion * delta_dna_multiplier
  delta_only_counterfactual_4n <- reference_mu_cen * R_4n * reference_distortion * delta_dna_multiplier
  delta_only_counterfactual_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) {
    reference_mu_cen * R_8n * reference_distortion * delta_dna_multiplier
  } else {
    rep(NA_real_, nrow(input_tbl))
  }

  fit_tbl <- input_tbl %>%
    mutate(
      ablation_name = spec$name,
      counterfactual_mode = counterfactual_mode,
      phi_date = unname(phi_date[date_label]),
      R_2n = R_2n,
      R_4n = R_4n,
      R_8n = R_8n,
      eta = eta,
      ell = eta,
      log_u = log_u,
      delta_dna = delta_dna,
      delta_dna_multiplier = delta_dna_multiplier,
      u = u,
      mu_cen = M_cen * u / (1 + u),
      distortion = (1 + u) / (1 + rho * u),
      mu_g1_2n = mu_cen * phi_date * R_2n * distortion * delta_dna_multiplier,
      mu_g1_4n = mu_cen * phi_date * R_4n * distortion * delta_dna_multiplier
    )

  if (isTRUE(spec$two_peak_adjacent_8n)) {
    fit_tbl <- fit_tbl %>%
      mutate(
        mu_g1_8n = mu_cen * phi_date * R_8n * distortion * delta_dna_multiplier,
        lower_peak_observed = modeled_peak_1,
        upper_peak_observed = modeled_peak_2,
        has_upper_peak = is.finite(modeled_peak_2),
        log_w4 = log(p_4n) +
          dlnorm(lower_peak_observed, meanlog = log(mu_g1_4n), sdlog = sigma_g1, log = TRUE) +
          ifelse(
            has_upper_peak,
            dlnorm(upper_peak_observed, meanlog = log(mu_g1_8n), sdlog = sigma_g1, log = TRUE),
            0
          ),
        log_w2 = log1p(-p_4n) +
          dlnorm(lower_peak_observed, meanlog = log(mu_g1_2n), sdlog = sigma_g1, log = TRUE) +
          ifelse(
            has_upper_peak,
            dlnorm(upper_peak_observed, meanlog = log(mu_g1_4n), sdlog = sigma_g1, log = TRUE),
            0
          ),
        prob_4n = plogis(log_w4 - log_w2),
        assigned_state = ifelse(prob_4n >= 0.5, "4N", "2N"),
        assigned_adjacent_state = ifelse(prob_4n >= 0.5, "8N", "4N"),
        fitted_g1 = ifelse(prob_4n >= 0.5, mu_g1_4n, mu_g1_2n),
        fitted_adjacent_peak = ifelse(prob_4n >= 0.5, mu_g1_8n, mu_g1_4n),
        phi1_only_counterfactual_peak = ifelse(assigned_state == "4N", phi1_only_counterfactual_4n, phi1_only_counterfactual_2n),
        delta_only_counterfactual_peak = ifelse(assigned_state == "4N", delta_only_counterfactual_4n, delta_only_counterfactual_2n),
        phi1_only_counterfactual_adjacent_peak = ifelse(assigned_state == "4N", (M_cen * u / (1 + u)) * R_8n * ((1 + u) / (1 + rho * u)), (M_cen * u / (1 + u)) * R_4n * ((1 + u) / (1 + rho * u))),
        delta_only_counterfactual_adjacent_peak = ifelse(assigned_state == "4N", delta_only_counterfactual_8n, delta_only_counterfactual_4n),
        phi1_counterfactual_peak = ifelse(counterfactual_mode == "delta_dna_only", delta_only_counterfactual_peak, phi1_only_counterfactual_peak),
        phi1_counterfactual_adjacent_peak = ifelse(counterfactual_mode == "delta_dna_only", delta_only_counterfactual_adjacent_peak, phi1_only_counterfactual_adjacent_peak)
      )
  } else {
    fit_tbl <- fit_tbl %>%
      mutate(
        log_w4 = log(p_4n) + dlnorm(g1_peak, meanlog = log(mu_g1_4n), sdlog = sigma_g1, log = TRUE),
        log_w2 = log1p(-p_4n) + dlnorm(g1_peak, meanlog = log(mu_g1_2n), sdlog = sigma_g1, log = TRUE),
        prob_4n = plogis(log_w4 - log_w2),
        assigned_state = ifelse(prob_4n >= 0.5, "4N", "2N"),
        assigned_adjacent_state = NA_character_,
        fitted_g1 = ifelse(prob_4n >= 0.5, mu_g1_4n, mu_g1_2n),
        fitted_adjacent_peak = NA_real_,
        phi1_only_counterfactual_peak = ifelse(assigned_state == "4N", phi1_only_counterfactual_4n, phi1_only_counterfactual_2n),
        delta_only_counterfactual_peak = ifelse(assigned_state == "4N", delta_only_counterfactual_4n, delta_only_counterfactual_2n),
        phi1_only_counterfactual_adjacent_peak = NA_real_,
        delta_only_counterfactual_adjacent_peak = NA_real_,
        phi1_counterfactual_peak = ifelse(counterfactual_mode == "delta_dna_only", delta_only_counterfactual_peak, phi1_only_counterfactual_peak),
        phi1_counterfactual_adjacent_peak = NA_real_
      )
  }

  list(
    fit_tbl = fit_tbl,
    date_tbl = tibble(
      date_label = date_levels,
      log_phi = unname(log_phi),
      phi_date = unname(phi_date)
    ),
    param_tbl = tibble(
      ablation_name = spec$name,
      alpha = alpha,
      beta = beta,
      M_cen = M_cen,
      R_2n = R_2n,
      R_4n = R_4n,
      R_8n = R_8n,
      rho = rho,
      sigma_u = sigma_u,
      sigma_delta_dna = sigma_delta_dna,
      sigma_cen = opt_row$sigma_cen[[1]],
      sigma_g1 = sigma_g1,
      p_4n = p_4n
    )
  )
}

get_ablation_registry <- function() {
  list(
    baseline = list(
      name = "baseline",
      stan_file = "hypoxia_peak_reduced.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = FALSE,
      latent_u = FALSE,
      sample_delta_dna = FALSE
    ),
    no_date_effect = list(
      name = "no_date_effect",
      stan_file = "hypoxia_peak_reduced_no_date.stan",
      include_date_effect = FALSE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = FALSE,
      latent_u = FALSE,
      sample_delta_dna = FALSE
    ),
    no_x_ratio_effect = list(
      name = "no_x_ratio_effect",
      stan_file = "hypoxia_peak_reduced_no_xratio.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = FALSE,
      rho_fixed_one = FALSE,
      latent_u = FALSE,
      sample_delta_dna = FALSE
    ),
    rho_fixed_one = list(
      name = "rho_fixed_one",
      stan_file = "hypoxia_peak_reduced_rho1.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = TRUE,
      latent_u = FALSE,
      sample_delta_dna = FALSE
    )
  )
}

get_extension_registry <- function() {
  list(
    latent_u = list(
      name = "latent_u",
      stan_file = "hypoxia_peak_reduced_latent_u.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = FALSE,
      latent_u = TRUE,
      sample_delta_dna = FALSE,
      two_peak_adjacent_8n = FALSE
    ),
    delta_dna = list(
      name = "delta_dna",
      stan_file = "hypoxia_peak_reduced_delta_dna.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = FALSE,
      latent_u = FALSE,
      sample_delta_dna = TRUE,
      two_peak_adjacent_8n = FALSE
    ),
    adjacent_8n = list(
      name = "adjacent_8n",
      stan_file = "hypoxia_peak_reduced_adjacent_8n.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = FALSE,
      latent_u = FALSE,
      sample_delta_dna = FALSE,
      two_peak_adjacent_8n = TRUE
    ),
    adjacent_8n_delta_dna = list(
      name = "adjacent_8n_delta_dna",
      stan_file = "hypoxia_peak_reduced_adjacent_8n_delta_dna.stan",
      include_date_effect = TRUE,
      include_x_ratio_effect = TRUE,
      rho_fixed_one = FALSE,
      latent_u = FALSE,
      sample_delta_dna = TRUE,
      two_peak_adjacent_8n = TRUE
    )
  )
}

build_ablation_init <- function(input_tbl, date_levels, spec) {
  alpha_init <- 0
  beta_init <- if (isTRUE(spec$include_x_ratio_effect)) 0.6 else 0
  eta_init <- alpha_init + beta_init * input_tbl$ratio_scaled
  init <- list(
    alpha = alpha_init,
    log_M_cen = log(max(input_tbl$cen_peak) * 2),
    log_R_2n = log(8),
    log_R_4n_over_2n = log(2),
    sigma_cen = 0.08,
    sigma_g1 = 0.08,
    p_4n = mean(grepl("^4N", input_tbl$condition))
  )

  if (isTRUE(spec$two_peak_adjacent_8n)) {
    init$log_R_8n_over_4n <- log(2)
  }

  if (isTRUE(spec$latent_u)) {
    init$log_u <- eta_init
    init$sigma_u <- 0.15
  }
  if (isTRUE(spec$sample_delta_dna)) {
    init$delta_dna <- rep(0, nrow(input_tbl))
    init$sigma_delta_dna <- 0.03
  }

  if (isTRUE(spec$include_x_ratio_effect)) {
    init$beta <- beta_init
  }
  if (!isTRUE(spec$rho_fixed_one)) {
    init$log_rho <- 0
  }
  if (isTRUE(spec$include_date_effect)) {
    init$log_phi_raw <- rep(0, length(date_levels))
  }

  init
}

run_ablation_fit <- function(project_root, prepared, spec, seed = 123, iter = 4000, use_primary_output = FALSE, output_category = "ablations") {
  dev_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced")
  ablation_out_dir <- if (isTRUE(use_primary_output) && identical(spec$name, "baseline")) {
    file.path(dev_dir, "output")
  } else {
    file.path(dev_dir, "output", output_category, spec$name)
  }
  dir.create(ablation_out_dir, recursive = TRUE, showWarnings = FALSE)

  input_tbl <- prepared$input_tbl
  date_levels <- prepared$date_levels
  stan_data <- build_stan_data_for_spec(prepared, spec)
  init <- build_ablation_init(input_tbl, date_levels, spec)

  model_stub <- tools::file_path_sans_ext(spec$stan_file)
  mod <- cmdstanr::cmdstan_model(
    stan_file = file.path(dev_dir, spec$stan_file),
    exe_file = file.path(ablation_out_dir, paste0(model_stub, ".exe")),
    compile = TRUE,
    force_recompile = FALSE
  )

  fit <- mod$optimize(
    data = stan_data,
    init = list(init),
    seed = seed,
    iter = iter,
    algorithm = "lbfgs",
    output_dir = ablation_out_dir,
    output_basename = paste0(model_stub, "_opt")
  )

  opt_row <- read.csv(fit$output_files(), comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
  fit_parts <- compute_fit_table_from_spec(
    input_tbl = input_tbl,
    opt_row = opt_row,
    date_levels = date_levels,
    spec = spec
  )

  list(
    spec = spec,
    out_dir = ablation_out_dir,
    fit = fit,
    opt_row = opt_row,
    fit_parts = fit_parts
  )
}

run_ablation_nuts <- function(
  project_root,
  prepared,
  spec,
  seed = 123,
  chains = 4,
  parallel_chains = chains,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.9,
  max_treedepth = 12,
  output_category = "extensions"
) {
  dev_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced")
  fit_out_dir <- file.path(dev_dir, "output", output_category, spec$name)
  nuts_out_dir <- file.path(fit_out_dir, "nuts")
  dir.create(nuts_out_dir, recursive = TRUE, showWarnings = FALSE)

  input_tbl <- prepared$input_tbl
  date_levels <- prepared$date_levels
  stan_data <- build_stan_data_for_spec(prepared, spec)
  init <- build_ablation_init(input_tbl, date_levels, spec)

  model_stub <- tools::file_path_sans_ext(spec$stan_file)
  mod <- cmdstanr::cmdstan_model(
    stan_file = file.path(dev_dir, spec$stan_file),
    exe_file = file.path(fit_out_dir, paste0(model_stub, "_nuts.exe")),
    compile = TRUE,
    force_recompile = FALSE
  )

  fit <- mod$sample(
    data = stan_data,
    init = replicate(chains, init, simplify = FALSE),
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    refresh = 100,
    output_dir = nuts_out_dir,
    output_basename = paste0(model_stub, "_nuts")
  )

  draws_matrix <- fit$draws(format = "matrix")
  summary_tbl <- fit$summary()
  metadata <- list(
    generated_at = Sys.time(),
    spec = spec,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    stan_csv_files = normalizePath(fit$output_files(), winslash = "/", mustWork = TRUE),
    draws_matrix_rds = normalizePath(file.path(nuts_out_dir, "draws_matrix.rds"), winslash = "/", mustWork = FALSE),
    summary_csv = normalizePath(file.path(nuts_out_dir, "draws_summary.csv"), winslash = "/", mustWork = FALSE)
  )

  saveRDS(draws_matrix, file.path(nuts_out_dir, "draws_matrix.rds"))
  write.csv(summary_tbl, file.path(nuts_out_dir, "draws_summary.csv"), row.names = FALSE)
  saveRDS(metadata, file.path(nuts_out_dir, "nuts_manifest.rds"))

  list(
    spec = spec,
    out_dir = fit_out_dir,
    nuts_out_dir = nuts_out_dir,
    fit = fit,
    draws_matrix = draws_matrix,
    summary_tbl = summary_tbl,
    metadata = metadata
  )
}

extract_indexed_draw_matrix <- function(draws_matrix, prefix, n) {
  col_tbl <- vapply(seq_len(n), function(i) sprintf("%s[%d]", prefix, i), character(1))
  missing_cols <- setdiff(col_tbl, colnames(draws_matrix))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing draw columns: %s", paste(missing_cols, collapse = ", ")))
  }
  as.matrix(draws_matrix[, col_tbl, drop = FALSE])
}

build_delta_only_counterfactual_draws <- function(draws_matrix, prepared, spec) {
  if (!isTRUE(spec$sample_delta_dna)) {
    stop("Delta-only counterfactual reconstruction requires `sample_delta_dna = TRUE`.")
  }

  input_tbl <- prepared$input_tbl
  n_draw <- nrow(draws_matrix)
  n_sample <- nrow(input_tbl)

  alpha <- draws_matrix[, "alpha"]
  beta <- if (isTRUE(spec$include_x_ratio_effect) && "beta" %in% colnames(draws_matrix)) draws_matrix[, "beta"] else rep(0, n_draw)
  log_M_cen <- draws_matrix[, "log_M_cen"]
  log_R_2n <- draws_matrix[, "log_R_2n"]
  log_R_4n_over_2n <- draws_matrix[, "log_R_4n_over_2n"]
  log_R_8n_over_4n <- if (isTRUE(spec$two_peak_adjacent_8n) && "log_R_8n_over_4n" %in% colnames(draws_matrix)) draws_matrix[, "log_R_8n_over_4n"] else rep(NA_real_, n_draw)
  log_rho <- if (!isTRUE(spec$rho_fixed_one) && "log_rho" %in% colnames(draws_matrix)) draws_matrix[, "log_rho"] else rep(0, n_draw)
  delta_dna <- extract_indexed_draw_matrix(draws_matrix, "delta_dna", n_sample)
  prob_4n <- extract_indexed_draw_matrix(draws_matrix, "prob_4n", n_sample)

  M_cen <- exp(log_M_cen)
  R_2n <- exp(log_R_2n)
  R_4n <- R_2n * exp(log_R_4n_over_2n)
  R_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) R_4n * exp(log_R_8n_over_4n) else rep(NA_real_, n_draw)
  rho <- if (isTRUE(spec$rho_fixed_one)) rep(1, n_draw) else exp(log_rho)

  reference_eta <- alpha + beta * 0
  reference_u <- exp(-reference_eta)
  reference_mu_cen <- M_cen * reference_u / (1 + reference_u)
  reference_distortion <- (1 + reference_u) / (1 + rho * reference_u)
  base_2n <- reference_mu_cen * R_2n * reference_distortion
  base_4n <- reference_mu_cen * R_4n * reference_distortion
  base_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) reference_mu_cen * R_8n * reference_distortion else rep(NA_real_, n_draw)

  delta_dna_multiplier <- exp(delta_dna)
  delta_only_counterfactual_2n <- delta_dna_multiplier * base_2n
  delta_only_counterfactual_4n <- delta_dna_multiplier * base_4n
  delta_only_counterfactual_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) delta_dna_multiplier * base_8n else matrix(NA_real_, nrow = n_draw, ncol = n_sample)
  assigned_state <- ifelse(prob_4n >= 0.5, "4N", "2N")
  delta_only_counterfactual_peak <- ifelse(assigned_state == "4N", delta_only_counterfactual_4n, delta_only_counterfactual_2n)
  delta_only_counterfactual_peak_mixture <- prob_4n * delta_only_counterfactual_4n + (1 - prob_4n) * delta_only_counterfactual_2n
  delta_only_counterfactual_adjacent_peak <- if (isTRUE(spec$two_peak_adjacent_8n)) {
    ifelse(assigned_state == "4N", delta_only_counterfactual_8n, delta_only_counterfactual_4n)
  } else {
    matrix(NA_real_, nrow = n_draw, ncol = n_sample)
  }
  delta_only_counterfactual_adjacent_peak_mixture <- if (isTRUE(spec$two_peak_adjacent_8n)) {
    prob_4n * delta_only_counterfactual_8n + (1 - prob_4n) * delta_only_counterfactual_4n
  } else {
    matrix(NA_real_, nrow = n_draw, ncol = n_sample)
  }

  draw_tbl <- tibble(
    draw_id = rep(seq_len(n_draw), times = n_sample),
    sample_index = rep(seq_len(n_sample), each = n_draw),
    sample_name = rep(input_tbl$sample_name, each = n_draw),
    condition = rep(input_tbl$condition, each = n_draw),
    latest_match_date = rep(as.character(input_tbl$latest_match_date), each = n_draw),
    relative_day = rep(input_tbl$relative_day, each = n_draw),
    ratio_above_below = rep(input_tbl$ratio_above_below, each = n_draw),
    delta_dna = as.vector(delta_dna),
    delta_dna_multiplier = as.vector(delta_dna_multiplier),
    prob_4n = as.vector(prob_4n),
    assigned_state = as.vector(assigned_state),
    delta_only_counterfactual_2n = as.vector(delta_only_counterfactual_2n),
    delta_only_counterfactual_4n = as.vector(delta_only_counterfactual_4n),
    delta_only_counterfactual_8n = as.vector(delta_only_counterfactual_8n),
    delta_only_counterfactual_peak = as.vector(delta_only_counterfactual_peak),
    delta_only_counterfactual_peak_mixture = as.vector(delta_only_counterfactual_peak_mixture),
    delta_only_counterfactual_adjacent_peak = as.vector(delta_only_counterfactual_adjacent_peak),
    delta_only_counterfactual_adjacent_peak_mixture = as.vector(delta_only_counterfactual_adjacent_peak_mixture)
  )

  summary_tbl <- draw_tbl %>%
    group_by(sample_index, sample_name, condition, latest_match_date, relative_day, ratio_above_below) %>%
    summarise(
      posterior_prob_4n = mean(prob_4n, na.rm = TRUE),
      assigned_state_map = ifelse(mean(prob_4n, na.rm = TRUE) >= 0.5, "4N", "2N"),
      mean_peak = mean(delta_only_counterfactual_peak, na.rm = TRUE),
      median_peak = stats::median(delta_only_counterfactual_peak, na.rm = TRUE),
      q05_peak = stats::quantile(delta_only_counterfactual_peak, probs = 0.05, na.rm = TRUE),
      q25_peak = stats::quantile(delta_only_counterfactual_peak, probs = 0.25, na.rm = TRUE),
      q75_peak = stats::quantile(delta_only_counterfactual_peak, probs = 0.75, na.rm = TRUE),
      q95_peak = stats::quantile(delta_only_counterfactual_peak, probs = 0.95, na.rm = TRUE),
      mean_peak_mixture = mean(delta_only_counterfactual_peak_mixture, na.rm = TRUE),
      median_peak_mixture = stats::median(delta_only_counterfactual_peak_mixture, na.rm = TRUE),
      mean_adjacent_peak = mean(delta_only_counterfactual_adjacent_peak, na.rm = TRUE),
      median_adjacent_peak = stats::median(delta_only_counterfactual_adjacent_peak, na.rm = TRUE),
      q05_adjacent_peak = stats::quantile(delta_only_counterfactual_adjacent_peak, probs = 0.05, na.rm = TRUE),
      q25_adjacent_peak = stats::quantile(delta_only_counterfactual_adjacent_peak, probs = 0.25, na.rm = TRUE),
      q75_adjacent_peak = stats::quantile(delta_only_counterfactual_adjacent_peak, probs = 0.75, na.rm = TRUE),
      q95_adjacent_peak = stats::quantile(delta_only_counterfactual_adjacent_peak, probs = 0.95, na.rm = TRUE),
      mean_adjacent_peak_mixture = mean(delta_only_counterfactual_adjacent_peak_mixture, na.rm = TRUE),
      median_adjacent_peak_mixture = stats::median(delta_only_counterfactual_adjacent_peak_mixture, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(median_peak, median_adjacent_peak) %>%
    mutate(sample_label = factor(sample_name, levels = sample_name))

  list(
    draw_tbl = draw_tbl,
    summary_tbl = summary_tbl
  )
}

plot_delta_only_counterfactual_posterior <- function(draw_tbl, summary_tbl) {
  plot_df <- draw_tbl %>%
    inner_join(summary_tbl %>% select(sample_name, sample_label, assigned_state_map), by = "sample_name") %>%
    mutate(sample_label = factor(sample_label, levels = levels(summary_tbl$sample_label)))

  lower_plot_df <- plot_df %>%
    transmute(
      sample_label,
      assigned_state_map,
      peak_role = "Lower peak",
      peak_value = delta_only_counterfactual_peak
    )
  adjacent_plot_df <- plot_df %>%
    transmute(
      sample_label,
      assigned_state_map,
      peak_role = "Adjacent upper peak",
      peak_value = delta_only_counterfactual_adjacent_peak
    ) %>%
    filter(is.finite(peak_value))
  combined_plot_df <- bind_rows(lower_plot_df, adjacent_plot_df)

  lower_summary_df <- summary_tbl %>%
    transmute(
      sample_label,
      peak_role = "Lower peak",
      median_value = median_peak,
      q05_value = q05_peak,
      q95_value = q95_peak,
      point_shape = "Lower peak"
    )
  adjacent_summary_df <- summary_tbl %>%
    transmute(
      sample_label,
      peak_role = "Adjacent upper peak",
      median_value = median_adjacent_peak,
      q05_value = q05_adjacent_peak,
      q95_value = q95_adjacent_peak,
      point_shape = "Adjacent upper peak"
    ) %>%
    filter(is.finite(median_value))
  combined_summary_df <- bind_rows(lower_summary_df, adjacent_summary_df)

  ggplot(combined_plot_df, aes(x = peak_value, y = sample_label)) +
    geom_violin(
      aes(fill = peak_role),
      scale = "width",
      trim = FALSE,
      color = "grey35",
      linewidth = 0.2,
      alpha = 0.45,
      position = position_identity()
    ) +
    geom_point(
      data = combined_summary_df,
      aes(x = median_value, y = sample_label, shape = point_shape),
      inherit.aes = FALSE,
      color = "black",
      size = 1.2
    ) +
    geom_segment(
      data = combined_summary_df,
      aes(x = q05_value, xend = q95_value, y = sample_label, yend = sample_label),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.35
    ) +
    scale_fill_manual(values = c("Lower peak" = "#3182bd", "Adjacent upper peak" = "#de2d26")) +
    scale_shape_manual(values = c("Lower peak" = 16, "Adjacent upper peak" = 17)) +
    labs(
      title = "Posterior distributions for the delta-only counterfactual peaks",
      subtitle = "Lower and adjacent upper counterfactual peaks are overlaid within each sample.",
      x = "Counterfactual peak location",
      y = NULL,
      fill = "Peak role",
      shape = "Peak role"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

write_delta_only_counterfactual_outputs <- function(project_root, spec_name = "delta_dna", output_category = "extensions") {
  spec_registry <- c(get_ablation_registry(), get_extension_registry())
  if (!spec_name %in% names(spec_registry)) {
    stop(sprintf("Unknown spec `%s`.", spec_name))
  }

  spec <- spec_registry[[spec_name]]
  nuts_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", output_category, spec$name, "nuts")
  draws_path <- file.path(nuts_dir, "draws_matrix.rds")
  prepared_path <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", "prepared_input.rds")

  draws_matrix <- readRDS(draws_path)
  prepared <- readRDS(prepared_path)
  posterior_parts <- build_delta_only_counterfactual_draws(draws_matrix, prepared, spec)

  draw_rds_path <- file.path(nuts_dir, "delta_only_counterfactual_draws.rds")
  summary_csv_path <- file.path(nuts_dir, "delta_only_counterfactual_summary.csv")
  plot_path <- file.path(nuts_dir, "delta_only_counterfactual_posterior.png")

  saveRDS(posterior_parts$draw_tbl, draw_rds_path)
  write.csv(posterior_parts$summary_tbl, summary_csv_path, row.names = FALSE)
  ggplot2::ggsave(plot_path, plot_delta_only_counterfactual_posterior(posterior_parts$draw_tbl, posterior_parts$summary_tbl), width = 10, height = 8, dpi = 220)

  list(
    spec = spec,
    nuts_dir = nuts_dir,
    draw_rds_path = draw_rds_path,
    summary_csv_path = summary_csv_path,
    plot_path = plot_path,
    draw_tbl = posterior_parts$draw_tbl,
    summary_tbl = posterior_parts$summary_tbl
  )
}

write_fit_outputs <- function(run_parts) {
  out_dir <- run_parts$out_dir
  fit_tbl <- run_parts$fit_parts$fit_tbl
  date_tbl <- run_parts$fit_parts$date_tbl
  param_tbl <- run_parts$fit_parts$param_tbl %>%
    mutate(
      objective_lp = if ("lp__" %in% names(run_parts$opt_row)) run_parts$opt_row$lp__[[1]] else NA_real_
    )

  write.csv(fit_tbl, file.path(out_dir, "sample_predictions.csv"), row.names = FALSE)
  write.csv(date_tbl, file.path(out_dir, "date_phi.csv"), row.names = FALSE)
  write.csv(param_tbl, file.path(out_dir, "fit_parameters.csv"), row.names = FALSE)

  invisible(out_dir)
}

render_fit_outputs <- function(out_dir) {
  fit_tbl <- read.csv(file.path(out_dir, "sample_predictions.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  summary_tbl <- build_sample_peak_summary(fit_tbl)

  write.csv(summary_tbl, file.path(out_dir, "sample_peak_summary.csv"), row.names = FALSE)
  ggplot2::ggsave(file.path(out_dir, "phi_counterfactual.png"), plot_phi_counterfactual(fit_tbl), width = 12, height = 8, dpi = 220)

  summary_tbl
}

summarise_ablation_run <- function(run_parts) {
  fit_tbl <- run_parts$fit_parts$fit_tbl
  spec <- run_parts$spec
  observed_lower_peak <- if ("lower_peak_observed" %in% names(fit_tbl)) fit_tbl$lower_peak_observed else fit_tbl$g1_peak
  n_date_levels <- dplyr::n_distinct(fit_tbl$date_label)
  free_parameter_count <- 7 +
    (if (isTRUE(spec$include_x_ratio_effect)) 1 else 0) +
    (if (!isTRUE(spec$rho_fixed_one)) 1 else 0) +
    (if (isTRUE(spec$include_date_effect)) n_date_levels else 0) +
    (if (isTRUE(spec$latent_u)) nrow(fit_tbl) + 1 else 0) +
    (if (isTRUE(spec$sample_delta_dna)) nrow(fit_tbl) + 1 else 0) +
    (if (isTRUE(spec$two_peak_adjacent_8n)) 1 else 0)

  tibble(
    ablation_name = spec$name,
    stan_file = spec$stan_file,
    include_date_effect = spec$include_date_effect,
    include_x_ratio_effect = spec$include_x_ratio_effect,
    rho_fixed_one = spec$rho_fixed_one,
    free_parameter_count = free_parameter_count,
    n_samples = nrow(fit_tbl),
    objective_lp = if ("lp__" %in% names(run_parts$opt_row)) run_parts$opt_row$lp__[[1]] else NA_real_,
    rmse_cen = sqrt(mean((fit_tbl$cen_peak - fit_tbl$mu_cen)^2, na.rm = TRUE)),
    rmse_g1 = sqrt(mean((observed_lower_peak - fit_tbl$fitted_g1)^2, na.rm = TRUE)),
    cor_counterfactual_vs_burden = suppressWarnings(stats::cor(fit_tbl$phi1_counterfactual_peak, fit_tbl$ratio_above_below, use = "complete.obs")),
    mean_phi_date = mean(fit_tbl$phi_date, na.rm = TRUE),
    sd_phi_date = stats::sd(fit_tbl$phi_date, na.rm = TRUE)
  )
}

build_sample_peak_summary <- function(fit_tbl) {
  lower_peak_observed <- if ("lower_peak_observed" %in% names(fit_tbl)) fit_tbl$lower_peak_observed else fit_tbl$g1_peak
  upper_peak_observed <- if ("upper_peak_observed" %in% names(fit_tbl)) fit_tbl$upper_peak_observed else fit_tbl$tumor_peak_2
  summary_tbl <- fit_tbl %>%
    transmute(
      sample_name,
      condition,
      latest_match_date,
      relative_day,
      assigned_state,
      cen_peak,
      observed_g1_peak = lower_peak_observed,
      observed_second_peak = upper_peak_observed,
      fitted_g1_peak = fitted_g1,
      phi1_counterfactual_peak,
      phi_date,
      ratio_above_below,
      count_ratio_above_below,
      gated_n_below = n_below,
      gated_n_above = n_above,
      ratio_source,
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
  has_adjacent_peak <- all(c("upper_peak_observed", "fitted_adjacent_peak", "assigned_adjacent_state") %in% names(fit_tbl))
  lower_peak_observed <- if ("lower_peak_observed" %in% names(fit_tbl)) fit_tbl$lower_peak_observed else fit_tbl$g1_peak
  upper_peak_observed <- if ("upper_peak_observed" %in% names(fit_tbl)) fit_tbl$upper_peak_observed else rep(NA_real_, nrow(fit_tbl))
  fitted_adjacent_peak <- if ("fitted_adjacent_peak" %in% names(fit_tbl)) fit_tbl$fitted_adjacent_peak else rep(NA_real_, nrow(fit_tbl))
  assigned_adjacent_state <- if ("assigned_adjacent_state" %in% names(fit_tbl)) fit_tbl$assigned_adjacent_state else rep(NA_character_, nrow(fit_tbl))
  counterfactual_mode <- if ("counterfactual_mode" %in% names(fit_tbl)) fit_tbl$counterfactual_mode else rep("phi1_only", nrow(fit_tbl))
  counterfactual_series_label <- if (all(counterfactual_mode == "delta_dna_only")) {
    "Predicted with only delta_dna varying"
  } else {
    "Predicted if phiF = 1"
  }
  phi1_adjacent_peak <- if ("phi1_counterfactual_adjacent_peak" %in% names(fit_tbl)) fit_tbl$phi1_counterfactual_adjacent_peak else rep(NA_real_, nrow(fit_tbl))

  plot_df <- fit_tbl %>%
    transmute(
      sample_name,
      latest_match_date = as.character(latest_match_date),
      assigned_state = as.character(assigned_state),
      assigned_adjacent_state = as.character(assigned_adjacent_state),
      observed_g1 = lower_peak_observed,
      observed_adjacent = upper_peak_observed,
      fitted_g1 = fitted_g1,
      fitted_adjacent = fitted_adjacent_peak,
      fitted_phi1_g1 = phi1_counterfactual_peak,
      fitted_phi1_adjacent = phi1_adjacent_peak,
      observed_cen = cen_peak,
      fitted_cen = mu_cen
    )

  sample_levels <- rev(unique(plot_df$sample_name))
  group_levels <- unique(plot_df$latest_match_date)
  group_colors <- setNames(grDevices::hcl.colors(length(group_levels), palette = "Dark 3"), group_levels)
  axis_label_map <- plot_df %>%
    distinct(sample_name, latest_match_date) %>%
    mutate(
      sample_name = factor(sample_name, levels = sample_levels),
      axis_label = sprintf("<span style='color:%s;'>%s</span>", group_colors[latest_match_date], sample_name)
    )
  axis_labels <- setNames(as.character(axis_label_map$axis_label), as.character(axis_label_map$sample_name))

  long_df <- bind_rows(
    plot_df %>% transmute(sample_name, peak_role = "CEN", state_label = "CEN", series = "Measured", peak = observed_cen),
    plot_df %>% transmute(sample_name, peak_role = "CEN", state_label = "CEN", series = "Predicted at fitted phiF", peak = fitted_cen),
    plot_df %>% transmute(sample_name, peak_role = "Lower tumor peak", state_label = assigned_state, series = "Measured", peak = observed_g1),
    plot_df %>% transmute(sample_name, peak_role = "Lower tumor peak", state_label = assigned_state, series = "Predicted at fitted phiF", peak = fitted_g1),
    plot_df %>% transmute(sample_name, peak_role = "Lower tumor peak", state_label = assigned_state, series = counterfactual_series_label, peak = fitted_phi1_g1)
  ) %>%
    bind_rows(
      plot_df %>%
        filter(is.finite(observed_adjacent), !is.na(assigned_adjacent_state)) %>%
        transmute(sample_name, peak_role = "Upper tumor peak", state_label = assigned_adjacent_state, series = "Measured", peak = observed_adjacent),
      plot_df %>%
        filter(is.finite(fitted_adjacent), !is.na(assigned_adjacent_state)) %>%
        transmute(sample_name, peak_role = "Upper tumor peak", state_label = assigned_adjacent_state, series = "Predicted at fitted phiF", peak = fitted_adjacent),
      plot_df %>%
        filter(is.finite(fitted_phi1_adjacent), !is.na(assigned_adjacent_state)) %>%
        transmute(sample_name, peak_role = "Upper tumor peak", state_label = assigned_adjacent_state, series = counterfactual_series_label, peak = fitted_phi1_adjacent)
    ) %>%
    filter(is.finite(peak)) %>%
    mutate(
      peak_type = ifelse(state_label == "CEN", "CEN", sprintf("%s (%s)", peak_role, state_label)),
      peak_type = factor(
        peak_type,
        levels = c(
          "CEN",
          "Lower tumor peak (2N)",
          "Lower tumor peak (4N)",
          "Upper tumor peak (4N)",
          "Upper tumor peak (8N)"
        )
      )
    ) %>%
    mutate(sample_name = factor(sample_name, levels = sample_levels))

  ggplot(long_df, aes(x = peak, y = sample_name, color = series, shape = peak_type)) +
    geom_point(position = position_jitter(width = 0, height = 0.05), size = 2.4, alpha = 0.9) +
    labs(
      title = "Measured and predicted peak locations across fitted and counterfactual settings",
      subtitle = "All modeled tumor peaks are shown; axis-label colors mark samples collected in the same date group.",
      x = "Peak location",
      y = NULL,
      color = NULL,
      shape = NULL
    ) +
    scale_color_manual(
      values = c(
        "Measured" = "black",
        "Predicted at fitted phiF" = "#1b9e77",
        "Predicted if phiF = 1" = "#d95f02",
        "Predicted with only delta_dna varying" = "#d95f02"
      )
    ) +
    scale_shape_manual(
      values = c(
        "CEN" = 17,
        "Lower tumor peak (2N)" = 16,
        "Lower tumor peak (4N)" = 15,
        "Upper tumor peak (4N)" = 1,
        "Upper tumor peak (8N)" = 2
      ),
      drop = FALSE
    ) +
    scale_y_discrete(labels = axis_labels) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.y = ggtext::element_markdown()
    )
}

parse_series_id <- function(x) {
  lineage_match <- regexec("([24]N)(?:_(C|O[12]))?(?:_A([0-9]+))?", x)
  lineage_parts <- regmatches(x, lineage_match)

  tibble(
    raw_id = x,
    ploidy_label = vapply(lineage_parts, function(v) if (length(v) >= 2) v[[2]] else NA_character_, character(1)),
    lineage_label = vapply(lineage_parts, function(v) if (length(v) >= 3 && nzchar(v[[3]])) v[[3]] else NA_character_, character(1)),
    passage = suppressWarnings(as.integer(vapply(lineage_parts, function(v) if (length(v) >= 4) v[[4]] else NA_character_, character(1))))
  )
}

compute_karyotype_ploidy_metrics <- function(karyotype_vec) {
  autosome_mbp <- c(
    chr1 = 249.0, chr2 = 242.2, chr3 = 198.3, chr4 = 190.2,
    chr5 = 181.5, chr6 = 170.8, chr7 = 159.3, chr8 = 145.1,
    chr9 = 138.4, chr10 = 133.8, chr11 = 135.1, chr12 = 133.3,
    chr13 = 114.4, chr14 = 107.0, chr15 = 102.0, chr16 = 90.3,
    chr17 = 83.3, chr18 = 80.4, chr19 = 58.6, chr20 = 64.4,
    chr21 = 46.7, chr22 = 50.8
  )
  unidentified_mbp <- mean(autosome_mbp)

  chrom_labels <- sub(":.*$", "", names(karyotype_vec))
  known_mask <- chrom_labels %in% as.character(seq_len(22))
  unidentified_mask <- chrom_labels == "999"
  supplied_unidentified_mask <- chrom_labels == "9999"

  known_copy_numbers <- as.numeric(karyotype_vec[known_mask])
  known_weights <- autosome_mbp[paste0("chr", chrom_labels[known_mask])]
  known_contrib <- sum(known_copy_numbers * known_weights, na.rm = TRUE)
  known_weight_total <- sum(known_weights, na.rm = TRUE)

  unidentified_copy_number <- if (any(unidentified_mask)) as.numeric(karyotype_vec[unidentified_mask][[1]]) else 0
  supplied_unidentified_dna_ploidy <- if (any(supplied_unidentified_mask)) as.numeric(karyotype_vec[supplied_unidentified_mask][[1]]) else 0
  has_unidentified_component <- any(unidentified_mask) || any(supplied_unidentified_mask)
  total_weight <- known_weight_total + if (has_unidentified_component) unidentified_mbp else 0

  dna_ploidy_avg_unidentified <- if (total_weight > 0) {
    (known_contrib + unidentified_copy_number * unidentified_mbp) / total_weight
  } else {
    NA_real_
  }

  dna_ploidy_supplied_unidentified <- if (total_weight > 0) {
    known_contrib / total_weight + supplied_unidentified_dna_ploidy
  } else {
    NA_real_
  }

  tibble(
    numerical_ploidy_excluding_unidentified = sum(known_copy_numbers, na.rm = TRUE),
    numerical_ploidy_including_unidentified = sum(known_copy_numbers, na.rm = TRUE) + unidentified_copy_number,
    dna_ploidy_using_supplied_unidentified = dna_ploidy_supplied_unidentified,
    dna_ploidy_using_average_unidentified = dna_ploidy_avg_unidentified
  )
}

summarise_karyotype_ploidy <- function(karyo_tbl) {
  karyo_tbl %>%
    bind_cols(bind_rows(lapply(karyo_tbl$karyotype, compute_karyotype_ploidy_metrics))) %>%
    group_by(id) %>%
    summarise(
      n_karyotypes = dplyr::n(),
      numerical_ploidy_excluding_unidentified = mean(numerical_ploidy_excluding_unidentified, na.rm = TRUE),
      numerical_ploidy_including_unidentified = mean(numerical_ploidy_including_unidentified, na.rm = TRUE),
      mean_dna_ploidy_using_supplied_unidentified = mean(dna_ploidy_using_supplied_unidentified, na.rm = TRUE),
      sd_dna_ploidy_using_supplied_unidentified = stats::sd(dna_ploidy_using_supplied_unidentified, na.rm = TRUE),
      mean_dna_ploidy_using_average_unidentified = mean(dna_ploidy_using_average_unidentified, na.rm = TRUE),
      sd_dna_ploidy_using_average_unidentified = stats::sd(dna_ploidy_using_average_unidentified, na.rm = TRUE),
      .groups = "drop"
    )
}

build_karyotype_flow_lookup <- function(project_root) {
  karyo_path <- file.path(project_root, "data", "Hypoxia_SUM159_karyotyping.Rds")
  metadata_path <- file.path(project_root, "processed_data", "hypoxia-sum159", "sample_metadata.csv")

  karyo_tbl <- readRDS(karyo_path)
  metadata_tbl <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

  flow_tbl <- metadata_tbl %>%
    transmute(
      flow_sample_name = sample_id,
      condition,
      latest_match_date,
      relative_day,
      flow_key = sub("^Sample_SUM159_", "", sample_id)
    ) %>%
    bind_cols(parse_series_id(.$flow_key) %>% select(ploidy_label, lineage_label, passage))

  karyo_id_tbl <- tibble(karyo_id = sort(unique(karyo_tbl$id))) %>%
    bind_cols(parse_series_id(.$karyo_id) %>% select(ploidy_label, lineage_label, passage)) %>%
    mutate(
      lineage_label = case_when(
        karyo_id %in% c(
          "SUM-159_NLS_2N_A7M_K_harvest",
          "SUM-159_NLS_2N_MRCA_harvest"
        ) ~ "C",
        karyo_id %in% c(
          "SUM-159_NLS_4N_A5M_K_harvest",
          "SUM-159_NLS_4N_MRCA_harvest"
        ) ~ "C",
        TRUE ~ lineage_label
      )
    )

  ploidy_summary <- summarise_karyotype_ploidy(karyo_tbl)

  lookup_tbl <- bind_rows(lapply(seq_len(nrow(karyo_id_tbl)), function(i) {
    karyo_row <- karyo_id_tbl[i, , drop = FALSE]
    candidates <- flow_tbl %>%
      filter(
        ploidy_label == karyo_row$ploidy_label[[1]],
        lineage_label == karyo_row$lineage_label[[1]]
      ) %>%
      mutate(
        passage_diff = if (is.na(karyo_row$passage[[1]])) Inf else abs(passage - karyo_row$passage[[1]]),
        missing_passage = is.na(karyo_row$passage[[1]])
      ) %>%
      arrange(passage_diff, passage, latest_match_date, flow_sample_name)

    if (nrow(candidates) == 0) {
      return(tibble(
        karyo_id = karyo_row$karyo_id[[1]],
        karyo_ploidy_label = karyo_row$ploidy_label[[1]],
        karyo_lineage_label = karyo_row$lineage_label[[1]],
        karyo_passage = karyo_row$passage[[1]],
        flow_sample_name = NA_character_,
        flow_condition = NA_character_,
        flow_passage = NA_integer_,
        latest_match_date = NA_character_,
        relative_day = NA_real_,
        passage_diff = NA_real_,
        mapping_note = "no matching flow lineage"
      ))
    }

    best <- candidates[1, , drop = FALSE]
    mapping_note <- if (is.na(karyo_row$passage[[1]])) "control lineage fallback: earliest available passage" else "closest passage within matched lineage"

    tibble(
      karyo_id = karyo_row$karyo_id[[1]],
      karyo_ploidy_label = karyo_row$ploidy_label[[1]],
      karyo_lineage_label = karyo_row$lineage_label[[1]],
      karyo_passage = karyo_row$passage[[1]],
      flow_sample_name = best$flow_sample_name[[1]],
      flow_condition = best$condition[[1]],
      flow_passage = best$passage[[1]],
      latest_match_date = best$latest_match_date[[1]],
      relative_day = best$relative_day[[1]],
      passage_diff = if (is.finite(best$passage_diff[[1]])) best$passage_diff[[1]] else NA_real_,
      mapping_note = mapping_note
    )
  })) %>%
    left_join(ploidy_summary, by = c("karyo_id" = "id"))

  list(
    lookup_tbl = lookup_tbl,
    ploidy_summary = ploidy_summary
  )
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

trapz_integral <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) {
    return(NA_real_)
  }
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

map_dna_to_ploidy_counterfactual <- function(x, anchor_2n, anchor_4n, anchor_8n) {
  if (!all(is.finite(c(anchor_2n, anchor_4n, anchor_8n))) || !(anchor_2n < anchor_4n && anchor_4n < anchor_8n)) {
    return(rep(NA_real_, length(x)))
  }

  out <- numeric(length(x))
  idx_mid <- x <= anchor_4n
  out[idx_mid] <- 2 + 2 * (x[idx_mid] - anchor_2n) / (anchor_4n - anchor_2n)
  out[!idx_mid] <- 4 + 4 * (x[!idx_mid] - anchor_4n) / (anchor_8n - anchor_4n)
  out
}

select_evenly_spaced_draw_indices <- function(n_draw, max_draws) {
  max_draws <- as.integer(max_draws[[1]])
  if (!is.finite(max_draws) || max_draws < 1) {
    stop("`max_draws` must be a positive integer.")
  }
  if (max_draws >= n_draw) {
    return(seq_len(n_draw))
  }

  idx <- unique(as.integer(round(seq(1, n_draw, length.out = max_draws))))
  if (length(idx) < max_draws) {
    remaining <- setdiff(seq_len(n_draw), idx)
    idx <- sort(c(idx, remaining[seq_len(min(length(remaining), max_draws - length(idx)))]))
  }
  idx[seq_len(min(length(idx), max_draws))]
}

infer_histogram_bin_count <- function(values, min_bins = 60, max_bins = 160) {
  values <- values[is.finite(values)]
  if (length(values) < 2) {
    return(min_bins)
  }
  fd_bins <- suppressWarnings(grDevices::nclass.FD(values))
  if (!is.finite(fd_bins) || fd_bins < 1) {
    fd_bins <- ceiling(sqrt(length(values)))
  }
  as.integer(max(min_bins, min(max_bins, fd_bins)))
}

build_flow_event_histogram <- function(values, threshold, min_bins = 60, max_bins = 160) {
  x <- values[is.finite(values) & values > threshold]
  if (length(x) < 50 || length(unique(x)) < 2) {
    stop("Not enough above-threshold events to build a histogram.")
  }

  n_bins <- infer_histogram_bin_count(x, min_bins = min_bins, max_bins = max_bins)
  hist_obj <- graphics::hist(x, breaks = n_bins, plot = FALSE, include.lowest = TRUE, right = FALSE)
  binwidth <- diff(hist_obj$breaks)
  density <- hist_obj$counts / sum(hist_obj$counts) / binwidth

  tibble(
    bin_left = head(hist_obj$breaks, -1),
    bin_right = tail(hist_obj$breaks, -1),
    bin_mid = hist_obj$mids,
    bin_width = binwidth,
    count = hist_obj$counts,
    density = density
  )
}

build_flow_anchor_draws <- function(draws_matrix, prepared, spec) {
  input_tbl <- prepared$input_tbl
  n_sample <- nrow(input_tbl)
  n_draw <- nrow(draws_matrix)

  alpha <- draws_matrix[, "alpha"]
  beta <- if (isTRUE(spec$include_x_ratio_effect) && "beta" %in% colnames(draws_matrix)) draws_matrix[, "beta"] else rep(0, n_draw)
  log_M_cen <- draws_matrix[, "log_M_cen"]
  log_R_2n <- draws_matrix[, "log_R_2n"]
  log_R_4n_over_2n <- draws_matrix[, "log_R_4n_over_2n"]
  log_R_8n_over_4n <- if (isTRUE(spec$two_peak_adjacent_8n) && "log_R_8n_over_4n" %in% colnames(draws_matrix)) draws_matrix[, "log_R_8n_over_4n"] else rep(NA_real_, n_draw)
  log_rho <- if (!isTRUE(spec$rho_fixed_one) && "log_rho" %in% colnames(draws_matrix)) draws_matrix[, "log_rho"] else rep(0, n_draw)
  prob_4n <- extract_indexed_draw_matrix(draws_matrix, "prob_4n", n_sample)
  delta_dna <- if (isTRUE(spec$sample_delta_dna)) extract_indexed_draw_matrix(draws_matrix, "delta_dna", n_sample) else matrix(0, nrow = n_draw, ncol = n_sample)
  anchor_2n <- extract_indexed_draw_matrix(draws_matrix, "mu_peak_2n", n_sample)
  anchor_4n <- extract_indexed_draw_matrix(draws_matrix, "mu_peak_4n", n_sample)
  anchor_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) extract_indexed_draw_matrix(draws_matrix, "mu_peak_8n", n_sample) else matrix(NA_real_, nrow = n_draw, ncol = n_sample)

  if ("log_R_4n_over_2n" %in% colnames(draws_matrix)) {
    ratio_4n_over_2n <- exp(draws_matrix[, "log_R_4n_over_2n"])
  } else if (all(c("R_4n", "R_2n") %in% colnames(draws_matrix))) {
    ratio_4n_over_2n <- draws_matrix[, "R_4n"] / draws_matrix[, "R_2n"]
  } else {
    stop("Could not reconstruct the posterior 4N/2N peak ratio draws.")
  }

  if (isTRUE(spec$two_peak_adjacent_8n)) {
    if ("log_R_8n_over_4n" %in% colnames(draws_matrix)) {
      ratio_8n_over_4n <- exp(draws_matrix[, "log_R_8n_over_4n"])
    } else if (all(c("R_8n", "R_4n") %in% colnames(draws_matrix))) {
      ratio_8n_over_4n <- draws_matrix[, "R_8n"] / draws_matrix[, "R_4n"]
    } else {
      stop("Could not reconstruct the posterior 8N/4N peak ratio draws.")
    }
  } else {
    ratio_8n_over_4n <- rep(NA_real_, n_draw)
  }

  M_cen <- exp(log_M_cen)
  R_2n <- exp(log_R_2n)
  R_4n <- R_2n * exp(log_R_4n_over_2n)
  R_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) R_4n * exp(log_R_8n_over_4n) else rep(NA_real_, n_draw)
  rho <- if (isTRUE(spec$rho_fixed_one)) rep(1, n_draw) else exp(log_rho)
  reference_eta <- alpha + beta * 0
  reference_u <- exp(-reference_eta)
  reference_mu_cen <- M_cen * reference_u / (1 + reference_u)
  reference_distortion <- (1 + reference_u) / (1 + rho * reference_u)
  delta_dna_multiplier <- exp(delta_dna)
  cf_anchor_2n <- delta_dna_multiplier * (reference_mu_cen * R_2n * reference_distortion)
  cf_anchor_4n <- delta_dna_multiplier * (reference_mu_cen * R_4n * reference_distortion)
  cf_anchor_8n <- if (isTRUE(spec$two_peak_adjacent_8n)) {
    delta_dna_multiplier * (reference_mu_cen * R_8n * reference_distortion)
  } else {
    matrix(NA_real_, nrow = n_draw, ncol = n_sample)
  }

  assigned_state <- ifelse(prob_4n >= 0.5, "4N", "2N")
  lower_anchor <- ifelse(assigned_state == "4N", anchor_4n, anchor_2n)
  upper_anchor <- ifelse(assigned_state == "4N", anchor_8n, anchor_4n)

  tibble(
    draw_id = rep(seq_len(n_draw), times = n_sample),
    sample_index = rep(seq_len(n_sample), each = n_draw),
    sample_name = rep(input_tbl$sample_name, each = n_draw),
    condition = rep(input_tbl$condition, each = n_draw),
    latest_match_date = rep(as.character(input_tbl$latest_match_date), each = n_draw),
    relative_day = rep(input_tbl$relative_day, each = n_draw),
    prob_4n = as.vector(prob_4n),
    assigned_state = as.vector(assigned_state),
    delta_dna = as.vector(delta_dna),
    anchor_2n = as.vector(anchor_2n),
    anchor_4n = as.vector(anchor_4n),
    anchor_8n = as.vector(anchor_8n),
    cf_anchor_2n = as.vector(cf_anchor_2n),
    cf_anchor_4n = as.vector(cf_anchor_4n),
    cf_anchor_8n = as.vector(cf_anchor_8n),
    lower_anchor = as.vector(lower_anchor),
    upper_anchor = as.vector(upper_anchor),
    ratio_4n_over_2n = rep(ratio_4n_over_2n, times = n_sample),
    ratio_8n_over_4n = rep(ratio_8n_over_4n, times = n_sample)
  )
}

build_flow_component_spec <- function(anchor_row) {
  if (identical(anchor_row$assigned_state[[1]], "2N")) {
    anchor_2n <- anchor_row$anchor_2n[[1]]
    anchor_4n <- anchor_row$anchor_4n[[1]]
    anchor_8n <- anchor_row$anchor_8n[[1]]
    midpoint_24 <- (anchor_2n + anchor_4n) / 2
    midpoint_48 <- (anchor_4n + anchor_8n) / 2
    spacing_24 <- max(anchor_4n - anchor_2n, 1)
    spacing_48 <- max(anchor_8n - anchor_4n, 1)
    means <- c(anchor_2n, midpoint_24, anchor_4n, anchor_4n, midpoint_48, anchor_8n)
    labels <- c("g1_2n", "s_2n_4n", "g2m_2n", "g1_4n", "s_4n_8n", "g2m_4n")
    lineage_group <- c("2N", "2N", "2N", "4N", "4N", "4N")
    phase_group <- c("g1", "s", "g2m", "g1", "s", "g2m")
    sigma_init <- c(spacing_24, spacing_24, spacing_24, spacing_48, spacing_48, spacing_48) * c(0.08, 0.16, 0.08, 0.08, 0.16, 0.08)
    sigma_max <- c(spacing_24, spacing_24, spacing_24, spacing_48, spacing_48, spacing_48) * c(0.30, 0.45, 0.30, 0.30, 0.45, 0.30)
  } else {
    anchor_4n <- anchor_row$anchor_4n[[1]]
    anchor_8n <- anchor_row$anchor_8n[[1]]
    midpoint_48 <- (anchor_4n + anchor_8n) / 2
    spacing_48 <- max(anchor_8n - anchor_4n, 1)
    means <- c(anchor_4n, midpoint_48, anchor_8n)
    labels <- c("g1_4n", "s_4n_8n", "g2m_4n")
    lineage_group <- c("4N", "4N", "4N")
    phase_group <- c("g1", "s", "g2m")
    sigma_init <- c(spacing_48, spacing_48, spacing_48) * c(0.08, 0.16, 0.08)
    sigma_max <- c(spacing_48, spacing_48, spacing_48) * c(0.30, 0.45, 0.30)
  }

  tibble(
    component_index = seq_along(means),
    component_label = labels,
    family = "gaussian",
    lineage_group = lineage_group,
    phase_group = phase_group,
    mean = means,
    sigma_init = pmax(sigma_init, 400),
    sigma_min = pmax(0.35 * sigma_init, 200),
    sigma_max = pmax(sigma_max, pmax(2.5 * sigma_init, 1200))
  ) %>%
    bind_rows(
      tibble(
        component_index = length(means) + 1L,
        component_label = "noise",
        family = "uniform",
        lineage_group = "noise",
        phase_group = "noise",
        mean = mean(range(means)),
        sigma_init = NA_real_,
        sigma_min = NA_real_,
        sigma_max = NA_real_
      )
    )
}

fit_fixed_mean_histogram_mixture <- function(hist_tbl, component_tbl, max_iter = 25, sigma_prior_strength = 10, dirichlet_alpha = 1.05) {
  x <- hist_tbl$bin_mid
  counts <- hist_tbl$count
  means <- component_tbl$mean
  sigma_init <- component_tbl$sigma_init
  sigma_min <- component_tbl$sigma_min
  sigma_max <- component_tbl$sigma_max
  family <- component_tbl$family
  lineage_group <- component_tbl$lineage_group
  phase_group <- component_tbl$phase_group
  k <- length(means)
  gaussian_idx <- which(family == "gaussian")
  uniform_idx <- which(family == "uniform")
  biological_idx <- which(family != "uniform")
  lineage_levels <- unique(lineage_group[biological_idx])
  phase_levels <- c("g1", "s", "g2m")
  x_range <- range(c(hist_tbl$bin_left, hist_tbl$bin_right))
  uniform_density <- 1 / diff(x_range)

  weights <- rep(1 / k, k)
  sigmas <- sigma_init
  loglik <- NA_real_

  for (iter in seq_len(max_iter)) {
    dens_mat <- vapply(seq_len(k), function(j) {
      if (identical(family[[j]], "uniform")) {
        rep(uniform_density, length(x))
      } else {
        stats::dnorm(x, mean = means[[j]], sd = sigmas[[j]])
      }
    }, numeric(length(x)))
    dens_mat <- pmax(dens_mat, 1e-300)
    weighted_mat <- sweep(dens_mat, 2, weights, `*`)
    row_total <- pmax(rowSums(weighted_mat), 1e-300)
    resp <- weighted_mat / row_total
    weighted_resp <- resp * counts
    nk <- colSums(weighted_resp)

    noise_weight <- if (length(uniform_idx)) {
      pmax(nk[uniform_idx] + dirichlet_alpha - 1, 1e-8)
    } else {
      0
    }
    noise_weight <- noise_weight / sum(c(noise_weight, pmax(sum(nk[biological_idx]), 1e-8)))
    biological_total <- 1 - noise_weight

    lineage_counts <- vapply(lineage_levels, function(lineage_name) {
      sum(nk[lineage_group == lineage_name], na.rm = TRUE)
    }, numeric(1))
    lineage_weights <- pmax(lineage_counts + dirichlet_alpha - 1, 1e-8)
    lineage_weights <- lineage_weights / sum(lineage_weights)
    names(lineage_weights) <- lineage_levels

    phase_counts <- vapply(phase_levels, function(phase_name) {
      sum(nk[phase_group == phase_name], na.rm = TRUE)
    }, numeric(1))
    phase_weights <- pmax(phase_counts + dirichlet_alpha - 1, 1e-8)
    phase_weights <- phase_weights / sum(phase_weights)
    names(phase_weights) <- phase_levels

    weights[] <- 0
    if (length(biological_idx)) {
      weights[biological_idx] <- biological_total *
        lineage_weights[lineage_group[biological_idx]] *
        phase_weights[phase_group[biological_idx]]
    }
    if (length(uniform_idx)) {
      weights[uniform_idx] <- noise_weight
    }

    if (length(gaussian_idx)) {
      for (phase_name in phase_levels) {
        phase_idx <- gaussian_idx[phase_group[gaussian_idx] == phase_name]
        if (!length(phase_idx)) {
          next
        }
        phase_resp_total <- 0
        phase_var_total <- 0
        prior_var <- mean(sigma_init[phase_idx]^2, na.rm = TRUE)
        for (j in phase_idx) {
          centered_sq <- (x - means[[j]])^2
          phase_resp_total <- phase_resp_total + sum(weighted_resp[, j] * centered_sq)
          phase_var_total <- phase_var_total + nk[[j]]
        }
        sigma_sq <- (phase_resp_total + sigma_prior_strength * prior_var) / pmax(phase_var_total + sigma_prior_strength, 1e-8)
        shared_sigma <- sqrt(sigma_sq)
        shared_sigma <- min(max(shared_sigma, min(sigma_min[phase_idx], na.rm = TRUE)), max(sigma_max[phase_idx], na.rm = TRUE))
        sigmas[phase_idx] <- shared_sigma
      }
    }

    loglik <- sum(counts * log(row_total))
  }

  bin_width <- stats::median(hist_tbl$bin_width)
  fitted_density <- rowSums(vapply(seq_len(k), function(j) {
    if (identical(family[[j]], "uniform")) {
      weights[[j]] * rep(uniform_density, length(x))
    } else {
      weights[[j]] * stats::dnorm(x, mean = means[[j]], sd = sigmas[[j]])
    }
  }, numeric(length(x))))

  list(
    component_tbl = component_tbl %>%
      mutate(weight = weights, sigma = sigmas),
    overlay_tbl = hist_tbl %>%
      transmute(
        bin_mid,
        observed_count = count,
        observed_density = density,
        fitted_density = fitted_density,
        fitted_count = fitted_density * sum(counts) * bin_width
      ),
    loglik = loglik
  )
}

fit_flow_posterior_histograms <- function(
  project_root,
  spec_name = "adjacent_8n_delta_dna",
  output_category = "extensions",
  max_draws = NULL,
  threshold = NULL,
  min_bins = 60,
  max_bins = 160
) {
  spec_registry <- c(get_ablation_registry(), get_extension_registry())
  if (!spec_name %in% names(spec_registry)) {
    stop(sprintf("Unknown spec `%s`.", spec_name))
  }

  spec <- spec_registry[[spec_name]]
  if (!isTRUE(spec$two_peak_adjacent_8n) || !isTRUE(spec$sample_delta_dna)) {
    stop("This posterior histogram fitter currently expects the adjacent_8n_delta_dna-style extension.")
  }

  prepared_path <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", "prepared_input.rds")
  prepared <- readRDS(prepared_path)
  threshold <- if (is.null(threshold)) prepared$threshold else threshold

  nuts_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", output_category, spec$name, "nuts")
  draws_path <- file.path(nuts_dir, "draws_matrix.rds")
  draws_matrix <- readRDS(draws_path)
  if (!is.null(max_draws)) {
    draw_idx <- select_evenly_spaced_draw_indices(nrow(draws_matrix), max_draws)
    draws_matrix <- draws_matrix[draw_idx, , drop = FALSE]
  }
  anchor_draw_tbl <- build_flow_anchor_draws(draws_matrix, prepared, spec)

  dna_path <- file.path(project_root, "processed_data", "hypoxia-sum159", "filtered_dna_area_vectors.rds")
  raw_vectors <- readRDS(dna_path)$raw

  sample_histograms <- lapply(prepared$input_tbl$sample_name, function(sample_name) {
    if (!sample_name %in% names(raw_vectors)) {
      stop(sprintf("Missing raw DNA vector for sample `%s`.", sample_name))
    }
    build_flow_event_histogram(raw_vectors[[sample_name]], threshold = threshold, min_bins = min_bins, max_bins = max_bins)
  })
  names(sample_histograms) <- prepared$input_tbl$sample_name

  sample_objects <- vector("list", length(sample_histograms))
  names(sample_objects) <- names(sample_histograms)
  component_draw_parts <- vector("list", length(sample_histograms))
  draw_summary_parts <- vector("list", length(sample_histograms))
  overlay_parts <- vector("list", length(sample_histograms))

  for (sample_i in seq_along(sample_histograms)) {
    sample_name <- names(sample_histograms)[[sample_i]]
    hist_tbl <- sample_histograms[[sample_i]]
    sample_anchor_tbl <- anchor_draw_tbl %>%
      filter(sample_name == !!sample_name)

    fit_parts <- lapply(seq_len(nrow(sample_anchor_tbl)), function(draw_i) {
      anchor_row <- sample_anchor_tbl[draw_i, , drop = FALSE]
      component_spec <- build_flow_component_spec(anchor_row)
      fit_fixed_mean_histogram_mixture(hist_tbl, component_spec)
    })

    component_draw_tbl <- bind_rows(lapply(seq_len(nrow(sample_anchor_tbl)), function(draw_i) {
      anchor_row <- sample_anchor_tbl[draw_i, , drop = FALSE]
      fit_obj <- fit_parts[[draw_i]]
      fit_obj$component_tbl %>%
        transmute(
          draw_id = anchor_row$draw_id[[1]],
          sample_index = anchor_row$sample_index[[1]],
          sample_name = anchor_row$sample_name[[1]],
          assigned_state = anchor_row$assigned_state[[1]],
          component_index,
          component_label,
          constrained_location = mean,
          weight,
          sigma
        )
    }))

    draw_summary_tbl <- bind_rows(lapply(seq_len(nrow(sample_anchor_tbl)), function(draw_i) {
      anchor_row <- sample_anchor_tbl[draw_i, , drop = FALSE]
      fit_obj <- fit_parts[[draw_i]]
      tibble(
        draw_id = anchor_row$draw_id[[1]],
        sample_index = anchor_row$sample_index[[1]],
        sample_name = anchor_row$sample_name[[1]],
        condition = anchor_row$condition[[1]],
        latest_match_date = anchor_row$latest_match_date[[1]],
        relative_day = anchor_row$relative_day[[1]],
        assigned_state = anchor_row$assigned_state[[1]],
        prob_4n = anchor_row$prob_4n[[1]],
        delta_dna = anchor_row$delta_dna[[1]],
        anchor_2n = anchor_row$anchor_2n[[1]],
        anchor_4n = anchor_row$anchor_4n[[1]],
        anchor_8n = anchor_row$anchor_8n[[1]],
        cf_anchor_2n = anchor_row$cf_anchor_2n[[1]],
        cf_anchor_4n = anchor_row$cf_anchor_4n[[1]],
        cf_anchor_8n = anchor_row$cf_anchor_8n[[1]],
        lower_anchor = anchor_row$lower_anchor[[1]],
        upper_anchor = anchor_row$upper_anchor[[1]],
        loglik = fit_obj$loglik
      )
    }))

    overlay_array <- simplify2array(lapply(fit_parts, function(fit_obj) fit_obj$overlay_tbl$fitted_count))
    overlay_density_array <- simplify2array(lapply(fit_parts, function(fit_obj) fit_obj$overlay_tbl$fitted_density))
    overlay_quant <- t(apply(overlay_array, 1, safe_quantile, probs = c(0.05, 0.5, 0.95)))
    overlay_density_quant <- t(apply(overlay_density_array, 1, safe_quantile, probs = c(0.05, 0.5, 0.95)))
    overlay_tbl <- hist_tbl %>%
      transmute(
        sample_name = sample_name,
        bin_left,
        bin_right,
        bin_mid,
        observed_count = count,
        observed_density = density,
        fitted_count_q05 = overlay_quant[, 1],
        fitted_count_median = overlay_quant[, 2],
        fitted_count_q95 = overlay_quant[, 3],
        fitted_density_q05 = overlay_density_quant[, 1],
        fitted_density_median = overlay_density_quant[, 2],
        fitted_density_q95 = overlay_density_quant[, 3]
      )

    component_summary_tbl <- component_draw_tbl %>%
      group_by(sample_name, assigned_state, component_index, component_label) %>%
      summarise(
        constrained_location_mean = mean(constrained_location, na.rm = TRUE),
        constrained_location_median = stats::median(constrained_location, na.rm = TRUE),
        weight_mean = mean(weight, na.rm = TRUE),
        weight_median = stats::median(weight, na.rm = TRUE),
        weight_q05 = safe_quantile(weight, 0.05),
        weight_q95 = safe_quantile(weight, 0.95),
        sigma_mean = mean(sigma, na.rm = TRUE),
        sigma_median = stats::median(sigma, na.rm = TRUE),
        sigma_q05 = safe_quantile(sigma, 0.05),
        sigma_q95 = safe_quantile(sigma, 0.95),
        .groups = "drop"
      )

    sample_summary_tbl <- draw_summary_tbl %>%
      summarise(
        sample_name = first(sample_name),
        condition = first(condition),
        latest_match_date = first(latest_match_date),
        relative_day = first(relative_day),
        posterior_prob_4n = mean(prob_4n, na.rm = TRUE),
        assigned_state_map = ifelse(mean(prob_4n, na.rm = TRUE) >= 0.5, "4N", "2N"),
        lower_anchor_mean = mean(lower_anchor, na.rm = TRUE),
        lower_anchor_median = stats::median(lower_anchor, na.rm = TRUE),
        lower_anchor_q05 = safe_quantile(lower_anchor, 0.05),
        lower_anchor_q95 = safe_quantile(lower_anchor, 0.95),
        upper_anchor_mean = mean(upper_anchor, na.rm = TRUE),
        upper_anchor_median = stats::median(upper_anchor, na.rm = TRUE),
        upper_anchor_q05 = safe_quantile(upper_anchor, 0.05),
        upper_anchor_q95 = safe_quantile(upper_anchor, 0.95),
        mean_loglik = mean(loglik, na.rm = TRUE),
        gated_n_above_threshold = sum(hist_tbl$observed_count)
      )

    sample_objects[[sample_i]] <- list(
      hist_tbl = hist_tbl,
      draw_summary_tbl = draw_summary_tbl,
      component_draw_tbl = component_draw_tbl,
      component_summary_tbl = component_summary_tbl,
      overlay_tbl = overlay_tbl,
      sample_summary_tbl = sample_summary_tbl
    )
    component_draw_parts[[sample_i]] <- component_draw_tbl
    draw_summary_parts[[sample_i]] <- draw_summary_tbl
    overlay_parts[[sample_i]] <- overlay_tbl
  }

  component_draw_tbl <- bind_rows(component_draw_parts)
  draw_summary_tbl <- bind_rows(draw_summary_parts)
  overlay_tbl <- bind_rows(overlay_parts)
  component_summary_tbl <- bind_rows(lapply(sample_objects, `[[`, "component_summary_tbl"))
  sample_summary_tbl <- bind_rows(lapply(sample_objects, `[[`, "sample_summary_tbl")) %>%
    arrange(lower_anchor_median)

  distribution_object <- list(
    metadata = list(
      project_root = normalizePath(project_root, winslash = "/", mustWork = TRUE),
      spec_name = spec$name,
      output_category = output_category,
      threshold = threshold,
      n_draw = nrow(draws_matrix),
      n_sample = nrow(prepared$input_tbl)
    ),
    samples = sample_objects,
    sample_summary_tbl = sample_summary_tbl,
    component_summary_tbl = component_summary_tbl,
    draw_summary_tbl = draw_summary_tbl,
    component_draw_tbl = component_draw_tbl,
    overlay_tbl = overlay_tbl
  )

  list(
    spec = spec,
    nuts_dir = nuts_dir,
    threshold = threshold,
    distribution_object = distribution_object
  )
}

plot_flow_histogram_overlays <- function(overlay_tbl, sample_summary_tbl) {
  sample_levels <- unique(sample_summary_tbl$sample_name)
  plot_tbl <- overlay_tbl %>%
    inner_join(sample_summary_tbl %>% select(sample_name, lower_anchor_median, assigned_state_map), by = "sample_name") %>%
    mutate(
      sample_label = factor(sample_name, levels = sample_levels)
    )

  ggplot(plot_tbl, aes(x = bin_mid)) +
    geom_rect(aes(xmin = bin_left, xmax = bin_right, ymin = 0, ymax = observed_density), fill = "grey82", color = "white", linewidth = 0.12) +
    geom_ribbon(aes(ymin = fitted_density_q05, ymax = fitted_density_q95), fill = "#9ecae1", alpha = 0.45) +
    geom_line(aes(y = fitted_density_median), color = "#08519c", linewidth = 0.5) +
    geom_vline(aes(xintercept = lower_anchor_median), color = "#cb181d", linewidth = 0.35, linetype = "dashed") +
    facet_wrap(~ sample_label, scales = "free_y", ncol = 4) +
    labs(
      title = "Flow-only posterior histogram overlays",
      subtitle = "Observed above-threshold DNA-A histograms with draw-propagated constrained fits.",
      x = "DNA-A",
      y = "Density"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

build_flow_state_distribution_outputs <- function(distribution_object, g0g1_samples_per_draw = 200, predictive_seed = 1L) {
  component_draw_tbl <- distribution_object$component_draw_tbl
  draw_anchor_tbl <- distribution_object$draw_summary_tbl %>%
    distinct(sample_name, draw_id, anchor_2n, anchor_4n, anchor_8n, cf_anchor_2n, cf_anchor_4n, cf_anchor_8n)
  sample_summary_tbl <- distribution_object$sample_summary_tbl %>%
    distinct(sample_name, .keep_all = TRUE) %>%
    mutate(sample_name = as.character(sample_name))
  sample_levels <- unique(sample_summary_tbl$sample_name)

  biological_draw_tbl <- component_draw_tbl %>%
    filter(component_label != "noise") %>%
    group_by(sample_name, draw_id) %>%
    mutate(
      biological_weight_total = sum(weight, na.rm = TRUE),
      state_fraction_biological = weight / pmax(biological_weight_total, 1e-8)
    ) %>%
    ungroup()

  state_summary_tbl <- biological_draw_tbl %>%
    group_by(sample_name, component_label) %>%
    summarise(
      mean_state_weight = mean(weight, na.rm = TRUE),
      median_state_weight = stats::median(weight, na.rm = TRUE),
      q05_state_weight = safe_quantile(weight, 0.05),
      q95_state_weight = safe_quantile(weight, 0.95),
      mean_state_fraction_biological = mean(state_fraction_biological, na.rm = TRUE),
      median_state_fraction_biological = stats::median(state_fraction_biological, na.rm = TRUE),
      q05_state_fraction_biological = safe_quantile(state_fraction_biological, 0.05),
      q95_state_fraction_biological = safe_quantile(state_fraction_biological, 0.95),
      .groups = "drop"
    ) %>%
    mutate(sample_name = factor(sample_name, levels = sample_levels))

  g0g1_draw_tbl <- component_draw_tbl %>%
    filter(component_label %in% c("g1_2n", "g1_4n")) %>%
    group_by(sample_name, draw_id) %>%
    mutate(
      g0g1_weight_total = sum(weight, na.rm = TRUE),
      g0g1_fraction = weight / pmax(g0g1_weight_total, 1e-8)
    ) %>%
    ungroup() %>%
    select(sample_name, draw_id, component_label, constrained_location, sigma, weight, g0g1_fraction) %>%
    tidyr::pivot_wider(
      names_from = component_label,
      values_from = c(constrained_location, sigma, weight, g0g1_fraction),
      values_fill = 0
    ) %>%
    mutate(
      sample_name = as.character(sample_name),
      g1_2n_location = ifelse(weight_g1_2n > 0, constrained_location_g1_2n, NA_real_),
      g1_4n_location = ifelse(weight_g1_4n > 0, constrained_location_g1_4n, NA_real_),
      g1_2n_sigma = ifelse(weight_g1_2n > 0, sigma_g1_2n, NA_real_),
      g1_4n_sigma = ifelse(weight_g1_4n > 0, sigma_g1_4n, NA_real_),
      g1_2n_weight = weight_g1_2n,
      g1_4n_weight = weight_g1_4n,
      g1_2n_fraction = g0g1_fraction_g1_2n,
      g1_4n_fraction = g0g1_fraction_g1_4n
    ) %>%
    select(sample_name, draw_id, g1_2n_location, g1_4n_location, g1_2n_sigma, g1_4n_sigma, g1_2n_weight, g1_4n_weight, g1_2n_fraction, g1_4n_fraction)

  g0g1_draw_tbl <- g0g1_draw_tbl %>%
    left_join(draw_anchor_tbl, by = c("sample_name", "draw_id"))

  g0g1_summary_tbl <- g0g1_draw_tbl %>%
    group_by(sample_name) %>%
    summarise(
      g1_2n_fraction_mean = mean(g1_2n_fraction, na.rm = TRUE),
      g1_2n_fraction_median = stats::median(g1_2n_fraction, na.rm = TRUE),
      g1_2n_fraction_q05 = safe_quantile(g1_2n_fraction, 0.05),
      g1_2n_fraction_q95 = safe_quantile(g1_2n_fraction, 0.95),
      g1_4n_fraction_mean = mean(g1_4n_fraction, na.rm = TRUE),
      g1_4n_fraction_median = stats::median(g1_4n_fraction, na.rm = TRUE),
      g1_4n_fraction_q05 = safe_quantile(g1_4n_fraction, 0.05),
      g1_4n_fraction_q95 = safe_quantile(g1_4n_fraction, 0.95),
      g1_2n_location_mean = mean(g1_2n_location, na.rm = TRUE),
      g1_2n_location_median = stats::median(g1_2n_location, na.rm = TRUE),
      g1_2n_location_q05 = safe_quantile(g1_2n_location, 0.05),
      g1_2n_location_q95 = safe_quantile(g1_2n_location, 0.95),
      g1_4n_location_mean = mean(g1_4n_location, na.rm = TRUE),
      g1_4n_location_median = stats::median(g1_4n_location, na.rm = TRUE),
      g1_4n_location_q05 = safe_quantile(g1_4n_location, 0.05),
      g1_4n_location_q95 = safe_quantile(g1_4n_location, 0.95),
      .groups = "drop"
    ) %>%
    left_join(sample_summary_tbl %>% select(sample_name, condition, latest_match_date, relative_day, assigned_state_map), by = "sample_name") %>%
    mutate(sample_name = factor(sample_name, levels = sample_levels))

  set.seed(as.integer(predictive_seed))
  g0g1_predictive_samples_tbl <- bind_rows(lapply(split(g0g1_draw_tbl, g0g1_draw_tbl$sample_name), function(sample_draws) {
    sample_name <- unique(sample_draws$sample_name)[[1]]
    bind_rows(lapply(seq_len(nrow(sample_draws)), function(i) {
      row_i <- sample_draws[i, , drop = FALSE]
      a2 <- row_i$cf_anchor_2n[[1]]
      a4 <- row_i$cf_anchor_4n[[1]]
      a8 <- row_i$cf_anchor_8n[[1]]
      if (!all(is.finite(c(a2, a4, a8))) || !(a2 < a4 && a4 < a8)) {
        return(NULL)
      }

      component_probs <- c(row_i$g1_2n_fraction[[1]], row_i$g1_4n_fraction[[1]])
      component_probs[!is.finite(component_probs)] <- 0
      if (sum(component_probs) <= 0) {
        return(NULL)
      }
      component_probs <- component_probs / sum(component_probs)
      sampled_component <- sample(c("g1_2n", "g1_4n"), size = g0g1_samples_per_draw, replace = TRUE, prob = component_probs)
      dna_values <- numeric(g0g1_samples_per_draw)

      idx_2n <- sampled_component == "g1_2n"
      if (any(idx_2n)) {
        dna_values[idx_2n] <- stats::rnorm(sum(idx_2n), mean = row_i$g1_2n_location[[1]], sd = row_i$g1_2n_sigma[[1]])
      }
      idx_4n <- sampled_component == "g1_4n"
      if (any(idx_4n)) {
        dna_values[idx_4n] <- stats::rnorm(sum(idx_4n), mean = row_i$g1_4n_location[[1]], sd = row_i$g1_4n_sigma[[1]])
      }

      tibble(
        sample_name = sample_name,
        draw_id = row_i$draw_id[[1]],
        sampled_component = sampled_component,
        dna_value = dna_values,
        ploidy_value = map_dna_to_ploidy_counterfactual(dna_values, a2, a4, a8)
      )
    }))
  }))

  g0g1_predictive_density_tbl <- bind_rows(lapply(split(g0g1_predictive_samples_tbl, g0g1_predictive_samples_tbl$sample_name), function(sample_tbl) {
    ploidy_vals <- sample_tbl$ploidy_value[is.finite(sample_tbl$ploidy_value)]
    if (length(ploidy_vals) < 10 || length(unique(ploidy_vals)) < 2) {
      return(tibble(sample_name = unique(sample_tbl$sample_name)[[1]], ploidy = numeric(), predictive_density = numeric()))
    }
    dens <- stats::density(ploidy_vals, from = 0.5, to = 8.5, n = 400, na.rm = TRUE)
    tibble(
      sample_name = unique(sample_tbl$sample_name)[[1]],
      ploidy = dens$x,
      predictive_density = dens$y
    )
  }))

  g0g1_density_grid_tbl <- bind_rows(lapply(split(g0g1_predictive_samples_tbl, g0g1_predictive_samples_tbl$sample_name), function(sample_tbl) {
    sample_name <- unique(sample_tbl$sample_name)[[1]]
    ploidy_vals <- sample_tbl$ploidy_value[is.finite(sample_tbl$ploidy_value)]
    if (length(ploidy_vals) < 10 || length(unique(ploidy_vals)) < 2) {
      return(tibble(
        sample_name = sample_name,
        grid_index = integer(),
        ploidy = numeric(),
        log_density = numeric()
      ))
    }

    dens <- stats::density(ploidy_vals, from = 0.5, to = 8.5, n = 200, na.rm = TRUE)
    density_floored <- pmax(dens$y, 1e-12)
    density_norm <- density_floored / trapz_integral(dens$x, density_floored)

    tibble(
      sample_name = sample_name,
      grid_index = seq_along(dens$x),
      ploidy = dens$x,
      log_density = log(density_norm)
    )
  }))

  list(
    state_summary_tbl = state_summary_tbl,
    g0g1_draw_tbl = g0g1_draw_tbl,
    g0g1_summary_tbl = g0g1_summary_tbl,
    g0g1_predictive_samples_tbl = g0g1_predictive_samples_tbl,
    g0g1_predictive_density_tbl = g0g1_predictive_density_tbl,
    g0g1_density_grid_tbl = g0g1_density_grid_tbl
  )
}

plot_flow_state_distribution <- function(state_summary_tbl) {
  plot_tbl <- state_summary_tbl %>%
    mutate(
      sample_name = factor(as.character(sample_name), levels = unique(as.character(sample_name))),
      component_label = factor(
        component_label,
        levels = c("g1_2n", "s_2n_4n", "g2m_2n", "g1_4n", "s_4n_8n", "g2m_4n")
      )
    )

  ggplot(plot_tbl, aes(x = sample_name, y = mean_state_fraction_biological, fill = component_label)) +
    geom_col(width = 0.82) +
    scale_fill_manual(
      values = c(
        "g1_2n" = "#2166ac",
        "s_2n_4n" = "#67a9cf",
        "g2m_2n" = "#d1e5f0",
        "g1_4n" = "#b2182b",
        "s_4n_8n" = "#ef8a62",
        "g2m_4n" = "#fddbc7"
      )
    ) +
    labs(
      title = "Posterior mean biological state composition",
      subtitle = "Stacked fractions are normalized across the six biological states only; noise is excluded.",
      x = NULL,
      y = "Fraction of biological mass",
      fill = "State"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
    )
}

plot_g0g1_ploidy_stacked_bar <- function(g0g1_summary_tbl) {
  plot_tbl <- g0g1_summary_tbl %>%
    mutate(sample_name = factor(as.character(sample_name), levels = unique(as.character(sample_name)))) %>%
    select(sample_name, g1_2n_fraction_mean, g1_4n_fraction_mean) %>%
    tidyr::pivot_longer(
      cols = c(g1_2n_fraction_mean, g1_4n_fraction_mean),
      names_to = "g0g1_state",
      values_to = "fraction"
    ) %>%
    mutate(
      g0g1_state = recode(
        g0g1_state,
        g1_2n_fraction_mean = "G0/G1 2N",
        g1_4n_fraction_mean = "G0/G1 4N"
      )
    )

  ggplot(plot_tbl, aes(x = sample_name, y = fraction, fill = g0g1_state)) +
    geom_col(width = 0.82) +
    scale_fill_manual(values = c("G0/G1 2N" = "#2166ac", "G0/G1 4N" = "#b2182b")) +
    labs(
      title = "Posterior mean G0/G1 ploidy composition",
      subtitle = "Fractions are normalized within the G0/G1 subpopulation.",
      x = NULL,
      y = "Fraction of G0/G1 mass",
      fill = "G0/G1 state"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)
    )
}

plot_g0g1_ploidy_density <- function(g0g1_predictive_density_tbl, g0g1_summary_tbl) {
  sample_levels <- unique(as.character(g0g1_summary_tbl$sample_name))
  plot_tbl <- g0g1_predictive_density_tbl %>%
    mutate(sample_name = factor(sample_name, levels = sample_levels))

  ggplot(plot_tbl, aes(x = ploidy)) +
    geom_line(aes(y = predictive_density), color = "#08519c", linewidth = 0.6) +
    facet_wrap(~ sample_name, scales = "free_y", ncol = 4) +
    scale_x_continuous(
      breaks = c(2, 4, 8),
      labels = c("2N", "4N", "8N")
    ) +
    labs(
      title = "Posterior predictive G0/G1 ploidy density",
      subtitle = "Each sample pools predictive ploidy samples drawn from the fitted G0/G1 mixtures under the draw-specific delta-only counterfactual DNA-to-ploidy map.",
      x = "Ploidy",
      y = "Density"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

write_flow_posterior_histogram_outputs <- function(
  project_root,
  spec_name = "adjacent_8n_delta_dna",
  output_category = "extensions",
  max_draws = NULL,
  threshold = NULL,
  min_bins = 60,
  max_bins = 160
) {
  fit_parts <- fit_flow_posterior_histograms(
    project_root = project_root,
    spec_name = spec_name,
    output_category = output_category,
    max_draws = max_draws,
    threshold = threshold,
    min_bins = min_bins,
    max_bins = max_bins
  )

  out_dir <- file.path(fit_parts$nuts_dir, "flow_only_distribution")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  distribution_rds_path <- file.path(out_dir, "flow_only_distribution.rds")
  sample_summary_csv_path <- file.path(out_dir, "sample_flow_distribution_summary.csv")
  component_summary_csv_path <- file.path(out_dir, "component_flow_distribution_summary.csv")
  draw_summary_csv_path <- file.path(out_dir, "draw_flow_distribution_summary.csv")
  component_draw_csv_path <- file.path(out_dir, "component_flow_distribution_draws.csv")
  overlay_csv_path <- file.path(out_dir, "flow_distribution_overlay_summary.csv")
  overlay_plot_path <- file.path(out_dir, "flow_distribution_overlays.png")
  state_summary_csv_path <- file.path(out_dir, "state_distribution_summary.csv")
  state_plot_path <- file.path(out_dir, "state_distribution_stacked_bar.png")
  g0g1_draw_csv_path <- file.path(out_dir, "g0g1_ploidy_distribution_draws.csv")
  g0g1_summary_csv_path <- file.path(out_dir, "g0g1_ploidy_distribution_summary.csv")
  g0g1_predictive_samples_csv_path <- file.path(out_dir, "g0g1_ploidy_predictive_samples.csv")
  g0g1_density_grid_csv_path <- file.path(out_dir, "g0g1_ploidy_density_grid.csv")
  g0g1_bar_plot_path <- file.path(out_dir, "g0g1_ploidy_stacked_bar.png")
  g0g1_density_plot_path <- file.path(out_dir, "g0g1_ploidy_density.png")

  saveRDS(fit_parts$distribution_object, distribution_rds_path)
  write.csv(fit_parts$distribution_object$sample_summary_tbl, sample_summary_csv_path, row.names = FALSE)
  write.csv(fit_parts$distribution_object$component_summary_tbl, component_summary_csv_path, row.names = FALSE)
  write.csv(fit_parts$distribution_object$draw_summary_tbl, draw_summary_csv_path, row.names = FALSE)
  write.csv(fit_parts$distribution_object$component_draw_tbl, component_draw_csv_path, row.names = FALSE)
  write.csv(fit_parts$distribution_object$overlay_tbl, overlay_csv_path, row.names = FALSE)

  state_parts <- build_flow_state_distribution_outputs(fit_parts$distribution_object)
  write.csv(state_parts$state_summary_tbl, state_summary_csv_path, row.names = FALSE)
  write.csv(state_parts$g0g1_draw_tbl, g0g1_draw_csv_path, row.names = FALSE)
  write.csv(state_parts$g0g1_summary_tbl, g0g1_summary_csv_path, row.names = FALSE)
  write.csv(state_parts$g0g1_predictive_samples_tbl, g0g1_predictive_samples_csv_path, row.names = FALSE)
  write.csv(state_parts$g0g1_density_grid_tbl, g0g1_density_grid_csv_path, row.names = FALSE)

  ggplot2::ggsave(
    overlay_plot_path,
    plot_flow_histogram_overlays(
      fit_parts$distribution_object$overlay_tbl,
      fit_parts$distribution_object$sample_summary_tbl
    ),
    width = 14,
    height = 10,
    dpi = 220
  )
  ggplot2::ggsave(
    state_plot_path,
    plot_flow_state_distribution(state_parts$state_summary_tbl),
    width = 13,
    height = 7,
    dpi = 220
  )
  ggplot2::ggsave(
    g0g1_bar_plot_path,
    plot_g0g1_ploidy_stacked_bar(state_parts$g0g1_summary_tbl),
    width = 13,
    height = 7,
    dpi = 220
  )
  ggplot2::ggsave(
    g0g1_density_plot_path,
    plot_g0g1_ploidy_density(state_parts$g0g1_predictive_density_tbl, state_parts$g0g1_summary_tbl),
    width = 14,
    height = 10,
    dpi = 220
  )

  c(
    fit_parts,
    list(
      out_dir = out_dir,
      distribution_rds_path = distribution_rds_path,
      sample_summary_csv_path = sample_summary_csv_path,
      component_summary_csv_path = component_summary_csv_path,
      draw_summary_csv_path = draw_summary_csv_path,
      component_draw_csv_path = component_draw_csv_path,
      overlay_csv_path = overlay_csv_path,
      overlay_plot_path = overlay_plot_path,
      state_summary_csv_path = state_summary_csv_path,
      state_plot_path = state_plot_path,
      g0g1_draw_csv_path = g0g1_draw_csv_path,
      g0g1_summary_csv_path = g0g1_summary_csv_path,
      g0g1_predictive_samples_csv_path = g0g1_predictive_samples_csv_path,
      g0g1_density_grid_csv_path = g0g1_density_grid_csv_path,
      g0g1_bar_plot_path = g0g1_bar_plot_path,
      g0g1_density_plot_path = g0g1_density_plot_path
    )
  )
}

plot_flow_histogram_draw_overlays <- function(draw_overlay_tbl) {
  plot_tbl <- draw_overlay_tbl %>%
    mutate(
      sample_label = factor(sample_name, levels = unique(sample_name)),
      draw_label = paste0("Draw ", draw_id, " (", assigned_state, ")")
    )

  ggplot(plot_tbl, aes(x = bin_mid)) +
    geom_rect(
      aes(xmin = bin_left, xmax = bin_right, ymin = 0, ymax = observed_density),
      fill = "grey82",
      color = "white",
      linewidth = 0.12
    ) +
    geom_line(
      aes(y = fitted_density, color = draw_label, group = draw_id),
      linewidth = 0.55,
      alpha = 0.9
    ) +
    facet_wrap(~ sample_label, scales = "free_y", ncol = 2) +
    labs(
      title = "Trial flow-only histogram overlays",
      subtitle = "Observed above-threshold DNA-A histograms with one fitted curve per posterior draw.",
      x = "DNA-A",
      y = "Density",
      color = "Posterior draw"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
}

write_flow_posterior_histogram_trial_outputs <- function(
  project_root,
  spec_name = "adjacent_8n_delta_dna",
  output_category = "extensions",
  sample_names,
  draw_ids = 1:3,
  threshold = NULL,
  min_bins = 60,
  max_bins = 160
) {
  spec_registry <- c(get_ablation_registry(), get_extension_registry())
  if (!spec_name %in% names(spec_registry)) {
    stop(sprintf("Unknown spec `%s`.", spec_name))
  }

  spec <- spec_registry[[spec_name]]
  prepared_path <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", "prepared_input.rds")
  prepared <- readRDS(prepared_path)
  threshold <- if (is.null(threshold)) prepared$threshold else threshold

  if (missing(sample_names) || !length(sample_names)) {
    stop("Please provide at least one sample name for the trial fit.")
  }

  missing_samples <- setdiff(sample_names, prepared$input_tbl$sample_name)
  if (length(missing_samples)) {
    stop(sprintf("Unknown sample(s): %s", paste(missing_samples, collapse = ", ")))
  }

  nuts_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", output_category, spec$name, "nuts")
  draws_matrix <- readRDS(file.path(nuts_dir, "draws_matrix.rds"))
  anchor_draw_tbl <- build_flow_anchor_draws(draws_matrix, prepared, spec) %>%
    filter(sample_name %in% sample_names, draw_id %in% draw_ids)

  dna_path <- file.path(project_root, "processed_data", "hypoxia-sum159", "filtered_dna_area_vectors.rds")
  raw_vectors <- readRDS(dna_path)$raw
  sample_histograms <- lapply(sample_names, function(sample_name) {
    build_flow_event_histogram(raw_vectors[[sample_name]], threshold = threshold, min_bins = min_bins, max_bins = max_bins)
  })
  names(sample_histograms) <- sample_names

  draw_overlay_tbl <- bind_rows(lapply(sample_names, function(sample_name) {
    hist_tbl <- sample_histograms[[sample_name]]
    sample_anchor_tbl <- anchor_draw_tbl %>% filter(sample_name == !!sample_name) %>% arrange(draw_id)
    bind_rows(lapply(seq_len(nrow(sample_anchor_tbl)), function(i) {
      anchor_row <- sample_anchor_tbl[i, , drop = FALSE]
      fit_obj <- fit_fixed_mean_histogram_mixture(hist_tbl, build_flow_component_spec(anchor_row))
      fit_obj$overlay_tbl %>%
        transmute(
          sample_name = sample_name,
          draw_id = anchor_row$draw_id[[1]],
          assigned_state = anchor_row$assigned_state[[1]],
          bin_left = hist_tbl$bin_left,
          bin_right = hist_tbl$bin_right,
          bin_mid,
          observed_count,
          observed_density,
          fitted_count,
          fitted_density
        )
    }))
  }))

  component_draw_tbl <- bind_rows(lapply(sample_names, function(sample_name) {
    hist_tbl <- sample_histograms[[sample_name]]
    sample_anchor_tbl <- anchor_draw_tbl %>% filter(sample_name == !!sample_name) %>% arrange(draw_id)
    bind_rows(lapply(seq_len(nrow(sample_anchor_tbl)), function(i) {
      anchor_row <- sample_anchor_tbl[i, , drop = FALSE]
      fit_obj <- fit_fixed_mean_histogram_mixture(hist_tbl, build_flow_component_spec(anchor_row))
      fit_obj$component_tbl %>%
        transmute(
          sample_name = sample_name,
          draw_id = anchor_row$draw_id[[1]],
          assigned_state = anchor_row$assigned_state[[1]],
          component_index,
          component_label,
          constrained_location = mean,
          weight,
          sigma
        )
    }))
  }))

  sample_summary_tbl <- anchor_draw_tbl %>%
    group_by(sample_name, condition, latest_match_date, relative_day) %>%
    summarise(
      posterior_prob_4n_subset = mean(prob_4n, na.rm = TRUE),
      lower_anchor_median_subset = stats::median(lower_anchor, na.rm = TRUE),
      upper_anchor_median_subset = stats::median(upper_anchor, na.rm = TRUE),
      draw_ids = paste(sort(unique(draw_id)), collapse = ","),
      .groups = "drop"
    )

  out_dir <- file.path(nuts_dir, "flow_only_distribution_trial")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  draw_overlay_csv_path <- file.path(out_dir, "trial_draw_overlay_counts.csv")
  component_draw_csv_path <- file.path(out_dir, "trial_component_draws.csv")
  sample_summary_csv_path <- file.path(out_dir, "trial_sample_summary.csv")
  overlay_plot_path <- file.path(out_dir, "trial_flow_distribution_overlays.png")

  write.csv(draw_overlay_tbl, draw_overlay_csv_path, row.names = FALSE)
  write.csv(component_draw_tbl, component_draw_csv_path, row.names = FALSE)
  write.csv(sample_summary_tbl, sample_summary_csv_path, row.names = FALSE)
  ggplot2::ggsave(
    overlay_plot_path,
    plot_flow_histogram_draw_overlays(draw_overlay_tbl),
    width = 12,
    height = 7,
    dpi = 220
  )

  list(
    spec = spec,
    nuts_dir = nuts_dir,
    out_dir = out_dir,
    draw_overlay_tbl = draw_overlay_tbl,
    component_draw_tbl = component_draw_tbl,
    sample_summary_tbl = sample_summary_tbl,
    draw_overlay_csv_path = draw_overlay_csv_path,
    component_draw_csv_path = component_draw_csv_path,
    sample_summary_csv_path = sample_summary_csv_path,
    overlay_plot_path = overlay_plot_path
  )
}
