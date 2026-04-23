source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

estimate_window_peaks <- function(
  values,
  threshold = 12000,
  tumor_peak_min_separation = 35000,
  tumor_peak_min_prominence = 0
) {
  x <- values[is.finite(values)]
  if (length(x) < 10 || length(unique(x)) < 2) {
    return(tibble(
      cen_peak = NA_real_,
      tumor_peak_1 = NA_real_,
      tumor_peak_2 = NA_real_,
      n_window = length(x),
      n_below = sum(x < threshold, na.rm = TRUE),
      n_above = sum(x > threshold, na.rm = TRUE)
    ))
  }

  below_peak <- summarise_threshold_peak(x, side = "below", threshold = threshold)
  tumor_candidates <- find_threshold_peak_candidates(x, side = "above", threshold = threshold)
  tumor_selected <- select_separated_peaks(
    tumor_candidates,
    n_keep = 2,
    min_separation = tumor_peak_min_separation,
    min_prominence = tumor_peak_min_prominence
  )

  tibble(
    cen_peak = below_peak$peak_location,
    tumor_peak_1 = if (nrow(tumor_selected) >= 1) tumor_selected$peak_location[[1]] else NA_real_,
    tumor_peak_2 = if (nrow(tumor_selected) >= 2) tumor_selected$peak_location[[2]] else NA_real_,
    n_window = length(x),
    n_below = sum(x < threshold, na.rm = TRUE),
    n_above = sum(x > threshold, na.rm = TRUE)
  )
}

build_index_windows <- function(n, window_size, step_size) {
  stopifnot(n >= 1, window_size >= 1, step_size >= 1)
  if (n <= window_size) {
    return(tibble(window_id = 1L, start_idx = 1L, end_idx = as.integer(n)))
  }

  starts <- seq.int(1L, n - window_size + 1L, by = step_size)
  ends <- starts + window_size - 1L
  out <- tibble(
    window_id = seq_along(starts),
    start_idx = as.integer(starts),
    end_idx = as.integer(ends)
  )

  if (tail(out$end_idx, 1) < n) {
    out <- bind_rows(
      out,
      tibble(
        window_id = nrow(out) + 1L,
        start_idx = as.integer(n - window_size + 1L),
        end_idx = as.integer(n)
      )
    ) %>%
      distinct(start_idx, end_idx, .keep_all = TRUE) %>%
      mutate(window_id = row_number())
  }

  out
}

analyze_event_windows <- function(
  sample_tbl,
  threshold,
  event_window_size,
  event_step_size,
  tumor_peak_min_separation,
  tumor_peak_min_prominence
) {
  bind_rows(lapply(seq_len(nrow(sample_tbl)), function(i) {
    sample_id <- sample_tbl$sample_id[[i]]
    x <- sample_tbl$raw_dna[[i]]
    n <- length(x)
    windows <- build_index_windows(n, event_window_size, event_step_size)

    bind_rows(lapply(seq_len(nrow(windows)), function(j) {
      window_row <- windows[j, , drop = FALSE]
      idx <- window_row$start_idx[[1]]:window_row$end_idx[[1]]
      peak_tbl <- estimate_window_peaks(
        values = x[idx],
        threshold = threshold,
        tumor_peak_min_separation = tumor_peak_min_separation,
        tumor_peak_min_prominence = tumor_peak_min_prominence
      )

      tibble(
        sample_id = sample_id,
        window_mode = "event_index",
        window_id = window_row$window_id[[1]],
        start_idx = window_row$start_idx[[1]],
        end_idx = window_row$end_idx[[1]],
        window_mid_idx = mean(c(window_row$start_idx[[1]], window_row$end_idx[[1]])),
        window_mid_fraction = mean(c(window_row$start_idx[[1]], window_row$end_idx[[1]])) / n,
        total_events = n
      ) %>%
        bind_cols(peak_tbl)
    }))
  }))
}

discover_hypoxia_flow_files <- function(project_root, sample_ids) {
  flow_dir <- file.path(project_root, "data", "Hypoxia_SUM159")
  if (!dir.exists(flow_dir)) {
    return(NULL)
  }

  flow_files <- list.files(flow_dir, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
  key <- tools::file_path_sans_ext(basename(flow_files))
  flow_map <- setNames(normalizePath(flow_files, winslash = "/", mustWork = TRUE), key)
  sample_ids <- sample_ids[sample_ids %in% names(flow_map)]
  if (!length(sample_ids)) {
    return(NULL)
  }

  flow_map[sample_ids]
}

try_load_time_vectors <- function(project_root, sample_tbl) {
  if (!requireNamespace("flowWorkspace", quietly = TRUE) || !requireNamespace("flowCore", quietly = TRUE)) {
    message("Skipping time-window analysis: flowWorkspace/flowCore not available.")
    return(NULL)
  }

  gs_dir <- file.path(project_root, "processed_data", "hypoxia-sum159", "gating_gs")
  if (!dir.exists(gs_dir)) {
    message("Skipping time-window analysis: gating set cache not found.")
    return(NULL)
  }

  flow_map <- discover_hypoxia_flow_files(project_root, sample_tbl$sample_id)
  if (is.null(flow_map) || !length(flow_map)) {
    message("Skipping time-window analysis: source .fcs files not found.")
    return(NULL)
  }

  gs <- flowWorkspace::load_gs(gs_dir)
  gs_names <- flowWorkspace::sampleNames(gs)
  gs_keys <- tools::file_path_sans_ext(basename(gs_names))
  names(gs_keys) <- gs_names

  time_parts <- lapply(seq_len(nrow(sample_tbl)), function(i) {
    sample_id <- sample_tbl$sample_id[[i]]
    gs_match <- names(gs_keys)[gs_keys %in% sample_id]
    flow_path <- flow_map[[sample_id]]
    if (!length(gs_match) || is.na(flow_path) || !file.exists(flow_path)) {
      return(NULL)
    }

    sn <- gs_match[[1]]
    keep_idx <- flowWorkspace::gh_pop_get_indices(gs[[sn]], "Singlets")
    ff <- flowCore::read.FCS(
      flow_path,
      transformation = FALSE,
      truncate_max_range = FALSE
    )
    expr <- flowCore::exprs(ff)
    if (!"Time" %in% colnames(expr)) {
      return(NULL)
    }

    dna_vec <- as.numeric(expr[keep_idx, "450-Violet C-A"])
    time_vec <- as.numeric(expr[keep_idx, "Time"])
    keep <- is.finite(dna_vec) & is.finite(time_vec)
    dna_vec <- dna_vec[keep]
    time_vec <- time_vec[keep]

    if (!length(dna_vec) || length(dna_vec) != length(sample_tbl$raw_dna[[i]])) {
      return(NULL)
    }

    # Convert from instrument ticks to seconds when the FCS metadata provides a timestep.
    kw <- flowCore::keyword(ff)
    time_step <- suppressWarnings(as.numeric(kw[["$TIMESTEP"]]))
    if (is.finite(time_step) && time_step > 0) {
      time_vec <- time_vec * time_step
    }

    tibble(sample_id = sample_id, time_seconds = time_vec)
  })

  time_tbl <- bind_rows(time_parts)
  if (!nrow(time_tbl)) {
    return(NULL)
  }

  split(time_tbl$time_seconds, time_tbl$sample_id)
}

build_time_windows <- function(time_seconds, approx_events_per_window, approx_events_per_step) {
  if (!length(time_seconds) || any(!is.finite(time_seconds))) {
    return(tibble())
  }

  n <- length(time_seconds)
  total_duration <- max(time_seconds) - min(time_seconds)
  if (!is.finite(total_duration) || total_duration <= 0) {
    return(tibble())
  }

  time_window_width <- total_duration * min(approx_events_per_window, n) / n
  time_step_width <- total_duration * min(approx_events_per_step, n) / n
  time_start <- min(time_seconds)
  latest_start <- max(time_seconds) - time_window_width

  starts <- seq(from = time_start, to = latest_start, by = max(time_step_width, .Machine$double.eps))
  if (!length(starts)) {
    starts <- time_start
  }

  ends <- pmin(starts + time_window_width, max(time_seconds))
  out <- tibble(
    window_id = seq_along(starts),
    time_start_s = starts,
    time_end_s = ends
  )

  if (tail(out$time_end_s, 1) < max(time_seconds)) {
    out <- bind_rows(
      out,
      tibble(
        window_id = nrow(out) + 1L,
        time_start_s = max(time_seconds) - time_window_width,
        time_end_s = max(time_seconds)
      )
    ) %>%
      distinct(time_start_s, time_end_s, .keep_all = TRUE) %>%
      mutate(window_id = row_number())
  }

  out
}

analyze_time_windows <- function(
  sample_tbl,
  time_vectors,
  threshold,
  event_window_size,
  event_step_size,
  tumor_peak_min_separation,
  tumor_peak_min_prominence
) {
  if (is.null(time_vectors) || !length(time_vectors)) {
    return(tibble())
  }

  bind_rows(lapply(seq_len(nrow(sample_tbl)), function(i) {
    sample_id <- sample_tbl$sample_id[[i]]
    x <- sample_tbl$raw_dna[[i]]
    time_seconds <- time_vectors[[sample_id]]
    if (is.null(time_seconds) || length(time_seconds) != length(x)) {
      return(NULL)
    }

    windows <- build_time_windows(
      time_seconds = time_seconds,
      approx_events_per_window = event_window_size,
      approx_events_per_step = event_step_size
    )
    if (!nrow(windows)) {
      return(NULL)
    }

    time_origin <- min(time_seconds)
    total_duration <- max(time_seconds) - time_origin

    bind_rows(lapply(seq_len(nrow(windows)), function(j) {
      window_row <- windows[j, , drop = FALSE]
      keep <- time_seconds >= window_row$time_start_s[[1]] & time_seconds <= window_row$time_end_s[[1]]
      idx <- which(keep)
      if (!length(idx)) {
        return(NULL)
      }

      peak_tbl <- estimate_window_peaks(
        values = x[idx],
        threshold = threshold,
        tumor_peak_min_separation = tumor_peak_min_separation,
        tumor_peak_min_prominence = tumor_peak_min_prominence
      )

      tibble(
        sample_id = sample_id,
        window_mode = "time_seconds",
        window_id = window_row$window_id[[1]],
        start_idx = min(idx),
        end_idx = max(idx),
        window_mid_idx = mean(range(idx)),
        window_mid_fraction = if (total_duration > 0) {
          (mean(c(window_row$time_start_s[[1]], window_row$time_end_s[[1]])) - time_origin) / total_duration
        } else {
          NA_real_
        },
        total_events = length(x),
        time_start_s = window_row$time_start_s[[1]],
        time_end_s = window_row$time_end_s[[1]],
        time_mid_s = mean(c(window_row$time_start_s[[1]], window_row$time_end_s[[1]])),
        total_duration_s = total_duration
      ) %>%
        bind_cols(peak_tbl)
    }))
  }))
}

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
event_window_size <- if (length(args) >= 2) as.integer(args[[2]]) else 10000L
event_step_size <- if (length(args) >= 3) as.integer(args[[3]]) else 2500L
threshold <- if (length(args) >= 4) as.numeric(args[[4]]) else 12000
tumor_peak_min_separation <- if (length(args) >= 5) as.numeric(args[[5]]) else 35000
tumor_peak_min_prominence <- if (length(args) >= 6) as.numeric(args[[6]]) else 0

out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", "peak_drift")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_tbl <- read.csv(
  file.path(project_root, "processed_data", "hypoxia-sum159", "sample_metadata.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
dna_vectors <- readRDS(
  file.path(project_root, "processed_data", "hypoxia-sum159", "filtered_dna_area_vectors.rds")
)$raw

sample_tbl <- tibble(
  sample_id = names(dna_vectors),
  raw_dna = unname(dna_vectors)
) %>%
  inner_join(metadata_tbl, by = "sample_id") %>%
  mutate(sample_name = sample_id) %>%
  select(sample_id, sample_name, cell_line, condition, latest_match_date, relative_day, raw_dna)

overall_peak_tbl <- bind_rows(lapply(seq_len(nrow(sample_tbl)), function(i) {
  estimate_window_peaks(
    values = sample_tbl$raw_dna[[i]],
    threshold = threshold,
    tumor_peak_min_separation = tumor_peak_min_separation,
    tumor_peak_min_prominence = tumor_peak_min_prominence
  ) %>%
    mutate(sample_id = sample_tbl$sample_id[[i]])
})) %>%
  relocate(sample_id)

event_tbl <- analyze_event_windows(
  sample_tbl = sample_tbl,
  threshold = threshold,
  event_window_size = event_window_size,
  event_step_size = event_step_size,
  tumor_peak_min_separation = tumor_peak_min_separation,
  tumor_peak_min_prominence = tumor_peak_min_prominence
) %>%
  left_join(
    sample_tbl %>% select(sample_id, sample_name, cell_line, condition, latest_match_date, relative_day),
    by = "sample_id"
  ) %>%
  left_join(
    overall_peak_tbl %>%
      rename(
        overall_cen_peak = cen_peak,
        overall_tumor_peak_1 = tumor_peak_1,
        overall_tumor_peak_2 = tumor_peak_2,
        overall_n_window = n_window,
        overall_n_below = n_below,
        overall_n_above = n_above
      ),
    by = "sample_id"
  )

time_vectors <- try_load_time_vectors(project_root, sample_tbl)
time_tbl <- analyze_time_windows(
  sample_tbl = sample_tbl,
  time_vectors = time_vectors,
  threshold = threshold,
  event_window_size = event_window_size,
  event_step_size = event_step_size,
  tumor_peak_min_separation = tumor_peak_min_separation,
  tumor_peak_min_prominence = tumor_peak_min_prominence
) %>%
  left_join(
    sample_tbl %>% select(sample_id, sample_name, cell_line, condition, latest_match_date, relative_day),
    by = "sample_id"
  ) %>%
  left_join(
    overall_peak_tbl %>%
      rename(
        overall_cen_peak = cen_peak,
        overall_tumor_peak_1 = tumor_peak_1,
        overall_tumor_peak_2 = tumor_peak_2,
        overall_n_window = n_window,
        overall_n_below = n_below,
        overall_n_above = n_above
      ),
    by = "sample_id"
  )

write.csv(event_tbl, file.path(out_dir, "event_window_peak_drift.csv"), row.names = FALSE)
write.csv(time_tbl, file.path(out_dir, "time_window_peak_drift.csv"), row.names = FALSE)
write.csv(overall_peak_tbl, file.path(out_dir, "overall_peak_reference.csv"), row.names = FALSE)
write.csv(
  tibble(
    threshold = threshold,
    event_window_size = event_window_size,
    event_step_size = event_step_size,
    tumor_peak_min_separation = tumor_peak_min_separation,
    tumor_peak_min_prominence = tumor_peak_min_prominence,
    has_time_windows = nrow(time_tbl) > 0
  ),
  file.path(out_dir, "peak_drift_config.csv"),
  row.names = FALSE
)

cat("Saved peak drift outputs under:", out_dir, "\n")
cat("Event-window rows:", nrow(event_tbl), "\n")
cat("Time-window rows:", nrow(time_tbl), "\n")
