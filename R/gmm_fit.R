require_mclust <- function() {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for GMM fitting. Install it with install.packages('mclust').")
  }
  suppressWarnings(
    suppressPackageStartupMessages(
      library("mclust", character.only = TRUE)
    )
  )
}

require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for GMM fit plots. Install it with install.packages('ggplot2').")
  }
}

clean_univariate_values <- function(values, min_n = 10L) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]

  if (length(values) < min_n) {
    stop(sprintf("Need at least %d finite values for GMM fitting; found %d.", min_n, length(values)))
  }
  if (length(unique(values)) < 2L) {
    stop("Need at least two distinct finite values for GMM fitting.")
  }

  values
}

normalize_filtered_dna_export <- function(payload, export_path = "<in-memory>") {
  if (!is.list(payload) || is.null(payload$raw) || !is.list(payload$raw)) {
    stop(sprintf("Expected '%s' to contain a list with a named '$raw' sample list.", export_path))
  }

  sample_names <- names(payload$raw)
  if (is.null(sample_names) || any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop(sprintf("Invalid or duplicated sample names in '$raw' for '%s'.", export_path))
  }

  transformed <- payload$transformed
  if (is.null(transformed)) {
    transformed <- payload$raw
  }
  if (!is.list(transformed) || !identical(sample_names, names(transformed))) {
    stop(sprintf("Raw and transformed sample names do not align in '%s'.", export_path))
  }

  list(raw = payload$raw, transformed = transformed)
}

read_filtered_dna_export <- function(export_path) {
  export_path <- normalizePath(export_path, winslash = "/", mustWork = TRUE)
  normalize_filtered_dna_export(readRDS(export_path), export_path)
}

find_filtered_dna_exports <- function(processed_root = "processed_data") {
  files <- list.files(
    processed_root,
    pattern = "^filtered_dna_area_vectors\\.rds$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  sort(normalizePath(files, winslash = "/", mustWork = TRUE))
}

get_filtered_dna_sample_values <- function(export_payload, sample_name, scale = c("raw", "transformed")) {
  scale <- match.arg(scale)
  export_payload <- normalize_filtered_dna_export(export_payload)

  sample_names <- names(export_payload[[scale]])
  if (!sample_name %in% sample_names) {
    stop(sprintf(
      "Sample '%s' is not present. Available samples: %s",
      sample_name,
      paste(sample_names, collapse = ", ")
    ))
  }

  export_payload[[scale]][[sample_name]]
}

criterion_grid_table <- function(values, value_name) {
  value_matrix <- as.matrix(values)
  grid <- as.data.frame(as.table(value_matrix), stringsAsFactors = FALSE)
  names(grid) <- c("G", "modelName", value_name)
  grid$G <- as.integer(as.character(grid$G))
  grid$modelName <- as.character(grid$modelName)
  grid[[value_name]] <- as.numeric(grid[[value_name]])
  grid[order(grid$modelName, grid$G), , drop = FALSE]
}

mclust_selection_grid <- function(bic_values, icl_values) {
  bic_grid <- criterion_grid_table(bic_values, "BIC")
  icl_grid <- criterion_grid_table(icl_values, "ICL")
  merge(bic_grid, icl_grid, by = c("G", "modelName"), all = TRUE, sort = FALSE)
}

mclust_model_icl <- function(model) {
  z <- model$z
  if (is.null(z)) {
    return(model$bic[[1]])
  }
  if (!is.matrix(z) || any(!is.finite(z))) {
    return(NA_real_)
  }
  class_idx <- max.col(z, ties.method = "first")
  if (any(!is.finite(class_idx))) {
    return(NA_real_)
  }
  hard <- matrix(0, nrow = nrow(z), ncol = ncol(z))
  hard[cbind(seq_len(nrow(z)), class_idx)] <- 1
  model$bic[[1]] + 2 * sum(hard * ifelse(z > 0, log(z), 0))
}

compute_mclust_icl_grid <- function(values, bic_values, G, modelNames) {
  bic_grid <- criterion_grid_table(bic_values, "BIC")
  bic_matrix <- as.matrix(bic_values)
  icl_matrix <- matrix(
    NA_real_,
    nrow = nrow(bic_matrix),
    ncol = ncol(bic_matrix),
    dimnames = dimnames(bic_matrix)
  )

  finite_grid <- bic_grid[is.finite(bic_grid$BIC), , drop = FALSE]
  for (i in seq_len(nrow(finite_grid))) {
    g <- finite_grid$G[[i]]
    model_name <- finite_grid$modelName[[i]]
    icl_value <- tryCatch({
      model <- mclust::mclustModel(
        data = values,
        BICvalues = bic_values,
        G = g,
        modelNames = model_name
      )
      mclust_model_icl(model)
    }, error = function(e) NA_real_)
    icl_matrix[as.character(g), model_name] <- icl_value
  }

  structure(
    icl_matrix,
    G = G,
    modelNames = modelNames,
    criterion = "ICL",
    class = "mclustICL"
  )
}

best_criterion_row <- function(selection_grid, criterion = c("ICL", "BIC")) {
  criterion <- match.arg(criterion)
  finite_grid <- selection_grid[is.finite(selection_grid[[criterion]]), , drop = FALSE]
  if (nrow(finite_grid) == 0L) {
    stop(sprintf("No finite %s values were produced by mclust.", criterion))
  }
  finite_grid[which.max(finite_grid[[criterion]]), , drop = FALSE]
}

deterministic_hc_subset <- function(n, max_n = 2000L) {
  max_n <- as.integer(max_n)
  if (!is.finite(max_n) || max_n < 1L || n <= max_n) {
    return(NULL)
  }
  unique(as.integer(round(seq(1, n, length.out = max_n))))
}

build_mclust_hc_initialization <- function(values, modelNames = "E", hc_subset = 2000L, hc_modelName = "E") {
  hc_modelName <- match.arg(hc_modelName, choices = c("E", "V"))

  subset_idx <- deterministic_hc_subset(length(values), max_n = hc_subset)
  hc_values <- if (is.null(subset_idx)) values else values[subset_idx]
  hc_pairs <- mclust::hc(hc_values, modelName = hc_modelName)

  if (is.null(subset_idx)) {
    list(hcPairs = hc_pairs)
  } else {
    list(hcPairs = hc_pairs, subset = subset_idx)
  }
}

fit_univariate_gmm_bic <- function(
    values,
    G = 1:16,
    modelNames = "E",
    criterion = c("ICL", "BIC"),
    min_n = 10L,
    control = mclust::emControl(),
    initialization = NULL,
    hc_subset = 2000L,
    hc_modelName = "E",
    verbose = FALSE) {
  require_mclust()

  values <- clean_univariate_values(values, min_n = min_n)
  G <- sort(unique(as.integer(G)))
  G <- G[is.finite(G) & G >= 1L]
  if (length(G) == 0L) {
    stop("'G' must contain at least one positive integer.")
  }
  modelNames <- match.arg(modelNames, choices = c("E", "V"), several.ok = TRUE)
  criterion <- match.arg(criterion)

  if (is.null(initialization)) {
    initialization <- build_mclust_hc_initialization(
      values = values,
      modelNames = modelNames,
      hc_subset = hc_subset,
      hc_modelName = hc_modelName
    )
  }

  bic_values <- mclust::mclustBIC(
    data = values,
    G = G,
    modelNames = modelNames,
    control = control,
    initialization = initialization,
    verbose = verbose
  )
  icl_values <- compute_mclust_icl_grid(values, bic_values, G = G, modelNames = modelNames)
  selection_grid <- mclust_selection_grid(bic_values, icl_values)
  best <- best_criterion_row(selection_grid, criterion = criterion)
  model <- mclust::mclustModel(
    data = values,
    BICvalues = bic_values,
    G = best$G[[1]],
    modelNames = best$modelName[[1]]
  )

  list(
    values = values,
    bic = bic_values,
    icl = icl_values,
    selection_grid = selection_grid,
    bic_grid = criterion_grid_table(bic_values, "BIC"),
    icl_grid = criterion_grid_table(icl_values, "ICL"),
    best = best,
    criterion = criterion,
    model = model,
    settings = list(
      G = G,
      modelNames = modelNames,
      criterion = criterion,
      min_n = min_n,
      hc_subset = hc_subset,
      hc_modelName = hc_modelName,
      initialization = initialization
    )
  )
}

component_summary <- function(fit) {
  model <- fit$model
  params <- model$parameters
  variance <- params$variance$sigmasq
  if (length(variance) == 1L) {
    variance <- rep(variance, model$G)
  }

  out <- data.frame(
    component = seq_len(model$G),
    weight = as.numeric(params$pro),
    mean = as.numeric(params$mean),
    variance = as.numeric(variance),
    stringsAsFactors = FALSE
  )
  out$sd <- sqrt(out$variance)
  out[order(out$mean), , drop = FALSE]
}

classification_table <- function(fit) {
  z <- fit$model$z
  component <- max.col(z, ties.method = "first")
  probability <- z[cbind(seq_len(nrow(z)), component)]

  data.frame(
    event_index = seq_along(fit$values),
    value = fit$values,
    component = component,
    probability = probability,
    uncertainty = 1 - probability,
    stringsAsFactors = FALSE
  )
}

gmm_fit_summary <- function(fit, sample_name = NA_character_, dataset_id = NA_character_, scale = NA_character_) {
  values <- fit$values
  model <- fit$model
  data.frame(
    dataset_id = dataset_id,
    sample_name = sample_name,
    scale = scale,
    n = length(values),
    min = min(values),
    median = stats::median(values),
    max = max(values),
    modelName = model$modelName[[1]],
    G = model$G[[1]],
    bic = model$bic[[1]],
    icl = fit$best$ICL[[1]],
    selection_criterion = fit$criterion %||% NA_character_,
    selection_value = fit$best[[fit$criterion %||% "BIC"]][[1]],
    loglik = model$loglik[[1]],
    stringsAsFactors = FALSE
  )
}

fit_filtered_dna_sample_gmm <- function(
    export_path,
    sample_name,
    scale = c("raw", "transformed"),
    G = 1:16,
    modelNames = "E",
    criterion = c("ICL", "BIC"),
    hc_subset = 2000L,
    hc_modelName = "E",
    min_n = 10L,
    verbose = FALSE) {
  scale <- match.arg(scale)
  criterion <- match.arg(criterion)
  payload <- read_filtered_dna_export(export_path)
  values <- get_filtered_dna_sample_values(payload, sample_name = sample_name, scale = scale)
  fit <- fit_univariate_gmm_bic(
    values = values,
    G = G,
    modelNames = modelNames,
    criterion = criterion,
    hc_subset = hc_subset,
    hc_modelName = hc_modelName,
    min_n = min_n,
    verbose = verbose
  )
  fit$source <- list(
    export_path = normalizePath(export_path, winslash = "/", mustWork = TRUE),
    sample_name = sample_name,
    scale = scale
  )
  fit
}

fit_filtered_dna_export_gmms <- function(
    export_path,
    sample_names = NULL,
    scale = c("raw", "transformed"),
    G = 1:16,
    modelNames = "E",
    criterion = c("ICL", "BIC"),
    hc_subset = 2000L,
    hc_modelName = "E",
    min_n = 10L,
    verbose = FALSE,
    progress = interactive()) {
  scale <- match.arg(scale)
  criterion <- match.arg(criterion)
  payload <- read_filtered_dna_export(export_path)
  available_samples <- names(payload[[scale]])

  if (is.null(sample_names)) {
    sample_names <- available_samples
  }
  missing_samples <- setdiff(sample_names, available_samples)
  if (length(missing_samples) > 0L) {
    stop(sprintf("Samples are not present in '%s': %s", export_path, paste(missing_samples, collapse = ", ")))
  }

  fits <- setNames(vector("list", length(sample_names)), sample_names)
  for (sample_name in sample_names) {
    if (isTRUE(progress)) {
      message(sprintf("Fitting %s", sample_name))
    }
    values <- get_filtered_dna_sample_values(payload, sample_name = sample_name, scale = scale)
    fits[[sample_name]] <- fit_univariate_gmm_bic(
      values = values,
      G = G,
      modelNames = modelNames,
      criterion = criterion,
      hc_subset = hc_subset,
      hc_modelName = hc_modelName,
      min_n = min_n,
      verbose = verbose
    )
    fits[[sample_name]]$source <- list(
      export_path = normalizePath(export_path, winslash = "/", mustWork = TRUE),
      sample_name = sample_name,
      scale = scale
    )
  }

  fits
}

summarize_gmm_fits <- function(fits, dataset_id = NA_character_, scale = NA_character_) {
  parts <- Map(
    f = function(fit, sample_name) {
      fit_scale <- fit$source$scale %||% scale
      gmm_fit_summary(fit, sample_name = sample_name, dataset_id = dataset_id, scale = fit_scale)
    },
    fit = fits,
    sample_name = names(fits)
  )
  do.call(rbind, parts)
}

gmm_density_table <- function(fit, n = 2048L) {
  values <- fit$values
  model <- fit$model
  params <- model$parameters
  variance <- params$variance$sigmasq
  if (length(variance) == 1L) {
    variance <- rep(variance, model$G)
  }

  x <- seq(min(values), max(values), length.out = n)
  parts <- lapply(seq_len(model$G), function(component) {
    data.frame(
      x = x,
      component = component,
      density = params$pro[[component]] * stats::dnorm(
        x,
        mean = params$mean[[component]],
        sd = sqrt(variance[[component]])
      ),
      stringsAsFactors = FALSE
    )
  })
  component_tbl <- do.call(rbind, parts)
  mixture_tbl <- stats::aggregate(density ~ x, data = component_tbl, FUN = sum)
  mixture_tbl$component <- "mixture"

  list(
    components = component_tbl,
    mixture = mixture_tbl
  )
}

plot_gmm_fit <- function(
    fit,
    sample_name = fit$source$sample_name %||% NA_character_,
    dataset_id = NA_character_,
    scale = fit$source$scale %||% NA_character_,
    bins = 256L,
    density_n = 2048L,
    component_color = "#2C7FB8",
    mixture_color = "#D95F02") {
  require_ggplot2()

  values_tbl <- data.frame(value = fit$values)
  densities <- gmm_density_table(fit, n = density_n)
  title <- if (is.na(sample_name)) "GMM fit" else sprintf("GMM fit: %s", sample_name)
  subtitle <- sprintf(
    "%s scale; selected %s, G=%d, %s=%.1f",
    scale,
    fit$model$modelName[[1]],
    fit$model$G[[1]],
    fit$criterion %||% "BIC",
    fit$best[[fit$criterion %||% "BIC"]]
  )
  if (!is.na(dataset_id)) {
    subtitle <- sprintf("%s; %s", dataset_id, subtitle)
  }

  ggplot2::ggplot(values_tbl, ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = bins,
      fill = "grey82",
      color = "grey55",
      linewidth = 0.15
    ) +
    ggplot2::geom_line(
      data = densities$components,
      ggplot2::aes(x = x, y = density, group = component),
      color = component_color,
      linewidth = 0.45,
      alpha = 0.65,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = densities$mixture,
      ggplot2::aes(x = x, y = density),
      color = mixture_color,
      linewidth = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Filtered DNA-A",
      y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

write_gmm_fit_plot <- function(
    fit,
    output_path,
    sample_name = fit$source$sample_name %||% NA_character_,
    dataset_id = NA_character_,
    scale = fit$source$scale %||% NA_character_,
    bins = 256L,
    width = 8,
    height = 5,
    dpi = 150) {
  require_ggplot2()
  plot <- plot_gmm_fit(
    fit = fit,
    sample_name = sample_name,
    dataset_id = dataset_id,
    scale = scale,
    bins = bins
  )
  ggplot2::ggsave(
    filename = output_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  invisible(output_path)
}

write_gmm_fit_outputs <- function(
    fit,
    output_dir,
    sample_name = fit$source$sample_name %||% NA_character_,
    dataset_id = NA_character_,
    scale = fit$source$scale %||% NA_character_,
    include_classification = TRUE,
    include_plot = TRUE,
    plot_bins = 256L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  summary_tbl <- gmm_fit_summary(fit, sample_name = sample_name, dataset_id = dataset_id, scale = scale)
  components_tbl <- component_summary(fit)

  saveRDS(fit, file.path(output_dir, "gmm_fit.rds"))
  utils::write.csv(summary_tbl, file.path(output_dir, "summary.csv"), row.names = FALSE)
  utils::write.csv(fit$selection_grid, file.path(output_dir, "selection_grid.csv"), row.names = FALSE)
  utils::write.csv(fit$bic_grid, file.path(output_dir, "bic_grid.csv"), row.names = FALSE)
  utils::write.csv(fit$icl_grid, file.path(output_dir, "icl_grid.csv"), row.names = FALSE)
  utils::write.csv(components_tbl, file.path(output_dir, "components.csv"), row.names = FALSE)

  if (isTRUE(include_classification)) {
    utils::write.csv(classification_table(fit), file.path(output_dir, "classification.csv"), row.names = FALSE)
  }
  if (isTRUE(include_plot)) {
    write_gmm_fit_plot(
      fit = fit,
      output_path = file.path(output_dir, "fit_plot.png"),
      sample_name = sample_name,
      dataset_id = dataset_id,
      scale = scale,
      bins = plot_bins
    )
  }

  invisible(list(
    fit = file.path(output_dir, "gmm_fit.rds"),
    summary = file.path(output_dir, "summary.csv"),
    selection_grid = file.path(output_dir, "selection_grid.csv"),
    bic_grid = file.path(output_dir, "bic_grid.csv"),
    icl_grid = file.path(output_dir, "icl_grid.csv"),
    components = file.path(output_dir, "components.csv"),
    classification = if (isTRUE(include_classification)) file.path(output_dir, "classification.csv") else NA_character_,
    plot = if (isTRUE(include_plot)) file.path(output_dir, "fit_plot.png") else NA_character_
  ))
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}
