residual_band_score <- function(residuals, clusters) {
  centers <- vapply(clusters, mean, numeric(1))
  spreads <- vapply(clusters, stats::mad, numeric(1), constant = 1)
  sizes <- vapply(clusters, length, integer(1))

  if (length(centers) < 2L || any(sizes < 5L)) {
    return(-Inf)
  }

  separation <- stats::median(diff(sort(centers)))
  weighted_spread <- stats::weighted.mean(pmax(spreads, .Machine$double.eps), sizes)
  separation / weighted_spread
}

residual_band_fit <- function(
  x,
  y,
  slope_min,
  slope_max,
  slope_steps,
  x_trim,
  K,
  central_prob,
  nstart
) {
  keep <- stats::complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < K * 10L) {
    stop("Not enough complete events to estimate residualBandGate.")
  }

  x_limits <- stats::quantile(x, probs = c(x_trim, 1 - x_trim), names = FALSE)
  fit_keep <- x >= x_limits[[1]] & x <= x_limits[[2]]
  x_fit <- x[fit_keep]
  y_fit <- y[fit_keep]

  if (length(x_fit) < K * 10L) {
    stop("Not enough x-trimmed events to estimate residualBandGate.")
  }

  slopes <- seq(slope_min, slope_max, length.out = slope_steps)
  best <- NULL

  for (slope in slopes) {
    residuals <- y_fit - slope * x_fit
    km <- stats::kmeans(residuals, centers = K, nstart = nstart)
    clusters <- split(residuals, km$cluster)
    score <- residual_band_score(residuals, clusters)

    if (is.null(best) || score > best$score) {
      best <- list(slope = slope, kmeans = km, score = score)
    }
  }

  residuals <- y - best$slope * x
  cluster_centers <- best$kmeans$centers[, 1]
  assigned <- vapply(residuals, function(r) {
    which.min(abs(r - cluster_centers))
  }, integer(1))
  cluster_sizes <- tabulate(assigned, nbins = length(cluster_centers))
  singlet_cluster <- which.max(cluster_sizes)
  singlet_residuals <- residuals[assigned == singlet_cluster]
  alpha <- (1 - central_prob) / 2
  residual_limits <- stats::quantile(
    singlet_residuals,
    probs = c(alpha, 1 - alpha),
    names = FALSE
  )

  list(
    slope = best$slope,
    residual_limits = residual_limits,
    x_limits = range(x, finite = TRUE),
    score = best$score,
    singlet_count = length(singlet_residuals)
  )
}

.residualBandGate <- function(
  fr,
  pp_res = NULL,
  channels,
  slope_min = 0.6,
  slope_max = 1.4,
  slope_steps = 81,
  x_trim = 0.01,
  K = 3,
  central_prob = 0.95,
  nstart = 10,
  ...
) {
  if (length(channels) != 2L) {
    stop("residualBandGate requires exactly two channels: x then y.")
  }

  expr <- flowCore::exprs(fr)
  x <- expr[, channels[[1]]]
  y <- expr[, channels[[2]]]
  fit <- residual_band_fit(
    x = x,
    y = y,
    slope_min = slope_min,
    slope_max = slope_max,
    slope_steps = slope_steps,
    x_trim = x_trim,
    K = K,
    central_prob = central_prob,
    nstart = nstart
  )

  x_min <- fit$x_limits[[1]]
  x_max <- fit$x_limits[[2]]
  r_lo <- fit$residual_limits[[1]]
  r_hi <- fit$residual_limits[[2]]
  slope <- fit$slope

  boundaries <- matrix(
    c(
      x_min, slope * x_min + r_lo,
      x_max, slope * x_max + r_lo,
      x_max, slope * x_max + r_hi,
      x_min, slope * x_min + r_hi
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(boundaries) <- channels

  flowCore::polygonGate(.gate = boundaries, filterId = "residualBandGate")
}

register_custom_gates <- function() {
  invisible(openCyto::register_plugins(
    fun = .residualBandGate,
    methodName = "residualBandGate"
  ))
}
