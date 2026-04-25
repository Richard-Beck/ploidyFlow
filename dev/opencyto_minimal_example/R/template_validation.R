load_template_table <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

template_dims <- function(template_table) {
  dims <- template_table$dims[template_table$dims != ""]
  dims <- unique(unlist(strsplit(dims, ",", fixed = TRUE)))
  trimws(dims[nzchar(dims)])
}

validate_template_dims <- function(template_table, available_channels, allowed_transforms) {
  dims <- template_dims(template_table)

  for (dim_name in dims) {
    parsed <- split_transformed_dim(dim_name, allowed_transforms)
    if (!parsed$base_dim %in% available_channels) {
      stop(
        "Template dim '",
        dim_name,
        "' refers to unknown channel '",
        parsed$base_dim,
        "'."
      )
    }
  }

  invisible(TRUE)
}

load_and_validate_template <- function(path, fs, config) {
  template_table <- load_template_table(path)
  validate_template_dims(template_table, colnames(fs), config$allowed_transforms)

  list(
    table = template_table,
    gt = openCyto::gatingTemplate(path)
  )
}
