normalize_metadata_string <- function(x) {
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

parse_sum159_hypoxia_sample_signature <- function(sample_name) {
  sample_id <- tools::file_path_sans_ext(basename(sample_name))
  normalized <- normalize_metadata_string(sample_id)

  parsed <- stringr::str_match(
    normalized,
    "^SAMPLE_SUM159_(2N|4N)_(C|O1|O2)_A([0-9]+)$"
  )

  tibble::tibble(
    sample_name = basename(sample_name),
    sample_id = sample_id,
    signature_ok = !is.na(parsed[, 1]),
    cell_line_key = ifelse(!is.na(parsed[, 1]), "SUM159", NA_character_),
    ploidy_key = parsed[, 2],
    oxygen_key = parsed[, 3],
    address_key = ifelse(!is.na(parsed[, 4]), paste0("A", parsed[, 4]), NA_character_)
  )
}

parse_sum159_passaging_signature <- function(passaging_tbl, id_col = "id", time_col = "date") {
  stopifnot(id_col %in% names(passaging_tbl))
  stopifnot(time_col %in% names(passaging_tbl))

  parse_match_time <- function(x) {
    parse_one <- function(value) {
      if (is.na(value)) {
        return(NA_real_)
      }

      value <- trimws(as.character(value))
      if (!nzchar(value)) {
        return(NA_real_)
      }

      formats <- c(
        "%Y-%m-%dT%H:%M:%OSZ",
        "%Y-%m-%d %H:%M:%OS",
        "%Y-%m-%d"
      )

      for (fmt in formats) {
        parsed <- as.POSIXct(value, format = fmt, tz = "UTC")
        if (!is.na(parsed)) {
          return(as.numeric(parsed))
        }
      }

      NA_real_
    }

    seconds <- vapply(x, parse_one, numeric(1))
    as.POSIXct(seconds, origin = "1970-01-01", tz = "UTC")
  }

  normalized_id <- normalize_metadata_string(passaging_tbl[[id_col]])
  parsed <- stringr::str_match(
    normalized_id,
    "^SUM_?159(?:_NLS)?_(2N|4N)_(C|O1|O2)_A([0-9]+)(?:_|$)"
  )

  dplyr::bind_cols(
    passaging_tbl,
    tibble::tibble(
      normalized_id = normalized_id,
      signature_ok = !is.na(parsed[, 1]),
      cell_line_key = ifelse(!is.na(parsed[, 1]), "SUM159", NA_character_),
      ploidy_key = parsed[, 2],
      oxygen_key = parsed[, 3],
      address_key = ifelse(!is.na(parsed[, 4]), paste0("A", parsed[, 4]), NA_character_),
      match_time = parse_match_time(passaging_tbl[[time_col]])
    )
  )
}

match_metadata_signatures <- function(
    sample_tbl,
    passaging_tbl,
    sample_col = "sample_name",
    passaging_id_col = "id",
    passaging_time_col = "date",
    sample_parser = parse_sum159_hypoxia_sample_signature,
    passaging_parser = parse_sum159_passaging_signature) {
  stopifnot(sample_col %in% names(sample_tbl))

  sample_sig <- sample_parser(sample_tbl[[sample_col]])
  passaging_sig <- passaging_parser(
    passaging_tbl = passaging_tbl,
    id_col = passaging_id_col,
    time_col = passaging_time_col
  )
  passaging_sig <- dplyr::filter(passaging_sig, signature_ok)

  sample_sig <- dplyr::filter(sample_sig, signature_ok)
  joined <- dplyr::left_join(
    sample_sig,
    passaging_sig,
    by = c("cell_line_key", "ploidy_key", "oxygen_key", "address_key"),
    relationship = "many-to-many"
  )

  dplyr::select(
    dplyr::rename(
      joined,
      matched_id = !!passaging_id_col,
      matched_time = match_time
    ),
    sample_name,
    sample_id,
    cell_line_key,
    ploidy_key,
    oxygen_key,
    address_key,
    matched_id,
    matched_time,
    dplyr::everything()
  )
}

summarize_latest_metadata_matches <- function(match_tbl) {
  first_non_missing <- function(x) {
    idx <- which(!is.na(x))[1]
    if (is.na(idx)) {
      return(x[NA_integer_])
    }
    x[idx]
  }

  match_tbl <- dplyr::distinct(
    match_tbl,
    sample_name,
    sample_id,
    cell_line_key,
    ploidy_key,
    oxygen_key,
    address_key,
    matched_id,
    matched_time,
    .keep_all = TRUE
  )

  out <- dplyr::group_by(
    match_tbl,
    sample_name,
    sample_id,
    cell_line_key,
    ploidy_key,
    oxygen_key,
    address_key
  )
  out <- dplyr::arrange(out, dplyr::desc(matched_time), matched_id, .by_group = TRUE)
  out <- dplyr::summarise(
    out,
    match_count = sum(!is.na(matched_id)),
    matched_ids = paste(stats::na.omit(matched_id), collapse = " | "),
    latest_matching_id = first_non_missing(matched_id),
    latest_match_time = first_non_missing(matched_time),
    .groups = "drop"
  )
  dplyr::mutate(
    out,
    latest_match_date = as.Date(latest_match_time),
    relative_day = dplyr::if_else(
      is.na(latest_match_date),
      NA_integer_,
      as.integer(latest_match_date - min(latest_match_date, na.rm = TRUE))
    )
  )
}

build_sum159_hypoxia_match_table <- function(sample_names, passaging_tbl) {
  sample_tbl <- tibble::tibble(sample_name = basename(sample_names))
  match_tbl <- match_metadata_signatures(sample_tbl, passaging_tbl = passaging_tbl)
  summary_tbl <- summarize_latest_metadata_matches(match_tbl)
  summary_tbl <- dplyr::distinct(summary_tbl, sample_name, .keep_all = TRUE)
  dplyr::arrange(summary_tbl, relative_day, latest_match_date, sample_name)
}
