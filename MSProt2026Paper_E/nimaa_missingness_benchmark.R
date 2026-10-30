#!/usr/bin/env Rscript

# NIMAA missingness benchmark for proteomics


required <- "NIMAA"
if (!requireNamespace(required, quietly = TRUE)) {
  stop(
    "Install NIMAA 0.2.2 before running, e.g. install.packages('NIMAA').",
    call. = FALSE
  )
}
if (utils::packageVersion("NIMAA") != "0.2.2") {
  warning(
    "This benchmark was recorded with NIMAA 0.2.2; installed version is ",
    utils::packageVersion("NIMAA"), ".",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(root, "nimaa_missingness_results")
root <- normalizePath(root, mustWork = TRUE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)
dir.create(file.path(out_dir, "rectangles"), showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Required input was not found: ", path, call. = FALSE)
  path
}

as_numeric_matrix <- function(x, label) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(label, " must be a numeric matrix (or numeric data frame).", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(label, " needs feature row names and run column names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(x))) stop(label, " has duplicated feature identifiers.", call. = FALSE)
  x
}

# RDS files may store one matrix directly, a named list of matrices, or a named
# list one level deeper.  Extract only numeric, named matrices.
read_matrix_set <- function(path, default_name) {
  object <- readRDS(stop_if_missing(path))
  collect <- function(x, prefix = "") {
    if (is.matrix(x) || is.data.frame(x)) {
      name <- if (nzchar(prefix)) prefix else default_name
      return(stats::setNames(list(as_numeric_matrix(x, name)), name))
    }
    if (!is.list(x)) return(list())
    nms <- names(x)
    if (is.null(nms)) nms <- as.character(seq_along(x))
    pieces <- Map(
      function(item, name) {
        next_prefix <- if (nzchar(prefix)) paste(prefix, name, sep = "_") else name
        collect(item, next_prefix)
      },
      x, nms
    )
    unlist(pieces, recursive = FALSE)
  }
  result <- collect(object)
  if (!length(result)) stop("No numeric matrices found in ", path, call. = FALSE)
  result
}

normalise_run_id <- function(x) {
  x <- basename(trimws(as.character(x)))
  sub("\\.raw$", "", x, ignore.case = TRUE)
}

find_column <- function(x, candidates, purpose) {
  hit <- candidates[tolower(candidates) %in% tolower(names(x))]
  if (!length(hit)) {
    stop("Could not find a ", purpose, " column. Expected one of: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  names(x)[match(tolower(hit[[1L]]), tolower(names(x)))]
}

read_metadata <- function(path) {
  x <- utils::read.csv(stop_if_missing(path), check.names = FALSE, stringsAsFactors = FALSE)
  run_col <- find_column(x, c("run_id", "run", "raw_file", "rawfile", "file_name"), "run ID")
  sample_col <- find_column(
    x,
    c("biological_sample", "biological_sample_id", "sample_id", "sample"),
    "biological-sample"
  )
  x$run_id <- normalise_run_id(x[[run_col]])
  x$biological_sample <- as.character(x[[sample_col]])
  if (anyDuplicated(x$run_id)) stop("Metadata has duplicated run IDs.", call. = FALSE)
  x
}

align_to_runs <- function(mat, run_ids, label) {
  run_names <- normalise_run_id(colnames(mat))
  if (anyDuplicated(run_names)) stop(label, " has duplicated run IDs after harmonisation.", call. = FALSE)
  colnames(mat) <- run_names
  missing <- setdiff(run_ids, run_names)
  extra <- setdiff(run_names, run_ids)
  if (length(missing) || length(extra)) {
    stop(
      label, " does not match the harmonised run set. Missing: ",
      paste(missing, collapse = ", "), "; extra: ", paste(extra, collapse = ", "),
      call. = FALSE
    )
  }
  mat[, run_ids, drop = FALSE]
}

to_detection <- function(x) {
  present <- is.finite(x) & x > 0
  out <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
  out[present] <- 1
  out
}

write_tsv <- function(x, filename) {
  utils::write.table(x, filename, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

matrix_summary <- function(det, level, workflow, calculate_temperature = TRUE) {
  observed <- !is.na(det)
  temperature <- if (calculate_temperature) {
    unname(bipartite::nested(observed * 1))
  } else {
    NA_real_
  }
  data.frame(
    level = level,
    workflow = workflow,
    features = nrow(det),
    runs = ncol(det),
    observed_cells = sum(observed),
    observed_cell_fraction = mean(observed),
    strict_complete_features = sum(rowSums(observed) == ncol(det)),
    nestedness_temperature = temperature,
    stringsAsFactors = FALSE
  )
}

run_summary <- function(det, level, workflow) {
  data.frame(
    level = level,
    workflow = workflow,
    run_id = colnames(det),
    observed_features = colSums(!is.na(det)),
    observed_feature_fraction = colMeans(!is.na(det)),
    stringsAsFactors = FALSE
  )
}

common_identifier_summary <- function(detection_set, level) {
  shared_ids <- Reduce(intersect, lapply(detection_set, rownames))
  if (!length(shared_ids)) stop("No identifiers are shared across ", level, " matrices.", call. = FALSE)
  do.call(rbind, lapply(names(detection_set), function(workflow) {
    det <- detection_set[[workflow]][shared_ids, , drop = FALSE]
    cbind(
      data.frame(level = level, workflow = workflow, common_features = length(shared_ids)),
      matrix_summary(det, level, workflow, calculate_temperature = FALSE)[
        , c("observed_cells", "observed_cell_fraction", "strict_complete_features")
      ]
    )
  }))
}

technical_consensus <- function(det, metadata, level, workflow) {
  groups <- split(seq_len(ncol(det)), metadata$biological_sample)
  if (any(lengths(groups) != 3L)) {
    warning(level, "/", workflow, ": biological samples are not all triplicates; ",
            "the 1/2/3-of-replicate results use the available runs.", call. = FALSE)
  }
  observed <- !is.na(det)
  do.call(rbind, lapply(1:3, function(k) {
    retained <- Reduce("&", lapply(groups, function(i) rowSums(observed[, i, drop = FALSE]) >= k))
    data.frame(
      level = level, workflow = workflow, replicate_rule = paste0(k, "_of_3"),
      retained_features = sum(retained), stringsAsFactors = FALSE
    )
  }))
}

extract_nimaa_no_plot <- function(x, bar = 1) {
  # NIMAA::extractSubMatrix additionally materialises two all-cell plotly
  # widgets. Reuse its internal arranging and rectangle-selection functions,
  # but omit that display-only step for scalable peptide analysis.
  arrange <- getFromNamespace("threeStepArrange", "NIMAA")
  element_max <- getFromNamespace("findElementMaxSubmatrix", "NIMAA")
  row_cutoff <- getFromNamespace("findRowwiseCutoffSubmatrix", "NIMAA")
  col_cutoff <- getFromNamespace("findColwiseCutoffSubmatrix", "NIMAA")

  x <- as.data.frame(x)
  x_arranged <- arrange(x)
  x_transposed <- as.data.frame(t(x))
  x_transposed_arranged <- arrange(x_transposed)
  choose_larger <- function(a, b) if (prod(dim(a)) < prod(dim(b))) b else a

  row_candidate <- choose_larger(
    row_cutoff(x_arranged, bar),
    t(col_cutoff(x_transposed_arranged, bar))
  )
  col_candidate <- choose_larger(
    col_cutoff(x_arranged, bar),
    t(row_cutoff(x_transposed_arranged, bar))
  )
  element_candidate <- choose_larger(
    element_max(x_arranged),
    t(element_max(x_transposed_arranged))
  )
  restore <- function(candidate) x[rownames(candidate), colnames(candidate), drop = FALSE]
  list(
    rectangles = list(
      Rectangular_element_max = restore(element_candidate),
      Rectangular_row = restore(row_candidate),
      Rectangular_col = restore(col_candidate)
    ),
    feature_order = rownames(x_arranged),
    run_order = colnames(x_arranged)
  )
}

run_nimaa <- function(det, level, workflow) {
  shapes <- c("Rectangular_element_max", "Rectangular_row", "Rectangular_col")
  extracted <- extract_nimaa_no_plot(det, bar = 1)
  selected <- extracted$rectangles
  map_stem <- paste(level, safe_name(workflow), "NIMAA_arrangement", sep = "_")
  write_tsv(
    data.frame(feature_id = extracted$feature_order, arranged_rank = seq_along(extracted$feature_order)),
    file.path(out_dir, "rectangles", paste0(map_stem, "_features.tsv"))
  )
  write_tsv(
    data.frame(run_id = extracted$run_order, arranged_rank = seq_along(extracted$run_order)),
    file.path(out_dir, "rectangles", paste0(map_stem, "_runs.tsv"))
  )
  do.call(rbind, lapply(shapes, function(shape) {
    rect <- selected[[shape]]
    if (is.null(rect)) stop("NIMAA did not return ", shape, " for ", workflow, call. = FALSE)
    if (anyNA(rect)) stop("NIMAA returned a non-complete ", shape, " rectangle.", call. = FALSE)
    stem <- paste(level, safe_name(workflow), safe_name(shape), sep = "_")
    write_tsv(data.frame(feature_id = rownames(rect)), file.path(out_dir, "rectangles", paste0(stem, "_features.tsv")))
    write_tsv(data.frame(run_id = colnames(rect)), file.path(out_dir, "rectangles", paste0(stem, "_runs.tsv")))
    data.frame(
      level = level, workflow = workflow, shape = shape,
      selected_features = nrow(rect), selected_runs = ncol(rect),
      selected_cell_fraction = prod(dim(rect)) / prod(dim(det)),
      stringsAsFactors = FALSE
    )
  }))
}

protein_rds <- file.path(root, "NIMAA_benchmark", "matrices", "primary_matrices.rds")
peptide_dir <- file.path(root, "NIMAA_peptide_benchmark", "matrices")
metadata_csv <- file.path(root, "NIMAA_benchmark", "results", "harmonized_samples.csv")

protein <- read_matrix_set(protein_rds, "protein")
peptide_files <- list.files(peptide_dir, pattern = "_primary\\.rds$", full.names = TRUE, ignore.case = TRUE)
if (!length(peptide_files)) stop("No *_primary.rds peptide matrices found in ", peptide_dir, call. = FALSE)
peptide <- do.call(c, unname(lapply(peptide_files, function(path) {
  workflow <- sub("_primary\\.rds$", "", basename(path), ignore.case = TRUE)
  result <- read_matrix_set(path, workflow)
  if (length(result) == 1L) names(result) <- workflow
  result
})))
if (anyDuplicated(names(peptide))) stop("Peptide workflow names are not unique.", call. = FALSE)

metadata <- read_metadata(metadata_csv)
run_ids <- metadata$run_id
process_level <- function(matrix_set, level) {
  matrix_names <- names(matrix_set)
  matrix_set <- lapply(matrix_names, function(name) {
    to_detection(align_to_runs(matrix_set[[name]], run_ids, paste(level, name)))
  })
  names(matrix_set) <- matrix_names

  summaries <- do.call(rbind, Map(function(x, name) matrix_summary(x, level, name), matrix_set, names(matrix_set)))
  runs <- do.call(rbind, Map(function(x, name) run_summary(x, level, name), matrix_set, names(matrix_set)))
  shared <- common_identifier_summary(matrix_set, level)
  consensus <- do.call(rbind, Map(
    function(x, name) technical_consensus(x, metadata, level, name), matrix_set, names(matrix_set)
  ))
  rectangles <- do.call(rbind, Map(function(x, name) run_nimaa(x, level, name), matrix_set, names(matrix_set)))

  saveRDS(matrix_set, file.path(out_dir, paste0(level, "_detection_1_NA.rds")))
  list(summary = summaries, run = runs, common = shared, consensus = consensus, rectangle = rectangles)
}

protein_results <- process_level(protein, "protein")
peptide_results <- process_level(peptide, "peptide")

write_tsv(rbind(protein_results$summary, peptide_results$summary), file.path(out_dir, "matrix_summary.tsv"))
write_tsv(rbind(protein_results$run, peptide_results$run), file.path(out_dir, "run_missingness_map.tsv"))
write_tsv(rbind(protein_results$common, peptide_results$common), file.path(out_dir, "common_identifier_summary.tsv"))
write_tsv(rbind(protein_results$consensus, peptide_results$consensus), file.path(out_dir, "technical_replicate_consensus.tsv"))
write_tsv(rbind(protein_results$rectangle, peptide_results$rectangle), file.path(out_dir, "nimaa_rectangle_summary.tsv"))

capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
message("Finished. Results written to: ", out_dir)
