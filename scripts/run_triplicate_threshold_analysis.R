# Triplicate threshold-dependent network analysis -----------------------------

suppressPackageStartupMessages({
  library(mgnet)
  library(igraph)
  library(netkit)
  
  library(dplyr)
  library(purrr)
  library(tibble)
  
  library(future)
  library(furrr)
  library(here)
})


# Helpers --------------------------------------------------------------------

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  key <- paste0("--", name)
  pos <- match(key, args)
  
  if (is.na(pos)) {
    return(default)
  }
  
  if (pos == length(args)) {
    stop("Missing value for argument: ", key)
  }
  
  args[[pos + 1]]
}


as_number <- function(x, name) {
  value <- suppressWarnings(as.numeric(x))
  
  if (is.na(value)) {
    stop("Argument `", name, "` must be numeric.")
  }
  
  value
}


as_integer <- function(x, name) {
  value <- suppressWarnings(as.integer(x))
  
  if (is.na(value)) {
    stop("Argument `", name, "` must be an integer.")
  }
  
  value
}


format_param <- function(x) {
  gsub("\\.", "p", as.character(x))
}


make_output_name <- function(prevalence, min_count, seed, method) {
  paste0(
    "triplicate",
    "_prev", format_param(prevalence),
    "_mincount", format_param(min_count),
    "_seed", seed,
    "_", method,
    ".rds"
  )
}


# Main function ---------------------------------------------------------------

run_triplicate_threshold_analysis <- function(prevalence = 0.25,
                                              min_count = 10,
                                              threshold_min = 0.50,
                                              threshold_max = 0.95,
                                              threshold_by = 0.01,
                                              method = "pearson",
                                              workers = max(2, parallel::detectCores() - 1),
                                              seed = 42) {
  input_rds <- here::here("data", "processed", "sewage_cph_mgnet.rds")
  functions_dir <- here::here("functions")
  
  output_name <- make_output_name(
    prevalence = prevalence,
    min_count = min_count,
    seed = seed,
    method = method
  )
  
  output_rds <- here::here("cache", output_name)
  
  if (!file.exists(input_rds)) {
    stop("Input file not found: ", input_rds)
  }
  
  if (!dir.exists(functions_dir)) {
    stop("Functions directory not found: ", functions_dir)
  }
  
  threshold_fun_file <- file.path(functions_dir, "network_threshold_diagnosis.R")
  
  if (!file.exists(threshold_fun_file)) {
    stop("Function file not found: ", threshold_fun_file)
  }
  
  if (!method %in% c("pearson", "spearman", "kendall")) {
    stop("`method` must be one of: pearson, spearman, kendall.")
  }
  
  if (workers <= 1) {
    stop("`workers` must be greater than 1.")
  }
  
  output_dir <- dirname(output_rds)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  netfun <- new.env(parent = globalenv())
  
  source(
    threshold_fun_file,
    local = netfun
  )
  
  message("Reading input object: ", input_rds)
  
  plants_cph <- readRDS(input_rds) |>
    split_meta(plant)
  
  message("Filtering taxa...")
  message("Prevalence threshold: ", prevalence)
  message("Minimum count: ", min_count)
  
  plants_cph_filt <- plants_cph |>
    update_mgnet(rela = abun / rowSums(abun)) |>
    mutate_meta(sample_sum = sum(abun)) |>
    filter_taxa(sum(abun >= min_count) / n() >= prevalence) |>
    mutate_meta(lost_counts = 1 - (sum(abun) / sample_sum))
  
  common_taxa <- taxa_id(plants_cph_filt) |>
    purrr::reduce(intersect)
  
  message("Number of common taxa after filtering: ", length(common_taxa))
  
  plants_cph_common <- plants_cph_filt |>
    mutate_meta(sample_sum_filt = sum(abun)) |>
    filter_taxa(taxa_id %in% common_taxa) |>
    mutate_meta(lost_counts_common = 1 - (sum(abun) / sample_sum_filt))
  
  message("CLR normalization...")
  
  plants_cph_common <- plants_cph_common |>
    update_mgnet(
      norm = abun |>
        zCompositions::cmultRepl(
          method = "CZM",
          z.delete = FALSE,
          z.warning = 1
        ) |>
        vegan::decostand(method = "clr")
    )
  
  message("Computing association matrices...")
  
  assoc_list <- norm(plants_cph_common) |>
    purrr::map(\(x) stats::cor(
      x,
      method = method,
      use = "pairwise.complete.obs"
    ))
  
  thresholds <- seq(
    from = threshold_min,
    to = threshold_max,
    by = threshold_by
  )
  
  message("Running threshold-dependent analysis...")
  message("Threshold range: ", threshold_min, " - ", threshold_max, " by ", threshold_by)
  message("Number of thresholds: ", length(thresholds))
  message("Correlation method: ", method)
  message("Workers: ", workers)
  message("Seed: ", seed)
  
  threshold_results <- netfun$threshold_network_analysis(
    assoc_list = assoc_list,
    thresholds = thresholds,
    method = method,
    cluster_fun = netkit::cluster_signed_louvain,
    workers = workers,
    seed = seed
  )
  
  output <- list(
    diagnosis = threshold_results$diagnosis,
    similarity = threshold_results$similarity
  )
  
  saveRDS(output, output_rds)
  
  message("Saved results to: ", output_rds)
  message("Done.")
  
  invisible(output)
}


# Command-line arguments ------------------------------------------------------

prevalence <- as_number(
  get_arg("prevalence", "0.25"),
  "prevalence"
)

min_count <- as_number(
  get_arg("min-count", "10"),
  "min-count"
)

threshold_min <- as_number(
  get_arg("threshold-min", "0.50"),
  "threshold-min"
)

threshold_max <- as_number(
  get_arg("threshold-max", "0.95"),
  "threshold-max"
)

threshold_by <- as_number(
  get_arg("threshold-by", "0.01"),
  "threshold-by"
)

method <- get_arg(
  name = "method",
  default = "pearson"
)

workers <- as_integer(
  get_arg(
    name = "workers",
    default = as.character(max(2, parallel::detectCores() - 1))
  ),
  "workers"
)

seed <- as_integer(
  get_arg("seed", "42"),
  "seed"
)


# Run -------------------------------------------------------------------------

run_triplicate_threshold_analysis(
  prevalence = prevalence,
  min_count = min_count,
  threshold_min = threshold_min,
  threshold_max = threshold_max,
  threshold_by = threshold_by,
  method = method,
  workers = workers,
  seed = seed
)