source(here::here("functions", "network_similarity_measures.R"), local = TRUE)

assoc_to_graph <- function(assoc, thr) {
  adj <- assoc
  
  diag(adj) <- 0
  adj[abs(adj) <= thr] <- 0
  adj <- (adj + t(adj)) / 2
  
  igraph::graph_from_adjacency_matrix(
    adjmatrix = adj,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
}

# Network threshold diagnosis ------------------------------------------------

similarity_from_graphs <- function(graphs, thr, method = "pearson",
                                   cluster_fun = netkit::cluster_signed_louvain) {
  graph_names <- names(graphs)
  
  pairs <- utils::combn(graph_names, 2, simplify = FALSE)
  
  purrr::map_dfr(
    pairs,
    \(pair) {
      g1 <- graphs[[pair[1]]]
      g2 <- graphs[[pair[2]]]
      
      tibble::tibble(
        threshold = thr,
        network1 = pair[1],
        network2 = pair[2],
        edge_jaccard_signed = edge_jaccard_signed(g1, g2),
        degree_cor = degree_correlation(g1, g2, method = method),
        signed_strength_cor = signed_strength_correlation(g1, g2, method = method),
        community_ari = community_ari_from_graphs(g1 = g1, 
                                                  g2 = g2,
                                                  cluster_fun = cluster_fun))
    })
}

diagnosis_from_graph <- function(g, thr, network_name = NA_character_) {
  deg <- igraph::degree(g)
  comp <- igraph::components(g)
  comp_sizes <- sort(comp$csize, decreasing = TRUE)
  
  n_nodes <- igraph::gorder(g)
  n_edges <- igraph::gsize(g)
  
  tibble::tibble(
    network = network_name,
    threshold = thr,
    
    n_nodes = n_nodes,
    n_edges = n_edges,
    density = igraph::edge_density(g, loops = FALSE),
    avg_degree = mean(deg),
    frac_isolated = mean(deg == 0),
    
    n_components = comp$no,
    n_components_gt5 = sum(comp$csize > 5),
    lcc_size = max(comp$csize),
    lcc_fraction = max(comp$csize) / n_nodes,
    second_lcc_size = if (length(comp_sizes) >= 2) comp_sizes[2] else 0L,
    second_lcc_fraction = if (length(comp_sizes) >= 2) comp_sizes[2] / n_nodes else 0,
    component_entropy = vegan::diversity(comp$csize),
    
    n_positive_edges = if (n_edges > 0) sum(igraph::E(g)$weight > 0) else 0L,
    n_negative_edges = if (n_edges > 0) sum(igraph::E(g)$weight < 0) else 0L,
    frac_positive = if (n_edges > 0) mean(igraph::E(g)$weight > 0) else NA_real_,
    frac_negative = if (n_edges > 0) mean(igraph::E(g)$weight < 0) else NA_real_,
    
    components_obj = list(comp)
  )
}


threshold_network_analysis <- function(
    assoc_list,
    thresholds,
    method = "pearson",
    cluster_fun = netkit::cluster_signed_louvain,
    workers = max(2, parallel::detectCores() - 1),
    seed = 42
) {
  
  if (workers <= 1) {
    stop("`workers` must be greater than 1.")
  }
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  future::plan(
    future::multisession,
    workers = workers
  )
  
  results <- furrr::future_map(
    thresholds,
    \(thr) {
      graphs <- purrr::imap(
        assoc_list,
        \(assoc, network_name) assoc_to_graph(assoc, thr)
      )
      
      diagnosis_tbl <- purrr::imap_dfr(
        graphs,
        \(g, network_name) diagnosis_from_graph(
          g = g,
          thr = thr,
          network_name = network_name
        )
      )
      
      similarity_tbl <- similarity_from_graphs(
        graphs = graphs,
        thr = thr,
        method = method,
        cluster_fun = cluster_fun
      )
      
      list(
        diagnosis = diagnosis_tbl,
        similarity = similarity_tbl
      )
    },
    .options = furrr::furrr_options(seed = seed),
    .progress = TRUE
  )
  
  list(
    diagnosis = purrr::map(results, \(x) x$diagnosis) %>% purrr::list_rbind(),
    similarity = purrr::map(results, \(x) x$similarity) %>% purrr::list_rbind()
  )
}