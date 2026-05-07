# Network similarity measures -------------------------------------------------

edge_ids_signed <- function(g) {
  ed <- igraph::as_edgelist(g, names = TRUE)
  
  if (nrow(ed) == 0) {
    return(character(0))
  }
  
  w <- igraph::E(g)$weight
  
  tibble::tibble(
    from = pmin(ed[, 1], ed[, 2]),
    to = pmax(ed[, 1], ed[, 2]),
    sign = sign(w)
  ) |>
    dplyr::mutate(edge_id = paste(from, to, sign, sep = "||")) |>
    dplyr::pull(edge_id)
}


edge_jaccard_signed <- function(g1, g2) {
  e1 <- edge_ids_signed(g1)
  e2 <- edge_ids_signed(g2)
  
  edge_union <- union(e1, e2)
  
  if (length(edge_union) == 0) {
    return(NA_real_)
  }
  
  length(intersect(e1, e2)) / length(edge_union)
}


degree_correlation <- function(g1, g2, method = "pearson") {
  nodes <- union(igraph::V(g1)$name, igraph::V(g2)$name)
  
  d1 <- stats::setNames(rep(0, length(nodes)), nodes)
  d2 <- stats::setNames(rep(0, length(nodes)), nodes)
  
  d1[igraph::V(g1)$name] <- igraph::degree(g1)
  d2[igraph::V(g2)$name] <- igraph::degree(g2)
  
  if (stats::sd(d1) == 0 || stats::sd(d2) == 0) {
    return(NA_real_)
  }
  
  stats::cor(d1, d2, method = method)
}


signed_strength_correlation <- function(g1, g2, method = "pearson") {
  nodes <- union(igraph::V(g1)$name, igraph::V(g2)$name)
  
  s1 <- stats::setNames(rep(0, length(nodes)), nodes)
  s2 <- stats::setNames(rep(0, length(nodes)), nodes)
  
  s1[igraph::V(g1)$name] <- igraph::strength(g1)
  s2[igraph::V(g2)$name] <- igraph::strength(g2)
  
  if (stats::sd(s1) == 0 || stats::sd(s2) == 0) {
    return(NA_real_)
  }
  
  stats::cor(s1, s2, method = method)
}

community_ari_from_graphs <- function(g1, g2, cluster_fun = netkit::cluster_signed_louvain) {
  
  comm1 <- tryCatch(
    cluster_fun(g1),
    error = function(e) {
      warning("Community detection failed for g1: ", conditionMessage(e))
      NULL
    }
  )

  comm2 <- tryCatch(
    cluster_fun(g2),
    error = function(e) {
      warning("Community detection failed for g2: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(comm1) || is.null(comm2)) {
    return(NA_real_)
  }

  memb1 <- igraph::membership(comm1)
  memb2 <- igraph::membership(comm2)

  names(memb1) <- igraph::V(g1)$name
  names(memb2) <- igraph::V(g2)$name

  common_nodes <- intersect(names(memb1), names(memb2))

  if (length(common_nodes) == 0) {
    return(NA_real_)
  }

  mclust::adjustedRandIndex(
    memb1[common_nodes],
    memb2[common_nodes]
  )
}
