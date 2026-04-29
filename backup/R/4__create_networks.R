library(psych)
library(readr)
library(igraph)

copenhagen_abundances <- readRDS("..\\data\\processed\\abun_plants.rds")

common_nodes <- read.csv("..\\data\\processed\\sp_name_idx.csv")
common_nodes$X <- NULL
common_nodes <- common_nodes$sp_name

#-------------------------------------------------------------------------------

abundances_RA <- copenhagen_abundances$RA[ , names(copenhagen_abundances$RA) %in% common_nodes]
abundances_RD <- copenhagen_abundances$RD[ , names(copenhagen_abundances$RD) %in% common_nodes]
abundances_RL <- copenhagen_abundances$RL[ , names(copenhagen_abundances$RL) %in% common_nodes]

# Zeroes correction
correct_dataset <- function(site_abundances){
    
    correct_counts <- site_abundances
    for(sample_name in rownames(site_abundances)){
        
        sample <- site_abundances[sample_name, ]
        
        K  = length(sample)
        Z = length(sample[sample==0])
        
        delta <- 1/(K^2)
        
        nonzero_correction <- 1 - Z*delta/sum(sample)
        
        sample[sample==0] <- delta
        sample[sample!=0] <- nonzero_correction*sample[sample!=0]
        
        correct_counts[sample_name, ] <- sample        
        
    }
    
    return(correct_counts)  
    
}

RA_correct <- correct_dataset(abundances_RA)
RD_correct <- correct_dataset(abundances_RD)
RL_correct <- correct_dataset(abundances_RL)

# Central Log Ratio
apply_CLR_to_dataset <- function(counts){
    
    # MARGIN = 1 for rows, = 2 for columns
    ref <- apply(counts, 1, function(x) mean(log(x)))
    transformed <- sweep(log(counts), 1, ref)
    return(transformed)
    
}

RA_clr <- apply_CLR_to_dataset(RA_correct)
RD_clr <- apply_CLR_to_dataset(RD_correct)
RL_clr <- apply_CLR_to_dataset(RL_correct)

#-------------------------------------------------------------------------------

# Correlations and p-values
RA_corr_test <- corr.test(
    x=RA_clr, 
    use="pairwise", method="pearson",
    adjust="bonferroni",
    ci=FALSE
)

RD_corr_test <- corr.test(
    x=RD_clr, 
    use="pairwise", method="pearson",
    adjust="bonferroni",
    ci=FALSE
)

RL_corr_test <- corr.test(
    x=RL_clr, 
    use="pairwise", method="pearson",
    adjust="bonferroni",
    ci=FALSE
)

RA_corr <- RA_corr_test$r
RD_corr <- RD_corr_test$r
RL_corr <- RL_corr_test$r
diag(RA_corr) <- 0
diag(RD_corr) <- 0
diag(RL_corr) <- 0

names(RA_corr) <- common_nodes
names(RD_corr) <- common_nodes
names(RL_corr) <- common_nodes

RA_pval <- RA_corr_test$p
RD_pval <- RD_corr_test$p
RL_pval <- RL_corr_test$p

# p-values returned by the function are weird
#     Upper triangular are the corrected
#     Lower triangular are the not corrected
RA_pval[lower.tri(RA_pval)] <- t(RA_pval)[lower.tri(RA_pval)]
RD_pval[lower.tri(RD_pval)] <- t(RD_pval)[lower.tri(RD_pval)]
RL_pval[lower.tri(RL_pval)] <- t(RL_pval)[lower.tri(RL_pval)]
diag(RA_pval) <- 1
diag(RD_pval) <- 1
diag(RL_pval) <- 1

# Check if 
#isSymmetric(RA_pval)

# Function which returns the correlation value which leaves out a given 
#     fraction of nodes
density_to_correlation_threshold <- function(site_corr, site_pval, density_thr){
    
    N_links <- length(site_corr)
    index_thr <- 1 + as.integer(density_thr*N_links)
    
    # Set to 0 not significative corrs, so they are not chosen
    significative_corr <- site_corr*apply(site_pval<0.05, c(1,2), as.integer)
    
    corr_thr <- abs(
        site_corr[
            order(abs(
                significative_corr
            ), decreasing=TRUE)
        ][[index_thr]]
    )
    
    return(corr_thr)
    
}

density_to_correlation_threshold(RA_corr, RA_pval, 5.6/100) # 0.77
density_to_correlation_threshold(RD_corr, RD_pval, 5.6/100) # 0.82
density_to_correlation_threshold(RL_corr, RL_pval, 5.6/100) # 0.79

# Create adjacency matrices for all sites with specified density
density_to_adjacencies <- function(density_thr){
    
    RA_r_thr <- density_to_correlation_threshold(RA_corr, RA_pval, density_thr)
    RD_r_thr <- density_to_correlation_threshold(RD_corr, RD_pval, density_thr)
    RL_r_thr <- density_to_correlation_threshold(RL_corr, RL_pval, density_thr)

    RA_adj <- RA_corr
    RA_adj[(abs(RA_corr) < RA_r_thr)] <- 0
    RA_adj[RA_pval > 0.05] <- 0
    RA_adj[RA_adj > 0] <-  1
    RA_adj[RA_adj < 0] <- -1
    
    RD_adj <- RD_corr
    RD_adj[(abs(RD_corr) < RD_r_thr)] <- 0
    RD_adj[RD_pval > 0.05] <- 0
    RD_adj[RD_adj > 0] <-  1
    RD_adj[RD_adj < 0] <- -1
    
    RL_adj <- RL_corr
    RL_adj[(abs(RL_corr) < RL_r_thr)] <- 0
    RL_adj[RL_pval > 0.05] <- 0
    RL_adj[RL_adj > 0] <-  1
    RL_adj[RL_adj < 0] <- -1
    
    adjacencies <- list("RA" = RA_adj, "RD" = RD_adj, "RL" = RL_adj)

    return(adjacencies)
    
}

# This value has been selected through the table computed below
adjacencies <- density_to_adjacencies(5.6/100)
write_rds(adjacencies, file="..\\data\\processed\\adjacencies.rds")

#-------------------------------------------------------------------------------
# Density search
#   Computation of the table used for choosing the density threshold
density_to_components <- function(density_thr){
    
    adjacencies <- density_to_adjacencies(density_thr)
    
    g_RA <- graph_from_adjacency_matrix(
        adjacencies$RA,
        mode = "undirected",
        weighted = TRUE
    )
    g_RD <- graph_from_adjacency_matrix(
        adjacencies$RD,
        mode = "undirected",
        weighted = TRUE
    )
    g_RL <- graph_from_adjacency_matrix(
        adjacencies$RL,
        mode = "undirected",
        weighted = TRUE
    )
    
    comp_RA <- components(g_RA)
    comp_RD <- components(g_RD)
    comp_RL <- components(g_RL)
    
    first_comp_RA <- sort(comp_RA$csize, decreasing=TRUE)[[1]]
    first_comp_RD <- sort(comp_RD$csize, decreasing=TRUE)[[1]]
    first_comp_RL <- sort(comp_RL$csize, decreasing=TRUE)[[1]]
    
    second_comp_RA <- sort(comp_RA$csize, decreasing=TRUE)[[2]]
    second_comp_RD <- sort(comp_RD$csize, decreasing=TRUE)[[2]]
    second_comp_RL <- sort(comp_RL$csize, decreasing=TRUE)[[2]]
    
    return(
        list(
            "first_comp" = list(
                "RA"=first_comp_RA, "RD"=first_comp_RD, "RL"=first_comp_RL
            ),
            "second_comp" = list(
                "RA"=second_comp_RA, "RD"=second_comp_RD, "RL"=second_comp_RL
            )
        )
    )
    
}

density_thresholds <- c(
    1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 5.6, 5.7, 5.8, 6.0, 7.0, 8.0, 9.0
)
first_comp <- list("RA"=c(), "RD"=c(), "RL"=c())
second_comp <- list("RA"=c(), "RD"=c(), "RL"=c())

log_i <- 1
for(d in density_thresholds){
    message(log_i, "/", length(density_thresholds))
    log_i <- log_i+1
    results <- density_to_components(d/100)
    first_comp$RA <- c(first_comp$RA, results$first_comp$RA)
    first_comp$RD <- c(first_comp$RD, results$first_comp$RD)
    first_comp$RL <- c(first_comp$RL, results$first_comp$RL)
    second_comp$RA <- c(second_comp$RA, results$second_comp$RA)
    second_comp$RD <- c(second_comp$RD, results$second_comp$RD)
    second_comp$RL <- c(second_comp$RL, results$second_comp$RL)
}

write_rds(
    list("first_comp"=first_comp, "second_comp"=second_comp, "densities"=density_thresholds),
    file="..\\data\\processed\\density_search.rds"
)
