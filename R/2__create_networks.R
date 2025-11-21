library(psych)

copenhagen_abundances <- readRDS("..\\data\\processed\\abun_plants.rds")
metadata <- read.csv("..\\data\\raw\\cph_sewage_2020_metadata.csv", sep=';')

#-------------------------------------------------------------------------------

get_species_stats <- function(df, site){

    nonzero_medians <- c()
    prevalence_and_numerosity <- c()

    species = names(df)
    
    for(bact in species){
        
        bact_abundances <-  df[ ,bact]
        
        prevalence_and_numerosity <- c(
            prevalence_and_numerosity, 
            length(bact_abundances[bact_abundances >= 10])/length(bact_abundances)
        )

        nonzero_filter <- bact_abundances!=0
        if(length(bact_abundances[nonzero_filter])!=0){
            nonzero_medians <- c(nonzero_medians, median(bact_abundances[nonzero_filter]))
        }else{
            nonzero_medians <- c(nonzero_medians, 0)
        }
        
    }
    
    return(
        data.frame(
            "name" = unlist(species),
            "nonzero_medians" = nonzero_medians,
            "prevalence_and_numerosity" = prevalence_and_numerosity,
            "site" = site
        )
    )
    
}

species_stats <- rbind(
    get_species_stats(copenhagen_abundances$RA, "RA"),
    get_species_stats(copenhagen_abundances$RD, "RD"),
    get_species_stats(copenhagen_abundances$RL, "RL")
)

# Threshold on prevalence/abundance
RA_prevnum_filter <- species_stats[species_stats$site=="RA", ]$prevalence_and_numerosity > 0.25
RD_prevnum_filter <- species_stats[species_stats$site=="RD", ]$prevalence_and_numerosity > 0.25
RL_prevnum_filter <- species_stats[species_stats$site=="RL", ]$prevalence_and_numerosity > 0.25

# Threshold on median
RA_median_filter <- species_stats[species_stats$site=="RA", ]$nonzero_medians > 5
RD_median_filter <- species_stats[species_stats$site=="RD", ]$nonzero_medians > 5
RL_median_filter <- species_stats[species_stats$site=="RL", ]$nonzero_medians > 5

nodes_RA <- names(copenhagen_abundances$RA)[RA_prevnum_filter & RA_median_filter]
nodes_RD <- names(copenhagen_abundances$RD)[RD_prevnum_filter & RD_median_filter]
nodes_RL <- names(copenhagen_abundances$RL)[RL_prevnum_filter & RL_median_filter]

common_nodes <- intersect(intersect(nodes_RA, nodes_RD), nodes_RL)
# table(table(c(nodes_RA, nodes_RD, nodes_RL)))

#-------------------------------------------------------------------------------

# Only a subset of these nodes has been kept, corresponding to the species.
# This filtering step has been done in python, through the NCBI's API.
# The resulting nodes are stored in the file loaded below.
common_nodes <- read.csv("..\\data\\processed\\common_nodes.csv")
common_nodes <- common_nodes$common_nodes

#-------------------------------------------------------------------------------

abundances_RA <- copenhagen_abundances$RA[ , names(copenhagen_abundances$RA) %in% nodes_RA]
abundances_RD <- copenhagen_abundances$RD[ , names(copenhagen_abundances$RD) %in% nodes_RD]
abundances_RL <- copenhagen_abundances$RL[ , names(copenhagen_abundances$RL) %in% nodes_RL]

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

RA_clr_df <- data.frame(RA_clr)
names(RA_clr_df) <- names(RA_clr)
RD_clr_df <- data.frame(RD_clr)
names(RD_clr_df) <- names(RD_clr)
RL_clr_df <- data.frame(RL_clr)
names(RL_clr_df) <- names(RL_clr)

write.csv(RA_clr_df, "..\\data\\processed\\RA_clr.csv")
write.csv(RD_clr_df, "..\\data\\processed\\RD_clr.csv")
write.csv(RL_clr_df, "..\\data\\processed\\RL_clr.csv")

#-------------------------------------------------------------------------------

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

names(RA_corr) <- nodes_RA
names(RD_corr) <- nodes_RD
names(RL_corr) <- nodes_RL

RA_corr <- RA_corr[common_nodes, common_nodes]
RD_corr <- RD_corr[common_nodes, common_nodes]
RL_corr <- RL_corr[common_nodes, common_nodes]

RA_pval <- RA_corr_test$p[common_nodes, common_nodes]
RD_pval <- RD_corr_test$p[common_nodes, common_nodes]
RL_pval <- RL_corr_test$p[common_nodes, common_nodes]
diag(RA_pval) <- 1
diag(RD_pval) <- 1
diag(RL_pval) <- 1

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

RA_r_thr <- density_to_correlation_threshold(RA_corr, RA_pval, 1/100)
RD_r_thr <- density_to_correlation_threshold(RD_corr, RD_pval, 1/100)
RL_r_thr <- density_to_correlation_threshold(RL_corr, RL_pval, 1/100)

#-------------------------------------------------------------------------------

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

adjacencies <- list(
    "RA" = RA_adj,
    "RD" = RD_adj,
    "RL" = RL_adj
)
write_rds(adjacencies, file="..\\data\\processed\\adjacencies.rds")
