library(psych)
library(readr)
library(igraph)

copenhagen_abundances <- readRDS("..\\data\\processed\\abun_plants.rds")

common_nodes <- read.csv("..\\data\\processed\\sp_name_idx.csv")
common_nodes$X <- NULL
common_nodes <- common_nodes$sp_name

metadata <- read.csv("..\\data\\raw\\cph_sewage_2020_metadata.csv", sep=';')

#-------------------------------------------------------------------------------
# Separate the data into two seasons

split_by_period <- function(site_data, start_date, end_date){

    site_dates <- c()
    for(sample_name in rownames(site_data)){
        sample_date <- metadata$Date.of.isolation[metadata$Complete.name == sample_name]
        site_dates <- c(site_dates, as.Date(sample_date, "%Y.%m.%d"))
    }
        
    site_data <- site_data[
        ((site_dates>=start_date) & (site_dates<end_date)), 
    ]
    
    return(site_data)
    
}

cold_start <- as.Date("28.10.2019", "%d.%m.%Y")
cold_end <- as.Date("25.04.2020", "%d.%m.%Y")
hot_start <- as.Date("25.04.2020", "%d.%m.%Y")
hot_end <- as.Date("01.10.2020", "%d.%m.%Y")

RA_cold <- split_by_period(copenhagen_abundances$RA, cold_start, cold_end)
RD_cold <- split_by_period(copenhagen_abundances$RD, cold_start, cold_end)
RL_cold <- split_by_period(copenhagen_abundances$RL, cold_start, cold_end)
cold_abundances <- list("RA"=RA_cold, "RD"=RD_cold, "RL"=RL_cold)

RA_hot <- split_by_period(copenhagen_abundances$RA, hot_start, hot_end)
RD_hot <- split_by_period(copenhagen_abundances$RD, hot_start, hot_end)
RL_hot <- split_by_period(copenhagen_abundances$RL, hot_start, hot_end)
hot_abundances <- list("RA"=RA_hot, "RD"=RD_hot, "RL"=RL_hot)

#-------------------------------------------------------------------------------

compute_CLR_corr_and_pval <- function(abundances){

    abundances_RA <- abundances$RA[ , names(abundances$RA) %in% common_nodes]
    abundances_RD <- abundances$RD[ , names(abundances$RD) %in% common_nodes]
    abundances_RL <- abundances$RL[ , names(abundances$RL) %in% common_nodes]

    message("correct dataset")    
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

    message("Correct RA...")    
    RA_correct <- correct_dataset(abundances_RA)
    message("Correct RD...")    
    RD_correct <- correct_dataset(abundances_RD)
    message("Correct RL...")    
    RL_correct <- correct_dataset(abundances_RL)
    
    message("clr computation")    
    apply_CLR_to_dataset <- function(counts){
        
        # MARGIN = 1 for rows, = 2 for columns
        ref <- apply(counts, 1, function(x) mean(log(x)))
        transformed <- sweep(log(counts), 1, ref)
        return(transformed)
        
    }
    
    RA_clr <- apply_CLR_to_dataset(RA_correct)
    RD_clr <- apply_CLR_to_dataset(RD_correct)
    RL_clr <- apply_CLR_to_dataset(RL_correct)

    message("correlation RA")    
    RA_corr_test <- corr.test(
        x=RA_clr, 
        use="pairwise", method="pearson",
        adjust="bonferroni",
        ci=FALSE
    )
    
    message("correlation RD")    
    RD_corr_test <- corr.test(
        x=RD_clr, 
        use="pairwise", method="pearson",
        adjust="bonferroni",
        ci=FALSE
    )
    
    message("correlation RL")    
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
    
    RA_pval <- RA_corr_test$p
    RD_pval <- RD_corr_test$p
    RL_pval <- RL_corr_test$p
    RA_pval[lower.tri(RA_pval)] <- t(RA_pval)[lower.tri(RA_pval)]
    RD_pval[lower.tri(RD_pval)] <- t(RD_pval)[lower.tri(RD_pval)]
    RL_pval[lower.tri(RL_pval)] <- t(RL_pval)[lower.tri(RL_pval)]
    diag(RA_pval) <- 1
    diag(RD_pval) <- 1
    diag(RL_pval) <- 1

    return(list(
        "RA_corr"=RA_corr, "RA_pval"=RA_pval,
        "RD_corr"=RD_corr, "RD_pval"=RD_pval,
        "RL_corr"=RL_corr, "RL_pval"=RL_pval
    ))

}

cold_clr_corr_pval <- compute_CLR_corr_and_pval(cold_abundances)
hot_clr_corr_pval <- compute_CLR_corr_and_pval(hot_abundances)

# Significative pvals are much lower in this case...
# If we were to take only significative p-vals, the densities would be as
#     shown below
100*sum(cold_clr_corr_pval$RA_pval<0.05)/(2170*2170) # 1.1 %
100*sum(cold_clr_corr_pval$RD_pval<0.05)/(2170*2170) # 5.1 %
100*sum(cold_clr_corr_pval$RL_pval<0.05)/(2170*2170) # 1.2 %

100*sum(hot_clr_corr_pval$RA_pval<0.05)/(2170*2170) # 0.96 %
100*sum(hot_clr_corr_pval$RD_pval<0.05)/(2170*2170) # 0.86 %
100*sum(hot_clr_corr_pval$RL_pval<0.05)/(2170*2170) # 1.28 %

density_to_correlation_threshold <- function(site_corr, site_pval, density_thr){
    
    N_links <- length(site_corr)
    index_thr <- 1 + as.integer(density_thr*N_links)
    
    # All correlations are treated as significative
    significative_corr <- site_corr
    
    corr_thr <- abs(
        site_corr[
            order(abs(
                significative_corr
            ), decreasing=TRUE)
        ][[index_thr]]
    )
    
    return(corr_thr)
    
}

density_to_adjacencies <- function(clr_corr_pval, density_thr){
    
    RA_corr <- clr_corr_pval$RA_corr
    RD_corr <- clr_corr_pval$RD_corr
    RL_corr <- clr_corr_pval$RL_corr
    RA_pval <- clr_corr_pval$RA_pval
    RD_pval <- clr_corr_pval$RD_pval
    RL_pval <- clr_corr_pval$RL_pval
    
    RA_r_thr <- density_to_correlation_threshold(RA_corr, RA_pval, density_thr)
    RD_r_thr <- density_to_correlation_threshold(RD_corr, RD_pval, density_thr)
    RL_r_thr <- density_to_correlation_threshold(RL_corr, RL_pval, density_thr)
    message(RA_r_thr)
    message(RD_r_thr)
    message(RL_r_thr)
    
    RA_adj <- RA_corr
    RA_adj[(abs(RA_corr) < RA_r_thr)] <- 0
    #RA_adj[RA_pval > 0.05] <- 0
    RA_adj[RA_adj > 0] <-  1
    RA_adj[RA_adj < 0] <- -1
    
    RD_adj <- RD_corr
    RD_adj[(abs(RD_corr) < RD_r_thr)] <- 0
    #RD_adj[RD_pval > 0.05] <- 0
    RD_adj[RD_adj > 0] <-  1
    RD_adj[RD_adj < 0] <- -1
    
    RL_adj <- RL_corr
    RL_adj[(abs(RL_corr) < RL_r_thr)] <- 0
    #RL_adj[RL_pval > 0.05] <- 0
    RL_adj[RL_adj > 0] <-  1
    RL_adj[RL_adj < 0] <- -1
    
    adjacencies <- list("RA" = RA_adj, "RD" = RD_adj, "RL" = RL_adj)
    
    return(adjacencies)
    
}

# Same value of the year network
adjacencies_cold <- density_to_adjacencies(cold_clr_corr_pval, 5.6/100)
adjacencies_hot <- density_to_adjacencies(hot_clr_corr_pval, 5.6/100)
write_rds(
    list("cold" = adjacencies_cold, "hot" = adjacencies_hot), 
    file="..\\data\\processed\\seasonal_adjacencies.rds"
)