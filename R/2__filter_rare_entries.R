library(readr)

copenhagen_abundances <- readRDS("..\\data\\processed\\abun_plants.rds")

#-------------------------------------------------------------------------------

get_species_prevalence <- function(df, site){
    # To each entry attach the fraction of samples in which 
    #   it is detected with more than 10 counts
    
    prevalence <- c()
    entries = names(df)
    
    for(e in entries){
        
        abundances <-  df[ ,e]
        
        prevalence <- c(
            prevalence, 
            length(abundances[abundances >= 10])/length(abundances)
        )

    }
    
    return(
        data.frame(
            "name" = unlist(entries),
            "prevalence" = prevalence,
            "site" = site
        )
    )
    
}

results <- rbind(
    get_species_prevalence(copenhagen_abundances$RA, "RA"),
    get_species_prevalence(copenhagen_abundances$RD, "RD"),
    get_species_prevalence(copenhagen_abundances$RL, "RL")
)

# Threshold on prevalence
RA_prev_filter <- results[results$site=="RA", ]$prevalence > 0.25
RD_prev_filter <- results[results$site=="RD", ]$prevalence > 0.25
RL_prev_filter <- results[results$site=="RL", ]$prevalence > 0.25

filtered_RA <- names(copenhagen_abundances$RA)[RA_prev_filter]
filtered_RD <- names(copenhagen_abundances$RD)[RD_prev_filter]
filtered_RL <- names(copenhagen_abundances$RL)[RL_prev_filter]

site_filtered <- list("RA" = filtered_RA, "RD" = filtered_RD, "RL" = filtered_RL)
write_rds(site_filtered, file="..\\data\\processed\\site_filtered.rds")

common_nodes <- intersect(intersect(filtered_RA, filtered_RD), filtered_RL)
table(table(c(filtered_RA, filtered_RD, filtered_RL)))
write.csv(common_nodes, "..\\data\\processed\\common_nodes.csv")
