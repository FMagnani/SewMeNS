library(ggplot2)
library(ggalluvial)
library(stringr)
library(mclust)
library(gridExtra)

communities <- readRDS("..\\data\\processed\\seasonal_communities.rds")

#-------------------------------------------------------------------------------

print_ARI <- function(comm1, comm2){
    ari <- adjustedRandIndex(comm1, comm2)
    message(round(ari, 2))

    rm1 <- (comm1 != '0')
    rm2 <- (comm2 != '0')
    comm1_nozero <- comm1[rm1&rm2]
    comm2_nozero <- comm2[rm1&rm2]

    ari_nozero <- adjustedRandIndex(comm1_nozero, comm2_nozero)
    message(round(ari_nozero, 2))
}

# Hot - Hot
# Avg: 0.29 0.45                                    zero no-zero
print_ARI(communities$Hot_RA, communities$Hot_RD) # 0.26 0.47
print_ARI(communities$Hot_RA, communities$Hot_RL) # 0.28 0.43
print_ARI(communities$Hot_RD, communities$Hot_RL) # 0.33 0.45

# Cold - Cold
# Avg: 0.31 0.41                                      zero no-zero
print_ARI(communities$Cold_RA, communities$Cold_RD) # 0.33 0.49
print_ARI(communities$Cold_RA, communities$Cold_RL) # 0.26 0.35
print_ARI(communities$Cold_RD, communities$Cold_RL) # 0.34 0.40

# Hot - Cold
# Avg: 0.21 0.36                                     zero no-zero
print_ARI(communities$Hot_RA, communities$Cold_RA) # 0.17 0.34
print_ARI(communities$Hot_RD, communities$Cold_RD) # 0.24 0.38
print_ARI(communities$Hot_RL, communities$Cold_RL) # 0.22 0.36

#-------------------------------------------------------------------------------

# Compute flux
alluvial_data_long <- data.frame(
    "H_RA" = communities$Hot_RA,
    "H_RD" = communities$Hot_RD,
    "H_RL" = communities$Hot_RL,
    "C_RA" = communities$Cold_RA,
    "C_RD" = communities$Cold_RD,
    "C_RL" = communities$Cold_RL
)

alluvial_data_long$flow_intersite_hot <- paste(
    alluvial_data_long$H_RA, alluvial_data_long$H_RD, alluvial_data_long$H_RL,
    sep='_'
)
alluvial_data_long$flow_intersite_cold <- paste(
    alluvial_data_long$C_RA, alluvial_data_long$C_RD, alluvial_data_long$C_RL,
    sep='_'
)

#-------------------------------------------------------------------------------

get_alluvial_data_freq <- function(flow, colnames, n_cols){

    flow_table <- table(flow)
    alluvial_data_freq <- data.frame(
        "flow" = names(flow_table),
        "numerosity" = unname(flow_table)
    )
    alluvial_data_freq$numerosity.Var1 <- NULL
    names(alluvial_data_freq) <- c("flow", "numerosity")
    
    flow_df <- as.data.frame(str_split_fixed(alluvial_data_freq$flow, '_', n_cols))
    names(flow_df) <- colnames
    alluvial_data_freq[colnames] <- flow_df
    
    alluvial_data_freq <- alluvial_data_freq[order(-alluvial_data_freq$numerosity), ]

    return(alluvial_data_freq)
        
}

intersite_hot_freq <- get_alluvial_data_freq(
    alluvial_data_long$flow_intersite_hot, c("RA", "RD", "RL"), 3
)

intersite_cold_freq <- get_alluvial_data_freq(
    alluvial_data_long$flow_intersite_cold, c("RA", "RD", "RL"), 3
)

#-------------------------------------------------------------------------------

name_metacommunities <- function(alluvial_data_freq, colnames){
    
    # Take largest flux
    mcf1 <- (alluvial_data_freq$numerosity>9)

    mcf2 <- (alluvial_data_freq[[colnames[[1]]]] != "0")
    mcf3 <- (alluvial_data_freq[[colnames[[2]]]] != "0")
    mcf4 <- (alluvial_data_freq[[colnames[[3]]]] != "0")
    metacomm_filter <- mcf1 & mcf2 & mcf3 & mcf4

    alluvial_data_freq$metacomm <- alluvial_data_freq$flow
    alluvial_data_freq$metacomm[ !metacomm_filter ] <- 'x'
    
    return(alluvial_data_freq)    
    
}

intersite_hot_freq <- name_metacommunities(intersite_hot_freq, c("RA", "RD", "RL"), 3)
intersite_cold_freq <- name_metacommunities(intersite_cold_freq, c("RA", "RD", "RL"), 3)

#-------------------------------------------------------------------------------

plot_meta_communities <- function(alluvial_data_freq, title){

    largest_flows = union(
        union(
            alluvial_data_freq$RA[alluvial_data_freq$metacomm != 'x'],
            alluvial_data_freq$RD[alluvial_data_freq$metacomm != 'x']
        ),
        alluvial_data_freq$RL[alluvial_data_freq$metacomm != 'x']
    )
    
    rm1 <- (alluvial_data_freq$RA %in% largest_flows)
    rm2 <- (alluvial_data_freq$RD %in% largest_flows)
    rm3 <- (alluvial_data_freq$RL %in% largest_flows)
    rm4 <- (alluvial_data_freq$metacomm != 'x')
    #rm4 <- (alluvial_data_freq$metacomm != '0_0_0')
    
    plt_meta <- ggplot(
        data = alluvial_data_freq[rm1&rm2&rm3&rm4,],
        aes(
            axis1=metacomm, axis2=RA, axis3=RD, axis4=RL, y=numerosity
        )
    ) + scale_x_discrete(
        limits = c("Meta", "RA", "RD", "RL"), 
        expand = c(.2, .05)
    )
        
    plt_meta <- plt_meta + ggtitle(
        title
    ) + xlab("") + ylab("") + geom_alluvium(
        aes(fill=metacomm), width = 1/3.5
    ) + geom_stratum(
        width = 1/3.5, alpha=0.5, reverse = TRUE, color = "black", linewidth = 0.7
    ) + geom_text(
        stat = "stratum", aes(label = after_stat(stratum), fontface="bold"), size=4
    ) + theme_minimal() + guides(
        fill=guide_legend(title="Community")
    ) + theme(
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        plot.title = element_text(size=18)
    ) 
    
    return(plt_meta)
        
}

plt_hot <- plot_meta_communities(intersite_hot_freq, "Meta communities - Hot season", 3)
plt_cold <- plot_meta_communities(intersite_cold_freq, "Meta communities - Cold season", 3)
grid.arrange(plt_hot, plt_cold, ncol = 1)

#-------------------------------------------------------------------------------

get_seasonal_metacommunity <- function(metacomm_name, season){
    
    ra <- strsplit(metacomm_name, '_')[[1]][1]
    rd <- strsplit(metacomm_name, '_')[[1]][2]
    rl <- strsplit(metacomm_name, '_')[[1]][3]

    if(season=="Hot"){
        comm_ra = names(communities$Hot_RA)[(communities$Hot_RA==ra)]
        comm_rd = names(communities$Hot_RD)[(communities$Hot_RD==rd)]
        comm_rl = names(communities$Hot_RL)[(communities$Hot_RL==rl)]
    }else if(season=="Cold"){
        comm_ra = names(communities$Cold_RA)[(communities$Cold_RA==ra)]
        comm_rd = names(communities$Cold_RD)[(communities$Cold_RD==rd)]
        comm_rl = names(communities$Cold_RL)[(communities$Cold_RL==rl)]
    }    
    
    comm = intersect(intersect(comm_rl, comm_ra), comm_rd)
    
    return(comm)  
}

get_seasonal_metacommunities <- function(alluvial_data_freq, season){
    species_list <- c()
    comm_names_list <- c()
    for(comm_name in alluvial_data_freq$metacomm[alluvial_data_freq$metacomm != 'x']){
        comm_species <- get_seasonal_metacommunity(comm_name, season)
        species_list <- c(species_list, unlist(comm_species))
        comm_names_list <- c(comm_names_list, rep(comm_name, times=length(comm_species)))
    }
    meta_communities_df <- data.frame(species=unlist(species_list), comm=unlist(comm_names_list))
    return(meta_communities_df)
}

hot_hot_metacommunities <- get_seasonal_metacommunities(intersite_hot_freq, "Hot")
cold_cold_metacommunities <- get_seasonal_metacommunities(intersite_cold_freq, "Cold")

write.csv(hot_hot_metacommunities, "..\\data\\processed\\hot_hot_metacommunities.csv")
write.csv(cold_cold_metacommunities, "..\\data\\processed\\cold_cold_metacommunities.csv")
