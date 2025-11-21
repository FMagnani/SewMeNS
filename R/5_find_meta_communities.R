library(ggplot2)
library(ggalluvial)
library(stringr)
library(mclust)

communities <- readRDS("..\\data\\processed\\communities.rds")

# Compute flux
alluvial_data_long <- data.frame(
    "col1" = communities$RA,
    "col2" = communities$RD,
    "col3" = communities$RL
)

alluvial_data_long$flow <- paste(
    alluvial_data_long$col1, alluvial_data_long$col2, alluvial_data_long$col3,
    sep='_'
)

flow_table <- table(alluvial_data_long$flow)
alluvial_data_freq <- data.frame(
    "flow" = names(flow_table),
    "numerosity" = unname(flow_table)
)
alluvial_data_freq$numerosity.Var1 <- NULL
names(alluvial_data_freq) <- c("flow", "numerosity")

flow_df <- as.data.frame(str_split_fixed(alluvial_data_freq$flow, '_', 3))
names(flow_df) <- c("col1", "col2", "col3")
alluvial_data_freq[c("col1", "col2", "col3")] <- flow_df

alluvial_data_freq <- alluvial_data_freq[order(-alluvial_data_freq$numerosity), ]

# Take largest flux
mcf1 <- (alluvial_data_freq$numerosity>9)
mcf2 <- (alluvial_data_freq$flow != "0_0_0")
mcf3 <- (alluvial_data_freq$col1 != "0")
mcf4 <- (alluvial_data_freq$col2 != "0")
mcf5 <- (alluvial_data_freq$col3 != "0")
metacomm_filter <- mcf1 & mcf2 & mcf3 & mcf4 & mcf5

metacomm_names <- c()
for(fl in alluvial_data_freq$flow){
    if(fl %in% alluvial_data_freq$flow[metacomm_filter]){
        metacomm_names <- c(metacomm_names, fl)
    }else{
        metacomm_names <- c(metacomm_names, 'x')
    }
}
alluvial_data_freq$metacomm <- metacomm_names

#-------------------------------------------------------------------------------

# NETWORK COMMUNITIES

selected_comms <- c("0","1","2","3","4","5","6","7","8","9","10")
RA_largest_comms <- selected_comms
RD_largest_comms <- selected_comms
RL_largest_comms <- selected_comms

alluvial_data_freq$col1[!(alluvial_data_freq$col1 %in% RA_largest_comms)] <- 's'
alluvial_data_freq$col2[!(alluvial_data_freq$col2 %in% RD_largest_comms)] <- 's'
alluvial_data_freq$col3[!(alluvial_data_freq$col3 %in% RL_largest_comms)] <- 's'


# Remove zero community
rm1 <- (alluvial_data_freq$col1 != '0')
rm2 <- (alluvial_data_freq$col2 != '0')
rm3 <- (alluvial_data_freq$col3 != '0')

# Remove small communities
rm4 <- (alluvial_data_freq$col1 != 's')
rm5 <- (alluvial_data_freq$col2 != 's')
rm6 <- (alluvial_data_freq$col3 != 's')

rm_filter_comm <- rm1&rm2&rm3&rm4&rm5&rm6


# Remove zero-zero, small-small and zero-small flux
z1 <- (alluvial_data_freq$col1 == "0")
z2 <- (alluvial_data_freq$col2 == "0")
z3 <- (alluvial_data_freq$col3 == "0")
s1 <- (alluvial_data_freq$col1 == "s")
s2 <- (alluvial_data_freq$col2 == "s")
s3 <- (alluvial_data_freq$col3 == "s")
rm_filter_flux <- !((z1|s1) & (z2|s2) & (z3|s3))

#-------------------------------------------------------------------------------

plt_communities <- ggplot(
    data = alluvial_data_freq[rm_filter_flux,],
    aes(
        axis1=col1, axis2=col2, axis3=col3, y=numerosity
    )
) + scale_x_discrete(
    limits = c("RA", "RD", "RL"), 
    expand = c(.2, .05)
) + xlab("") + ylab("") + geom_alluvium(
    aes(fill=col1), width = 1/3.5
) + geom_stratum(
    width = 1/3.5, alpha=0.5, reverse = TRUE, color = "black", linewidth = 0.7
) + geom_text(
    stat = "stratum", aes(label = after_stat(stratum), fontface="bold"), size=4
) + theme_minimal() + ggtitle(
    "Largest network communities in the three sites"
) + guides(
    fill=guide_legend(title="Community")
) + theme(
    legend.title=element_text(size=15),
    legend.text=element_text(size=15),
    axis.text.x = element_text(size=15), 
    axis.text.y = element_text(size=15),
    plot.title = element_text(size=18)
) #+ geom_text(stat = "stratum", aes(label = after_stat(stratum)), size= 5, color = "black")

plt_communities

#-------------------------------------------------------------------------------

# ADJUSTED RAND INDEX

ari_RL_RA <- adjustedRandIndex(communities$RL, communities$RA)
ari_RL_RD <- adjustedRandIndex(communities$RL, communities$RD)
ari_RA_RD <- adjustedRandIndex(communities$RA, communities$RD)
#message(ari_RL_RA) # 0.38
#message(ari_RL_RD) # 0.46
#message(ari_RA_RD) # 0.40

rm1 <- (communities$RL != '0')
rm2 <- (communities$RA != '0')
rm3 <- (communities$RD != '0')
communities_RL_nozero <- communities$RL[rm1&rm2&rm3]
communities_RA_nozero <- communities$RA[rm1&rm2&rm3]
communities_RD_nozero <- communities$RD[rm1&rm2&rm3]

ari_RL_RA <- adjustedRandIndex(communities_RL_nozero, communities_RA_nozero)
ari_RL_RD <- adjustedRandIndex(communities_RL_nozero, communities_RD_nozero)
ari_RA_RD <- adjustedRandIndex(communities_RA_nozero, communities_RD_nozero)
#message(ari_RL_RA) # 0.63
#message(ari_RL_RD) # 0.58
#message(ari_RA_RD) # 0.67

#-------------------------------------------------------------------------------

largest_flows = union(
    union(
        alluvial_data_freq$col1[alluvial_data_freq$metacomm != 'x'],
        alluvial_data_freq$col2[alluvial_data_freq$metacomm != 'x']
    ),
    alluvial_data_freq$col3[alluvial_data_freq$metacomm != 'x']
)

rm1 <- (alluvial_data_freq$col1 %in% largest_flows)
rm2 <- (alluvial_data_freq$col2 %in% largest_flows)
rm3 <- (alluvial_data_freq$col3 %in% largest_flows)
rm4 <- (alluvial_data_freq$metacomm != 'x')

plt_meta <- ggplot(
    data = alluvial_data_freq[rm1&rm2&rm3&rm4,],
    aes(
        axis1=metacomm, axis2=col1, axis3=col2, axis4=col3, y=numerosity
    )
) + scale_x_discrete(
    limits = c("Meta", "RA", "RD", "RL"), 
    expand = c(.2, .05)
) + xlab("") + ylab("") + geom_alluvium(
    aes(fill=metacomm), width = 1/3.5
) + geom_stratum(
    width = 1/3.5, alpha=0.5, reverse = TRUE, color = "black", linewidth = 0.7
) + geom_text(
    stat = "stratum", aes(label = after_stat(stratum), fontface="bold"), size=4
) + theme_minimal() + ggtitle(
    "Network meta-communities"
) + guides(
    fill=guide_legend(title="Community")
) + theme(
    legend.title=element_text(size=15),
    legend.text=element_text(size=15),
    axis.text.x = element_text(size=15), 
    axis.text.y = element_text(size=15),
    plot.title = element_text(size=18)
) 

plt_meta

#-------------------------------------------------------------------------------

get_meta_community <- function(metacomm_name){
    
    ra <- strsplit(metacomm_name, '_')[[1]][1]
    rd <- strsplit(metacomm_name, '_')[[1]][2]
    rl <- strsplit(metacomm_name, '_')[[1]][3]
    
    comm_ra = names(communities$RA)[(communities$RA==ra)]
    comm_rd = names(communities$RD)[(communities$RD==rd)]
    comm_rl = names(communities$RL)[(communities$RL==rl)]
    
    comm = intersect(intersect(comm_rl, comm_ra), comm_rd)
    
    return(comm)  
    
}

species_list <- c()
comm_names_list <- c()
for(comm_name in alluvial_data_freq$metacomm[alluvial_data_freq$metacomm != 'x']){
    comm_species <- get_meta_community(comm_name)
    species_list <- c(species_list, unlist(comm_species))
    comm_names_list <- c(comm_names_list, rep(comm_name, times=length(comm_species)))
}

meta_communities_df <- data.frame(species=unlist(species_list), comm=unlist(comm_names_list))
write.csv(meta_communities_df, "..\\data\\processed\\meta_communities.csv")
