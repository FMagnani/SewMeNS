library(ggplot2)
library(data.table)
library(ggalluvial)
library(stringr)
library(mclust)
library(gridExtra)
library(VennDiagram)
library(grid)

# Name of variable: communities
load("..\\data\\processed\\seasonalCommunities_res_0.RData")

taxonomy <- read.csv("..\\data\\processed\\sp_name_idx.csv")
taxonomy$X <- NULL

#-------------------------------------------------------------------------------

# List the communities with more than 9 species

large_Hot_RA <- names(table(communities$Hot_RA))[table(communities$Hot_RA)>9]
large_Hot_RD <- names(table(communities$Hot_RD))[table(communities$Hot_RD)>9]
large_Hot_RL <- names(table(communities$Hot_RL))[table(communities$Hot_RL)>9]

large_Cold_RA <- names(table(communities$Cold_RA))[table(communities$Cold_RA)>9]
large_Cold_RD <- names(table(communities$Cold_RD))[table(communities$Cold_RD)>9]
large_Cold_RL <- names(table(communities$Cold_RL))[table(communities$Cold_RL)>9]

#-------------------------------------------------------------------------------

# In this section we plot one flowchart with 4 columns
#     Column 1: Taxonomy (es. the class)
#     Column 2: Community (RA)
#     Column 3: Community (RD)
#     Column 4: Community (RL)

alluvial_data_hot <- data.frame(
    "RA" = communities$Hot_RA,
    "RD" = communities$Hot_RD,
    "RL" = communities$Hot_RL
)

alluvial_data_cold <- data.frame(
    "RA" = communities$Cold_RA,
    "RD" = communities$Cold_RD,
    "RL" = communities$Cold_RL
)

# Argument "season" is only to set the title of the plot
# You may want to change: 
#     - The taxonomic level used 
#     - Filters on the numerosity of the communities to be displayed
get_flux_plot <- function(alluvial_data, season){

    title <- paste("Taxonomic composition of network communities (", season, " season)", sep='')
    
    # TAXONOMIC LEVEL
    tax_lvl <- "cl_name"
    taxonomic_info <- c()
    for(sp in names(communities$Hot_RA)){
        taxonomic_info <- c(taxonomic_info, taxonomy[taxonomy[["sp_name"]] == sp, tax_lvl])
    }
    alluvial_data$tax <- taxonomic_info
    
    alluvial_data$flow <- paste(
        alluvial_data$tax,
        alluvial_data$RA, alluvial_data$RD, alluvial_data$RL,
        sep='_'
    )
    
    get_alluvial_data_freq <- function(flow, colnames, n_cols, separator='_'){
        
        flow_table <- table(flow)
        alluvial_data_freq <- data.frame(
            "flow" = names(flow_table),
            "numerosity" = unname(flow_table)
        )
        alluvial_data_freq$numerosity.Var1 <- NULL
        names(alluvial_data_freq) <- c("flow", "numerosity")
        
        flow_df <- as.data.frame(str_split_fixed(alluvial_data_freq$flow, separator, n_cols))
        names(flow_df) <- colnames
        alluvial_data_freq[colnames] <- flow_df
        
        alluvial_data_freq <- alluvial_data_freq[order(-alluvial_data_freq$numerosity), ]
        
        return(alluvial_data_freq)
        
    }
    
    cols <- c("Class", "RA", "RD", "RL")
    alluvial_freq <- get_alluvial_data_freq(
        alluvial_data$flow, 
        cols, length(cols), separator="_"
    )
    
    # DISPLAY FILTER

    # f1: Numerosity of community: at least 10 species
    # f2, f3, f4: Exclude community '0' from the plot (isolated nodes)
    f1 <- (alluvial_freq$numerosity>9)
    f2 <- (alluvial_freq[["RA"]] != "0")
    f3 <- (alluvial_freq[["RD"]] != "0")
    f4 <- (alluvial_freq[["RL"]] != "0")
    filt <- f1&f2&f3&f4
    
    plt <- ggplot(
        data = alluvial_freq[filt, ],
        aes(
            y=numerosity,
            axis1=Class, axis2=RA, axis3=RD, axis4=RL
        )
    ) + scale_x_discrete(
        limits = cols,
        expand = c(.2, .05)
    ) + ggtitle(
        title
    ) + xlab("") + ylab("") + geom_alluvium(
        aes(fill=Class), width = 0.6
    ) + geom_stratum(
        width = 0.6, 
        alpha=0.1, reverse = TRUE, color = "black", linewidth = 0.7
    ) + geom_text(
        stat = "stratum", aes(label = after_stat(stratum), fontface="bold"), size=5
    ) + theme_minimal() + guides(
        fill=guide_legend(title="Taxonomic class")
    ) + theme(
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        plot.title = element_text(size=18)
    )
    
    return(plt)

}
    
get_flux_plot(alluvial_data_hot, "Hot")
get_flux_plot(alluvial_data_cold, "Cold")

#-------------------------------------------------------------------------------

# In this section we plot three separated flowcharts with 2 columns each:
#     Column 1: Community Hot
#     Column 2: Community Cold
# One of these for each site

plot_site_interseason <- function(site){

    alluvial_data <- data.frame(
        "Hot" = communities[[paste("Hot", site, sep='_')]],
        "Cold" = communities[[paste("Cold", site, sep='_')]]
    )
    
    alluvial_data$flow <- paste(
        alluvial_data$Hot, alluvial_data$Cold, sep='_'
    )
    
    cols <- c("Hot", "Cold")
    alluvial_freq <- get_alluvial_data_freq(
        alluvial_data$flow, 
        cols, length(cols), separator="_"
    )
    
    # Take largest flux
    f1 <- (alluvial_freq$numerosity>4)
    f2 <- (alluvial_freq[["Hot"]] != "0")
    f3 <- (alluvial_freq[["Cold"]] != "0")
    filt <- f1&f2&f3
    
    plt <- ggplot(
        data = alluvial_freq[filt,],
        aes(
            y=numerosity, axis1=Hot, axis2=Cold
        )
    ) + scale_x_discrete(
        limits = cols,
        expand = c(.2, .05)
    ) + ggtitle(
        site
    ) + xlab("") + ylab("") + geom_alluvium(
        aes(fill=Hot), width = 0.25
    ) + geom_stratum(
        width = 0.25, 
        alpha=0.1, reverse = TRUE, color = "black", linewidth = 0.7
    ) + geom_text(
        stat = "stratum", aes(label = after_stat(stratum), fontface="bold"), size=5
    ) + theme_minimal() + guides(
        fill=guide_legend(title="Community")
    ) + theme(
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15),
        plot.title = element_text(size=18),
        legend.position = "none"
    )
    
    return(plt)

}

RA_plt <- plot_site_interseason("RA")
RD_plt <- plot_site_interseason("RD")
RL_plt <- plot_site_interseason("RL")

title_grob <- textGrob(
    "Correspondence between seasonal communities",
    gp = gpar(fontsize = 20)
)
grid.arrange(RA_plt, RD_plt, RL_plt, ncol=3, top=title_grob)
