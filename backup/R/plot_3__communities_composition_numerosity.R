library(ggplot2)
library(data.table)
library(ggalluvial)
library(stringr)
library(mclust)
library(gridExtra)
library(VennDiagram)

load("..\\data\\processed\\communities_res_0.RData")

taxonomy <- read.csv("..\\data\\processed\\sp_name_idx.csv")
taxonomy$X <- NULL

large_RA <- names(table(communities$RA))[table(communities$RA)>9]
large_RD <- names(table(communities$RD))[table(communities$RD)>9]
large_RL <- names(table(communities$RL))[table(communities$RL)>9]

alluvial_data <- data.frame(
    "RA" = communities$RA,
    "RD" = communities$RD,
    "RL" = communities$RL
)

tax_lvl <- "cl_name"
tax <- c()
for(sp in names(communities$RA)){
    tax <- c(tax, taxonomy[taxonomy[["sp_name"]] == sp, tax_lvl])
}
alluvial_data[["taxa"]] <- tax

plot_site_composition <- function(site_communities, site){
    
    community_col <- c()
    taxa_col <- c()
    N_col <- c()
    for(comm_name in site_communities[2:length(site_communities)]){
        composition <- table(alluvial_data$taxa[alluvial_data[[site]]==comm_name])
        n_others <- sum(composition[composition<5])
        composition <- composition[composition>5]
        composition["Others"] <- n_others
        
        community_col <- c(community_col, rep(comm_name, length(composition)))
        taxa_col <- c(taxa_col, names(composition))
        N_col <- c(N_col, unname(composition))
    }
    site_composition_df <- data.frame(
        "Community" = community_col,
        "Taxa" = taxa_col,
        "N" = N_col
    )
    site_composition_df$Community <- factor(
        site_composition_df$Community,
        levels = sort(as.numeric(unique(site_composition_df$Community)))
    )
    
    my_colors <- scales::hue_pal()(length(unique(site_composition_df$Taxa)))
    names(my_colors) <- unique(site_composition_df$Taxa)
    my_colors["Others"] <- "black"
    
    plt <- ggplot(
        site_composition_df, aes(fill=Taxa, y=N, x=Community)
    ) + geom_bar(
        position='stack', stat='identity'
    ) + scale_fill_manual(values = my_colors) + ggtitle(
        site
    ) + theme_minimal() + guides(
        fill=guide_legend(title="Taxonomic class")
    ) + xlab("Community") + ylab("Number of species") + theme(
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text = element_text(size=15), 
        axis.title = element_text(size=15), 
        #    axis.text.y = element_text(size=15),
        plot.title = element_text(size=18),
    ) + theme(legend.position = "none")

    return(plt)
    
}

RA_plot <- plot_site_composition(large_RA, "RA")
RD_plot <- plot_site_composition(large_RD, "RD")
RL_plot <- plot_site_composition(large_RL, "RL")

title_grob <- textGrob(
    "Composition and numerosity of largest communities",
    gp = gpar(fontsize = 20)
)
grid.arrange(
    RA_plot, RD_plot, RL_plot, ncol=3, top=title_grob
)


