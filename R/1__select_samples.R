library(readr)
genomic_species_abundance <- readRDS("..\\data\\processed\\abun.rds")
metadata <- read.csv("..\\data\\raw\\cph_sewage_2020_metadata.csv", sep=';')

samples <- rownames(genomic_species_abundance)
metadata <- metadata[metadata$Complete.name %in% samples, ]

site_RA_filter <- (metadata$Plant == "RA")
site_RD_filter <- (metadata$Plant == "RD")
site_RL_filter <- (metadata$Plant == "RL")

RA_metadata <- metadata[site_RA_filter, ]
RD_metadata <- metadata[site_RD_filter, ]
RL_metadata <- metadata[site_RL_filter, ]

site_RA_sample_names <- metadata$Complete.name[site_RA_filter]
site_RD_sample_names <- metadata$Complete.name[site_RD_filter]
site_RL_sample_names <- metadata$Complete.name[site_RL_filter]

RA_abundance <- genomic_species_abundance[
    rownames(genomic_species_abundance) %in% site_RA_sample_names, 
]
RD_abundance <- genomic_species_abundance[
    rownames(genomic_species_abundance) %in% site_RD_sample_names, 
]
RL_abundance <- genomic_species_abundance[
    rownames(genomic_species_abundance) %in% site_RL_sample_names, 
]

RA_metadata$sample_abundance <- rowSums(RA_abundance)
RD_metadata$sample_abundance <- rowSums(RD_abundance)
RL_metadata$sample_abundance <- rowSums(RL_abundance)

collection_dates_RA <- unique(metadata$Date.of.isolation[site_RA_filter])
collection_dates_RD <- unique(metadata$Date.of.isolation[site_RD_filter])
collection_dates_RL <- unique(metadata$Date.of.isolation[site_RL_filter])

RA_selected_samples <- c()
for(date in collection_dates_RA){
    date_metadata <- RA_metadata[RA_metadata$Date.of.isolation==date, ]
    if(nrow(date_metadata) == 1){
        selected_sample <- date_metadata$Complete.name
    }else{
        selected_sample <- date_metadata$Complete.name[
            date_metadata$sample_abundance == max(date_metadata$sample_abundance)
        ]
    }
    RA_selected_samples <- c(RA_selected_samples, selected_sample)
}

RD_selected_samples <- c()
for(date in collection_dates_RD){
    date_metadata <- RD_metadata[RD_metadata$Date.of.isolation==date, ]
    if(nrow(date_metadata) == 1){
        selected_sample <- date_metadata$Complete.name
    }else{
        selected_sample <- date_metadata$Complete.name[
            date_metadata$sample_abundance == max(date_metadata$sample_abundance)
        ]
    }
    RD_selected_samples <- c(RD_selected_samples, selected_sample)
}

RL_selected_samples <- c()
for(date in collection_dates_RL){
    date_metadata <- RL_metadata[RL_metadata$Date.of.isolation==date, ]
    if(nrow(date_metadata) == 1){
        selected_sample <- date_metadata$Complete.name
    }else{
        selected_sample <- date_metadata$Complete.name[
            date_metadata$sample_abundance == max(date_metadata$sample_abundance)
        ]
    }
    RL_selected_samples <- c(RL_selected_samples, selected_sample)
}

RA_abundance <- RA_abundance[
    rownames(RA_abundance) %in% RA_selected_samples, 
]
RD_abundance <- RD_abundance[
    rownames(RD_abundance) %in% RD_selected_samples, 
]
RL_abundance <- RL_abundance[
    rownames(RL_abundance) %in% RL_selected_samples, 
]

RA_metadata <- RA_metadata[RA_metadata$complete_name %in% RA_selected_samples, ]
RD_metadata <- RD_metadata[RD_metadata$complete_name %in% RD_selected_samples, ]
RL_metadata <- RL_metadata[RL_metadata$complete_name %in% RL_selected_samples, ]

abun_plants <- list(
    "RA" = RA_abundance,
    "RD" = RD_abundance,
    "RL" = RL_abundance
)
write_rds(abun_plants, file="..\\data\\processed\\abun_plants.rds")
