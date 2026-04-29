library(igraph)

seasonal_adjacencies <- readRDS("..\\data\\processed\\seasonal_adjacencies.rds")

adjacencies <- seasonal_adjacencies$cold
species <- rownames(adjacencies$RA)

#-------------------------------------------------------------------------------
# Links in common

# Multiply adjacencies: common links are entries equal to 1
common_links_RA_RD <- sum(adjacencies$RA*adjacencies$RD == 1)
common_links_RA_RL <- sum(adjacencies$RA*adjacencies$RL == 1)
common_links_RD_RL <- sum(adjacencies$RD*adjacencies$RL == 1)
tot_links <- sum(abs(adjacencies$RA))
message(100*common_links_RA_RD/tot_links) # Cold: 37.5 %    Hot: 40.1 %
message(100*common_links_RA_RL/tot_links) # Cold: 36.5 %    Hot: 38.8 %
message(100*common_links_RD_RL/tot_links) # Cold: 43.6 %    Hot: 42.8 %

# Average
(40.1 + 38.8 + 42.8)/3 # Hot:  40.6 %
(37.5 + 36.5 + 43.6)/3 # Cold: 39.2 %

common_links_RA_RA <- sum(seasonal_adjacencies$hot$RA*seasonal_adjacencies$cold$RA == 1)
common_links_RD_RD <- sum(seasonal_adjacencies$hot$RD*seasonal_adjacencies$cold$RD == 1)
common_links_RL_RL <- sum(seasonal_adjacencies$hot$RL*seasonal_adjacencies$cold$RL == 1)
tot_links <- sum(abs(adjacencies$RA))
message(100*common_links_RA_RA/tot_links) # 24.6 %
message(100*common_links_RD_RD/tot_links) # 33.7 %
message(100*common_links_RL_RL/tot_links) # 39.8 %

# Average
(24.6+33.7+39.8)/3 # 32.7 %

#-------------------------------------------------------------------------------
# Isolated nodes

isolated_RA <- species[rowSums(abs(adjacencies$RA)) == 0]
isolated_RD <- species[rowSums(abs(adjacencies$RD)) == 0]
isolated_RL <- species[rowSums(abs(adjacencies$RL)) == 0]
all_isolated <- length(unique(c(isolated_RA, isolated_RD, isolated_RL)))
common_isolated <- table(table(c(isolated_RA, isolated_RD, isolated_RL)))[[3]]
message(100*common_isolated/all_isolated) # Cold: 15.7 %    Hot: 13.6 %

#-------------------------------------------------------------------------------
# Components

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

size_lcc_RA <- sort(table(comp_RA$membership), decreasing=TRUE)[[1]]
size_lcc_RD <- sort(table(comp_RD$membership), decreasing=TRUE)[[1]]
size_lcc_RL <- sort(table(comp_RL$membership), decreasing=TRUE)[[1]]

idx_lcc_RA <- names(sort(table(comp_RA$membership), decreasing=TRUE)[1])
idx_lcc_RD <- names(sort(table(comp_RD$membership), decreasing=TRUE)[1])
idx_lcc_RL <- names(sort(table(comp_RL$membership), decreasing=TRUE)[1])

lcc_RA <- names(comp_RA$membership)[comp_RA$membership == idx_lcc_RA]
lcc_RD <- names(comp_RD$membership)[comp_RD$membership == idx_lcc_RD]
lcc_RL <- names(comp_RL$membership)[comp_RL$membership == idx_lcc_RL]

common_lcc <- table(table(c(lcc_RA, lcc_RD, lcc_RL)))[[3]]
message(100*common_lcc/length(species)) # Cold: 58.1 %    Hot: 65.5 %

#-------------------------------------------------------------------------------
# Correlation of node degrees

RA_degrees_hot <- rowSums(abs(seasonal_adjacencies$hot$RA))
RD_degrees_hot <- rowSums(abs(seasonal_adjacencies$hot$RD))
RL_degrees_hot <- rowSums(abs(seasonal_adjacencies$hot$RL))

RA_degrees_cold <- rowSums(abs(seasonal_adjacencies$cold$RA))
RD_degrees_cold <- rowSums(abs(seasonal_adjacencies$cold$RD))
RL_degrees_cold <- rowSums(abs(seasonal_adjacencies$cold$RL))

cor(RA_degrees_hot, RD_degrees_hot) # Hot: 0.62
cor(RA_degrees_hot, RL_degrees_hot) # Hot: 0.62
cor(RD_degrees_hot, RL_degrees_hot) # Hot: 0.63
# Average 0.62

cor(RA_degrees_cold, RD_degrees_cold) # Cold: 0.63
cor(RA_degrees_cold, RL_degrees_cold) # Cold: 0.57
cor(RD_degrees_cold, RL_degrees_cold) # Cold: 0.70
# Average 0.63

cor(RA_degrees_hot, RA_degrees_cold) # 0.40
cor(RD_degrees_hot, RD_degrees_cold) # 0.56
cor(RL_degrees_hot, RL_degrees_cold) # 0.69
# Average 0.55

#-------------------------------------------------------------------------------
# Common first quantile



#-------------------------------------------------------------------------------
# Connectivity

# N. isolated nodes
length(isolated_RA) # Cold: 416    Hot: 303
length(isolated_RD) # Cold: 447    Hot: 268
length(isolated_RL) # Cold: 301    Hot: 427

# N. connected components (no isolated nodes)
comp_RA$no - length(isolated_RA) # Cold: 42    Hot: 36
comp_RD$no - length(isolated_RD) # Cold: 52    Hot: 39
comp_RL$no - length(isolated_RL) # Cold: 34    Hot: 35

# Size of lcc
size_lcc_RA # Cold: 1570    Hot: 1773
size_lcc_RD # Cold: 1497    Hot: 1784
size_lcc_RL # Cold: 1760    Hot: 1626
