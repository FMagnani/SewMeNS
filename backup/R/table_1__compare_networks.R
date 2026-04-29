library(igraph)

adjacencies <- readRDS("..\\data\\processed\\adjacencies.rds")
species <- rownames(adjacencies$RA)

#-------------------------------------------------------------------------------
# Links in common

# Multiply adjacencies: common links are entries equal to 1
common_links_RA_RD <- sum(adjacencies$RA*adjacencies$RD == 1)
common_links_RA_RL <- sum(adjacencies$RA*adjacencies$RL == 1)
common_links_RD_RL <- sum(adjacencies$RD*adjacencies$RL == 1)
tot_links <- sum(abs(adjacencies$RA))
message(100*common_links_RA_RD/tot_links) # 48.1 %
message(100*common_links_RA_RL/tot_links) # 56.3 %
message(100*common_links_RD_RL/tot_links) # 49.4 %

#-------------------------------------------------------------------------------
# Isolated nodes

isolated_RA <- species[rowSums(abs(adjacencies$RA)) == 0]
isolated_RD <- species[rowSums(abs(adjacencies$RD)) == 0]
isolated_RL <- species[rowSums(abs(adjacencies$RL)) == 0]
all_isolated <- length(unique(c(isolated_RA, isolated_RD, isolated_RL)))
common_isolated <- table(table(c(isolated_RA, isolated_RD, isolated_RL)))[[3]]
message(100*common_isolated/all_isolated) # 27.1 %

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
message(100*common_lcc/length(species)) # 62.5 %

#-------------------------------------------------------------------------------
# Correlation of node degrees

RA_degrees <- rowSums(abs(adjacencies$RA))
RD_degrees <- rowSums(abs(adjacencies$RD))
RL_degrees <- rowSums(abs(adjacencies$RL))

cor(RA_degrees, RD_degrees) # 0.69
cor(RA_degrees, RL_degrees) # 0.83
cor(RD_degrees, RL_degrees) # 0.72

most_connected_quartile_RA <- names(RA_degrees)[order(RA_degrees, decreasing = TRUE)][1:542]
most_connected_quartile_RD <- names(RD_degrees)[order(RD_degrees, decreasing = TRUE)][1:542]
most_connected_quartile_RL <- names(RL_degrees)[order(RL_degrees, decreasing = TRUE)][1:542]

100* length(intersect(most_connected_quartile_RA, most_connected_quartile_RD)) /542 # 61 %
100* length(intersect(most_connected_quartile_RA, most_connected_quartile_RL)) /542 # 77 %
100* length(intersect(most_connected_quartile_RD, most_connected_quartile_RL)) /542 # 62 %

#-------------------------------------------------------------------------------
# Connectivity

# N. isolated nodes
length(isolated_RA) # 360
length(isolated_RD) # 443
length(isolated_RL) # 418

# N. connected components (no isolated nodes)
comp_RA$no - length(isolated_RA) # 43
comp_RD$no - length(isolated_RD) # 43
comp_RL$no - length(isolated_RL) # 37

# Size of lcc
size_lcc_RA # 1670
size_lcc_RD # 1571
size_lcc_RL # 1615
