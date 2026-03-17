library(mclust)

load("..\\data\\processed\\communities_res_0.RData")
adjustedRandIndex(communities$RA, communities$RD) # 0.31
adjustedRandIndex(communities$RA, communities$RL) # 0.37
adjustedRandIndex(communities$RD, communities$RL) # 0.30

# Name of variable: communities
load("..\\data\\processed\\seasonalCommunities_res_0.RData")

adjustedRandIndex(communities$Hot_RA, communities$Hot_RD) # 0.20
adjustedRandIndex(communities$Hot_RA, communities$Hot_RL) # 0.23
adjustedRandIndex(communities$Hot_RD, communities$Hot_RL) # 0.26

adjustedRandIndex(communities$Cold_RA, communities$Cold_RD) # 0.25
adjustedRandIndex(communities$Cold_RA, communities$Cold_RL) # 0.31
adjustedRandIndex(communities$Cold_RD, communities$Cold_RL) # 0.26

adjustedRandIndex(communities$Hot_RA, communities$Cold_RA) # 0.14
adjustedRandIndex(communities$Hot_RA, communities$Cold_RD) # 0.16
adjustedRandIndex(communities$Hot_RD, communities$Cold_RL) # 0.19