library(ggplot2)

density_search <- readRDS("..\\data\\processed\\density_search.rds")

data <- data.frame(
    "d" = density_search$densities, 
    "first_comp_RA" = density_search$first_comp$RA/21.70, #To have %: *100/2170
    "first_comp_RD" = density_search$first_comp$RD/21.70,
    "first_comp_RL" = density_search$first_comp$RL/21.70,
    "second_comp_RA" = density_search$second_comp$RA/21.70,
    "second_comp_RD" = density_search$second_comp$RD/21.70,
    "second_comp_RL" = density_search$second_comp$RL/21.70
)

# First component (% of total nodes)
plt1 <- ggplot(data) + geom_line(
    aes(x=d, y=first_comp_RA), color='red'
) + geom_line(
    aes(x=d, y=first_comp_RD), color='green'
) + geom_line(
    aes(x=d, y=first_comp_RL), color='blue'
) + geom_point(
    aes(x=d, y=first_comp_RA), color='red'
) + geom_point(
    aes(x=d, y=first_comp_RD), color='green'
) + geom_point(
    aes(x=d, y=first_comp_RL), color='blue'
)

# Second component (% of total nodes)
plt2 <- ggplot(data) + geom_line(
    aes(x=d, y=second_comp_RA), color='red'
) + geom_line(
    aes(x=d, y=second_comp_RD), color='green'
) + geom_line(
    aes(x=d, y=second_comp_RL), color='blue'
) + geom_point(
    aes(x=d, y=second_comp_RA), color='red'
) + geom_point(
    aes(x=d, y=second_comp_RD), color='green'
) + geom_point(
    aes(x=d, y=second_comp_RL), color='blue'
)
