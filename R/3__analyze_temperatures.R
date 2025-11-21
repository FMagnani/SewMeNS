library(ggplot2)    

#-------------------------------------------------------------------------------

load("..\\data\\processed\\copenhagen_metadata.RData")
RA_order <- order(copenhagen_metadata$RA$collection_date)
RD_order <- order(copenhagen_metadata$RD$collection_date)
RL_order <- order(copenhagen_metadata$RL$collection_date)
RA_dates <- copenhagen_metadata$RA$collection_date[RA_order]
RD_dates <- copenhagen_metadata$RD$collection_date[RD_order]
RL_dates <- copenhagen_metadata$RL$collection_date[RL_order]
RA_abundance <- copenhagen_metadata$RA$sample_abundance[RA_order]
RD_abundance <- copenhagen_metadata$RD$sample_abundance[RD_order]
RL_abundance <- copenhagen_metadata$RL$sample_abundance[RL_order]

dates_df <- rbind(
  data.frame(
    "dates"=RA_dates, "frags"=RA_abundance, "col"="RA", "pos"=-9
  ),
  data.frame(
      "dates"=RD_dates, "frags"=RD_abundance, "col"="RD", "pos"=-10
  ),
  data.frame(
      "dates"=RL_dates, "frags"=RL_abundance, "col"="RL", "pos"=-11
  )
)

#-------------------------------------------------------------------------------

temperatures <- read.csv("..\\data\\raw\\temperatures.csv")

temperatures <- temperatures[!is.na(temperatures$TMAX),]
temperatures <- temperatures[!is.na(temperatures$TMIN),]

temperatures$DATE <- as.Date(temperatures$DATE, "%Y-%m-%d")

start_date <- as.Date("2019.01.01", "%Y.%m.%d")
end_date <- as.Date("2020.10.01", "%Y.%m.%d")

vline <- as.Date("2019.10.01", "%Y.%m.%d")

temperatures$TMEAN <- (temperatures$TMAX+temperatures$TMIN)/2

mean_of_period <- mean((temperatures[temperatures$DATE>vline, ]$TMAX+temperatures[temperatures$DATE>vline, ]$TMIN)/2)
std_of_period <- sd((temperatures[temperatures$DATE>vline, ]$TMAX+temperatures[temperatures$DATE>vline, ]$TMIN)/2)

temperatures$BELOW <- temperatures$TMEAN<mean_of_period

window_1 <- temperatures[
  temperatures$DATE<as.Date("2020.01.01", "%Y.%m.%d") & temperatures$DATE>vline, 
]
window_1_start <- window_1[window_1$BELOW==TRUE, ]$DATE[[1]]
window_1_end <- window_1[window_1$BELOW==FALSE, ]$DATE[[length(window_1[window_1$BELOW==FALSE, ]$DATE)]]

window_2 <- temperatures[
  temperatures$DATE>as.Date("2020.01.01", "%Y.%m.%d") & temperatures$DATE<end_date,
]
window_2_start <- window_2[window_2$BELOW==FALSE, ]$DATE[[1]]
window_2_end <- window_2[window_2$BELOW==TRUE, ]$DATE[[length(window_2[window_2$BELOW==TRUE, ]$DATE)]]

cold_period_start <- mean(c(window_1_start, window_1_end))
cold_period_end <- mean(c(window_2_start, window_2_end))
hot_period_start <- cold_period_end
hot_period_end <- end_date

temperatures$lower_for_blue_shading <- pmin(temperatures$TMIN, mean_of_period)
temperatures$upper_for_blue_shading <- pmin(temperatures$TMAX, mean_of_period)
temperatures$lower_for_red_shading <- pmax(temperatures$TMIN, mean_of_period)
temperatures$upper_for_red_shading <- pmax(temperatures$TMAX, mean_of_period)

#-------------------------------------------------------------------------------

ggplot(
  temperatures[
    (temperatures$DATE>start_date)&(temperatures$DATE<end_date),
    ]
) + geom_vline(
  xintercept = c(
    as.Date("2019.10.01", "%Y.%m.%d"), as.Date("2019.11.01", "%Y.%m.%d"), 
    as.Date("2019.12.01", "%Y.%m.%d"),as.Date("2020.01.01", "%Y.%m.%d"), 
    as.Date("2020.02.01", "%Y.%m.%d"), as.Date("2020.03.01", "%Y.%m.%d"), 
    as.Date("2020.04.01", "%Y.%m.%d"), as.Date("2020.05.01", "%Y.%m.%d"), 
    as.Date("2020.06.01", "%Y.%m.%d"), as.Date("2020.07.01", "%Y.%m.%d"),
    as.Date("2020.08.01", "%Y.%m.%d"), as.Date("2020.09.01", "%Y.%m.%d"), 
    as.Date("2020.10.01", "%Y.%m.%d")
  ), lty='dashed', linewidth=1, alpha=.6
) + annotate(
  "rect", fill = "darkblue", alpha = 0.5, 
  xmin = cold_period_start, xmax = cold_period_end, ymin = -8, ymax =Inf
) + annotate(
  "rect", fill = "darkred", alpha = 0.5, 
  xmin = hot_period_start, xmax = hot_period_end, ymin = -8, ymax = Inf
) + geom_ribbon(
  data=temperatures[(temperatures$DATE>start_date)&(temperatures$DATE<cold_period_start),], 
  aes(x=DATE, ymin=TMIN, ymax=TMAX), fill="black", alpha=0.5
) + geom_ribbon(
  data=temperatures[(temperatures$DATE>vline)&(temperatures$DATE<end_date),], 
  aes(x=DATE, ymin=lower_for_blue_shading, ymax=upper_for_blue_shading), fill="lightblue", alpha=0.7
) + geom_ribbon(
  data=temperatures[(temperatures$DATE>vline)&(temperatures$DATE<end_date),], 
  aes(x=DATE, ymin=lower_for_red_shading, ymax=upper_for_red_shading), fill="red", alpha=0.7
) + geom_line(
  aes(x=DATE, y=TMAX), linewidth=1, alpha=1
) + geom_line(
  aes(x=DATE, y=TMIN), linewidth=1, alpha=1
) + geom_point(
  data = dates_df, 
  aes(x=as.Date(dates, "%Y-%m-%d"), y=pos, fill=col),
  size=7, shape=21, alpha=1, show.legend=TRUE, col="black"
) + labs(
  title="Temperature in Copenhagen and dates of sampling",
  y="Temperature", x="Date"
) + geom_line(
  data = data.frame("x"=c(vline, end_date), "y"=mean_of_period),
  aes(x=x, y=y)
) + geom_text(
  data = data.frame(
    "Date" = c(
      as.Date("2019.10.15", "%Y.%m.%d"), as.Date("2019.11.15", "%Y.%m.%d"), 
      as.Date("2019.12.15", "%Y.%m.%d"), as.Date("2020.01.15", "%Y.%m.%d"), 
      as.Date("2020.02.15", "%Y.%m.%d"), as.Date("2020.03.15", "%Y.%m.%d"), 
      as.Date("2020.04.15", "%Y.%m.%d"), as.Date("2020.05.15", "%Y.%m.%d"), 
      as.Date("2020.06.15", "%Y.%m.%d"), as.Date("2020.07.15", "%Y.%m.%d"),
      as.Date("2020.08.15", "%Y.%m.%d"), as.Date("2020.09.15", "%Y.%m.%d")
    ), 
    "Text" = c(
      "Oct", "Nov", "Dic", "Gen", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
      "Aug", "Sep"
    )  
  ),
  aes(x=Date, y=35, label=Text), fontface="bold"
) + annotate(
  geom="text", x=as.Date("2018.12.28", "%Y.%m.%d"), y=-10, label="Sampling", fontface="bold"
)
