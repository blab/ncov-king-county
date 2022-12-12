library(tidyverse)
dates <- read_csv("~/Desktop/gitrepos/ncov-king-county/analysis/data-files//dates.csv")

immunity <- data.frame(date = as.Date(c("2020-07-20","2020-08-24","2020-09-14", "2020-10-10", "2020-11-07", "2020-12-15", "2021-01-13", "2021-02-17","2021-03-18", "2021-04-14", "2021-05-18", "2021-06-17", "2021-07-15", "2021-08-13", "2021-09-18", "2021-10-10", "2021-11-17", "2021-12-12")), 
                       imm_perc = c(0.018, 0.018, 0.018, 0.018, 0.02, 0.045, 0.083, 0.229, 0.38, 0.637, 0.836, 0.93, 0.941,0.937, 0.945, 0.965,0.973,0.966 ))

#need to find dates and then interpolate based on those dates
approxData <- data.frame(
  with(immunity, 
       approx(date, imm_perc, xout = dates$Date, method = "linear", rule = 2)
  ),
  method = "approx()"
)


splineData <- data.frame(
  with(immunity, 
       spline(date, imm_perc, xout = dates$Date)
  ),
  method = "spline()"
)


spline_dup <-splineData[rep(seq_len(nrow(splineData)), each = 2), ]
spline_sub <- spline_dup %>% select(y)
# calculating reverse
rev_data_frame <- apply(spline_sub, 2, rev)

# converting the result to dataframe
rev_data_frame <- as.data.frame(rev_data_frame)
write.table(rev_data_frame, "~/Desktop/glm_multi/immunity.tsv", sep = " ",eol = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

plot(splineData)

library(ggplot2)
ggplot(rbind(approxData, splineData), aes(x = as.numeric(date ), y =imm_perc)) + 
  geom_point(immunity = immunity, aes(as.numeric(date), imm_perc), alpha = 0.2, col = "red") +
  geom_line(col = "blue") +
  facet_wrap(~method) +
  ggtitle("Interpolation and smoothing functions in R") +
  theme_bw(16)

f <- with(immunity, approxfun(date, imm_perc, rule = 2))
curve(f(x), as.Date("2020-01-20"), as.Date("2021-12-12"), col = "green2")
with(immunity, points(date, imm_perc))
is.function(fc <- with(immunity, approxfun(date, imm_perc, method = "const"))) # TRUE
curve(fc(x), as.Date("2020-01-20"), as.Date("2021-12-12"), col = "darkblue", add = TRUE)
## different extrapolation on left and right side :
plot(with(immunity, approxfun(date, imm_perc, rule = 2),
     col = "tomato", add = TRUE, lty = 3, lwd = 2))

