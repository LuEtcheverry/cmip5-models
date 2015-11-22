# Comparing the CMIP5 suite of climate models to observed global temperatures.

Let's start by readying the necessary packages that we'll be using to analyse the data.



```r
rm(list=ls()) ## Clear data

require(readr) ## For reading in data
require(tidyr) ## For data tidying (gather, etc.)
require(dplyr) ## For data munging and manipulation (filter, mutate, etc.)
require(quantreg) ## For efficient recursive regression estimates
require(Hmisc) ## For additional minor ticks in plot
require(alr3) ## For calculating confidence intervals for regression estimates
```

Next, we read in the data. Note that I have already uploaded the CMIP5 ensemble data to this GitHub repository, since that saves you from having to go through the [Climate Explorer](http://climexp.knmi.nl/) interface.[^1] The HadCRUT4 data, on the other hand, is already in a convenient form and so we'll download this data directly from the Met Office website. 


```r
### CMIP5 ensemble ###
cmip5 <- read_csv("cmip5.csv")

### HadCRUT4 ###          
had <- read_table("http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.4.0.0.annual_ns_avg.txt",
                         col_names = c("Year", "HadCRUT4", 3:12)
                         )
```

Now we need to tidy the data and get the annual temperature means. We'll be using the extremely useful `tidyr` and `dplyr` packages from [Hadley](http://had.co.nz/), and also be making liberal use of the awesome pipe operator (`%>%`).


```r
## First decide on dates for the recursive sample.
end_year <- 2014 ## M&K use 2014. This is the fixed point for the recursive estimates.
start_year <- end_year - 63 ## To coincide with's M&K's choice of 1951 when using 2014 as an end date.

### CMIP5 Ensemble ###
## Already done. See: https://github.com/grantmcdermott/cmip5-models/blob/master/cmip5_download/cmip5_download.md

### HadCRUT4 ###
had <- 
  had %>%
  select(Year, HadCRUT4) %>%
  filter(Year >= cmip5$Year[1] & Year<=end_year) %>%
  group_by(Year) %>%
  summarise(HadCRUT4 = mean(HadCRUT4)) %>% ## Collapse into all annual averages
  gather(Series, Temp, HadCRUT4)
```

Using the annual data, we can now calculate the trends and recursive estimates for both the CMIP5 ensemble and the observed (i.e. HadCRUT4) temperature data. Note that the CMIP5 recursive regressions are grouped by model and that we multiply the end result by 10 to obtain the decadal trend, so as to be consistent with M&K's paper.


```r
### CMIP5 ensemble ###
cmip5_trend <- 
  cmip5 %>%
  group_by(Model) %>% 
  arrange(desc(Year)) %>%
  mutate(Length = max(Year) - Year + 1) %>% 
  mutate(Time = Year - min(Year) + 1) %>% 
  filter(Year >= start_year) %>%
  mutate(Trend = 10 * as.numeric(lm.fit.recursive(as.matrix(Time), Temp, int = T)[2, ]))

## Summarise the indiv. model results i.t.o. ensemble mean and C.I.s.
cmip5_trend_summ <- 
  cmip5_trend %>% 
  group_by(Year) %>%
  summarise(Trend.mean = mean(Trend),
            Trend_025 = quantile(Trend, p = 0.025),
            Trend_05 = quantile(Trend, p = 0.05),
            Trend_95 = quantile(Trend, p = 0.95),
            Trend_975 = quantile(Trend, p = 0.975)
            ) %>%
  mutate(Length = max(Year) - Year + 1) %>%
  filter(Length >= 10)

### HadCRUT4 ###
had_trend <- 
  had %>%
  arrange(desc(Year)) %>%
  mutate(Length = max(Year) - Year + 1) %>% 
  mutate(Time = Year - min(Year) + 1) %>% 
  filter(Year >= start_year) %>%
  mutate(Trend = 10 * as.numeric(lm.fit.recursive(as.matrix(Time), Temp, int = T)[2, ]) ) %>%
  filter(Length >= 10)
```

Nearly there. We now combine our recursive trend estimates in a single figure that matches [**this one**](http://object.cato.org/sites/cato.org/files/wp-content/uploads/gsr_12_18_14_fig2.png) from M&K. I'd normally do this in `ggplot2`, but I'm using R's base plotting features (with a number of little tweaks) to try and get it looking as close to M&K's original as possible.


```r
par(las = 1) # Rotate y-axis labels sideways
plot(cmip5_trend_summ$Length, cmip5_trend_summ$Trend.mean,
     type="l", lwd = 1.5,
     main = "", 
     xlab = paste("Beginning Year of Trend Length (all trends end in ", end_year, ")", sep = ""),
     ylab = "Trend (°C/decade)",
     xlim = rev(range(cmip5_trend_summ$Length)),
     ylim = c(-.1, .7),
     cex = 0.8
     )
minor.tick(nx = 10, ny = 10) ## Optional
points(cmip5_trend_summ$Length, cmip5_trend_summ$Trend.mean, pch = 20)
lines(cmip5_trend_summ$Length, cmip5_trend_summ$Trend_025, lwd = 1, lty = 3)
lines(cmip5_trend_summ$Length, cmip5_trend_summ$Trend_05, lwd = 1, lty = 1)
lines(cmip5_trend_summ$Length, cmip5_trend_summ$Trend_95, lwd = 1, lty = 1)
lines(cmip5_trend_summ$Length, cmip5_trend_summ$Trend_975, lwd = 1, lty = 3)
lines(had_trend$Length, had_trend$Trend, col = "limegreen", lwd = 1.5, lty = 1)
points(had_trend$Length, had_trend$Trend, col = "limegreen", pch = 20)
legend("topleft", 
       c("Black circles = Multi-model mean trend", 
         "Thin Black Lines = 5th and 95th percentiles of model trends", 
         "Dotted Black Lines = 2.5th and 97.5th percentiles of model trends",
         "", 
         "Green Circles = Observed trend consistent with modeled trend", 
         "Yellow Circles = Observed trend at or below the 5th percentile of modeled trends", 
         "Red Circles = Observed trend at or below the 2.5th percentile of modeled trends"
         ), 
       text.col = c("black", "black", "black", "black", "limegreen", "gold", "red"),
       bty = "n", cex = 0.8)
## Some additional lines of code to match colour coding of M&K.
had_05 <- had_trend$Trend < cmip5_trend_summ$Trend_05
had_025 <- had_trend$Trend < cmip5_trend_summ$Trend_025
had.col <- ifelse(had_025, "red", ifelse(had_05, "gold", "limegreen"))
had_05_line <- ifelse(had_05, had_trend$Trend, NA)
had_025_line <- ifelse(had_025, had_trend$Trend, NA)
lines(had_trend$Length, had_05_line, col = "gold", lwd=1.5, lty = 1)
lines(had_trend$Length, had_025_line, col = "red", lwd = 1.5, lty = 1)
points(had_trend$Length, had_trend$Trend, col = had.col, pch = 20)
```

![](modelsVSobs_files/figure-html/unnamed-chunk-5-1.png) 

Good news: We have successfully replicated M&K's [key figure](http://object.cato.org/sites/cato.org/files/wp-content/uploads/gsr_12_18_14_fig2.png)! Now comes the interesting part, however, because M&K use this figure [to claim](http://www.cato.org/blog/agu-2014-quantifying-lack-consistency-between-climate-model-projections-observations-evolution) that: *"[A]t the global scale, this suite of climate models has failed. Treating them as mathematical hypotheses, which they are, means that it is the duty of scientists to reject their predictions in lieu of those with a lower climate sensitivity."* 

As, we're about to see, neither claim makes much sense at all because M&K haven't actually done any proper hypothesis testing... Indeed, they have completely forgotten to calculate confidence intervals!

This is a really important point. The "spread" that we see in the above graph is just a pseudo-interval generated by the full CMIP5 ensemble (i.e. the individual trend means across all 108 climate models). In other words, it is a measure of *model uncertainty*, **not** *statistical uncertainty*. We haven't said anything yet about the confidence intervals attached to each trend estimate, $\hat{\beta_1}$. And, as any first-year econometrics student will tell you, regression coefficients must come with standard errors and implied confidence intervals. Indeed, accounting for these uncertainties is primarily what hypothesis testing is all about: Can we confidently rule out that a parameter of interest doesn't overlap with some value or range? Absent these uncertainty measures, one cannot talk meaningfully about rejecting a hypothesis. 

With these points in mind, let us add the 95% error bars to each of the "observed" HadCRUT4 trend estimates. First we, calculate we need to calculate them. A somewhat technical point is that the `lm.fit.recursive` function from the `quantreg` package (which we used above to get the recursive regression coefficient means) doesn't actually produce standard errors or confidence intervals. Nonetheless, we are able to obtain these manually by calling a loop of sequential (i.e recursive) regressions using the `do()` function in combination with various `apply` functions to extract the relevant parameters. We could also have used a standard "for loop"", but this way is much more efficient.



```r
had_rec <- 
  had %>%
  arrange(desc(Year)) %>%
  mutate(Length = max(Year) - Year + 1) %>% 
  mutate(Time = Year - min(Year) + 1) %>%
  filter(Year >= start_year) %>%
  do(mod = lapply(seq(10, nrow(.)), function(x) lm(Temp ~ Time, data = .[1:x, ])))
  
## Extract recursive trend (already calculated earlier with the "lm.fit.recursive" call) and 95% C.I.
## Note: Again, we multiply by ten to match the decadal time scale.
trend_rec <- 10*unlist(sapply(1:nrow(had_trend), function(j) had_rec$mod[[1]][[j]]$coefficients["Time"])) ## Already calculated, but just illustrating the point.
ci_rec <- 10*unlist(sapply( 1:nrow(had_trend), function(j) confint(had_rec$mod[[1]][[j]], level = .95)["Time",]))
```

Now we add these confidence intervals to the previous plot.


```r
## Add to plot
arrows(10:max(had_trend$Length), ci_rec[1,], 10:max(had_trend$Length), ci_rec[2, ], length = 0.05, angle = 90, code=3)
text(62, -0.08, "Note: 95% confidence intervals added to observed trends.",
     adj = c(0,0), cex = 0.85)
```

![](modelsVSobs_files/figure-html/unnamed-chunk-7-1.png) 

As we can see from the updated graph, the error bars comfortably overlap the model ensemble at every point in time. M&K's bold assertion that we need to reject the climate models as "failed hypotheses" simply does not hold water.

[^1]: If you *really* can't help yourself, click [here](https://github.com/grantmcdermott/cmip5-models/blob/master/cmip5_download/cmip5_download.md) for instructions on how to download and clean the CMIP5 ensemble data.
