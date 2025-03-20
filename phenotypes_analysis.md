Stabilizing_selection_analysis
================
Siraj
2025-03-19

\#Loading the data I got from Neda from for both the SS and control

``` bash
cd /Users/selmasri/projects/neda_data/data
ls *.csv
```

    ## sieving_count_data.csv
    ## sieving_count_data_controls.csv

``` r
neda_data_ss <- read.csv("/Users/selmasri/projects/neda_data/data/sieving_count_data.csv", header = TRUE, sep = ";")
str(neda_data_ss)
```

    ## 'data.frame':    1230 obs. of  9 variables:
    ##  $ treatment : chr  "SS" "SS" "SS" "SS" ...
    ##  $ replicate : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ sieve_num : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ sieve_size: int  1700 1600 1526 1454 1400 1319 1250 1180 1120 1085 ...
    ##  $ female_num: int  0 0 0 0 2 9 2 5 17 321 ...
    ##  $ male_num  : int  0 0 0 0 3 6 0 1 5 20 ...
    ##  $ who       : chr  "samaneh" "samaneh" "samaneh" "samaneh" ...
    ##  $ date      : chr  "04.11.24" "04.11.24" "04.11.24" "04.11.24" ...
    ##  $ generation: int  1 1 1 1 1 1 1 1 1 1 ...

``` r
neda_data_control <- read.csv("/Users/selmasri/projects/neda_data/data/sieving_count_data_controls.csv", 
                              header = TRUE, sep = ";")
str(neda_data_control)
```

    ## 'data.frame':    150 obs. of  9 variables:
    ##  $ treatment : chr  "SS" "SS" "SS" "SS" ...
    ##  $ replicate : chr  "control1" "control1" "control1" "control1" ...
    ##  $ sieve_num : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ sieve_size: int  1700 1600 1526 1454 1400 1319 1250 1180 1120 1085 ...
    ##  $ female_num: int  0 0 0 0 2 1 5 3 16 203 ...
    ##  $ male_num  : int  0 0 0 0 1 0 3 3 4 48 ...
    ##  $ who       : chr  "samaneh" "samaneh" "samaneh" "samaneh" ...
    ##  $ date      : chr  "07.01.25" "07.01.25" "07.01.25" "07.01.25" ...
    ##  $ generation: int  5 5 5 5 5 5 5 5 5 5 ...

``` r
str(neda_data_ss)
```

    ## 'data.frame':    1230 obs. of  9 variables:
    ##  $ treatment : chr  "SS" "SS" "SS" "SS" ...
    ##  $ replicate : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ sieve_num : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ sieve_size: int  1700 1600 1526 1454 1400 1319 1250 1180 1120 1085 ...
    ##  $ female_num: int  0 0 0 0 2 9 2 5 17 321 ...
    ##  $ male_num  : int  0 0 0 0 3 6 0 1 5 20 ...
    ##  $ who       : chr  "samaneh" "samaneh" "samaneh" "samaneh" ...
    ##  $ date      : chr  "04.11.24" "04.11.24" "04.11.24" "04.11.24" ...
    ##  $ generation: int  1 1 1 1 1 1 1 1 1 1 ...

## function for data cleaning and transformation for both data

``` r
# Convert columns to formats useful for later:

clean_data <- function(data) {
  data$treatment   <- factor(data$treatment)
  data$replicate   <- factor(data$replicate)
  data$sieve_num   <- factor(data$sieve_num)
  data$sieve_size  <- as.numeric(data$sieve_size)
  data$who         <- factor(data$who)
  data$date        <- as.Date(data$date, format = "%d.%m.%y") 
  data$generation  <- factor(data$generation)
  return(data)
}
neda_data_ss <- clean_data(neda_data_ss)
neda_data_control <- clean_data(neda_data_control)
str(neda_data_control)
```

    ## 'data.frame':    150 obs. of  9 variables:
    ##  $ treatment : Factor w/ 1 level "SS": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ replicate : Factor w/ 5 levels "control1","control2",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ sieve_num : Factor w/ 15 levels "1","2","3","4",..: 1 2 3 4 5 6 7 8 9 10 ...
    ##  $ sieve_size: num  1700 1600 1526 1454 1400 ...
    ##  $ female_num: int  0 0 0 0 2 1 5 3 16 203 ...
    ##  $ male_num  : int  0 0 0 0 1 0 3 3 4 48 ...
    ##  $ who       : Factor w/ 2 levels "davide","samaneh": 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ date      : Date, format: "2025-01-07" "2025-01-07" ...
    ##  $ generation: Factor w/ 2 levels "5","10": 1 1 1 1 1 1 1 1 1 1 ...

\#Tranforming female data into pseudo continous by asssinging a row with
the sieve size value for each fly in a replicate

``` r
create_pseudo_continuous_female <- function(data) {
  library(dplyr)
  library(tidyr)
  data %>% 
    uncount(female_num, .remove = FALSE) %>%
    mutate(sieve_expanded = sieve_size) %>%
    select(generation, replicate, sieve_expanded, who)
}
pseudo_continuous_female_ss <- create_pseudo_continuous_female(neda_data_ss)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
pseudo_continuous_female_control <- create_pseudo_continuous_female(neda_data_control)
```

\###Computing pop size for each replicate and including it in the pseudo
continuous, in case needed for later

``` r
pop_size_data_ss <- neda_data_ss %>%
  group_by(generation, replicate) %>%
  summarise(pop_size = sum(female_num, na.rm = TRUE) + sum(male_num, na.rm = TRUE), .groups = "drop")

##including population size in the pseudo continous data for later analysis
pseudo_continuous_female_ss <- pseudo_continuous_female_ss %>%
  left_join(pop_size_data_ss, by = c("generation", "replicate"))
```

\####plot the pop size per generation

``` r
library(ggplot2)
 ggplot(pop_size_data_ss, aes(x = as.factor(generation), y = pop_size, color = replicate, group = replicate)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Population Size by Replicate and Generation",
    x = "Generation",
    y = "Population Size"
  )
```

![](phenotypes_analysis_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

\#Visualisation of replicates across generations in two ways, proportion
and pseudo continuous

``` r
library(Cairo)
library(RColorBrewer)


palette.12 <- brewer.pal(n = 12, name = "Paired")

pseudo_continuous_ss <- ggplot(pseudo_continuous_female_ss, aes(x = factor(generation), y = sieve_expanded, fill = factor(generation))) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(alpha = 0.7, width = 0.1) +
  labs(title = "Violin and Boxplot of Sieve Distributions Across Generations",
       x = "Generation", y = "Sieve Bin (Pseudo-Continuous)", fill = "Generation") +
  scale_fill_manual(values = palette.12) +
  facet_wrap(~ replicate, ncol = 3) +
  theme_minimal()
png(filename = "/Users/selmasri/projects/neda_data/figs/pseudo_continuous_ss.png",
    width = 8, height = 6, units = "in", res = 300, bg = "white")
print(pseudo_continuous_ss)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# ---- Boxplot of Female Proportions Across Generations ----

informative_bins_ss <- neda_data_ss %>%
  group_by(generation, replicate) %>%
  mutate(total_females = sum(female_num, na.rm = TRUE)) %>%  # Compute total females per generation & replicate
  ungroup() %>%
  group_by(generation, replicate, sieve_size) %>%
  summarise(female_prop = sum(female_num, na.rm = TRUE) / unique(total_females), .groups = "drop")



 proportions <- ggplot(informative_bins_ss, aes(x = factor(generation), y = sieve_size, fill = factor(generation))) +
  geom_boxplot(aes(weight = female_prop), alpha = 0.7, outlier.shape = NA) +
  labs(title = "Distribution of Female Proportions Across Generations (Faceted by Replicate)",
       x = "Generation", y = "Sieve Size (µm)", fill = "Generation") +
  scale_fill_manual(values = palette.12) +
  facet_wrap(~ replicate, ncol = 3) +
  scale_y_continuous(breaks = c(937,1085, 1250)) +  # Only key sieve sizes
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"))

png(filename = "/Users/selmasri/projects/neda_data/figs/proportions_ss.png",
    width = 8, height = 6, units = "in", res = 300, bg = "white")
print(proportions)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Computing an average fly size per replicate per generation as a summary statistic based on pseudo-continuous transformation of the data

``` r
###I will change 800 to 850 in the psuedo continuous data so the flies 800 dont drag the mean down although there is rarely females in 800 and then calculate a mean fly size
calc_avg_fly_size <- function(pseudo_data) {
  pseudo_data %>%
    mutate(sieve_expanded = ifelse(sieve_expanded == 800, 850, sieve_expanded)) %>%
    group_by(generation, replicate) %>%
    summarise(avg_fly_size = mean(sieve_expanded, na.rm = TRUE), .groups = "drop")
}

##Apply for SS
avg_fly_size_ss <- calc_avg_fly_size(pseudo_continuous_female_ss)
head(avg_fly_size_ss)
```

    ## # A tibble: 6 × 3
    ##   generation replicate avg_fly_size
    ##   <fct>      <fct>            <dbl>
    ## 1 1          1                1044.
    ## 2 1          2                1114.
    ## 3 1          3                1050.
    ## 4 1          4                1071.
    ## 5 1          5                1039.
    ## 6 1          6                1048.

``` r
#Apply for control
avg_fly_size_control <- calc_avg_fly_size(pseudo_continuous_female_control)
head(avg_fly_size_control)
```

    ## # A tibble: 6 × 3
    ##   generation replicate avg_fly_size
    ##   <fct>      <fct>            <dbl>
    ## 1 5          control1         1030.
    ## 2 5          control2         1028.
    ## 3 5          control3         1007.
    ## 4 5          control4          984.
    ## 5 5          control5         1007.
    ## 6 10         control1          960.

``` r
pseudo_continuous_female_ss_850 <- pseudo_continuous_female_ss %>%
  mutate(sieve_expanded = ifelse(sieve_expanded == 800, 850, sieve_expanded))

avg_fly_size_ss <- pseudo_continuous_female_ss_850 %>%
  group_by(generation, replicate) %>%
  summarise(avg_fly_size = mean(sieve_expanded, na.rm = TRUE), .groups = "drop")  # Compute mean sieve size

head(avg_fly_size_ss)
```

    ## # A tibble: 6 × 3
    ##   generation replicate avg_fly_size
    ##   <fct>      <fct>            <dbl>
    ## 1 1          1                1044.
    ## 2 1          2                1114.
    ## 3 1          3                1050.
    ## 4 1          4                1071.
    ## 5 1          5                1039.
    ## 6 1          6                1048.

\##Visualisation

``` r
avg_fly_plot_ss <- ggplot(avg_fly_size_ss, aes(x = factor(generation), y = avg_fly_size)) +
  geom_boxplot(fill = "orange", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, color = "black") +
  labs(title = "Average Fly Size per Generation", 
       x = "Generation", 
       y = "Mean Fly Size (µm)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),        # Bold, larger title
    axis.title.x = element_text(size = 14, face = "bold"),        # Bold, larger x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),        # Bold, larger y-axis title
    legend.title = element_text(size = 12, face = "bold"),        # Bold legend title
    legend.text = element_text(size = 10)                         # Larger legend text
  )

ggsave(
  filename = "/Users/selmasri/projects/neda_data/figs/mean_size_ss.png",
  plot = avg_fly_plot_ss,
  width = 8,
  height = 6,
  dpi = 300
)


avg_fly_ss_labeled <- ggplot(avg_fly_size_ss, aes(x = factor(generation), y = avg_fly_size)) +
  geom_boxplot(fill = "orange", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, aes(color = replicate)) +
  geom_text(aes(label = replicate), 
            position = position_jitter(width = 0.2, height = 0), 
            vjust = -1.5, size = 3.5) +
  labs(title = "Average Fly Size per Generation",
       x = "Generation", y = "Mean Fly Size (µm)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave(
  filename = "/Users/selmasri/projects/neda_data/figs/mean_size_ss_labels.png",
  plot = avg_fly_ss_labeled,
  width = 8,
  height = 6,
  dpi = 300
)
```

\##Control visualisation

``` r
 library(ggrepel)
mean_size_controls <- ggplot(avg_fly_size_control, aes(x = factor(generation), y = avg_fly_size)) +
  geom_boxplot(fill = "green", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, aes(color = replicate)) +
  # Use geom_text_repel instead of geom_text for non-overlapping labels
  geom_text_repel(aes(label = replicate), size = 3.5) +
  # Add some vertical padding so the boxes and labels aren't jammed at the edges
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(
    title = "Average Fly Size per Generation (Increasing Selection Replicates 90->82%)",
    x = "Generation", 
    y = "Mean Fly Size (µm)"
  ) +
  theme_bw() +
  theme(
    plot.title    = element_text(size = 16, face = "bold"),
    axis.title.x  = element_text(size = 14, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold"),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 10),
    # Optional: add extra margin around the entire plot
    plot.margin   = margin(t = 20, r = 20, b = 20, l = 20)
  )
ggsave(
  filename = "/Users/selmasri/projects/neda_data/figs/mean_size_controls.png",
  plot = mean_size_controls ,
  width = 8,
  height = 6,
  dpi = 300
)
```

\#Anova with generation and replicate as predictor and mean fly size as
response

``` r
#the summarize function above produced a tibble and the anova needs a df

avg_fly_size_ss <- as.data.frame(avg_fly_size_ss)

anova_result_ss <- aov(avg_fly_size ~ generation + replicate, data = avg_fly_size_ss)
summary(anova_result_ss)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## generation   9  58200    6467   6.335 2.44e-06 ***
    ## replicate    9   2042     227   0.222     0.99    
    ## Residuals   63  64308    1021                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
shapiro.test(residuals(anova_result_ss))  # Should be > 0.05 for normality, which is the case
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(anova_result_ss)
    ## W = 0.98739, p-value = 0.6075

``` r
qqnorm(residuals(anova_result_ss))
qqline(residuals(anova_result_ss), col = "red")
```

![](phenotypes_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
residuals_ss <- residuals(anova_result_ss)

library(car)
```

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

``` r
leveneTest(residuals_ss ~ avg_fly_size_ss$generation) #p> 0.05 assumption met of homogenity of variance.
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  9  1.4314  0.191
    ##       72

``` r
TukeyHSD(anova_result_ss)
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = avg_fly_size ~ generation + replicate, data = avg_fly_size_ss)
    ## 
    ## $generation
    ##            diff          lwr        upr     p adj
    ## 2-1  -15.613678  -62.4725427  31.245187 0.9837199
    ## 3-1  -34.045854  -80.9047186  12.813011 0.3534701
    ## 4-1  -17.867413  -71.9753695  36.240544 0.9846930
    ## 5-1  -42.622204 -104.6106554  19.366248 0.4319442
    ## 6-1  -90.933551 -172.0954858  -9.771616 0.0165316
    ## 7-1  -58.279912 -105.1387767 -11.421047 0.0047375
    ## 8-1  -82.857331 -129.7161962 -35.998466 0.0000100
    ## 9-1  -30.676650  -77.5355150  16.182215 0.5028033
    ## 10-1 -11.545366  -58.4042312  35.313499 0.9982140
    ## 3-2  -18.432176  -65.2910408  28.426689 0.9523089
    ## 4-2   -2.253735  -56.3616918  51.854221 1.0000000
    ## 5-2  -27.008526  -88.9969777  34.979926 0.9133317
    ## 6-2  -75.319873 -156.4818080   5.842062 0.0909049
    ## 7-2  -42.666234  -89.5250990   4.192631 0.1042707
    ## 8-2  -67.243654 -114.1025185 -20.384789 0.0005694
    ## 9-2  -15.062972  -61.9218373  31.795893 0.9872885
    ## 10-2   4.068312  -42.7905534  50.927177 0.9999997
    ## 4-3   16.178441  -37.9295159  70.286397 0.9923889
    ## 5-3   -8.576350  -70.5648018  53.412102 0.9999846
    ## 6-3  -56.887697 -138.0496322  24.274238 0.4043173
    ## 7-3  -24.234058  -71.0929231  22.624807 0.7934168
    ## 8-3  -48.811478  -95.6703426  -1.952613 0.0345755
    ## 9-3    3.369204  -43.4896614  50.228069 1.0000000
    ## 10-3  22.500487  -24.3583776  69.359352 0.8550329
    ## 5-4  -24.754791  -92.3897366  42.880155 0.9696505
    ## 6-4  -73.066138 -158.6183291  12.486053 0.1586082
    ## 7-4  -40.412499  -94.5204555  13.695458 0.3156481
    ## 8-4  -64.989918 -119.0978750 -10.881962 0.0073827
    ## 9-4  -12.809237  -66.9171938  41.298719 0.9986917
    ## 10-4   6.322047  -47.7859099  60.430003 0.9999964
    ## 6-5  -48.311347 -139.0531490  42.430455 0.7650063
    ## 7-5  -15.657708  -77.6461598  46.330744 0.9978362
    ## 8-5  -40.235128 -102.2235793  21.753324 0.5150969
    ## 9-5   11.945554  -50.0428981  73.934005 0.9997496
    ## 10-5  31.076837  -30.9116143  93.065289 0.8211535
    ## 7-6   32.653639  -48.5082958 113.815574 0.9453412
    ## 8-6    8.076220  -73.0857153  89.238155 0.9999991
    ## 9-6   60.256901  -20.9050341 141.418836 0.3237426
    ## 10-6  79.388185   -1.7737503 160.550120 0.0602792
    ## 8-7  -24.577420  -71.4362845  22.281445 0.7799514
    ## 9-7   27.603262  -19.2556033  74.462127 0.6479197
    ## 10-7  46.734546   -0.1243195  93.593410 0.0511603
    ## 9-8   52.180681    5.3218162  99.039546 0.0176588
    ## 10-8  71.311965   24.4531001 118.170830 0.0002056
    ## 10-9  19.131284  -27.7275811  65.990149 0.9403508
    ## 
    ## $replicate
    ##             diff       lwr      upr     p adj
    ## 2-1    6.0308460 -44.88289 56.94458 0.9999960
    ## 3-1   -6.3878139 -55.78139 43.00577 0.9999914
    ## 4-1   -0.0253061 -50.93904 50.88843 1.0000000
    ## 5-1   -3.2448388 -54.15858 47.66890 1.0000000
    ## 6-1   -0.2078449 -51.12158 50.70589 1.0000000
    ## 7-1    2.5714839 -48.34225 53.48522 1.0000000
    ## 8-1  -11.9529188 -62.86666 38.96082 0.9987746
    ## 9-1   -7.5384351 -58.45217 43.37530 0.9999725
    ## 10-1  -2.5277678 -53.44151 48.38597 1.0000000
    ## 3-2  -12.4186599 -63.33240 38.49508 0.9983482
    ## 4-2   -6.0561521 -58.44596 46.33365 0.9999967
    ## 5-2   -9.2756849 -61.66549 43.11412 0.9998757
    ## 6-2   -6.2386910 -58.62849 46.15111 0.9999958
    ## 7-2   -3.4593622 -55.84917 48.93044 1.0000000
    ## 8-2  -17.9837649 -70.37357 34.40604 0.9800904
    ## 9-2  -13.5692811 -65.95908 38.82052 0.9973784
    ## 10-2  -8.5586138 -60.94842 43.83119 0.9999367
    ## 4-3    6.3625078 -44.55123 57.27625 0.9999936
    ## 5-3    3.1429750 -47.77076 54.05671 1.0000000
    ## 6-3    6.1799689 -44.73377 57.09371 0.9999950
    ## 7-3    8.9592978 -41.95444 59.87304 0.9998819
    ## 8-3   -5.5651049 -56.47884 45.34863 0.9999980
    ## 9-3   -1.1506212 -52.06436 49.76312 1.0000000
    ## 10-3   3.8600461 -47.05369 54.77378 0.9999999
    ## 5-4   -3.2195327 -55.60934 49.17027 1.0000000
    ## 6-4   -0.1825388 -52.57234 52.20726 1.0000000
    ## 7-4    2.5967900 -49.79301 54.98659 1.0000000
    ## 8-4  -11.9276127 -64.31742 40.46219 0.9990380
    ## 9-4   -7.5131290 -59.90293 44.87667 0.9999791
    ## 10-4  -2.5024617 -54.89227 49.88734 1.0000000
    ## 6-5    3.0369939 -49.35281 55.42680 1.0000000
    ## 7-5    5.8163227 -46.57348 58.20613 0.9999977
    ## 8-5   -8.7080800 -61.09788 43.68172 0.9999268
    ## 9-5   -4.2935962 -56.68340 48.09621 0.9999998
    ## 10-5   0.7170710 -51.67273 53.10687 1.0000000
    ## 7-6    2.7793288 -49.61047 55.16913 1.0000000
    ## 8-6  -11.7450739 -64.13488 40.64473 0.9991488
    ## 9-6   -7.3305901 -59.72039 45.05921 0.9999830
    ## 10-6  -2.3199229 -54.70973 50.06988 1.0000000
    ## 8-7  -14.5244027 -66.91421 37.86540 0.9956243
    ## 9-7  -10.1099189 -62.49972 42.27988 0.9997467
    ## 10-7  -5.0992517 -57.48906 47.29055 0.9999993
    ## 9-8    4.4144838 -47.97532 56.80429 0.9999998
    ## 10-8   9.4251510 -42.96465 61.81495 0.9998581
    ## 10-9   5.0106673 -47.37914 57.40047 0.9999994
