1_data_prep
================
Siraj
2025-03-13

# Phenotypic data preparation

We start with our phenotypic data stored in the following csv files:

``` bash
cd /Users/selmasri/projects/trunc_sel_ee/data/phenotypic/raw
ls *.csv
```

    ## E&R_Phenotypic_Data_D.sim - phenotype_data.csv
    ## Thorax_measurments - Thorax_measurments_unbiased.csv
    ## Thorax_measurments - Variance_estimation.csv
    ## Thorax_measurments_Variance_estimation2days.csv

Next, we process this into

``` r
library(readr)
library(dplyr)
library(ggsignif)
library(stringr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(lme4)
library(emmeans)
library(RColorBrewer)
library(MuMIn)

# ---- Load Data ----
datap <- read.csv("/Users/selmasri/projects/trunc_sel_ee/data/phenotypic/raw/E&R_Phenotypic_Data_D.sim\ -\ phenotype_data.csv", header = T, sep = ",")
str(datap)
```

    ## 'data.frame':    1885 obs. of  12 variables:
    ##  $ generation      : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ regime          : chr  "constant" "constant" "constant" "constant" ...
    ##  $ date            : chr  "04/11/2024" "04/11/2024" "04/11/2024" "04/11/2024" ...
    ##  $ replicate       : chr  "K1" "K1" "K1" "K1" ...
    ##  $ sieve_num       : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ sieve_size      : int  1700 1600 1526 1454 1400 1319 1250 1180 1120 1085 ...
    ##  $ female_num      : chr  NA NA "0" "4" ...
    ##  $ male_num        : int  NA NA 0 2 0 4 3 3 36 202 ...
    ##  $ dead            : chr  "" "" "" "1" ...
    ##  $ sexing_counting : chr  "misa_siraj" "misa_siraj" "misa_siraj" "misa_siraj" ...
    ##  $ siever          : chr  "Siraj" "Siraj" "Siraj" "Siraj" ...
    ##  $ selected_females: chr  "0" "0" "0" "4" ...

``` r
# ---- Data Preprocessing ----
# Replace NA in `selected_females` with 0 for rows where replicate == "control"
datap <- datap %>%
  mutate(
    # For "control" selected_females is NA, I set it to 0 to avoid numeric conversion issues
    selected_females = if_else(
      replicate == "control" & is.na(selected_females),
      "0",  # use "0" here so it stays a string for now
      selected_females
    )
  ) %>%
  # convert columns to numeric and chr to factors.
  mutate(
    across(
      c(female_num, male_num, selected_females),
      ~ as.numeric(trimws(.))
    ),
    date     = dmy(date),
    siever   = as.factor(str_to_title(siever)),
    generation   = as.factor(generation),
    replicate = as.factor(replicate)
  ) %>%
  filter(!is.na(female_num) & !is.na(male_num)) %>%  ##filter rows in female_num and male_num where there is NA
  dplyr :: select(-dead) #remove uneeded rows


# Check if transformation worked
str(datap)
```

    ## 'data.frame':    1142 obs. of  11 variables:
    ##  $ generation      : Factor w/ 9 levels "1","2","3","4",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ regime          : chr  "constant" "constant" "constant" "constant" ...
    ##  $ date            : Date, format: "2024-11-04" "2024-11-04" ...
    ##  $ replicate       : Factor w/ 21 levels "C10","C6","C7",..: 6 6 6 6 6 6 6 6 6 6 ...
    ##  $ sieve_num       : num  3 4 5 6 7 8 9 10 11 12 ...
    ##  $ sieve_size      : int  1526 1454 1400 1319 1250 1180 1120 1085 1033 1000 ...
    ##  $ female_num      : num  0 4 4 26 25 143 201 381 203 29 ...
    ##  $ male_num        : num  0 2 0 4 3 3 36 202 299 132 ...
    ##  $ sexing_counting : chr  "misa_siraj" "misa_siraj" "misa_siraj" "misa_siraj" ...
    ##  $ siever          : Factor w/ 3 levels "Kati","Misa",..: 3 3 3 3 3 3 3 3 3 3 ...
    ##  $ selected_females: num  0 4 4 26 25 143 16 0 0 0 ...

``` r
datap$sieve_num <- factor(datap$sieve_num)
sum(is.na(datap$date))  # Checking missing dates
```

    ## [1] 0

Next, we transform the data into pseudo-continuous:

``` r
# ---- Transforming Data into Pseudo-Continuous ----
pseudo_continuous_female <- datap %>%
  filter(!is.na(female_num) & female_num > 0) %>%  # Remove NAs and zero counts
  uncount(female_num, .remove = FALSE) %>%  # Expand rows based on fly count
  mutate(sieve_expanded = sieve_size) %>%  # Assign sieve size
  dplyr :: select(generation, replicate, sieve_expanded, siever)  # Ensure siever is included

# Filtering Optimized Replicates
pseudo_continuous_female_O <- pseudo_continuous_female %>%
  filter(grepl("^O", replicate))  # Select only optimized replicates

# Merging 9.5 Bin into 10, to compare after F6 to before F6 but this approach should be changed for later generations
pseudo_continuous_female_O <- pseudo_continuous_female_O %>%
  mutate(sieve_expanded = ifelse(sieve_expanded == 9.5, 10, sieve_expanded))
```

Then we can compute the population size for each replicate:

``` r
# ---- Computing Population Size ----
pop_size_data <- datap %>%
  group_by(generation, replicate) %>%
  summarise(pop_size = sum(female_num, na.rm = TRUE) + sum(male_num, na.rm = TRUE), .groups = "drop")

##plotting population size filtering for less than 200 because this is the incomplete data from K 

pop_size_data_filtered <- pop_size_data %>%
  filter(pop_size >= 200)

##including population size in the pseudo continous data for later analysis
pseudo_continuous_female <- pseudo_continuous_female %>%
  left_join(pop_size_data, by = c("generation", "replicate"))
```

## Population size visualisation

``` r
# 2. Plot pop_size vs. generation, colored by replicate

#replicate = mutate()
pop_size <- ggplot(pop_size_data_filtered, aes(x = as.factor(generation), y = pop_size, color = replicate, group = replicate)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Population Size by Replicate and Generation",
    x = "Generation",
    y = "Population Size"
  )

ggsave("figs/pop_size_across_time.png", pop_size, width = 12, height = 8, dpi = 400)
ggsave("figs/pop_size_across_time.svg", pop_size, width = 12, height = 8)

knitr::include_graphics("figs/pop_size_across_time.png")
```

<img src="figs/pop_size_across_time.png" width="4800" />

## Visualisation of the increasing regimes replicates

``` r
# ---- Visualization: Violin and Boxplots Faceted by Replicate ----
palette.12 <- brewer.pal(n = 12, name = "Paired")

ggplot(pseudo_continuous_female_O, aes(x = factor(generation), y = sieve_expanded, fill = factor(generation))) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(alpha = 0.7, width = 0.1) +
  labs(title = "Violin and Boxplot of Sieve Distributions Across Generations",
       x = "Generation", y = "Sieve Bin (Pseudo-Continuous)", fill = "Generation") +
  scale_fill_manual(values = palette.12) +
  facet_wrap(~ replicate, ncol = 3) +
  theme_minimal()
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# ---- Boxplot of Female Proportions Across Generations ----

informative_bins <- datap %>%
  filter(replicate %in% paste0("O", 1:8)) %>%  # Keep only optimized replicates
  group_by(generation, replicate) %>%
  mutate(total_females = sum(female_num, na.rm = TRUE)) %>%  # Compute total females per generation & replicate
  ungroup() %>%
  group_by(generation, replicate, sieve_size) %>%
  summarise(female_prop = sum(female_num, na.rm = TRUE) / unique(total_females), .groups = "drop")


# Merge 1102 into 1085
informative_bins <- informative_bins %>% 
  mutate(sieve_size = ifelse(sieve_size == 1102, 1085, sieve_size)) %>%
  group_by(replicate, generation, sieve_size) %>%
  summarize(female_prop = mean(female_prop, na.rm = TRUE), .groups = "drop")


poportion_O_repl <- ggplot(informative_bins, aes(x = factor(generation), y = sieve_size, fill = factor(generation))) +
  geom_boxplot(aes(weight = female_prop), alpha = 0.7, outlier.shape = NA) +
  labs(title = "Distribution of Female Proportions Across Generations (Faceted by Replicate)",
       x = "Generation", y = "Sieve Size (µm)", fill = "Generation") +
  scale_fill_manual(values = palette.12) +
  facet_wrap(~ replicate, ncol = 3) +
  scale_y_continuous(breaks = c(937,1085, 1250)) +  # Only key sieve sizes
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"))
poportion_O_repl
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggsave(filename = "/Users/selmasri/projects/trunc_sel_ee/figs/proportion_repl_O.png",
       plot = poportion_O_repl,
       width = 20,    # less wide, adjust as needed
       height = 10,   # adjust height as desired
       units = "in")
```

## Calculations and Visualisation of constant and control replicates proportions

``` r
constant_data_prop <- datap %>%
  filter((replicate == "K7" & generation %in% c(1, 2, 7)) | 
           (replicate == "K1" & generation %in% c(1, 2, 8))) %>%
  mutate(sieve_size = ifelse(sieve_size == 1102, 1085, sieve_size)) %>%  # Merge 1102 into 1085
  group_by(replicate, generation, sieve_size) %>%
  summarise(female_count = sum(female_num, na.rm = TRUE), .groups = "drop") %>%
  group_by(replicate, generation) %>%
  mutate(female_prop = female_count / sum(female_count))  # Compute proportion

head(constant_data_prop)
```

    ## # A tibble: 6 × 5
    ## # Groups:   replicate, generation [1]
    ##   replicate generation sieve_size female_count female_prop
    ##   <fct>     <fct>           <dbl>        <dbl>       <dbl>
    ## 1 K1        1                 800            0     0      
    ## 2 K1        1                 900            3     0.00289
    ## 3 K1        1                 937           18     0.0174 
    ## 4 K1        1                1000           29     0.0280 
    ## 5 K1        1                1033          203     0.196  
    ## 6 K1        1                1085          381     0.367

``` r
control_data <- datap %>%
  filter((replicate == "C7" & generation %in% c(2, 8)) | 
           (replicate == "C10" & generation %in% c(2))) %>% ##i can add another generation when i sex C10
  mutate(sieve_size = ifelse(sieve_size == 1102, 1085, sieve_size)) %>%  # Merge 1102 into 1085
  group_by(replicate, generation, sieve_size) %>%
  summarise(female_count = sum(female_num, na.rm = TRUE), .groups = "drop") %>%
  group_by(replicate, generation) %>%
  mutate(female_prop = female_count / sum(female_count))
head(control_data)
```

    ## # A tibble: 6 × 5
    ## # Groups:   replicate, generation [1]
    ##   replicate generation sieve_size female_count female_prop
    ##   <fct>     <fct>           <dbl>        <dbl>       <dbl>
    ## 1 C10       2                 800            0     0      
    ## 2 C10       2                 900            1     0.00108
    ## 3 C10       2                 937           97     0.105  
    ## 4 C10       2                1000          128     0.138  
    ## 5 C10       2                1033          440     0.474  
    ## 6 C10       2                1085          245     0.264

``` r
ggplot(constant_data_prop %>% filter(replicate == "K1"),
       aes(x = factor(generation), y = sieve_size, fill = factor(generation))) +
  geom_violin(alpha = 0.5, trim = FALSE) +  # Add violin plot with semi-transparency
  geom_boxplot(aes(weight = female_prop), width = 0.2, alpha = 0.7, outlier.shape = NA) +
  labs(title = "Comparison of Female Distributions in K1 (Gen 1, 2, 8)",
       x = "Generation", y = "Sieve Size (µm)", fill = "Generation") +
  scale_fill_manual(values = c("red", "blue", "purple")) +
  scale_y_continuous(breaks = c(937,1033,1085,1400)) +
  theme_bw()
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
##plot the proportions in different generations in some replicates after filtering
constant_data_prop_K1 <- constant_data_prop %>%
  filter(replicate == "K1", generation %in% c(1, 2, 8)) %>%
  mutate(type = "Constant")

control_data_C7 <- control_data %>%
  filter(replicate == "C7", generation %in% c(2, 8)) %>%
  mutate(type = "Control")

informative_bins_O1 <- informative_bins %>%
  filter(replicate == "O1", generation %in% c(1, 2, 5, 8)) %>%
  mutate(type = "Optimized")

# Combine the three datasets into one data frame
combined_data <- bind_rows(constant_data_prop_K1, control_data_C7, informative_bins_O1)
facet_labeller <- as_labeller(c("Constant" = "Constant", "Control" = "Control", "Optimized" = "Increasing"))
# Create the faceted plot
plot_save <- ggplot(combined_data, aes(x = sieve_size, y = female_prop, color = factor(generation))) +
  geom_point() +
  geom_line(aes(group = generation)) +
  facet_wrap(~ type, ncol = 1, labeller = facet_labeller) +
  labs(title = "Female Proportions Across Sieves by Generation",
       x = "Sieve Size (µm)",
       y = "Female Proportion",
       color = "Generation") +
  scale_x_continuous(breaks = c(937, 1085, 1400)) +
  theme_bw() +
  theme(strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_blank(),           # Remove the title
        legend.title = element_text(size = 14),   # Increase legend title font size
        legend.text = element_text(size = 12))      # Increase legend text font size
print(plot_save)
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
ggsave(filename = "figs/my_plot_K1_O1_C7.png",
       plot = plot_save,
       width = 20,    # less wide, adjust as needed
       height = 10,   # adjust height as desired
       units = "in")
ggsave(filename = "figs/my_plot_K1_O1_C7.svg",
       plot = plot_save,
       width = 20,    # less wide, adjust as needed
       height = 10,   # adjust height as desired
       units = "in")
```

# Computing an average ly size per replicate per generation as a summary statistic based on pseudo-continuous transformation of the data

``` r
###I will change 800 to 850 in the psuedo continuous data so the flies 800 dont drag the mean down although there is rarely females in 800
pseudo_continuous_female_O_850 <- pseudo_continuous_female_O %>%
  mutate(sieve_expanded = ifelse(sieve_expanded == 800, 850, sieve_expanded))

gen9_data <- pseudo_continuous_female_O %>%
  filter(generation == 9)

# Compute the average fly size for each replicate and siever
avg_fly_size_gen9 <- gen9_data %>%
  group_by(replicate, siever) %>%
  summarise(avg_fly_size = mean(sieve_expanded, na.rm = TRUE), .groups = "drop")

# Print the new dataset to verify
print(avg_fly_size_gen9)
```

    ## # A tibble: 7 × 3
    ##   replicate siever avg_fly_size
    ##   <fct>     <fct>         <dbl>
    ## 1 O1        Kati           989.
    ## 2 O2        Misa          1087.
    ## 3 O3        Misa          1064.
    ## 4 O4        Misa          1097.
    ## 5 O5        Kati          1006.
    ## 6 O6        Kati          1055.
    ## 7 O7        Kati          1053.

``` r
##this for all other generations
avg_fly_size <- pseudo_continuous_female_O_850 %>%
  group_by(generation, replicate) %>%
  summarise(avg_fly_size = mean(sieve_expanded, na.rm = TRUE), .groups = "drop")  # Compute mean sieve size

print(avg_fly_size)
```

    ## # A tibble: 54 × 3
    ##    generation replicate avg_fly_size
    ##    <fct>      <fct>            <dbl>
    ##  1 1          O1               1077.
    ##  2 1          O2               1082.
    ##  3 1          O3               1134.
    ##  4 1          O4               1037.
    ##  5 1          O5               1046.
    ##  6 1          O6               1078.
    ##  7 1          O7               1016.
    ##  8 2          O1               1051.
    ##  9 2          O2               1072.
    ## 10 2          O3               1083.
    ## # ℹ 44 more rows

## Visualizatoin of the avg fly size for the optimized regime

``` r
avg_fly_optimized <- ggplot(avg_fly_size, aes(x = factor(generation), y = avg_fly_size)) +
  geom_boxplot(fill = "blue", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, color = "black") +
  labs(title = "Average Fly Size per Generation (Increasing selection Replicates 90->82%)", 
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
print(avg_fly_optimized)
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave(filename = "figs/Avg_rep_gen.png",
       plot = avg_fly_optimized,
       width = 10,    # less wide, adjust as needed
       height = 8,   # adjust height as desired
       units = "in")

avg_fly_optimzed_labeled <- ggplot(avg_fly_size, aes(x = factor(generation), y = avg_fly_size)) +
  geom_boxplot(fill = "blue", alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, aes(color = replicate)) +
  geom_text(aes(label = replicate), 
            position = position_jitter(width = 0.2, height = 0), 
            vjust = -1.5, size = 3.5) +
  labs(title = "Average Fly Size per Generation (Increasing Selection Replicates 90->82%)",
       x = "Generation", y = "Mean Fly Size (µm)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

print(avg_fly_optimzed_labeled)
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

``` r
ggsave(filename = "figs/Avg_rep_gen_labeled_rep.png",
       plot = avg_fly_optimzed_labeled,
       width = 10,    # less wide, adjust as needed
       height = 8,   # adjust height as desired
       units = "in")
```

## Pseudo continuous data for controls and constant then coputing average size

``` r
# ---- Pseudo-continuous data for Constant Replicates (K) ----
pseudo_continuous_female_K <- datap %>%
  filter(replicate %in% c("K1", "K2", "K3", "K6", "K7", "K8"), generation %in% c(1, 2, 8)) %>%
  filter(!is.na(female_num) & female_num > 0) %>%  
  uncount(female_num, .remove = FALSE) %>%
  mutate(sieve_expanded = ifelse(sieve_size == 1102, 1085, sieve_size))  # Merge 1102 into 1085

# Compute average fly size per replicate and generation
avg_fly_size_K <- pseudo_continuous_female_K %>%
  group_by(generation, replicate) %>%
  summarise(avg_fly_size = mean(sieve_expanded, na.rm = TRUE), .groups = "drop")

print(avg_fly_size_K)
```

    ## # A tibble: 15 × 3
    ##    generation replicate avg_fly_size
    ##    <fct>      <fct>            <dbl>
    ##  1 1          K1               1102.
    ##  2 1          K2               1059.
    ##  3 1          K3               1066.
    ##  4 1          K6               1052.
    ##  5 1          K7               1027.
    ##  6 1          K8               1049.
    ##  7 2          K1               1049.
    ##  8 2          K2               1061.
    ##  9 2          K3               1056.
    ## 10 2          K6               1041.
    ## 11 2          K7               1045.
    ## 12 2          K8               1030.
    ## 13 8          K1               1049.
    ## 14 8          K6               1093.
    ## 15 8          K8               1167

``` r
# ---- Pseudo-continuous data for Control Replicates (C) ----

pseudo_continuous_female_C <- datap %>%
  filter(replicate %in% c("C6", "C7", "C8", "C9", "C10")) %>%
  filter((replicate %in% c("C6", "C7", "C8", "C9", "C10") & generation == 2) |  
           (replicate == "C7" & generation == 8)) %>%
  filter(!is.na(female_num) & female_num > 0) %>%  
  uncount(female_num, .remove = FALSE) %>%
  mutate(sieve_expanded = ifelse(sieve_size == 1102, 1085, sieve_size))
pseudo_continuous_female_C %>%
  filter(replicate == "C7", generation == 2)  %>% 
   head(10)
```

    ##    generation  regime       date replicate sieve_num sieve_size female_num
    ## 1           2 control 2024-11-22        C7         7       1250          5
    ## 2           2 control 2024-11-22        C7         7       1250          5
    ## 3           2 control 2024-11-22        C7         7       1250          5
    ## 4           2 control 2024-11-22        C7         7       1250          5
    ## 5           2 control 2024-11-22        C7         7       1250          5
    ## 6           2 control 2024-11-22        C7         8       1180          3
    ## 7           2 control 2024-11-22        C7         8       1180          3
    ## 8           2 control 2024-11-22        C7         8       1180          3
    ## 9           2 control 2024-11-22        C7         9       1120          7
    ## 10          2 control 2024-11-22        C7         9       1120          7
    ##    male_num sexing_counting siever selected_females sieve_expanded
    ## 1         0       misa_anna  Siraj               NA           1250
    ## 2         0       misa_anna  Siraj               NA           1250
    ## 3         0       misa_anna  Siraj               NA           1250
    ## 4         0       misa_anna  Siraj               NA           1250
    ## 5         0       misa_anna  Siraj               NA           1250
    ## 6         3       misa_anna  Siraj               NA           1180
    ## 7         3       misa_anna  Siraj               NA           1180
    ## 8         3       misa_anna  Siraj               NA           1180
    ## 9         2       misa_anna  Siraj               NA           1120
    ## 10        2       misa_anna  Siraj               NA           1120

``` r
# Compute average fly size per replicate and generation
avg_fly_size_C <- pseudo_continuous_female_C %>%
  group_by(generation, replicate) %>%
  summarise(avg_fly_size = mean(sieve_expanded, na.rm = TRUE), .groups = "drop")
```

\###Visualisation of constant and control mean size

``` r
# ---- Visualization: Boxplot of Mean Fly Size Across Generations ----
mean_size_K <- ggplot(avg_fly_size_K, aes(x = factor(generation), y = avg_fly_size)) +
  geom_boxplot(data = avg_fly_size_K %>% filter(generation != 8),  
               fill = "red", alpha = 0.7) +  # Exclude Gen 8 from boxplot
  geom_jitter(width = 0, aes(color = replicate), size = 2) +  # Missing "+"
  geom_text(aes(label = replicate), position = position_jitter(width = 0.2, height = 0),  
            vjust = -1.5, size = 3.5) +  # Missing "+"
  labs(title = "Average Fly Size per Generation (Constant Selection Replicates)",
       x = "Generation", y = "Mean Fly Size (µm)") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),        # Bold, larger title
    axis.title.x = element_text(size = 14, face = "bold"),        # Bold, larger x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),        # Bold, larger y-axis title
    legend.title = element_text(size = 12, face = "bold"),        # Bold legend title
    legend.text = element_text(size = 10)                         # Larger legend text
  )
ggsave(filename = "figs/mean_size_K.png",
       plot = mean_size_K,
       width = 10,    # less wide, adjust as needed
       height = 8,   # adjust height as desired
       units = "in")

 
mean_control <- ggplot(avg_fly_size_C, aes(x = factor(generation), y = avg_fly_size, color = replicate)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = replicate), vjust = -1.2, size = 4) +
  labs(title = "Mean Fly Size Across Generations (Control Replicates)",
       x = "Generation", y = "Mean Fly Size (µm)") +
  scale_color_manual(values = c("red", "blue", "green", "purple")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),        # Bold, larger title
    axis.title.x = element_text(size = 14, face = "bold"),        # Bold, larger x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),        # Bold, larger y-axis title
    legend.title = element_text(size = 12, face = "bold"),        # Bold legend title
    legend.text = element_text(size = 10)                         # Larger legend text
  )


ggsave(filename = "figs/mean_control.png",
       plot = mean_control,
       width = 20,    # less wide, adjust as needed
       height = 10,   # adjust height as desired
       units = "in")
```

# Linear mixed model for the increasing regime replicates

``` r
# ---- Linear Mixed Model for Optimized Replicates ----

pseudo_continuous_female_O_copy <- pseudo_continuous_female_O

# First, ensure 'generation' is a factor and set the levels in the order you want.
pseudo_continuous_female_O_copy$generation <- factor(pseudo_continuous_female_O_copy$generation, 
                                                levels = c("1", "2", "3", "5", "6", "7", "8", "9"))        

# Create dummy variables for all non-reference levels (all except "1")
levels_gen <- levels(pseudo_continuous_female_O_copy$generation) 
for (lev in levels_gen[-1]) {   
  dummy_name <- paste0("gen_dummy", lev)   
  pseudo_continuous_female_O_copy[[dummy_name]] <- as.numeric(pseudo_continuous_female_O_copy$generation == lev)    
  
  # center the dummy variable
  pseudo_continuous_female_O_copy[[dummy_name]] <- pseudo_continuous_female_O_copy[[dummy_name]] -      
    mean(pseudo_continuous_female_O_copy[[dummy_name]]) 
}   

# Fit the linear mixed model (adding gen_dummy9)
lmm_result_O <- lme4::lmer(sieve_expanded ~ generation + siever +  
                             (1 + gen_dummy2 + gen_dummy3 + gen_dummy5 +  
                                gen_dummy6 + gen_dummy7 + gen_dummy8 + gen_dummy9 | replicate),  
                           data = pseudo_continuous_female_O_copy)  

summary(lmm_result_O)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: sieve_expanded ~ generation + siever + (1 + gen_dummy2 + gen_dummy3 +  
    ##     gen_dummy5 + gen_dummy6 + gen_dummy7 + gen_dummy8 + gen_dummy9 |  
    ##     replicate)
    ##    Data: pseudo_continuous_female_O_copy
    ## 
    ## REML criterion at convergence: 574392.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -5.6146 -0.7031  0.0956  0.6445  8.9198 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev. Corr                               
    ##  replicate (Intercept)  218.6   14.79                                       
    ##            gen_dummy2   806.9   28.41     0.54                              
    ##            gen_dummy3  4968.2   70.49     0.10  0.62                        
    ##            gen_dummy5   707.3   26.60     0.61  0.81  0.44                  
    ##            gen_dummy6  3228.9   56.82     0.55  0.42 -0.16  0.73            
    ##            gen_dummy7  3312.7   57.56    -0.05  0.37  0.87  0.07 -0.59      
    ##            gen_dummy8  1997.1   44.69     0.55  0.95  0.54  0.74  0.27  0.42
    ##            gen_dummy9  2917.9   54.02     0.04  0.79  0.87  0.51 -0.10  0.76
    ##  Residual              3533.5   59.44                                       
    ##       
    ##       
    ##       
    ##       
    ##       
    ##       
    ##       
    ##       
    ##   0.78
    ##       
    ## Number of obs: 52167, groups:  replicate, 7
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error t value
    ## (Intercept) 1042.461     11.044  94.388
    ## generation2   -9.025     10.789  -0.836
    ## generation3  -74.167     26.681  -2.780
    ## generation5  -39.111     10.105  -3.871
    ## generation6  -62.010     21.556  -2.877
    ## generation7  -20.256     21.785  -0.930
    ## generation8  -57.307     16.941  -3.383
    ## generation9  -24.024     20.471  -1.174
    ## sieverMisa    73.699      1.400  52.640
    ## sieverSiraj    4.689      1.582   2.965
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) gnrtn2 gnrtn3 gnrtn5 gnrtn6 gnrtn7 gnrtn8 gnrtn9 sivrMs
    ## generation2 -0.675                                                        
    ## generation3 -0.883  0.617                                                 
    ## generation5 -0.468  0.806  0.439                                          
    ## generation6  0.109  0.420 -0.154  0.722                                   
    ## generation7 -0.744  0.369  0.866  0.072 -0.585                            
    ## generation8 -0.625  0.941  0.537  0.736  0.268  0.422                     
    ## generation9 -0.965  0.791  0.866  0.508 -0.099  0.763  0.777              
    ## sieverMisa  -0.078  0.002 -0.015  0.010  0.003  0.023  0.019  0.013       
    ## sieverSiraj -0.118  0.001  0.025 -0.014  0.059  0.019  0.036  0.053  0.401
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
r_squared <- r.squaredGLMM(lmm_result_O)
print(r_squared)
```

    ##            R2m       R2c
    ## [1,] 0.1903245 0.4286759

``` r
isSingular(lmm_result_O)
```

    ## [1] TRUE

## Visualisation of the LMM effects and diagnostics

``` r
library(effects)
library(lattice)

# EFFECTS PLOT
# FIRST: Plot inline in Markdown
plot(allEffects(lmm_result_O))
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
# SECOND: Open file, plot again, and save 
png("figs/Mixed_model_effects.png", width = 12, height = 8, units = "in", res = 400)
plot(allEffects(lmm_result_O))
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("figs/Mixed_model_effects.svg", width = 12, height = 8)
plot(allEffects(lmm_result_O))
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# RANDOM EFFECTS PLOT
ranef_plot <- ranef(lmm_result_O, condVar = TRUE)
raneff <- dotplot(ranef_plot)

# FIRST: Inline plot (displayed in markdown):
print(raneff)
```

    ## $replicate

![](1_data_prep_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

``` r
# SECOND: explicitly save to files 
png("figs/Mixed_model_raneff.png", width = 12, height = 8, units = "in", res = 400)
print(raneff)
```

    ## $replicate

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("figs/Mixed_model_raneff.svg", width = 12, height = 8)
print(raneff)
```

    ## $replicate

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# MODEL DIAGNOSTICS PLOT
source("/Users/selmasri/Desktop/scripts/mundry_mixedmodelscourse/functions/diagnostic_fcns.r")

# FIRST: Show inline plot clearly in Markdown:
diagnostics.plot(lmm_result_O)
```

![](1_data_prep_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->

``` r
# Then SAVE plots explicitly as files:
png("figs/Mixed_model_diagnostics.png", width = 12, height = 8, units = "in", res = 400)
diagnostics.plot(lmm_result_O)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
svg("figs/Mixed_model_diagnostics.svg", width = 12, height = 8)
diagnostics.plot(lmm_result_O)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Extract residuals and fitted values 
residuals_O <- residuals(lmm_result_O)
fitted_values_O <- fitted(lmm_result_O)
```

Should i give kendal test a try?
