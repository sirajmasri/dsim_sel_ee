# ========================== R script for analysis of phenotypic response # ==========================
setwd("/Users/selmasri/Phenotype_data_analysis_feb2025")
getwd()
# ---- Load Libraries ----
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(lubridate)

# ---- Load Data ----
datap <- read.csv("E&R_Phenotypic_Data_D.sim_phenotype_data.csv", header = T, sep = ";")
str(datap)
# ---- mutating variables as they should be ----
datap <- datap %>%
  mutate( female_num = as.numeric(female_num),  
    male_num = as.numeric(male_num), 
    selected_females = as.numeric(selected_females) 
  ) %>%
  select(-dead)  
str(datap)
sum(is.na(datap$date))
datap <- datap %>%
  mutate(date = lubridate::dmy(date))    ###converting the date column from chr to date
str(datap)
sum(is.na(datap$date))  #checking if the conversion worked

# ---- Phenotypic response per replicate ----
dataO1 <- datap %>% filter(replicate == "O1", generation %in% 1:6) #filtering the data for one replicate 
dataO2 <- datap %>% filter(replicate == "O2", generation %in% 1:6)
dataO3 <- datap %>% filter(replicate == "O3", generation %in% 1:6)

## calculating proportions since pop size is variable without collapsing 9 and 9.5
prop_dataO1 <- dataO1 %>% group_by(generation, sieve_num) %>%
  summarise(total_females = sum(female_num, na.rm = TRUE),
            total_males = sum(male_num, na.rm = TRUE),
            .groups = 'drop')%>% 
  group_by(generation)%>%
  mutate(female_prop = total_females/sum(total_females),
         male_prop = total_males / sum(total_males))
#### visualisation
informative_binsO1 <- prop_dataO1 %>%
  filter(sieve_num >= 4 & sieve_num <= 13) ##since i dont care about all bins or the ones that i am not using

ggplot(informative_binsO1, aes(x = sieve_num, y = female_prop, color = factor(generation))) +
  geom_line(size = 1.2) +  # Connects points with lines
  geom_point(size = 2) +   # Marks each value with a dot
  labs(title = "Proportion of Females Across Sieve Bins",
       x = "Sieve Bin",
       y = "Proportion of Females",
       color = "Generation") +
  scale_x_continuous(breaks = unique(informative_binsO1$sieve_num), trans = "reverse") +   # Reverse sieve bin order
  theme_minimal()

ggplot(informative_binsO1, aes(x = generation, y = female_prop, color = factor(sieve_num), group = sieve_num)) +
  geom_line(size = 1.2) +  # Connect points over generations
  geom_point(size = 2) +   # Add points for clarity
  labs(title = "Proportion of Females Across Generations (Sieve Bins 4-14)",
       x = "Generation",
       y = "Proportion of Females",
       color = "Sieve Bin") +
  scale_x_continuous(breaks = unique(informative_binsO1$generation)) +  # Ensure all generations are shown
  scale_color_viridis_d(option = "plasma") +  # Better color scale
  theme_minimal()


###peak bins and population size
peak_bins <- prop_dataO1 %>%
  group_by(generation) %>%
  slice_max(female_prop, n = 1, with_ties = FALSE) %>%  ##slice_max or slice_min select rows with largest or smallest values of a variable
  select(generation, sieve_num, female_prop) #select() is used to get specific columns from a data frame.
pop_size <- dataO1 %>%
  group_by(generation) %>%
  summarise(total_population = sum(female_num + male_num, na.rm = TRUE))
pop_sizeO2 <- dataO2 %>%
  group_by(generation) %>%
  summarise(total_population = sum(female_num + male_num, na.rm = TRUE))
pop_sizeO3 <- dataO2 %>%
  group_by(generation) %>%
  summarise(total_population = sum(female_num + male_num, na.rm = TRUE))
##i can transform the data into pseudo-continuous
pseudo_femaleO1 <- dataO1 %>%
  filter(!is.na(female_num) & female_num > 0) %>%  # Remove NA or zero counts
  group_by(generation) %>%
  summarise(sieve_expanded = list(rep(sieve_num, times = female_num))) %>%
  unnest(cols = sieve_expanded)

summary(pseudo_femaleO1)
str(pseudo_femaleO1)
# Histogram for each generation
ggplot(pseudo_femaleO1, aes(x = sieve_expanded, fill = factor(generation))) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 10) +
  facet_wrap(~generation, scales = "free_y") +
  labs(title = "Distribution of Sieve Values Across Generations",
       x = "Sieve Bin",
       y = "Count",
       fill = "Generation") +
  theme_minimal()

#check normality with Shapiro-Wilk test
by(pseudo_femaleO1$sieve_expanded, pseudo_femaleO1$generation, shapiro.test)
##all p values are highly significant the data is not normally distributed
##i can use wilcoxon rank sum to compare between generations pairwise and kruska wallis if comparing all generations
# convert the dataframe into separate vectors per generation
pseudo_female_vectors <- pseudo_femaleO1 %>%
  group_by(generation) %>%
  summarise(sieve_values = list(sieve_expanded)) %>%
  pull(sieve_values)  # Extract list of numeric vectors
#pairwise wilcox test between generations
str(pseudo_femaleO1)
###i should change the sieve_num to sieve size for the pseudo continous and make a loop for wilcoxon that increases by 1 so i can use it for later generations
##ican also plug it into an annova(lm for pseudo size for example lm pseudo ~generation and it should be a facotr
###i can have per replicate a column of pseudo continuous data, and then i will have everything in one data frame, so ican keep track of who sieved, create a column pop size and check 
lm_result<-lm(pseudolength~generation.f, pseudo_continuous_female)
summary(lm_result)
emmeans::emmeans(lm_result,pairwise~generation.f)


lm_result<-lm(pseudolength~generation.f*replicate + siever, pseudo_continuous_female)

wilcox.test(pseudo_femaleO1$sieve_expanded [generation =1], pseudo_femaleO1$sieve_expanded [generation =2])
wilcox.test(pseudo_female_vectors[[1]], pseudo_female_vectors[[2]])  #  1 vs 2
wilcox.test(pseudo_female_vectors[[2]], pseudo_female_vectors[[3]])  #  2 vs  3
wilcox.test(pseudo_female_vectors[[3]], pseudo_female_vectors[[4]])  # 3 vs 5
wilcox.test(pseudo_female_vectors[[4]], pseudo_female_vectors[[5]])  # 5 vs 6
###all pairwise comparison indicate a significant shift in distribution, now i can check the direction
pseudo_femaleO1 %>%
  group_by(generation) %>%
  summarise(median_sieve = median(sieve_expanded))
##plot
ggplot(pseudo_femaleO1, aes(x = factor(generation), y = sieve_expanded, fill = factor(generation))) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparison of Sieve Distributions Across Generations",
       x = "Generation",
       y = "Sieve Bin (Pseudo-Continuous)",
       fill = "Generation") +
  scale_fill_manual(values = c("green", "red", "black", "yellow", "pink")) +  # Custom colors
  scale_y_reverse(breaks = seq(4, 15, by = 1)) +  # Reverse Y-axis direction, excluding bin 3
  theme_minimal()







#####i can collapse 10 and 9.5
pseudo_femaleO1_10 <- pseudo_femaleO1 %>%
  mutate(sieve_expanded = ifelse(sieve_expanded == 9.5, 10, sieve_expanded))
ggplot(pseudo_femaleO1_10, aes(x = sieve_expanded, fill = factor(generation))) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 10) +
  facet_wrap(~generation, scales = "free_y") +
  labs(title = "Distribution of Sieve Values Across Generations",
       x = "Sieve Bin",
       y = "Count",
       fill = "Generation") +
  theme_minimal()
ggplot(pseudo_femaleO1_10, aes(x = factor(generation), y = sieve_expanded, fill = factor(generation))) +
  geom_boxplot(alpha = 0.7) +
  labs(title = "Comparison of Sieve Distributions Across Generations",
       x = "Generation",
       y = "Sieve Bin (Pseudo-Continuous)",
       fill = "Generation") +
  scale_fill_manual(values = c("green", "red", "black", "yellow", "pink")) +  # Custom colors
  scale_y_reverse(breaks = seq(4, 15, by = 1)) +  # Reverse Y-axis direction, excluding bin 3
  theme_minimal()

###create a pseudo continous data for all optimized replicats
pseudo_continuous_female <- datap %>%  
  filter(replicate %in% paste0("O", 1:7)) %>%  # Select only replicates O1 to O7
  mutate(
    sieve_expanded = sieve_num * female_num  # Multiply bin size by female count
  )
unique(pseudo_continuous_female$sieve_num)
unique(informative_bins$sieve_num)
informative_bins <- pseudo_continuous_female %>%
  filter(sieve_num >= 4 & sieve_num <= 14, replicate %in% paste0("O", 1:7))
unique(pseudo_continuous_female$sieve_num)
unique(informative_bins$sieve_num)
informative_bins %>%
  filter(sieve_num %in% c(10, 11)) %>%
  select(replicate, generation, sieve_num)

##I can compact 9.5 and 10 
informative_bins <- informative_bins %>%
  mutate(sieve_num = ifelse(sieve_num == 9.5, 10, sieve_num))
informative_bins <- informative_bins %>%
  group_by(generation, replicate, sieve_num) %>%
  summarise(female_prop = sum(female_prop, na.rm = TRUE), .groups = "drop")


###plotting all replicates

custom_colors <- c("red", "pink", "green", "black", "yellow", "purple", "orange", "blue", "brown", "cyan", "darkgreen")

ggplot(informative_bins, aes(x = generation, y = female_prop, color = factor(sieve_num), group = sieve_num)) +
  geom_line(size = 1.2) +  
  geom_point(size = 2) +  
  labs(title = "Proportion of Females Across Generations (Sieve Bins 4-14)",
       x = "Generation",
       y = "Proportion of Females",
       color = "Sieve Bin") +
  scale_x_continuous(breaks = unique(informative_bins$generation)) +  
  scale_fill_manual(values = custom_colors) +  # Uses a better, clearer color scale
  facet_wrap(~ replicate, ncol = 2) +  # Creates separate plots for each replicate
  theme_minimal(base_size = 14) +  # Increases font size for better readability
  theme(legend.position = "right")  # Moves legend for better spacing



###plotting all repolicates

ggplot(informative_bins, aes(x = factor(generation), y = sieve_num, fill = factor(generation))) +
  geom_boxplot(aes(weight = female_prop), alpha = 0.7, outlier.shape = NA) +  
  labs(title = "Distribution of Female Proportions Across Generations (Faceted by Replicate)",
       x = "Generation",
       y = "Sieve Bin",
       fill = "Generation") +
  scale_fill_manual(values = c("red", "pink", "green", "black", "yellow", "purple", "orange")) +  
  facet_wrap(~ replicate, ncol = 2) +  
  scale_y_reverse(breaks = c(10, 11, 13)) +  # Show only sieve bins 10, 11, and 13
  theme_minimal(base_size = 14) +  
  theme(legend.position = "right",
        strip.text = element_text(size = 12, face = "bold"))



###hypothesis testing

####transform into pseudo continuos and do analysis
##visualize
##hypothesis testing, siever etc

##all replicates comaprison
