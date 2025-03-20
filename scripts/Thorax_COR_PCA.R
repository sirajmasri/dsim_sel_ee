# ========================== R script for analysis of images phenotype data ==========================
setwd("/Users/selmasri/Phenotype_data_analysis_feb2025")
getwd()
install.packages("psych")
install.packages("corrplot")
install.packages("factoextra") 
library(emmeans)
library(factoextra)
library(corrplot)
library(psych)
library(dplyr)
library(ggplot2)
library(car)
library(tidyr)
library(dplyr)
install.packages("GGally") 
library(GGally)
library(lmtest)

main_df <- read.csv("Thorax_measurments - Thorax_measurments_unbiased.csv", header =  TRUE, sep = ",")
str(main_df)
main_df$Date_measured <- as.Date(main_df$Date_measured, format = "%d/%m/%Y")
str(main_df)
# Aggregating replicate measurements by Fly_ID (and Date_measured if needed)
main_df_tidy <- main_df %>%
  group_by(Fly_ID, Date_measured) %>%
  summarise(
    Thorax_length_avg = mean(Thorax_length.mm., na.rm = TRUE),
    Full_body_length_avg = mean(Full_body_length.mm., na.rm = TRUE),
    Thorax_width_avg = mean(Thorax_width.mm., na.rm = TRUE),
    Wing_length_avg = mean(wing.length.mm., na.rm = TRUE)
  )

##now i have to load the rest of flies from the measurments
#the upcoming code is given by chatgpt since I dind't know anyway to bring those scattered,
  #individual flies files together in one file while each one has 11 measurments
## the code below given by chatgpt is a bit dangerous though so i had to check all files measurments because it 
#averages the first 3 rows as thorax length then 3 rows after that as body length then 3 rows as thor width than 2 as wing length
##if in any individual fly the rows are mixed up so i did my measruments not in order the average will be totally wrong
##for the futre I should avoid this type of organizing files, and especially with different unlabeled measurments in the csv files
# i could modify the macros to label inside the csv file or by awk or something
#now i was able to check everyone of them by eye because they are only 28


# ---- Loading individual csv files Data by chatgpt ----

# Path to your folder of individual CSVs
path <- "/Users/selmasri/Desktop/Project/measurments_results_gwas1"

# List all CSV files in that directory
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

# A small helper function to parse Fly_ID from a filename like "fly_0189_results.csv"
get_fly_id <- function(filename) {
  base <- basename(filename)               # e.g. "fly_0189_results.csv"
  no_ext <- sub("\\.csv$", "", base)       # "fly_0189_results"
  fly_id <- sub("_results$", "", no_ext)   # "fly_0189"
  return(fly_id)
}

# A helper function that safely averages the "Length" column for a given set of row indices.
# If the file doesn't have enough rows, it returns NA.
get_mean <- function(df, idx, col = "Length") {
  if (all(idx <= nrow(df))) {
    mean(df[idx, col], na.rm = TRUE)
  } else {
    NA
  }
}

list_of_dfs <- lapply(files, function(file) {
  # 1) Parse the Fly_ID from the filename
  fly_id <- get_fly_id(file)
  
  # 2) Read the CSV
  df <- read.csv(file)
  
  # 3) For each trait, average the correct rows if they exist
  #    (Adjust these if your data is 1-indexed differently)
  thorax_length_avg   <- get_mean(df, 1:3,  "Length")
  body_length_avg     <- get_mean(df, 4:6,  "Length")
  thorax_width_avg    <- get_mean(df, 7:9,  "Length")
  wing_length_avg     <- get_mean(df, 10:11,"Length")  # may be NA if only 9 rows
  
  # 4) Create a single-row data frame with all traits
  out <- data.frame(
    Fly_ID             = fly_id,
    Date_measured      = as.Date("2025-01-27"),  # Hard-code this date
    Thorax_length_avg  = thorax_length_avg,
    Full_body_length_avg = body_length_avg,
    Thorax_width_avg   = thorax_width_avg,
    Wing_length_avg    = wing_length_avg
  )
  
  return(out)
})

# Combine all mini data frames for flies 188–215 into one
other_df <- bind_rows(list_of_dfs)

# Check the result
str(other_df)
head(other_df)

# ---- Merging the two data by chatgpt ----

# 1) Remove flies 188 through 215 from the main data frame(note that original one had empty rows from 188 to 215)
#    Here, we're filtering out any Fly_ID that matches "fly_0188" through "fly_0215".
main_df_no_placeholders <- main_df_tidy %>%
  filter(
    !(Fly_ID %in% paste0("fly_0", 188:215))
  )

# 2) Stack the main data (with placeholders removed) and the other data
all_flies_df <- bind_rows(main_df_no_placeholders, other_df)

# 3) Inspect the final combined data
str(all_flies_df)
head(all_flies_df)
# change the outlier 
# ---- Correlation between measurments ----
##correlation between measurments
#creating a vector with column of interest
measure_cols <- c("Thorax_length_avg", "Full_body_length_avg", 
                  "Thorax_width_avg", "Wing_length_avg")
#subseting for columns of interest
sub_df <- all_flies_df[, measure_cols]

# Perform pairwise correlation with corr.test() storing in res_psych to use,
#the psych pacakge beacuse many rows have only two obersvations and some have 3 some have 4, and i dont want to omit any data
#bonferroni p value correction
res_psych <- corr.test(sub_df, use = "pairwise", method = "pearson", adjust = "bonferroni")

# View the results
res_psych$r   # Correlation matrix
res_psych$n   # Number of observations used for each pair
res_psych$p   # p-values for each correlation
##plotting



pairsplot <- ggpairs(sub_df,
                     # Lower panels: scatterplots with a linear smoothing line and transparent points
                     lower = list(continuous = wrap("smooth", method = "lm", se = FALSE, na.rm = TRUE, alpha = 0.4)),
                     # Diagonal panels: density plots or histograms
                     diag = list(continuous = "densityDiag"),
                     # Upper panels: display correlation coefficients using pairwise complete observations
                     upper = list(continuous = wrap("cor", use = "pairwise.complete.obs", method= "pearson")),
                     title = "Pairs Plot of Fly Traits"
) + theme_bw()

# Display the plot
pairsplot

# Save the plot to a file
ggsave(
  filename = "/Users/selmasri/Phenotype_data_analysis_feb2025/pairsplotphenotypes.png",
  plot = pairsplot,
  width = 10,
  height = 10,
  units = "in"
)



summary(sub_df)
range(sub_df$Thorax_length_avg, na.rm = TRUE)
range(sub_df$Thorax_width_avg, na.rm = TRUE)
range(sub_df$Wing_length_avg, na.rm = TRUE)

# ---- Which traits correlates the most with bin size ----
# load the mapping file 
mapping <- read.csv("/Users/selmasri/Desktop/Project/mapping.csv", stringsAsFactors = FALSE)

# Check the first rows
head(mapping)

### I will merge the data so i know later which fly has which id

# adding .tif extension to Fly_ID so the merging works
all_flies_df$Fly_ID <- ifelse(grepl("\\.tif$", all_flies_df$Fly_ID),
                              all_flies_df$Fly_ID,
                              paste0(all_flies_df$Fly_ID, ".tif"))

# merge the data 
merged_data <- merge(all_flies_df, mapping, 
                     by.x = "Fly_ID", by.y = "Anonymized_ID", 
                     all.x = TRUE)

# Check the merged data again
head(merged_data)
# Remove the ".tif" extension from Original_File
merged_data$original_id_clean <- gsub("\\.tif$", "", merged_data$Original_File)

# Split the cleaned original_id into parts using the underscore "_" as the delimiter
meta_list <- strsplit(merged_data$original_id_clean, "_")

# Convert the list to a data frame; we expect 5 parts: Batch, SieveSize, Plate, Column, and Row.
meta_df <- do.call(rbind, meta_list)
colnames(meta_df) <- c("batch", "sieve_size", "plate", "column", "row")

# Convert meta_df to a proper data frame
meta_df <- as.data.frame(meta_df, stringsAsFactors = FALSE)

# Combine the metadata with merged_data
merged_data <- cbind(merged_data, meta_df)

# Check the merged data to see the new metadata columns
head(merged_data)


long_data <- merged_data %>%
  pivot_longer(cols = c("Thorax_length_avg", "Full_body_length_avg", 
                        "Thorax_width_avg", "Wing_length_avg"),
               names_to = "Trait",
               values_to = "Value")



long_data$sieve_size <- factor(long_data$sieve_size, levels = sort(unique(long_data$sieve_size)))

traits_dis_by_sieve_size <- ggplot(long_data, aes(x = factor(sieve_size, levels = sort(unique(as.numeric(as.character(sieve_size))))), 
                                                  y = Value, fill = factor(sieve_size))) +
  geom_boxplot() +
  facet_wrap(~ Trait, scales = "free_y") +
  labs(x = "Sieve Size", y = "Measurement",
       title = "Fly Traits by Sieve Size") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold")  # Make facet labels bold & bigger
  )

ggsave(
  filename = "/Users/selmasri/Phenotype_data_analysis_feb2025/traits_by_sieve_size.png",
  plot = traits_dis_by_sieve_size,
  width = 15,
  height = 10,
  units = "in"
)


ggplot(long_data, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  facet_wrap(~ Trait, scales = "free") +
  labs(
    title = "Distribution of Fly Traits",
    x = "Measurement Value (mm)",
    y = "Count"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.y  = element_text(face = "bold", size = 12),
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5)
  )




###PCA
##omit NAs for the pca analysis to have complete data for the 4 traits
pca_df <- na.omit(sub_df)
pca_res <- prcomp(pca_df, center = TRUE, scale. = TRUE)
summary(pca_res)
pca_res$rotation #gives the loadings (i.e., how each original variable contributes to each PC).
pca_res$x #gives the coordinates of each observation (fly) in the new PCA space.
##
# Retrieve the sieve size for the rows used in the PCA
#grouping <- merged_data[rownames(pca_df), "sieve_size"]
# i can do the grouping for date, plate or whatever and rerun the pca
##all my flies come from same batch
# Calculate percentage of variance for PC1 and PC2 then plot
eig.val <- pca_res$sdev^2
var.percent <- eig.val / sum(eig.val) * 100
pc1 <- round(var.percent[1], 1)
pc2 <- round(var.percent[2], 1)
#fviz_pca_biplot(pca_res,
#label = "var",      # label only variables (arrows)
                #repel = TRUE,       # reduce label overlap
                #geom = "point",     # plot individuals as points
                #habillage = grouping,   # color points by sieve size
                #addEllipses = TRUE,      # optional: add confidence ellipses
                #palette = "jco") +       # or use any palette you prefer
  # labs(x = paste0("PC1 (", pc1, "%)"),
#y = paste0("PC2 (", pc2, "%)"),
#title = "PCA Biplot Colored by Sieve Size") +
  # theme_classic() +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  # geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

 fviz_pca_biplot(pca_res,
                label = "var",
                repel = TRUE,
                geom = "point") +
  labs(
    x = paste0("PC1 (", pc1, "%)"),
    y = paste0("PC2 (", pc2, "%)")
  ) +
  theme_classic() +                                # Use a classic theme for clearer axes
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Add horizontal axis line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey")      # Add vertical axis line


 
 
 

# ---- linear models for different traits as response sieves as predictors ----

# 1. Ensure sieve_size is a factor (adjust levels if needed)
# convert to a factor
merged_data$sieve_size <- factor(merged_data$sieve_size, 
                                 levels = sort(unique(merged_data$sieve_size)))
 levels(merged_data$sieve_size)
 
 # Replace the level "800" with "850" so the value doesnt drive significance
 levels(merged_data$sieve_size)[levels(merged_data$sieve_size) == "800"] <- "850"
 
 
 levels(merged_data$sieve_size)
 # (Optional) Reorder levels if needed
 merged_data$sieve_size <- factor(merged_data$sieve_size, levels = sort(unique(merged_data$sieve_size)))
 
########fitting separate linear models to see which trait is explained best by sieve size

# Fit models for each trait
lm_thorax <- lm(Thorax_length_avg ~ sieve_size, data = merged_data)
lm_body   <- lm(Full_body_length_avg ~ sieve_size, data = merged_data)
lm_width  <- lm(Thorax_width_avg ~ sieve_size, data = merged_data)
lm_wing <- lm(Wing_length_avg ~ sieve_size, data = merged_data)
lm_wing_noextremes   <- lm(Wing_length_avg ~ sieve_size, data = merged_data_noextremes)
lm_thorax_date <- lm(Thorax_length_avg ~ sieve_size + as.numeric(Date_measured), data = merged_data)
summary(lm_thorax_date)
anova(lm_thorax, lm_thorax_date)
summary(lm_thorax)
unique_dates <- unique(merged_data$Date_measured)
length(unique_dates)
unique_dates

merged_data$DateFactor <- factor(merged_data$Date_measured)
lm_thorax_date_factor <- lm(Thorax_length_avg ~ sieve_size + DateFactor, data = merged_data)
summary(lm_thorax_date_factor)
table_counts <- table(merged_data$Date_measured, merged_data$sieve_size)

print(table_counts)

chisq_result <- chisq.test(table_counts)

print(chisq_result).  ####conclusion the high inflation in R2 result with date measured included is a bias caused by dependance(OR CORR) between sive ize and date measured in this data

# Extract adjusted R-squared values for each model.
adj_r2_thorax <- summary(lm_thorax)$adj.r.squared
adj_r2_body   <- summary(lm_body)$adj.r.squared
adj_r2_width  <- summary(lm_width)$adj.r.squared
adj_r2_wing   <- summary(lm_wing)$adj.r.squared

# Print the values
print(paste("Adjusted R² for Thorax:", adj_r2_thorax))
print(paste("Adjusted R² for Body:", adj_r2_body))
print(paste("Adjusted R² for Width:", adj_r2_width))
print(paste("Adjusted R² for Wing:", adj_r2_wing))

adj_r2_values <- data.frame(
  Trait = c("Thorax", "Body", "Width", "Wing"),
  Adj_R2 = c(adj_r2_thorax, adj_r2_body, adj_r2_width, adj_r2_wing)
)

## the last sieve was showing a decrease for most flies due to techincal issue because we notice it has little flies and mostly males so i rerun the association by removing extremes
# Remove rows where sieve_size is 800, 1319, or 1400
merged_data_noextremes <- merged_data[!(merged_data$sieve_size %in% c(800, 1319, 1400)), ]

# Ensure sieve_size is still treated as a factor with updated levels
merged_data_noextremes$sieve_size <- factor(merged_data_noextremes$sieve_size)

# Check the unique sieve sizes remaining
unique(merged_data_noextremes$sieve_size)

# Fit linear models with filtered data
lm_thorax_filtered <- lm(Thorax_length_avg ~ sieve_size, data = merged_data_noextremes)
lm_body_filtered  <- lm(Full_body_length_avg ~ sieve_size, data = merged_data_noextremes)
lm_width_filtered  <- lm(Thorax_width_avg ~ sieve_size, data = merged_data_noextremes)
lm_wing_filtered   <- lm(Wing_length_avg ~ sieve_size, data = merged_data_noextremes)

# Model summaries
summary(lm_thorax_filtered)
summary(lm_body_filtered)
summary(lm_width_filtered)
summary(lm_wing_filtered)
#####COCNLUSION FILTERING THE DATA DOES NOT CHANGE THE MOST ASSOCIATED TRAIT SO IT IS NOT USEFUL THE CODE ABOVE OF FILTERING ETC

ggplot(adj_r2_values, aes(x = Trait, y = Adj_R2)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  labs(title = "Adjusted R² for Each Trait", x = "Trait", y = "Adjusted R²") +
 theme_bw()
######Model duagnostics and emmeans 
# 3. Model diagnostics
# Plot diagnostic plots to check for normality, homoscedasticity, and influential points.
plot(lm_thorax, which = 1)  # Residuals vs. Fitted
plot(lm_thorax, which = 2)
plot(lm_wing, which = 1)
plot(lm_wing, which = 2)
# Normal Q-Q plot
merged_data[28, ]
##fucntion to test model assumptions
test_model_assumptions <- function(model) {
  # 1. Normality of residuals: Shapiro–Wilk
  shapiro_res <- shapiro.test(residuals(model))
  
  # 2. Homoscedasticity: Breusch–Pagan
  bptest_res <- bptest(model)
  
  cat("\nShapiro–Wilk Test for Normality:\n")
  print(shapiro_res)
  
  cat("\nBreusch–Pagan Test for Homoscedasticity:\n")
  print(bptest_res)
}

# Usage example (for your four models):
test_model_assumptions(lm_thorax)
test_model_assumptions(lm_body)
test_model_assumptions(lm_width)
test_model_assumptions(lm_wing)
test_model_assumptions(lm_wing_noextremes)
# ---- Summary of Model Assumptions ----
# 1. lm_thorax
#    - SW p-value = 0.01537 -> Residuals are NOT normal (p < 0.05)
#    - BP p-value = 0.2188  -> Residuals are homoscedastic (p > 0.05)
#
# 2. lm_body
#    - SW p-value = 0.5119  -> Residuals are normal (p > 0.05)
#    - BP p-value = 0.01303 -> Residuals are NOT homoscedastic (p < 0.05)
#
# 3. lm_width
#    - SW p-value = 0.6716  -> Residuals are normal (p > 0.05)
#    - BP p-value = 0.0365  -> Residuals are NOT homoscedastic (p < 0.05)
#
# 4. lm_wing
#    - SW p-value = 0.003986 -> Residuals are NOT normal (p < 0.05)
#    - BP p-value = 0.3059   -> Residuals are homoscedastic (p > 0.05)


## Here i will try to log transform the data
merged_data$Wing_length_log <- log(merged_data$Wing_length_avg)
merged_data$Thorax_length_log <- log(merged_data$Thorax_length_avg)
lm_wing_log <- lm(Wing_length_log ~ sieve_size, data = merged_data)
lm_thorax_log <- lm(Thorax_length_log ~ sieve_size, data = merged_data)

test_model_assumptions(lm_wing_log) ##did not improve anything
test_model_assumptions(lm_thorax_log)## the residuals are normally distriubted now

adj_r2_thorax_log <- summary(lm_thorax_log)$adj.r.squared
adj_r2_thorax_log ## now the adjusted r squared is a little bit higher than before but stilll lower than wing etc


emm <- emmeans(lm_model, specs = "sieve_size")
print(emm)
# Pairwise comparisons among levels of sieve_size:
pairwise_comp <- pairs(emm)
print(pairwise_comp)



emm_thorax <- emmeans(lm_thorax, ~ sieve_size)
plot(emm_thorax)

# Function to compare null vs full model for a given trait, I still dont get exactly whats hapenning here, so dismiss it
run_model_comparison <- function(trait, data) {
  # Build formulas for the null model and full model
  full_formula <- as.formula(paste(trait, "~ sieve_size"))
  null_formula <- as.formula(paste(trait, "~ 1"))
  
  # Fit the models
  model_null <- lm(null_formula, data = data)
  model_full <- lm(full_formula, data = data)
  
  # Compare the models using ANOVA
  comparison <- anova(model_null, model_full)
  return(comparison)
}







# 1) Extract the unique sieve sizes as character
# 2) Convert them to numeric
# 3) Sort them
sorted_levels <- sort(as.numeric(unique(as.character(merged_data$sieve_size))))

# 4) Re-factor sieve_size with these sorted levels
merged_data$sieve_size <- factor(
  merged_data$sieve_size, 
  levels = sorted_levels
)

# Check that the levels are in ascending order
levels(merged_data$sieve_size)


# Create a subset for wing length (flies with non-missing wing length)
wing_data <- subset(merged_data, !is.na(Wing_length_avg))

# Drop unused factor levels in the subset (so "800" won't appear if it has no observations)
wing_data$sieve_size <- droplevels(wing_data$sieve_size)

# Check the levels for wing_data:
levels(wing_data$sieve_size)

ggplot(wing_data, aes(x = sieve_size, y = Wing_length_avg, fill = factor(sieve_size)) +
  geom_boxplot() +
  labs(title = "Wing Length by Sieve Size",
       x = "Sieve Size",
       y = "Wing Length (mm)") +
  theme_minimal() +
  theme(legend.position = "none")



ggplot(wing_data, aes(x = factor(sieve_size), y = Wing_length_avg, fill = factor(sieve_size))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "Wing Length by Sieve Size", 
       x = "Sieve Size", 
       y = "Wing Length (mm)") +
  theme_bw() +
  theme(legend.position = "none")







sum(merged_data$sieve_size == "1250")
sum(merged_data$sieve_size == "1250" & !is.na(merged_data$Wing_length_avg))




##later I can remove any fly id of 800 and substitue with 850 and do the analysis again of the linear models



# ---- linear models for different traits with each other ----







