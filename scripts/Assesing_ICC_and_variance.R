# ========================== R script for analysis of image fly measurments ==========================
setwd("/Users/selmasri/Phenotype_data_analysis_feb2025")
getwd()
library(dplyr)
library(ggplot2)
install.packages("irr")
datam <- read.csv("Thorax_measurments - Variance_estimation.csv", header = TRUE, sep = ",")
str(datam)
as.Date(datam$Date_measured, format ="%d/%m/%Y")
datam$Date_measured <- as.Date(datam$Date_measured, format="%d/%m/%Y")
str(datam)
# ---- Variance calculation ----
variance_summary <- dplyr::group_by(datam, Fly_ID) %>%
  dplyr::summarise(
    var_thorax   = stats::var(Thorax_length.mm.),
    sd_thorax    = stats::sd(Thorax_length.mm.),
    cv_thorax    = sd_thorax / base::mean(Thorax_length.mm.),
    
    var_full_body = stats::var(Full_body_length.mm.),
    sd_full_body  = stats::sd(Full_body_length.mm.),
    cv_full_body  = sd_full_body / base::mean(Full_body_length.mm.),
    
    var_wing     = stats::var(Wing.length),
    sd_wing      = stats::sd(Wing.length),
    cv_wing      = sd_wing / base::mean(Wing.length)
  )

print(variance_summary)
##VISUAL INSPECTION
plot_measurements <- ggplot2::ggplot(datam, aes(x = Fly_ID, y = Thorax_length.mm.)) +
  geom_boxplot() +
  scale_x_discrete(labels = function(x) gsub("fly_", "", x)) +
  labs(title = "Thorax Length Variation per Fly", 
       x = "Fly ID", 
       y = "Thorax Length (mm)") +
  theme_bw()

print(plot_measurements)

plot_measurements_body <- ggplot2::ggplot(datam, aes(x = Fly_ID, y = Full_body_length.mm.)) +
  geom_boxplot() +
  scale_x_discrete(labels = function(x) gsub("fly_", "", x)) +
  labs(title = "Body Length Variation per Fly", 
       x = "Fly ID", 
       y = "Body Length (mm)") +
  theme_bw()

print(plot_measurements_body)

plot_measurements_wing <- ggplot2::ggplot(datam, aes(x = Fly_ID, y = Wing.length)) +
  geom_boxplot() +
  scale_x_discrete(labels = function(x) gsub("fly_", "", x)) +
  labs(title = "Wing Length Variation per Fly", 
       x = "Fly ID", 
       y = "Wing Length (mm)") +
  theme_bw()

print(plot_measurements_wing)


ggsave(filename = "Variance_Thorax_Same_Day.png", plot = plot_measurements, width = 10, height = 6, dpi = 300)
ggsave(filename = "Variance_Body_Same_Day.png", plot = plot_measurements_body, width = 10, height = 6, dpi = 300)


# ---- Intra class calculation ----
datam <- dplyr::group_by(datam, Fly_ID) %>%            ###ADDING A MEASURMENT IDENTIFIER 1,2,3
  dplyr::mutate(measurement = dplyr::row_number()) %>%
  dplyr::ungroup()
#####making the data in a wide format for thorax length and assiging i to thorax_wide in this way each fly is a row with column are different measurment for that trait
thorax_wide <- tidyr::pivot_wider(
  datam,
  id_cols = Fly_ID,
  names_from = measurement,
  values_from = Thorax_length.mm.,
  names_prefix = "thorax_"
)

###calculating ICC using irr for thorax, exclude the first column since its not a measurment
icc_resultthorax <- irr::icc(as.data.frame(thorax_wide[ , -1]), model = "oneway") 
icc_resultthorax
icc_func<-function(phenotype.df,onetwoway,agreementconsistency)
{
  print(irr::icc(as.data.frame(phenotype.df[ , -1]),
                 model = onetwoway, 
                  type = agreementconsistency)) 
}


icc_func(thorax_wide,"oneway","consistency")
icc_func(thorax_wide,"twoway","agreement")
  
  (as.data.frame(thorax_wide[ , -1]), model = "twoway", 
         type = c("agreement"))
print(icc_resultthorax)

####for full body
full_body_wide <- tidyr::pivot_wider(
  datam,
  id_cols = Fly_ID,
  names_from = measurement,
  values_from = Full_body_length.mm.,
  names_prefix = "Full_body_length_"
)

###calculating ICC using irr for full body, exclude the first column since its not a measurment
icc_resultbody <- irr::icc(as.data.frame(full_body_wide[ , -1]), model = "oneway", type = "agreement")
print(icc_resultbody)
####for wing length
wing_length_wide <- tidyr::pivot_wider(
  datam,
  id_cols = Fly_ID,
  names_from = measurement,
  values_from = Wing.length,
  names_prefix = "Wing_length_"
)

###calculating ICC using irr for full body, exclude the first column since its not a measurment
icc_resultwing <- irr::icc(as.data.frame(wing_length_wide[ , -1]))
print(icc_resultwing)
##Visualization
# ICC results for the 3-measurement data:
icc_data1 <- data.frame(
  Trait = c("Thorax Length", "Full Body Length", "Wing Length"),
  ICC   = c(0.774, 0.955, 0.971),     
  Lower = c(0.594, 0.908, 0.94),      # Replace 0.80 with the lower CI for Wing Length
  Upper = c(0.894, 0.98, 0.987)        # Replace 0.90 with the upper CI for Wing Length
)


# Create the forest plot for ICC values
plot_icc1 <- ggplot(icc_data1, aes(x = Trait, y = ICC)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  ylim(0, 1) +
  labs(title = "ICC for 3 Measurements (Same Day) for Each Trait",
       x = "Trait",
       y = "Intraclass Correlation Coefficient (ICC)") +
  theme_bw()

print(plot_icc1)

# Save the plot as a PNG file
ggsave(filename = "icc1_forest_plot.png", plot = plot_icc1, width = 10, height = 6, dpi = 300)


####all the previous analysis is based on flies measured same day(within fly within day effect, but different days between flies)
# now i can include a fourth measurment which was taken on a different day within fly to assess the overall variability,
##although this might confound day effect with fly effect, assessing the repeatiality of my measurments under a condition of measurment of a fly on the same day is much more logical since 
##it is very unlikely that i will measure the fly different measurments on different days
#but in this way I'll be able to quantify how much additional variability the different day introduces
data2 <- read.csv("Thorax_measurments_Variance_estimation2days.csv", header = TRUE, sep = ",")
# ---- Variance calculation 2----
variance_summary2 <- dplyr::group_by(data2, Fly_ID) %>%
  dplyr::summarise(
    var_thorax   = stats::var(Thorax_length.mm., na.rm = TRUE),
    sd_thorax    = stats::sd(Thorax_length.mm., na.rm = TRUE),
    cv_thorax    = stats::sd(Thorax_length.mm., na.rm = TRUE) / base::mean(Thorax_length.mm., na.rm = TRUE),
    
    var_full_body = stats::var(Full_body_length.mm., na.rm = TRUE),
    sd_full_body  = stats::sd(Full_body_length.mm., na.rm = TRUE),
    cv_full_body  = stats::sd(Full_body_length.mm., na.rm = TRUE) / base::mean(Full_body_length.mm., na.rm = TRUE),
    
    var_wing     = stats::var(Wing.length, na.rm = TRUE),
    sd_wing      = stats::sd(Wing.length, na.rm = TRUE),
    cv_wing      = stats::sd(Wing.length, na.rm = TRUE) / base::mean(Wing.length, na.rm = TRUE)
  )

print(variance_summary2)

ggplot2::ggplot(data2, aes(x = Fly_ID, y = Thorax_length.mm.)) + ##I CAN CHANGE Y AXIS TO A DIFFERENT TRAIT
 geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
labs(title = "Thorax Length Variation per Fly", x = "Fly ID", y = "Thorax Length (mm)")

# ---- Intra class calculation 2 ----
data2 <- dplyr::group_by(data2, Fly_ID) %>%            ###ADDING A MEASURMENT IDENTIFIER 1,2,3
  dplyr::mutate(measurement = dplyr::row_number()) %>%
  dplyr::ungroup()
#####making the data in a wide format for thorax length and assiging i to thorax_wide in this way each fly is a row with column are different measurment for that trait
thorax_wide2 <- tidyr::pivot_wider(
  data2,
  id_cols = Fly_ID,
  names_from = measurement,
  values_from = Thorax_length.mm.,
  names_prefix = "thorax_"
)

###calculating ICC using irr for thorax, exclude the first column since its not a measurment
icc_resulthorax2 <- irr::icc(as.data.frame(thorax_wide2[ , -1]))
print(icc_resulthorax2)

####for full body
full_body_wide2 <- tidyr::pivot_wider(
  data2,
  id_cols = Fly_ID,
  names_from = measurement,
  values_from = Full_body_length.mm.,
  names_prefix = "Full_body_length_"
)

###calculating ICC using irr for full body, exclude the first column since its not a measurment
icc_resultbody2 <- irr::icc(as.data.frame(full_body_wide[ , -1]))
print(icc_resultbody2)
#In my thorax measurements, the ICC dropped from 0.774 to 0.742 when I included the fourth measurement from a different day, so the extra day adds some additional variability day‐to‐day differences,reducing repeatability.

#in my full body measurements i got the same ICC (0.955) regardless of whether I included the fourth measurement or not. So full body length more stable or repeatable across days
# ---- Visualization ----
# creating a data frame with my ICC estimates and confidence intervals
icc_data <- data.frame(
  Trait = rep(c("Thorax Length", "Full Body Length"), each = 2),
  MeasurementSet = rep(c("3 Measurements_Same_day", "Adding +1 measurment_dayeffect"), times = 2),
  ICC   = c(0.774, 0.742, 0.955, 0.955),
  Lower = c(0.594, 0.573, 0.908, 0.908),
  Upper = c(0.894, 0.873, 0.98, 0.98)
)

# Load ggplot2 and create the forest plot
plot_icc <- ggplot(icc_data, aes(x = MeasurementSet, y = ICC, color = Trait)) +
  geom_point(size = 3) +
geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2) +
 facet_wrap(~ Trait) +
  ylim(0, 1) +
 labs(title = "Comparison of ICC Values by Trait and Measurement Set",
                x = "Measurement Set",
                y = "Intraclass Correlation Coefficient (ICC)") +
  theme_bw()

# Print the plot
print(plot_icc)
ggsave(filename = "icc_forest_plot.png", plot = plot_icc, width = 10, height = 6, dpi = 300)



