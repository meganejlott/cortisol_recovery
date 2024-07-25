################################################################################
################################SET UP##########################################
################################################################################

#Load libraries
library(tidyverse)
library(ggpubr)

#Load in data
data = read_csv("./data/raw_data/solids_removal.csv")
metadata = read_csv("./data/raw_data/sample_list.csv")


#Make area numeric
data = data %>% mutate(area = as.numeric(area))

#Separate standards from samples
stds = data %>% filter(sample_description == "Std Bracket Sample") 
samples = data %>% filter(sample_description != "Std Bracket Sample")

################################################################################
###########################STANDARD CURVES######################################
################################################################################

#Take a quick look to know what we're expecting
stds %>%
  ggplot(aes(x= amount, y = area)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE, formula = y ~ x) +
  facet_wrap(~analyte) + 
  stat_cor(size = 4, label.x = 250, label.y = 4E6, r.digits = 4) +
  stat_regline_equation(size = 4, label.x = 250, label.y = 5E6) 

##The equation for cortisol should be around y = 50000 + 5700x
##The equation for cortisone should be around y = 43000 + 4400x

#Calculate standard curves by hand

##What are the LoDs for these? Based on the standard curves
#for cortisol, the LoD is ~ 0.49 ug/L
##for cortisone, the LoD is ~ 0.98 ug/L

##Pull out the names of the analytes
analytes = data.frame(unique(data$analyte))

##For now, we only have standard curves for cortisol and cortisone, pull those out
analytes =  analytes %>% filter(unique.data.analyte. == "cortisol" | unique.data.analyte. == "cortisone")

##Choose only the observations greater than the LoD 
stds_cleaned = 
  stds %>%
  filter(amount >= 0.49 & analyte == "cortisol" | amount >= 0.98 & analyte == "cortisone") %>%
  group_by(analyte, amount) %>%
  summarize(area = mean(area))

##Set up an empty data frame to store the standard curve data
std_curves = data.frame(matrix(ncol = 4, nrow = 0))

##Provide column names to the data frame
colnames(std_curves) = c('analyte', 'slope', 'intercept', 'r.squared')

##Calculate the standard curves, using area (later we will use area ratio)
##This loop will read the analyte name from the df "anaytes," and find that data in the df "stds". 
##Then, it will run a regression where y is the area and x is the specified amount.
##The code will save the slope, intercept, and the R^2 value from the linear regression into the new df "std_curves"
for(i in 1:nrow(analytes)){
  std_curves[i,1] = analytes[i,1]
  std_curves[i,2] = summary(lm(area ~ amount, na.action=na.exclude, stds_cleaned[stds_cleaned$analyte == unlist(analytes[i,]),]))$coefficients[2,1]
  std_curves[i,3] = summary(lm(area ~ amount, na.action=na.exclude, stds_cleaned[stds_cleaned$analyte == unlist(analytes[i,]),]))$coefficients[1,1]
  std_curves[i,4] = summary(lm(area ~ amount, na.action=na.exclude, stds_cleaned[stds_cleaned$analyte == unlist(analytes[i,]),]))$r.squared
}


##For now, we are going to extrapolate our values for the deuterated standards from the 
##standard curves of the cortisol and cortisone analytical standards

##Here is a code where we make a duplicate of cortisol standard curve data, and we call it cortisol-d4, 
##And make a duplicate of cortisone standard data and call it cortisol-d8
index = rep(1:nrow(std_curves), 2)
std_curves = std_curves[index,]
std_curves[1,1] = "cortisone_d8"
std_curves[2,1] = "cortisol_d4"


################################################################################
##############################PROCESS DATA######################################
################################################################################

#Calculate observed values for spike solutions and samples

spike_solutions = 
  samples %>%
  filter(sample_name == "spike_solution_1" | sample_name == "spike_solution_2") %>%
  left_join(std_curves, by = "analyte") %>% #bind the standard curves to the df "samples" 
  mutate(spike_solution.ug_L = (area-intercept)/slope) %>% #use the standard curve to calculate the observed value 
  filter(analyte == "cortisol_d4" | analyte == "cortisone_d8") %>% 
  select(sample_name, analyte, spike_solution.ug_L)

samples_cleaned =
  samples %>%
  left_join(std_curves, by = "analyte") %>% #bind the standard curves to the df "samples" 
  mutate(eluate.ug_L = (area-intercept)/slope) %>% #use the standard curve to calculate the observed value 
  filter(analyte == "cortisol_d4" | analyte == "cortisone_d8") %>% #filter out just cortisol_d8 and
  mutate(eluate.ug_L = case_when(analyte == "cortisol_d4" & eluate.ug_L < 1 ~ 0.48, TRUE ~ eluate.ug_L)) %>% #impute LoD
  mutate(eluate.ug_L = case_when(analyte == "cortisone_d8" & eluate.ug_L < 1 ~ 0.98, TRUE ~ eluate.ug_L)) %>% #impute LoD
  right_join(metadata, by = "sample_name") %>% #keep only the samples for which we have metadata
  mutate(observed.concentration.ug_L = eluate.ug_L*volume.eluate.L/volume.processed.L) %>% #calculate the observed concentration 
  left_join(spike_solutions, by = c("spike_solution" = "sample_name", "analyte")) %>% #add in spike solution data 
  mutate(expected.concentration.ug_L = spike_solution.ug_L*0.002/0.2) %>% #calculate the expected concentration for each sample
  mutate(recovery.perc = observed.concentration.ug_L/expected.concentration.ug_L*100) %>%
  mutate(sample_name = stringr::str_replace(sample_name, "positive_control", "positive.control")) %>% 
  mutate(sample_name = stringr::str_replace(sample_name, "negative_control", "negative.control")) %>% 
  separate(sample_name, into = c("treatment", "replicate"), sep = "_") %>%
  select(treatment, replicate, analyte, recovery.perc, observed.concentration.ug_L, expected.concentration.ug_L, 
         area, slope, intercept, r.squared, eluate.ug_L, volume.eluate.L, volume.processed.L,
         spike_solution.ug_L)

#The many-to-many relationship error is due to the fact that each sample has multiple analytes

#save cleaned data 
write_rds(samples_cleaned, "./data/processed_data/my.data.RDS")

################################################################################
###########################PLOT STANDARD CURVES#################################
################################################################################

#Plot the standard curves for Supplemental 
std_curves = 
  std_curves %>%
  mutate(analyte = as.factor(analyte), 
         analyte = fct_recode(analyte, "Cortisol-d4" = "cortisol_d4", "Cortisone-d8" = "cortisone_d8"))


standard_curves_plot = stds_cleaned %>%
  mutate(analyte = as.factor(analyte), 
         analyte = fct_recode(analyte, "Cortisol-d4" = "cortisol", "Cortisone-d8" = "cortisone")) %>%
  ggplot(aes(x = amount, y = area)) + 
  geom_point(aes(color = analyte)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", size = 0.5) + 
  facet_wrap(~analyte) +
  geom_text(data = std_curves, 
            aes(x = 250, y = 5E6, 
                label = paste0("y = ", round(slope),"x +", round(intercept)))) +
  geom_text(data = std_curves, 
            aes(x = 250, y = 4E6, 
                label = paste0(expression(R^2), " = ", round(r.squared, digits = 3)))) +
  xlab("Amount (µg/L)") +
  ylab("Area") +
  theme_prism() + 
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")


tiff('./figures/standard_curves_plot.tiff', units="in", width = 9, height = 6, res=600, compression = 'lzw', pointsize = 12)
plot(standard_curves_plot)
dev.off()