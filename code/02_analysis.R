################################################################################
################################SET UP##########################################
################################################################################

#Load libraries
library(tidyverse)
#install.packages("ggprism")
library(ggprism)
library(ggpubr)
library(rstatix)

#Load in data
data = readRDS("./data/processed_data/my.data.RDS")

################################################################################
##############################VISUALIZE DATA####################################
################################################################################

#Make sure that the data is normally distributed 
shapiro.test(data$recovery.perc) #p-value = 0.1678

#Try to recreate Abby's figure

##Make the dataframe aesthetically pleasing

data = 
  data %>%
  dplyr::filter(treatment != "negative.control") %>%
  mutate(treatment = as.factor(treatment), 
         treatment = fct_recode(treatment, "Control" = "positive.control", 
                                "Centrifuged" = "centrifuged", 
                                "Filtered" = "filtered", 
                                "Unfiltered" = "unfiltered"), 
         analyte = as.factor(analyte), 
         analyte = fct_recode(analyte, "Cortisol-d4" = "cortisol_d4", "Cortisone-d8" = "cortisone_d8"))

##Run the statistics
stat.test = 
  data %>%
  dplyr::filter(treatment != "negative.control") %>%
  dplyr::group_by(analyte) %>%
  rstatix::t_test(recovery.perc ~ treatment,  p.adjust.method = "none") %>%
  add_xy_position(x = "treatment") %>%
  dplyr::filter(p <= 0.05) %>%
  mutate(y.position = c(100, 100, 95))

##Plot the data
fig_1 = data %>%
  left_join(data %>% 
              group_by(treatment, analyte) %>%
              summarize(recovery.mean = mean(recovery.perc), recovery.std = sd(recovery.perc))) %>%
  ggplot() +
  geom_bar(aes(x = reorder(treatment, -recovery.mean), y = recovery.mean, fill = treatment), position = 'identity', stat = "identity") + 
  geom_point(aes(x = reorder(treatment, -recovery.mean), y = recovery.perc)) +
  geom_errorbar(aes(x = treatment, ymin = recovery.mean - recovery.std, ymax = recovery.mean + recovery.std), 
                position = 'identity', stat = "identity", size = 0.01) +
  facet_wrap(~analyte) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)+
  xlab("Treatment") + 
  ylab("Recovery (%)") +
  theme_prism() + 
  scale_fill_prism() + 
  scale_fill_brewer(palette = "Dark2") +
  #theme(axis.text.x = element_text(angle = 45)) + 
  theme(legend.position = "none") 

##Save the plot
tiff('./figures/fig_1.tiff', units="in", width = 9, height = 6, res=600, compression = 'lzw', pointsize = 12)
plot(fig_1)
dev.off()
