# Chromatography multivariate data analysis (MVDA)

#clear the environment to save on driver space and run faster
rm(list = ls())
setwd("/Users/mphomafata/Documents/GitHub/Chromatographic-data")

# import the necessary libraries
req_packages <- c("tidyverse", # to wrangle data frames
                  "glue",
                  "ggplot2", # to plot the spectra
                  "scales", # for formating numerical data in plots
                  "FactoMineR",
                  "factoextra")
lapply(req_packages, require, character.only = TRUE)

# Read-in the data tables
dataset1 <- data.frame(read.csv(file = "/Users/mphomafata/Documents/GitHub/Chromatographic-data/dataset1.csv"))
dataset1$retention_time = round(dataset1$retention_time, 2)
# dataset1 <- dataset1 %>% column_to_rownames(., var = 'retention_time')
# dataset1[-1] = mutate_all(dataset1[-1], function(x) as.numeric(as.character(x)))

dataset2 <- data.frame(read.csv(file = "/Users/mphomafata/Documents/GitHub/Chromatographic-data/dataset2.csv"))
dataset2$retention_time = round(dataset2$retention_time, 2)
# dataset2 <- dataset2 %>% column_to_rownames(., var = 'retention_time')
# dataset2[-1] = mutate_all(dataset2[-1], function(x) as.numeric(as.character(x)))

# dataset3 <- data.frame(read.csv(file = "/Users/mphomafata/Documents/GitHub/Chromatographic-data/dataset3.csv"))
# dataset3$retention_time = round(dataset3$retention_time, 2)
# dataset3 <- dataset3 %>% column_to_rownames(., var = 'retention_time')
# dataset3[-1] = mutate_all(dataset3[-1], function(x) as.numeric(as.character(x)))


merged_dataset = merge(dataset1, dataset2, by= 'retention_time')
merged_dataset <- merged_dataset %>% column_to_rownames(., var = 'retention_time')
merged_dataset[-1] = mutate_all(merged_dataset[-1], function(x) as.numeric(as.character(x)))

# run the MFAs
print("Starting MFA calculation")
mfa_plot <- MFA(
  merged_dataset,
  group = c(length(dataset1)-1,
            length(dataset2)-1),
  type = c(rep("s", 2)),
  ncp = 4,# length(merged_dataset),
  name.group = c('dataset1', 'dataset2'),
  graph = FALSE
)
print("MFA calculation complete!")

mfa_groups <- plot.MFA(mfa_plot, choix = "group")
ggsave(
  filename = "mfa_groups.jpg",
  plot = mfa_groups,
  width = 30,
  height = 15,
  units = 'cm',
  dpi = 300
)
print("MFA individual plot is complete!")

plot.MFA(mfa_plot, choix = "var")