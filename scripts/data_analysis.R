#Importing packages
library(readxl)
library(tidyverse)
library(rio)


#importing data_set and inspection
AMR <- import("data/AMR-genes_amrfinder_22032024 - Sheet1.csv")
AMR <- janitor::clean_names(AMR)
glimpse(AMR)

# data cleaning
AMR <- AMR %>% mutate_at(vars(gene_symbol, element_type, element_subtype,
                              class, subclass, method), as.factor) # converting variables into factors


#exploratory analysis ----
##Summary of class and subclass
AMR_class <- table(AMR$sample_id, AMR$class)
AMR_subclass <- table(AMR$sample_id, AMR$subclass)

#creating stalked bar plots
ggplot(AMR, aes(x= sample_id, fill = class)) +
  geom_bar(position = "stack") + 
  theme(axis.text = element_text(angle = 45 , hjust = 1)) +
  labs(x = "Sample ID", y = "COunt", title = "AMR classes by sample ID")



