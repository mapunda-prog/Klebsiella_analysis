#Importing packages
library(readxl)
library(tidyverse)
library(rio)
library(cowplot)

#importing data_set and inspection
AMR <- import("data/new_amr.csv")
AMR <- janitor::clean_names(AMR)
glimpse(AMR)

# data cleaning
AMR <- AMR %>% mutate_at(vars(gene_symbol, element_type, element_subtype,
                              class, subclass, method), as.factor) # converting variables into factors


#exploratory analysis ----
##Summary of class and subclass
AMR_clean <- AMR %>% select(-c ("protein_identifier", "start", "gene_symbol",
                                          "stop", "strand", "alignment_length", 
                                          "accession_of_closest_sequence", 
                                          "name_of_closest_sequence", "hmm_id", 
                                          "hmm_description")) # selecting columns to use

#separating virulence genes from AMR genes

AMR_data <- AMR_clean %>% filter(element_type == "AMR") 
Virulence_data <- AMR_clean %>% filter(element_type != "AMR") 

#Filter AMR genes were coverage is less the 90%
AMR_data <- AMR_data %>% filter(percent_coverage_of_reference_sequence > 90) 

#creating stalked bar plots
plot_A <- ggplot(AMR_data, aes(x= name, fill = class)) +
  geom_bar(position = "stack") + 
  theme(axis.text = element_text(angle = 45 , hjust = 1)) +
  labs(x = "Sample ID", y = "Count", title = "AMR classes by sample ID") + theme_classic() + 
  scale_fill_manual(name='name',values=c('goldenrod2','darkorange3','red3', 'hotpink3','thistle3','magenta4',                                                                                                              'dodgerblue1','dodgerblue3','slategray3','darkseagreen3','darkgreen',
                                                                            'tan','cornsilk3','bisque4','rosybrown3',                                                                                                                                'grey90'))+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_text(size=10))+
  theme(legend.position = 'bottom',legend.direction = 'horizontal')+
  xlab('Sample IDs')+
  ylab('Available Resistance Genes')

plot_B <- ggplot(AMR, aes(x= sample_id, fill = gene_symbol)) +
  geom_bar(position = "stack") + 
  theme(axis.text = element_text(angle = 45 , hjust = 1)) +
  labs(x = "Sample ID", y = "COunt", title = "AMR classes by sample ID") 

#phenotype data analysis
library(AMR)
AMR_phenotypes <- import("data/phenotype_data.csv")
AMR_phenotypes <- janitor::clean_names(AMR_phenotypes) # cleaning variable names
names(AMR_phenotypes)

#Fomatting data for AMR package
AMR_phenotypes$bacteria <- as.mo(AMR_phenotypes$organism, info= TRUE) #checking for Names

AMR_phenotypes <- AMR_phenotypes %>% 
  mutate_at(vars(amk:cip), as.sir) # converting amr data to SIR

summary(AMR_phenotypes) # quick check of data distribution

#antibiogram by region by phenotype data
Phenotype_abm_by_sample <- antibiogram(AMR_phenotypes, syndromic_group = "organism", minimum = 1)


bact_susceptibility <- AMR_phenotypes %>% 
  summarise(ampicillin = susceptibility(ampicillin),
            chlorampenicol = susceptibility(chloramphenicol),
            Ciprofloxacin = susceptibility(ciprofloxacin),
            trimethoprim = susceptibility(trimethoprim),
            Azithromycin = susceptibility(azithromycin),
            doxycycline = susceptibility(doxycycline),
            Tetracyclines= susceptibility(tetracycline))


bact_Resistance <- AMR_phenotypes %>%
  summarise(ampicillin = resistance(ampicillin),
            chlorampenicol = resistance(chloramphenicol),
            Ciprofloxacin = resistance(ciprofloxacin),
            trimethoprim = resistance(trimethoprim),
            Azithromycin = resistance(azithromycin),
            doxycycline = resistance(doxycycline),
            Tetracyclines= resistance(tetracycline))

bact_susce_region <- AMR_phenotypes %>% group_by(region) %>% 
  summarise(chlorampenicol = susceptibility(chloramphenicol),
            Ciprofloxacin = susceptibility(ciprofloxacin),
            trimethoprim = susceptibility(trimethoprim),
            Azithromycin = susceptibility(azithromycin),
            doxycycline = susceptibility(doxycycline),
            Tetracyclines= susceptibility(tetracycline))


AMR_phenotypes %>%
  select(region, antiseptic, chloramphenicol, ciprofloxacin, trimethoprim, azithromycin, doxycycline,
         tetracycline) %>%
  ggplot_sir(translate_ab = "ab", datalabels = FALSE) 


AMR_phenotypes %>% group_by(region) %>% susceptibility(antiseptic)






