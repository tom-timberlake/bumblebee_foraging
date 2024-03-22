#######################################################################################################################################################
################################## TIMBERLAKE ET AL. - R SCRIPT FOR ALL ANALYSES AND FIGURE PLOTTING  #################################################
#######################################################################################################################################################

library(dplyr)
library(tidyr)
library(lattice)
library("ggplot2")
#install.packages("econullnetr")
library("econullnetr")
library(nlme)
library(lme4)
library(MuMIn)
library(emmeans)
library(multcomp)
library(contrast)
library(plotrix)
library(gridExtra)
library(sjPlot)
library(broom.mixed)
library(vegan)
library(car)
library(readxl)
library(data.table)
library(tibble)
library(glmmTMB)

##################################
##### Loading in data files
##################################

# Paths (Edit these before running)
input.path <- "C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/csv input files/"
output.path <- "C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Results/"

#Input files
barcoding_data_quantitative <- read.csv(file.path(input.path, "All farms_all seasons_relative read abun.csv"))

plant_characterisation_data <- read.csv(file.path(input.path,"barcoding_data_plant_characterisation.csv"))

#Create copy of barcoding data with binary data format (i.e. all non-zero values are converted to 1 to provide presence-absence data only)
barcoding_data_binary <- barcoding_data_quantitative
barcoding_data_binary[] <- lapply(barcoding_data_quantitative, function(x) ifelse(is.numeric(x) & x > 0, 1, x))


##########################################################################################################################
#######  Background: Visualising differences in pollen load composition between species, sites and seasons ###############
##########################################################################################################################

#### Running and Plotting MDS on Barcoding Data 

##Setting data out as a matrix
barcoding_matrix <- barcoding_data_quantitative[,5:182]
rownames(barcoding_matrix) <- barcoding_data_quantitative[,1]
month <- barcoding_data_quantitative[,3]
farm <- barcoding_data_quantitative[,2]
insect <- barcoding_data_quantitative[,4]

####Performing the MDS 
Barcode_MDS=metaMDS(barcoding_matrix,distance = "bray",k=2,trymax=100) # K=The number of reduced dimensions, try max = number of iterations

Barcode_MDS # Shows results of MDS. 

stressplot(Barcode_MDS) # Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (distances between each pair of communities) against their original dissimilarities
# Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions

MDS_coordinates <- Barcode_MDS$points
MDS_dataframe <- as.data.frame(MDS_coordinates)
MDS_dataframe$sample <- rownames(MDS_dataframe)
rownames(MDS_dataframe) <- NULL

#Merge MDS coordinates with original barcoding data
barcoding_nmds_data <- merge(x=barcoding_data_quantitative, y=nmds_dataframe,by="sample",all.x=TRUE)


###### Plotting the MDS ######

### Grouping by month
hull_month <- 
  barcoding_nmds_data %>%
  group_by(month) %>% 
  slice(chull(MDS1, MDS2))

MDS_plot_month <- ggplot(barcoding_nmds_data, aes(x=MDS1, y=MDS2, color=month, shape = month, fill=month)) + 
  geom_point(size=3, show.legend = FALSE) +
  geom_polygon(data = hull_month,
               aes(fill = month,
                   colour = month),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("darkorange", "deepskyblue4", "firebrick")) +
  scale_fill_manual(values = c("darkorange", "deepskyblue4", "firebrick")) +
  theme_classic(base_size = 22, ) +  #sets plot theme - there are various other themes available
  labs(x=bquote("MDS 1"),
       y=bquote("MDS 2")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

### Grouping by farm
hull_farm <- 
  barcoding_nmds_data %>%
  group_by(farm) %>% 
  slice(chull(MDS1, MDS2))

MDS_plot_farm <- ggplot(barcoding_nmds_data, aes(x=MDS1, y=MDS2, color=farm, shape = farm, fill=farm)) + 
  geom_point(size=3, show.legend = FALSE) +
  geom_polygon(data = hull_farm,
               aes(fill = farm,
                   colour = farm),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("orangered1", "turquoise4",  "red4")) +
  scale_fill_manual(values = c("orangered1", "turquoise4",  "red4")) +
  theme_classic(base_size = 22, ) +  #sets plot theme - there are various other themes available
  labs(x=bquote("MDS 1"),
       y=bquote("MDS 2")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

### Grouping by species
hull_species <- 
  barcoding_nmds_data %>%
  group_by(species) %>% 
  slice(chull(MDS1, MDS2))

MDS_plot_species <- ggplot(barcoding_nmds_data, aes(x=MDS1, y=MDS2, color=species, shape = species, fill=species)) + 
  geom_point(size=3, show.legend = FALSE) +
  geom_polygon(data = hull_species,
               aes(fill = species,
                   colour = species),
               alpha = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("darkgreen", "darkorange3",  "red3", "deeppink3", "deepskyblue4")) +
  scale_fill_manual(values = c("darkgreen", "darkorange3",  "red3", "deeppink3", "deepskyblue4")) +
  theme_classic(base_size = 22, ) +  #sets plot theme - there are various other themes available
  labs(x=bquote("MDS 1"),
       y=bquote("MDS 2")) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))


all_MDS_plots <- grid.arrange(MDS_plot_month, MDS_plot_farm, MDS_plot_species, ncol = 3, nrow = 1)

ggsave(plot=all_MDS_plots, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Barcoding_MDS_plots.svg", width=12, height=4, dpi=500, bg="white")
ggsave(plot=all_MDS_plots, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Barcoding_MDS_plots.png", width=12, height=4, dpi=500, bg="white")



############################################################################################################################################
#### Question 1: How many plant species do individual bumblebees utilise and does this vary amongst species, sites and periods of the year?
############################################################################################################################################

barcoding_data_binary_counts <- barcoding_data_binary # Create copy of dataset to maintain structure of original dataset

#Create new variable in binary dataset to calculate number of plant taxa in the pollen load of each insect
binary_data_only <- barcoding_data_binary[, 5:ncol(barcoding_data_binary)] # Subset your data to include only binary variables
row_sums <- rowSums(binary_data_only) # Calculate row sums
barcoding_data_binary_counts$total_plant_taxa <- row_sums # Add row sums to your original dataset

#Summarise means +-Sd for each month
resource_richness_month <- barcoding_data_binary_counts %>%   group_by(month) %>%  
  dplyr::summarise(mean_plant_taxa = mean(total_plant_taxa),
                   sd_plant_taxa = sd(total_plant_taxa),
                   se_plant_taxa = sd(total_plant_taxa) / sqrt(n()))

print(resource_richness_month)

#Summarise means +-Sd for all specimens
resource_richness_all <- barcoding_data_binary_counts %>% 
  dplyr::summarise(median_plant_taxa = median(total_plant_taxa),
                   mean_plant_taxa = mean(total_plant_taxa),
                   sd_plant_taxa = sd(total_plant_taxa),
                   se_plant_taxa = sd(total_plant_taxa) / sqrt(n()))

print(resource_richness_all)

###### Comparing diet breadth between sampling period, species and sites using a GLMM

# Convert variables into factors
barcoding_data_binary_counts$month <- factor(barcoding_data_binary_counts$month)
barcoding_data_binary_counts$species <- factor(barcoding_data_binary_counts$species)
barcoding_data_binary_counts$farm <- factor(barcoding_data_binary_counts$farm)

# Calculate the grand mean of the response variable
grand_mean <- mean(barcoding_data_binary_counts$total_plant_taxa)

# Center the response variable
barcoding_data_binary_counts$total_plant_taxa_centered <- barcoding_data_binary_counts$total_plant_taxa - grand_mean

#Check normality of response variable
hist(barcoding_data_binary_counts$total_plant_taxa_centered)

model_glm <- glm(total_plant_taxa_centered ~ month + species + farm, data = barcoding_data_binary_counts, family = gaussian())

# View the summary of the model
summary(model_glm)

#Post-hoc test for months
month_test <- glht(model_glm, linfct = mcp(month = "Tukey"))
summary(month_test)

######## Plot model coefficients as forest plot ######## 

# Extract coefficients and confidence intervals
coef_ci <- confint(model_glm)

# Prepare data for plotting
coef_data <- data.frame(
  Coefficients = coef(model_glm),
  Lower_CI = coef_ci[, 1],
  Upper_CI = coef_ci[, 2]
)

# Rename the row names of coef_data
new_row_names <- c("Intercept", "Month: July", "Month: Sept", "Species: B.hypnorum", "Species: B.lapidarius", "Species: B.pascuorum", "Species: B.terrestris", "Farm: EM", "Farm: ET")
rownames(coef_data) <- new_row_names

#Remove intercept
coef_data_no_intercept <- coef_data[-which(rownames(coef_data) == "Intercept"), ]

# Plot
resource_richness_forest_plot <- ggplot(coef_data_no_intercept, aes(x = Coefficients, y = rownames(coef_data_no_intercept))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0) +
  labs(title = "",
       x = "Coefficient Estimate",
       y = "Factor Levels") +
  theme_minimal()
  

ggsave(plot=resource_richness_forest_plot, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_richness_coefficients.svg", width=5, height=5, dpi=500, bg="white")
ggsave(plot=resource_richness_forest_plot, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_richness_coefficients.png", width=5, height=5, dpi=500, bg="white")


##########################################
## Plotting resource richness density plot
##########################################

month_colors <- c("April" = "darkorange", 
                  "July" = "deepskyblue4",
                  "Sept" = "firebrick")

# Plot kernel density plot with separate curves for each month
density_plot <- ggplot(barcoding_data_binary_counts, aes(x = total_plant_taxa, color = month, fill = month)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = month_colors) +
  scale_fill_manual(values = month_colors) +
  labs(title = "",
       x = "Total plant taxa in pollen loads",
       y = "Density")

# Calculate mean total_plant_taxa for each month
mean_total_plant_taxa <- aggregate(total_plant_taxa ~ month, data = barcoding_data_binary_counts, FUN = mean)

# Add vertical dashed lines for mean value of total_plant_taxa for each month
density_plot_ann <- density_plot + geom_vline(data = mean_total_plant_taxa, aes(xintercept = total_plant_taxa, color = month), 
                            linetype = "dashed", size = 1) +
  geom_text(data = mean_total_plant_taxa, aes(x = total_plant_taxa, y = 0, label = round(total_plant_taxa, 2)), 
            vjust = 1, hjust = 0, color = "black") +
  theme_minimal()

ggsave(plot=density_plot_ann, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_richness_density.svg", width=5, height=3, dpi=500, bg="white")
ggsave(plot=density_plot_ann, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_richness_density.png", width=5, height=3, dpi=500, bg="white")


##############################################################################################################################
## Question 2: Do bumblebees utilise a greater diversity of plant species than expected under neutral mechanisms?
##############################################################################################################################


###########################################################################
######### Generating the null distribution for pooled species data #########


### Step 1: import the trimmed 'resource' and 'consumer' datasets for each sampling period (earl, mid & late) showing only species for which both resource and consumer data is available
#Resource data shows the estimates of floral resource availability in the landscape for each species (subdivided by farm and season and trimmed to match consumer data)
#Consumer data shows the quantitative pollen barcoding data (subdivided by farm and season and trimmed to match resource data)


##Consumer data
pollinators_early_pooled <- data.table(read_excel(file.path(input.path, "Consumer_data.xlsx"), sheet = "early", na = c("", "---", NA))) 
pollinators_mid_pooled <- data.table(read_excel(file.path(input.path, "Consumer_data.xlsx"), sheet = "mid", na = c("", "---", NA))) 
pollinators_late_pooled <- data.table(read_excel(file.path(input.path, "Consumer_data.xlsx"), sheet = "late", na = c("", "---", NA))) 
                                              

##Resource data
plants_early_FU <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "early_FU", na = c("", "---", NA)))  
plants_early_nectar <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "early_nectar", na = c("", "---", NA)))   
plants_early_pollen <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "early_pollen", na = c("", "---", NA)))  

plants_mid_FU <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "mid_FU", na = c("", "---", NA)))  
plants_mid_nectar <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "mid_nectar", na = c("", "---", NA)))   
plants_mid_pollen <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "mid_pollen", na = c("", "---", NA)))  

plants_late_FU <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "late_FU", na = c("", "---", NA)))  
plants_late_nectar <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "late_nectar", na = c("", "---", NA)))   
plants_late_pollen <- data.table(read_excel(file.path(input.path, "Resource_data.xlsx"), sheet = "late_pollen", na = c("", "---", NA)))  

### Step 2: Generate null networks of resource used based on relative availability of each floral resource measure (floral abundance, nectar and pollen)
#The model redistributes interactions between plants and pollinators and reallocated them in direct proportion to reseource availability of each plant species
#NOTE that the data is subdivided by field site and sampling month so the redistribution is done within these sub-categories of the data
#NOTE that a separate null model is generated for each floral resource measure (nectar, pollen and floral abundance)


#Set seed
set.seed(1234)

###### Early Season ######

##Floral unit data
null.net_early_pooled_FU<-generate_null_net(pollinators_early_pooled[, 3:20], plants_early_FU[, 2:18], sims=10000,
                                            data.type = "quantities",
                                            c.samples=pollinators_early_pooled[, 2],
                                            r.samples=plants_early_FU[, 1])

##Nectar data
null.net_early_pooled_nectar<-generate_null_net(pollinators_early_pooled[, 3:20], plants_early_nectar[, 2:18], sims=10000,
                                                data.type = "quantities",
                                                c.samples=pollinators_early_pooled[, 2],
                                                r.samples=plants_early_nectar[, 1])
##Pollen data
null.net_early_pooled_pollen<-generate_null_net(pollinators_early_pooled[, 3:20], plants_early_pollen[, 2:18], sims=10000,
                                                data.type = "quantities",
                                                c.samples=pollinators_early_pooled[, 2],
                                                r.samples=plants_early_pollen[, 1])

###### Mid Season ######

##Floral unit data
null.net_mid_pooled_FU<-generate_null_net(pollinators_mid_pooled[, 3:28], plants_mid_FU[, 2:26], sims=10000,
                                          data.type = "quantities",
                                          c.samples=pollinators_mid_pooled[, 2],
                                          r.samples=plants_mid_FU[, 1])

##Nectar data
null.net_mid_pooled_nectar<-generate_null_net(pollinators_mid_pooled[, 3:28], plants_mid_nectar[, 2:26], sims=10000,
                                              data.type = "quantities",
                                              c.samples=pollinators_mid_pooled[, 2],
                                              r.samples=plants_mid_nectar[, 1])
##Pollen data
null.net_mid_pooled_pollen<-generate_null_net(pollinators_mid_pooled[, 3:28], plants_mid_pollen[, 2:26], sims=10000,
                                              data.type = "quantities",
                                              c.samples=pollinators_mid_pooled[, 2],
                                              r.samples=plants_mid_pollen[, 1])

###### Late Season ######

##Floral unit data
null.net_late_pooled_FU<-generate_null_net(pollinators_late_pooled[, 3:23], plants_late_FU[, 2:21], sims=10000,
                                           data.type = "quantities",
                                           c.samples=pollinators_late_pooled[, 2],
                                           r.samples=plants_late_FU[, 1])

##Nectar data
null.net_late_pooled_nectar<-generate_null_net(pollinators_late_pooled[, 3:23], plants_late_nectar[, 2:21], sims=10000,
                                               data.type = "quantities",
                                               c.samples=pollinators_late_pooled[, 2],
                                               r.samples=plants_late_nectar[, 1])
##Pollen data
null.net_late_pooled_pollen<-generate_null_net(pollinators_late_pooled[, 3:23], plants_late_pollen[, 2:21], sims=10000,
                                               data.type = "quantities",
                                               c.samples=pollinators_late_pooled[, 2],
                                               r.samples=plants_late_pollen[, 1])


### Step 3: Calculate Shannon diversity of the observed network as well as each of the 3 null networks (based on floral abundance, nectar and pollen)
#NOTE: Shannon diversity is calculated separately for each sampling period (early, mid and late) and for each null network, as well as observed network 

#Early season
net.stats_early_FU <- bipartite_stats(null.net_early_pooled_FU, index.type = "networklevel",
                                          indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

net.stats_early_nectar <- bipartite_stats(null.net_early_pooled_nectar, index.type = "networklevel",
                                          indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

net.stats_early_pollen <- bipartite_stats(null.net_early_pooled_pollen, index.type = "networklevel",
                                          indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

#mid season
net.stats_mid_FU <- bipartite_stats(null.net_mid_pooled_FU, index.type = "networklevel",
                                        indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

net.stats_mid_nectar <- bipartite_stats(null.net_mid_pooled_nectar, index.type = "networklevel",
                                        indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

net.stats_mid_pollen <- bipartite_stats(null.net_mid_pooled_pollen, index.type = "networklevel",
                                        indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

#late season
net.stats_late_FU <- bipartite_stats(null.net_late_pooled_FU, index.type = "networklevel",
                                         indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

net.stats_late_nectar <- bipartite_stats(null.net_late_pooled_nectar, index.type = "networklevel",
                                         indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

net.stats_late_pollen <- bipartite_stats(null.net_late_pooled_pollen, index.type = "networklevel",
                                         indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")

### Step 4: Restructure and bind all results so that they can be plotted together

#Early results
net.stats_early_FU$null_model <- "FU"
net.stats_early_FU$statistic <- rownames(net.stats_early_FU)
rownames(net.stats_early_FU) <- NULL

net.stats_early_nectar$null_model <- "nectar"
net.stats_early_nectar$statistic <- rownames(net.stats_early_nectar)
rownames(net.stats_early_nectar) <- NULL

net.stats_early_pollen$null_model <- "pollen"
net.stats_early_pollen$statistic <- rownames(net.stats_early_pollen)
rownames(net.stats_early_pollen) <- NULL

early_results <- rbind(net.stats_early_FU, net.stats_early_nectar, net.stats_early_pollen)
early_results$season <- "early"

#Mid results
net.stats_mid_FU$null_model <- "FU"
net.stats_mid_FU$statistic <- rownames(net.stats_mid_FU)
rownames(net.stats_mid_FU) <- NULL

net.stats_mid_nectar$null_model <- "nectar"
net.stats_mid_nectar$statistic <- rownames(net.stats_mid_nectar)
rownames(net.stats_mid_nectar) <- NULL

net.stats_mid_pollen$null_model <- "pollen"
net.stats_mid_pollen$statistic <- rownames(net.stats_mid_pollen)
rownames(net.stats_mid_pollen) <- NULL

mid_results <- rbind(net.stats_mid_FU, net.stats_mid_nectar, net.stats_mid_pollen)
mid_results$season <- "mid"

#Late results
net.stats_late_FU$null_model <- "FU"
net.stats_late_FU$statistic <- rownames(net.stats_late_FU)
rownames(net.stats_late_FU) <- NULL

net.stats_late_nectar$null_model <- "nectar"
net.stats_late_nectar$statistic <- rownames(net.stats_late_nectar)
rownames(net.stats_late_nectar) <- NULL

net.stats_late_pollen$null_model <- "pollen"
net.stats_late_pollen$statistic <- rownames(net.stats_late_pollen)
rownames(net.stats_late_pollen) <- NULL

late_results <- rbind(net.stats_late_FU, net.stats_late_nectar, net.stats_late_pollen)
late_results$season <- "late"

all_results <- rbind(early_results, mid_results, late_results)

#Shannon results only
shannon_results <- all_results[all_results$statistic == 'Shannon diversity', ]

shannon_results <- shannon_results %>%
  mutate(season = case_when(season == "early" ~ "April",
                            season == "mid" ~ "July",
                            season == "late" ~ "Sept",
                            TRUE ~ season))
    

### Step 5: Plot the results

# Convert 'season' to a factor with the defined order
season_order <- c("Sept", "July", "April")
shannon_results$season <- factor(shannon_results$season, levels = season_order)


diversity_plot <- ggplot(shannon_results, aes(x = Null, y = season, xmin = Lower.CL, xmax = Upper.CL)) +
  geom_point(aes(color = null_model), size = 3, position = position_dodge(width = 0.2)) +
  geom_errorbarh(aes(xmin = Lower.CL, xmax = Upper.CL, color = null_model), height = 0.1, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Observed), size = 5, shape = 9, stroke = 1) +  # Add a point for observed values
  labs(title = "",
       x = "Shannon diversity",
       y = "",
       color = "null_model") +
  scale_color_manual(values = c("FU" = "#EB8252", "nectar" = "#8FA33F", "pollen" = "#008080")) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 10),  # Adjust x-axis label size
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"))  # Adjust y-axis label size

ggsave(plot=diversity_plot, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Diversity_plot.svg", width=4, height=5, dpi=500, bg="white")
ggsave(plot=diversity_plot, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Diversity_plot.png", width=4, height=5, dpi=500, bg="white")


#Save the results as a csv file
#write.csv(all_results,file = file.path(output.path, "Network-metrics_all_results.csv"))


##############################################################################################################################
## Question 3: Which plant species do bumblebees preferentially utilise?
##############################################################################################################################

### Step 1: Test whether each plant species is utilised significantly more or less than expected from the null models
#NOTE these tests are run separately for each null model (floral abundance, nectar and pollen) and for each sampling period (earl, mid and late)
#Results are tested at the 95% confidence level - if observed values fall outside the 95% confidence interval of the null model, they are considered significantly different

#FU data
Early_interactions_test_pooled_FU <- test_interactions(null.net_early_pooled_FU, 0.95)
Early_interactions_test_pooled_FU$null_model <- "FU"
Early_interactions_test_pooled_FU$season <- "April"

Mid_interactions_test_pooled_FU <- test_interactions(null.net_mid_pooled_FU, 0.95)
Mid_interactions_test_pooled_FU$null_model <- "FU"
Mid_interactions_test_pooled_FU$season <- "July"

Late_interactions_test_pooled_FU <- test_interactions(null.net_late_pooled_FU, 0.95)
Late_interactions_test_pooled_FU$null_model <- "FU"
Late_interactions_test_pooled_FU$season <- "Sept"


#Nectar data
Early_interactions_test_pooled_nectar <- test_interactions(null.net_early_pooled_nectar, 0.95)
Early_interactions_test_pooled_nectar$null_model <- "nectar"
Early_interactions_test_pooled_nectar$season <- "April"

Mid_interactions_test_pooled_nectar <- test_interactions(null.net_mid_pooled_nectar, 0.95)
Mid_interactions_test_pooled_nectar$null_model <- "nectar"
Mid_interactions_test_pooled_nectar$season <- "July"

Late_interactions_test_pooled_nectar <- test_interactions(null.net_late_pooled_nectar, 0.95)
Late_interactions_test_pooled_nectar$null_model <- "nectar"
Late_interactions_test_pooled_nectar$season <- "Sept"


#Pollen data
Early_interactions_test_pooled_pollen <- test_interactions(null.net_early_pooled_pollen, 0.95)
Early_interactions_test_pooled_pollen$null_model <- "pollen"
Early_interactions_test_pooled_pollen$season <- "April"

Mid_interactions_test_pooled_pollen <- test_interactions(null.net_mid_pooled_pollen, 0.95)
Mid_interactions_test_pooled_pollen$null_model <- "pollen"
Mid_interactions_test_pooled_pollen$season <- "July"

Late_interactions_test_pooled_pollen <- test_interactions(null.net_late_pooled_pollen, 0.95)
Late_interactions_test_pooled_pollen$null_model <- "pollen"
Late_interactions_test_pooled_pollen$season <- "Sept"

#Bind together the results for each season
All_results_April <- rbind(Early_interactions_test_pooled_FU, Early_interactions_test_pooled_nectar, Early_interactions_test_pooled_pollen)
All_results_July <- rbind(Mid_interactions_test_pooled_FU, Mid_interactions_test_pooled_nectar, Mid_interactions_test_pooled_pollen)
All_results_Sept <- rbind(Late_interactions_test_pooled_FU, Late_interactions_test_pooled_nectar, Late_interactions_test_pooled_pollen)

#Bind together all results
Resource_use_all_results <- rbind(All_results_April, All_results_July, All_results_Sept)

#Save the results as a csv file
#write.csv(Resource_use_all_results,file = file.path(output.path, "Resource_use_all_results.csv"))


##Create figure showing all results - restructure data and then plot all results


##################
## April results
##################

# Reorder Plant species based on their 'Observed' resource use values
# Group by 'Resource' and calculate the mean of 'Observed'
resource_summary_April <- All_results_April %>%
  group_by(Resource) %>%
  summarise(mean_observed = mean(Observed, na.rm = TRUE)) %>%
  arrange(mean_observed)  # Arrange by descending mean_observed

# Reorder 'Resource' based on the mean_observed
All_results_April$Resource <- factor(All_results_April$Resource, levels = resource_summary_April$Resource)

resource_use_April <- ggplot(All_results_April, aes(x = Null, y = Resource, xmin = Lower.95.CL, xmax = Upper.95.CL)) +
  geom_point(aes(color = null_model), size = 3, position = position_dodge(width = 0.2)) +
  geom_errorbarh(aes(xmin = Lower.95.CL, xmax = Upper.95.CL, color = null_model), height = 0.1, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Observed), size = 3, shape = 9, stroke = 0.6) +  # Add a point for observed values
  labs(title = "April",
       x = "Number of bees utilising this resource",
       y = "",
       color = "null_model") +
  scale_color_manual(values = c("FU" = "#EB8252", "nectar" = "#8FA33F", "pollen" = "#008080")) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 12),  # Adjust x-axis label size
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold.italic"),
        legend.text = element_text(size = 12))  # Adjust y-axis label size

ggsave(plot=resource_use_April, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_use_plot_April_P-A_data.svg", width=6, height=7, dpi=500, bg="white")
ggsave(plot=resource_use_April, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_use_plot_April_P-A_data.png", width=6, height=7, dpi=500, bg="white")


##################
## July results
##################

# Reorder Plant species based on their 'Observed' resource use values
# Group by 'Resource' and calculate the mean of 'Observed'
resource_summary_July <- All_results_July %>%
  group_by(Resource) %>%
  summarise(mean_observed = mean(Observed, na.rm = TRUE)) %>%
  arrange(mean_observed)  # Arrange by descending mean_observed

# Reorder 'Resource' based on the mean_observed
All_results_July$Resource <- factor(All_results_July$Resource, levels = resource_summary_July$Resource)

resource_use_July <- ggplot(All_results_July, aes(x = Null, y = Resource, xmin = Lower.95.CL, xmax = Upper.95.CL)) +
  geom_point(aes(color = null_model), size = 3, position = position_dodge(width = 0.2)) +
  geom_errorbarh(aes(xmin = Lower.95.CL, xmax = Upper.95.CL, color = null_model), height = 0.1, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Observed), size = 3, shape = 9, stroke = 0.6) +  # Add a point for observed values
  labs(title = "July",
       x = "Number of bees utilising this resource",
       y = "",
       color = "null_model") +
  scale_color_manual(values = c("FU" = "#EB8252", "nectar" = "#8FA33F", "pollen" = "#008080")) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 12),  # Adjust x-axis label size
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold.italic"),
        legend.text = element_text(size = 12))  # Adjust y-axis label size

ggsave(plot=resource_use_July, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_use_plot_July_P-A_data.svg", width=6, height=7, dpi=500, bg="white")
ggsave(plot=resource_use_July, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_use_plot_July_P-A_data.png", width=6, height=7, dpi=500, bg="white")

##################
## Sept results
##################

# Reorder Plant species based on their 'Observed' resource use values
# Group by 'Resource' and calculate the mean of 'Observed'
resource_summary_Sept <- All_results_Sept %>%
  group_by(Resource) %>%
  summarise(mean_observed = mean(Observed, na.rm = TRUE)) %>%
  arrange(mean_observed)  # Arrange by descending mean_observed

# Reorder 'Resource' based on the mean_observed
All_results_Sept$Resource <- factor(All_results_Sept$Resource, levels = resource_summary_Sept$Resource)

resource_use_Sept <- ggplot(All_results_Sept, aes(x = Null, y = Resource, xmin = Lower.95.CL, xmax = Upper.95.CL)) +
  geom_point(aes(color = null_model), size = 3, position = position_dodge(width = 0.2)) +
  geom_errorbarh(aes(xmin = Lower.95.CL, xmax = Upper.95.CL, color = null_model), height = 0.1, position = position_dodge(width = 0.2)) +
  geom_point(aes(x = Observed), size = 3, shape = 9, stroke = 0.6) +  # Add a point for observed values
  labs(title = "Sept",
       x = "Number of bees utilising this resource",
       y = "",
       color = "null_model") +
  scale_color_manual(values = c("FU" = "#EB8252", "nectar" = "#8FA33F", "pollen" = "#008080")) +
  theme_minimal()+
  theme(axis.text.x = element_text(size = 12),  # Adjust x-axis label size
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold.italic"),
        legend.text = element_text(size = 12))  # Adjust y-axis label size

ggsave(plot=resource_use_Sept, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_use_plot_Sept_P-A_data.svg", width=6, height=7, dpi=500, bg="white")
ggsave(plot=resource_use_Sept, filename="C:/Users/tt15117/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/AA_Paper writing/Re-analysis/Plots/Resource_use_plot_Sept_P-A_data.png", width=6, height=7, dpi=500, bg="white")



####################################################
## Characterisation of plant types in pollen loads
####################################################

#We will be using the plant characterisation dataset which categorises each plant recorded in the barcoding data as either garden, native or NA, as well as separate characterisation of the plant's growth form/habit

# Step 1: Summarise the plant characterisation data to calculate the number of garden vs non-garden plants recorded in the pollen load of each bee sample
garden_plant_counts <- plant_characterisation_data %>%
  dplyr:: group_by(garden_wild_invasive) %>%
  dplyr::summarize(across(where(is.numeric), sum, na.rm = TRUE)) 

#Pivot dataframe
garden_data_long <- as.data.frame(t(garden_plant_counts))
colnames(garden_data_long) <- garden_data_long[1, ]  # Set the column names to the values of the first row
garden_data_long <- garden_data_long[-1, ]  # Remove the first row which is now the column names
garden_data_long$sample <- rownames(garden_data_long)
rownames(garden_data_long) <- NULL

numeric_columns <- c("garden", "invasive", "unknown", "wild")  # Specify columns to be numeric
garden_data_long[numeric_columns] <- lapply(garden_data_long[numeric_columns], as.numeric)


#Calculate proportion of garden plants in pollen basket of each bee
garden_data_long$proportion_garden <- garden_data_long$garden / rowSums(garden_data_long[, -which(names(garden_data_long) == "sample")])

#Merge in other sample information (season, site, species etc.)

# Replace '.' with '-' in the 'sample' column to get in matching format with barcoding data
garden_data_long$sample <- gsub("\\.", "-", garden_data_long$sample)

#Merge the two dataframes
barcoding_data_garden_info <- merge(x=barcoding_data_binary , y=garden_data_long,by="sample",all.x=TRUE)

#Summarise data
garden_plant_summary <- barcoding_data_garden_info %>% group_by(month) %>%
    dplyr::summarise(mean_garden_prop = mean(proportion_garden, na.rm = TRUE),
                   sd_garden_prop = sd(proportion_garden, na.rm = TRUE),
                   se_garden_prop = sd(proportion_garden, na.rm = TRUE) / sqrt(n()))

# Calculate the grand mean of garden plant proportions
mean(barcoding_data_garden_info$proportion_garden, na.rm = TRUE) # Report grant mean of garden plants in pollen loads
sd(barcoding_data_garden_info$proportion_garden, na.rm = TRUE) # Report grant mean of garden plants in pollen loads


#######################################################################################################################
############# Floral resource analysis - modeling smooth, non-linear trend in floral abundance over time  #############
#######################################################################################################################


#The code below provides an example of the model structure and selection process used to model the phenology of nectar for each plant species on each farm
#Note that this process was completed in Timberlake et al. 2019. J.Appl.Ecol and the model predictions from this paper (same farms and sampling year) were used to provide resource data to inform the null models in this study
#The code below shows the model selection process for just one species (Rubus fruticosus) on one farm (Eastwood Manor)

Species_nectar <- read.csv(file.path(input.path, "All farms_2017_nectar by species_per km2_grams.csv"))

# Using week as a non-linear covariate (continuous variable)
Species_nectar$day <- as.numeric(Species_nectar$day)

#Split the data by farm
Species_nectar_Birches<-Species_nectar[Species_nectar$farm == "Birches",]
Species_nectar_Eastwood<-Species_nectar[Species_nectar$farm == "Eastwood",]
Species_nectar_Elmtree<-Species_nectar[Species_nectar$farm == "Elmtree",]

summary(Species_nectar_Eastwood) #summarise data

str(Species_nectar_Eastwood) #check data structure


#### GAM model - one model constructed for each k value

library(mgcv)

Rubus.fruticosus_model1 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 4), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model2 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 5), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model3 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 6), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model4 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 7), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model5 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 8), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model6 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 9), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model7 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 10), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model8 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 11), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model9 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 12), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model10 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 13), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model11 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 14), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model12 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 15), family = Gamma(link = "log"),data = Species_nectar_Eastwood)
Rubus.fruticosus_model13 <- gam(Rubus.fruticosus + 0.0001 ~s(day, fx = T, k= 4), family = Gamma(link = "identity"),data = Species_nectar_Eastwood)


AIC(Rubus.fruticosus_model1)
AIC(Rubus.fruticosus_model2)
AIC(Rubus.fruticosus_model3)
AIC(Rubus.fruticosus_model4)
AIC(Rubus.fruticosus_model5)
AIC(Rubus.fruticosus_model6)
AIC(Rubus.fruticosus_model7)
AIC(Rubus.fruticosus_model8)
AIC(Rubus.fruticosus_model9)
AIC(Rubus.fruticosus_model10)
AIC(Rubus.fruticosus_model11)
AIC(Rubus.fruticosus_model12)
AIC(Rubus.fruticosus_model13)


#### Plotting the data with model overlain
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))

plot(jitter(Rubus.fruticosus, amount = 0.5) ~ day, # plotting daily sugar of given species per day with a minor adjustment Jitter shifts points slightly so that they don't overlap and get hidden
     data = Species_nectar,
     ylab = "Sugar/km2/day (grams)",    #labels y axis      
     xlab = "Date", # labels x axis
     cex.lab = 2, # size of y and x axis labels
     las = 1,  # Ensure that all the axis labels are horizontal
     
     col = "red", # colour the dots in different colours by habitat type, based on the alphabetical order of habitats
     pch = 19, # changes the type of dot for data
     xlim=c(1,260),
     ylim=c(0,1000),#fixes the range of x values to set the x axis
     xaxt="n",
     bty = "L") # changes the box type to an L rather than square

axis(side = 1, at = c(1,32,62,93,123,154,185,215,246), 
     labels = c("March","April","May","June","July","Aug","Sept","Oct","Nov"),las=1)

##Using pred function to plot model predictions over the data
pdat <- expand.grid(day = seq(-1,260,1)) 
pred <- predict (Rubus.fruticosus_model5, newdata = pdat, na.rm = T, type= "response", se.fit = TRUE)
predframe <- data.frame (pdat,level=0, preds = pred$fit, se = pred$se.fit)
lines(predframe$preds+predframe$se~predframe$day, lwd=1, lty=2, col="red")
lines(predframe$preds~predframe$day, lwd=3, col="red")
lines(predframe$preds-predframe$se~predframe$day, lwd=1, lty=2, col="red")


write.csv(pred, file="Chosen_file_name.csv") #write daily predictions of nectar to csv file


########################################################################################################################################
#######################################################   END OF SCRIPT  ###############################################################
########################################################################################################################################
