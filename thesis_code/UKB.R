# 
####June.9th-June15th####
# # load the haven package
# library(haven)
# 
# # read the .dta file
# data <- read_dta("UKB_<application>_extract.dta")
# 
# # view the data
# View(data)
# 
# #write csv 
# write.csv(data,"UKB.csv")

data<-read.csv("UKB.csv")

#select protein data into new csv file 
protein <- data[, 34:761]

####1. Find NA cutoff####
#count how many NAs, mean, SD of each column 
summary_protein <- data.frame(
  variable = colnames(protein),  # column names of the protein,
  mean = sapply(protein, function(x) mean(x, na.rm = TRUE)),  # mean 
  sd = sapply(protein, function(x) sd(x, na.rm = TRUE)),  # sd
  n_NA = sapply(protein, function(x) sum(is.na(x))),  # amount of NAs
  na_percent = sapply(protein, function(x) sum(is.na(x)) / length(x)*100) #percentage of NAs 
)

####1.1draw a histogram to see the distribution####
NA_h<-hist(summary_protein$na_percent, breaks = 100, main = "NA Percentage Distribution", xlab = "NA Percentage")

#only one data is apart from around 25% 
#narrow down the area
hist(summary_protein$na_percent, 
     breaks = 200,
     main = "NA Percentage Distribution",
     xlab = "NA Percentage",
     xlim = c(0,30),
     ylim = c(0, max(NA_h$counts) * 1.15))
# add labels (counts shown above each bar)

keep <- NA_h$counts > 0 #get rid of all the zeros

text(x = NA_h$mids[keep],        # middle position of each bar
     y = NA_h$counts[keep],      # height of each bar (count)
     labels = NA_h$counts[keep], # text to display
     pos = 3,           # place the text above the bar
     cex = 0.7)         # font size (adjustable)

#only one data is around 25 
#keep shrinking the targeting area 
hist(summary_protein$na_percent, 
     breaks = 100,
     main = "NA Percentage Distribution",
     xlab = "NA Percentage",
     xlim = c(0,20),
     ylim = c(0, max(NA_h$counts) * 1.15)) 
# add labels (counts shown above each bar)

keep <- NA_h$counts > 0 #get rid of all the zeros

text(x = NA_h$mids[keep],        # middle position of each bar
     y = NA_h$counts[keep],      # height of each bar (count)
     labels = NA_h$counts[keep], # text to display
     pos = 3,           # place the text above the bar
     cex = 0.7)         # font size (adjustable)

#cutoff at 5%
hist(summary_protein$na_percent,
     breaks = 100,
     main = "NA Percentage Distribution",
     xlab = "NA Percentage",
     xlim = c(0, 25),
     ylim = c(0, max(NA_h$counts) * 1.15))
abline(v = 6, col = "red", lty = 2, lwd = 2)
keep <- NA_h$counts > 0 #get rid of all the zeros

text(x = NA_h$mids[keep],        # middle position of each bar
     y = NA_h$counts[keep],      # height of each bar (count)
     labels = NA_h$counts[keep], # text to display
     pos = 3,           # place the text above the bar
     cex = 0.7)         # font size (adjustable)

#decide the cutoff 

table(summary_protein$na_percent <= 6)
370/(370+358)*100
#percentage of proteins kept at this cutoff

table(summary_protein$na_percent <= 5)

252/(476+252)*100
#percentage of proteins kept at this cutoff

coverage_proteins <- summary_protein[summary_protein$na_percent <= 5, ]
 
####June16th-June22nd####
####1.2cut off go with 10% ####

#next task 
#which of these protin has 10% missing has spearman's correlated with PRS? correlation (each protein with the PRS)
#what p value threhold define significant is? 0.05? or something else 
#writing: intro & method 
#polygenic score

# keep proteins with <= 10% missing values
selected_proteins <- summary_protein$variable[summary_protein$na_percent <= 10]

# select the chosen protein columns from the protein data
filtered_protein <- protein[, selected_proteins]

library(dplyr)
prs <- data %>% select(prs_bp)
#combine prs with filtered_protein 
filtered_protein <- cbind(filtered_protein,prs)

####2.1 Spearman rank correlation ####

# create an empty matrix to store the correlation coefficients (rho)
cor_matrix <- matrix(NA, nrow = ncol(filtered_protein), ncol = ncol(prs))
rownames(cor_matrix) <- colnames(filtered_protein)  # row names are the protein names
colnames(cor_matrix) <- colnames(prs)               # column names are the PRS names

# create an empty matrix to store the p-values
pval_matrix <- matrix(NA, nrow = ncol(filtered_protein), ncol = ncol(prs))
rownames(pval_matrix) <- colnames(filtered_protein)
colnames(pval_matrix) <- colnames(prs)

# this double for loop goes through every protein x PRS combination
for (i in seq_along(colnames(filtered_protein))) {
  for (j in seq_along(colnames(prs))) {
    
    # Spearman correlation test between protein i and PRS j
    result <- cor.test(
      filtered_protein[[i]],     # expression vector of protein i
      prs[[j]],                  # score vector of PRS j
      method = "spearman",       # use Spearman rank correlation
      use = "pairwise.complete.obs"  # if there are NAs, use only complete pairs
    )
    
    # extract the correlation coefficient rho from the result
    cor_matrix[i, j] <- result$estimate
    
    # extract the p-value from the result
    pval_matrix[i, j] <- result$p.value
  }
}

####data observation####

# 1. turn the two matrices into long format (one row per protein x PRS pair)

# convert the matrices to data frames
cor_df <- as.data.frame(cor_matrix)
pval_df <- as.data.frame(pval_matrix)

# turn both matrices into long format (one row per protein x PRS pair)
library(tidyr)
cor_long <- pivot_longer(cor_df, cols = everything(), names_to = "PRS", values_to = "rho", names_transform = list(PRS = as.character))
cor_long$Protein <- rep(rownames(cor_df), times = ncol(cor_df))

pval_long <- pivot_longer(pval_df, cols = everything(), names_to = "PRS", values_to = "p_value", names_transform = list(PRS = as.character))
pval_long$Protein <- rep(rownames(pval_df), times = ncol(pval_df))

# 2. merge rho and p-value
spearman <- merge(cor_long, pval_long, by = c("Protein", "PRS"))

# 3. save as a CSV file
write.csv(spearman, "spearman_correlation_results.csv", row.names = FALSE)

#find the data range 
min_vals <- sapply(filtered_protein, function(x) min(x, na.rm = TRUE))  # minimum of each column
max_vals <- sapply(filtered_protein, function(x) max(x, na.rm = TRUE))  # maximum of each column
range_table <- data.frame(
  filtered_protein = colnames(filtered_protein),
  min = min_vals,
  max = max_vals
)

# the whole row with the maximum value
range_table[which.max(range_table$max), ]

# the whole row with the minimum value
range_table[which.min(range_table$min), ]

####2.2. Spearman cutoff####
#without any cutoff
# step 1: save the histogram object without plotting
p_hist <- hist(spearman$p_value, breaks = 200, plot = FALSE)

# step 2: use it to draw the plot and set ylim
hist(spearman$p_value,
     breaks = 200,
     main = "P value Distribution",
     xlab = "P value",
     col = "gray",
     ylim = c(0, max(p_hist$counts) * 1.15))  # set the y axis from the p_hist object

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour

#cutoff at 0.2 
p_hist <- hist(spearman$p_value,
               breaks = 200,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.2),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour

#cutoff at 0.05 
p_hist <- hist(spearman$p_value,
               breaks = 200,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.05),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour

#cutoff at 0.01
p_hist <- hist(spearman$p_value,
               breaks = 200,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.01),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour
#probably cutoff at 0.005

#use rho value to see the effect size 
significant <- spearman[spearman$p_value < 0.005, ]
summary(significant$rho)

#histogram plot 
# hist(significant$rho,
#      breaks = 30,
#      main = "Distribution of Spearman Rho (p < 0.005)",
#      xlab = "Spearman Rho",
#      col = "lightblue")

#density plot 
plot(density(significant$rho),
     main = "Density Plot of Spearman Rho (p < 0.005)",
     xlab = "Spearman Rho",
     col = "steelblue",
     lwd = 2)
abline(v = 0, col = "red", lty = 2)  # reference line

#cutoff at 0.2

plot(
  density(significant$rho,  from = -0.2, to = 0.2),  # compute the density and set its range
  main = "Density Plot of Spearman Rho (p < 0.005)",  # plot title
  xlab = "Spearman Rho",      # x axis title
  col = "steelblue",          # line colour
  lwd = 2                     # line width
)

# add a red vertical line at rho = 0 to separate negative from positive correlations
abline(v = 0, col = "red", lty = 2)

#check the direction of the correlations

#p value divide by how many proteins 

####June 23rd-June 30th####

####Use Bonferonni correction to decide thes significant protein cutoff####
#significance threshold = 0.05/365
0.05/363
#Bonferroni-corrected significance level
significant <- spearman[spearman$p_value < 0.0001 & spearman$p_value !=0 & spearman$p_value !=1,]


plot(density(significant$rho),
     main = "Density Plot of Spearman Rho (p < 0.0001)",
     xlab = "Spearman Rho",
     col = "steelblue",
     lwd = 2)
abline(v = 0, col = "red", lty = 2)  # reference line

#negative and positive both are good 

plot(density(significant$rho, from = -0.2, to = 0.2),
     main = "Density Plot of Spearman Rho (p < 0.0001)",
     xlab = "Spearman Rho",
     col = "steelblue",
     lwd = 2)
abline(v = 0, col = "red", lty = 2)  # reference line

#check the direction of the correlations
#

####3.regression analysis #####
#linear regression, what effect prs has on the protein level?
#age, sex (sex plink), population structure (PGC1-PGC8), genotyping chips (chip)
#multivariate regression analysis  (prs and all the other covariates)
#make sure what we finding is reflecting genetic and protein level, not just across the population 

#Task 
#take the significant proteins
#run one multivariate linear regression model for each protein 
#protein ~ PRS-BD-std age sex pgc1 pgc2 pgc3 pgc4 pgc5 pgc6 pgc7 pgc8 chip

#presentation + report (two week before ddl :D)

#opportunity for the conference. :D 

####Basic model####

# first take the protein names to select (assuming the column is called protein)
selected_proteins2 <- significant$Protein

# then use them to select columns from the protein data
filtered_protein2 <- protein[, selected_proteins2]

covaritate<- data %>% select(prs_bp_std,age_0_0,sex_plink,PGC1,PGC2,PGC3,PGC4,PGC5,PGC6,PGC7,PGC8,chip)
#combine prs with filtered_protein 
filtered_protein2 <- cbind(filtered_protein2,covaritate)

# take the column names of the proteins selected earlier
protein_cols <- selected_proteins2  # these names must be found in full_data

# create an empty list to store the regression results
results_list <- list()

# loop over each protein and run a regression
for (protein_name in protein_cols) {
  
  # build the regression formula
  formula <- as.formula(paste(protein_name, "~ prs_bp_std + age_0_0 + sex_plink +",
                              "PGC1 + PGC2 + PGC3 + PGC4 + PGC5 + PGC6 + PGC7 + PGC8 + chip"))
  
  # fit the model
  model <- lm(formula, data = filtered_protein2)
  model_summary <- summary(model)
  
  # extract the coefficient for prs_bp_std
  if ("prs_bp_std" %in% rownames(model_summary$coefficients)) {
    coef <- model_summary$coefficients["prs_bp_std", ]
    
    # store the result
    results_list[[protein_name]] <- data.frame(
      protein = protein_name,
      estimate = coef["Estimate"],
      std_error = coef["Std. Error"],
      t_value = coef["t value"],
      p_value = coef["Pr(>|t|)"]
    )
  }
}

# combine all results
regression_results <- do.call(rbind, results_list)


hist(regression_results$p_value, breaks = 50, main = "P-value Distribution", xlab = "P-value")


#cutoff to 0.0010->0.0020

p_hist <- hist(regression_results$p_value,
               breaks = 50,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.002),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour


#Bonferroni
0.05/24
#Bonferroni-corrected significance level
significant2 <- regression_results[regression_results$p_value < 0.0021 & regression_results$p_value !=0,]

####June 30th-July 6th####
####extended model1####
#add more covariates 
#recode alcohol and smoking
#townsend+n_210001_0_0+whr_0+alcoholstatus+smokestatu+basic regression 

#how should tell the story:
#method ->result 
#introduction: we don't know biology of bipolar, study genetics we might now the blood, covariant to bipolar disorder, biomarker who help 
# why we bother this study and what's the point 

#recode alcohol status and smoke status 

# merge 0 and 1 into "non-current"; 2 is "current"
data$alcoholstatus <- factor(ifelse(data$alcoholstatus == 2, "current", "non_current"),
                                  levels = c("non_current", "current"))

data$smokestatus <- factor(ifelse(data$smokestatus == 2, "current", "non_current"),
                                levels = c("non_current", "current"))


# first take the protein names to select (assuming the column is called protein)
selected_proteins3 <- significant2$protein

# then use them to select columns from the protein data
filtered_protein3 <- protein[, selected_proteins3]
#townsend+n_210001_0_0+whr_0+alcoholstatus+smokestatu+basic regression 

covaritate2<- data %>% select(townsend,n_21001_0_0,whr_0,alcoholstatus,smokestatus,prs_bp_std,age_0_0,sex_plink,PGC1,PGC2,PGC3,PGC4,PGC5,PGC6,PGC7,PGC8,chip)
#combine prs with filtered_protein 
filtered_protein3 <- cbind(filtered_protein3,covaritate2)

# take the column names of the proteins selected earlier
protein_cols <- selected_proteins3  # these names must be found in full_data

# create an empty list to store the regression results
results_list <- list()

# loop over each protein and run a regression
for (protein_name in protein_cols) {
  
  # build the regression formula
  formula <- as.formula(paste(
    protein_name, "~ townsend + n_21001_0_0 + whr_0 + alcoholstatus + smokestatus + prs_bp_std + age_0_0 + sex_plink +",
    "PGC1 + PGC2 + PGC3 + PGC4 + PGC5 + PGC6 + PGC7 + PGC8 + chip"
  ))
  
  # fit the model
  model <- lm(formula, data = filtered_protein3)
  model_summary <- summary(model)
  
  # extract the coefficient for prs_bp_std
  if ("prs_bp_std" %in% rownames(model_summary$coefficients)) {
    coef <- model_summary$coefficients["prs_bp_std", ]
    
    # store the result
    results_list[[protein_name]] <- data.frame(
      protein = protein_name,
      estimate = coef["Estimate"],
      std_error = coef["Std. Error"],
      t_value = coef["t value"],
      p_value = coef["Pr(>|t|)"]
    )
  }
}

# combine all results
regression_results2 <- do.call(rbind, results_list)


hist(regression_results2$p_value, breaks = 50, main = "P-value Distribution", xlab = "P-value")


#cutoff to 0.0010->0.0020

p_hist <- hist(regression_results2$p_value,
               breaks = 50,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.0005),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour


#Bonferroni
0.05/20
#Bonferroni-corrected significance level
significant3 <- regression_results2[regression_results2$p_value < 0.0025 & regression_results2$p_value !=0,]

####July 7th to July 13th####

####extension model 2####
# first take the protein names to select (assuming the column is called protein)
selected_proteins4 <- significant3$protein

# then use them to select columns from the protein data
filtered_protein4 <- protein[, selected_proteins4]
#townsend+n_210001_0_0+whr_0+alcoholstatus+smokestatu+basic regression 

covaritate3<- data %>% select(sbp_0,dbp_0,stroke_0,ht_0,ish_0,lipidlowering_0,htnmed_0,adj_sbp_0,adj_dbp_0,adj_probposs_t2dm_0,townsend,n_21001_0_0,whr_0,alcoholstatus,smokestatus,prs_bp_std,age_0_0,sex_plink,PGC1,PGC2,PGC3,PGC4,PGC5,PGC6,PGC7,PGC8,chip)
#combine prs with filtered_protein 
filtered_protein4 <- cbind(filtered_protein4,covaritate3)

# take the column names of the proteins selected earlier
protein_cols <- selected_proteins4  # these names must be found in full_data

# create an empty list to store the regression results
results_list <- list()

# loop over each protein and run a regression
for (protein_name in protein_cols) {
  
  # build the regression formula
  formula <- as.formula(paste(
    protein_name,"~ sbp_0 + dbp_0 + stroke_0 + ht_0 + ish_0 + lipidlowering_0 + htnmed_0 +
      adj_sbp_0 + adj_dbp_0 + adj_probposs_t2dm_0 + townsend + n_21001_0_0 + whr_0 +
      alcoholstatus + smokestatus + prs_bp_std + age_0_0 + sex_plink + 
      PGC1 + PGC2 + PGC3 + PGC4 + PGC5 + PGC6 + PGC7 + PGC8 + chip"
  ))
  
  # fit the model
  model <- lm(formula, data = filtered_protein4)
  model_summary <- summary(model)
  
  # extract the coefficient for prs_bp_std
  if ("prs_bp_std" %in% rownames(model_summary$coefficients)) {
    coef <- model_summary$coefficients["prs_bp_std", ]
    
    # store the result
    results_list[[protein_name]] <- data.frame(
      protein = protein_name,
      estimate = coef["Estimate"],
      std_error = coef["Std. Error"],
      t_value = coef["t value"],
      p_value = coef["Pr(>|t|)"]
    )
  }
}

# combine all results
regression_results3 <- do.call(rbind, results_list)


hist(regression_results3$p_value, breaks = 50, main = "P-value Distribution", xlab = "P-value")


#cutoff to 0.0010->0.0020

p_hist <- hist(regression_results3$p_value,
               breaks = 50,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.0010),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour


#Bonferroni
0.05/19
#Bonferroni-corrected significance level
significant4 <- regression_results3[regression_results3$p_value < 0.002631579 & regression_results3$p_value !=0,]

####July 14th to July 21st####

####extended model 3####
# first take the protein names to select (assuming the column is called protein)
selected_proteins5 <- significant4$protein

# then use them to select columns from the protein data
filtered_protein5 <- protein[, selected_proteins5]
#townsend+n_210001_0_0+whr_0+alcoholstatus+smokestatu+basic regression 

covaritate4<- data %>% select(any_antidepBL,any_antipsychBL,sbp_0,dbp_0,stroke_0,ht_0,ish_0,lipidlowering_0,htnmed_0,adj_sbp_0,adj_dbp_0,adj_probposs_t2dm_0,townsend,n_21001_0_0,whr_0,alcoholstatus,smokestatus,prs_bp_std,age_0_0,sex_plink,PGC1,PGC2,PGC3,PGC4,PGC5,PGC6,PGC7,PGC8,chip)
#combine prs with filtered_protein 
filtered_protein5 <- cbind(filtered_protein5,covaritate4)

# take the column names of the proteins selected earlier
protein_cols <- selected_proteins5 # these names must be found in full_data

# create an empty list to store the regression results
results_list <- list()

# loop over each protein and run a regression
for (protein_name in protein_cols) {
  
  # build the regression formula
  formula <- as.formula(paste(
    protein_name, "~ any_antidepBL + any_antipsychBL + sbp_0 + dbp_0 + stroke_0 + ht_0 + ish_0 + lipidlowering_0 + htnmed_0 +",
    "adj_sbp_0 + adj_dbp_0 + adj_probposs_t2dm_0 + townsend + n_21001_0_0 + whr_0 +",
    "alcoholstatus + smokestatus + prs_bp_std + age_0_0 + sex_plink +",
    "PGC1 + PGC2 + PGC3 + PGC4 + PGC5 + PGC6 + PGC7 + PGC8 + chip"
  ))
  
  # rest of the code unchanged...
# fit the model
  model <- lm(formula, data = filtered_protein5)
  model_summary <- summary(model)
  
  # extract the coefficient for prs_bp_std
  if ("prs_bp_std" %in% rownames(model_summary$coefficients)) {
    coef <- model_summary$coefficients["prs_bp_std", ]
    
    # store the result
    results_list[[protein_name]] <- data.frame(
      protein = protein_name,
      estimate = coef["Estimate"],
      std_error = coef["Std. Error"],
      t_value = coef["t value"],
      p_value = coef["Pr(>|t|)"]
    )
  }
}

# combine all results
regression_results4 <- do.call(rbind, results_list)


hist(regression_results4$p_value, breaks = 50, main = "P-value Distribution", xlab = "P-value")


#cutoff to 0.0010->0.0020

p_hist <- hist(regression_results4$p_value,
               breaks = 50,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.01),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour


#Bonferroni
0.05/11
#Bonferroni-corrected significance level
significant5 <- regression_results4[regression_results4$p_value < 0.004545455 & regression_results4$p_value !=0,]

####July 22nd- July27thh####
                                                                                                                                                                  
#take the more stable proteins
#have a look of what these proteins do 
#search: protein bipolar disorder 
#omen database

# # load the haven package
# library(haven)
# 
# # read the .dta file
# data2 <- read_dta("<psychiatric history file>.dta")
# 
# # view the data
# View(data2)
# 
# #write csv 
# write.csv(data2,"supplement.csv")
# 
# #remove all negatives
# data2 <-data2 %>% filter(n_eid>0)
# 
# #merge two data together 
# data3<-merge(data,data2,by = "n_eid")
# 
# write.csv(data3,"sensitve.csv")
 
 
 data3 <- read.csv("sensitive.csv")
 
 #remove anyMentalIllness 
 data4 <- data3 %>% filter(anyMentalIllness != 1 | is.na(anyMentalIllness))
 

#townsend+n_210001_0_0+whr_0+alcoholstatus+smokestatu+basic regression 

selected_data <- data4 %>% select(
  any_antidepBL, any_antipsychBL, sbp_0, dbp_0, stroke_0, ht_0, ish_0, lipidlowering_0,
  htnmed_0, adj_sbp_0, adj_dbp_0, adj_probposs_t2dm_0, townsend, n_21001_0_0, whr_0,
  alcoholstatus, smokestatus, prs_bp_std, age_0_0, sex_plink, PGC1, PGC2, PGC3, PGC4,
  PGC5, PGC6, PGC7, PGC8, chip,
  all_of(tolower(significant5$protein))  # PUBLIC COPY: the original typed the protein names here (withheld until publication)
)

#select_protein 
selected_proteins6 <-data4%>% select(all_of(tolower(significant5$protein)))  # PUBLIC COPY: names withheld, as above

# take the column names of the proteins selected earlier
# take the column names of selected_proteins6 as a vector of variable names (so the protein names appear as values, not as data):
protein_cols <- names(selected_proteins6) # these names must be found in full_data

# create an empty list to store the regression results
results_list <- list()

# loop over each protein and run a regression
for (protein_name in protein_cols) {
  
  # build the regression formula
  formula <- as.formula(paste(
    protein_name, "~ any_antidepBL + any_antipsychBL + sbp_0 + dbp_0 + stroke_0 + ht_0 + ish_0 + lipidlowering_0 + htnmed_0 +",
    "adj_sbp_0 + adj_dbp_0 + adj_probposs_t2dm_0 + townsend + n_21001_0_0 + whr_0 +",
    "alcoholstatus + smokestatus + prs_bp_std + age_0_0 + sex_plink +",
    "PGC1 + PGC2 + PGC3 + PGC4 + PGC5 + PGC6 + PGC7 + PGC8 + chip"
  ))
  
  # rest of the code unchanged...
  # fit the model
  model <- lm(formula, data = selected_data)
  model_summary <- summary(model)
  
  # extract the coefficient for prs_bp_std
  if ("prs_bp_std" %in% rownames(model_summary$coefficients)) {
    coef <- model_summary$coefficients["prs_bp_std", ]
    
    # store the result
    results_list[[protein_name]] <- data.frame(
      protein = protein_name,
      estimate = coef["Estimate"],
      std_error = coef["Std. Error"],
      t_value = coef["t value"],
      p_value = coef["Pr(>|t|)"]
    )
  }
}

# combine all results
regression_results5 <- do.call(rbind, results_list)


hist(regression_results5$p_value, breaks = 50, main = "P-value Distribution", xlab = "P-value")


#cutoff to 0.0010->0.0020

p_hist <- hist(regression_results5$p_value,
               breaks = 50,
               main = "P value Distribution",
               xlab = "P value",
               xlim = c(0, 0.005),
               col = "gray",
               ylim = c(0, max(p_hist$counts) * 1.15))  # leave room on the y axis for the labels

# add a label to each bar
text(x = p_hist$mids,               # center of each bar
     y = p_hist$counts,             # height of each bar
     labels = ifelse(p_hist$counts == 0, "", p_hist$counts),  # show non-zero counts only
     pos = 3,                       # show above the bar
     cex = 0.6,                     # font size
     col = "black")                # label colour


#Bonferroni
0.05/11
#Bonferroni-corrected significance level
significant6 <- regression_results5[regression_results5$p_value < 0.004545455 & regression_results4$p_value !=0,]

#####July 28th-Aug4th####


# 1. main analysis results (extended model3)
extended_model3_main <- significant5 %>%
  rename(
    beta_main = estimate,
    p_main = p_value
  )

# 2. sensitivity analysis results
extended_model3_sens <- significant6 %>%
  rename(
    beta_sens = estimate,
    p_sens = p_value
  )

# 3. merge the two result tables, matched by protein
comparison_df <- left_join(extended_model3_main, extended_model3_sens, by = "protein")

# 4. add the change in estimate
comparison_df <- comparison_df %>%
  mutate(
    delta_beta = beta_sens - beta_main,
    delta_p = p_sens - p_main
  )

# 5. view the result
print(comparison_df)

library(ggplot2)

# scatter plot of the beta values
ggplot(comparison_df, aes(x = beta_main, y = beta_sens, label = protein)) +
  geom_point(size = 3, color = "steelblue") +  # points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # y = x line
  geom_text(vjust = -1, size = 3) +  # add protein labels
  labs(
    title = "Comparison of Beta Estimates: Main vs Sensitivity Analysis",
    x = "Beta (Main Analysis)",
    y = "Beta (Sensitivity Analysis)"
  ) +
  theme_minimal()

# rename "Protein" to "protein" in the column names of significant
names(significant)[names(significant) == "Protein"] <- "protein"

# create a list containing all the data frames
df_list <- list(significant, significant2, significant3, significant4, significant5, significant6)

# extract the protein column and sort it
protein_lists <- lapply(df_list, function(df) sort(df$protein))

# pad to the same length (fill with NA)
max_len <- max(sapply(protein_lists, length))
protein_lists_padded <- lapply(protein_lists, function(vec) {
  length(vec) <- max_len  # pad with NA
  return(vec)
})

# combine into a data frame
protein_df <- as.data.frame(do.call(cbind, protein_lists_padded))

# set the column names
colnames(protein_df) <- c("Original", "Basic", "Extended1", "Extended2", "Extended3", "Sensitive")

# view the result
print(protein_df)

# Step 1: main list, from the original model (significant)
main_proteins <- sort(unique(significant$protein))  # sort alphabetically

# Step 2: put all the other data into a list
df_list <- list(Basic = significant2,
                Extended1 = significant3,
                Extended2 = significant4,
                Extended3 = significant5,
                Sensitive = significant6)

# Step 3: create an empty result data frame with the main list as row names
result_df <- data.frame(Original = main_proteins, stringsAsFactors = FALSE)

# Step 4: for each other data frame, check whether each protein is in the main list
for (name in names(df_list)) {
  proteins <- df_list[[name]]$protein
  matched <- ifelse(main_proteins %in% proteins, main_proteins, "")
  result_df[[name]] <- matched
}

# Step 5: view the result
print(result_df)

print(result_df$Original)

write.csv(result_df,"proteinlist.csv")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

BiocManager::install("org.Hs.eg.db")


library(clusterProfiler)
library(org.Hs.eg.db)

# assume the protein names are gene symbols
genes <- tolower(significant6$protein)  # PUBLIC COPY: the original typed the final protein list here (withheld until publication)

genes_upper <- toupper(genes)

# convert to Entrez IDs (required by clusterProfiler)
gene_ids <- bitr(genes_upper,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)
# GO enrichment
ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
head(ego)

# KEGG enrichment
ekegg <- enrichKEGG(gene = gene_ids$ENTREZID, organism = 'hsa')
head(ekegg)

library(enrichplot)
#enrichplot v1.28.4 Learn more at https://yulab-smu.top/contribution-knowledge-mining/

# Please cite:
#   
#   Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for reactome pathway analysis and
# visualization. Molecular BioSystems. 2016, 12(2):477-479

# GO bubble plot
dotplot(ego, showCategory = 15, title = "GO Enrichment") 


# extract the GO enrichment results
go_results <- ego@result

write.csv(go_results,"goenrichment.csv")


#Data structure 
# select numeric columns only
numeric_data <- data[sapply(data, is.numeric)]
total_rows = nrow(data)

# compute mean, sd, max, min and NA percentage of each column
mean_values <- sapply(numeric_data, mean, na.rm = TRUE)
sd_values <- sapply(numeric_data, sd, na.rm = TRUE)
max_values <-sapply(numeric_data,max,na.rm = TRUE)
min_values <-sapply(numeric_data,min,na.rm = TRUE)
na_percentage <- sapply(numeric_data, function(x) sum(is.na(x)) / total_rows * 100)

# combine into a new table
summary_df <- data.frame(
  Variable = names(mean_values),
  Mean = round(mean_values,2),
  SD =  round(sd_values,2),
  Max = round(max_values,2),
  Min = round( min_values,2),
 NAPercent = round(na_percentage,2)
)

print(summary_df)

write.csv(summary_df,"EDA.csv")


# assume spearman is the original table, which contains non-numeric columns
Spearman_rank <- spearman

# round the numeric columns and leave the others unchanged
Spearman_rank[sapply(Spearman_rank, is.numeric)] <- 
  round(Spearman_rank[sapply(Spearman_rank, is.numeric)], 2)

write.csv(Spearman_rank,"Spearman.csv")

# use the protein names as row names
rownames(protein_df) <- protein_df$Original
protein_df$Original <- NULL

# list of proteins to keep
sensitive_proteins <- tolower(significant6$protein)  # PUBLIC COPY: the original typed the final protein list here
                                                     # (withheld until publication)

# keep the rows of the spearman data whose Protein is in sensitive_proteins
Sensitive_Spearman <- spearman %>%
  filter(Protein %in% significant6$protein)  # note the $protein

write.csv(Sensitive_Spearman,"Senstive_spearman.csv")

library(dplyr)

Analysis <- data %>%
  subset(select = c(n_21001_0_0, age_0_0, sex_plink, alcoholstatus, smokestatus,townsend)) %>%
  na.omit()

#BMI

ggplot(Analysis, aes(x = n_21001_0_0)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(x = "BMI", y = "Count") +
  theme_minimal()

#Age

ggplot(Analysis, aes(x = age_0_0)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(x = "Age", y = "Count") +
  theme_minimal()
#Sex

# convert sex to a factor with labels
Analysis <- Analysis %>%
  mutate(sex = factor(sex_plink, levels = c(1, 2), labels = c("Male", "Female")))

# draw a bar chart
ggplot(Analysis, aes(x = sex)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(x = "Sex", y = "Count") +
  theme_minimal()

#Alchohol

ggplot(Analysis, aes(x = alcoholstatus)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(x = "Alcohol Status", y = "Count") +
  theme_minimal()


#Smoke

ggplot(Analysis, aes(x = smokestatus)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(x = "Smoke Status", y = "Count") +
  theme_minimal()

#townsend 
ggplot(Analysis, aes(x = age_0_0)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(x = "Townsend", y = "Count") +
  theme_minimal()






