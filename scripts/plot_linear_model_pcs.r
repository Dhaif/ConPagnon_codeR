# Plotting interresting results from the linear analysis for the NAIS cohort.

require(ggplot2)
require(ggrepel)
require(readxl)

root_directory <- "/media/dhaif/Samsung_T5/Work/Neurospin/AVCnn/AVCnn_Dhaif/ConPagnon_data/patients_behavior_ACM"
score_folder <- "lexExp_zscore_081019"
score_name <- "lexExp_zscore"


# Load the dataframe containing behavioral attribute to each subject
cohort_behavioral_dataframe <- as.data.frame(read_excel(
  path = "/media/dhaif/Samsung_T5/Work/Neurospin/AVCnn/AVCnn_Dhaif/ConPagnon_data/regression_data/regression_data_2.xlsx"))
rownames(cohort_behavioral_dataframe) <- cohort_behavioral_dataframe$subjects
cohort_behavioral_dataframe$subjects <- NULL
# fetch the patient only dataframe
patients_data <- cohort_behavioral_dataframe[cohort_behavioral_dataframe$Groupe == "P", ]


model <- "mean_contra"

save_in <- file.path(root_directory, score_folder)

network_name <- c( "Primary_Visual")
kind <- "tangent"


# Graph parameter
group_by <- patients_data$Syntaxe_exp
legend_title <- "Syntax (expression)"
legend_filename <- "expression_profile"
 # Lesion side color
group_by <- patients_data$Lesion
legend_title <- "Lesion side"
legend_filename <- "lesion_side"


# Plot the results

#Network model: Fetch the design matrix in the network folder corresponding to 
# the chosen model
for (network in network_name) {
  model_design_matrix <- file.path(root_directory, score_folder, "regression_analysis", kind,
                                   network, paste(model, "_design_matrix.csv", sep=""))
  design_matrix_data <- read.csv(model_design_matrix)
  
  # Make subjects column the new index
  rownames(design_matrix_data) <- design_matrix_data$subjects
  design_matrix_data$subjects <- NULL
  
  # Change row order of the design to match the row order
  # of the patients dataframe
  design_matrix_data <- design_matrix_data[rownames(patients_data), ]
  
  # plot the results
  #dev.new()
 # pdf(file = file.path(save_in, paste(score_name,"_", network, "_", model,"_",legend_filename,".pdf", sep="")))
  ggplot(data = design_matrix_data) + 
    geom_point(aes(x = design_matrix_data[, paste("intra_", network, "_connectivity", sep="")],
                   y = design_matrix_data[, score_name],
                   color = factor(group_by))) +
    geom_smooth(method = "lm", aes(x = design_matrix_data[, paste("intra_", network, "_connectivity", sep="")],
                                  y = design_matrix_data[, score_name]), colour = "black") +
    labs(color = legend_title)  +
    xlab(paste(model, "_connectivity", sep="")) +
    ylab(score_name) + 
    ggtitle(label = paste("Relationship between ", score_name, " and", " intra_", network, "_connectivity", sep="")) +
    theme_classic() +
    theme(legend.position = "bottom")
  ggsave(filename = file.path(save_in, paste(score_name,"_", network, "_", model,"_",legend_filename,".jpg", sep="")))
}


# Plot global model
global_model_design_matrix <- file.path(root_directory, score_folder, "regression_analysis", kind,
                                 paste(model, "_design_matrix.csv", sep=""))
global_design_matrix_data <- read.csv(global_model_design_matrix)

# Make subjects column the new index
rownames(global_design_matrix_data) <- global_design_matrix_data$subjects
global_design_matrix_data$subjects <- NULL

# Change row order of the design to match the row order
# of the patients dataframe
global_design_matrix_data <- global_design_matrix_data[rownames(patients_data), ]

ggplot(data = global_design_matrix_data) + 
  geom_point(aes(x = global_design_matrix_data[, model],
                 y = global_design_matrix_data[, score_name],
                 color = factor(group_by))) +
  geom_smooth(method = "lm", aes(x = global_design_matrix_data[, model],
                                 y = global_design_matrix_data[, score_name]), colour = "black") +
  labs(color = legend_title)  +
  xlab(paste(model, "_connectivity", sep="")) +
  ylab(score_name) + 
  ggtitle(label = paste("Relationship between ", score_name, " and ", model, "_connectivity", sep="")) +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(filename = file.path(save_in, paste(score_name, "_", model,"_",legend_filename,".jpg", sep="")))


