library(googlesheets)

# Statistical analysis of principal components, raw scores
root_data_directory <- "/media/db242421/db242421_data/ConPagnon_data"
save_results_directory <- "/media/db242421/db242421_data/ConPagnon_reports/resultsPCA/additional_stats_analyses"

# Create directory for saving results
dir.create(save_results_directory, showWarnings = TRUE)

# Load behavioral dataframe.
data_gs_title <- gs_title("Resting State AVCnn: cohort data")
data <- gs_read(ss = data_gs_title, ws="Middle Cerebral Artery resting state cohort")


#data <- as.data.frame(read_excel(file.path(root_data_directory,"regression_data/regression_data.xlsx")))
patients_data <- as.data.frame(data[(data$Groupe == 'P'), ])
rownames(patients_data) <- patients_data$X1
patients_data$X1 <- NULL
patients_data <- as.data.frame(patients_data)

# Raw variable, and principal components results per
# behavioral domain
language_neel <- c("uni_deno", "plu_deno", "uni_rep", "plu_rep", 
                   "empan", "phono", "elision_i", "invers",
                   "ajout", "elision_f", "morpho", "listea", 
                   "listeb", "topo", "voc1", "voc2",
                   "voc1_ebauche", "voc2_ebauche", "abstrait_diff", "abstrait_pos",
                   "lex1", "lex2")

pca_language_neel <- c("pc1_language", "pc2_language", "pc3_language")

wisc_tests <- c("wisc_sim", "wisc_voca",  "wisc_comp", 
                "wisc_cube", "wisc_idc", "wisc_mat", 
                "wisc_memo", "wisc_seq", 
                "wisc_code", "wisc_sym")

pca_wisc <- c("pc1_wisc", "pc2_wisc", "pc3_wisc")


motor <- c("bbt_left_hand", "bbt_right_hand", "nhpt_left", 
           "nhpt_right")

pca_motor <- c("motor_pc1", "motor_pc2", "motor_pc3")

lexical_decoding <- c("alou_m", "alou_e", "alou_c",
                      "alou_cm", "alou_ctl")

pca_lexical_decoding <- c("pc1_lexical_decoding")

executive_functions <- c("rey_copie", "rey_dessin")

pca_executive <- c("executive_pc1")


all_scores <- c(executive_functions, lexical_decoding, motor,
                wisc_tests, language_neel)

all_pca <- c(pca_wisc, pca_motor, pca_language_neel, pca_executive, pca_lexical_decoding)


# Raw variables: link with normalized lesion volume
# Per domain ?
domain <- language_neel
domain_name <- "Language"

# Subsetting the dataframe to the clinical domain of interest
domain_patients_data <-patients_data[, domain, drop = FALSE]
# Drop rows containing missing values
domain_patients_data <- domain_patients_data[complete.cases(domain_patients_data), ]
# Make sure dataframe contain numeric values only
domain_patients_data_ <- sapply(domain_patients_data, as.numeric)
rownames(domain_patients_data_) <- rownames(domain_patients_data)
domain_patients_data <- as.data.frame(domain_patients_data_)
# Add normalized lesion size
variable <- c("lesion_normalized")
domain_patients_data_w_lesion <-merge(domain_patients_data , patients_data[variable ], by = "row.names")
rownames(domain_patients_data_w_lesion) <- domain_patients_data_w_lesion$Row.names
domain_patients_data_w_lesion$Row.names <- NULL

# Clean dataframe, column containing ","
domain_patients_data_w_lesion <- sapply(domain_patients_data_w_lesion[variable], as.numeric)
sapply(domain_patients_data_w_lesion, class)
# Correlation with lesion volume
correlation_w_lesion <- cor(domain_patients_data_w_lesion)
# Correlation without lesion volume
correlation_wo_lesion <- cor(domain_patients_data)

# Color palette
col3 <- colorRampPalette(c("blue", "white", "red"))

pdf(file = file.path(save_results_directory, paste("correlation_matrix_", domain_name, "_without_lesion.pdf", sep = "")))
corrplot(correlation_wo_lesion, method = 'color',
         title = paste('Correlation matrix between behavioral tests for ', 
                       domain_name, sep = " "),
         number.cex = 0.5, col = col3(100), 
         addgrid.col = 'black', tl.col = 'black', tl.cex = 0.75,
         mar=c(0,0,1,0))
dev.off()

# Compute correlation matrix between raw behavioral scores regressing out 
# lesion volume effect
correlation_w_lesion_volume_regressed <- partial.r(data = domain_patients_data_w_lesion, 
                                                   x = domain, 
                                                   y = variable)

pdf(file = file.path(save_results_directory, 
                     paste("correlation_matrix_", domain_name, "_with_lesion_regressed.pdf", sep = "")),
    height = 25,
    width = 20)
corrplot(correlation_w_lesion_volume_regressed, method = 'color',
         title = paste('Correlation matrix between behavioral tests for ', 
                       domain_name, "\n with lesion volume regressed out",sep = " "),
         number.cex = 0.5, col = col3(100), 
         addgrid.col = 'black', tl.col = 'black', tl.cex = 0.75,
         mar=c(0,0,1,0))
dev.off()

# Compute difference with and without volume regression
correlation_difference <- correlation_wo_lesion - correlation_w_lesion_volume_regressed



