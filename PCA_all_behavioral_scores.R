library(factoextra)
library(FactoMineR)
library(readxl)
library(psych)
library(corrplot)
library(plyr)
library(easyGgplot2)
library(ggrepel)
library(pracma)
library(NbClust)
library(randomcoloR)
library(sparsepca)
library(dplyr)
library(reshape2)
library(GPArotation)
library(paran)


wd = "/media/db242421/db242421_data/ConPagnon_data/regression_data"
setwd(wd)

# Domain of interest
domain <- "Executive"

# Save results
save_results_directory <- paste("/media/db242421/db242421_data/ConPagnon_reports/resultsPCA/", domain, sep = "")
dir.create(path = save_results_directory, showWarnings = TRUE)

sheet <- "ACM_resting_state_cohort"
data <- read_excel(path = 'regression_data.xlsx', sheet = sheet)
patients_data <- as.data.frame(data[(data$Groupe == 'P'), ])
rownames(patients_data) <- patients_data$X__1
patients_data$X__1 <- NULL
patients_data <- as.data.frame(patients_data)

# Columns name of behavioral domain of interest
# TODO: Check with Lucie
language_neel <- c("uni_deno", "plu_deno", "uni_rep", "plu_rep", 
                    "empan", "phono", "elision_i", "invers",
                    "ajout", "elision_f", "morpho", "listea", 
                    "listeb", "topo", "voc1", "voc2",
                    "voc1_ebauche", "voc2_ebauche", "abstrait_diff", "abstrait_pos",
                    "lex1", "lex2")

wisc_tests <- c("wisc_sim", "wisc_voca",  "wisc_comp", 
                "wisc_cube", "wisc_idc", "wisc_mat", 
                "wisc_memo", "wisc_seq", 
                "wisc_code", "wisc_sym")

motor <- c("bbt_left_hand", "bbt_right_hand", "nhpt_left", 
           "nhpt_right")

lexical_decoding <- c("alou_tl", "alou_m", "alou_e", "alou_c",
                      "alou_cm", "alou_ctl")

executive_functions <- c("rey_copie", "rey_dessin")


# Subsetting the dataframe to the clinical domain of interest 
domain_patients_data <-patients_data[, executive_functions, drop = FALSE]
# Drop rows containing missing values
domain_patients_data <- domain_patients_data[complete.cases(domain_patients_data), ]
# Make sure dataframe contain numeric values only
domain_patients_data_ <- sapply(domain_patients_data, as.numeric)
rownames(domain_patients_data_) <- rownames(domain_patients_data)
domain_patients_data <- domain_patients_data_

# Correlation between scores
cormat <- cor(domain_patients_data)
col3 <- colorRampPalette(c("blue", "white", "red"))

#dev.new()
pdf(file = file.path(save_results_directory, paste("correlation_matrix_", domain, ".pdf", sep = "")))
corrplot(cormat, method = 'color',
         title = paste('Correlation matrix between behavioral tests for ', 
                       domain, sep = " "),
         number.cex = 0.5, col = col3(100), 
         addgrid.col = 'black', tl.col = 'black', tl.cex = 0.75,
         mar=c(0,0,1,0))
dev.off()

# Choose the number of components: scree plot method,
# percentage of variance explained, Parralel Analysis, 
# Broken Stick model.

# Screeplot to choose number of components
#dev.new()
pdf(file = file.path(save_results_directory, paste("scree_plot_", domain, ".pdf", sep = "")))
scree(rx = domain_patients_data, 
      main = 'eigenvalues - principal component plot', 
      factors = FALSE)
dev.off()

# Number of components to keep based on Kaiser Rule: eigenvalues superior to 1, at least
# 70 % of variance explained
ncomp <- 1

# Perform PCA
res.pca <- PCA(X = domain_patients_data, scale.unit = TRUE,
              graph = FALSE, ncp = ncomp)

res.pca.scores <- res.pca$ind$coord
res.pca.scores.std <- scale(res.pca.scores, scale = TRUE, center = TRUE)

# Save eigenvalues (no rotation) and percentage of 
# variances explained 
res.pca.eig <- res.pca$eig
write.file.csv(x = res.pca.eig, 
               f = file.path(save_results_directory, "no_rotation_eig.csv"),
               row.names = TRUE)

# Show percentage of variance explained
#dev.new()
pdf(file = file.path(save_results_directory, paste("barplot_eigenvalues_", domain, ".pdf", sep = "")))
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80),
         xlab = "Principal Components", ylab = "Percentage of variance explained",
         main = "Variance explained (%) by each principal components",
         barfill = "orange", barcolor = "black")
dev.off()


# Parallel Analysis
#dev.new()
pdf(file = file.path(save_results_directory, paste("paralel_analysis_", domain, ".pdf", sep = "")))
parralel_analysis <- paran(domain_patients_data, iterations = 5000, centile = 0, quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)
dev.off()

# Broken Stick model
# Quelles sont les valeurs propres plus grandes que la moyenne?
ev <- res.pca$eig[,"eigenvalue"]
ev[ev > mean(ev)]

# Modele du baton brise (broken stick model)

# Distribute in a random fashion the amount of variance 
# along all the components
n = length(ev)
bsm = data.frame(j=seq(1:n), p=0)
bsm$p[1] = 1/n
for (i in 2:n) {
  bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
}
bsm$p = 100*bsm$p/n

# Dessiner les valeurs propres et le % de variance de chaque axe
#dev.new()
pdf(file = file.path(save_results_directory, paste("broken_stick_model_", domain, ".pdf", sep = "")))
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")	# mean of eigenvalues
legend("topright", "Mean of eigenvalues", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev),bsm$p[n:1])), beside=T, 
        main="% variance", col=c("bisque",2), las=2)
legend("topright", c("% of variance explained", "Broken stick model"), 
       pch=15, col=c("bisque",2), bty="n")
dev.off()

# Correlation between variable and principal components: loadings
#res.pca.loadings <- res.pca$var$coord

#dev.new()
#corrplot(res.pca.loadings, method = 'number', 
#         title = paste('Correlation between principal components and variables for ', 
#               domain, sep = " "), 
#         tl.cex = 1, col = col3(100), 
#         tl.col = 'black', 
#         number.cex = 0.75,
#         mar=c(0,0,1,0))

# interpretation of loadings: Try a rotation of principal components if desired
# in the psych packages, computed scores are standardized
rotation <- "none"

# Perform PCA with psych package
# Note: The scores are standardized
res.pca.withRotation <- principal(r = domain_patients_data, 
                            nfactors = ncomp, 
                            scores = TRUE, 
                            rotate = rotation)
if (rotation == "none") {
  rotation_colnames <- c(paste0("PC", 1:ncomp))
  
}
if (rotation == "varimax") {
  rotation_colnames <- c(paste0("RC", 1:ncomp))
  
}
if (rotation == "oblimin") {
  rotation_colnames <- c(paste0("TC", 1:ncomp))  
}

# Loadings matrix
res.pca.withRotation.loading <- res.pca.withRotation$loadings[, rotation_colnames]

# plot the loading matrix

# Reshape the data before plotting loading 
# as heatmap with ggplot
melted_loading <- melt(res.pca.withRotation.loading)

# Plot the loading as heatmap matrix with ggplot
#dev.new()
pdf(file = file.path(save_results_directory, paste("loadings_heatmap_", domain, ".pdf", sep = "")))
ggplot(melted_loading, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 3))) +
  scale_fill_distiller(palette ="RdBu", direction = -1) +
  theme_gray(base_size = 10) +
  coord_flip() +
  ylab("Principal Components") + 
  xlab("Variables") + 
  ggtitle(label = paste("Loadings of variables on principal components for", domain, sep = " "),
          subtitle = paste("Rotation:", rotation)) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(face = "bold", colour = "grey50"),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 8, face = "bold", colour = "black"),
        plot.title = element_text(size = 12, colour = "grey50", face = "bold")) 
dev.off()

# Asses the significance of the Loadings
# Compute the Peres-Neto methods.
# See: "Giving meaningul interpretation to ordination axes: Assesing
# loading significance in principal components analysis", 2003,
# Ecology.

# Squared the loadings matrix
n_var <- dim(x = res.pca.withRotation.loading)[1]
squared_loading <- res.pca.withRotation.loading**2

# Compute the loadings distribution under the broken stick 
# model
b_k <-  data.frame(j=seq(1:n_var), p=0)
b_k$p[1] = (1/n_var)

for (i in 2:n_var) {
  b_k$p[i] = b_k$p[i-1] + (1/(n_var + 1 - i))
}
b_k$p = b_k$p/n_var
# Reverse element: the biggest proportion
# for PC1, etc...
b_k$p = rev(b_k$p)

# Rank squared loading and assess their significance according 
# to the borken stick model in descending order
squared_loading_ranking <- t(apply(-squared_loading,1,rank))
significant_loading_boolean <- cbind(squared_loading_ranking)

for(var in rownames(squared_loading_ranking)){
  # compare squared loading to expected
  # proportion of variance under broken stick model
  cmp <- squared_loading[var, ] > b_k$p[squared_loading_ranking[var, ]]
  # Fill the significant matrix of loading for each variables
  significant_loading_boolean[var, ] <- cmp
  
}

# Finaly, the raw loading under the broken
# stick model
significant_loading <- res.pca.withRotation.loading * significant_loading_boolean

# Plot the heatmap of significant loadings
# Reshape the data before plotting loading 
# as heatmap with ggplot
melted_significant_loading <- melt(significant_loading)

# Plot the loading as heatmap matrix with ggplot
#dev.new()
pdf(file = file.path(save_results_directory, paste("significant_loadings_", domain, ".pdf", sep = "")))
ggplot(melted_significant_loading, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 3))) +
  scale_fill_distiller(palette ="RdBu", direction = -1) +
  theme_gray(base_size = 10) +
  coord_flip() +
  ylab("Principal Components") + 
  xlab("Variables") + 
  ggtitle(label = paste("Significant loadings under the broken stick models for ", domain, sep = " "),
          subtitle = paste("Rotation:", rotation)) +
  theme(axis.ticks = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(face = "bold", colour = "grey50"),
        axis.text.x = element_text(size = 10, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 8, face = "bold", colour = "black"),
        plot.title = element_text(size = 12, colour = "grey50", face = "bold")) 
dev.off()


# Write significant loadings for easier interpretation
write.file.csv(x = significant_loading_boolean,
               f = file.path(save_results_directory, 
                             paste(rotation, "_significant_loadings_interpretation.csv", sep = "")),
               row.names = TRUE)


# Get Correlation between scores: depending
# rotation used, scores can be correlated if the
# rotatation was not orthogonal.

# If rotation is none save the raw scores
# Else: save standardized output of psych packages
scores.correlation <- res.pca.withRotation$r.scores
write.file.csv(x = scores.correlation, f = file.path(save_results_directory, 
                                                     paste(rotation, "_correlations_scores.csv", sep = "")),
               row.names = TRUE)

# Save the scores: the projection of individuals
# in the principal components space (scores are standardized)
if (rotation == "none") {
  individuals.coord <- res.pca$ind$coord
  
}else{
  individuals.coord <- res.pca.withRotation$scores
}
write.file.csv(x = individuals.coord, f = file.path(save_results_directory, 
                                                     paste(rotation, "_individuals_coordinates.csv", sep = "")),
               row.names = TRUE)

# Save the loadings: the correlation
# between raw variables and principal components.
write.file.csv(x = res.pca.withRotation.loading, f = file.path(save_results_directory, 
                                                    paste(rotation, "_loadings.csv", sep = "")),
               row.names = TRUE)


# Plot the results in the desired plane
# of chosen principal components

# Create a dataframe with the principal components
# clinical variables for plotting purposes.
clinical_variables <- c("Sexe", "Lesion", "langage_clinique", "cerebral_palsy", "Parole",
                        "Lexique_comp", "Lexique_exp") 
plotting_dataframe <- merge(individuals.coord, patients_data[, clinical_variables], by = 0, all = TRUE)
# Drop rows containing missing values
plotting_dataframe <- as.data.frame(plotting_dataframe[complete.cases(plotting_dataframe), ])
# Make the newly created Row.names the rownames of the dataframe
rownames(plotting_dataframe) <- plotting_dataframe$Row.names
plotting_dataframe$Row.names <- NULL

# Build composite factor variable for plotting purpose


plotting_dataframe$langage_clinique <- factor(plotting_dataframe$langage_clinique)
plotting_dataframe$Parole <- factor(plotting_dataframe$Parole)
plotting_dataframe$speech_language_profile <- factor(as.numeric(with(plotting_dataframe, 
                                                                     interaction(langage_clinique, 
                                                                                 Parole)))-1)

plotting_dataframe$lexical_comprehension_speech <- factor(as.numeric(with(plotting_dataframe, 
                                                                          interaction(Lexique_comp, 
                                                                                      Parole)))-1)


plotting_dataframe$lexical_expression_speech <- factor(as.numeric(with(plotting_dataframe, 
                                                                          interaction(Lexique_exp, 
                                                                                      Parole)))-1)

plotting_dataframe$lexical_comprehension_expression_speech <- factor(as.numeric(with(plotting_dataframe, 
                                                                                     interaction(Lexique_exp, 
                                                                                                 Parole,
                                                                                                 Lexique_comp)))-1)
# Speech and Language
# groups <- mapvalues(plotting_dataframe$speech_language_profile, from = c("0","2","3"), 
#                     to = c("Impaired Language and Impaired Speech",
#                            "Impaired language and Non Impaired Speech",
#                            "Non Impaired language and Non Impaired Speech"))


# Speech and lexical comprehension
# groups <- mapvalues(plotting_dataframe$speech_language_profile, from = c("0","1","2","3"),
#                     to = c("Impaired lexical comprehension and Impaired Speech",
#                            "Non Impaired lexical comprehension and Impaired Speech",
#                            "Impaired lexical comprehension and Non Impaired Speech",
#                            "Non Impaired lexical comprehension and Non Impaired Speech"))


# Speech and lexical expression
# groups <- mapvalues(plotting_dataframe$lexical_expression_speech, from = c("0","1","2","3"),
#                     to = c("Impaired lexical expression and Impaired Speech",
#                            "Non Impaired lexical expression and Impaired Speech",
#                            "Impaired lexical expression and Non Impaired Speech",
#                            "Non Impaired lexical expression and Non Impaired Speech"))


# Speech and lexical expression and comprehension
# groups <- mapvalues(plotting_dataframe$lexical_comprehension_expression_speech, from = c("0","1","2","3", "4", "6", "7"),
#                     to = c("Impaired lexical exp/comp and Impaired Speech",
#                            "Impaired lexical comp and Speech and Non Impaired lexical exp",
#                            "Non Impaired Speech and Impaired lexical exp/comp",
#                            "Non Impaired Speech and lexical exp and Impaired lexical comp",
#                            "Impaired Speech and lexical exp and Non Impaired lexical comp",
#                            "Non Impairesd Speech and lexical comp and Impaired lexical exp",
#                            "Non Impaired Speech, lexical comp/exp"
#                            ))


groups <- mapvalues(plotting_dataframe$langage_clinique, from = c("A","N"),
                    to = c("Impaired Language",
                           "Non Impaired Language"))


# groups <- plotting_dataframe$Lesion

groups <- as.factor(groups)

# Figures parameters
dim_on_x <- 1
dim_on_y <- 2
points_labels <- rownames(plotting_dataframe)
size_of_points <- 3
size_of_points_labels <- 4
legend_title <- "Legend: "
legend_labels_size <- 10
points_labels_color <- groups

x_label <- "PC2"
y_label <- "PC3"
figure_title <- paste(domain, ": Projection of subjects in the PC2 and PC3 plan", sep = "")

figure_width <- 20
figure_heigth <- 10

# Draw the plot
#dev.new()
pdf(file = file.path(save_results_directory, 
                     paste(domain,"_",x_label,"_",y_label,"language.pdf")),
    width = figure_width,
    height = figure_heigth)
print(ggplot(data = plotting_dataframe) + 
  geom_point(mapping = aes(x = plotting_dataframe[, colnames(plotting_dataframe)[dim_on_x]], 
                           y = plotting_dataframe[, colnames(plotting_dataframe)[dim_on_y]],
                           color = groups), 
             size =size_of_points) +
  geom_vline(xintercept = 0, linetype ='dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_text_repel(mapping = aes(x = plotting_dataframe[, colnames(plotting_dataframe)[dim_on_x]], 
                                y = plotting_dataframe[, colnames(plotting_dataframe)[dim_on_y]],
                                label = points_labels,
                                color = points_labels_color),
                  size = size_of_points_labels,
                  fontface = 'bold',
                  box.padding = 0.5
  ) +
  xlab(x_label) + 
  ylab(y_label)+
  ggtitle(figure_title, subtitle = paste(' ')) + 
  theme_classic() + 
  theme(plot.title = element_text(size=15, face="bold"), 
        legend.position="bottom", legend.text = element_text(size=legend_labels_size),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size=5), title = legend_title)))
dev.off()

# Save plotting dataframe
write.file.csv(x = plotting_dataframe,
               file = file.path(save_results_directory, paste(domain,"_pca_results_dataframe.csv", sep = "")),
               row.names = TRUE)


