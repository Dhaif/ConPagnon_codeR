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
domain <- "WISC"

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

wisc_tests <- c("wisc_sim", "wisc_voca",  "wisc_comp", "wisc_irp", 
                "wisc_cube", "wisc_idc", "wisc_mat", "wisc_imt", 
                "wisc_memo", "wisc_seq",  "wisc_arith", "wisc_ivt", 
                "wisc_code", "wisc_sym")

motor <- c("bbt_left_hand", "bbt_right_hand", "t9c_drtt_g", 
           "t9c_drtt_d", "EHI")

lexical_decoding <- c("alou_tl", "alou_m", "alou_e", "alou_c",
                      "alou_cm", "alou_ctl")

excutive_functions <- c("rey_copie", "rey_dessin")


# Subsetting the dataframe to the clinical domain of interest 
domain_patients_data <-patients_data[, wisc_tests, drop = FALSE]
# Drop rows containing missing values
domain_patients_data <- domain_patients_data[complete.cases(domain_patients_data), ]
# Make sure dataframe contain numeric values only
domain_patients_data_ <- sapply(domain_patients_data, as.numeric)
rownames(domain_patients_data_) <- rownames(domain_patients_data)
domain_patients_data <- domain_patients_data_



# Screeplot to choose number of components
dev.new()
scree(rx = domain_patients_data, 
      main = 'eigenvalues - principal component plot', 
      factors = FALSE)

# Correlation between scores
cormat <- cor(domain_patients_data)
col3 <- colorRampPalette(c("blue", "white", "red"))
dev.new()
corrplot(cormat, method = 'color',
         title = paste('Correlation matrix between behavioral tests for ', 
                       domain, sep = " "),
         number.cex = 0.5, col = col3(100), 
         addgrid.col = 'black', tl.col = 'black', tl.cex = 0.75,
         mar=c(0,0,1,0))

# Number of components to keep based on Kaiser Rule: eigenvalues superior to 1, at least
# 70 % of variance explained
ncomp <- 3

# Perform PCA
res.pca <- PCA(X = domain_patients_data, scale.unit = TRUE,
              graph = FALSE, ncp = ncomp)

res.pca.scores <- res.pca$ind$coord
res.pca.scores.std <- scale(res.pca.scores, scale = TRUE, center = TRUE)

# Show percentage of variance explained
dev.new()
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 80),
         xlab = "Principal Components", ylab = "Percentage of variance explained",
         main = "Variance explained (%) by each principal components",
         barfill = "orange", barcolor = "black")

# Parallel Analysis
paran(domain_patients_data, iterations = 5000, centile = 0, quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)



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
dev.new()
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
  
# Get Correlation between scores: depending
# rotation used, scores can be correlated if the
# rotatation was not orthogonal.
scores.correlation <- res.pca.withRotation$r.scores
write.file.csv(x = scores.correlation, f = file.path(save_results_directory, 
                                                     paste(rotation, "_correlations_scores.csv", sep = "")),
               row.names = TRUE)

# Save the scores: the projection of individuals
# in the principal components space (scores are standardized)
individuals.coord <- res.pca.withRotation$scores
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

clinical_variables <- c("Sexe", "Lesion", "langage_clinique", "cerebral_palsy", "Parole") 
plotting_dataframe <- merge(individuals.coord, patients_data[, clinical_variables], by = 0, all = TRUE)
# Drop rows containing missing values
plotting_dataframe <- as.data.frame(plotting_dataframe[complete.cases(plotting_dataframe), ])
# Make the newly created Row.names the rownames of the dataframe
rownames(plotting_dataframe) <- plotting_dataframe$Row.names
plotting_dataframe$Row.names <- NULL





# Draw the plot
dev.new()
ggplot(data = plotting_dataframe) + 
  geom_point(mapping = aes(x = PC1, 
                           y = PC2,
                           color = Parole), size =3) +
  geom_vline(xintercept = 0, linetype ='dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_text_repel(mapping = aes(x = PC1, 
                                y = PC2,
                                label = rownames(plotting_dataframe),
                                color = Parole),
                  size = 4,
                  fontface = 'bold',
                  box.padding = 0.5
  ) +
  xlab('PC1') + 
  ylab('PC2')+
  ggtitle("Ind. Coord", subtitle = paste(' ')) + 
  theme_classic() + 
  theme(plot.title = element_text(size=15), legend.position="bottom", legend.text = element_text(size=10)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))




