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

wd = "/media/db242421/db242421_data/ConPagnon_data/regression_data"
setwd(wd)

sheet <- "ACM_resting_state_cohort"
data <- read_excel(path = 'regression_data.xlsx', sheet = sheet)
patients_data <- data[(data$Groupe == 'P'), ]
rownames(patients_data) <- patients_data$X__1
patients_data$X__1 <- NULL
patients_data <- as.data.frame(patients_data)

# Columns name of behavioral domain of interest
language_tests <- c("uni_deno", "plu_deno", "uni_rep", "plu_rep", 
                    "empan", "phono", "elision_i", "invers",
                    "ajout", "elision_f", "morpho", "listea", 
                    "listeb", "topo", "voc1", "voc2",
                    "voc1_ebauche", "voc2_ebauche", "abstrait_diff", "abstrait_pos",
                    "lex1", "lex2")
wisc_tests <- c("wisc_sim", "wisc_voca",  "wisc_comp", "wisc_irp", 
                "wisc_cube", "wisc_idc", "wisc_mat", "wisc_imt", 
                "wisc_memo", "wisc_seq",  "wisc_arith", "wisc_ivt", 
                "wisc_code", "wisc_sym")

# Subsetting the dataframe to the clinical domain of interest 
domain_patients_data <-patients_data[, language_tests, drop = FALSE]

# Domain of interest
domain <- "Language"


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


# Screeplot to choose number of components
dev.new()
scree(rx = domain_patients_data, 
      main = 'eigenvalues - principal component plot', 
      factors = FALSE)


# Number of components to keep based on Kaiser Rule: eigenvalues superior to 1.
ncomp <- 5

# Perform PCA
res.pca <- PCA(X = domain_patients_data, scale.unit = TRUE, 
               graph = FALSE, ncp = ncomp)

# Show percentage of variance explained
res.pca.eig <- res.pca$eig
barplot(res.pca.eig[, "percentage of variance"])

# Correlation between variable and principal components: loadings
res.pca.loadings <- res.pca$var$coord

dev.new()
corrplot(res.pca.loadings, method = 'number', 
         title = paste('Correlation between principal components and variables for ', 
               domain, sep = " "), 
         tl.cex = 1, col = col3(100), 
         tl.col = 'black', 
         number.cex = 0.75,
         mar=c(0,0,1,0))

# interpretation of loadings
rotation <- "oblimin"
res.pca.withRotation <- pca(r = domain_patients_data, 
                            nfactors = ncomp, 
                            scores = TRUE, 
                            rotate = rotation)
res.pca.withRotation.loading <- res.pca.withRotation$loadings


