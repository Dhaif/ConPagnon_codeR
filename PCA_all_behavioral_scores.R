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

wd = "/media/db242421/db242421_data/ConPagnon_data/regression_data"
setwd(wd)

sheet <- "ACM_resting_state_cohort"
data <- read_excel(path = 'regression_data.xlsx', sheet = sheet)
patients_data <- data[(data$Groupe == 'P'), ]
rownames(patients_data) <- patients_data$X__1
patients_data$X__1 <- NULL


# Subset the dataframe
data_subset <- data[(data$Groupe == 'P'),  48:61]

# Patients names for plotting purpose
patients_names <- data[(data$Groupe == 'P'),]$X__1
rownames(data_subset) <- patients_names

# Correlation between scores
cormat <- cor(data_subset)
col3 <- colorRampPalette(c("blue", "white", "red"))
dev.new()
corrplot(cormat, method = 'color',title = 'Correlation matrix between behavioral tests.',
         number.cex = 0.5, col = col3(100), 
         addgrid.col = 'black', tl.col = 'black', tl.cex = 0.75)


# Screeplot to choose number of components
dev.new()
scree(rx = data_subset, 
      main = 'eigenvalues - principal component plot', 
      factors = FALSE)


# Number of components to keep based on Kaiser Rule: eigenvalues superior to 1.
ncomp <- 3

# Perform PCA
res.pca <- PCA(X = data_subset, scale.unit = TRUE, 
               graph = FALSE, ncp = ncomp)

# Correlation between variable and principal components: loadings
res.pca.loadings <- res.pca$var$coord

dev.new()
corrplot(res.pca.loadings, method = 'number', 
         title = 'Raw loadings', 
         tl.cex = 1, col = col3(100), 
         tl.col = 'black', number.cex = 0.75)




