library(factoextra)
library(FactoMineR)
library(readxl)
library(psych)
library(corrplot)
library(plyr)
library(easyGgplot2)
library(ggrepel)
library(fpc)
source("cluster_on_pca.R")
source("plot_prcomp_cluster.R")
library(pracma)
library(NbClust)
library(randomcoloR)

##### Preparation of DATA #####

# This routine aims to perfom a principal components analysis on the raw scores of language test (NEEL).
# The PCA is performed using the FactoMineR packages, and for the sake of interpretation I perfom a varimax rotation
# on raw loadings.
wd = "/media/db242421/db242421_data/vlsmpca/pca"
setwd(wd)
# Close all open graphics windows
graphics.off()
#### DATA LECTURE AND PREPARATION

rawdatafile = '/media/db242421/db242421_data/Document_utile/SCORES_NEEL.xlsx'
langdataxls = as.data.frame(read_excel(rawdatafile, sheet = 'Scores_brut'))

# Drop some data giving row index, comment if you don't
# DropIndexList <- c(22)
# langdataxls <- langdataxls[-DropIndexList, ]

# Drop subjects with non complete sub test : select only the column with the sub-test
langdataxls <- langdataxls[!rowSums(is.na(langdataxls[c(15:36)])), ]

# Fetch numerical data only 
nums <- sapply(langdataxls, is.numeric)
langdata = langdataxls[, nums]
nsub = nrow(langdata)
nvar = ncol(langdata)

# Drops all subtest referring to metaphonology, comment if you don't
drops <- c('elision_i','invers','ajout','elision_f','empan','phono')
langdata <- langdata[, !(names(langdata) %in% drops )  ]

# Replace missing value with the corresponding column mean
newLangdata <- cbind(langdata)
for(i in 1:ncol(langdata))
{
  newLangdata[is.na(langdata[,i]), i] <- mean(langdata[,i], na.rm = TRUE)
}

# Copie of newLangdata, containing only numerical values.
newLangdata_onlynum <- cbind(newLangdata)

##### Correlation matrices and SCREE plot

# Correlational structure of the data and variables
cormat <- cor(newLangdata_onlynum)
col3 <- colorRampPalette(c("blue", "white", "red"))
dev.new()
#tiff(filename = './Illustration/Correlation_matrix_no_phono.tiff', width = 20, height = 10, units = 'in', res = 360)
#tiff(filename = './Illustration/Correlation_matrix_no_phono_lowres.tiff', width = 1920, height = 1080, units = 'px')
corrplot(cormat, method = 'color',title = 'Correlational structure of NEEL raw scores',number.cex = 0.5, col = col3(100), 
         addgrid.col = 'black', tl.col = 'black', tl.cex = 2)
#dev.off()

# Screeplot : Use of the Kaiser rules, keeping all components with eigenvalue > 1.0
dev.new()
#tiff(filename = './Illustration/Scree_plot.tiff', width = 20, height = 10, units = 'in', res = 360)
scree(rx = newLangdata_onlynum, main = 'eigenvalues - principal component plot', factors = FALSE)
#dev.off()
# Number of components to keep
ncomp <- 3

# With prcomp
# pca.prcomp <- prcomp(x = newLangdata, scale. = TRUE, center = TRUE )
# rawLoadings     <- pca.prcomp$rotation[,1:ncomp] %*% diag(pca.prcomp$sdev, ncomp, ncomp)
# rotatedLoadings <- varimax(rawLoadings)$loadings

##### PERFORM PCA ANALYSIS ####

# With principal
pca.psych <- principal(r = newLangdata_onlynum, nfactors = ncomp, rotate = 'varimax')

# Add a subjects identifier for illustration purposes
#rownames(newLangdata) <- langdataxls$monogramme
rownames(newLangdata) <- langdataxls$monogramme
# Performs PCA with FactoMineR 

# OPERATION ON FACTOR FOR INTERPRETATION PURPOSE
# SOUS COMPOSANTE: syntaxe expression, syntaxe compréhension, lexique expression, lexique compréhension , Parole, Profil Langage
#newLangdata <- cbind(Lexique_C=langdataxls$Lexique_C, newLangdata )
#newLangdata <- cbind(Lexique_E=langdataxls$Lexique_E, newLangdata )
#newLangdata <- cbind(Syntaxe_C=langdataxls$Syntaxe_C, newLangdata )
#newLangdata <- cbind(Syntaxe_E=langdataxls$Syntaxe_E, newLangdata )
#newLangdata <- cbind(Parole=langdataxls$Parole, newLangdata)
#newLangdata <- cbind(Profil_Langagier=langdataxls$Profil_langagier, newLangdata)
#newLangdata <- cbind(Territoire_lesionel=langdataxls$Territoire_simple, newLangdata)
#newLangdata <- cbind(Phonologie=langdataxls$Phonologie, newLangdata)
#newLangdata <- cbind(Lesion=langdataxls$Lésion, newLangdata)


# Composante Parole + Langage : good separation.
#newLangdata$Profil_Langagier_Parole <- factor(as.numeric(with(newLangdata, interaction(Profil_Langagier, Parole)))-1)
#groups <- mapvalues(newLangdata$Profil_Langagier_Parole, from = c("0","1","2","3"), to = c("Impaired Language and Impaired Speech", 
#                                                                                           "Non Impaired language and Impaired Speech", 
#                                                                                           "Impaired language and Non Impaired Speech", 
#                                                                                           "Non Impaired language and Non Impaired Speech"))
#groups <- as.factor(groups)

# Lesion
#groups <- mapvalues(newLangdata$Lesion, from = c("D", "G", "Bilatérale"), to = c("Right lesion", "Left Lesion", "Bilateral Lesion"))
#groups <- as.factor(groups)

# Parole
# groups <- mapvalues(newLangdata$Parole, from = c("A","T"), to = c("Parole Atypique", "Parole Typique"))
# groups <- as.factor(groups)

# Syntax + Lexicon comphrension
# newLangdata$Comphrension <- factor(as.numeric(with(newLangdata, interaction(Lexique_C, Syntaxe_C)))-1)
# groups <- as.factor(mapvalues(newLangdata$Comphrension, from = c("0","1","2","3"), to = c("Lexique A et Syntaxe A", "Lexique A et Syntaxe T", "Lexique T et Syntaxe A","Lexique T et Syntaxe T")))

# Phonologie + Parole
# newLangdata$Phonologie_Parole <-  factor(as.numeric(with(newLangdata, interaction(Phonologie, Parole)))-1)
# groups <- as.factor(mapvalues(newLangdata$Phonologie_Parole, from = c("0","1","2","3"), to = c("Phonologie A et Parole A", "Phonologie T et Parole A", "Phonologie A et Parole T","Phonologie T et Parole T")))

# Phonologie + Parole + Production lexicale
# newLangdata$Phonologie_Parole_Lexique <-  factor(as.numeric(with(newLangdata, interaction(Phonologie, Parole, Lexique_E)))-1)
# groups <- as.factor(newLangdata$Phonologie_Parole_Lexique)
 
# Parole + Lexique_C
# newLangdata$Parole_Lexique_C <-  factor(as.numeric(with(newLangdata, interaction(Parole, Lexique_C)))-1)
# groups <- as.factor(mapvalues(newLangdata$Parole_Lexique_C, from = c("0","1","2","3"), to = c("Parole Atypique Lexique Atypique","Parole Typique Lexique Atypique","Parole Atypique Lexique Typique", "Parole Typique Lexique Typique")))

# Lesion area label
#groups <- as.factor(newLangdata$Territoire_lesionel)

# Langage + Lexique_C
# newLangdata$Langage_Lexique_C <-  factor(as.numeric(with(newLangdata, interaction(Profil_Langagier, Lexique_C)))-1)
# groups <- as.factor(mapvalues(newLangdata$Langage_Lexique_C, from = c("0","2","3"), to = c("Langage Atypique Atypique Lexique Atypique","Langage Atypique Lexique Typique","Langage Typique Lexique Typique")))

# Parole + Lexique_C + Syntaxe_C
# newLangdata$Parole_Lexique_C_Syntaxe_C <-  factor(as.numeric(with(newLangdata, interaction(Parole, Lexique_C, Syntaxe_C)))-1)
# groups <- as.factor(mapvalues(newLangdata$Parole_Lexique_C, from = c("1","2","3","4","5","6","7","0"), 
#                               to = c("Parole Typique Reception lexicale Atypique Syntaxique Atypique",
#                                      "Parole ATypique Reception lexicale Typique Syntaxique Atypique",
#                                      "Parole Typique Reception lexicale Typique Syntaxique ATypique", 
#                                      "Parole ATypique Reception lexicale Atypique Syntaxique Typique",
#                                      "Parole Typique Reception lexicale Atypique Syntaxique Typique",
#                                      "Parole Atypique Reception lexicale Typique Syntaxique Typique",
#                                      "Parole Typique Reception lexicale Typique Syntaxique Typique",
#                                      "Parole ATypique Reception lexicale ATypique Syntaxique ATypique")))
# 

# Parole + Lexique_E
# newLangdata$Parole_Lexique_E <-  factor(as.numeric(with(newLangdata, interaction(Parole, Lexique_E)))-1)
# groups <- as.factor(mapvalues(newLangdata$Parole_Lexique_E, from = c("0","1","2","3"), to = c("Parole Atypique Lexique Atypique","Parole Typique Lexique Atypique","Parole Atypique Lexique Typique", "Parole Typique Lexique Typique")))


# Drops some subjects for the PCA
remove_subjects <-  'no'
remove_subjects_for_PCA <- c(62)
if(remove_subjects == 'yes') {
dropSubjects <- remove_subjects_for_PCA  
groups <- groups[-dropSubjects]
} else {
  dropSubjects <- NULL
}

res.pca <- PCA(X = newLangdata, ncp = ncomp ,scale.unit = TRUE, quali.sup = NULL , graph = FALSE, ind.sup = dropSubjects)

# The loadings, defined as the sqrt of eigenvalues multplity by the eigenvectors are in *$var$coord
res.pca.loadings <- res.pca$var$coord

# For easier interpretation of principal directions, a varimax is performed
rotatedLoadings <- varimax(x = res.pca.loadings, normalize = TRUE)$loadings

##### ILLUSTATION PLOT ####

# Correlation plot of the loadings matrix
dev.new()
#tiff(filename = './Illustration/Raw_loadings.tiff', width = 20, height = 10, units = 'in', res = 360)
#tiff(filename = './Illustration/Raw_loadings_matrix_no_phono_lowres.tiff', width = 1920, height = 1080, units = 'px')
corrplot(res.pca.loadings, method = 'number', title = 'Raw loadings', tl.cex = 2, col = col3(100), tl.col = 'black', number.cex = 2)
#dev.off()

dev.new()
corrplot(rotatedLoadings, method = 'number', title = 'Loadings after varimax rotation')
# Print a summary 
summary(res.pca)

# Variables contribution
dev.new()
fviz_contrib(X = res.pca, choice = 'var', axes = 1)
dev.new()
fviz_contrib(X = res.pca, choice = 'var', axes = 2)
dev.new()
fviz_contrib(X = res.pca, choice = 'var', axes = 3)


# Projection of individuals on retained components
dev.new()
fviz_pca_ind(res.pca,
             col.ind = groups, # colors by groups
             palette = randomColor(length(unique(groups)), hue = "random", luminosity = "bright"),
             addEllipses = FALSE,
             legend.title = "Groups",
             ellipse.type = 'convex',
             repel = TRUE,
             axes = c(1,2),
             title = 'Projection of subjects on PC1 and PC2'
)


#### CLUSTER ANALYSIS ON PRINCIPAL COMPONENTS FOR INTERPRETATION PURPOSE ONLY ####
axes_to_cluster <- c(1,2)
xvar <- "PC1"
yvar <- "PC2"
# Determines the optimal number of cluster
PC_to_cluster <- cluster_on_pca(pca.object = res.pca,axes = axes_to_cluster, center=FALSE, scl=FALSE)$PCdata[,1:2] # Fetch the component to cluster

# Elbow method
dev.new()
fviz_nbclust(PC_to_cluster, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
# Silhouette method
dev.new()
fviz_nbclust(PC_to_cluster, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
dev.new()
fviz_nbclust(PC_to_cluster, kmeans, nstart = 3,  method = "gap_stat", nboot = 500)+
  labs(subtitle = "Gap statistic method")

# With NbClut lib
dev.new()
nb <- NbClust(PC_to_cluster, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans")
dev.new()
fviz_nbclust(nb)

# Cluster analysis
ncluster <- 2
nstart <- 1
PCcluster <- cluster_on_pca(pca.object = res.pca, axes = axes_to_cluster, nstart = nstart, ncluster = ncluster, scl = FALSE, center =FALSE)$PCdata
PCclusterPlot <- plot_prcomp_cluster(PC.data.frame = PCcluster, xvar = xvar, yvar = xvar,luminosity = "bright")

dev.new()
plot_prcomp_cluster(PC.data.frame = PCcluster, xvar = xvar , yvar = yvar,luminosity = "bright")


dev.new()
fviz_pca_ind(res.pca,
             col.ind = as.factor(PCcluster$V3), # colors by groups
             palette = randomColor(6),
             addEllipses = TRUE,
             legend.title = "Groups",
             ellipse.type = 'convex',
             repel = TRUE,
             axes = axes_to_cluster,
             title = 'Projection of subjects on PC1 and PC2 according to cluster groups'
)





###### Exploration of data #########
PC_ <- data.frame(PC2=res.pca$ind$coord[,1], PC3=res.pca$ind$coord[,2],Comphrension=groups)
#PC_ <- data.frame(PC2=pca.psych$scores[,1], PC3=pca.psych$scores[,3],Comphrension=groups)
#rownames(PC_) <- rownames(res.pca$ind$coord)
#dev.new()
#tiff(filename = './Illustration/Parole_Langage_PC1_PC2_no_phono.tiff', width = 20, height = 10, units = 'in', res = 360)
pdf('/home/db242421/CPM_results_23_05_2018/LG_LDFLIP_gender_lesion_0.02/Lesion_side_PC1_PC2_no_phono_lowres.pdf')
ggplot(data = PC_, fill=groups) + 
  geom_point(mapping = aes(x =PC2, y = PC3, color = factor(groups)), size =3) +
  geom_vline(xintercept = 0, linetype ='dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_text_repel(mapping = aes(x = PC2, y = PC3,
                                color = factor(groups),
                                label = rownames(PC_)),
                  size = 4,
                  fontface = 'bold',
                  box.padding = 0.5
                  ) +
  scale_color_discrete(name = 'Lesion side') +
  xlab("Language profil scores") + 
  ylab("Speech profil scores") +
  ggtitle("Projection of subjects in the speech and language profil plan", subtitle = "legend according to lesion side") + 
  theme_classic() + 
  theme(plot.title = element_text(size=15), legend.position="bottom", legend.text = element_text(size=10)) + 
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()  

# Distribution of PC scores
hist(res.pca$ind)
lines(density(res.pca$ind$coord[,1]), col="blue", lwd=2)

# Linear Model fit 
monograme <- langdataxls$monogramme[1:61]
WISC <- as.numeric(langdataxls$WISC[1:61])
PC1 <- res.pca$ind$coord[,1]
PC2 <- res.pca$ind$coord[,2]
# Construct dataframe for convenience
Profil_Langagier <- factor(langdataxls$Profil_langagier[1:61])
wisc.dataframe <- as.data.frame(cbind(WISC, PC1, PC2, Profil_Langagier))

# Basic scatter plot of WISC according to PC1
dev.new()
ggplot2.scatterplot(data = wisc.dataframe, xName ="PC1", yName ="WISC",
            groupName="Profil_Langagier", mainTitle = paste("WISC score according to first principal components scores", ", Pearson correlation = ",round(cor(WISC,PC1),digits = 2))) +
                     scale_color_manual(labels = c("Atypique", "Typique"), values = c("#ff6666","#996600")) + stat_smooth(method = "lm", se = FALSE, colour = 'purple')


lm.WISC <- lm(data = wisc.dataframe,formula = 'WISC~PC1')
dev.new()
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(lm.WISC, las = 1)



library(pracma)
# Manual computation of the scores after rotation
pinvLoadings <- t(pinv(rotatedLoadings)) # transpose of pseudo inverse of rotated loadings
rotatedScores <- scale(newLangdata_onlynum) %*% pinvLoadings # the rotated scores are standardized scores times the transpose pinv of
# rotated loading.
ScoresVerif <- scale(res.pca$ind$coord) %*% varimax(res.pca.loadings)$rotmat # Or simpler, there are the stantardized scores times the
# rotation matrix.


res.pca$ind$coord[,1] <- rotatedScores[,1]
res.pca$ind$coord[,2] <- rotatedScores[,2]
res.pca$ind$coord[,3] <- rotatedScores[,3]

dev.new()
fviz_pca_ind(res.pca,
              col.ind = groups, # colors by groups
              palette = c("#00AFBB", "#E7B800","#551A8B"),
              addEllipses = TRUE,
              legend.title = "Groups",
              repel = TRUE,
              axes = c(1,2),
              ellipse.type = 'convex',
              title = 'Scores after varimax rotation')



