# Relation between lesion and behavior
library(googlesheets)
library(ggplot2)
library(ggrepel)

wd <- "/media/db242421/db242421_data/ConPagnon_reports/resultsPCA"
setwd(wd)
# Get clinical data on goole drive
# Load behavioral dataframe.
data_gs_title <- gs_title("Resting State AVCnn: cohort data")
data <- as.data.frame(gs_read(ss = data_gs_title, ws="Middle Cerebral Artery resting state cohort"))
rownames(data) <- data$X1
data$X1 <- NULL
# Clean dataframe format
for(i in 1:ncol(x = data)){
  data[, i] <- gsub(",", ".", data[ , i])
}
data <- as.data.frame(data)
# Check everything is ok.
head(data)

# Read PCA on lesion map dataframe
pca_lesion_maps <- read.csv("pca_lesion_map.csv", row.names = "subjects")
scaled_pca_lesion_maps <- as.data.frame(scale(x = pca_lesion_maps, center = TRUE, scale = TRUE))


# Merge with overall dataframe 
new_data <- as.data.frame(merge(data, scaled_pca_lesion_maps, by = "row.names"))
rownames(new_data) <- new_data$Row.names
new_data$Row.names <- NULL

# Plot
library(RColorBrewer)
groups <- as.numeric(new_data$lesion_normalized)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(0, 0.05))
dev.new()
ggplot(data = new_data) + 
  geom_point(mapping = aes(x = PC1, 
                           y = PC2,
                           color = groups), 
             size = 5) + sc +
  geom_vline(xintercept = 0, linetype ='dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_text_repel(mapping = aes(x = PC1, 
                                y = PC2,
                                label = rownames(new_data),
                                color = groups),
                  size = 5,
                  fontface = 'bold',
                  box.padding = 0.5
  ) +
  xlab("PC1 (lesion map)") + 
  ylab("PC2 (lesion map)")+
  ggtitle("Projection of subjects in the first principal components plan from lesion map", 
          subtitle = paste(' ')) + 
  theme_classic() + 
  theme(plot.title = element_text(size=15, face="bold"), 
        legend.position="bottom", legend.text = element_text(size=10),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")) + 
  guides(colour = guide_legend(override.aes = list(size=5), title = "Legend"))


