# Additional plot from whitin domain PCA analysis

wd <- "/media/db242421/db242421_data/ConPagnon_reports/resultsPCA"
domain <- "WISC"
setwd(file.path(wd,domain))
filename <- "WISC_pca.xlsx"
domain_data_path <- file.path(wd, domain, filename)

# Sheetname
sheet_name <- "scores"
# Read the file
domain_data <- as.data.frame(read_excel(domain_data_path, sheet = sheet_name))
rownames(domain_data) <- domain_data$subjects
domain_data$subjects <- NULL

# Read dataframe containing behavioral data
behavioral_dataframe <- as.data.frame(read.csv(file.path(wd, domain, paste(domain, "_pca_results_dataframe.csv", sep = ""))))
rownames(behavioral_dataframe) <- behavioral_dataframe$subjects
behavioral_dataframe$subjects <- NULL
column_to_merge <- c("Sexe", "Lesion", "cerebral_palsy", "Parole", "langage_clinique")

# Merge the dataframe to the domain dataframe
domain_data <- merge(x = domain_data, y = behavioral_dataframe[, column_to_merge],by = "row.names")
rownames(domain_data) <- domain_data$Row.names
domain_data$Row.names <- NULL


# Make the plot
# Figures parameters

groups <- domain_data$Lesion
dim_on_x <- 3
dim_on_y <- 13
points_labels <- rownames(domain_data)
size_of_points <- 3
size_of_points_labels <- 4
legend_title <- "Legend: "
legend_labels_size <- 10
points_labels_color <- groups

x_label <- "WISC sequences test socres"
y_label <- "PC1"
figure_title <- paste(domain, ": Relation between wisc sequence scores and PC1", sep = "")

figure_width <- 20
figure_heigth <- 10

save_results_directory <- file.path(wd, domain)

pdf(file = file.path(save_results_directory, 
                     paste(domain,"_",x_label,"_",y_label,"lesion_side.pdf")),
    width = figure_width,
    height = figure_heigth)
print(ggplot(data = domain_data) + 
        geom_point(mapping = aes(x = domain_data[, colnames(domain_data)[dim_on_x]], 
                                 y = domain_data[, colnames(domain_data)[dim_on_y]],
                                 color = groups), 
                   size =size_of_points) +
        geom_text_repel(mapping = aes(x = domain_data[, colnames(domain_data)[dim_on_x]], 
                                      y = domain_data[, colnames(domain_data)[dim_on_y]],
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

