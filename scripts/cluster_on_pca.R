cluster_on_pca <- function(pca.object ,axes = c(1,2), nstart = 1, ncluster = 1, center = TRUE, scl = TRUE){
  # Fetch the chosen PC
  PCplan <- scale(cbind(pca.object$ind$coord[, axes[1]], pca.object$ind$coord[, axes[2]]), center = center, scale = scl)
  # Performs a cluster analysis with kmeans
  km <- kmeans(x = PCplan, nstart = nstart, centers = ncluster)
  # Fetch the label cluster for each subject
  km.cluster <- km$cluster
  # Construct the dataframe with a "Cluster" columns containing the cluster
  PC.data.frame <- as.data.frame(cbind(PCplan,km$cluster))
  
  output <- list(PCdata=PC.data.frame, ClusterList=km$cluster)
  return(output)
}