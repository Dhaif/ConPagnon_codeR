plot_prcomp_cluster <- function(PC.data.frame, xvar,yvar, title = "Cluster analysis on Principal components" ,luminosity="light"){

  # Plot the results in the Principal component plan with the label corresponding to the attributed cluster number
  #dev.new()
  # generate random set of color
  randCol <- randomColor(length(unique(PC.data.frame$V3)), luminosity = luminosity)
  p <- ggplot(data = PC.data.frame, aes(x = V1, y = V2, group=factor(V3)))+ geom_point(aes(colour=factor(V3), fill=factor(V3)),shape=21)+scale_fill_manual(values = randCol)  + 
    geom_vline(xintercept=c(0), linetype="dotted",size=1) + geom_hline(yintercept = c(0), linetype="dotted",size=1) +geom_label_repel(aes(x = V1 , y = V2 , fill = factor(V3),
                                                                                                                                          label = rownames(PC.data.frame)),
                                                                                                                                      fontface = 'bold', color = 'white',
                                                                                                                                      box.padding = 0.35, point.padding = 0.5,
                                                                                                                                      segment.color = 'grey50') + theme_classic(base_size = 16) + xlab(xvar) + ylab(yvar) + labs(title = title)
  
  
  output <- list(ClusterPlot=p, colorPlot=randCol)
  return(output)
}