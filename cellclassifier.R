setGeneric("ClassifyCells", function(object, classifier, training.genes = NULL, training.classes = NULL, new.data = NULL, ... ) standardGeneric("ClassifyCells"))
#' @export
setMethod("ClassifyCells", "seurat",
          function(object, classifier, training.genes = NULL, training.classes = NULL, new.data = NULL, ...) {
            # build the classifier
            if (missing(classifier)){
              classifier <- BuildRFClassifier(object, training.genes = training.genes, training.classes = training.classes,...)
            }
            # run the classifier on the new data
            features <- classifier$forest$independent.variable.names
            genes.to.add <- setdiff(features, rownames(new.data))
            data.to.add <- matrix(0, nrow = length(genes.to.add), ncol = ncol(new.data))
            rownames(data.to.add) <- genes.to.add
	    #####*******Wen added the below 1 sentence**********##############
	    colnames(data.to.add)<-colnames(new.data)
            new.data <- rbind(new.data, data.to.add)
            new.data <- new.data[features, ]
            new.data <- as.matrix(t(new.data))
            cat("Running Classifier ...", file = stderr())
            prediction <- predict(classifier, new.data)
            new.classes <- prediction$predictions
            return(new.classes)
          }
)

#' Build Random Forest Classifier
#'
#' Train the random forest classifier 
#'
#'
#' @param object Seurat object on which to train the classifier
#' @param training.genes Vector of genes to build the classifier on
#' @param training.classes Vector of classes to build the classifier on
#' @param verbose Additional progress print statements
#' @param ... additional parameters passed to ranger
#' @return Returns the random forest classifier
#' @import Matrix
#' @importFrom ranger ranger 
#' @importFrom plyr mapvalues
#' @export
setGeneric("BuildRFClassifier", function(object, training.genes = NULL, training.classes = NULL, verbose = TRUE, ... ) standardGeneric("BuildRFClassifier"))
#' @export
setMethod("BuildRFClassifier", "seurat",
          function(object, training.genes = NULL, training.classes = NULL, verbose = TRUE, ...) {
            training.classes <- as.vector(training.classes)
            training.genes <- set.ifnull(training.genes, rownames(object@data))
            training.data <- as.data.frame(as.matrix(t(object@data[training.genes, ])))
            training.data$class <- factor(training.classes)
            if (verbose) cat("Training Classifier ...", file = stderr())
            classifier <- ranger(data = training.data, dependent.variable.name = "class", classification = T, 
                                 write.forest = T, ...)
            return(classifier)
          }
)



#' Highlight classification results
#'
#' This function is useful to view where proportionally the clusters returned from 
#' classification map to the clusters present in the given object. Utilizes the FeaturePlot()
#' function to color clusters in object.
#'
#'
#' @param object Seurat object on which the classifier was trained and 
#' onto which the classification results will be highlighted
#' @param clusters vector of cluster ids (output of ClassifyCells)
#' @param ... additional parameters to pass to FeaturePlot()
#' @return Returns a feature plot with clusters highlighted by proportion of cells
#' mapping to that cluster
#' @export

############ATTENTION!!!!!!!!!!! the object "out" is just the parameter "clusters",so, the code has an error here, you shoule assign "out" with the results of ClassifyCells()###############
setGeneric("VizClassification", function(object, clusters, ... ) standardGeneric("VizClassification"))
#' @export
setMethod("VizClassification", "seurat",
          function(object, clusters, ...) {
            #####********Wen added the below one sentence*********#########
   	    out<-clusters
	    cluster.dist <- prop.table(table(out))
            object@data.info$Classification <- numeric(nrow(object@data.info))
            for (cluster in 1:length(cluster.dist)) {
              cells.to.highlight <- WhichCells(object, names(cluster.dist[cluster]))
              if(length(cells.to.highlight) > 0){
                object@data.info[cells.to.highlight, ]$Classification <- cluster.dist[cluster]
              }
            }
            if(any(grepl("cols.use", deparse(match.call())))){
              return(FeaturePlot(object, "Classification", ...))
            }
            cols.use = c("#f6f6f6", "black")
            return(FeaturePlot(object, "Classification", cols.use = cols.use, ...))
          }
)


calc.drop.prob=function(x,a,b) {
  return(exp(a+b*x)/(1+exp(a+b*x)))
}

