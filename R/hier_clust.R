#' Implements a very basic single linkage hclust algorithm, return the same value as  cutree(hclust(dist(dat), method="single"), k=k)
#'
#' @param dat A dataset and desired number of cluster
#' @param k An integer to select desired number of clusters
#'
#'
#' @return Cluster assignments
#'
#' @import dplyr
#' @export
#'
#'

hier_clust<- function(dat, k){
N <- nrow(dat)
distMAT<- as.matrix(dist(dat))


 rownames(distMAT) <- -(1:N)
 colnames(distMAT) <- -(1:N)


 merge <- matrix(0, N-1, ncol=2)
 height <- rep(0,N-1)

 for (i in 1:(N-1)) {
     diag(distMAT) <- Inf
     cols <- colnames(distMAT)

     minD <- which(distMAT == min(distMAT), arr.ind = TRUE)[1,,drop=FALSE]

     height[i] <- min(distMAT)
     merge[i,] <- as.numeric(cols[minD])

     clusters <- c(minD, which(cols %in% cols[minD[1, cols[minD] > 0]]))
     colnames(distMAT)[clusters] <- i

     pairs <- apply(distMAT[minD,], 2, min)    #single linkage  method:mininum
     distMAT[,min(minD)] <- pairs
     distMAT[min(minD),] <- pairs


     diag(distMAT) <- Inf
     distMAT[,max(minD)] <- Inf
     distMAT[max(minD),] <- Inf

 }
 merge<-Thermimage::mirror.matrix(merge)

 res <- list("merged"=merge,"height"=height,"clusters"=clusters, "labels"=rownames(dat))
 cutree(res, k)



}
