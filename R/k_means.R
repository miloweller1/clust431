#' Implements a very basic k means algorithm using Lloyd's algorithm
#'
#' @param dat A dataset
#' @param k Number of cluster
#' @param pca default=F If true, automatically performs a PCA before doing k means clustering

#'
#' @return Cluster assignments, TOTSS, and a dataframe of cluster means
#'
#' @import dplyr
#' @export
#'
#'
k_means<- function(dat, K, pca=F) {

    #Write one function called k_means() that implements a very basic k-means algorithm.
    #Choose k random observations in the data as your starting points.
    #Do not do any fancy adjustments to balance cluster sizes and so forth.
    #Include an option in k_means() to automatically perform PCA before doing the k_means() clustering, using only the first 2 dimensions. (You may use built-in functions like princomp() for this.)
    #At a minimum, your function should output the cluster assignments and total sum of squares.

    if (isTRUE(pca)){

      pca_df = prcomp(dat, center = TRUE, scale = TRUE)

        X = as.matrix(pca_df$x[,1:2])
    }
    else{X=as.matrix(dat)}


    initialCentroids<-matrix(0,K,ncol(X))
    centroids<-X[sample(1:nrow(X),K),]
    set<-rep(0, K)
    Cluster<-rep(0, nrow(X))

    converged=F
    while(converged==F) {

         for(i in 1:dim(X)[1]){

               for(j in 1:dim(centroids)[1]){

                set[j]=(X[i,]-centroids[j,])%*%(X[i,]-centroids[j,]);
               }

            Cluster[i]=which.min(set);
        }

        initialCentroids=centroids;

        for(k in 1:K){
            centroids[k,]=colMeans(X[which(Cluster==k),]);
        }
        if (isTRUE(all.equal(centroids, initialCentroids))) {converged=T}
        else{converged=F}

    }
    #print(Cluster)
    obs<-dat
    sstotal<-sum((obs-centroids)^2)
    print(sstotal)
    results<-data.frame(rbind("clusters"=Cluster))
    colnames(results)<-rownames(dat)
    return(results)
}

