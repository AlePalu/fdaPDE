compute_h <- function(mesh){
    h_min = 1e3
    h_max = 1-16
    
    for(e in 1:nrow(mesh$edges)){
        edge = mesh$edges[e,]
        h = norm(mesh$nodes[edge[1],]-mesh$nodes[edge[2],], "2")
        h_min = min(h_min, h)
        h_max = max(h_max, h)
    }
    return(list(h_min = h_min, h_max = h_max) )
}

sparseNeigh <- function(mesh) {
    nnodes <- nrow(mesh$nodes)
    neigh = Matrix(0, nrow = nnodes, ncol=nnodes)
    for(i in 1:nrow(mesh$edges)){ neigh[ mesh$edges[i,1], mesh$edges[i,2] ] = 1; neigh[ mesh$edges[i,2], mesh$edges[i,1] ] = 1}
    neigh <- as(neigh, "TsparseMatrix")

    neigh.TMP = matrix(0, nrow=length(neigh@i), ncol=3)
    for(k in 1:length(neigh@i)){
        neigh.TMP[k, ] = c( neigh@i[k], neigh@j[k], neigh@x[k] )
    }

    neigh.TMP[1:nrow(neigh.TMP), 1:2] = neigh.TMP[1:nrow(neigh.TMP), 1:2] + matrix(1, nrow=nrow(neigh.TMP), ncol=2)
    neigh.TMP

    neigh = neigh.TMP
    storage.mode(neigh) <- "integer"
    return(neigh)
}

neigh1D <- function(mesh){
    
    neigh = matrix(0, nrow=nrow(mesh$edges), ncol=2)
    for(i in 1:nrow(mesh$edges))
        for(j in 1:2)
            if( (i == 1 && j ==2) || ( i == nrow(mesh$edges) && j == 1))
                neigh[i,j] = -1
            else    
                neigh[i,j] = mesh$neighbors[[i,j]]

    return(neigh)
}