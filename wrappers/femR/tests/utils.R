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
