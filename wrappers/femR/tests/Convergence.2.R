cat("############ ORDER 2 #############\n\n")
library(Rcpp)
library(roxygen2)
roxygenise()

library(fdaPDE)

exact_solution <- function(points){
    return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]) )
}

forcing <- function(points){
    return(8.*pi^2* sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) ) 
}

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)

N = 3
errors.L2 <- rep(0, times = N)

mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T), order = 2)
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.05)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.mesh.2D(mesh = mesh, minimum_angle=30, maximum_area=0.05/(i+1))
    
    square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
    
    PDE <- new(PDE_2D_isotropic_ORDER_2, square)
    PDE$set_PDEparameters(PDE_parameters)
    
    #_Map_base::at
    #dirichletBC <- as.matrix(rep(0., times = dim(square$nodes)[1]))
    #PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("############ ", i , " ############\n")
   
    f <- forcing(quadrature_nodes)
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(square$nodes))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")
cat("##################################\n\n")
