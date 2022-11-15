cat("######## ORDER 1 #######\n\n")
library(Rcpp)
library(roxygen2)
roxygenise()

exact_solution <- function(points){
    return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]) )
}

forcing <- function(points){
    return(8.*pi^2* sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) ) 
}

PDE <- new(PDE_2D_isotropic_ORDER_1, femR::unit_square)

# attenzione al tipo di "transport"
PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- as.matrix(rep(0., times = dim(femR::unit_square$nodes)[1]))
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes)
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(femR::unit_square$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

cat("########################\n\n")

cat("######## ORDER 2 #######\n\n")
library(fdaPDE)
mesh <- create.mesh.2D(nodes = femR::unit_square$nodes, order = 2)
square <-list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)

# FAILS !!!
PDE <- new(PDE_2D_isotropic_ORDER_2, square)

# attenzione al tipo di "transport"
PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- as.matrix(rep(0., times = dim(femR::unit_square$nodes)[1]))
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes)
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(femR::unit_square$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")
cat("########################\n\n")

