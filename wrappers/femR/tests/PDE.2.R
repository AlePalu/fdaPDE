cat("######## ORDER 2 #######\n\n")
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

mesh <- create.mesh.2D(nodes = femR::unit_square$nodes, order = 2)

PDE <- new(PDE_2D_isotropic_ORDER_2, femR::unit_square)

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)
PDE$set_PDEparameters(PDE_parameters)

#_Map_base::at
dirichletBC <- as.matrix(rep(0., times = dim(femR::unit_square$nodes)[1]))
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes)
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(mesh$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

FEMbasis = create.FEM.basis(mesh = mesh)
FEM.1 = fdaPDE::FEM(coeff = result$solution, FEMbasis = FEMbasis)
FEM.2 = fdaPDE::FEM(coeff = u_ex, FEMbasis = FEMbasis)
plot(FEM.1)
plot(FEM.2)

cat("########################\n\n")

