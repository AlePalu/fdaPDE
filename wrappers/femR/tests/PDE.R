library(Rcpp)
library(roxygen2)
roxygenise()

#library(femR)

PDE <- new(PDE_2D_isotropic_ORDER_1, femR::unit_square)

forcing <- function(points){
    return(4.*(points[,1]^2+points[,2]^2)*sin(2.*points[,1]*points[,2])*(points[,1]-1.0)*(points[,2]-1.) - 4.*cos(2.*points[,1]*points[,2])*(points[,1]*(points[,1]-1.) + points[,2]*(points[,2]-1.)))
}

# attenzione al tipo di "transport"
PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- as.matrix(rep(0., times = sum(femR::unit_square$boundary)))
PDE$set_dirichletBC(dirichletBC)


quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes)
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

print(dim(result$solution))
print(dim(result$Mass))


exact_solution <- function(points){
    return(sin(2.*points[,1]*points[,2])*(points[,1]-1.)*(points[,2]-1.))
}

u_ex <- as.matrix(exact_solution(femR::unit_square$nodes))

error.L2 <- sum(result$Mass %*% (u_ex - result$solution)^2)
print(paste("L2 error = ", error.L2, sep=""))
