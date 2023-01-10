cat("######## Advection-Diffusion ########\n")
cat("########       ORDER 1       ########\n\n")
rm(list=ls())
library(Rcpp)
library(roxygen2)
library(latex2exp)
roxygenise()

library(fdaPDE)
rm(list=ls())
source("tests/utils.R")

W_ <- 1.
R_ <- 1.
H_ <- 1.
beta_ <- 1.

alpha_ <- H_ * beta_ / R_
gamma_ <- pi * W_ / R_  

lambda1 <- -alpha_/2 - sqrt((alpha_/2)^2 + pi^2)
lambda2 <- -alpha_/2 + sqrt((alpha_/2)^2 + pi^2)

p_ <- (1-exp(lambda2))/(exp(lambda1)-exp(lambda2))

exact_solution <- function(points){
    return( -gamma_/pi^2 * ( p_ * exp(lambda1 * points[,1]) + (1 - p_) * exp(lambda2 * points[,1]) - 1. ) * sin(pi * points[,2] ) )
}

forcing <- function(points){
    return(gamma_ * sin(pi * points[,2])) 
}

PDE <- new(PDE_2D_isotropic_ORDER_1, femR::unit_square)

PDE_parameters <- list("diffusion" = 1., "transport" = rbind(-alpha_, 0.), "reaction" = 0.)
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

cat("#####################################\n\n")

cat("############ Convergence ############\n\n")
N <- 5
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T))
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.1)
for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.2D(mesh = mesh)
    
    square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
    h[i] <- compute_h(mesh)$h_max
    
    PDE <- new(PDE_2D_isotropic_ORDER_1, square)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- as.matrix(rep(0., times = dim(square$nodes)[1]))
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("################ ", i , " ################\n")
   
    f <- forcing(quadrature_nodes)
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(mesh$nodes))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)

pdf(paste(imgdir_,"adv_diff_rates_order_1.pdf",sep=""))
plot(log2(h), log2(errors.L2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h^2), log2(errors.L2)), max(log2(h), log2(errors.L2))+2),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^2$")), 
        col=c("red", "black"), 
        lty = 2, 
        cex=1.25)
dev.off()

cat("#####################################\n\n")

cat("####### Advection-dominanted ########\n\n")
W_ <- 1.
R_ <- 1.
H_ <- 1.
beta_ <- 2e2
alpha_ <- H_ * beta_ / R_

gamma_ <- pi * W_ / R_  
lambda1 <- -alpha_/2 - sqrt((alpha_/2)^2 + pi^2)
lambda2 <- -alpha_/2 + sqrt((alpha_/2)^2 + pi^2)
p_ <- (1-exp(lambda2))/(exp(lambda1)-exp(lambda2))

PDE_parameters <- list("diffusion" = 1., "transport" = rbind(-alpha_, 0.), "reaction" = 0.)

N <- 3
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T))
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.1)
for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.2D(mesh = mesh)
    
    square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
    h[i] <- 1/nrow(square$elements)
    
    PDE <- new(PDE_2D_isotropic_ORDER_1, square)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- as.matrix(rep(0., times = dim(square$nodes)[1]))
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("################ ", i , " ################\n")
   
    f <- forcing(quadrature_nodes)
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(mesh$nodes))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
    
    FEMbasis = create.FEM.basis(mesh = mesh)
    FEM.1 = fdaPDE::FEM(coeff = result$solution, FEMbasis = FEMbasis)
    FEM.2 = fdaPDE::FEM(coeff = u_ex, FEMbasis = FEMbasis)
    plot(FEM.1)
    plot(FEM.2)
}
     
cat("#####################################\n\n")

