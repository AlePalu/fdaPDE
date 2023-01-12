cat("######## Reaction-Diffusion ########\n")
cat("########       ORDER 2      ########\n\n")
rm(list=ls())
library(Rcpp)
library(roxygen2)
library(latex2exp)
roxygenise()

library(fdaPDE)
rm(list=ls())
source("tests/utils.R")

exact_solution <- function(points){
    return( exp(points[,1] + points[,2]) )
}

forcing <- function(points){
    return(0.) 
}

PDE <- new(PDE_2D_isotropic_ORDER_2, femR::unit_square)

PDE_parameters <- list("diffusion" = 1., "transport" = rbind(0., 0.), "reaction" = 2.)
PDE$set_PDEparameters(PDE_parameters)

nodes <- PDE$get_dofs_coordinates()
dirichletBC <- as.matrix(exact_solution(nodes))
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- rep(0., times= nrow(quadrature_nodes))
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(nodes))
error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

cat("#####################################\n\n")

cat("############ Convergence ############\n\n")
N <- 4
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T))
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.1)
for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.2D(mesh = mesh)
    
    square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
    h[i] <- compute_h(mesh)$h_max
    
    PDE <- new(PDE_2D_isotropic_ORDER_2, square)
    PDE$set_PDEparameters(PDE_parameters)
    
    nodes <- PDE$get_dofs_coordinates()
    dirichletBC <- as.matrix(exact_solution(nodes))
    PDE$set_dirichletBC(dirichletBC)

    cat("################ ", i , " ################\n")
    
    quadrature_nodes <- PDE$get_quadrature_nodes()
    f <- rep(0., times= nrow(quadrature_nodes))
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(nodes))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)

pdf(paste(imgdir_,"reac_diff_rates_order_2.pdf",sep=""))
plot(log2(h), log2(errors.L2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h^3), log2(errors.L2)), max(log2(h^3), log2(errors.L2))+2),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^3), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^3$")), 
        col=c("red", "black"), 
        lty = 2, 
        cex=1.25)
dev.off()
