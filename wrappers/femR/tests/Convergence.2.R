cat("############ ORDER 2 #############\n\n")
library(Rcpp)
library(roxygen2)
library(latex2exp)
roxygenise()

library(fdaPDE)
rm(list=ls())
source("tests/utils.R")

exact_solution <- function(points){
    return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]) )
}

forcing <- function(points){
    return(8.*pi^2* sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) ) 
}

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)

N = 4
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T), order = 1)
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.05)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.2D(mesh = mesh)
    
    square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
    h[i] <- compute_h(mesh)$h_max
    
    PDE <- new(PDE_2D_isotropic_ORDER_2, square)
    PDE$set_PDEparameters(PDE_parameters)
    
    #dirichletBC <- as.matrix(rep(0., times = dim(square$nodes)[1]))
    #PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("############ ", i , " ############\n")
   
    f <- forcing(quadrature_nodes)
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    nodes <- PDE$get_dofs_coordinates()
    u_ex <- as.matrix(exact_solution(nodes))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)

pdf(paste(imgdir_,"diffusion_rates_order_2.pdf",sep=""))
plot(log(h), log(errors.L2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h^3), log2(errors.L2)), max(log2(h^3), log2(errors.L2))+2),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log(h), log(h^3), col = "black", type = "b", pch = 16, lwd = 3, lty = 2, cex = 2,)
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^3$")), 
        col=c("red", "black"), 
        lty = 2)
dev.off()

png(paste(imgdir_,"exact_solution.png",sep=""))
xx <- seq(from=0, to=1., length.out=100)
yy <- xx
grid <- expand.grid(xx,yy)
zz <- matrix(0, nrow=100, ncol=100)
for(i in 1:100)
    for( j in 1:100)
        zz[i,j] = exact_solution(matrix(c(xx[i], yy[j]), nrow=1, ncol=2) )
filled.contour(z=zz)
dev.off()
cat("##################################\n\n")
