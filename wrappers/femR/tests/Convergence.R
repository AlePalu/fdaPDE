library(latex2exp)
library(femR)
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

N <- 5
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)
nnodes <- rep(0, times = N)
mesh <- fdaPDE::create.mesh.2D(nodes = matrix(c(0,0,1,0,1,1,0,1), nrow=4,ncol=2, byrow=T))
mesh <- fdaPDE::refine.mesh.2D(mesh=mesh, minimum_angle=30, maximum_area=0.1)

cat("############ ORDER 1 #############\n\n")
for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.2D(mesh = mesh)
    
    square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
    h[i] <- compute_h(mesh)$h_max #1/nrow(square$elements)
    nnodes[i] <- nrow(mesh$nodes)
    PDE <- new(PDE_2D_isotropic_ORDER_1, square)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- as.matrix(rep(0., times = dim(square$nodes)[1]))
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("############ ", i , " ############\n")
   
    f <- forcing(quadrature_nodes)
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(mesh$nodes))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}


q = log2(errors.L2[1:(N-1)]/errors.L2[2:N]) # /log2(h[1:(N-1)]/h[2:N]) 
cat("order = ", q, "\n")

p = log10(errors.L2[1:(N-1)]/errors.L2[2:N])/log10(nnodes[2:N]/nnodes[1:(N-1)])
cat("order (nodes) = ", p, "\n")
imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)

domain_ = "2D/"
imgdir_ = paste(imgdir_,domain_,sep="")

if(!dir.exists(imgdir_))
    dir.create(imgdir_)
    
pdf(paste(imgdir_,"diffusion_rates_order_1.pdf",sep=""))
plot(log2(h), log2(errors.L2), col="red", type="b", pch =16, lwd = 4, lty = 2, cex = 3,
        ylim = c(min(log2(h^2), log2(errors.L2)), max(log2(h), log2(errors.L2))),
        xlab = TeX("$log_2(h)$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 4, lty =2, cex = 3 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^2$")), 
        col=c("red", "black"), 
        lty = 2, 
        lwd = 6,
        cex=1.25)
dev.off()

pdf(paste(imgdir_,"diffusion_rates_order_1_nodes.pdf",sep=""))
plot(log10(nnodes), log10(errors.L2), col="red", type="b", pch =16, lwd = 4, lty = 2, cex = 3,
        ylim = c(min(log10(nnodes^(-1)), log10(errors.L2)), max(log10(nnodes^(-1)), log10(errors.L2))),
        xlab = TeX("$log_{10}(nodes)$"), ylab="", cex.lab=1.25)
lines(log10(nnodes), log10(nnodes^(-1)), col = "black", type = "b", pch = 16, lwd = 4, lty =2, cex = 3 )
legend("topright", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$nodes^{-1}$")), 
        col=c("red", "black"), 
        lty = 2,
        lwd=6, 
        cex=1.25)
dev.off()

cat("##################################\n\n")

FEMbasis = create.FEM.basis(mesh = mesh)
FEM.1 = fdaPDE::FEM(coeff = result$solution, FEMbasis = FEMbasis)
FEM.2 = fdaPDE::FEM(coeff = u_ex, FEMbasis = FEMbasis)
plot(FEM.1)
plot(FEM.2)
