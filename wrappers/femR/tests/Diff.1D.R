cat("###  Diff-Adv-Reac (exact) ###\n\n")
cat("###         ORDER 1        ###\n\n")

library(femR)

rm(list=ls())

source("tests/utils.R")
data("unit_interval", package="femR")

exact_solution <- function(x){
    return( (1. - x))
}

forcing <- function(x){
    return(x) 
}

PDE <- new(PDE_1D_isotropic_ORDER_1, unit_interval)

PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(-1.), "reaction" = -1.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- exact_solution(PDE$get_dofs_coordinates())
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes)
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(PDE$get_dofs_coordinates()))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

library(latex2exp)
{
x11()
plot(PDE$get_dofs_coordinates(), result$solution, col="blue", lwd=1, pch=16)
points(PDE$get_dofs_coordinates(), u_ex, col="black", lwd=3, pch=21)
legend("topleft", legend=c(TeX("$u_{h}$"), TeX("$u_{ex}$")), 
        col=c("blue", "black"), 
        pch = c(16,21), 
        cex=1.25)
}


cat("########   Diff  #######\n\n")
cat("######## ORDER 1 #######\n\n")

library(femR)

rm(list=ls())

source("tests/utils.R")
data("unit_interval", package="femR")

exact_solution <- function(x){
    return( sin(x))
}

forcing <- function(x){
    return( sin(x) ) 
}

PDE <- new(PDE_1D_isotropic_ORDER_1, unit_interval)

PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(0.), "reaction" = 0.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- exact_solution(PDE$get_dofs_coordinates())
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes[,1])
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(unit_interval$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

library(latex2exp)
{
x11()
plot(unit_interval$nodes[,1], result$solution[,1], col="blue", lwd = 1, pch=16)
points(unit_interval$nodes[,1], u_ex, col="black", lwd=3, pch=21)
legend("topleft", legend=c(TeX("$u_{h}$"), TeX("$u_{ex}$")), 
        col=c("blue", "black"), 
        pch = c(16,21), 
        cex=1.25)
}

cat("############ Convergence ############\n\n")
library(femR)
library(fdaPDE)

rm(list=ls())
source("tests/utils.R")

exact_solution <- function(x){
    return( sin(x))
}

forcing <- function(x){
    return( sin(x) ) 
}

N <- 4
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

x = seq(from=0, to=1, length.out=3)
y = rep(0, times=3)
nodes = cbind(x,y)
edges = matrix(c(1,2,2,3), nrow=2, ncol=2, byrow=T)
mesh = create.mesh.1.5D(nodes=nodes, edges=edges)

delta = 0.25
mesh <- refine.mesh.1.5D(mesh, delta)
PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(0.), "reaction" = 0.)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.1.5D(mesh = mesh)
    
    interval <- list(nodes= as.matrix(mesh$nodes[,1]), edges= mesh$edges, elements= mesh$edges, neigh= neigh1D(mesh), boundary= mesh$nodesmarkers)
    h[i] <- delta/(2^i)
    
    PDE <- new(PDE_1D_isotropic_ORDER_1, interval)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- exact_solution(interval$nodes)
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("################ ", i , " ################\n")
   
    f <- forcing(quadrature_nodes[,1])
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(mesh$nodes[,1]))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

cat("####  Diff-Adv-Reac ###\n\n")
cat("####      ORDER 1   ###\n\n")
library(femR)

rm(list=ls())
data("unit_interval", package="femR")
source("tests/utils.R")

exact_solution <- function(x){
    return( exp(x)*sin(x))
}

forcing <- function(x){
    return(  exp(x)*sin(x) ) 
}

PDE <- new(PDE_1D_isotropic_ORDER_1, unit_interval)

PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(2.), "reaction" = -1.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- exact_solution(PDE$get_dofs_coordinates())
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes[,1])
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(unit_interval$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

library(latex2exp)
{
x11()
plot(unit_interval$nodes[,1], result$solution[,1], col="blue", lwd = 1, pch=16)
points(unit_interval$nodes[,1], u_ex, col="black", lwd=3, pch=21)
legend("topleft", legend=c(TeX("$u_{h}$"), TeX("$u_{ex}$")), 
        col=c("blue", "black"), 
        pch = c(16,21), 
        cex=1.25)
}

cat("############ Convergence ############\n\n")
library(femR)
library(fdaPDE)

rm(list=ls())

source("tests/utils.R")

exact_solution <- function(x){
    return( exp(x)*sin(x))
}

forcing <- function(x){
    return( exp(x)*sin(x) ) 
}

N <- 4
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

x = seq(from=0, to=1, length.out=3)
y = rep(0, times=3)
nodes = cbind(x,y)
edges = matrix(c(1,2,2,3), nrow=2, ncol=2, byrow=T)
mesh = create.mesh.1.5D(nodes=nodes, edges=edges)

delta = 0.25
mesh <- refine.mesh.1.5D(mesh, delta)
PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(2.), "reaction" = -1.)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.1.5D(mesh = mesh)
    
    interval <- list(nodes= as.matrix(mesh$nodes[,1]), edges= mesh$edges, elements= mesh$edges, neigh= neigh1D(mesh), boundary= mesh$nodesmarkers)
    h[i] <- delta/(2^i)
    
    PDE <- new(PDE_1D_isotropic_ORDER_1, interval)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- exact_solution(interval$nodes)
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("################ ", i , " ################\n")
   
    f <- forcing(quadrature_nodes[,1])
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(mesh$nodes[,1]))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)
    
imgdir_ = paste(imgdir_, "1D/", sep = "")
if(!dir.exists(imgdir_))
    dir.create(imgdir_)

pdf(paste(imgdir_,"reac_diff_rates_order_1.pdf",sep=""))
plot(log2(h), log2(errors.L2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h^2), log2(errors.L2)), max(log2(h^2), log2(errors.L2))+1),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^2$")), 
        col=c("red", "black"), 
        lty = 2, 
        cex=1.25)
dev.off()


cat("####  Diff-Adv-Reac ###\n\n")
cat("####      ORDER 2   ###\n\n")
library(femR)
rm(list=ls())

data("unit_interval", package="femR")
source("tests/utils.R")

exact_solution <- function(x){
    return( exp(x)*sin(x))
}

forcing <- function(x){
    return(  exp(x)*sin(x) ) 
}

PDE <- new(PDE_1D_isotropic_ORDER_2, unit_interval)

PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(2.), "reaction" = -1.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- exact_solution(PDE$get_dofs_coordinates())
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes[,1])
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(unit_interval$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

cat("############ Convergence ############\n\n")
library(femR)
library(fdaPDE)

rm(list=ls())
source("tests/utils.R")

exact_solution <- function(x){
    return( exp(x)*sin(x))
}

forcing <- function(x){
    return( exp(x)*sin(x) ) 
}

N <- 4
errors.L2 <- rep(0, times = N)
h <- rep(0, times = N)

x = seq(from=0, to=1, length.out=3)
y = rep(0, times=3)
nodes = cbind(x,y)
edges = matrix(c(1,2,2,3), nrow=2, ncol=2, byrow=T)
mesh = create.mesh.1.5D(nodes=nodes, edges=edges)

delta = 0.25
mesh <- refine.mesh.1.5D(mesh, delta)
PDE_parameters <- list("diffusion" = 1., "transport" = as.matrix(2.), "reaction" = -1.)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.1.5D(mesh = mesh)
    
    interval <- list(nodes= as.matrix(mesh$nodes[,1]), edges= mesh$edges, elements= mesh$edges, neigh= neigh1D(mesh), boundary= mesh$nodesmarkers)
    h[i] <- delta/(2^i)
    
    PDE <- new(PDE_1D_isotropic_ORDER_2, interval)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- exact_solution(interval$nodes)
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("################ ", i , " ################\n")
   
    f <- forcing(quadrature_nodes[,1])
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(PDE$get_dofs_coordinates()))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)
    
imgdir_ = paste(imgdir_, "1D/", sep = "")
if(!dir.exists(imgdir_))
    dir.create(imgdir_)

pdf(paste(imgdir_,"reac_diff_rates_order_2.pdf",sep=""))
plot(log2(h), log2(errors.L2), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h^2), log2(errors.L2)), max(log2(h^2), log2(errors.L2))+1),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25)
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c(TeX("$\\| u - u_{ex} \\|_{2}$"), TeX("$h^2$")), 
        col=c("red", "black"), 
        lty = 2, 
        cex=1.25)
dev.off()

