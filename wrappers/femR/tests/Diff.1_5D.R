cat("########   Diff  #######\n\n")
cat("######## ORDER 1 #######\n\n")

library(femR)

rm(list=ls())

source("tests/utils.R")
data("unit_interval_graph", package="femR")

exact_solution <- function(x){
    return( sin(x))
}

forcing <- function(x){
    return( sin(x) ) 
}

PDE <- new(PDE_1_5D_isotropic_ORDER_2, unit_interval_graph)

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- exact_solution(PDE$get_dofs_coordinates())
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes[,1])
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(PDE$get_dofs_coordinates()[,1]))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

cat("############ Convergence ############\n\n")
rm(list=ls())

library(femR)
library(fdaPDE)
source("tests/utils.R")

data("unit_interval_graph", package="femR")
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
PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.1.5D(mesh = mesh)
    
    interval <- list(nodes= as.matrix(mesh$nodes[,1]), edges= mesh$edges, elements= mesh$edges, neigh= sparseNeigh(mesh), boundary= mesh$nodesmarkers)
    h[i] <- delta/(2^i)
    
    PDE <- new(PDE_1_5D_isotropic_ORDER_1, unit_interval_graph)
    PDE$set_PDEparameters(PDE_parameters)
    
    dirichletBC <- exact_solution(interval$nodes)
    PDE$set_dirichletBC(dirichletBC)

    quadrature_nodes <- PDE$get_quadrature_nodes()
    cat("################ ", i , " ################\n")
   
    f <- forcing(quadrature_nodes[,1])
    PDE$set_forcingTerm(as.matrix(f))
    result <- PDE$solve()
    
    u_ex <- as.matrix(exact_solution(PDE$get_dofs_coordinates()[,1]))

    errors.L2[i] <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
    cat("L2 error = ", errors.L2[i], "\n")
}
    
q = log2(errors.L2[1:(N-1)]/errors.L2[2:N])
cat("order = ", q, "\n")

cat("#####   Reac-Diff  #####\n\n")
cat("######## ORDER 1 #######\n\n")

library(femR)
library(fdaPDE)
source("tests/utils.R")

data("unit_interval_graph", package="femR")
C_ = 1/(exp(-1) - exp(1))
exact_solution <- function(x){
    return( C_ * exp(-x) - C_ * exp(x))
}

forcing <- function(x){
    return( 0.*x^0 ) 
}

PDE <- new(PDE_1_5D_isotropic_ORDER_1, unit_interval_graph)

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 1.)
PDE$set_PDEparameters(PDE_parameters)

dirichletBC <- exact_solution(PDE$get_dofs_coordinates())
PDE$set_dirichletBC(dirichletBC)

quadrature_nodes <- PDE$get_quadrature_nodes()
f <- forcing(quadrature_nodes[,1])
PDE$set_forcingTerm(as.matrix(f))

result <- PDE$solve()

u_ex <- as.matrix(exact_solution(unit_interval_graph$nodes))

error.L2 <- sqrt(sum(result$Mass %*% (u_ex - result$solution)^2))
cat("L2 error = ", error.L2, "\n")

library(latex2exp)
{
x11()
plot(unit_interval_graph$nodes[,1], result$solution[,1], col="blue", lwd = 1, pch=16)
points(unit_interval_graph$nodes[,1], u_ex, col="black", lwd=3, pch=21)
legend("topleft", legend=c(TeX("$u_{h}$"), TeX("$u_{ex}$")), 
        col=c("blue", "black"), 
        pch = c(16,21), 
        cex=1.25)
}

cat("############ Convergence ############\n\n")
rm(list=ls())
library(femR)
library(fdaPDE)
source("tests/utils.R")

C_ = 1/(exp(-1) - exp(1))
exact_solution <- function(x){
    return( C_ * exp(-x) - C_ * exp(x))
}

forcing <- function(x){
    return( 0.*x^0 ) 
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
PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 1.)

for(i in 1:N){
    
    mesh <- fdaPDE::refine.by.splitting.mesh.1.5D(mesh = mesh)
    
    interval <- list(nodes= as.matrix(mesh$nodes[,1]), edges= mesh$edges, elements= mesh$edges, neigh= sparseNeigh(mesh), boundary= mesh$nodesmarkers)
    h[i] <- delta/(2^i)
    
    PDE <- new(PDE_1_5D_isotropic_ORDER_1, interval)
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

imgdir_ = paste(imgdir_, "1_5D/", sep = "")

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

# cat("########################\n\n")
# library(fdaPDE)
# x = seq(from=0, to=1, length.out=3)
# y = rep(0, times=3)
# nodes = cbind(x,y)
# edges = matrix(c(1,2,2,3), nrow=2, ncol=2, byrow=T)
# mesh = create.mesh.1.5D(nodes=nodes, edges=edges)

# for(i in 1:5){
# mesh = refine.by.splitting.mesh.1.5D(mesh)
# }
# nodes = mesh$nodes
# edges = mesh$edges
# boundary = as.integer(mesh$nodesmarkers)
# nedges = nrow(mesh$nodes)
# neigh = Matrix(0, nrow = nedges, ncol=nedges)
# for(i in 1:nrow(mesh$edges)){ neigh[ mesh$edges[i,1], mesh$edges[i,2] ] = 1; neigh[ mesh$edges[i,2], mesh$edges[i,1] ] = 1}
# neigh <- as(neigh, "TsparseMatrix")

# neigh.TMP = matrix(0, nrow=length(neigh@i), ncol=3)
# for(k in 1:length(neigh@i)){
#     neigh.TMP[k, ] = c( neigh@i[k], neigh@j[k], neigh@x[k] )
# }

# neigh.TMP[1:nrow(neigh.TMP), 1:2] = neigh.TMP[1:nrow(neigh.TMP), 1:2] + matrix(1, nrow=nrow(neigh.TMP), ncol=2)
# neigh.TMP

# neigh = neigh.TMP
# storage.mode(neigh) <- "integer"
# #neigh

# unit_interval_graph = list( "nodes" = matrix(nodes[,1], nrow=nrow(mesh$nodes), ncol=1), 
#                       "edges" = edges,
#                       "elements" = edges, 
#                       "neigh" = neigh,
#                       "boundary" = matrix(boundary, nrow=length(boundary), ncol=1))
# save(unit_interval_graph, file = paste(getwd(), "/data/unit_interval_graph.RData",sep=""))

#############################
# library(fdaPDE)
# x = seq(from= 0, to=1, length.out= 61)
# nodes = cbind(x, rep(0, times= length(x)))
# edges = cbind( seq(from=1, to=60, by=1), seq(from=2, to= 61, by=1))
# mesh = create.mesh.1.5D(nodes= nodes, edges= edges)
# elements= edges
# boundary= mesh$nodesmarkers
# neigh = matrix(0, nrow= nrow(mesh$edges), ncol=2)

# mesh$neighbors[[1,2]] = -1
# mesh$neighbors[[60,1]] = -1

# for(i in 1:nrow(mesh$edges))
#     for(j in 1:2)
#         neigh[i,j] = mesh$neighbors[[i,j]]

# unit_interval = list(nodes= as.matrix(nodes[,1]), edges= edges, 
#                      elements= elements, neigh=neigh,
#                      boundary= boundary)

# save(unit_interval, file="data/unit_interval.RData")
