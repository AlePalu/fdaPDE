library(femR)

## a simple script for an anisotropic diffusion problem

## create square domain
n <- 10
x.2D <- seq(0, 1, length.out = n)
y.2D <- x.2D
locations.2D <- expand.grid(x.2D, y.2D)
mesh.2D <- fdaPDE::create.mesh.2D(locations.2D)

## domain definition
D <- list(
    "nodes"    = mesh.2D$nodes,
    "edges"    = mesh.2D$edges,
    "elements" = mesh.2D$triangles,
    "neigh"    = mesh.2D$neighbors,
    "boundary" = mesh.2D$nodesmarkers
)

## use order 1 finite elements
f <- fe_function(fe_order = 1)
## define differential operator
K = matrix(c(1, 2, 2, 1), nrow = 2) ## diffusion tensor
L = div(K*grad(f))
## forcing term
u <- function(points) {
    return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]))
}

## create pde Lf = u on D
PDE <- pde(D, L, u)
## add boundary conditions
zero_bc <- as.matrix(rep(0., times = nrow(D$nodes)))
PDE$set_dirichlet_bc(zero_bc)

## solve PDE and get solution 
PDE$solve()
PDE$solution()
