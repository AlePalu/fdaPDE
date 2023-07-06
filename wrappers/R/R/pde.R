## create pde object backed by Cpp_pde_module
pde <- function(D, L, u) {
    ## NB: D should provide informations about mesh dimensions,
    ## for now we assume a 2D planar domain only!

    f <- L@f
    fe_order <- f@fe_order

    ## set pde type
    pde_type <- 0
    pde_parameters <- NULL
    if (length(L@tokens) == 1 && L@tokens[1] == "diffusion" &&
        length(L@params$diffusion) == 1 && L@params$diffusion == 1.0) {
        ## specialized implementation for laplace operator
        pde_type <- 1 
    } else {
        ## general diffusion-transport-reaction problem, constant coefficients
        pde_type <- 2
        ## prepare pde_parameters list
        pde_parameters$diffusion <- matrix(0, nrow = 2, ncol = 2)
        pde_parameters$transport <- matrix(0, nrow = 2, ncol = 1)
        pde_parameters$reaction  <- 0.0

        for(i in 1:length(L@tokens)) {
            pde_parameters[[paste(L@tokens[i])]] <- L@params[[paste(L@tokens[i])]]
        }
    }

    ## define Rcpp module
    if (fe_order == 1) { ## linear finite elements
        pde_ <- new(PDE_2D_ORDER_1, D, pde_type, pde_parameters)
    }
    if (fe_order == 2) { ## quadratic finite elements
        pde_ <- new(PDE_2D_ORDER_2, D, pde_type, pde_parameters)
    }
        
    ## evaluate forcing term on quadrature nodes
    quad_nodes <- as.matrix(pde_$get_quadrature_nodes())
    pde_$set_force(as.matrix(u(quad_nodes)))
    ## initialize and return
    pde_$init()
    return(pde_)
}

## plot.pde <- function() {

## }

## set_dirichlet_bc <- function(pde, g) {
##     ## get boundary nodes and evaluate boundary condition
## }

## ## convinient wrapper for solve() method
## solve.pde <- function(pde) { pde$solve() }

## ## mesh D ...

## ## linear finite element expansion
## f = fe_function(order = 1)

## ## differential operator
## K <- matrix(c(1,1,1,1), nrow = 2, ncol = 2)
## L <- div(K*grad(f))
## ## forcing term
## u <- function(p) { return(p[0] + p[1] + 2) }

## ## Lf = u
## pde <- pde(L, u)
## set_dirichlet_bc(pde, ...) ## set boundary conditions

## ## solve and plot
## solve(pde)
## plot(f)

## implementa -div

## `-.div` <- function() {
    
## } 

## L = div(K * grad()) + 4 * I() + dot(c(1, 1), grad())



## triangulation - pacchetto R
