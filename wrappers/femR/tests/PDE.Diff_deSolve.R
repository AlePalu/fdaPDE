library(ReacTran)
library(latex2exp)
library(roxygen2)
roxygenise()

rm(list=ls())
graphics.off()

rmse <- function(x,y){return(sqrt(mean( (x-y)^2)))}

#2D
exact <- function(x,y) sin(2*pi*x)*sin(2*pi*y)
f <- function(x,y) 8*pi^2*sin(2*pi*x)*sin(2*pi*y)
Dy    <- Dx <- -1.   # diffusion coeff, X- and Y-direction

PDE_parameters <- list("diffusion" = 1., "transport" = matrix(0.,nrow=2,ncol=1), "reaction" = 0.)

N = c(20,40,80,160) # int. nodes N[i]^2

errors.l2 <- list("deSolve" = rep(0, times = length(N)),
                  "femR" = rep(0, times = length(N)))

times <- list("deSolve" = rep(0, times = length(N)),
                  "femR" = rep(0, times = length(N)))
h <- rep(0, times = length(N))

for(i in 1:length(N)){
  cat("############ ", i , " ############\n")
  x.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  y.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  grid2D    <- setup.grid.2D(x.grid, y.grid)
  
  h[i] = max(grid2D$dx, grid2D$dy)

  D.grid    <- setup.prop.2D(value = Dx, y.value = Dy, grid = grid2D)
  v.grid    <- setup.prop.2D(value = 0, grid = grid2D)
  A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
  VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)

  C.x.up    <- rep(0, times = N[i])
  C.y.up    <- rep(0, times = N[i])

  C.x.down  <- rep(0, times = N[i])
  C.y.down  <- rep(0, times = N[i])

  forcing = matrix(nrow=N[i], ncol=N[i], data=0)
  y.ex = matrix(nrow=N[i], ncol=N[i], data=0)
  for(k in 1:N[i]){
    for(l in 1:N[i]){
      forcing[k,l] = f(grid2D$x.mid[k], grid2D$y.mid[l])
      y.ex[k,l] = exact(grid2D$x.mid[k], grid2D$y.mid[l])  
    }
  }

  Diff2Db <- function (t, y, parms)  {
  
    Y  <- matrix(nrow = N[i], ncol = N[i], data = y)
  
    dY <- tran.2D(C = Y, 
                     C.x.up = C.x.up, C.x.down = C.x.down,
                     C.y.up = C.y.up, C.y.down = C.y.down,
                     grid = grid2D, 
                     D.grid = D.grid, 
                     A.grid = A.grid, 
                     VF.grid = VF.grid, 
                     v.grid = v.grid)$dC
  
    dY <- dY - forcing 
  
    return (list(dY))
  }

  y = matrix(data = rep(1,times=N[i]*N[i]))
  y = matrix(data = rnorm(N[i]*N[i], mean = 1, sd =0.5))

  start_ = Sys.time()
  Y <- steady.2D(y=y,
                   dimens = c(N[i],N[i]), 
                   time = 0, 
                   func = Diff2Db, parms=NULL, lrw=1e8)
  times$deSolve[i] = as.numeric(difftime(Sys.time(), start_, units="secs"))                 
  cat("deSolve ", times$deSolve[i], " s\n")

  y.hat <- matrix(Y$y, nrow=N[i], ncol=N[i])
  errors.l2$deSolve[i] = rmse(Y$y, as.vector(y.ex))
  
  ###### femR ######
  grid2D <- expand.grid(x.grid$x.int, y.grid$x.int)
  mesh <- fdaPDE::create.mesh.2D(grid2D)
  square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
  PDE <- new(PDE_2D_isotropic_ORDER_1, square)
  PDE$set_PDEparameters(PDE_parameters)
    
  dirichletBC <- as.matrix(rep(0., times = dim(square$nodes)[1]))
  PDE$set_dirichletBC(dirichletBC)

  quadrature_nodes <- PDE$get_quadrature_nodes()
   
  PDE$set_forcingTerm(as.matrix(f(quadrature_nodes[,1], quadrature_nodes[,2])))
  start_ <- Sys.time()
  result <- PDE$solve()
  times$femR[i] <- as.numeric(difftime(Sys.time(), start_, units="secs"))
  cat("femR ", times$femR[i] , " s\n")  
  u_ex <- as.matrix(exact(mesh$nodes[,1], mesh$nodes[,2]))

  errors.l2$femR[i] <- rmse(u_ex, result$solution)
}

q = list("deSolve" = 0, "femR"= 0)
q$deSolve = log2(errors.l2$deSolve[1:(length(N)-1)]/errors.l2$deSolve[2:length(N)])
q$femR = log2(errors.l2$femR[1:(length(N)-1)]/errors.l2$femR[2:length(N)])

cat("deSolve order = ", q$deSolve, "\n")
cat("femR order = ", q$femR, "\n")


imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
    dir.create(imgdir_)
    
pdf(paste(imgdir_,"diff_rates_order_1_deSolve.pdf",sep=""))
par(mfrow=c(1,1))
plot(log2(h), log2(errors.l2$femR), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        ylim = c(min(log2(h^2), log2(errors.l2$femR)), max(log2(h), log2(errors.l2$femR))+2),
        xlab = TeX("$h$"), ylab="", cex.lab=1.25, main=TeX("$\\| u - u_{ex} \\|_{2}$"))
lines(log2(h), log2(errors.l2$deSolve), col = "blue", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
lines(log2(h), log2(h^2), col = "black", type = "b", pch = 16, lwd = 3, lty =2, cex = 2 )
legend("topleft", legend=c("femR", "deSolve", TeX("$h^2$")), 
        col=c("red", "blue", "black"), 
        lty = 2, 
        cex=1.25)
dev.off()

pdf(paste(imgdir_,"diff_times_order_1_deSolve.pdf",sep=""))
plot(log2(h), log2(times$femR), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        xlab = TeX("$h$"), ylab="", cex.lab=1.25, main=TeX("$time [s]$"))
lines(log2(h), log2(times$deSolve), col = "blue", type = "b", pch = 16, lwd = 3, lty =2, cex = 2)
legend("topright", legend=c("femR", "deSolve"), 
        col=c("red", "blue"), 
        lty = 2, 
        cex=1.25)
dev.off()

png(paste(imgdir_,"diff_times_order_1_deSolve.png",sep=""))
plot(log2(h), log2(times$femR), col="red", type="b", pch =16, lwd = 3, lty = 2, cex = 2,
        xlab = TeX("$h$"), ylab="", cex.lab=1.25, main=TeX("$time [s]$"))
lines(log2(h), log2(times$deSolve), col = "blue", type = "b", pch = 16, lwd = 3, lty =2, cex = 2)
legend("topright", legend=c("femR", "deSolve"), 
        col=c("red", "blue"), 
        lty = 2, 
        cex=1.25)
dev.off()
