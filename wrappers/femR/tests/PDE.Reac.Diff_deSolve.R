library(ReacTran)
library(latex2exp)
library(femR)
library(fdaPDE)

rm(list=ls())
graphics.off()

rmse <- function(x,y){return(sqrt(mean( (x-y)^2)))}

#2D
exact <- function(x,y){
  return( exp(x + y) )
}

f <- function(x,y){
  return(0. * x^0) 
}

Dy    <- Dx <- -1.   # diffusion coeff, X- and Y-direction

PDE_parameters <- list("diffusion" = 1., "transport" = rbind(0.,0.), "reaction" = 2.)

N = 2*c(20,40,80,160) 

errors.l2 <- list("deSolve" = rep(0, times = length(N)),
                  "femR_1" = rep(0, times = length(N)),
                  "femR_2" = rep(0, times = length(N)))

times <- list("deSolve" = rep(0, times = length(N)),
              "femR_1" = rep(0, times = length(N)),
              "femR_2" = rep(0, times = length(N)))
h <- rep(0, times = length(N))
nnodes <- rep(0, times = length(N))

imgdir_ = "imgs/"
if(!dir.exists(imgdir_))
  dir.create(imgdir_)

domain_ = "2D/"
imgdir_ = paste(imgdir_,domain_,sep="")

if(!dir.exists(imgdir_))
  dir.create(imgdir_)    

imgdir_ = paste(imgdir_,"ReacDiff/",sep="")
if(!dir.exists(imgdir_))
  dir.create(imgdir_)    


CEX.axis = 1.25
CEX = 3
LWD = 3

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
  
  C.x.up    <- exact(rep(x.grid$x.up, times=N[i]), y.grid$x.mid)
  C.y.up    <- exact(x.grid$x.mid, rep(y.grid$x.up, times=N[i]))
  
  C.x.down  <- exact(rep(x.grid$x.down, times=N[i]), y.grid$x.mid)
  C.y.down  <- exact(x.grid$x.mid, rep(y.grid$x.down, times=N[i]))
  
  #forcing = matrix(nrow=N[i], ncol=N[i], data=0)
  y.ex = matrix(nrow=N[i], ncol=N[i], data=0)
  for(k in 1:N[i]){
    for(l in 1:N[i]){
      #forcing[k,l] = f(grid2D$x.mid[k], grid2D$y.mid[l])
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
    
    dY <- dY + 2*Y 
    
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
  error.deSolve = abs(y.hat-y.ex)
  # if(i == 1){
  #   png(paste(imgdir_, "error_deSolve.png", sep = ""))
  #   filled.contour(x = x.grid$x.mid,
  #                  y = y.grid$x.mid,
  #                  abs(y.hat-y.ex), main="deSolve",cex.axis = CEX.axis, cex.main=CEX.axis)
  #   dev.off()
  #   
  #   png(paste(imgdir_, "ReacDiff_deSolve.png", sep = ""))
  #   filled.contour(x = x.grid$x.mid,
  #                  y = x.grid$x.mid,
  #                  y.hat, main="deSolve",cex.axis = CEX.axis, cex.main=CEX.axis)
  #   dev.off()
  # }
  # ###### femR - ORDER 1 ######
  grid2D <- expand.grid(x.grid$x.int, y.grid$x.int)
  mesh <- fdaPDE::create.mesh.2D(grid2D)
  nnodes[i] <- nrow(mesh$nodes)
  square <- list(nodes= mesh$nodes, edges= mesh$segments, elements= mesh$triangles, neigh= mesh$neighbors, boundary= mesh$nodesmarkers)
  
  PDE <- new(PDE_2D_isotropic_ORDER_1, square)
  PDE$set_PDEparameters(PDE_parameters)
  
  dirichletBC <- as.matrix(exact(square$nodes[,1], square$nodes[,2]))
  PDE$set_dirichletBC(dirichletBC)
  
  quadrature_nodes <- PDE$get_quadrature_nodes()
  
  PDE$set_forcingTerm( matrix(0, nrow=nrow(quadrature_nodes), ncol=1))
  start_ <- Sys.time()
  result <- PDE$solve()
  times$femR_1[i] <- as.numeric(difftime(Sys.time(), start_, units="secs"))
  cat("femR order 1", times$femR_1[i] , " s\n")  
  u_ex <- as.matrix(exact(mesh$nodes[,1], mesh$nodes[,2]))
  
  errors.l2$femR_1[i] <- rmse(u_ex, result$solution)
  
  if(i == 1){
    
    FEMbasis = create.FEM.basis(mesh)
    u_h = FEM(result$solution, FEMbasis)
    x = x.grid$x.int
    y = x.grid$x.int
    U_h = matrix(0, nrow=(N[i]+1),ncol=(N[i]+1))
    y.ex = matrix(0,nrow=(N[i]+1), ncol=(N[i]+1))
    for(k in 1:(N[i]+1)){
      for(l in 1:(N[i]+1)){
        U_h[k,l] = eval.FEM(u_h, locations= cbind(x[k],y[l]))
        y.ex[k,l] =  exact(x[k], y[l]) 
      }
    }
    
    error.femR = abs(U_h-y.ex)
    min.err = min(error.femR, error.deSolve)
    max.err = max(error.femR, error.deSolve)
    
    min.sol = min(min(y.ex, U_h), y.hat)
    max.sol = max(max(y.ex, U_h), y.hat)
    png(paste(imgdir_, "error_femR_1.png", sep = ""))
    filled.contour(error.femR, main=TeX("$|u - u_h|$"), zlim = c(min.err, max.err),
                   cex.axis = CEX.axis, cex.main=CEX.axis)
    dev.off()
    
    png(paste(imgdir_, "ReacDiff_femR_1.png", sep = ""))
    filled.contour(U_h, main=TeX("$u_h$"), zlim = c(min.sol, max.sol),
                   cex.axis = CEX.axis, cex.main=CEX.axis)
    dev.off()
    
    png(paste(imgdir_, "ReacDiff.png", sep = ""))
    filled.contour(y.ex, main="", zlim = c(min.sol, max.sol),
                   cex.axis = CEX.axis, cex.main=CEX.axis)
    dev.off()
    
    png(paste(imgdir_, "error_deSolve.png", sep = ""))
    filled.contour(x = x.grid$x.mid,
                   y = y.grid$x.mid, 
                   error.deSolve, main=TeX("$|u - u_h|$"), zlim = c(min.err, max.err),
                   cex.axis = CEX.axis, cex.main=CEX.axis)
    dev.off()
    
    png(paste(imgdir_, "ReacDiff_deSolve.png", sep = ""))
    filled.contour(x = x.grid$x.mid,
                   y = x.grid$x.mid,
                   y.hat, main=TeX("$u_h$"), zlim = c(min.sol, max.sol),
                   cex.axis = CEX.axis, cex.main=CEX.axis)
    dev.off()
    
    # png(paste(imgdir_, "error_femR_1.png", sep = ""))
    # filled.contour(abs(U_h-y.ex), main="femR order 1",cex.axis = CEX.axis, cex.main=CEX.axis)
    # dev.off()
    # 
    # png(paste(imgdir_, "ReacDiff_femR_1.png", sep = ""))
    # filled.contour(U_h, main="femR order 1",cex.axis = CEX.axis, cex.main=CEX.axis)
    # dev.off()
  }
  
  ###### femR - ORDER 2 ######
  
  PDE <- new(PDE_2D_isotropic_ORDER_2, square)
  PDE$set_PDEparameters(PDE_parameters)
  
  nodes <- PDE$get_dofs_coordinates()  
  dirichletBC <- as.matrix(exact(nodes[,1], nodes[,2]))
  PDE$set_dirichletBC(dirichletBC)
  
  quadrature_nodes <- PDE$get_quadrature_nodes()
  
  PDE$set_forcingTerm(as.matrix(f(quadrature_nodes[,1], quadrature_nodes[,2])))
  start_ <- Sys.time()
  result <- PDE$solve()
  times$femR_2[i] <- as.numeric(difftime(Sys.time(), start_, units="secs"))
  cat("femR order 2", times$femR_2[i] , " s\n")  
  u_ex <- as.matrix(exact(nodes[,1], nodes[,2]))
  
  errors.l2$femR_2[i] <- rmse(u_ex, result$solution)
}

q = list("deSolve" = 0, "femR_1"= 0, "femR_2"= 0)
q$deSolve = log2(errors.l2$deSolve[1:(length(N)-1)]/errors.l2$deSolve[2:length(N)])
q$femR_1 = log2(errors.l2$femR_1[1:(length(N)-1)]/errors.l2$femR_1[2:length(N)])
q$femR_2 = log2(errors.l2$femR_2[1:(length(N)-1)]/errors.l2$femR_2[2:length(N)])

cat("deSolve order = ", q$deSolve, "\n")
cat("femR 1 order = ", q$femR_1, "\n")
cat("femR 2 order = ", q$femR_2, "\n")

p = list("deSolve" = 0, "femR_1"= 0, "femR_2"= 0)
p$deSolve = log10(errors.l2$deSolve[1:(length(N)-1)]/errors.l2$deSolve[2:length(N)])/log10(nnodes[2:length(N)]/nnodes[1:(length(N)-1)])
p$femR_1 = log10(errors.l2$femR_1[1:(length(N)-1)]/errors.l2$femR_1[2:length(N)])/log10(nnodes[2:length(N)]/nnodes[1:(length(N)-1)])
p$femR_2 = log10(errors.l2$femR_2[1:(length(N)-1)]/errors.l2$femR_2[2:length(N)])/log10(nnodes[2:length(N)]/nnodes[1:(length(N)-1)])
cat("deSolve order (nodes)  = ", p$deSolve, "\n")
cat("femR 1 order (nodes) = ", p$femR_1, "\n")
cat("femR 2 order (nodes) = ", p$femR_2, "\n")

pdf(paste(imgdir_,"ReacDiff_rates_deSolve.pdf",sep=""))
par(mai=c(1,1,0.5,0.5))
plot(log2(h), log2(errors.l2$femR_1), col="red", type="b", pch =16, lwd = LWD, lty = 1, cex = CEX,
     ylim = c( min( min(log2(h^3), log2(errors.l2$femR_1)), min(log2(errors.l2$deSolve), log2(errors.l2$femR_2)) ), 
               max( max(log2(h^2), log2(errors.l2$femR_1)), max(log2(errors.l2$deSolve), log2(errors.l2$femR_2))) ),
     xlab = TeX("$log_2(h)$"), ylab = TeX("$log_2(\\| u - u_{ex} \\|_{2})$"),
     cex.lab=CEX.axis, cex.axis=CEX.axis, cex.main=CEX.axis)
lines(log2(h), log2(errors.l2$femR_2), col = "orange", type = "b", pch = 16, lwd = LWD, lty =1, cex = CEX)
lines(log2(h), log2(errors.l2$deSolve), col = "blue", type = "b", pch = 16, lwd = LWD, lty =1, cex = CEX)
lines(log2(h), log2(h^2), col = "black", lwd = 3, lty =2)
lines(log2(h), log2(h^3)-2, col = "black", lwd = 3, lty =2)
text(log2(h[3])+0.65,(log2(h[3]^2)+2.), TeX("$h^2$"),cex=CEX.axis)
text(log2(h[3]),(log2(h[3]^3)-2+0.75), TeX("$h^3$"),cex=CEX.axis)
legend("topleft", legend=c("deSolve", "femR order 1", "femR order 2"), 
       col=c("blue", "red", "orange"), 
       lty = 1, 
       lwd = LWD+2, 
       cex=CEX.axis)

dev.off()

pdf(paste(imgdir_,"ReacDiff_times_deSolve.pdf",sep=""))
par(mai=c(1,1,0.5,0.5))
plot(log10(nnodes), log10(times$femR_1), col="red", type="b", pch =16, lwd = LWD, lty = 1, cex = CEX,
     ylim = c(min( min(log10(times$femR_1),log10(times$deSolve)), log10(times$femR_2)), 
              max( max(log10(times$femR_1), log10(times$deSolve)), log10(times$femR_2))),
     xlab = TeX("$log_{10}(nodes)$"), ylab=TeX("$log_{10}(time)$"), 
     cex.lab=CEX.axis, cex.axis=CEX.axis, cex.main=CEX.axis)
lines(log10(nnodes), log10(times$deSolve), col = "blue", type = "b", pch = 16, lwd = LWD, lty =1, cex = CEX)
lines(log10(nnodes), log10(times$femR_2), col = "orange", type = "b", pch = 16, lwd = LWD, lty =1, cex = CEX)
legend("topleft", legend=c("deSolve", "femR order 1", "femR order 2"), 
       col=c("blue", "red", "orange"), 
       lty = 1, 
       lwd = LWD+2, 
       cex=CEX.axis)
dev.off()