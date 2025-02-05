# Covariance matrix
# ---------------------------------------------------------------
covariance = function(coords, sigma2, gamma, phi){
  distM = crossdist(coords)/phi
  Rm    = distM*besselK(distM, 1)
  diag(Rm) = 1
  cov   = gamma*sigma2*Rm + (1-gamma)*sigma2*diag(nrow(coords))
  return((cov + t(cov))/2)
}

# Random number generation
# ---------------------------------------------------------------




# Matrices used in the SPDE approximation
# ---------------------------------------------------------------
matrix.Q = function(coord, iedge, plot=FALSE) {
  coord = as.matrix(coord)
  dmat  = crossdist(coord)
  mesh  = INLA::inla.mesh.2d(loc.domain=coord,
                             max.edge=c(iedge*max(range(dmat)), max(range(dmat))))

  if (plot) {
    plot(mesh)
    points(coord, pch = 20, col = "seagreen", cex = .75)
  }
  mat = INLA::inla.mesh.fem(mesh, order = 2)
  A   = INLA::inla.spde.make.A(mesh, loc = coord)
  D  = mat[["c0"]]
  G1 = mat[["g1"]]
  G2 = mat[["g2"]]
  return(list(D=D, G1=G1, G2=G2, A=A, mesh=mesh))
}

# Log-likelihood function
# -----------------------------------------------------------------------------
log_likelihood = function(y, cc, lower, upper, x, beta, sigma2, Psi){
  mean = x%*%beta
  Variance = sigma2*Psi

  if (sum(cc)==0) {
    logver = dmvnorm(c(y), c(mean), Variance, TRUE, TRUE)
  } else {
    invObs = solve(Variance[cc==0,cc==0])
    meanC  = mean[cc==1] + Variance[cc==1,cc==0]%*%invObs%*%(y[cc==0]-mean[cc==0])
    meanC  = as.vector(meanC)
    varC   = Variance[cc==1,cc==1] - Variance[cc==1,cc==0]%*%invObs%*%Variance[cc==0,cc==1]
    varC   = 0.5*(varC + t(varC))

    logden2 = log(pmvnorm(lower[cc==1], upper[cc==1], meanC, NULL, varC)[1])
    logver = dmvnorm(c(y[cc==0]), c(mean[cc==0]), as.matrix(Variance[cc==0,cc==0]), TRUE, TRUE) + logden2
  }
  results = list(loglik=logver, AIC=-2*logver+2*(length(beta)+3), BIC=-2*logver+(length(beta)+3)*log(length(y)))
  return (results)
}

# Estimation
# ---------------------------------------------------------------
SPDE_Spatial = function(y, x, ci, lcl, ucl, coords, phi0, gamma0, lower, upper,
                      iedge, MaxIter, M, pc, error, show_se){

  if (sum(ci)==0) {
    t1 = Sys.time()
    #output = Spatial_model(y, x, coords, phi0, gamma0, lower, upper, MaxIter, error, show_se)
    t2 = Sys.time()
    ptime = t2 - t1

  } else {
    t1 = Sys.time()

    mats = matrix.Q(coords, iedge=iedge, plot=FALSE)
    output = saemEsp(y, x, ci, coords, phi0, gamma0, lcl, ucl, lower, upper, M, MaxIter,
                     pc, error, mats, show_se)

    t2 = Sys.time()
    ptime = t2 - t1
  }
  likeli = log_likelihood(y, ci, lcl, ucl, x, output$beta, output$sigma2, output$Psi)
  if (show_se){
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, gamma=output$gamma, EY=output$EY, EYY=output$EYY, SE=output$SE,
                  InfMat=output$InfMat, loglik=likeli$loglik, AIC=likeli$AIC, BIC=likeli$BIC,
                  Iter=output$Iterations, time=ptime, X=x, coord=coords, nodes=ncol(mats$A), mats=mats)
  } else {
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, gamma=output$gamma, EY=output$EY, EYY=output$EYY, loglik=likeli$loglik,
                  AIC=likeli$AIC, BIC=likeli$BIC, Iter=output$Iterations, time=ptime, X=x,
                  coord=coords, nodes=ncol(mats$A), mats=mats)
  }
  return (output)
}


# Convergence plot
# -----------------------------------------------------------------------------
plotconvergence = function(model){
  Theta = model$Theta
  X = model$X
  q = length(model$beta)
  Iter = model$Iter

  myplot = vector("list",q+3)
  if (all(X[,1]==1)){ namesE = c(seq(0,(q-1)),0,0,0) } else { namesE = c(seq(1,q),0,0,0) }
  listabeta = rbind(namesE,Theta)
  listabeta = as.list(data.frame(listabeta[,1:q]))
  myplot[1:q] = lapply(listabeta, function(.x) ggplot(data.frame(.x[-1]),aes(x=seq(0,Iter),y=.x[-1])) +
                         geom_line() + labs(x="Iteration", y=bquote(beta[.(.x[1])])) + theme_bw())

  myplot[[q+1]] = ggplot(data.frame(Theta),aes(x=seq(0,Iter),y=Theta[,q+1])) + geom_line() + labs(x="Iteration", y=bquote(sigma^2)) + theme_bw()
  myplot[[q+2]] = ggplot(data.frame(Theta),aes(x=seq(0,Iter),y=Theta[,q+2])) + geom_line() + labs(x="Iteration", y=bquote(phi)) + theme_bw()
  myplot[[q+3]] = ggplot(data.frame(Theta),aes(x=seq(0,Iter),y=Theta[,q+3])) + geom_line() + labs(x="Iteration", y=bquote(gamma)) + theme_bw()

  grid.arrange(grobs=myplot, ncol=3)
}
