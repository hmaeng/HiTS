library(xtable)
library(hdbinseg)
library(InspectChangepoint)


###############################################################################
###############################################################################
############################### R functions ###################################
###############################################################################
###############################################################################

#########################################################
### decomposition
hitguh.decomp <- function(x, p = .01) { # x is a matrix
  
  d <- dim(x)[1]
  n <- dim(x)[2]
  noe <- n-1
  
  edges <- matrix(0, noe, 2)
  
  edges[,1] <- 1:(n-1)
  edges[,2] <- 2:n
  
  decomp.hist <- array(0, dim=c(4*d, 2, n-1))
  
  tguh.coeffs <- x
  vec.weights <- rep(1, n)
  
  steps.left <- n-1
  current.step <- 0
  
  while (dim(edges)[1]) {
    
    max.current.steps <- ceiling(p * steps.left)
    
    removable.nodes <- rep(1, max(edges))
    
    a <- vec.weights[edges[,1]]
    b <- vec.weights[edges[,2]]
    
    h1 <- 1/sqrt(1 + (a/b)^2)
    h2 <- -1/sqrt(1 + (b/a)^2)
    l1 <- -h2
    l2 <- h1
    
    Dmat <- sweep(tguh.coeffs[, edges[,1], drop=F], MARGIN=2, h1,`*`) + sweep(tguh.coeffs[, edges[,2], drop=F], MARGIN=2, h2,`*`)
    
    details <- apply(abs(Dmat), 2, max) # L infinity
    
    ord.det <- order(abs(details))
    
    edge.indices.2b.removed <- 1
    traverse.edges.index <- 1
    removable.nodes[edges[ord.det[1],1]] <- removable.nodes[edges[ord.det[1],2]] <- 0
    
    while  ( (length(edge.indices.2b.removed) < max.current.steps) & (traverse.edges.index < noe) ) {
      traverse.edges.index <- traverse.edges.index + 1
      if (removable.nodes[edges[ord.det[traverse.edges.index],1]] & removable.nodes[edges[ord.det[traverse.edges.index],2]]) {
        edge.indices.2b.removed <- c(edge.indices.2b.removed, traverse.edges.index)
        removable.nodes[edges[ord.det[traverse.edges.index],1]] <- removable.nodes[edges[ord.det[traverse.edges.index],2]] <- 0
      }
    }
    
    details.min.ind <- ord.det[edge.indices.2b.removed]
    
    no.of.current.steps <- length(edge.indices.2b.removed)
    
    for(j in 1:d){
      decomp.hist[(j-1)*4+1,,(current.step+1):(current.step+no.of.current.steps)] <- t(edges[details.min.ind,])
      decomp.hist[(j-1)*4+2,1,(current.step+1):(current.step+no.of.current.steps)] <- h1[details.min.ind]
      decomp.hist[(j-1)*4+2,2,(current.step+1):(current.step+no.of.current.steps)] <- h2[details.min.ind]
      
      decomp.hist[(j-1)*4+3,1,(current.step+1):(current.step+no.of.current.steps)] <- Dmat[j, details.min.ind]
      
      smooth.at.min <- l1[details.min.ind] * tguh.coeffs[j, edges[details.min.ind,1]] +
        l2[details.min.ind] * tguh.coeffs[j, edges[details.min.ind,2]]
      decomp.hist[(j-1)*4+3,2,(current.step+1):(current.step+no.of.current.steps)] <- smooth.at.min
      
      decomp.hist[(j)*4,1,(current.step+1):(current.step+no.of.current.steps)] <- vec.weights[edges[details.min.ind,1]]^2
      decomp.hist[(j)*4,2,(current.step+1):(current.step+no.of.current.steps)] <- vec.weights[edges[details.min.ind,2]]^2
    }
    
    eating.up <- apply(matrix(edges[details.min.ind,], no.of.current.steps, 2), 1, min)
    eaten.up <- apply(matrix(edges[details.min.ind,], no.of.current.steps, 2), 1, max)
    
    tguh.coeffs[,eating.up] <- sweep(tguh.coeffs[,edges[details.min.ind,1], drop=F], MARGIN=2, l1[details.min.ind],`*`) + sweep(tguh.coeffs[,edges[details.min.ind,2], drop=F], MARGIN=2, l2[details.min.ind],`*`)
    
    tguh.coeffs[,eaten.up] <- Dmat[, details.min.ind]
    
    det.weight.at.min <- h1[details.min.ind] * vec.weights[edges[details.min.ind,1]] +
      h2[details.min.ind] * vec.weights[edges[details.min.ind,2]]
    sm.weight.at.min <- l1[details.min.ind] * vec.weights[edges[details.min.ind,1]] +
      l2[details.min.ind] * vec.weights[edges[details.min.ind,2]]
    
    vec.weights[eating.up] <- sm.weight.at.min
    vec.weights[eaten.up] <- det.weight.at.min
    
    edges <- plyr::mapvalues(edges, eaten.up, eating.up)
    
    edges <- edges[edges[,1] != edges[,2],]
    
    if (length(edges) == 2) edges <- matrix(edges, 1, 2)
    
    noe <- dim(edges)[1]
    steps.left <- steps.left - no.of.current.steps
    current.step <- current.step + no.of.current.steps
    
  }
  
  list(n = n, d = d, decomp.hist=decomp.hist, tguh.coeffs=tguh.coeffs)
  
}

#########################################################
### Thresholding
hitguh.denoise <- function(tguh.obj, lambda, minseglen = 1, bal = 1/20) {
  
  n <- tguh.obj$n
  d <- tguh.obj$d
  
  if(d>1){
    detail.all <- tguh.obj$decomp.hist[c(4*(1:d)-1),1,]
    details <- apply(abs(detail.all), 2, max) # columnwise maximum of Detail matrix
  } else {
    details <- abs(tguh.obj$decomp.hist[3,1,])
  }
  
  protected <- rep(0, n)
  
  for (i in 1:(n-1)) {
    
    if (!protected[tguh.obj$decomp.hist[1,1,i]] & !protected[tguh.obj$decomp.hist[1,2,i]])
      tguh.obj$decomp.hist[c(4*(1:d)-1),1,i] <- tguh.obj$decomp.hist[c(4*(1:d)-1),1,i] * (
        details[i] > lambda &
          (round(tguh.obj$decomp.hist[4,1,i]) >= minseglen) &
          (round(tguh.obj$decomp.hist[4,2,i]) >= minseglen) &
          (tguh.obj$decomp.hist[2,1,i]^2 > bal) &
          (tguh.obj$decomp.hist[2,2,i]^2 > bal)
      )
    
    if (abs(tguh.obj$decomp.hist[3,1,i]) > 0) protected[tguh.obj$decomp.hist[1,1,i]] <- 1
    
  }
  
  tguh.obj
  
}

#########################################################
### inverse transformation 
hitguh.reconstr <- function(tguh.obj) {
  
  n <- tguh.obj$n
  d <- tguh.obj$d
  
  for (i in (n-1):1) {
    
    ind <- tguh.obj$decomp.hist[1,,i]
    
    tguh.obj$decomp.hist[c(4*(1:d)-1),2,i] <- tguh.obj$tguh.coeffs[,min(ind)]
    
    tguh.obj$tguh.coeffs[,ind[2]] <- tguh.obj$decomp.hist[2,2,i]*tguh.obj$decomp.hist[c(4*(1:d)-1),1,i] + tguh.obj$decomp.hist[2,1,i]*tguh.obj$decomp.hist[c(4*(1:d)-1),2,i]
    
    tguh.obj$tguh.coeffs[,ind[1]] <- tguh.obj$decomp.hist[2,1,i]*tguh.obj$decomp.hist[c(4*(1:d)-1),1,i] - tguh.obj$decomp.hist[2,2,i]*tguh.obj$decomp.hist[c(4*(1:d)-1),2,i]
    
    
  }
  
  tguh.obj
  
}

#########################################################
### finding change-point index
find.cpt.ind <- function(tguh.obj, lambda){
  
  pp1fit <- tguh.obj$tguh.coeffs
  chp <- which(abs(diff(tguh.obj$tguh.coeffs[1, ])) > sqrt(.Machine$double.eps))
  
  if(length(chp) >0){
    chp <- c(1, chp+1, tguh.obj$n+1)
    pqr <- cbind(chp[1:(length(chp)-2)], chp[2:(length(chp)-1)]-1, chp[3:length(chp)]-1)
    a <- sqrt((pqr[,3]- pqr[,2]) / (pqr[,3]- pqr[,1] + 1))
    b <- sqrt((pqr[,2]- pqr[,1] + 1) / (pqr[,3]- pqr[,1] + 1))
    details <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(pqr)[1])
    
    for(j in 1:dim(pqr)[1]){
      details[, j] <- abs(a[j]*apply(pp1fit[, (pqr[j,1]:pqr[j,2]), drop=F], 1, sum)/sqrt(pqr[j,2]-pqr[j,1]+1) - b[j]*apply(pp1fit[, ((pqr[j,2]+1):pqr[j,3]), drop=F], 1, sum)/sqrt(pqr[j,3]-pqr[j,2]))
    }
    
    ### final change points
    tguh.obj$chp <- chp[-c(1, length(chp))]-1
    tguh.obj$details <- details
    tguh.obj$tguh.coeffs <- pp1fit
  } else{
    details <- matrix(0, nrow=dim(pp1fit)[1], ncol=2)
    
    tguh.obj$chp <- chp
    tguh.obj$details <- details
    tguh.obj$tguh.coeffs <- pp1fit
  }
  
  
  return(tguh.obj)
}

#########################################################
### all three steps (from inverse transform)
hitguh.cpt <- function(x, th.const = 1, p = .01, minseglen = 1, bal = 0) {
  
  n <- dim(x)[2]
  
  if (n == 1) {
    
    est <- x
    no.of.cpt <- 0
    cpt <- integer(0)
    
  } else {
    
    ### normalisation and letting sigma=1
    varest <- apply(x, 1, function(x){(mad(diff(x))/sqrt(2))}) # scale each time series with their MAD
    x <- x/varest
    lambda <- th.const * sqrt(2 * log(n*dim(x)[1])) # let sigma=1
    
    dcmp <- hitguh.decomp(x, p=p)
    
    dns <- hitguh.denoise(dcmp, lambda=lambda, minseglen, bal)
    
    inv <- hitguh.reconstr(dns)
    inv.end <- find.cpt.ind(inv, lambda=lambda)
    
    # same definition of cp with that of DC, INSPECT (last point of the constant subregion)
    cpt <- which(abs(diff(inv.end$tguh.coeffs[1, ])) > sqrt(.Machine$double.eps))
    no.of.cpt <- length(cpt)
    cptind <- ifelse(inv.end$details > lambda, 1, 0)
    
    if(exists("varest")){
      est <- inv$tguh.coeffs*varest # rescale (same with the above line)
    } else{
      est <- inv$tguh.coeffs
    }
  }
  
  list(est=est, no.of.cpt=no.of.cpt, cpt=cpt, cptind=cptind)
  
}




#########################################################
### misc functions
cpts.into.intervals <- function(cpt, n) {
  
  len.cpt <- length(cpt)
  
  if (len.cpt) cpt <- sort(cpt)
  iv <- matrix(0, 2, len.cpt+1)
  iv[1,1] <- 1
  iv[2,len.cpt+1] <- n
  
  if (len.cpt) {
    iv[1, 2:(len.cpt+1)] <- cpt+1
    iv[2, 1:len.cpt] <- cpt
  }
  
  iv
  
}

mean.from.cpts <- function(x, cpt) { # x is a matrix
  
  iv <- cpts.into.intervals(cpt, dim(x)[2])
  
  len.cpt <- length(cpt)
  
  means <- matrix(NA, nrow=dim(x)[1], ncol=len.cpt+1)
  for(j in 1:dim(x)[1]){
    for (i in 1:(len.cpt+1)){
      means[j, i] <- mean(x[j, iv[1,i]:iv[2,i]])
    }
  }
  
  estx <- matrix(NA, nrow=dim(x)[1], ncol=dim(x)[2])
  for(j in 1:dim(x)[1]){
    estx[j, ] <- rep(means[j,], iv[2,]-iv[1,]+1)
  }
  
  estx
}

assess <- function(object, modelnum, models, truex){
  qdiff <- length(which(object$cpts>0)) - length(models[[modelnum]]$cpt)
  mse <- sum((truex-object$fit)^2)/(dim(object$fit)[2]) # divided by (n) in high-dim case
  dh <- finding.dH(chp=object$cpts, modelnum=modelnum, models=models)
  return(list(qdiff=qdiff, mse=mse, dh=dh))
}

have.signal <- function(model, sparsity, d, ovl){ #ovl=c("complete", "half", "no")
  
  if(model$cpt.type == "pcwsConstMean"){
    
    signal <- rep(0, model$n)
    
    if(ovl=="complete" | (ovl=="half" & length(model$cpt)==1) | (ovl=="no" & length(model$cpt)==1) | d*sparsity==1){
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      signal[segments[1,1]:segments[1,2]] <- model$start[1]
      
      for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      for(j in 1:d1){
        x[j, ] <- signal}
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
    } else if(ovl=="half" & length(model$cpt) > 1 & d*sparsity>1){
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      ### set 1
      segments <- cbind(c(1, model$cpt[1:round(length(model$cpt)/2)]+1), c(model$cpt[1:round(length(model$cpt)/2)],model$n))
      signal[segments[1,1]:segments[1,2]] <- model$start[1]
      for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
      
      for(j in 1:round(d1/2)){x[j, ] <- signal}
      
      ### set 2
      segments <- cbind(c(1, model$cpt+1), c(model$cpt,model$n))
      signal[segments[1,1]:segments[1,2]] <- model$start[1]
      for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
      
      for(j in (round(d1/2)+1):d1){x[j, ] <- signal}
      
      ### set 3
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
      
    } else if(ovl=="no" & length(model$cpt) > 1 & d*sparsity>1){
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      ### set 1
      segments <- cbind(c(1, model$cpt[1:round(length(model$cpt)/2)]+1), c(model$cpt[1:round(length(model$cpt)/2)],model$n))
      signal[segments[1,1]:segments[1,2]] <- model$start[1]
      for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
      
      for(j in 1:round(d1/2)){x[j, ] <- signal}
      
      ### set 2
      segments <- cbind(c(1, model$cpt[(round(length(model$cpt)/2)+1):length(model$cpt)]+1), c(model$cpt[(round(length(model$cpt)/2)+1):length(model$cpt)], model$n))
      signal[segments[1,1]:segments[1,2]] <- sum(c(model$start[1], model$jump.size[1:round(length(model$cpt)/2)]))
      for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[round(length(model$cpt)/2)+j-1]
      
      for(j in (round(d1/2)+1):d1){x[j, ] <- signal}
      
      ### set 3
      segments <- cbind(c(1, model$cpt+1), c(model$cpt,model$n))
      signal[segments[1,1]:segments[1,2]] <- model$start[1]
      for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
      
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
      
    }
    
    
  } else if(model$cpt.type == "pcwsLinContMean"){
    
    signal <- rep(0, model$n)
    
    if(ovl=="complete"){
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[j-1]
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      for(j in 1:d1){x[j, ] <- signal}
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
    } else if(ovl=="half"){
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      ### set 1
      segments <- cbind(c(1,model$cpt[1:round(length(model$cpt)/2)]+1), c(model$cpt[1:round(length(model$cpt)/2)],model$n))
      
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      
      for(j in 2:max(2, nrow(segments)-1)) {
        slope <- slope +  model$jump.size[j-1]
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      for(k in segments[nrow(segments), 1]:segments[nrow(segments), 2]){
        signal[k] <- signal[k-1] + 0
      }
      
      
      for(j in 1:round(d1/2)){x[j, ] <- signal}
      
      ### set 2
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[j-1]
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      for(j in (round(d1/2)+1):d1){x[j, ] <- signal}
      
      ### set 3
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
    } else if(ovl=="no") {
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      ### set 1
      segments <- cbind(c(1,model$cpt[1:round(length(model$cpt)/2)]+1), c(model$cpt[1:round(length(model$cpt)/2)],model$n))
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      
      for(j in 2:max(2, nrow(segments)-1)) {
        slope <- slope +  model$jump.size[j-1]
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      for(k in segments[nrow(segments), 1]:segments[nrow(segments), 2]){
        signal[k] <- signal[k-1] + 0
      }
      
      for(j in 1:round(d1/2)){x[j, ] <- signal}
      
      ### set 2
      segments <- cbind(c(1, model$cpt[(round(length(model$cpt)/2)+1):length(model$cpt)]+1), c(model$cpt[(round(length(model$cpt)/2)+1):length(model$cpt)], model$n))
      slope <- sum(c(model$start[2], model$jump.size[1:round(length(model$cpt)/2)]))
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * 0
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[round(length(model$cpt)/2)+j-1]
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      for(j in (round(d1/2)+1):d1){x[j, ] <- signal}
      
      ### set 3
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[j-1]
        for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
      }
      
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
    }
    
    
  } else if(model$cpt.type == "pcwsLinMean"){
    
    signal <- rep(0, model$n)
    
    if(ovl=="complete"){
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * slope
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[j-1,2]
        signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[j-1,1]
        if(segments[j,1]+1 < segments[j,2]){
          for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
        }
      }
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      for(j in 1:d1){
        x[j, ] <- signal}
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
    } else if(ovl=="half"){
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      ### set 1
      segments <- cbind(c(1,model$cpt[1:round(length(model$cpt)/2)]+1), c(model$cpt[1:round(length(model$cpt)/2)],model$n))
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * slope
      
      for(j in 2:max(2, nrow(segments)-1)) {
        slope <- slope +  model$jump.size[j-1,2]
        signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[j-1,1]
        if(segments[j,1]+1 < segments[j,2]){
          for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
        }
      }
      for(k in (segments[nrow(segments),1]+1):segments[nrow(segments),2]){
        signal[k] <- signal[k-1] + 0
      }
      
      for(j in 1:round(d1/2)){
        x[j, ] <- signal}
      
      ### set 2
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * slope
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[j-1,2]
        signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[j-1,1]
        if(segments[j,1]+1 < segments[j,2]){
          for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
        }
      }
      
      for(j in (round(d1/2)+1):d1){
        x[j, ] <- signal}
      
      ### set 3
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
      
    } else if(ovl=="no") {
      
      ### data generating
      x <- matrix(0, nrow=d, ncol=model$n)
      d1 <- round(d*sparsity)
      
      ### set 1
      segments <- cbind(c(1,model$cpt[1:round(length(model$cpt)/2)]+1), c(model$cpt[1:round(length(model$cpt)/2)],model$n))
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * slope
      
      for(j in 2:max(2, nrow(segments)-1)) {
        slope <- slope +  model$jump.size[j-1,2]
        signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[j-1,1]
        if(segments[j,1]+1 < segments[j,2]){
          for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
        }
      }
      for(k in (segments[nrow(segments),1]+1):segments[nrow(segments),2]){
        signal[k] <- signal[k-1] + 0
      }
      
      for(j in 1:round(d1/2)){
        x[j, ] <- signal}
      
      ### set 2
      segments <- cbind(c(1, model$cpt[(round(length(model$cpt)/2)+1):length(model$cpt)]+1), c(model$cpt[(round(length(model$cpt)/2)+1):length(model$cpt)], model$n))
      slope <- sum(c(model$start[2], model$jump.size[1:round(length(model$cpt)/2),2]))
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * 0
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[round(length(model$cpt)/2)+j-1,2]
        signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[round(length(model$cpt)/2)+j-1,1]
        if(segments[j,1]+1 < segments[j,2]){
          for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
        }
      }
      
      for(j in (round(d1/2)+1):d1){x[j, ] <- signal}
      
      ### set 3
      segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
      slope <- model$start[2]
      signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * slope
      
      for(j in 2:nrow(segments)) {
        slope <- slope +  model$jump.size[j-1,2]
        signal[segments[j,1]] <-  signal[segments[j-1,2]] + model$jump.size[j-1,1]
        if(segments[j,1]+1 < segments[j,2]){
          for(k in (segments[j,1]+1):segments[j,2]) signal[k] <- signal[k-1] + slope
        }
      }
      if(d1<d){
        for(j in (d1+1):d){
          set.seed(j)
          x[j, ] <- rep(mean(signal), model$n)}
      }
    }
    
  }
  
  return(x)
}

finding.dH <- function(chp, modelnum, models){
  n <- models[[modelnum]]$n
  segments.endpoints.est <- sort(unique(c(0, chp, n)))
  segments.endpoints.true <- sort(unique(c(0, models[[modelnum]]$cpt, n)))
  
  distm <- abs(matrix(rep(segments.endpoints.est, length(segments.endpoints.true)), nrow=length(segments.endpoints.est))
    -matrix(rep(segments.endpoints.true, length(segments.endpoints.est)), nrow=length(segments.endpoints.est), byrow=T))
  
  screening.dist <- max(apply(distm, 2, min)) * 100 / n
  precision.dist <- max(apply(distm, 1, min)) * 100 / n
  dH <- mean((abs(screening.dist-precision.dist) + screening.dist + precision.dist)/2)
  
  return(dH)
}

###############################################################################
###############################################################################
################################## Methods ####################################
###############################################################################
###############################################################################

hits <- function(x, thr, bal, p){
  
  tic <- proc.time()
  object <- hitguh.cpt(x, th.const = thr, bal = bal, p = p, minseglen = 1)
  toc <- proc.time()
  
  list(fit = object$est, cpts=object$cpt, elapsed=(toc-tic)[3], cptind=object$cptind)
  
}

sbs <- function(x){
  
  tic <- proc.time()
  object <- sbs.alg(x, cp.type=1, thr=NULL, temporal=FALSE, do.parallel=4)
  cpts <- object$ecp
  fit <- mean.from.cpts(x, cpt=cpts)
  toc <- proc.time()
  
  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}

dc <- function(x){
  
  tic <- proc.time()
  object <- dcbs.alg(x, cp.type=1, phi=-1, thr=NULL, temporal=FALSE, do.parallel=4)
  cpts <- object$ecp
  fit <- mean.from.cpts(x, cpt=cpts)
  toc <- proc.time()
  
  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}

is <- function(x){
  
  tic <- proc.time()
  object <- inspect(x, threshold = compute.threshold(n=dim(x)[2], p=dim(x)[1]))
  cpts <- object$changepoints[,1]
  fit <- mean.from.cpts(x, cpt=cpts)
  toc <- proc.time()
  
  list(fit = fit, cpts=cpts, elapsed=(toc-tic)[3])
}

###############################################################################
###############################################################################
################################### Models ####################################
###############################################################################
###############################################################################

model.three <-  list(name = "three",
  cpt.type = "pcwsConstMean",
  jump.size = c(1,-1.5,2),
  #jump.size = c(1,-1.5,2),
  #jump.size = c(1.2,-1.8,2.4),
  start = 0,
  cpt = c(50,100,150),
  n = 200)

model.two <-  list(name = "two",
  cpt.type = "pcwsConstMean",
  #jump.size = 0.7*c(2,-3),
  jump.size = c(2,-2),
  start = -1,
  cpt = c(33,66),
  n = 100)

model.twosmall <-  list(name = "two",
  cpt.type = "pcwsConstMean",
  #jump.size = c(2,-3)/2,
  jump.size = c(2,-2)/1.5,
  start = -0.5,
  cpt = c(33,66),
  n = 100)

model.teeth <-  list(name = "teeth",
  cpt.type = "pcwsConstMean",
  jump.size = 2*(-1)^{1:9},
  start = 1,
  cpt = (1:9) * 30,
  n = 30 * 10)

model.extremeteeth <-  list(name = "teeth",
  cpt.type = "pcwsConstMean",
  jump.size = 3*(-1)^{1:39},
  start = 1.5,
  cpt = ceiling((1:39) * 12.5),
  n = 12.5 * 40)

model.blocks <-  list(name = "blocks",
  cpt.type = "pcwsConstMean",
  cpt = ceiling(c(205, 267, 308, 472, 512, 820, 902, 1332, 1557, 1598, 1659)/2),
  jump.size = diff(c(0, 14.64, -4.66, 8.32, -7.32, 10.98, -5.39, 6.29, 19.03, 3.68, 19.37, 0))/ 10,
  start = 0,
  n = 1000)


models <- list(model.two, model.twosmall, model.three, model.teeth, model.extremeteeth, model.blocks)



##### one-dimensional signal only 
par(mfrow=c(3,2),mar=c(4,3,2,2))
for(i in 1:length(models)){
  truex <- have.signal(models[[i]], sparsity=0.05, d=100, ovl="complete")
  x <- truex + rnorm(length(truex), sd=1)
  plot(x[1,], type="l", ylab="", xlab="t", ylim=range(x), col=8, lwd=2)
  lines(truex[1,], col=1, lwd=2)
}


##### data matrix image
modelnum=1
sgma=1 
sparsity=0.1 
d=100 
ovl="complete"

truex <- have.signal(model=models[[modelnum]], sparsity=sparsity, d=d, ovl=ovl) 
x <- truex + matrix(rnorm(dim(truex)[1]*dim(truex)[2]), nrow=dim(truex)[1])

par(mfrow=c(1,1),mar=c(2,2,2,2))
image(t(x), axes=F)


###############################################################################
###############################################################################
################################ SIMULATIONS ##################################
###############################################################################
###############################################################################


simall <- function(N, modelnum, sgma, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=10, ovl="complete"){ # ovl=c("complete", "half", "no")
  
  truex <- have.signal(model=models[[modelnum]], sparsity=sparsity, d=d, ovl=ovl) # whole matrix
  n <- dim(truex)[2]
  
  par(mfrow=c(3,1),mar=c(2,2,1.5,1.5))
  
  ### for assesement
  result <- list(qdiff=matrix(NA, nrow=N, ncol=4), mse=matrix(NA, nrow=N, ncol=4), dh=matrix(NA, nrow=N, ncol=4), 
    time=matrix(NA, nrow=N, ncol=4), ri=matrix(NA, nrow=N, ncol=4))
  result <- lapply(result, function(x) {colnames(x) <- c("HITS.pwc", "SBS", "DC", "INSPECT"); x})
  
  for(K in 1:N){
    
    
    set.seed(K)
    x <- truex + matrix(rnorm(d*n), nrow=d)
    
    #####################################################################
    ############################## hits ################################
    #####################################################################
    obj1 <- hits(x, thr=thr, bal=bal, p=p) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj1, modelnum=modelnum, models=models, truex=truex)
    
    g1 <- ifelse(apply(truex, 1, function(x){sum(abs(diff(x)) > sqrt(.Machine$double.eps))})!=0, 1, 0)
    g2 <- ifelse(rowSums(obj1$cptind)!=0, 1, 0)
    ri <- fossil::rand.index(g1, g2)
    
    result[[1]][K,1] <- a.obj$qdiff
    result[[2]][K,1] <- a.obj$mse
    result[[3]][K,1] <- a.obj$dh
    result[[4]][K,1] <- obj1$elapsed
    result[[5]][K,1] <- ri
    
    #####################################################################
    ################################# SBS ###############################
    #####################################################################
    obj2 <- sbs(x) ## cpts is turned "NA" if nothing is detected
    a.obj <- assess(object=obj2, modelnum=modelnum, models=models, truex=truex)
    result[[1]][K,2] <- a.obj$qdiff
    result[[2]][K,2] <- a.obj$mse
    result[[3]][K,2] <- a.obj$dh
    result[[4]][K,2] <- obj2$elapsed
    
    #####################################################################
    ################################# DC ################################
    #####################################################################
    obj3 <- dc(x) ## cpts is turned "0" if nothing is detected
    a.obj <- assess(object=obj3, modelnum=modelnum, models=models, truex=truex)
    result[[1]][K,3] <- a.obj$qdiff
    result[[2]][K,3] <- a.obj$mse
    result[[3]][K,3] <- a.obj$dh
    result[[4]][K,3] <- obj3$elapsed
    
    #####################################################################
    ############################## inspect ##############################
    #####################################################################
    obj4 <- is(x) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj4, modelnum=modelnum, models=models, truex=truex)
    result[[1]][K,4] <- a.obj$qdiff
    result[[2]][K,4] <- a.obj$mse
    result[[3]][K,4] <- a.obj$dh
    result[[4]][K,4] <- obj4$elapsed
    
    
    ### plot 1
    plot(x[1,], type="l", col=8, ylim=range(x), xlab="")
    lines(truex[1, ], col=7, lwd=2, lty=1)
    lines(obj1$fit[1, ], col=1, lwd=1, lty=2)
    lines(obj2$fit[1, ], col=2, lwd=1, lty=3)
    lines(obj3$fit[1, ], col=3, lwd=1, lty=4)
    lines(obj4$fit[1, ], col=4, lwd=1, lty=5)
    
    legend("topright", c("obs","true","HITS.pwc", "SBS", "DC", "INSPECT"), col=c(8,7,1:4), lty=c(1,1,2:5), lwd=c(1,2,rep(1,4)), cex=1, bty="n", ncol=2)
    
    ### plot 2
    plot(x[sparsity*d,], type="l", col=8, ylim=range(x), xlab="")
    lines(truex[sparsity*d, ], col=7, lwd=2, lty=1)
    lines(obj1$fit[sparsity*d, ], col=1, lwd=1, lty=2)
    lines(obj2$fit[sparsity*d, ], col=2, lwd=1, lty=3)
    lines(obj3$fit[sparsity*d, ], col=3, lwd=1, lty=4)
    lines(obj4$fit[sparsity*d, ], col=4, lwd=1, lty=5)
    
    legend("topright", c("obs","true","HITS.pwc", "SBS", "DC", "INSPECT"), col=c(8,7,1:4), lty=c(1,1,2:5), lwd=c(1,2,rep(1,4)), cex=1, bty="n", ncol=2)
    
    
    ### plot 3
    plot(x[d,], type="l", col=8, ylim=range(x), xlab="")
    lines(truex[d, ], col=7, lwd=2, lty=1)
    lines(obj1$fit[d, ], col=1, lwd=1, lty=2)
    lines(obj2$fit[d, ], col=2, lwd=1, lty=3)
    lines(obj3$fit[d, ], col=3, lwd=1, lty=4)
    lines(obj4$fit[d, ], col=4, lwd=1, lty=5)
    
    legend("topright", c("obs","true","HITS.pwc", "SBS", "DC", "INSPECT"), col=c(8,7,1:4), lty=c(1,1,2:5), lwd=c(1,2,rep(1,4)), cex=1, bty="n", ncol=2)
    
    print(K)
  }
  
  return(result)
  
}


###################################################################
######################## ovl="complete" ###########################
###################################################################


###################################################################
### N=100 / d=100
m1N100s1d100com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m1N100s1d100com, file="m1N100s1d100com.RData")
m1N100s2d100com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m1N100s2d100com, file="m1N100s2d100com.RData")
m1N100s3d100com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m1N100s3d100com, file="m1N100s3d100com.RData")

m2N100s1d100com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m2N100s1d100com, file="m2N100s1d100com.RData")
m2N100s2d100com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m2N100s2d100com, file="m2N100s2d100com.RData")
m2N100s3d100com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m2N100s3d100com, file="m2N100s3d100com.RData")

m3N100s1d100com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m3N100s1d100com, file="m3N100s1d100com.RData")
m3N100s2d100com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m3N100s2d100com, file="m3N100s2d100com.RData")
m3N100s3d100com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m3N100s3d100com, file="m3N100s3d100com.RData")

m4N100s1d100com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m4N100s1d100com, file="m4N100s1d100com.RData")
m4N100s2d100com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m4N100s2d100com, file="m4N100s2d100com.RData")
m4N100s3d100com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m4N100s3d100com, file="m4N100s3d100com.RData")

m5N100s1d100com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m5N100s1d100com, file="m5N100s1d100com.RData")
m5N100s2d100com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m5N100s2d100com, file="m5N100s2d100com.RData")
m5N100s3d100com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m5N100s3d100com, file="m5N100s3d100com.RData")

m6N100s1d100com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m6N100s1d100com, file="m6N100s1d100com.RData")
m6N100s2d100com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m6N100s2d100com, file="m6N100s2d100com.RData")
m6N100s3d100com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m6N100s3d100com, file="m6N100s3d100com.RData")

###################################################################
### N=100 / d=300
m1N100s1d300com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m1N100s1d300com, file="m1N100s1d300com.RData")
m1N100s2d300com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m1N100s2d300com, file="m1N100s2d300com.RData")
m1N100s3d300com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m1N100s3d300com, file="m1N100s3d300com.RData")

m2N100s1d300com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m2N100s1d300com, file="m2N100s1d300com.RData")
m2N100s2d300com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m2N100s2d300com, file="m2N100s2d300com.RData")
m2N100s3d300com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m2N100s3d300com, file="m2N100s3d300com.RData")

m3N100s1d300com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m3N100s1d300com, file="m3N100s1d300com.RData")
m3N100s2d300com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m3N100s2d300com, file="m3N100s2d300com.RData")
m3N100s3d300com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m3N100s3d300com, file="m3N100s3d300com.RData")

m4N100s1d300com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m4N100s1d300com, file="m4N100s1d300com.RData")
m4N100s2d300com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m4N100s2d300com, file="m4N100s2d300com.RData")
m4N100s3d300com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m4N100s3d300com, file="m4N100s3d300com.RData")

m5N100s1d300com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m5N100s1d300com, file="m5N100s1d300com.RData")
m5N100s2d300com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m5N100s2d300com, file="m5N100s2d300com.RData")
m5N100s3d300com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m5N100s3d300com, file="m5N100s3d300com.RData")

m6N100s1d300com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m6N100s1d300com, file="m6N100s1d300com.RData")
m6N100s2d300com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m6N100s2d300com, file="m6N100s2d300com.RData")
m6N100s3d300com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m6N100s3d300com, file="m6N100s3d300com.RData")

###################################################################
### N=100 / d=500
m1N100s1d500com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m1N100s1d500com, file="m1N100s1d500com.RData")
m1N100s2d500com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m1N100s2d500com, file="m1N100s2d500com.RData")
m1N100s3d500com <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m1N100s3d500com, file="m1N100s3d500com.RData")

m2N100s1d500com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m2N100s1d500com, file="m2N100s1d500com.RData")
m2N100s2d500com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m2N100s2d500com, file="m2N100s2d500com.RData")
m2N100s3d500com <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m2N100s3d500com, file="m2N100s3d500com.RData")

m3N100s1d500com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m3N100s1d500com, file="m3N100s1d500com.RData")
m3N100s2d500com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m3N100s2d500com, file="m3N100s2d500com.RData")
m3N100s3d500com <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m3N100s3d500com, file="m3N100s3d500com.RData")

m4N100s1d500com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m4N100s1d500com, file="m4N100s1d500com.RData")
m4N100s2d500com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m4N100s2d500com, file="m4N100s2d500com.RData")
m4N100s3d500com <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m4N100s3d500com, file="m4N100s3d500com.RData")

m5N100s1d500com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m5N100s1d500com, file="m5N100s1d500com.RData")
m5N100s2d500com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m5N100s2d500com, file="m5N100s2d500com.RData")
m5N100s3d500com <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m5N100s3d500com, file="m5N100s3d500com.RData")

m6N100s1d500com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m6N100s1d500com, file="m6N100s1d500com.RData")
m6N100s2d500com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m6N100s2d500com, file="m6N100s2d500com.RData")
m6N100s3d500com <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m6N100s3d500com, file="m6N100s3d500com.RData")





###################################################################
########################### ovl="half" ############################
###################################################################
## In the case of d=100/sparsity=0.01 with ovl="half" or ovl="no",
## only one coordinate has true change-points.

###################################################################
### N=100 / d=100
m1N100s1d100half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m1N100s1d100half, file="m1N100s1d100half.RData")
m1N100s2d100half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m1N100s2d100half, file="m1N100s2d100half.RData")
m1N100s3d100half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m1N100s3d100half, file="m1N100s3d100half.RData")

m2N100s1d100half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m2N100s1d100half, file="m2N100s1d100half.RData")
m2N100s2d100half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m2N100s2d100half, file="m2N100s2d100half.RData")
m2N100s3d100half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m2N100s3d100half, file="m2N100s3d100half.RData")

m3N100s1d100half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m3N100s1d100half, file="m3N100s1d100half.RData")
m3N100s2d100half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m3N100s2d100half, file="m3N100s2d100half.RData")
m3N100s3d100half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m3N100s3d100half, file="m3N100s3d100half.RData")

m4N100s1d100half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m4N100s1d100half, file="m4N100s1d100half.RData")
m4N100s2d100half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m4N100s2d100half, file="m4N100s2d100half.RData")
m4N100s3d100half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m4N100s3d100half, file="m4N100s3d100half.RData")

m5N100s1d100half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m5N100s1d100half, file="m5N100s1d100half.RData")
m5N100s2d100half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m5N100s2d100half, file="m5N100s2d100half.RData")
m5N100s3d100half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m5N100s3d100half, file="m5N100s3d100half.RData")

m6N100s1d100half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m6N100s1d100half, file="m6N100s1d100half.RData")
m6N100s2d100half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m6N100s2d100half, file="m6N100s2d100half.RData")
m6N100s3d100half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m6N100s3d100half, file="m6N100s3d100half.RData")

###################################################################
### N=100 / d=300
m1N100s1d300half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m1N100s1d300half, file="m1N100s1d300half.RData")
m1N100s2d300half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m1N100s2d300half, file="m1N100s2d300half.RData")
m1N100s3d300half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m1N100s3d300half, file="m1N100s3d300half.RData")

m2N100s1d300half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m2N100s1d300half, file="m2N100s1d300half.RData")
m2N100s2d300half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m2N100s2d300half, file="m2N100s2d300half.RData")
m2N100s3d300half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m2N100s3d300half, file="m2N100s3d300half.RData")

m3N100s1d300half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m3N100s1d300half, file="m3N100s1d300half.RData")
m3N100s2d300half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m3N100s2d300half, file="m3N100s2d300half.RData")
m3N100s3d300half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m3N100s3d300half, file="m3N100s3d300half.RData")

m4N100s1d300half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m4N100s1d300half, file="m4N100s1d300half.RData")
m4N100s2d300half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m4N100s2d300half, file="m4N100s2d300half.RData")
m4N100s3d300half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m4N100s3d300half, file="m4N100s3d300half.RData")

m5N100s1d300half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m5N100s1d300half, file="m5N100s1d300half.RData")
m5N100s2d300half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m5N100s2d300half, file="m5N100s2d300half.RData")
m5N100s3d300half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m5N100s3d300half, file="m5N100s3d300half.RData")

m6N100s1d300half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m6N100s1d300half, file="m6N100s1d300half.RData")
m6N100s2d300half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m6N100s2d300half, file="m6N100s2d300half.RData")
m6N100s3d300half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m6N100s3d300half, file="m6N100s3d300half.RData")

###################################################################
### N=100 / d=500
m1N100s1d500half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m1N100s1d500half, file="m1N100s1d500half.RData")
m1N100s2d500half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m1N100s2d500half, file="m1N100s2d500half.RData")
m1N100s3d500half <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m1N100s3d500half, file="m1N100s3d500half.RData")

m2N100s1d500half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m2N100s1d500half, file="m2N100s1d500half.RData")
m2N100s2d500half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m2N100s2d500half, file="m2N100s2d500half.RData")
m2N100s3d500half <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m2N100s3d500half, file="m2N100s3d500half.RData")

m3N100s1d500half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m3N100s1d500half, file="m3N100s1d500half.RData")
m3N100s2d500half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m3N100s2d500half, file="m3N100s2d500half.RData")
m3N100s3d500half <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m3N100s3d500half, file="m3N100s3d500half.RData")

m4N100s1d500half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m4N100s1d500half, file="m4N100s1d500half.RData")
m4N100s2d500half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m4N100s2d500half, file="m4N100s2d500half.RData")
m4N100s3d500half <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m4N100s3d500half, file="m4N100s3d500half.RData")

m5N100s1d500half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m5N100s1d500half, file="m5N100s1d500half.RData")
m5N100s2d500half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m5N100s2d500half, file="m5N100s2d500half.RData")
m5N100s3d500half <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m5N100s3d500half, file="m5N100s3d500half.RData")

m6N100s1d500half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m6N100s1d500half, file="m6N100s1d500half.RData")
m6N100s2d500half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m6N100s2d500half, file="m6N100s2d500half.RData")
m6N100s3d500half <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m6N100s3d500half, file="m6N100s3d500half.RData")



###################################################################
########################### ovl="no" ##############################
###################################################################

###################################################################
### N=100 / d=100
m1N100s1d100no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m1N100s1d100no, file="m1N100s1d100no.RData")
m1N100s2d100no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m1N100s2d100no, file="m1N100s2d100no.RData")
m1N100s3d100no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m1N100s3d100no, file="m1N100s3d100no.RData")

m2N100s1d100no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m2N100s1d100no, file="m2N100s1d100no.RData")
m2N100s2d100no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m2N100s2d100no, file="m2N100s2d100no.RData")
m2N100s3d100no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m2N100s3d100no, file="m2N100s3d100no.RData")

m3N100s1d100no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m3N100s1d100no, file="m3N100s1d100no.RData")
m3N100s2d100no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m3N100s2d100no, file="m3N100s2d100no.RData")
m3N100s3d100no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m3N100s3d100no, file="m3N100s3d100no.RData")

m4N100s1d100no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m4N100s1d100no, file="m4N100s1d100no.RData")
m4N100s2d100no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m4N100s2d100no, file="m4N100s2d100no.RData")
m4N100s3d100no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m4N100s3d100no, file="m4N100s3d100no.RData")

m5N100s1d100no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m5N100s1d100no, file="m5N100s1d100no.RData")
m5N100s2d100no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m5N100s2d100no, file="m5N100s2d100no.RData")
m5N100s3d100no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m5N100s3d100no, file="m5N100s3d100no.RData")

m6N100s1d100no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m6N100s1d100no, file="m6N100s1d100no.RData")
m6N100s2d100no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m6N100s2d100no, file="m6N100s2d100no.RData")
m6N100s3d100no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m6N100s3d100no, file="m6N100s3d100no.RData")

###################################################################
### N=100 / d=300
m1N100s1d300no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m1N100s1d300no, file="m1N100s1d300no.RData")
m1N100s2d300no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m1N100s2d300no, file="m1N100s2d300no.RData")
m1N100s3d300no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m1N100s3d300no, file="m1N100s3d300no.RData")

m2N100s1d300no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m2N100s1d300no, file="m2N100s1d300no.RData")
m2N100s2d300no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m2N100s2d300no, file="m2N100s2d300no.RData")
m2N100s3d300no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m2N100s3d300no, file="m2N100s3d300no.RData")

m3N100s1d300no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m3N100s1d300no, file="m3N100s1d300no.RData")
m3N100s2d300no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m3N100s2d300no, file="m3N100s2d300no.RData")
m3N100s3d300no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m3N100s3d300no, file="m3N100s3d300no.RData")

m4N100s1d300no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m4N100s1d300no, file="m4N100s1d300no.RData")
m4N100s2d300no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m4N100s2d300no, file="m4N100s2d300no.RData")
m4N100s3d300no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m4N100s3d300no, file="m4N100s3d300no.RData")

m5N100s1d300no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m5N100s1d300no, file="m5N100s1d300no.RData")
m5N100s2d300no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m5N100s2d300no, file="m5N100s2d300no.RData")
m5N100s3d300no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m5N100s3d300no, file="m5N100s3d300no.RData")

m6N100s1d300no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m6N100s1d300no, file="m6N100s1d300no.RData")
m6N100s2d300no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m6N100s2d300no, file="m6N100s2d300no.RData")
m6N100s3d300no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m6N100s3d300no, file="m6N100s3d300no.RData")

###################################################################
### N=100 / d=500
m1N100s1d500no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m1N100s1d500no, file="m1N100s1d500no.RData")
m1N100s2d500no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m1N100s2d500no, file="m1N100s2d500no.RData")
m1N100s3d500no <- simall(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m1N100s3d500no, file="m1N100s3d500no.RData")

m2N100s1d500no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m2N100s1d500no, file="m2N100s1d500no.RData")
m2N100s2d500no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m2N100s2d500no, file="m2N100s2d500no.RData")
m2N100s3d500no <- simall(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m2N100s3d500no, file="m2N100s3d500no.RData")

m3N100s1d500no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m3N100s1d500no, file="m3N100s1d500no.RData")
m3N100s2d500no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m3N100s2d500no, file="m3N100s2d500no.RData")
m3N100s3d500no <- simall(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m3N100s3d500no, file="m3N100s3d500no.RData")

m4N100s1d500no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m4N100s1d500no, file="m4N100s1d500no.RData")
m4N100s2d500no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m4N100s2d500no, file="m4N100s2d500no.RData")
m4N100s3d500no <- simall(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m4N100s3d500no, file="m4N100s3d500no.RData")

m5N100s1d500no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m5N100s1d500no, file="m5N100s1d500no.RData")
m5N100s2d500no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m5N100s2d500no, file="m5N100s2d500no.RData")
m5N100s3d500no <- simall(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m5N100s3d500no, file="m5N100s3d500no.RData")

m6N100s1d500no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m6N100s1d500no, file="m6N100s1d500no.RData")
m6N100s2d500no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m6N100s2d500no, file="m6N100s2d500no.RData")
m6N100s3d500no <- simall(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m6N100s3d500no, file="m6N100s3d500no.RData")





###############################################################################
###############################################################################
############################# DATA APPLICATION ################################
###############################################################################
###############################################################################


############################################################################
######################### TEMPERATURE DATA #################################
############################################################################
####### http://berkeleyearth.org/data/
####### https://www.kaggle.com/berkeleyearth/climate-change-earth-surface-temperature-data
############################################################################

library(tidyr)

MyData <- read.csv(file="GlobalLandTemperaturesByCity.csv", header=TRUE, sep=",")
head(MyData)
ctry <- levels(MyData$Country)
ctry

##########################################################
################### NUMBER OF CITIES #####################
##########################################################
no.city <- rep(0, length(ctry))
for(i in 1:length(ctry)){
  no.city[i] <- length(unique(MyData[MyData$Country==ctry[i],]$City))
}
cbind(c(ctry), no.city)

##########################################################
####################### COUNTRY ##########################
##########################################################
ctrs <- c(ctry)[which(no.city>30)]
ctrs
i <- 19 # SOUTH AFRICA / i=22 is UK 
CR <- MyData[MyData$Country==ctrs[i],]
unique(CR$City)
d <- length( unique(CR$City)) # Number of city
d
ctname <- as.character(unique(CR$City))

##########################################################
################ REARRANGE FOR EACH CITY #################
##########################################################
DAT <- cbind(separate(data.frame(date = CR[,1]), "date", c("Year", "Month", "Day"), sep = "-"), CR)
CT <- unique(CR$City)
head(DAT)

##### EXTRACT January only
DAT.1 <- DAT[DAT$Month=="01", 1:7, drop=F]
head(DAT.1)

##### get rid of the years which contain NA
DAT.1 <- DAT.1[!is.na(DAT.1$AverageTemperature),]
head(DAT.1)
sort(as.numeric(unique(DAT.1$Year)))[diff(sort(as.numeric(unique(DAT.1$Year))))!=1]
range(as.numeric(unique(DAT.1$Year)))

DAT.1 <- DAT.1[DAT.1$Year!="1745", 1:7, drop=F]
yr <- sort(as.numeric(unique(DAT.1$Year)))
yr

##### each column corresponds to each city
allcity <- unique(CR$City)
data <- matrix(0, nrow=length(yr), ncol=d+1)
data[,1] <- yr
colnames(data) <- c("year", as.character(allcity))
for(k in 1:d){
  data[,k+1] <- DAT.1[DAT.1$City==allcity[k], 5]
}
head(data)


###############################################
########## change-point detection #############
###############################################
X <- t(data[,-1,drop=F])
d <- dim(X)[1]
n <- dim(X)[2]
d;n

plot(data[,1], X[1,], type="l", ylim=range(X), ylab="temperature", xlab="Year")
for(k in 2:d){
  lines(data[,1], X[k,], type="l", col=k, lty=1)
}

### HITS_pwc
hitspwc <- hits(x=X, thr=1.2, p=0.04, bal=0)
cind <- hitspwc$cptind
hitspwc$cpts
yr[hitspwc$cpts]

# south africa
g10 <- which(hitspwc$cptind[,1]==1 & hitspwc$cptind[,2]==0)
g01 <- which(hitspwc$cptind[,1]==0 & hitspwc$cptind[,2]==1)
g11 <- which(hitspwc$cptind[,1]==1 & hitspwc$cptind[,2]==1)
g00 <- which(hitspwc$cptind[,1]==0 & hitspwc$cptind[,2]==0)
g10; g01; g11; g00
ctname[g10]
ctname[g01] 
ctname[g11]
ctname[g00]

### SBS 
sbsfit <- sbs(x=X)
yr[sbsfit$cpts]

### DC 
dcfit <- dc(x=X)
yr[dcfit$cpts]

### INSPECT 
isfit <- is(x=X)
yr[isfit$cpts]



### Cape Town
par(mfrow=c(1,1),mar=rep(2.5, 4))
smp <- 10 # 10=Cape town
for(k in smp){
  plot(yr, X[k,], type="b", col=8, ylab="", xlab="", lwd=2)
  
  lines(yr, hitspwc$fit[k,], type="l", col=1, lwd=2, lty=2)
  lines(yr, sbsfit$fit[k,], type="l", col=2, lwd=2, lty=3) 
  lines(yr, dcfit$fit[k,], type="l", col=4, lwd=1, lty=4) 
  lines(yr, isfit$fit[k,], type="l", col=6, lwd=1, lty=5)
  legend("topleft", c("obs","HiTS", "SBS", "DC", "IS"), col=c(8,1,2,4,6), lty=c(1,2:5), lwd=c(2,2,2,1,1), pch=c(1, rep(NA, 4)), cex=1.5, bty="n", ncol=3)
}


