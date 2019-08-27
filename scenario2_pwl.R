
library(tidyr)
library(dplyr)

###############################################################################
###############################################################################
############################### R functions ###################################
###############################################################################
###############################################################################

#########################################################
### decomposition
hd.bts.dcmp <- function(x, p = .01) { 
  
  d <- dim(x)[1]
  n <- dim(x)[2]
  noe <- n-2 
  
  weights.const <- rep(1, n)
  weights.lin <- 1:n
  idx <- 1:n
  paired <- c()
  
  edges <- matrix(0, noe, 3) 
  
  edges[,1] <- 1:(n-2)
  edges[,2] <- 2:(n-1)
  edges[,3] <- 3:n
  
  edges <- cbind(edges, 0)
  
  decomp.hist <- array(0, dim=c(4*d, 3, n-2))
  
  bts.coeffs <- x
  
  steps.left <- n-2 
  current.step <- 0 
  
  sameboat <- c()
  
  
  while (dim(edges)[1]) {
    
    max.current.steps <- ceiling(p * steps.left)
    removable.nodes <- rep(1, max(idx))
    
    Dmat <- matrix(0, nrow=d, ncol=dim(edges)[1])
    
    if(any(edges[, 4]>0)){
      
      pr <- matrix(which(edges[,4]!=0), nrow=2)
      for(j in 1:d){
        cd1 <- computeDET(edges=edges, edgerow=pr[1,], weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,])
        detcoef <- cd1$detcoef
        p1.details <- cd1$det
        
        M0 <- apply(detcoef, 2, orth.matrix)
        upd.wc <- cbind(rowSums(cd1$wc*t(M0[c(2,5,8),])), rowSums(cd1$wc*t(M0[c(3,6,9),])), weights.const[edges[pr[2,],3]])
        upd.wl <- cbind(rowSums(cd1$wl*t(M0[c(2,5,8),])), rowSums(cd1$wl*t(M0[c(3,6,9),])), weights.lin[edges[pr[2,],3]])
        upd.bts.coeffs <- cbind(rowSums(cd1$tc*t(M0[c(2,5,8),])), rowSums(cd1$tc*t(M0[c(3,6,9),])), bts.coeffs[edges[pr[2,],3]])
        p2.details <- colSums(apply(cbind(upd.wc, upd.wl), 1, filter.bts)*t(upd.bts.coeffs))
        
        p.detail <- rep(apply(cbind(abs(p1.details), abs(p2.details)), 1, max), each=2)
        
        if(length(pr)!=dim(edges)[1]){
          
          cd3 <- computeDET(edges=edges, edgerow=c(1:dim(edges)[1])[-c(pr)], weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,])
          np.details <- cd3$det
          Dmat[j, ] <- c(p.detail, np.details)
        } else{
          Dmat[j, ] <- p.detail
        }
      }
      
    } else{
      edgerow <- c(1:dim(edges)[1])
      sub.wc <- cbind(weights.const[edges[edgerow,1]], weights.const[edges[edgerow,2]], weights.const[edges[edgerow,3]])
      sub.wl <- cbind(weights.lin[edges[edgerow,1]], weights.lin[edges[edgerow,2]], weights.lin[edges[edgerow,3]])
      detcoef <- apply(cbind(sub.wc, sub.wl), 1, filter.bts)
      
      Dmat <- sweep(bts.coeffs[, edges[,1], drop=F], MARGIN=2, detcoef[1,],`*`) + sweep(bts.coeffs[, edges[,2], drop=F], MARGIN=2, detcoef[2,],`*`) +
        sweep(bts.coeffs[, edges[,3], drop=F], MARGIN=2, detcoef[3,],`*`)
      
    }
    
    colmaxD <- apply(abs(Dmat), 2, max) # L infinity
    
    ord.det <- order(colmaxD)
    
    cand <- rbind(ord.det, edges[ord.det,4])
    
    eitr <- 1 
    tei <- 1 
    if(cand[2,1] > 0){
      removable.nodes[edges[ord.det[1:2],1]] <- removable.nodes[edges[ord.det[1:2],2]] <- removable.nodes[edges[ord.det[1:2],3]] <- 0
      tei <- tei + 1
      eitr <- c(eitr, tei)
    } else{
      removable.nodes[edges[ord.det[1],1]] <- removable.nodes[edges[ord.det[1],2]] <- removable.nodes[edges[ord.det[1],3]] <- 0
    }
    
    while  ( ( length(eitr) < max.current.steps ) & ( tei < noe ) ) {
      
      tei <- tei + 1
      
      if(cand[2, tei] > 0){
        if(sum(removable.nodes[edges[ord.det[tei:(tei+1)],1]])==2 & sum(removable.nodes[edges[ord.det[tei:(tei+1)],2]])==2 & sum(removable.nodes[edges[ord.det[tei:(tei+1)],3]])==2){
          removable.nodes[edges[ord.det[tei:(tei+1)],1]] <- removable.nodes[edges[ord.det[tei:(tei+1)],2]] <- removable.nodes[edges[ord.det[tei:(tei+1)],3]] <- 0
          eitr <- c(eitr, tei, tei+1)
          tei <- tei + 1
          
        }
      } else{
        if(removable.nodes[edges[ord.det[tei],1]] & removable.nodes[edges[ord.det[tei],2]] & removable.nodes[edges[ord.det[tei],3]]){
          eitr <- c(eitr, tei)
          removable.nodes[edges[ord.det[tei],1]] <- removable.nodes[edges[ord.det[tei],2]] <- removable.nodes[edges[ord.det[tei],3]] <- 0
        }
      }
    }
    
    details.min.ind <- ord.det[eitr]
    
    no.of.current.steps <- length(eitr)
    
    ee <- matrix(edges[details.min.ind,], no.of.current.steps, 4)
    sameboat <- c(sameboat, c(ee[,4]))
    idx0 <- idx 
    
    if(sum(ee[,4]>0)==0){
      ee <- ee[, -4, drop=F]
      
      for(j in 1:d){
        udt <- updating(ee=ee, weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,], idx=idx0)
        bts.coeffs[j,] <- udt$bts.coeffs
        
        decomp.hist[(j-1)*4+1,,(current.step+1):(current.step+no.of.current.steps)] <- t(ee)
        decomp.hist[(j-1)*4+2,,(current.step+1):(current.step+no.of.current.steps)] <- udt$h
        decomp.hist[(j-1)*4+3,,(current.step+1):(current.step+no.of.current.steps)] <- t(udt$tc1)
        decomp.hist[(j)*4,,(current.step+1):(current.step+no.of.current.steps)] <- balance.np(paired=paired, ee=ee, idx=idx0, no.of.current.steps=no.of.current.steps, n=n)
      }
      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      idx <- udt$idx
      
    } else{
      
      pr <- matrix(which(ee[,4]!=0), nrow=2)
      ee[pr[2,], 1:2] <- ee[pr[1,], 1:2]
      ee <- ee[, -4, drop=F]
      ee.p1 <- ee[pr[1,],,drop=F]
      ee.p2 <- ee[pr[2,],,drop=F]
      ee.np <- ee[-c(pr),,drop=F]
      
      for(j in 1:d){
        udt <- updating(ee=ee.p1, weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs[j,], idx=idx0)
        bts.coeffs[j,] <- udt$bts.coeffs
        wg.c <- udt$weights.const
        wg.l <- udt$weights.lin
        idx <- udt$idx
        
        decomp.hist[(j-1)*4+1,,(current.step+1):(current.step+dim(ee.p1)[1])] <- t(ee.p1)
        decomp.hist[(j-1)*4+2,,(current.step+1):(current.step+dim(ee.p1)[1])] <- udt$h
        decomp.hist[(j-1)*4+3,,(current.step+1):(current.step+dim(ee.p1)[1])] <- t(udt$tc1)
        decomp.hist[(j-1)*4+4,,(current.step+1):(current.step+dim(ee.p1)[1])] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx0, ee.p2=ee.p2, n=n)
        
        udt <- updating(ee=ee.p2, weights.const=wg.c, weights.lin=wg.l, bts.coeffs=bts.coeffs[j,], idx=idx)
        bts.coeffs[j,] <- udt$bts.coeffs
        wg.c <- udt$weights.const
        wg.l <- udt$weights.lin
        idx <- udt$idx
        
        decomp.hist[(j-1)*4+1,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- t(ee.p2)
        decomp.hist[(j-1)*4+2,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- udt$h
        decomp.hist[(j-1)*4+3,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- t(udt$tc1)
        decomp.hist[(j-1)*4+4,,(current.step+dim(ee.p1)[1]+1):(current.step+length(pr))] <- balance.p(pr=pr, ee.p1=ee.p1, idx=idx, ee.p2=ee.p2, n=n)
        
        if(length(pr)!=dim(ee)[1]){
          
          udt <- updating(ee=ee.np, weights.const=wg.c, weights.lin=wg.l, bts.coeffs=bts.coeffs[j,], idx=idx)
          bts.coeffs[j,] <- udt$bts.coeffs
          
          decomp.hist[(j-1)*4+1,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- t(ee.np)
          decomp.hist[(j-1)*4+2,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- udt$h
          decomp.hist[(j-1)*4+3,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- t(udt$tc1)
          decomp.hist[(j-1)*4+4,,(current.step+length(pr)+1):(current.step+no.of.current.steps)] <- balance.np(paired=paired, ee=ee.np, idx=idx, 
            no.of.current.steps=no.of.current.steps-length(pr), n=n)
        }
      }
      
      weights.const <- udt$weights.const
      weights.lin <- udt$weights.lin
      idx <- udt$idx
      
    }
    
    paired <- sort(unique(c(paired, c(ee[,1:2]))))
    if(length(which(!is.na(match(paired, c(ee[,3])))))>0){
      paired <- sort(paired[-c(which(!is.na(match(paired, c(ee[,3])))))])
    }
    
    edges <- matrix(0, length(idx)-2, 3) 
    edges[,1] <- idx[1:(length(idx)-2)]
    edges[,2] <- idx[2:(length(idx)-1)]
    edges[,3] <- idx[3:(length(idx))]
    
    matchpair <- ifelse(is.na(matrix(match(edges, paired), ncol=3)), NA, 1)
    rs <- which(rowSums(matchpair, na.rm=T)==3)
    if(length(rs)>0){
      edges <- as.matrix(rbind(edges[rs,], edges[-rs,]))
      matchpair <- rbind(matchpair[rs,], matchpair[-rs,])
      edges <- as.matrix(cbind(edges, c(rep(1:(length(rs)/2), each=2), rep(0, dim(edges)[1]-length(rs)))))
    } else{
      edges <- as.matrix(cbind(edges, rep(0, dim(edges)[1])))
    }
    
    removed <- c(which(rowSums(matchpair, na.rm=T)==1), which(matchpair[,1]==1 & is.na(matchpair[,2]) & matchpair[,3]==1))
    if(length(removed)>0){
      edges <- edges[-removed,,drop=F]
    }
    
    noe <- dim(edges)[1]
    steps.left <- steps.left - no.of.current.steps
    current.step <- current.step + no.of.current.steps
    
    if(noe==1 & dim(edges)[1]==1){
      edges[,4] <- 0
    }
    
  }
  
  return(list(n = n, sameboat=sameboat, decomp.hist=decomp.hist, bts.coeffs=bts.coeffs))
  
}

#########################################################
### Thresholding
hd.bts.dns <- function(bts.obj, lambda, bal = 1/20) {
  
  if(dim(bts.obj$decomp.hist)[1]>4){
    d <- dim(bts.obj$decomp.hist)[1]/4
    detail.all <- bts.obj$decomp.hist[c(4*(1:d)-1),1,] 
    details <- apply(abs(detail.all), 2, max) 
  } else {
    d <- 1
    detail.all <- bts.obj$decomp.hist[c(4*(1:d)-1),1,,drop=F]
    details <- abs(bts.obj$decomp.hist[c(4*(1:d)-1),1,]) 
  }
  
  sameboat <- bts.obj$sameboat
  n <- bts.obj$n
  
  protected <- rep(0, n)
  
  for (i in 1:(n-2)) {
    
    if (!protected[bts.obj$decomp.hist[1,1,i]] & !protected[bts.obj$decomp.hist[1,2,i]] &
        !protected[bts.obj$decomp.hist[1,3,i]])
      
      bts.obj$decomp.hist[c(4*(1:d)-1),1,i] <- bts.obj$decomp.hist[c(4*(1:d)-1),1,i] * (
        (details[i] > lambda) &
          (bts.obj$decomp.hist[4,1,i] > bal) &
          (bts.obj$decomp.hist[4,2,i] > bal)
      )
    
    if (abs(bts.obj$decomp.hist[3,1,i]) > 0) protected[bts.obj$decomp.hist[1,1,i]] <- protected[bts.obj$decomp.hist[1,2,i]] <- 1
    
  }
  
  paired <- matrix(which(bts.obj$sameboat!=0), nrow=2)
  
  if(length(paired)>0){
    for(i in 1:ncol(paired)){
      
      overzero <- is.element(paired[,i], which(abs(bts.obj$decomp.hist[3,1,])>0))
      zero <- is.element(paired[,i], which(abs(bts.obj$decomp.hist[3,1,])==0))
      if(sum(overzero)==1 & sum(zero)==1){
        bts.obj$decomp.hist[c(4*(1:d)-1),1,paired[which(overzero==F),i]] <- detail.all[,paired[which(overzero==F),i]]
      }
    }
  }
  
  return(bts.obj)
  
}

#########################################################
### inverse transformation
hd.bts.inv <- function(bts.obj) {
  
  n <- bts.obj$n
  d <- dim(bts.obj$decomp.hist)[1]/4
  
  for (i in (n-2):1) {
    
    M.inv <- t(orth.matrix(bts.obj$decomp.hist[2,,i]))
    
    ind <- bts.obj$decomp.hist[1,,i]
    
    bts.obj$decomp.hist[c(4*(1:d)-1),2,i] <- bts.obj$bts.coeffs[,ind[1]]
    bts.obj$decomp.hist[c(4*(1:d)-1),3,i] <- bts.obj$bts.coeffs[,ind[2]]
    
    tmp <- bts.obj$decomp.hist[c(4*(1:d)-1),,i]
    
    if(d==1){
      rcstr.tmp <- M.inv %*% matrix(tmp, ncol=1)
    } else{
      rcstr.tmp <- M.inv %*% t(tmp)
    }
    
    bts.obj$bts.coeffs[,ind[1]] <- rcstr.tmp[1,]
    bts.obj$bts.coeffs[,ind[2]] <- rcstr.tmp[2,]
    bts.obj$bts.coeffs[,ind[3]] <- rcstr.tmp[3,]
    
  }
  
  return(bts.obj)
  
}

#########################################################
### post processing - stage 1 
hd.bts.pp1 <- function(bts.obj, lambda){ 
  
  wc <- rep(1, bts.obj$n)
  wl <- 1:(bts.obj$n)
  pp1fit <- bts.obj$bts.coeffs
  inicp <- finding.cp(bts.obj)
  
  if(length(inicp) >0){
    chp <- c(1, inicp+1, bts.obj$n+1)
    pqr <- cbind(chp[1:(length(chp)-2)], chp[2:(length(chp)-1)]-1, chp[3:length(chp)]-1)
    d.pqr <- t(apply(pqr, 1, diff)) 
    
    details <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1])
    detail.1 <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1]) 
    detail.2 <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1]) 
    
    while(length(chp)>2){
      ######################################
      ### 1) (x, x)
      c1 <- which(d.pqr[,1]==0 & d.pqr[,2]==1)
      if(length(c1)>0){
        for(i in 1:length(c1)){
          p <- pqr[c1[i],1]
          q <- pqr[c1[i],2]
          r <- pqr[c1[i],3]
          details[,c1[i]] <- abs((pp1fit[,p:q, drop=F]-pp1fit[,(q+1):r, drop=F])/sqrt(2))
        }
      }
      ######################################
      ### 23) (x, x, x)
      c23 <- which((d.pqr[,1]==0 & d.pqr[,2]==2) | (d.pqr[,1]==1 & d.pqr[,2]==1))
      if(length(c23)>0){
        for(i in 1:length(c23)){
          p <- pqr[c23[i],1]
          r <- pqr[c23[i],3]
          h <- filter.bts(c(wc[p:r], wl[p:r]))
          details[,c23[i]] <- abs(pp1fit[,p:r, drop=F]%*%matrix(h, ncol=1))
        }
      }
      ######################################
      ### 4) (xx, xx)
      c4 <- which(d.pqr[,1]==1 & d.pqr[,2]==2)
      if(length(c4)>0){
        for(i in 1:length(c4)){
          p1 <- pqr[c4[i],1]
          q1 <- pqr[c4[i],2]
          r1 <- q1+1
          h1 <- filter.bts(c(wc[p1:r1], wl[p1:r1]))
          detail.1[,c4[i]] <- abs(pp1fit[,p1:r1, drop=F]%*%matrix(h1, ncol=1))
          
          ### updating
          M <- orth.matrix(h1)
          
          r2 <- pqr[c4[i],3]
          h2 <- filter.bts(c((M%*%matrix(wc[p1:r1], ncol=1))[2:3,1], wc[r2], (M%*%matrix(wl[p1:r1], ncol=1))[2:3,1], wl[r2]))
          detail.2[,c4[i]] <- abs(cbind(t(M%*%t(pp1fit[,p1:r1, drop=F]))[,2:3,drop=F], pp1fit[, r2, drop=F])%*%matrix(h2, ncol=1))
          
          details[,c4[i]] <- apply(cbind(abs(detail.1[,c4[i]]), abs(detail.2[,c4[i]])), 1, max)
        }
      }
      ######################################
      ### 56) (x, xx from chunk) or (xx from chunk, x)
      c5 <- which(d.pqr[,1]==0 & d.pqr[,2]>2)
      if(length(c5)>0){
        for(i in 1:length(c5)){
          p <- pqr[c5[i],1]
          q <- pqr[c5[i],2]
          r <- pqr[c5[i],3]
          
          l12 <- L12(r-q)[[r-q-2]]
          wcL <- matrix(wc[(q+1):r], nrow=1)%*%l12
          wlL <- matrix(wl[(q+1):r], nrow=1)%*%l12
          xL <- pp1fit[,(q+1):r, drop=F]%*%l12
          
          h <- filter.bts(c(wc[p:q], wcL , wl[p:q], wlL))
          details[,c5[i]] <- abs(cbind(pp1fit[,p:q, drop=F], xL)%*%matrix(h, ncol=1))
        }
      }
      
      c6 <- which(d.pqr[,1]>1 & d.pqr[,2]==1)
      if(length(c6)>0){
        for(i in 1:length(c6)){
          p <- pqr[c6[i],1]
          q <- pqr[c6[i],2]
          r <- pqr[c6[i],3]
          
          l12 <- L12(q-p+1)[[q-p+1-2]]
          wcL <- matrix(wc[p:q], nrow=1)%*%l12
          wlL <- matrix(wl[p:q], nrow=1)%*%l12
          xL <- pp1fit[,p:q, drop=F]%*%l12
          
          h <- filter.bts(c(wcL, wc[(q+1):r], wlL, wl[(q+1):r]))
          details[,c6[i]] <- abs(cbind(xL, pp1fit[,(q+1):r, drop=F])%*%matrix(h, ncol=1))
        }
      }
      ######################################
      ### 78) (xx, xx from chunk) or (xx from chunk, xx)
      c7 <- which(d.pqr[,1]==1 & d.pqr[,2]>2)
      if(length(c7)>0){
        for(i in 1:length(c7)){
          p <- pqr[c7[i],1]
          q <- pqr[c7[i],2]
          r <- pqr[c7[i],3]
          
          l12 <- L12(r-q)[[r-q-2]]
          wcL <- matrix(wc[(q+1):r], nrow=1)%*%l12
          wlL <- matrix(wl[(q+1):r], nrow=1)%*%l12
          xL <- pp1fit[,(q+1):r, drop=F]%*%l12
          
          new.wc <- c(wc[p:q], wcL)
          new.wl <- c(wl[p:q], wlL)
          new.x <- cbind(pp1fit[,p:q, drop=F], xL)
          
          ### first edge
          h1 <- filter.bts(c(new.wc[1:3], new.wl[1:3]))
          detail.1[,c7[i]] <- abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1))
          
          ### updating
          M <- orth.matrix(h1)
          
          ### second edge
          h2 <- filter.bts(c( (M%*%matrix(new.wc[1:3], ncol=1))[2:3,1], new.wc[4], (M%*%matrix(new.wl[1:3], ncol=1))[2:3,1], new.wl[4]))
          detail.2[,c7[i]] <- abs(cbind(t(M%*%t(new.x[,1:3,drop=F]))[,2:3,drop=F], new.x[,4])%*%matrix(h2, ncol=1))
          
          details[,c7[i]] <-  apply(cbind(abs(detail.1[,c7[i]]), abs(detail.2[,c7[i]])), 1, max)
        }
      }
      
      c8 <- which(d.pqr[,1]>1 & d.pqr[,2]==2)
      if(length(c8)>0){
        for(i in 1:length(c8)){
          p <- pqr[c8[i],1]
          q <- pqr[c8[i],2]
          r <- pqr[c8[i],3]
          
          l12 <- L12(q-p+1)[[q-p+1-2]]
          wcL <- matrix(wc[p:q], nrow=1)%*%l12
          wlL <- matrix(wl[p:q], nrow=1)%*%l12
          xL <- pp1fit[,p:q, drop=F]%*%l12
          
          new.wc <- c(wcL, wc[(q+1):r])
          new.wl <- c(wlL, wl[(q+1):r])
          new.x <- cbind(xL, pp1fit[,(q+1):r, drop=F])
          
          ### first edge
          h1 <- filter.bts(c(new.wc[1:3], new.wl[1:3]))
          detail.1[,c8[i]] <- abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1))
          
          ### updating
          M <- orth.matrix(h1)
          
          ### second edge
          h2 <- filter.bts(c( (M%*%matrix(new.wc[1:3], ncol=1))[2:3,1], new.wc[4], (M%*%matrix(new.wl[1:3], ncol=1))[2:3,1], new.wl[4]))
          detail.2[,c8[i]] <- abs(cbind(t(M%*%t(new.x[,1:3,drop=F]))[,2:3,drop=F], new.x[,4])%*%matrix(h2, ncol=1))
          
          details[,c8[i]] <- apply(cbind(abs(detail.1[,c8[i]]), abs(detail.2[,c8[i]])), 1, max)
        }
      }
      ######################################
      ### 9) the others
      all.idx <- c(1:dim(d.pqr)[1])
      c9 <- all.idx[!all.idx %in% unique(c(c1, c23, c4, c5, c6, c7, c8))]
      
      if(length(c9)>0){
        for(i in 1:length(c9)){
          p <- pqr[c9[i],1]
          q <- pqr[c9[i],2]
          r <- pqr[c9[i],3]
          
          l12 <- L12(q-p+1)[[q-p+1-2]]
          wcL1 <- matrix(wc[p:q], nrow=1)%*%l12
          wlL1 <- matrix(wl[p:q], nrow=1)%*%l12
          xL1 <- pp1fit[,p:q, drop=F]%*%l12
          
          l12 <- L12(r-q)[[r-q-2]]
          wcL2 <- matrix(wc[(q+1):r], nrow=1)%*%l12
          wlL2 <- matrix(wl[(q+1):r], nrow=1)%*%l12
          xL2 <- pp1fit[,(q+1):r, drop=F]%*%l12
          
          new.wc <- c(wcL1, wcL2)
          new.wl <- c(wlL1, wlL2)
          new.x <- cbind(xL1, xL2)
          
          ### first edge
          h1 <- filter.bts(c(new.wc[1:3], new.wl[1:3]))
          detail.1[,c9[i]] <- abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1))
          #detail.1[c9[i]] <- max(abs(new.x[,1:3,drop=F]%*%matrix(h1, ncol=1)))
          
          ### updating
          M <- orth.matrix(h1)
          
          ### second edge
          h2 <- filter.bts(c( (M%*%matrix(new.wc[1:3], ncol=1))[2:3,1], new.wc[4], (M%*%matrix(new.wl[1:3], ncol=1))[2:3,1], new.wl[4]))
          detail.2[,c9[i]] <- abs(cbind(t(M%*%t(new.x[,1:3,drop=F]))[,2:3,drop=F], new.x[,4])%*%matrix(h2, ncol=1))
          
          details[,c9[i]] <- apply(cbind(abs(detail.1[,c9[i]]), abs(detail.2[,c9[i]])), 1, max)
          
        }
      }
      
      ##### Remove the change point having a smalleast |detail|
      maxdet <- apply(details, 2, max)
      if(maxdet[which.min(maxdet)] < lambda){
        ### update
        chp <- chp[-c(which.min(maxdet)+1)]
        
        for(k in 1:(length(chp)-1)){
          domain <- c(chp[k]:(chp[k+1]-1))
          pp1fit[,domain] <- t(apply(pp1fit[,domain, drop=F], 1, function(x){lm(x~domain)$fitted.values}))
        }
      } else{
        break
      }
      
      ### update
      pqr <- cbind(chp[1:(length(chp)-2)], chp[2:(length(chp)-1)]-1, chp[3:length(chp)]-1)
      d.pqr <- t(apply(pqr, 1, diff)) # difference of pqr
      details <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1])
      detail.1 <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1]) # first merging
      detail.2 <- matrix(NA, nrow=dim(pp1fit)[1], ncol=dim(d.pqr)[1]) # second merging
    }
    
    ### final change points
    bts.obj$chp <- chp[-c(1, length(chp))]
    bts.obj$details <- details
    bts.obj$bts.coeffs <- pp1fit
    
  }else{
    
    details <- matrix(0, nrow=dim(pp1fit)[1], ncol=2)
    
    bts.obj$chp <- inicp
    bts.obj$details <- details
    bts.obj$bts.coeffs <- pp1fit
  }
  
  return(bts.obj)
  
  
}

#########################################################
### all four steps (from post processing)
hd.bts.cpt <- function(x, th.const = 1, p = .01, bal = 0) {
  
  n <- dim(x)[2]
  
  if (n == 1) {
    
    est <- x
    no.of.cpt <- 0
    cpt <- integer(0)
    
  } else {
    
    varest <- apply(x, 1, function(x){(mad(diff(diff(x)))/sqrt(6))})
    x <- x/varest
    
    dcmp <- hd.bts.dcmp(x, p=p)
    lambda <- th.const * sqrt(2 * log(n*dim(x)[1])) 
    
    dns <- hd.bts.dns(dcmp, lambda=lambda, bal=bal)
    
    inv <- hd.bts.inv(dns)
    
    pp1 <- hd.bts.pp1(bts.obj=inv, lambda=lambda)
    
    cpt <- pp1$chp
    cptind <- ifelse(pp1$details > lambda, 1, 0)
    no.of.cpt <- length(cpt)
    
    est <- pp1$bts.coeffs*varest 
    
  }
  
  list(est=est, no.of.cpt=no.of.cpt, cpt=cpt, cptind=cptind)
}





#########################################################
### misc functions
finding.cp <- function(bts.obj){
  
  if(length(bts.obj$sameboat)==1){
    all.edges <- matrix(c(bts.obj$decomp.hist[1,,], bts.obj$sameboat), ncol=1) #sameboat shows pairs
  } else{
    all.edges <- rbind(bts.obj$decomp.hist[1,,], bts.obj$sameboat) #sameboat shows pairs
  }
  survived.edges <- all.edges[ ,which(abs(bts.obj$decomp.hist[3,1,])>sqrt(.Machine$double.eps)), drop=F]
  
  if(length(survived.edges)==0){
    cp <- c()
  } else if(length(survived.edges)>0 & dim(survived.edges)[2]>1){
    i <- 1
    cp <- c()
    while(i<dim(survived.edges)[2]){
      part.obj <- bts.obj$decomp.hist[1,c(1,2),] 
      matched <- which(!is.na(match(data.frame(part.obj), data.frame(matrix(survived.edges[c(2:3),i], ncol=1))))) 
      if((survived.edges[4,i]!=0) & diff(survived.edges[-4,i])[1]==1 & diff(survived.edges[-4,i+1])[1]==1){
        cp <- c(cp, survived.edges[c(1, 3), i])
        i <- i+2
      } else if((survived.edges[4,i]!=0) & diff(survived.edges[-4,i])[2]==1 & diff(survived.edges[-4,i+1])[1]==1){
        cp <- c(cp, survived.edges[c(1, 3), i+1])
        i <- i+2
      } else if((survived.edges[4,i]==0) & diff(survived.edges[-4,i])[1]==1 & diff(survived.edges[-4,i])[2]!=1){ 
        cp <-c(cp, survived.edges[c(1,3), i])
        i <- i+1
      } else if((survived.edges[4,i]==0) & diff(survived.edges[-4,i])[1]==1 & diff(survived.edges[-4,i])[2]==1 & length(matched)>0){
        cp <-c(cp, survived.edges[c(1,2), i])
        i <- i+1
      } else {
        cp <-c(cp, survived.edges[c(1,2,3), i]) 
        i <- i+1
      }
    }
    cp <- unique(sort(cp))
    cp <- c(cp, bts.obj$n+1)
  } else if(length(survived.edges)>0 & dim(survived.edges)[2]==1 & survived.edges[4,1]==0 & diff(survived.edges[-4,1])[1]==1 & diff(survived.edges[-4,1])[2]!=1){
    cp <- c()
    cp <- c(cp, survived.edges[c(1,3), 1])
  } else if(length(survived.edges)>0 & dim(survived.edges)[2]==1 & survived.edges[4,1]==0 & diff(survived.edges[-4,1])[1]==1 & diff(survived.edges[-4,1])[2]==1){ ### 3) (x, xx from a chunk)
    cp <- c()
    cp <- c(cp, survived.edges[c(1,2), 1])
  } else {
    cp <- c()
  }
  
  ### last adjustment
  if(bts.obj$n==3 & length(cp)>0 & dim(survived.edges)[2]==1){
    cp <- cp[-1]
    cp <- c(cp, bts.obj$n)
  } else {
    cp <- cp[which(cp<=bts.obj$n & cp>1)]
  }
  
  ### for comparing with other methods
  if(length(cp)>0){
    cp <- cp-1
  }
  return(cp)
}

computeDET <- function(edges = edges, edgerow = edgerow, weights.const = weights.const, weights.lin = weights.lin, bts.coeffs = bts.coeffs){
  
  sub.wc <- cbind(weights.const[edges[edgerow,1]], weights.const[edges[edgerow,2]], weights.const[edges[edgerow,3]])
  sub.wl <- cbind(weights.lin[edges[edgerow,1]], weights.lin[edges[edgerow,2]], weights.lin[edges[edgerow,3]])
  sub.tc <- cbind(bts.coeffs[edges[edgerow,1]], bts.coeffs[edges[edgerow,2]], bts.coeffs[edges[edgerow,3]])
  detcoef <- apply(cbind(sub.wc, sub.wl), 1, filter.bts)
  details <- colSums(detcoef*t(sub.tc))
  return(list(detcoef=detcoef, det=details, wc=sub.wc, wl=sub.wl, tc=sub.tc))
}

updating <- function(ee = ee, weights.const = weights.const, weights.lin = weights.lin, bts.coeffs = bts.coeffs, idx = idx){
  
  wc0 <- cbind(weights.const[ee[,1]], weights.const[ee[,2]], weights.const[ee[,3]])
  wl0 <- cbind(weights.lin[ee[,1]], weights.lin[ee[,2]], weights.lin[ee[,3]])
  tc0 <- cbind(bts.coeffs[ee[,1]], bts.coeffs[ee[,2]], bts.coeffs[ee[,3]])
  
  h <- apply(cbind(wc0,wl0), 1, filter.bts)
  M0 <- apply(h, 2, orth.matrix)
  
  wc1 <- cbind(rowSums(wc0*t(M0[c(1,4,7),])), rowSums(wc0*t(M0[c(2,5,8),])), rowSums(wc0*t(M0[c(3,6,9),])))
  wl1 <- cbind(rowSums(wl0*t(M0[c(1,4,7),])), rowSums(wl0*t(M0[c(2,5,8),])), rowSums(wl0*t(M0[c(3,6,9),])))
  tc1 <- cbind(rowSums(tc0*t(M0[c(1,4,7),])), rowSums(tc0*t(M0[c(2,5,8),])), rowSums(tc0*t(M0[c(3,6,9),])))
  
  eating.up0 <- ee[,1]
  eating.up1 <- ee[,2]
  eaten.up <- ee[,3]
  idx <- idx[-which(is.na(match(idx, c(eaten.up)))==F)]
  
  ### 1) updating X
  bts.coeffs[eating.up0] <- c(tc1[,2])
  bts.coeffs[eating.up1] <- c(tc1[,3])
  bts.coeffs[eaten.up] <- c(tc1[,1])
  
  ### 2) updating weight.const
  weights.const[eaten.up] <- c(wc1[,1])
  weights.const[eating.up0] <- c(wc1[,2])
  weights.const[eating.up1] <- c(wc1[,3])
  
  ### 3) updating weight.lin
  weights.lin[eaten.up] <- c(wl1[,1])
  weights.lin[eating.up0] <- c(wl1[,2])
  weights.lin[eating.up1] <- c(wl1[,3])
  
  return(list(weights.const=weights.const, weights.lin=weights.lin, bts.coeffs=bts.coeffs, idx=idx, h=h, tc1=tc1))
}

balance.np <- function(paired=paired, ee=ee, idx, no.of.current.steps=no.of.current.steps, n=n){
  
  prd <- !is.na(matrix(match(ee, paired), ncol=3))
  blnc <- matrix(NA, nrow=3, ncol=no.of.current.steps)
  firsttwo <- which(rowSums(prd[,1:2, drop=F])==2)
  lasttwo <- which(rowSums(prd[,2:3, drop=F])==2)
  
  if(length(c(firsttwo, lasttwo))>0){
    nopair <- c(1:dim(ee)[1])[-c(firsttwo, lasttwo)]
  } else{
    nopair <- c(1:dim(ee)[1])
  }
  
  if(length(firsttwo)>0){
    prtn <- ee[firsttwo,3] - ee[firsttwo,1]
    blnc[1:2, firsttwo] <- rbind(prtn/(prtn+1), 1/(prtn+1))
  }
  if(length(lasttwo)>0){
    prtn <- idx[match(ee[lasttwo, 3], idx)+1] - ee[lasttwo, 2]
    if(sum(is.na(prtn))>0){
      prtn[which(is.na(prtn))] <- n - ee[which(is.na(prtn)), 2] + 1
    }
    blnc[1:2, lasttwo] <- rbind(1/(prtn+1), prtn/(prtn+1))
  }
  if(length(nopair)>0){
    blnc[1:3, nopair] <- 1/3
  }
  return(blnc)
}

balance.p <- function(pr=pr, ee.p1=ee.p1, idx, ee.p2=ee.p2, n=n){
  
  blnc <- matrix(NA, nrow=3, ncol=dim(pr)[2])
  
  c1 <- ee.p1[,3] - ee.p1[,1]
  c2 <- idx[match(ee.p2[, 3], idx)+1] - ee.p1[, 3]
  if(sum(is.na(c2))>0){
    c2[which(is.na(c2))] <- n - ee.p1[which(is.na(c2)), 3] + 1
  }
  
  blnc[1:2,] <-  t(matrix(c(c1/(c1+c2), c2/(c1+c2)), nrow=2))
  
  return(blnc)
}

filter.bts <- function(a) {
  
  w <- -sqrt( (a[5]*a[1] - a[4]*a[2])^2 / ((a[2]*a[6] - a[3]*a[5])^2 + (a[4]*a[3] - a[1]*a[6])^2 + (a[5]*a[1] - a[4]*a[2])^2))
  
  u <- w * (a[2]*a[6] - a[3]*a[5])/(a[5]*a[1] - a[4]*a[2])
  
  v <- w * (a[3]*a[4] - a[1]*a[6])/(a[5]*a[1] - a[4]*a[2])
  
  df <- c(u, v, w)
  
  if (any(is.na(df))) {
    z <- filter.bts(a[6:1])
    df <- z[3:1]
  }
  return(df)
  
}

orth.matrix <- function(d) {
  
  M <- matrix(0, 3, 3)
  
  M[1,] <- d
  M[1,] <- M[1,] / sqrt(sum(M[1,]^2))
  u <- M[1, 1]
  v <- M[1, 2]
  w <- M[1, 3]
  
  M[2,] <- c(1-u^2, -u*v, -u*w)
  M[3,] <- c(0, -w, v)
  
  M[2,] <- M[2,] / sqrt(sum(M[2,]^2))
  M[3,] <- M[3,] / sqrt(sum(M[3,]^2))
  
  return(M)
  
}

L12 <- function(l){
  x <- 1:l
  n <- length(x)
  
  edges <- matrix(0, n-2, 3) # to be updated for each scale j
  for(i in 1:(n-2)){
    edges[i,] <- c((n-i-1):(n-i+1))
  }
  
  weights.const <- rep(1, n)
  weights.lin <- 1:n
  updatedS <- diag(n)
  L1L2 <- list()
  
  for(st in 1:dim(edges)[1]) {
    ee <- matrix(edges[st,], 1, 3)
    
    h <- filter.bts(c(weights.const[edges[st,]], weights.lin[edges[st,]] ))
    M <- orth.matrix(h)
    
    tmp <- matrix(0, nrow=3, ncol=2)
    tmp[,1] <- weights.const[ee]
    tmp[,2] <- weights.lin[ee]
    
    sm.det <- M %*% tmp
    updatedS[, ee[,1]] <- updatedS[, edges[st,]]%*%t(M[2,,drop=F])
    updatedS[, ee[,2]] <- updatedS[, edges[st,]]%*%t(M[3,,drop=F])
    updatedS[, ee[,3]] <- updatedS[, edges[st,]]%*%t(M[1,,drop=F])
    
    L1L2[[st]] <- updatedS[ee[,1]:n, ee[,1:2], drop=F]
    
    eating.up0 <- ee[,1]
    eating.up1 <- ee[,2]
    eaten.up <- ee[,3]
    
    weights.const[eaten.up] <- sm.det[1, 1]
    weights.const[eating.up0] <- sm.det[2, 1]
    weights.const[eating.up1] <- sm.det[3, 1]
    
    weights.lin[eaten.up] <- sm.det[1, 2]
    weights.lin[eating.up0] <- sm.det[2, 2]
    weights.lin[eating.up1] <- sm.det[3, 2]
  }
  
  return(L1L2)
}

assess <- function(object, modelnum, models, truex){
  qdiff <- length(which(object$cpts>0)) - length(models[[modelnum]]$cpt)
  mse <- mean((truex-object$fit)^2) 
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
    
    if(ovl=="complete" | (ovl=="half" & length(model$cpt)==1) | (ovl=="no" & length(model$cpt)==1) | d*sparsity==1){
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
    } else if(ovl=="half" & length(model$cpt) > 1 & d*sparsity>1){
      
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
        #slope +  mean(model$jump.size[(nrow(segments)-1):length(model$jump.size)])
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
    } else if(ovl=="no" & length(model$cpt) > 1 & d*sparsity>1){
      
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
        #slope +  mean(model$jump.size[(nrow(segments)-1):length(model$jump.size)])
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
    
    if(ovl=="complete" | (ovl=="half" & length(model$cpt)==1) | (ovl=="no" & length(model$cpt)==1) | d*sparsity==1){
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
    } else if(ovl=="half" & length(model$cpt) > 1 & d*sparsity>1){
      
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
      
    } else if(ovl=="no" & length(model$cpt) > 1 & d*sparsity>1){
      
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
################################## Method #####################################
###############################################################################
###############################################################################

hibts <- function(x, thr, bal, p){
  
  tic <- proc.time()
  object <- hd.bts.cpt(x = x, th.const = thr, bal = bal, p = p)
  toc <- proc.time()
  
  list(fit = object$est, cpts=object$cpt, elapsed=(toc-tic)[3], cptind=object$cptind)
  
}


###############################################################################
###############################################################################
################################### Models ####################################
###############################################################################
###############################################################################

model.one <-  list(name = "one",
  cpt.type = "pcwsLinMean",
  cpt = 50,
  jump.size = matrix(c(c(0), c(-2)/16), ncol=2),
  n = 100,
  start=c(-1,1/16))

model.wave <-  list(name = "wave2",
  cpt.type = "pcwsLinContMean",
  cpt = (1:9) * 20,
  jump.size = (-1)^{1:9} / 2,
  n = 20*10,
  start=c(-2,1/4))

model.mix1 <-  list(name = "mix1",
  cpt.type = "pcwsLinMean",
  cpt = (1:7) * 40,
  jump.size = matrix(c(c(1,-1,0,-1,1.5,-1,0), c(1,-1,-1,1,1,-2,2)/5), ncol=2),
  n = 320,
  start=c(-4,0))

model.mix2 <- list(name = "mix2",
  cpt.type = "pcwsLinMean",
  cpt = sort(c((1:7)*50, 100+5, 250+5)),
  jump.size = matrix(c(c(-4, 7,-5,-1,1,-6.5,5.5,2.5,0), c(0, 2,-3,2,-2,1,0,1,-2)/12), ncol=2),
  n = 400,
  start=c(4,0))

model.extremewave <-  list(name = "extremewave",
  cpt.type = "pcwsLinContMean",
  cpt = (1:24) * 20,
  jump.size = (-1)^{1:25} / 1.5,
  n = 20* 25,
  start=c(-4,1/3))

model.linsgmts <-  list(name = "linsgmts", 
  cpt.type = "pcwsLinMean",
  cpt = sort(c((1:4)*100, (1:4)*100+5)),
  jump.size = matrix(c(rep(c(6.5,-6.5-4/32), 4), rep(c(1,-1)/32,4)), ncol=2),
  n = 100*5,
  start=c(-2,0))

models <- list(model.one, model.wave, model.mix1, model.mix2, model.extremewave, model.linsgmts)


##### one-dimensional signal only 
par(mfrow=c(3,2),mar=c(2,2,2,2))

for(i in 1:length(models)){
  truex <- have.signal(models[[i]], sparsity=0.1, d=100, ovl="complete")
  x <- truex + rnorm(length(truex), sd=1)
  plot(x[1,], type="l", ylab="", xlab="t", lwd=2, ylim=range(x), col=8)
  lines(truex[1,], col=1, lwd=2)
}


##### data matrix image
modelnum=6
sgma=1 
sparsity=0.1
d=100 
ovl="complete"

truex <- have.signal(model=models[[modelnum]], sparsity, d, ovl) # whole matrix
x <- truex + matrix(rnorm(d*dim(truex)[2]), nrow=d)

par(mfrow=c(1,1),mar=c(2,2,2,2))
image(t(x), axes=F)


###############################################################################
###############################################################################
################################ SIMULATIONS ##################################
###############################################################################
###############################################################################

simhitslin <- function(N, modelnum, sgma, thr=1.2, bal=0, p=0.01, sparsity=0.1, d=10, ovl="complete"){
  
  truex <- have.signal(model=models[[modelnum]], sparsity, d, ovl) 
  n <- dim(truex)[2]
  
  par(mfrow=c(3,1),mar=c(2,2,1.5,1.5))
  
  result <- list(qdiff=matrix(NA, nrow=N, ncol=1), mse=matrix(NA, nrow=N, ncol=1), dh=matrix(NA, nrow=N, ncol=1), 
    ri=matrix(NA, nrow=N, ncol=1), time=matrix(NA, nrow=N, ncol=1))
  result <- lapply(result, function(x) {colnames(x) <- c("hibts"); x})
  
  for(K in 1:N){ 
    
    set.seed(K)
    x <- truex + matrix(rnorm(d*n), nrow=d)
    
    obj1 <- hibts(x=x, thr=thr, bal=bal, p=p) ## cpts is turned "integer(0)" if nothing is detected
    a.obj <- assess(object=obj1, modelnum=modelnum, models=models, truex=truex)
    
    g1 <- ifelse(apply(truex, 1, function(x){sum(abs(diff(diff(x))) > sqrt(.Machine$double.eps))})!=0, 1, 0)
    g2 <- ifelse(rowSums(obj1$cptind)!=0, 1, 0)
    ri <- fossil::rand.index(g1, g2)
    
    result[[1]][K,1] <- a.obj$qdiff
    result[[2]][K,1] <- a.obj$mse
    result[[3]][K,1] <- a.obj$dh
    result[[4]][K,1] <- ri
    result[[5]][K,1] <- obj1$elapsed
    
    ### plot 1
    plot(x[1,], type="l", col=8, ylim=range(x), xlab="")
    lines(truex[1, ], col=7, lwd=2, lty=1)
    lines(obj1$fit[1, ], col=1, lwd=1, lty=2)
    legend("topleft", c("obs","true","hibts"), col=c(8,7,1), lty=c(1,1,2), lwd=c(1,2,1), cex=1, bty="n", ncol=2)
    
    ### plot 2
    plot(x[sparsity*d,], type="l", col=8, ylim=range(x), xlab="")
    lines(truex[sparsity*d, ], col=7, lwd=2, lty=1)
    lines(obj1$fit[sparsity*d, ], col=1, lwd=1, lty=2)
    legend("topleft", c("obs","true","hibts"), col=c(8,7,1), lty=c(1,1,2), lwd=c(1,2,1), cex=1, bty="n", ncol=2)
    
    ### plot 3
    plot(x[d,], type="l", col=8, ylim=range(x), xlab="")
    lines(truex[d, ], col=7, lwd=2, lty=1)
    lines(obj1$fit[d, ], col=1, lwd=1, lty=2)
    legend("topleft", c("obs","true","hibts"), col=c(8,7,1), lty=c(1,1,2), lwd=c(1,2,1), cex=1, bty="n", ncol=2)
    
    print(K)
  }
  
  
  return(result)
  
}




###################################################################
######################## ovl="complete" ###########################
###################################################################

###################################################################
### N=100 / d=100
m1N100s1d100comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m1N100s1d100comlin, file="m1N100s1d100comlin.RData")
m1N100s2d100comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m1N100s2d100comlin, file="m1N100s2d100comlin.RData")
m1N100s3d100comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m1N100s3d100comlin, file="m1N100s3d100comlin.RData")

m2N100s1d100comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m2N100s1d100comlin, file="m2N100s1d100comlin.RData")
m2N100s2d100comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m2N100s2d100comlin, file="m2N100s2d100comlin.RData")
m2N100s3d100comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m2N100s3d100comlin, file="m2N100s3d100comlin.RData")

m3N100s1d100comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m3N100s1d100comlin, file="m3N100s1d100comlin.RData")
m3N100s2d100comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m3N100s2d100comlin, file="m3N100s2d100comlin.RData")
m3N100s3d100comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m3N100s3d100comlin, file="m3N100s3d100comlin.RData")

m4N100s1d100comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m4N100s1d100comlin, file="m4N100s1d100comlin.RData")
m4N100s2d100comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m4N100s2d100comlin, file="m4N100s2d100comlin.RData")
m4N100s3d100comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m4N100s3d100comlin, file="m4N100s3d100comlin.RData")

m5N100s1d100comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m5N100s1d100comlin, file="m5N100s1d100comlin.RData")
m5N100s2d100comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m5N100s2d100comlin, file="m5N100s2d100comlin.RData")
m5N100s3d100comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m5N100s3d100comlin, file="m5N100s3d100comlin.RData")

m6N100s1d100comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="complete")
save(m6N100s1d100comlin, file="m6N100s1d100comlin.RData")
m6N100s2d100comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="complete")
save(m6N100s2d100comlin, file="m6N100s2d100comlin.RData")
m6N100s3d100comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="complete")
save(m6N100s3d100comlin, file="m6N100s3d100comlin.RData")

###################################################################
### N=100 / d=300
m1N100s1d300comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m1N100s1d300comlin, file="m1N100s1d300comlin.RData")
m1N100s2d300comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m1N100s2d300comlin, file="m1N100s2d300comlin.RData")
m1N100s3d300comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m1N100s3d300comlin, file="m1N100s3d300comlin.RData")

m2N100s1d300comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m2N100s1d300comlin, file="m2N100s1d300comlin.RData")
m2N100s2d300comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m2N100s2d300comlin, file="m2N100s2d300comlin.RData")
m2N100s3d300comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m2N100s3d300comlin, file="m2N100s3d300comlin.RData")

m3N100s1d300comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m3N100s1d300comlin, file="m3N100s1d300comlin.RData")
m3N100s2d300comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m3N100s2d300comlin, file="m3N100s2d300comlin.RData")
m3N100s3d300comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m3N100s3d300comlin, file="m3N100s3d300comlin.RData")

m4N100s1d300comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m4N100s1d300comlin, file="m4N100s1d300comlin.RData")
m4N100s2d300comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m4N100s2d300comlin, file="m4N100s2d300comlin.RData")
m4N100s3d300comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m4N100s3d300comlin, file="m4N100s3d300comlin.RData")

m5N100s1d300comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m5N100s1d300comlin, file="m5N100s1d300comlin.RData")
m5N100s2d300comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m5N100s2d300comlin, file="m5N100s2d300comlin.RData")
m5N100s3d300comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m5N100s3d300comlin, file="m5N100s3d300comlin.RData")

m6N100s1d300comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="complete")
save(m6N100s1d300comlin, file="m6N100s1d300comlin.RData")
m6N100s2d300comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="complete")
save(m6N100s2d300comlin, file="m6N100s2d300comlin.RData")
m6N100s3d300comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="complete")
save(m6N100s3d300comlin, file="m6N100s3d300comlin.RData")

###################################################################
### N=100 / d=500
m1N100s1d500comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m1N100s1d500comlin, file="m1N100s1d500comlin.RData")
m1N100s2d500comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m1N100s2d500comlin, file="m1N100s2d500comlin.RData")
m1N100s3d500comlin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m1N100s3d500comlin, file="m1N100s3d500comlin.RData")

m2N100s1d500comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m2N100s1d500comlin, file="m2N100s1d500comlin.RData")
m2N100s2d500comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m2N100s2d500comlin, file="m2N100s2d500comlin.RData")
m2N100s3d500comlin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m2N100s3d500comlin, file="m2N100s3d500comlin.RData")

m3N100s1d500comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m3N100s1d500comlin, file="m3N100s1d500comlin.RData")
m3N100s2d500comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m3N100s2d500comlin, file="m3N100s2d500comlin.RData")
m3N100s3d500comlin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m3N100s3d500comlin, file="m3N100s3d500comlin.RData")

m4N100s1d500comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m4N100s1d500comlin, file="m4N100s1d500comlin.RData")
m4N100s2d500comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m4N100s2d500comlin, file="m4N100s2d500comlin.RData")
m4N100s3d500comlin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m4N100s3d500comlin, file="m4N100s3d500comlin.RData")

m5N100s1d500comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m5N100s1d500comlin, file="m5N100s1d500comlin.RData")
m5N100s2d500comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m5N100s2d500comlin, file="m5N100s2d500comlin.RData")
m5N100s3d500comlin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m5N100s3d500comlin, file="m5N100s3d500comlin.RData")

m6N100s1d500comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="complete")
save(m6N100s1d500comlin, file="m6N100s1d500comlin.RData")
m6N100s2d500comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="complete")
save(m6N100s2d500comlin, file="m6N100s2d500comlin.RData")
m6N100s3d500comlin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="complete")
save(m6N100s3d500comlin, file="m6N100s3d500comlin.RData")



###################################################################
########################### ovl="half" ############################
###################################################################

###################################################################
### N=100 / d=100
m1N100s1d100halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m1N100s1d100halflin, file="m1N100s1d100halflin.RData")
m1N100s2d100halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m1N100s2d100halflin, file="m1N100s2d100halflin.RData")
m1N100s3d100halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m1N100s3d100halflin, file="m1N100s3d100halflin.RData")

m2N100s1d100halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m2N100s1d100halflin, file="m2N100s1d100halflin.RData")
m2N100s2d100halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m2N100s2d100halflin, file="m2N100s2d100halflin.RData")
m2N100s3d100halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m2N100s3d100halflin, file="m2N100s3d100halflin.RData")

m3N100s1d100halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m3N100s1d100halflin, file="m3N100s1d100halflin.RData")
m3N100s2d100halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m3N100s2d100halflin, file="m3N100s2d100halflin.RData")
m3N100s3d100halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m3N100s3d100halflin, file="m3N100s3d100halflin.RData")

m4N100s1d100halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m4N100s1d100halflin, file="m4N100s1d100halflin.RData")
m4N100s2d100halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m4N100s2d100halflin, file="m4N100s2d100halflin.RData")
m4N100s3d100halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m4N100s3d100halflin, file="m4N100s3d100halflin.RData")

m5N100s1d100halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m5N100s1d100halflin, file="m5N100s1d100halflin.RData")
m5N100s2d100halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m5N100s2d100halflin, file="m5N100s2d100halflin.RData")
m5N100s3d100halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m5N100s3d100halflin, file="m5N100s3d100halflin.RData")

m6N100s1d100halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="half")
save(m6N100s1d100halflin, file="m6N100s1d100halflin.RData")
m6N100s2d100halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="half")
save(m6N100s2d100halflin, file="m6N100s2d100halflin.RData")
m6N100s3d100halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="half")
save(m6N100s3d100halflin, file="m6N100s3d100halflin.RData")

###################################################################
### N=100 / d=300
m1N100s1d300halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m1N100s1d300halflin, file="m1N100s1d300halflin.RData")
m1N100s2d300halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m1N100s2d300halflin, file="m1N100s2d300halflin.RData")
m1N100s3d300halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m1N100s3d300halflin, file="m1N100s3d300halflin.RData")

m2N100s1d300halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m2N100s1d300halflin, file="m2N100s1d300halflin.RData")
m2N100s2d300halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m2N100s2d300halflin, file="m2N100s2d300halflin.RData")
m2N100s3d300halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m2N100s3d300halflin, file="m2N100s3d300halflin.RData")

m3N100s1d300halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m3N100s1d300halflin, file="m3N100s1d300halflin.RData")
m3N100s2d300halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m3N100s2d300halflin, file="m3N100s2d300halflin.RData")
m3N100s3d300halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m3N100s3d300halflin, file="m3N100s3d300halflin.RData")

m4N100s1d300halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m4N100s1d300halflin, file="m4N100s1d300halflin.RData")
m4N100s2d300halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m4N100s2d300halflin, file="m4N100s2d300halflin.RData")
m4N100s3d300halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m4N100s3d300halflin, file="m4N100s3d300halflin.RData")

m5N100s1d300halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m5N100s1d300halflin, file="m5N100s1d300halflin.RData")
m5N100s2d300halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m5N100s2d300halflin, file="m5N100s2d300halflin.RData")
m5N100s3d300halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m5N100s3d300halflin, file="m5N100s3d300halflin.RData")

m6N100s1d300halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="half")
save(m6N100s1d300halflin, file="m6N100s1d300halflin.RData")
m6N100s2d300halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="half")
save(m6N100s2d300halflin, file="m6N100s2d300halflin.RData")
m6N100s3d300halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="half")
save(m6N100s3d300halflin, file="m6N100s3d300halflin.RData")

###################################################################
### N=100 / d=500
m1N100s1d500halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m1N100s1d500halflin, file="m1N100s1d500halflin.RData")
m1N100s2d500halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m1N100s2d500halflin, file="m1N100s2d500halflin.RData")
m1N100s3d500halflin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m1N100s3d500halflin, file="m1N100s3d500halflin.RData")

m2N100s1d500halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m2N100s1d500halflin, file="m2N100s1d500halflin.RData")
m2N100s2d500halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m2N100s2d500halflin, file="m2N100s2d500halflin.RData")
m2N100s3d500halflin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m2N100s3d500halflin, file="m2N100s3d500halflin.RData")

m3N100s1d500halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m3N100s1d500halflin, file="m3N100s1d500halflin.RData")
m3N100s2d500halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m3N100s2d500halflin, file="m3N100s2d500halflin.RData")
m3N100s3d500halflin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m3N100s3d500halflin, file="m3N100s3d500halflin.RData")

m4N100s1d500halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m4N100s1d500halflin, file="m4N100s1d500halflin.RData")
m4N100s2d500halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m4N100s2d500halflin, file="m4N100s2d500halflin.RData")
m4N100s3d500halflin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m4N100s3d500halflin, file="m4N100s3d500halflin.RData")

m5N100s1d500halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m5N100s1d500halflin, file="m5N100s1d500halflin.RData")
m5N100s2d500halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m5N100s2d500halflin, file="m5N100s2d500halflin.RData")
m5N100s3d500halflin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m5N100s3d500halflin, file="m5N100s3d500halflin.RData")

m6N100s1d500halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="half")
save(m6N100s1d500halflin, file="m6N100s1d500halflin.RData")
m6N100s2d500halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="half")
save(m6N100s2d500halflin, file="m6N100s2d500halflin.RData")
m6N100s3d500halflin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="half")
save(m6N100s3d500halflin, file="m6N100s3d500halflin.RData")





###################################################################
########################### ovl="no" ##############################
###################################################################


###################################################################
### N=100 / d=100
m1N100s1d100nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m1N100s1d100nolin, file="m1N100s1d100nolin.RData")
m1N100s2d100nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m1N100s2d100nolin, file="m1N100s2d100nolin.RData")
m1N100s3d100nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m1N100s3d100nolin, file="m1N100s3d100nolin.RData")

m2N100s1d100nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m2N100s1d100nolin, file="m2N100s1d100nolin.RData")
m2N100s2d100nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m2N100s2d100nolin, file="m2N100s2d100nolin.RData")
m2N100s3d100nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m2N100s3d100nolin, file="m2N100s3d100nolin.RData")

m3N100s1d100nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m3N100s1d100nolin, file="m3N100s1d100nolin.RData")
m3N100s2d100nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m3N100s2d100nolin, file="m3N100s2d100nolin.RData")
m3N100s3d100nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m3N100s3d100nolin, file="m3N100s3d100nolin.RData")

m4N100s1d100nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m4N100s1d100nolin, file="m4N100s1d100nolin.RData")
m4N100s2d100nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m4N100s2d100nolin, file="m4N100s2d100nolin.RData")
m4N100s3d100nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m4N100s3d100nolin, file="m4N100s3d100nolin.RData")

m5N100s1d100nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m5N100s1d100nolin, file="m5N100s1d100nolin.RData")
m5N100s2d100nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m5N100s2d100nolin, file="m5N100s2d100nolin.RData")
m5N100s3d100nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m5N100s3d100nolin, file="m5N100s3d100nolin.RData")

m6N100s1d100nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=100, ovl="no")
save(m6N100s1d100nolin, file="m6N100s1d100nolin.RData")
m6N100s2d100nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=100, ovl="no")
save(m6N100s2d100nolin, file="m6N100s2d100nolin.RData")
m6N100s3d100nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=100, ovl="no")
save(m6N100s3d100nolin, file="m6N100s3d100nolin.RData")

###################################################################
### N=100 / d=300
m1N100s1d300nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m1N100s1d300nolin, file="m1N100s1d300nolin.RData")
m1N100s2d300nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m1N100s2d300nolin, file="m1N100s2d300nolin.RData")
m1N100s3d300nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m1N100s3d300nolin, file="m1N100s3d300nolin.RData")

m2N100s1d300nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m2N100s1d300nolin, file="m2N100s1d300nolin.RData")
m2N100s2d300nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m2N100s2d300nolin, file="m2N100s2d300nolin.RData")
m2N100s3d300nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m2N100s3d300nolin, file="m2N100s3d300nolin.RData")

m3N100s1d300nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m3N100s1d300nolin, file="m3N100s1d300nolin.RData")
m3N100s2d300nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m3N100s2d300nolin, file="m3N100s2d300nolin.RData")
m3N100s3d300nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m3N100s3d300nolin, file="m3N100s3d300nolin.RData")

m4N100s1d300nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m4N100s1d300nolin, file="m4N100s1d300nolin.RData")
m4N100s2d300nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m4N100s2d300nolin, file="m4N100s2d300nolin.RData")
m4N100s3d300nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m4N100s3d300nolin, file="m4N100s3d300nolin.RData")

m5N100s1d300nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m5N100s1d300nolin, file="m5N100s1d300nolin.RData")
m5N100s2d300nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m5N100s2d300nolin, file="m5N100s2d300nolin.RData")
m5N100s3d300nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m5N100s3d300nolin, file="m5N100s3d300nolin.RData")

m6N100s1d300nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=300, ovl="no")
save(m6N100s1d300nolin, file="m6N100s1d300nolin.RData")
m6N100s2d300nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=300, ovl="no")
save(m6N100s2d300nolin, file="m6N100s2d300nolin.RData")
m6N100s3d300nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=300, ovl="no")
save(m6N100s3d300nolin, file="m6N100s3d300nolin.RData")

###################################################################
### N=100 / d=500
m1N100s1d500nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m1N100s1d500nolin, file="m1N100s1d500nolin.RData")
m1N100s2d500nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m1N100s2d500nolin, file="m1N100s2d500nolin.RData")
m1N100s3d500nolin <- simhitslin(N=100, modelnum=1, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m1N100s3d500nolin, file="m1N100s3d500nolin.RData")

m2N100s1d500nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m2N100s1d500nolin, file="m2N100s1d500nolin.RData")
m2N100s2d500nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m2N100s2d500nolin, file="m2N100s2d500nolin.RData")
m2N100s3d500nolin <- simhitslin(N=100, modelnum=2, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m2N100s3d500nolin, file="m2N100s3d500nolin.RData")

m3N100s1d500nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m3N100s1d500nolin, file="m3N100s1d500nolin.RData")
m3N100s2d500nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m3N100s2d500nolin, file="m3N100s2d500nolin.RData")
m3N100s3d500nolin <- simhitslin(N=100, modelnum=3, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m3N100s3d500nolin, file="m3N100s3d500nolin.RData")

m4N100s1d500nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m4N100s1d500nolin, file="m4N100s1d500nolin.RData")
m4N100s2d500nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m4N100s2d500nolin, file="m4N100s2d500nolin.RData")
m4N100s3d500nolin <- simhitslin(N=100, modelnum=4, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m4N100s3d500nolin, file="m4N100s3d500nolin.RData")

m5N100s1d500nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m5N100s1d500nolin, file="m5N100s1d500nolin.RData")
m5N100s2d500nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m5N100s2d500nolin, file="m5N100s2d500nolin.RData")
m5N100s3d500nolin <- simhitslin(N=100, modelnum=5, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m5N100s3d500nolin, file="m5N100s3d500nolin.RData")

m6N100s1d500nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.01, d=500, ovl="no")
save(m6N100s1d500nolin, file="m6N100s1d500nolin.RData")
m6N100s2d500nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.1, d=500, ovl="no")
save(m6N100s2d500nolin, file="m6N100s2d500nolin.RData")
m6N100s3d500nolin <- simhitslin(N=100, modelnum=6, sgma=1, thr=1.2, bal=0, p=0.04, sparsity=0.7, d=500, ovl="no")
save(m6N100s3d500nolin, file="m6N100s3d500nolin.RData")






###############################################################################
###############################################################################
############################# DATA APPLICATION ################################
###############################################################################
###############################################################################


############################################################################
############################## seaice DATA #################################
############################################################################
####### http://nsidc.org/data/nsidc-0051.html
####### https://www.kaggle.com/nsidcorg/daily-sea-ice-extent-data
####### ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/south/daily/data/
####### ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/
############################################################################

##### north
dat <- read.csv(file="N_seaice_extent_daily.csv", header=TRUE, sep=",")
dat <- dat[-1,-6, drop=F] # remove the website address column
# change the factor to numeric
dat$Year <- as.numeric(as.character(dat$Year))
dat$Month <- as.numeric(as.character(dat$Month))
dat$Day <- as.numeric(as.character(dat$Day))
dat$Extent <- as.numeric(as.character(dat$Extent))
head(dat)


##### south
dat <- read.csv(file="S_seaice_extent_daily.csv", header=TRUE, sep=",")
dat <- dat[-1,-6, drop=F] # remove the website address column
# change the factor to numeric
dat$Year <- as.numeric(as.character(dat$Year))
dat$Month <- as.numeric(as.character(dat$Month))
dat$Day <- as.numeric(as.character(dat$Day))
dat$Extent <- as.numeric(as.character(dat$Extent))
head(dat)


##########################################################
######## data matrix for each year and month #############
##########################################################
dat <- dat[dat$Year!=1978,]
dat <- dat[dat$Year!=2019,]

dat$ym <- paste(dat$Year, dat$Month)
dat.mean <- dat %>% group_by(ym) %>% summarise(mean.extent = mean(Extent))
dat.mean <- separate(dat.mean, ym, into = c('year', 'month'), sep=" ")

allmonth.dat <- list()
for(k in 1:12){
  dat.eachmonth <- dat.mean[dat.mean$month==k,]
  allmonth.dat[[k]] <- rbind(as.numeric(dat.eachmonth$year), dat.eachmonth$mean.extent)
}

yr <- unique(dat$Year)
commonyear <- matrix(0, nrow=length(yr), ncol=13)
commonyear[,1] <- yr
for(k in 1:12){
  t <- which(allmonth.dat[[k]][1,] %in% yr)
  commonyear[, k+1] <- allmonth.dat[[k]][2,t]
}

#### plotting
mth <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = T))
par(mfrow=c(1,1),mar=c(5,5,3,3))

plot(commonyear[,1], commonyear[,2], type="l", ylim=range(commonyear[,-1]), ylab="ice extent", xlab="Year", main="Antarctic")
for(k in 2:12){
  lines(commonyear[,1], commonyear[,k+1], type="l", col=k, lty=k)
  text(commonyear[40,1], commonyear[40,k+1],  mth[k], cex=0.65, pos=4)
}
legend("bottomleft", mth, col=c(1:12), lty=c(1:12), lwd=1, bty="n", ncol=6)


###############################################
########## change-point detection #############
###############################################

X <- t(commonyear[,-1,drop=F])
d <- dim(X)[1]
n <- dim(X)[2]
d;n

hitspwl <- hibts(x=X, thr=1.2, p=0.04, bal=0)
hitspwl$cpts
commonyear[hitspwl$cpts]
cind <- hitspwl$cptind
cind

mth <- c("January", "February", "March", "April", "May", "June",
  "July", "August", "September", "October", "November", "December")
### from the post processing
par(mfrow=c(4,3),mar=c(3,3,3,3))
for(k in 1:12){
  plot(yr, X[k,], type="b", col="darkgrey", ylab="", xlab="", main=paste("Antarctic in", mth[k]))
  lines(yr, hitspwl$fit[k,], type="l", col=1, lwd=2, lty=3) #pp-stage 1
  abline(v=yr[hitspwl$cpts], lty=1, col=3, lwd=1) #pp-stage 1
}


#### all months in one (Arctic)
mth <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = T))
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(yr, X[1,], type="l", col=1, ylab="ice extent", xlab="Year", ylim=range(commonyear[,-1]), xlim=c(yr[1], yr[length(yr)]+1))
text(commonyear[40,1], commonyear[40,1+1],  mth[1], cex=1.2, pos=4)

for(k in 2:12){
  lines(yr, X[k,], type="l", col=1)
  if(k==4 | k==12){
    text(commonyear[40,1], commonyear[40,k+1]-0.2,  mth[k], cex=1.2, pos=4)
  } else{
    text(commonyear[40,1], commonyear[40,k+1],  mth[k], cex=1.2, pos=4)
  }
}
abline(v=yr[hitspwl$cpts], lty=2, col=2, lwd=2) #pp-stage 1


#### all months in one (Antarctic)
mth <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
layout(matrix(c(1,1,2,3,4,5), 3, 2, byrow = T))
par(mfrow=c(1,1),mar=rep(2.5,4))
plot(yr, X[1,], type="l", col=1, ylab="ice extent", xlab="Year", ylim=range(commonyear[,-1]), xlim=c(yr[1], yr[length(yr)]+1))
text(commonyear[40,1], commonyear[40,1+1],  mth[1], cex=1.2, pos=4)

for(k in 2:12){
  lines(yr, X[k,], type="l", col=1)
  if(k==5 | k==9){
    text(commonyear[40,1], commonyear[40,k+1]+0.25,  mth[k], cex=1.2, pos=4)
  } else if(k==8){
    text(commonyear[40,1], commonyear[40,k+1]-0.2,  mth[k], cex=1.2, pos=4)
  } else{
    text(commonyear[40,1], commonyear[40,k+1],  mth[k], cex=1.2, pos=4) 
  }
}
abline(v=yr[hitspwl$cpts], lty=2, col=2, lwd=2) #pp-stage 1



### from the cpt index for each city
eachfit <- matrix(NA, nrow=d, ncol=n)
indc <- ifelse(cind==1, TRUE, FALSE)
for(i in 1:d){
  if(length(hitspwl$cpts[indc[i,]])==0){
    eachfit[i, ] <- lm(X[i,]~c(1:n))$fitted.values 
  }else{
    cp <- c(1, hitspwl$cpts[indc[i,]]+1, n+1)
    for(k in 1:(length(cp)-1)){
      domain <- c(cp[k]:(cp[k+1]-1))
      eachfit[i, domain] <- lm(X[i,domain]~domain)$fitted.values 
    }
  }
}

### Arctic
par(mfrow=c(4,3),mar=c(2,3,2,1))
for(k in 1:12){
  plot(yr, X[k,], type="b", col=8, lwd=2, ylab="", xlab="", main=paste("Arctic in", mth[k]))
  lines(yr, eachfit[k,], type="l", col=2, lwd=2, lty=2) #pp-stage 1
  if(k==1 | k==2){
    abline(v=yr[hitspwl$cpts][c(2,3)], lty=1, col=4, lwd=1) 
  } else if(k==12){
    abline(v=yr[hitspwl$cpts][1], lty=1, col=4, lwd=1) 
  }
}


### Antarctic
par(mfrow=c(4,3),mar=c(2,3,2,2))
for(k in 1:12){
  plot(yr, X[k,], type="b", col=8, lwd=2, ylab="", xlab="", main=paste("Antarctic in", mth[k]))
  lines(yr, eachfit[k,], type="l", col=2, lwd=2, lty=2) #pp-stage 1
  if(k==2 | k==3| k==5| k==6| k==7| k==8| k==9| k==11| k==12){
    abline(v=yr[hitspwl$cpts][c(3)], lty=1, col=4, lwd=1) 
  } else if(k==4){
    abline(v=yr[hitspwl$cpts][c(1,3)], lty=1, col=4, lwd=1) 
  } else if(k==10){
    abline(v=yr[hitspwl$cpts][c(2,3)], lty=1, col=4, lwd=1) 
  }
}
















