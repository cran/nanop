displacePart <- function(nanop, sigma=NA, sigmacore=NA, sigmashell=NA,
                         rcore=NA, rcenter=FALSE, center=c(0,0,0),
                         latticep=4.08) {
  onanop <- nanop 
  if(length(nanop) < 1)
    return() 
  if(rcenter)
    if(is.null(attributes(nanop)$center))
      center <- runif(3,max=latticep)
    else
      center <- attributes(nanop)$center
   
  
  if((is.na(sigmacore) || is.na(sigmashell)) && is.numeric(sigma))  {
    ## uniform
    #np <- rmvnorm(nrow(nanop), mean = rep(0,3), sigma = diag(rep(sigma,3)))
    np <- rnorm(nrow(nanop)*3, mean = 0, sd = sqrt(sigma))
    nanop <- nanop + np
  }
  else if(is.numeric(sigmacore)&&is.numeric(sigmashell)&&
          is.numeric(rcore)) {
    dist <- sqrt(rowSums((t( t(nanop)-center ))^2))
    rc <- which(dist<=rcore)

    if(length(rc) > 0) { ## if there are atoms in core
      nanop_rc <- nanop[rc,]
      if(length(rc) == 1) ## one atom in core 
        nanop_rc <- as.matrix(t(nanop_rc))
      #np_rc <-  rmvnorm(nrow(nanop_rc), mean =  rep(0,3),
      #                  sigma = diag(rep(sigmacore,3)))
      np_rc <- rnorm(nrow(nanop_rc)*3, mean = 0, sd=sqrt(sigmacore)) 
      nanop_rc <- nanop_rc + np_rc
    }
    else nanop_rc <- matrix(nrow=0, ncol=3)
    
    if(length(rc) < nrow(nanop)) { ## if there are atoms in shell
      if(length(rc) == 0)
        nanop_rs <- nanop
      else
        nanop_rs <- nanop[-rc,]
      if(length(nanop_rs) == 3) ## one atom in shell
        nanop_rs <- as.matrix(t(nanop_rs))
      
      #np_rs <-  rmvnorm(nrow(nanop_rs), mean =  rep(0,3),
      #                  sigma = diag(rep(sigmashell,3)))
      np_rs <- rnorm(nrow(nanop_rs)*3, mean=0, sd=sqrt(sigmashell))
      nanop_rs <- nanop_rs + np_rs
    }
    else nanop_rs <-  matrix(nrow=0, ncol=3)
    nanop <- rbind(nanop_rc,nanop_rs)
  }
  else 
    stop("Non-sensible input argument")
  attributes(nanop) <- attributes(onanop)
  nanop
}
