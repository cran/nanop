simPart <- function(sym = "fcc", latticep = 4.08, r=10, 
                    latticepshell = NA, rcore = NA, rcenter=FALSE,
                    center=c(0,0,0),foranalcs = FALSE) {
  
  if(rcenter)
    center <- runif(3,max=latticep)
    ##center <- rnorm(3,0,sqrt(latticep))
    
  if(sym=="fcc") {
    coords <- list(c(0,0,0), c(.5,.5,0), c(.5,0,.5), c(0,.5,.5),
                   c(0,1,0), c(1,1,0), c(1,0,0), c(1,1,1), c(1,0,1),
                   c(0,1,1), c(0,0,1), c(.5,1,.5), c(1,.5,.5), c(.5,.5,1))
  }
  unit.cell <- matrix(unlist(coords),byrow=TRUE,nrow=length(coords), ncol=3)
  base <- unit.cell * latticep
  a <- b <- c <- ceiling(r / latticep)+2
 
  res <- rep(0,a*b*c*20*length(base))
  
  nanop <- matrix(.C("simPart", res=as.double(res),
                     lenres=as.integer(length(res)),
                     base=as.double(as.vector(t(base))),
                     lenbase=as.integer(length(base)), 
                     shiftx=as.double(latticep),
                     shifty=as.double(latticep),
                     shiftz=as.double(latticep),
                     a=as.integer(a),b=as.integer(b),
                     c=as.integer(c), PACKAGE="nanop")$res,
                  byrow=TRUE, ncol=3)
  
  nanop <- unique(nanop) 
  dist <- sqrt(rowSums((t( t(nanop)-center ))^2)) 
  if(! is.na(rcore))
    rr <- rcore
  else
    rr <- r
  ans <- nanop[which(dist<=rr),]
  
  if(! is.na(rcore)) {
    if(is.na(latticepshell))
      latticepshell <- latticep  
    a <- b <- c <- ceiling(r / latticepshell)  
    res <- rep(0,a*b*c*20*length(base))
    nanops <- matrix(.C("simPart", res=as.double(res),
                        lenres=as.integer(length(res)),
                        base=as.double(as.vector(t(base))),
                        lenbase=as.integer(length(base)), 
                        shiftx=as.double(latticepshell),
                        shifty=as.double(latticepshell),
                        shiftz=as.double(latticepshell),
                        a=as.integer(a),b=as.integer(b),
                        c=as.integer(c), PACKAGE="nanop")$res,
                     byrow=TRUE, ncol=3)
    
    nanops <- unique(nanops) 
    dists <- sqrt(rowSums((t( t(nanops)-center ))^2))
    
    rd <- which(dists <= r)
    rs <- intersect(which(dists > rcore), rd)
       
    ans <- rbind(ans, nanops[rs,])
    attr(ans, "rowcore") <- nrow(ans) -  length(rs) 
    attr(ans, "rowshell") <- length(rs) 
    
  }
  attr(ans,"center") <- center
  ans  
}
