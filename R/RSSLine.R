RSSLine <- function(ref, rcore, sigc, sigs, r, 
                    latticep = 4.08, 
                    calpha=1, dr=.01, minR=1, maxR=20, p = 1,
                    dQ=.01, minQ=1,maxQ=20, a1 = 16.8819,
                    b1=.4611, a2=18.5913, b2=8.6216,
                    a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                    c=12.0658, calcPDF=TRUE, calcScatt=TRUE, norm=TRUE,
                    avRes = 1, rcenter=TRUE) {
 
  if(length(r) > 1)
    ans <- RSSLine_r(ref=ref, rcore=rcore, sigc=sigc, sigs=sigs, r=r, latticep
                     = latticep,  calpha=calpha,
                     dr=dr, minR=minR, maxR=maxR, p = p, dQ=dQ,
                     minQ=minQ,maxQ=maxQ, a1 = a1, b1=b1, a2=a2,
                     b2=b2, a3=a3, b3=b3, a4=a4, b4=b4, c=12.0658,
                     calcPDF=calcPDF, calcScatt=calcScatt, norm=norm,
                     avRes = avRes, rcenter=rcenter)
  else if(length(rcore) > 1)
    ans <- RSSLine_rcore(ref=ref,
                         rcore=rcore, sigc=sigc, sigs=sigs, r=r, latticep
                         = latticep,  calpha=calpha,
                         dr=dr, minR=minR, maxR=maxR, p = p, dQ=dQ,
                         minQ=minQ,maxQ=maxQ, a1 = a1, b1=b1, a2=a2,
                         b2=b2, a3=a3, b3=b3, a4=a4, b4=b4, c=12.0658,
                         calcPDF=calcPDF, calcScatt=calcScatt, norm=norm,
                         avRes = avRes, rcenter=rcenter)
  else if(length(latticep) > 1)
    ans <- RSSLine_latticep(ref=ref,
                            rcore=rcore, sigc=sigc, sigs=sigs, r=r, latticep
                            = latticep, calpha=calpha,
                            dr=dr, minR=minR, maxR=maxR, p = p, dQ=dQ,
                            minQ=minQ,maxQ=maxQ, a1 = a1, b1=b1, a2=a2,
                            b2=b2, a3=a3, b3=b3, a4=a4, b4=b4, c=12.0658,
                            calcPDF=calcPDF, calcScatt=calcScatt, norm=norm,
                            avRes = avRes, rcenter=rcenter)
  else stop("At least one of rcore, r, latticep must have length >1.")
  ans
}
RSSLine_rcore <- function(ref, rcore, sigc, sigs, r, 
                    latticep = 4.08, 
                    calpha=1, dr=.01, minR=1, maxR=20, p = 1,
                    dQ=.01, minQ=1,maxQ=20, a1 = 16.8819,
                    b1=.4611, a2=18.5913, b2=8.6216,
                    a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                    c=12.0658, calcPDF=TRUE, calcScatt=TRUE, norm=TRUE,
                    avRes = 1, rcenter=TRUE) {
  
  if(calcPDF) {
    pdfR <- 0 
    for(i in 1:avRes) {
      if(i == 1 || rcenter) 
        aa <- simPart(r=r, latticep = latticep, rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                         rcore=ref)
      pdfR <- pdfR + calcPDF(nanop=pd, calpha=calpha, dr=dr, minR=minR,
                             maxR=maxR, p=p)$gr
    }
    pdfR <- pdfR / avRes 
    if(norm)
      sr <- sum(pdfR)
  }
  if(calcScatt) {
    debQ <- 0
    for(i in 1:avRes) {
      if(i == 1 || rcenter) 
        aa <- simPart(r=r, latticep = latticep, rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                          rcore=ref)
      debQ <- debQ +
        calcTotalScatt(nanop=pd, dQ=dQ, minQ=minQ,maxQ=maxQ, a1=a1,
                          b1=b1, a2=a2, b2=b2, a3=a3, b3=b3, a4=a4, b4=b4,
                          c=c)$gQ
    }
    debQ <- debQ / avRes 
    if(norm)
      sq <- sum(debQ)
  }
  ansR <- ansQ <- vector(length=length(rcore))
  rssR <- rssQ <- list() 
  for(i in 1:length(rcore)) {
    pdf1 <- deb1 <- 0
    for(j in 1:avRes) {
      if(rcenter) 
        aa <- simPart(r=r, latticep = latticep, rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                         rcore=rcore[i])
      if(calcPDF) 
        pdf1 <- pdf1 + calcPDF(nanop=pd, calpha=calpha, dr=dr, minR=minR, maxR=
                               maxR, p=p)$gr
        
      if(calcScatt) 
        deb1 <- deb1 +calcTotalScatt(nanop=pd, dQ=dQ,
                                     minQ=minQ,maxQ=maxQ, a1=a1,
                                     b1=b1, a2=a2, b2=b2, a3=a3,
                                     b3=b3, a4=a4, b4=b4,c=c)$gQ
    }
    if(calcPDF) {
      pdf1 <- pdf1 / avRes
      ansR[i] <- sum( (pdfR - pdf1)^2 ) 
      if(norm) ansR[i] <-  ansR[i] / (sr*length(pdfR))
     
    }
    if(calcScatt) {
      deb1 <- deb1 / avRes
      ansQ[i] <- sum( (debQ - deb1)^2 )  
      if(norm) ansQ[i] <- ansQ[i] / (sq *length(debQ))
    }
    cat("Done with rcore=",rcore[i], "\n")
  }
  list(rcore=rcore, ansR=ansR, ansQ=ansQ)
}
RSSLine_r <- function(ref, r, sigc, sigs, rcore, 
                    latticep = 4.08, 
                    calpha=1, dr=.01, minR=1, maxR=20, p = 1,
                    dQ=.01, minQ=1,maxQ=20, a1 = 16.8819,
                    b1=.4611, a2=18.5913, b2=8.6216,
                    a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                    c=12.0658, calcPDF=TRUE, calcScatt=TRUE, norm=TRUE,
                    avRes = 1, rcenter=TRUE) {
  
  if(calcPDF) {
    pdfR <- 0 
    for(i in 1:avRes) {
      if(i == 1 || rcenter) 
        aa <- simPart(r=ref, latticep = latticep, 
                      rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                         rcore=rcore)
      pdfR <- pdfR + calcPDF(nanop=pd, calpha=calpha, dr=dr, minR=minR,
                             maxR=maxR, p=p)$gr
    }
    pdfR <- pdfR / avRes 
    if(norm)
      sr <- sum(pdfR)
  }
  if(calcScatt) {
    debQ <- 0
    for(i in 1:avRes) {
      if(I == 1 || rcenter) 
        aa <- simPart(r=ref, latticep = latticep, 
                      rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                          rcore=rcore)
      debQ <- debQ +
        calcTotalScatt(nanop=pd, dQ=dQ, minQ=minQ,maxQ=maxQ, a1=a1,
                          b1=b1, a2=a2, b2=b2, a3=a3, b3=b3, a4=a4, b4=b4,
                          c=c)$gQ
    }
    debQ <- debQ / avRes 
    if(norm)
      sq <- sum(debQ)
  }
  ansR <- ansQ <- vector(length=length(rcore))
  rssR <- rssQ <- list() 
  for(i in 1:length(r)) {
    pdf1 <- deb1 <- 0
    for(j in 1:avRes) {
      if(rcenter) 
        aa <- simPart(r=r[i], latticep = latticep, 
                      rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                         rcore=rcore)
      if(calcPDF) 
        pdf1 <- pdf1 + calcPDF(nanop=pd, calpha=calpha, dr=dr, minR=minR, maxR=
                               maxR, p=p)$gr
        
      if(calcScatt) 
        deb1 <- deb1 +calcTotalScatt(nanop=pd, dQ=dQ,
                                     minQ=minQ,maxQ=maxQ, a1=a1,
                                     b1=b1, a2=a2, b2=b2, a3=a3,
                                     b3=b3, a4=a4, b4=b4,c=c)$gQ
    
      
    }
    if(calcPDF) {
      pdf1 <- pdf1 / avRes
      ansR[i] <- sum( (pdfR - pdf1)^2 ) 
      if(norm) ansR[i] <-  ansR[i] / (sr*length(pdfR))
     
    }
    if(calcScatt) {
      deb1 <- deb1 / avRes
      ansQ[i] <- sum( (debQ - deb1)^2 )  
      if(norm) ansQ[i] <- ansQ[i] / (sq *length(debQ))
    }
    cat("Done with r=",r[i], "\n")
  }
  list(r=r, ansR=ansR, ansQ=ansQ)
}
RSSLine_latticep <- function(ref, latticep, sigc, sigs, rcore, r,
                    calpha=1, dr=.01, minR=1, maxR=20, p = 1,
                    dQ=.01, minQ=1,maxQ=20, a1 = 16.8819,
                    b1=.4611, a2=18.5913, b2=8.6216,
                    a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                    c=12.0658, calcPDF=TRUE, calcScatt=TRUE, norm=TRUE,
                    avRes = 1, rcenter=TRUE) {
  
  if(calcPDF) {
    pdfR <- 0 
    for(i in 1:avRes) {
      if(i == 1 || rcenter) 
        aa <- simPart(r=r, latticep = ref, rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                         rcore=rcore)
      pdfR <- pdfR + calcPDF(nanop=pd, calpha=calpha, dr=dr, minR=minR,
                             maxR=maxR, p=p)$gr
    }
    pdfR <- pdfR / avRes 
    if(norm)
      sr <- sum(pdfR)
  }
  if(calcScatt) {
    debQ <- 0
    for(i in 1:avRes) {
      if(i == 1 || rcenter) 
        aa <- simPart(r=r, latticep = ref, rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                          rcore=rcore)
      debQ <- debQ +
        calcTotalScatt(nanop=pd, dQ=dQ, minQ=minQ,maxQ=maxQ, a1=a1,
                          b1=b1, a2=a2, b2=b2, a3=a3, b3=b3, a4=a4, b4=b4,
                          c=c)$gQ
    }
    debQ <- debQ / avRes 
    if(norm)
      sq <- sum(debQ)
  }
  ansR <- ansQ <- vector(length=length(latticep))
  rssR <- rssQ <- list() 
  for(i in 1:length(latticep)) {
    pdf1 <- deb1 <- 0
    for(j in 1:avRes) {
      
      aa <- simPart(r=r, latticep = latticep[i], 
                    rcenter=TRUE)
      pd <- displacePart(aa, sigmacore=sigc, sigmashell=sigs,
                         rcore=rcore)
      if(calcPDF) 
        pdf1 <- pdf1 + calcPDF(nanop=pd, calpha=calpha, dr=dr, minR=minR, maxR=
                               maxR, p=p)$gr
        
      if(calcScatt) 
        deb1 <- deb1 +calcTotalScatt(nanop=pd, dQ=dQ,
                                     minQ=minQ,maxQ=maxQ, a1=a1,
                                     b1=b1, a2=a2, b2=b2, a3=a3,
                                     b3=b3, a4=a4, b4=b4,c=c)$gQ
    }
    if(calcPDF) {
      pdf1 <- pdf1 / avRes
      ansR[i] <- sum( (pdfR - pdf1)^2 ) 
      if(norm) ansR[i] <-  ansR[i] / (sr*length(pdfR))
     
    }
    if(calcScatt) {
      deb1 <- deb1 / avRes
      ansQ[i] <- sum( (debQ - deb1)^2 )  
      if(norm) ansQ[i] <- ansQ[i] / (sq *length(debQ))
    }
    cat("Done with latticep=", latticep[i], "\n")
  }
  list(r=latticep, ansR=ansR, ansQ=ansQ)
}
