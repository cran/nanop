calcQDepPDF <- function(nanop, dr=.1, minR=1, maxR=20, minQ=1, maxQ=20,
           a1 = 16.8819,
           b1=.4611, a2=18.5913, b2=8.6216,
           a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
           c=12.0658) {
 
    r <- seq(minR, maxR, by = dr)
    res <- vector(length=length(r))
    for(i in 1:length(r)) {
      
      res[i] <- integrate(f=calcQDepPDFAux, lower=minQ, upper=maxQ,
                          r=r[i], nanop=nanop, dr=dr, minR=minR, a1=a1,
                          b1=b1, a2=a2, b2=b2, a3=a3, b3=b3, a4=a4, b4=b4,
                          c=c, subdivisions=100)$value * (2/pi)
    }
    list(r=r,gr=res)
  }
calcQDepPDFAux <- function(Q, nanop, r, a1 = 16.8819,
                           b1=.4611, a2=18.5913, b2=8.6216,
                           a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                           c=12.0658, dr, minR) {
  
  .C("calcQDepPDF",
     res = as.double(Q),
     Q = as.double(Q),
     r = as.double(r), 
     len = as.integer(length(Q)),
     np = as.double(as.vector(t(nanop))),
     nrow = as.integer(nrow(nanop)),
     a1 = as.double(a1),
     b1 = as.double(b1),
     a2 = as.double(a2),
     b2 = as.double(b2),
     a3 = as.double(a3),
     b3 = as.double(b3),
     a4 = as.double(a4),
     b4 = as.double(b4),
     c = as.double(c),
     PACKAGE="nanop")$res
  
}
