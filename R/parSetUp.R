simPartPar <- function(sym = "fcc", latticep = 4.08, r=10, a=5, b=5, c=5,
                       rcenter=FALSE, latticepshell=NA, rcore=NA)
  list(sym=sym, latticep=latticep, r=r, a=a,b=b,c=c,
       latticepshell = latticepshell,
       rcore = rcore, rcenter=rcenter)

PDFPar <- function(calpha=1, dr=.01, minR=1, maxR=20, p = 1)
  list(calpha=calpha, dr=dr, minR=minR, maxR=maxR, p=p)

TotalScattPar <- function(dQ=.01, minQ=1,maxQ=20, a1 = 16.8819,
                          b1=.4611, a2=18.5913, b2=8.6216,
                          a3=25.5582, b3=1.48260, a4=5.86, b4=36.3956,
                          c=12.0658)
  list(dQ=dQ, minQ=minQ,maxQ=maxQ, a1 = a1, b1=b1, a2=a2, b2=b2,
       a3=a3, b3=b3, a4=a4, b4=b4, c=c)
  
