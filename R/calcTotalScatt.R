calcTotalScatt <- function(nanop, 
                           dQ=.01, minQ=1,maxQ=20, a1 =
                           16.8819, b1=.4611, a2=18.5913, b2=8.6216,
                           a3=25.5582, b3=1.48260,
                           a4=5.86, b4=36.3956, c=12.0658) {
  ## Q = 4*pi*sin(theta)/lambda
  Q <- seq(minQ,maxQ,by=dQ)
  list(Q=Q, gQ=.C("calcTotalScatt",
                res = as.double(Q),
                Q = as.double(Q), 
                len = as.integer(length(Q)),
                minQ = as.double(minQ),
                dQ = as.double(dQ),
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
                PACKAGE="nanop")$res)
 
}

