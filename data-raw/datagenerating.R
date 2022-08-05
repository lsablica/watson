library(microbenchmark)
library(watson)

ns = c(1, 3, 5, 10, 20, 50, 100, 500, 1000, 10000)
kappas <- c(-100, -50, -10, -1, 1, 10, 50, 100)
ds <- c( 3, 5, 10, 20, 50, 100, 200, 1000)
resultTinflex <- array(0, dim = c(length(ns), length(kappas), length(ds)), dimnames = list(ns,kappas,ds))
resultACG <- array(0, dim = c(length(ns), length(kappas), length(ds)), dimnames = list(ns,kappas,ds))

for(nn in seq_along(ns)){
   n = ns[nn]
   set.seed(1)
   print(n)
   for(dd in seq_along(ds)){
      d = ds[dd]
      for(kappaa in seq_along(kappas)){
         kappa = kappas[kappaa]
         mu = rep(1, times = d)
         mu = mu/sqrt(sum(mu^2))
         m = microbenchmark(watson:::rwatTinflex(n, kappa, mu, cT = 0, rho = 1.1),  watson:::rwatACG(n, kappa, mu) , times = 150L, unit = "ms")
         mm=summary(m)
         or = order(mm[,"expr"])
         resultTinflex[nn,kappaa,dd] <- mm[,"mean"][or][1]
         resultACG[nn,kappaa,dd] <- mm[,"mean"][or][2]
         print(c(d, kappa))
      }
   }
}
