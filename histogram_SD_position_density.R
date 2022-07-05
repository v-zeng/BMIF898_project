### Histogram and density of Standard Deviation positions for clock CpGs
library(data.table)

CpGs <- c(1,5,11,23,40,41,42,44,45,46,47,48,49,55,54,53,52,50,56,100)
hist(CpGs,breaks=50,xlim=c(0,100), freq = F) # try adding freq=F
lines(density(CpGs),
      lwd=2,
      col="blue")

# see https://stackoverflow.com/questions/67166083/density-curve-on-histogram-is-flat
d <- density(CpGs)
dx <- diff(d$x)[1]
sum(d$y)*dx # 1.000852
h <- hist(CpGs)
lines(x=d$x,y=max(h$counts)*d$y/dx)
