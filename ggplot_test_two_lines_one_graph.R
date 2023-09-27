library(tidyverse)
Time = 50
N =20
K =10
phi = 0.5

y = matrix(rep(0, N*Time), nrow = Time)

for (t in 2:Time){
  delta = rnorm(1)
  epsilons = rnorm(K)
  lam = phi*y[(t-1),] + rep(epsilons, each=2)
  noise = rnorm(20)
  y[t,] = delta + lam + noise
}


df <- as.data.frame(y)
df<-mutate(df, "time"= 1:Time)
melt(df, id="time") %>% ggplot(aes(x=time,y=value)) + geom_line(aes(color=variable))


library(reshape2) # install.packages("reshape2")

R = data.frame(Delta = c(1,2), UE = c(1,1), RE = c(3.8, 2.4))
meltR = melt(R, id = "Delta")

ggplot(meltR, aes(x = Delta, y = value, group = variable, colour = variable)) + 
  geom_line() 
