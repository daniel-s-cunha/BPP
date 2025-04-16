# BPP
This repository is associated with the paper, *Bernstein Polynomial Processes for Continuous Change Time Change Detection*. Here is example code to run the **BPP** model using an intercept only likelihood model,
```
#Set working directory to BPP first
setwd('~/Documents/BPP')
source('./src/functions.R')
set.seed(7)
#
n = 500; k_star = 3; #k_star is the true number of segments
soc = 0.5 #size of true changes
#
ti = c(0,sort(rbeta(n-2,0.5,0.5)),1) #continuous time points
chpts = sample_BPP(ti,n,k_star) #sample change points from a Bernstein polynomial process
#
y = c(); prev=0; cnt = 1
for(j in c(chpts,n)){
  y = c(y, sqrt(0.1)*rt(j-prev,df=3) + cnt*soc)
  prev = j
  cnt = cnt + 1
}
# Fit BPP model with intercept only t-distributed likelihood with nu=3 degrees of freedom
res = fit_EM(ti,y,K=5,intercept=T)
taus = res$res[[2]]$taus
#
#Plot data and true changes
plot(ti,y);abline(v=ti[chpts],col='orange',lwd=2);abline(v=ti[taus],col='purple',lty=2,lwd=2)
```
Here is example code to run the Interannually Varying Harmonics model on remote sensing data from a deforestation example in the Amazon rainforest,
```
#Set working directory first
setwd('~/Documents/BPP')
source('./src/functions.R')
#
df = read.csv('./case_study/rondonia/rondonia.csv')
df$ndvi = (df$SR_B4 - df$SR_B3) / (df$SR_B4 + df$SR_B3)
#
# Setup data to pass to the fit function:
df = na.omit(df); X = df[,13,drop=F]; y = as.numeric(df[,14])
X[,1] = as.POSIXct(X[,1], format="%Y-%m-%d")
colnames(X) = 'datetime'
y = y[order(X[,1])]
X[,1] = X[order(X[,1]),]
#
# Fit BPP model with intercept only t-distributed likelihood with nu=3 degrees of freedom
res = fit_EM(X,y,K=5,intercept=F)
taus = res$res[[2]]$taus
#
#Plot data and true changes
plot(X[,1],y);abline(v=X[taus,1],col='purple',lty=2,lwd=2)
```
