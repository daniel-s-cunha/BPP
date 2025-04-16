library(RColorBrewer)
pal = brewer.pal(8, 'Dark2')
#pal = brewer.pal(8, 'RdBu')
library(latex2exp)
source('~/Documents/BPP/src/functions.R')
setwd("~/Documents/BPP/plots")

pdf('fig_prior_on_number_segments_v3.pdf') #, width = 14, height = 5
set.seed(13)
n=20
K = 4
t = sort(runif(n))
logp = c()
ti = t[2:n]; tim1 = t[1:(n-1)]
for(k in 1:K){
    logp = c(logp, -k * sum(log( (1-ti) / (1-tim1) )) - (k/2)*log(2*pi))
    #logp = c(logp, n*sum(sapply(1:k, function(j) lfactorial(j))) + k^2*sum(log( (1-ti) / (1-tim1) )) - k*(k+1)/2 * sum(log( (ti - tim1)/(1 - tim1) )))
}

par(cex.lab = 1.7, cex.axis = 1.1, cex.main = 1.8,mgp=c(2.25,1.,0),#oma = c(6.5, 3.0, 2.25, 2.25),
  mai = c(1., 1.05, .75, 0.6)) #bottom,left,top,right
plot(1:K,exp(logp - lse(logp)),type='l',col=adjustcolor(pal[8], alpha = 0.2),ylab='',xlab='',xaxt='n',ylim=c(0,1))
axis(1, at=1:K)
title('Comparing priors on number of segments',line=1.,adj = 0.75)#0.275)
title(ylab=TeX(r'($\pi(k)$)'),mgp=c(2.5,1.6,0),cex.lab=2.) #cex.lab=2.2
title(xlab='Number of segments k',mgp=c(2.85,1.6,0),cex.lab=1.75) #cex.lab=2.2
lines(1:K,exp(-logp - lse(-logp)),col=adjustcolor(pal[4], alpha = 0.2))

for(i in 1:2000){
    t = sort(runif(n))
    logp = c()
    ti = t[2:n]; tim1 = t[1:(n-1)]
    for(k in 1:K){
        logp = c(logp, -k * sum(log( (1-ti) / (1-tim1) )) - (k/2)*log(2*pi))
        #logp = c(logp, n*sum(sapply(1:k, function(j) lfactorial(j))) + k^2*sum(log( (1-ti) / (1-tim1) )) - k*(k+1)/2 * sum(log( (ti - tim1)/(1 - tim1) )))
    }
    lines(1:K,exp(logp - lse(logp)),type='l',col=adjustcolor(pal[8], alpha = 0.2))
    lines(1:K,exp(-logp - lse(-logp)),col=adjustcolor(pal[4], alpha = 0.2))
    cat(exp(logp - lse(logp)),'\n')
}
dev.off()



