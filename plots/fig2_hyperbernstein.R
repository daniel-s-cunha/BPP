library(RColorBrewer)
library(latex2exp)

bernstein<-function(j,k,t){
 choose(k-1,j-1)*t^(j-1)*(1-t)^(k-j)
}

compute_p_z_cont <- function(t,k){
  n = length(t)
  p_z = array(0,dim=c(n,k))
  p_z[,1] = bernstein(1,k,t)
  for(j in 2:k){
    p_z[,j] = bernstein(j,k,t)
  }
  return(p_z)
}

compute_p_z_disc <- function(t,k,n_disc){
  #This n needs to be changed to include the total number of possible observations 365*20/16
  n = n_disc#floor(365*20/16) #length(t) 
  is = floor(t*n)
  p_z = sapply(1:k, function(j) choose(is,j-1)*choose(n-is,k-j)/choose(n,k-1))
  return(p_z)
}

k=4
n_disc = 10
t_disc = seq(1,n_disc,length.out=n_disc)/n_disc
p_z_disc = compute_p_z_disc(t_disc,k,n_disc)
t_disc = c(0,t_disc)
p_z_disc = rbind(c(1,rep(0,k-1)),p_z_disc)

n_cont = 1000
t_cont = seq(0,1,length.out=n_cont)
p_z_cont = compute_p_z_cont(t_cont,k)

pdf(file=paste0("figure2_v2.pdf"),family="ArialMT")
pal = brewer.pal(k, 'Dark2')
par(cex.lab = 1.75, cex.axis = 1.2, cex.main = 1.75)
par(mar = c(5, 4.75, 4.75, 2))
# par(cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.2)
# par(mar = c(5, 4.75, 4.75, 2))
plot(t_cont,p_z_cont[,1],type='l',col=pal[1],xlab='',xaxt='n',ylab='Marginal probability',lwd=1.7,family="ArialMT"); 
for(j in 2:k) lines(t_cont,p_z_cont[,j],col=pal[j],lwd=1.7)
title(TeX(r'(Discrete and Continuous Marginals $p(z_{t(i)}=j)$)'), adj = 0.1, line = 1.,family="ArialMT") #' for $j=1,\ldots,k$)')
title(xlab='Time',line=2.75,family="ArialMT")
axis(side=1,at=c(0,0.3,0.7,1.0))#seq(0,10,2)/10
axis(side=1,at=seq(0,10,1)/10,labels=F)
#
#lines(c(0,0), c(0, 1), type = "l", lwd = 1,col=pal[1])  # Line segments
pchi = 15
points(0, 1, pch = pchi, cex = .75,col=pal[1])  # Circles at the top
lines(c(0,t_disc), c(1,p_z_disc[,1]), type = "l", col=pal[1],lty=2,lwd=1.25)  # Line segments
for(j in 2:k) lines(t_disc, p_z_disc[,j], type = "l",col=pal[j],lty=2,lwd=1.25)  # Line segments
for(j in 1:k){
  for (i in 1:(n_disc+1)) {
    if(!((j==k & i<=2)|(j==(k-1) & i<=1))){
      points(t_disc[i], p_z_disc[i,j], pch = pchi, cex = .75,col=pal[j])  # Circles at the top
    }
  }
}
dev.off()
