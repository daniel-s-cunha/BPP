library(RColorBrewer)
library(latex2exp)

compute_p_z <- function(i,k,n){
  #This n needs to be changed to include the total number of possible observations 365*20/16
  p_z = sapply(1:k, function(j) choose(i,j-1)*choose(n-i,k-j)/choose(n,k-1))
  p_z = rbind(c(1,rep(0,k-1)),p_z)
  return(p_z)
}

compute_p_z_zm1 <- function(n,k,p_z){
  # The columns of p_z_zm1 are the conditions on zm1.
  # Since for each condition there can only be the current and next state,
  # we only need to store the probability for the current state for each condition.
  # Then the probability of the next state is p(z=c+1 | zm1=c) = 1 - p(z=c | zm1=c)  
  #
  # Only n-1 transitions; p(z_1 = 1) = 1, THERE ARE MORE CONSTRAINTS FOR i<k THAT I'M NOT IMPOSING
  p_z_zm1     = array(dim=c(n-1,k)) 
  p_z_zm1[,1] = p_z[2:n,1]/p_z[1:(n-1),1]; p_z_zm1[n-1,1]=0;
  if(k>2){
  for(j in 2:(k-1)){
    p_z_zm1[,j] = (apply(p_z[2:n,1:j],1,"sum") - apply(p_z[1:(n-1),1:j-1,drop=F],1,"sum"))/p_z[1:(n-1),j]  
    }
  }
  p_z_zm1[(k+1):(n-1),k] = 1
  return(p_z_zm1)
}

k=10
n_disc = 25
t_disc = seq(1,n_disc,length.out=n_disc)
p_z_disc = compute_p_z(t_disc,k,n_disc)
p_z_zm1_disc = compute_p_z_zm1(n_disc+1,k,p_z_disc)
t_disc = c(0,t_disc)
#
pdf(file=paste0("figure3_v2.pdf"),family="ArialMT")
pal = brewer.pal(k, 'Dark2')
palgen = colorRampPalette(c(pal[3],pal[4]))#c("#8DA0CB","#FC8D62")) #c(pal[5], pal[9]))
pal = palgen(k)
par(cex.lab = 1.75, cex.axis = 1.2, cex.main = 1.75)
par(mar = c(5, 4.75, 4.75, 2))
plot(t_disc[2:(n_disc+1)],p_z_zm1_disc[,1],type='l',lty=2,ylim=c(0,1),col=pal[1],xlab='',ylab='Transition probability',lwd=1.25)
points(t_disc[2:(n_disc+1)],p_z_zm1_disc[,1],pch=15,cex=0.85)
for(j in 2:k){
  lines(t_disc[2:(n_disc+1)],p_z_zm1_disc[,j],col=pal[j],lty=2,lwd=1.25)
  points(t_disc[2:(n_disc+1)],p_z_zm1_disc[,j],col=pal[j],pch=15,cex=0.85)
}
axis(side=1,at=c(1))
title(TeX(r'(Discrete Self-Transitions $p(z_{i+1}=j|z_{i}=j)$)'), adj = 0.1, line = 1.,family="ArialMT") #' for $j=1,\ldots,k$)')
title(xlab='Time', adj = .5, line = 2.75,family="ArialMT") #' for $j=1,\ldots,k$)')
dev.off()

