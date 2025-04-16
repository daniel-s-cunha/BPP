library(RColorBrewer)
library(latex2exp)

compute_ctmc<-function(n,k,t){
  ctmc = array(dim=c(n-1,k,k))
  for(j in 1:k){
    for(h in j:k){
      ctmc[,j,h] = choose(k-j,h-j)*((t[2:n]-t[1:(n-1)])/(1-t[1:(n-1)]))^(h-j) * ((1-t[2:n])/(1-t[1:(n-1)]))^(k-h) 
    }
  }
  return(ctmc)
}

compute_p_z <- function(i,k,n){
  #This n needs to be changed to include the total number of possible observations 365*20/16
  p_z = sapply(1:k, function(j) choose(i,j-1)*choose(n-i,k-j)/choose(n,k-1))
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

compute_p_z_zm1_cont <- function(n,k,t){
  # The columns of p_z_zm1 are the conditions on zm1.
  # Since for each condition there can only be the current and next state,
  # we only need to store the probability for the current state for each condition.
  # Then the probability of the next state is p(z=c+1 | zm1=c) = 1 - p(z=c | zm1=c)  
  #
  # Only n-1 transitions; p(z_1 = 1) = 1,
  p_z_zm1     = array(dim=c(n-1,k)) 
  for(j in 1:k){
    p_z_zm1[,j] = ((1-t[2:n])/(1-t[1:(n-1)]))^(k-j) #/ (((1-t[2:n])/(1-t[1:(n-1)]))^(k-j) + (k-j)*(t[2:n]-t[1:(n-1)])/(1-t[1:(n-1)]) * ((1-t[2:n])/(1-t[1:(n-1)]))^(k-j-1));
    #if(j>1){p_z_zm1[1:(j-1),j] = NA}
    #p_z_zm1[(n-1-(k-j)):(n-1),j] = 0
  }
  p_z_zm1[(n-1),] = 0
  return(p_z_zm1)
}

# # Compare Discrete and Continuous on the same times
# k=10
# n_disc = 25
# t_disc = seq(1,n_disc,length.out=n_disc)
# t_cont = t_disc/n_disc
# p_z_disc = compute_p_z(t_disc,k,n_disc)
# p_z_zm1_disc = compute_p_z_zm1(n_disc,k,p_z_disc)
# p_z_zm1_cont = compute_p_z_zm1_cont(n_disc,k,(t_disc-1)/(n_disc-1))
# #Use Set2: color3 to color2
# #pdf(file=paste0("figure4_v1.pdf"),family="ArialMT")
# #
# pal = brewer.pal(k, 'Dark2')
# palgen = colorRampPalette(c(pal[3],pal[4]))#c("#8DA0CB","#FC8D62")) #c(pal[5], pal[9]))
# pal1 = palgen(k)
# palgen = colorRampPalette(c(pal[5],pal[6]))#c("#8DA0CB","#FC8D62")) #c(pal[5], pal[9]))
# pal2 = palgen(k)
# #
# par(cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.2)
# par(mar = c(5, 4.75, 4.75, 2))
# plot(t_disc[2:n_disc],p_z_zm1_disc[,1],type='l',ylim=c(0,1),col=pal1[1],xlab='Time',ylab='Transition probability',lwd=4.5)
# for(j in 2:k) lines(t_disc[2:n_disc],p_z_zm1_disc[,j],col=pal1[j],lwd=4.5)
# lines(t_disc[2:n_disc],p_z_zm1_cont[,1],col=pal2[1],lwd=4.5)
# for(j in 2:k) lines(t_disc[2:n_disc],p_z_zm1_cont[,j],col=pal2[j],lwd=4.5)

pal = brewer.pal(8, 'Dark2')
palgen = colorRampPalette(c(pal[8],pal[3]))#c("#8DA0CB","#FC8D62")) #c(pal[5], pal[9]))
pal2 = palgen(k)
#
eff_seed <- sample(1:2^15, 1)
print(sprintf("Seed for session: %s", eff_seed))
set.seed(eff_seed)
set.seed(27851)
#
k=5
n_disc = 25
t_cont = c(0,sort(runif(n_disc-2)),1)
# t_disc = seq(1,n_disc,length.out=n_disc)
#
# Show k-1 normalized(self, or next) self-transition probabilities on non-uniform time 
pdf(file=paste0("figure4_v2.pdf"),family="ArialMT")
par(cex.lab = 1.75, cex.axis = 1.2, cex.main = 1.75)
par(mar = c(5, 4.75, 4.75, 2))
p_z_zm1_cont = compute_p_z_zm1_cont(n_disc,k,t_cont)
plot(t_cont[2:n_disc],p_z_zm1_cont[,1],type='n',ylim=c(0,1),col=pal[4],xaxt = "n",
      xlab='',ylab='Transition probability')
title(TeX(r'($p(z_{t(i+1)}=j|z_{t(i)}=j)$ for $j=1,\ldots,k-1$)'), adj = 0.1, line = 1.)
title(xlab="Time", line=1.75, cex.lab=1.75)
axis(side = 1, labels=F, at = t_cont,lwd=0.25,col='black')
axis(side = 1, labels=T, at = c(0.0,1.0))
for(i in t_cont){abline(v=i,col='lightgrey',lwd=0.75)}
lines(t_cont[2:n_disc],p_z_zm1_cont[,1],col=pal[4],lwd=1.5)
for(j in 2:(k-1)) lines(t_cont[2:n_disc],p_z_zm1_cont[,j],col=pal[8],lwd=1.5)
dev.off()

pdf(file=paste0("figure5_v2.pdf"),family="ArialMT")
par(cex.lab = 1.75, cex.axis = 1.2, cex.main = 1.75)
par(mar = c(5, 4.75, 4.75, 2))
pal = brewer.pal(8, 'Dark2')
palgen = colorRampPalette(c(pal[8],pal[1]))#c("#8DA0CB","#FC8D62")) #c(pal[5], pal[9]))
pal2 = palgen(k)
pal3 = c(pal[1],pal[6],pal[2],pal[7])
# Show probability of transitioning from 1 to 1:k
#pdf(file=paste0("figure4_v1.pdf"),family="ArialMT")
ctmc = compute_ctmc(n_disc,k,t_cont)
plot(t_cont[2:n_disc],ctmc[,1,1],type='n',ylim=c(0,1),col=pal[4],xaxt = "n",
      xlab='',ylab='Transition probability',lwd=3)
for(i in t_cont){abline(v=i,col='lightgray',lwd=0.75)}
title(TeX(r'($p(z_{t(i+1)}=h|z_{t(i)}=1)$ for $h=1,\ldots,k$)'), adj = 0.1, line = 1.)
title(xlab="Time", line=1.75, cex.lab=1.75)
axis(side = 1, labels=F, at = t_cont,lwd=0.25,col='black')
axis(side = 1, labels=T, at = c(0.0,1.0))
lines(t_cont[2:n_disc],ctmc[,1,1],col=pal[4],lwd=1.5)
for(j in 2:k) lines(t_cont[2:n_disc]+rnorm(n-1,0,0.001),ctmc[,1,j],col=pal3[j-1],lwd=1.5)
dev.off()
#interesting seeds:
# 27851





