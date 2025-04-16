library(RColorBrewer)
pal = brewer.pal(8, 'Set1')
source('~/Documents/BPP/src/functions.R')
locations = c('loving_county','rondonia','vineyards')
loc_names = c('Semi-arid Shrubland','Rondonia Deforestation','Crop Rotation')
loc_pids = c(6,16,3)
loc_Ks = c(6,6,6)
#
l = 1
setwd(paste0("~/Documents/BPP/case_study/",locations[l]))

#pixel_id
pid = loc_pids[l]

#Initialize data
df = read.csv(paste0('./',locations[l],'.csv'))
df$ndvi = (df$SR_B4 - df$SR_B3) / (df$SR_B4 + df$SR_B3)
#
#Permian Basin:
change_dates = list()
change_dates[[1]] = c()#use drought data
#
#Rondonia
change_dates[[2]] = c('12/1/2005')
#
#Vineyard
change_dates[[3]] = c('08/01/2006','08/01/2012','10/01/2016')

change_dates = change_dates[[l]]
change_dates = as.POSIXct(change_dates,format="%m/%d/%Y",origin='2000-01-01')
#
tau_true = data.frame(change_dates)
colnames(tau_true)= 'datetime'

j = unique(df$pixel_id)[pid]
h=2; psi=.1; lam=1; K=loc_Ks[l]; geometric=F 
#Permian Basin: pid=10,K=10
#Rondonia: pid=16,K=6
#CropRotation: pid=20,K=6

sr = na.omit(df[df$pixel_id==j,]); X = sr[,13,drop=F]; y = as.numeric(sr[,14])
X[,1] = as.POSIXct(X[,1], format="%Y-%m-%d")
colnames(X) = 'datetime'
tau_true = data.frame(c(X[1,1],tau_true[,1],X[dim(X)[1],1]))

colnames(tau_true) = 'datetime'
tt = design_setup(tau_true,1)[,2]
#
y = y[order(X[,1])]
X[,1] = X[order(X[,1]),]
quants=c(0.025,0.5,0.975)
re = fit_EM(X,y,K,h,psi,lam,quants=quants,nonInf_pk=F)
n = length(y)
res = re$res
taus.25 = res[[1]]$taus
taus = res[[2]]$taus
taus.75 = res[[3]]$taus
z = res[[2]]$zm
z_true = z_from_taus(n,floor(n*tt[-c(1,length(tt))]))
t = re$t
wl = whamming_loss(n,t,z,z_true)
lk = abs(max(z)-max(z_true))
theta = re$theta

mi = min(theta[3,]^2+theta[5,]^2);ma = max(theta[3,]^2+theta[5,]^2)
amp = (theta[3,]^2+theta[5,]^2 - mi)/(ma-mi)

mi = min(theta[1,]);ma = max(theta[1,])
int = (theta[1,] - mi)/(ma-mi)

mi = min(theta[2,]);ma = max(theta[2,])
slp = (theta[2,] - mi)/(ma-mi)

# pdf(paste0('./',locations[l],'_invPrior_impliedpk_pixel',pid,'.pdf'), width = 14, height = 5) #, width = 14, height = 5
# #bottom,left,top,right
# layout(matrix(c(1,1,2), 3, 1, byrow = T))
# par(cex.lab = 1.75, cex.axis = 1.1, cex.main = 1.85,mgp=c(2.25,1.,0),#oma = c(6.5, 3.0, 2.25, 2.25),
#   mai = c(0.05, 0.75, .65, 0.25))
# plot(t,y,xaxt='n',ylab='',xlab='',pch=20,col=adjustcolor('black',alpha.f=0.5))
# lines(re$t_pred,re$Ey_x)
# #abline(v=t[taus],col='black',lwd=2,lty=2) #1:length(res[[2]]$taus))
# abline(v = tt[-c(1,length(tt))],col=adjustcolor(pal[8],alpha.f=0.75),lwd=4,lty=2)
# title(paste0(loc_names[l],': Case Study Pixel'),line=.95,adj = 0.05)#0.275)
# title(ylab='NDVI',mgp=c(2.75,1.6,0),cex.lab=2.1) #cex.lab=2.2
# #legend('bottomright', legend = c('slp','amp','int'), col = pal[c(2,4,5)],lty=c(4,3,1),lwd=c(2.5,3.5,2))#"bottomright"
# #
# par(cex.lab = 1.75, cex.axis = 1.1, cex.main = 1.75,mgp=c(2.25,1.,0),
#   mai = c(.75, 0.75, 0.05, 0.25))
# plot(2,2,yaxt='n',xlim=c(0,1),ylim=c(0,1),ylab='',xlab='')#plot(t, res[[2]]$zm,type='l',lwd=2,col='gray') #zms
# #abline(v=t[taus.25],col='black',lwd=.9,lty=2) #1:length(res[[1]]$taus))
# abline(v=t[taus],col='black',lwd=4,lty=2) #1:length(res[[2]]$taus))
# #abline(v=t[taus.75],col='black',lwd=.9,lty=2) #1:length(res[[3]]$taus))
# for(j in 2:(length(taus)+2)){
# lines(t[c(1,taus,n)[c(j-1,j)]],c(slp[j-1],slp[j-1]),col=pal[2],lty=4,lwd=2.5)
# lines(t[c(1,taus,n)[c(j-1,j)]],c(amp[j-1],amp[j-1]),col=pal[4],lty=3,lwd=3.5)
# lines(t[c(1,taus,n)[c(j-1,j)]],c(int[j-1],int[j-1]),col=pal[5],lty=1,lwd=2)
# }
# #
# # if(!is.null(taus.25)){
# #   rect(xleft = t[taus.25], xright = t[taus.75], ybottom = 0, ytop = K+1, 
# #     border = NA, col = c(adjustcolor("blue", alpha = 0.2),adjustcolor("purple", alpha = 0.2)))
# # }
# title(ylab='Detected\n changes',mgp=c(1.2,1.6,0),cex.lab=2.1) #cex.lab=2.2
# title(xlab="Time",mgp=c(2.6,1.6,0),cex.lab=1.85) #cex.lab=2
# legend('topright', legend = c('slp','amp','int'), col = pal[c(2,4,5)],lty=c(4,3,1),lwd=c(2.5,3.5,2))#"bottomright"
# dev.off()


#Plot for loving county
pdf(paste0('./',locations[l],'_invPrior_Drought_impliedpk_pixel_',pid,'.pdf'), width = 14, height = 7.5) #, width = 14, height = 5
#bottom,left,top,right
layout(matrix(c(1,1,2,3,3), nrow=5, ncol=1, byrow = T))
par(cex.lab = 1.75, cex.axis = 1.1, cex.main = 1.85,mgp=c(2.25,1.,0),#oma = c(6.5, 3.0, 2.25, 2.25),
  mai = c(0.05, 0.75, .65, 0.25))
plot(t,y,xaxt='n',ylab='',xlab='',pch=20,col=adjustcolor('black',alpha.f=0.5))
lines(re$t_pred,re$Ey_x)
#abline(v = tt[-c(1,length(tt))],col=pal[1],lwd=2,lty=2)
title(paste0(loc_names[l],': Case Study Pixel'),line=.95,adj = 0.05)#0.275)
title(ylab='NDVI',mgp=c(2.75,1.6,0),cex.lab=2.1) #cex.lab=2.2
#
par(cex.lab = 1.75, cex.axis = 1.1, cex.main = 1.75,mgp=c(2.25,1.,0),
  mai = c(.05, 0.75, 0.05, 0.25))
plot(2,2,yaxt='n',xaxt='n',xlim=c(0,1),ylim=c(0,1),ylab='',xlab='')#plot(t, res[[2]]$zm,type='l',lwd=2,col='gray') #zms
abline(v=t[taus],col='black',lwd=4,lty=2) #1:length(res[[2]]$taus))
for(j in 2:(length(taus)+2)){
  lines(t[c(1,taus,n)[c(j-1,j)]],c(slp[j-1],slp[j-1]),col=pal[2],lty=4,lwd=2.5)
  lines(t[c(1,taus,n)[c(j-1,j)]],c(amp[j-1],amp[j-1]),col=pal[4],lty=3,lwd=3.5)
  lines(t[c(1,taus,n)[c(j-1,j)]],c(int[j-1],int[j-1]),col=pal[5],lty=1,lwd=2)
}
title(ylab='Detected\n changes',mgp=c(1.2,1.6,0),cex.lab=2.1) #cex.lab=2.2
#title(xlab="Time",mgp=c(2.6,1.6,0),cex.lab=1.85) #cex.lab=2
legend('bottomright', legend = c('slp','amp','int'), col = pal[c(2,4,5)],lty=c(4,3,1),lwd=c(2.5,3.5,2))#"bottomright"
#
dr = read.csv('./drought.csv')
dr$date = as.POSIXct(dr$Week,format="%Y-%m-%d",origin='2000-01-01')
dr = dr[dr$date<as.POSIXct('2023-01-01',format="%Y-%m-%d",origin='2000-01-01'),]
dr_dt = data.frame(dr$date)
colnames(dr_dt) = 'datetime'
dr$time = design_setup(dr_dt,1)[,2]
#
par(cex.lab = 1.75, cex.axis = 1.1, cex.main = 1.75,mgp=c(2.25,1.,0),
  mai = c(.75, 0.75, 0.05, 0.25))
plot(dr$time,dr$None,ylab='',xlab='',yaxt='n',pch=4,col=adjustcolor('black',alpha.f=0.4))
axis(2,c(0,100))
title(ylab='Absence of\n Drought',mgp=c(1.2,1.6,0),cex.lab=2.1) #cex.lab=2.2
title(xlab="Time",mgp=c(2.6,1.6,0),cex.lab=1.85) #cex.lab=2
dev.off()
