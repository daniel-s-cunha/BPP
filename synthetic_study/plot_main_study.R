library(latex2exp)
library(MASS)
library(RColorBrewer)
setwd('~/Documents/BPP/synthetic_study/')
source('../src/functions.R')

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Main plot for synthetic study
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


#Discrete time, normal, and BPP_impliedpk
res12 = read.csv('main_study_results_v1.csv')
res12[] <- lapply(res12, function(x) as.numeric(as.character(x)))
res12 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res12,FUN=sum)

res12$BPP_or = res12$BPP_omission/res12$len_tau_true
res12$BPP_cr = res12$BPP_commission/res12$BPP_len_tau

res12$BPP_impliedpk_or = res12$BPP_impliedpk_omission/res12$len_tau_true
res12$BPP_impliedpk_cr = res12$BPP_impliedpk_commission/res12$BPP_impliedpk_len_tau

res12$BinSeg_or = res12$BinSeg_omission/res12$len_tau_true
res12$BinSeg_cr = res12$BinSeg_commission/res12$BinSeg_len_tau

res12$PELT_or = res12$PELT_omission/res12$len_tau_true
res12$PELT_cr = res12$PELT_commission/res12$PELT_len_tau

res12[] <- lapply(res12, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}


res11 = res12[res12$robustness==3, ] #res1$time_dist==1 & 
#Main plot:
pal = brewer.pal(8, 'RdBu')#'PuOr'
pdf('main_study_full_nu3_v3.pdf', width = 12, height = 8)
layout(matrix(1:6, nrow = 2, ncol = 3,byrow=T))
#bottom,left,top,right
par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.675, 0.65, .6, 0.4))
x = res12$BPP_cr; y = 1-res12$BPP_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BPP: all data',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

x = res12$PELT_cr; y = 1-res12$PELT_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
	pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('PELT: all data',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

x = res12$BinSeg_cr; y = 1-res12$BinSeg_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
	pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BinSeg: all data',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

x = res11$BPP_cr; y = 1-res11$BPP_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BPP: nu=3 data',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

x = res11$PELT_cr; y = 1-res11$PELT_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
	pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('PELT: nu=3 data',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

x = res11$BinSeg_cr; y = 1-res11$BinSeg_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
	pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BinSeg: nu=3 data',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")


dev.off()



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Subset plots for appendix
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

res = read.csv('main_study_results_v1.csv')
res[] <- lapply(res, function(x) as.numeric(as.character(x)))

res1 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res,FUN=sum)

res1$BPP_or = res1$BPP_omission/res1$len_tau_true
res1$BPP_cr = res1$BPP_commission/res1$BPP_len_tau

res1$BPP_impliedpk_or = res1$BPP_impliedpk_omission/res1$len_tau_true
res1$BPP_impliedpk_cr = res1$BPP_impliedpk_commission/res1$BPP_impliedpk_len_tau

res1$BinSeg_or = res1$BinSeg_omission/res1$len_tau_true
res1$BinSeg_cr = res1$BinSeg_commission/res1$BinSeg_len_tau

res1$PELT_or = res1$PELT_omission/res1$len_tau_true
res1$PELT_cr = res1$PELT_commission/res1$PELT_len_tau

res1[] <- lapply(res1, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}

#Subplots
ev = c(0.1,0.2,0.3); rb = c(3,10,100)
res1_list = list()
for(j in 1:3) res1_list[[j]] = res1[res1$time_dist==j, ]
for(j in 1:3) res1_list[[j+3]] = res1[res1$error_var==ev[j], ]
for(j in 1:3) res1_list[[j+6]] = res1[res1$robustness==rb[j], ]
pal = brewer.pal(8, 'RdBu')#'PuOr'

for(ll in 1:3){
	pdf(paste0('main_study_subset',ll,'_v1.pdf'), width = 12, height = 12)
	layout(matrix(1:9, nrow = 3, ncol = 3,byrow=T))
	#bottom,left,top,right
	par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.75, 0.65, .5, 0.35))
	for(lll in 1:3){
		res11 = res1_list[[lll + 3*(ll-1)]]
		#	
		x = res11$BPP_cr; y = 1-res11$BPP_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
				pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('BPP',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		print(z)
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

		x = res11$PELT_cr; y = 1-res11$PELT_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('PELT',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		print(z)
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

		x = res11$BinSeg_cr; y = 1-res11$BinSeg_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('BinSeg',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		print(z)
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")
	}
	dev.off()
}


# Add subplot for k=2,3,4


# NUMBER OF SEGMENTS SUBPLOT
pal = brewer.pal(8, 'RdBu')#'PuOr'
pdf(paste0('main_study_numSeg_v0.pdf'), width = 12, height = 12)
layout(matrix(1:12, nrow = 4, ncol = 3,byrow=T))
#bottom,left,top,right
for(ll in 1:4){
	par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.75, 0.65, .5, 0.35))
	res11 = res1[res1$num_regimes==ll, ];
	#	
	x = res11$BPP_cr; y = 1-res11$BPP_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('BPP',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	print(z)
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

	x = res11$PELT_cr; y = 1-res11$PELT_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('PELT',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	print(z)
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

	x = res11$BinSeg_cr; y = 1-res11$BinSeg_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('BinSeg',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	print(z)
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")
}
dev.off()










####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Violin plot for main synthetic study
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

library(vioplot)
f1 <- function(cr,or){
	2*(1-cr)*(1-or)/(2-cr-or)
}

pal = brewer.pal(10, 'Paired')#'PuOr'
pdf('main_study_ComOmF1_v3.pdf', width = 12, height = 4)
layout(matrix(1:3, nrow = 1, ncol = 3,byrow=T))
#bottom,left,top,right
par(cex.lab = 2., cex.axis = 1.7, cex.main = 2.25, mgp=c(2.75,1.,0),mai = c(.55, 0.65, .5, 0.375))

data_BPP = data.frame(value = rep(res12$BPP_cr, 2),
					group = factor(c(rep(c('BPP vs PELT','BPP vs BinSeg'), each = 648))))
data_BPP$group <- factor(data_BPP$group, levels = unique(data_BPP$group))
data_PBin = data.frame(value = c(res12$PELT_cr,res12$BinSeg_cr),
					group = factor(c(rep(c('BPP vs PELT','BPP vs BinSeg'), each = 648))))
data_PBin$group <- factor(data_PBin$group, levels = unique(data_PBin$group))

histoplot(value~group, data=data_BPP, col = pal[2], plotCentre = "line", side = "left",xlab='',ylab='',yaxt='n')
histoplot(value~group, data=data_PBin, col = pal[c(6,5)], plotCentre = "line", side = "right", add = T)
title(xlab = "", ylab = "Commission Rate")
axis(side = 2, cex.axis = 1.5)
title('Commission Rate',line=.5,adj = 0.09)


data_BPP = data.frame(value = rep(res12$BPP_or, 2),
					group = factor(c(rep(c('BPP vs PELT','BPP vs BinSeg'), each = 648))))
data_BPP$group <- factor(data_BPP$group, levels = unique(data_BPP$group))
data_PBin = data.frame(value = c(res12$PELT_or,res12$BinSeg_or),
					group = factor(c(rep(c('BPP vs PELT','BPP vs BinSeg'), each = 648))))
data_PBin$group <- factor(data_PBin$group, levels = unique(data_PBin$group))

histoplot(value~group, data=data_BPP, col = pal[2], plotCentre = "line", side = "left",xlab='',ylab='',yaxt='n')
histoplot(value~group, data=data_PBin, col = pal[c(6,5)], plotCentre = "line", side = "right", add = T)
title(xlab = "", ylab = "Omission Rate")
axis(side = 2, cex.axis = 1.5)
title('Omission Rate',line=.5,adj = 0.09)


data_BPP = data.frame(value = rep(f1(res12$BPP_cr,res12$BPP_or), 2),
					group = factor(c(rep(c('BPP vs PELT','BPP vs BinSeg'), each = 648))))
data_BPP$group <- factor(data_BPP$group, levels = unique(data_BPP$group))
data_PBin = data.frame(value = c(f1(res12$PELT_cr,res12$PELT_or),f1(res12$BinSeg_cr,res12$BinSeg_or)),
					group = factor(c(rep(c('BPP vs PELT','BPP vs BinSeg'), each = 648))))
data_PBin$group <- factor(data_PBin$group, levels = unique(data_PBin$group))

histoplot(value~group, data=data_BPP, col = pal[2], plotCentre = "line", side = "left",xlab='',ylab='',yaxt='n')
histoplot(value~group, data=data_PBin, col = pal[c(6,5)], plotCentre = "line", side = "right", add = T)
title(xlab = "", ylab = "F1 Score")
axis(side = 2, cex.axis = 1.5)
title('F1 Score',line=.5,adj = 0.09)

dev.off()









####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Gibbs sampler plot
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################



#gibbs sampler plot
res = read.csv('main_study_results_gibbs_v2_nsamp500.csv')
res[] <- lapply(res, function(x) as.numeric(as.character(x)))

res1 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res,FUN=sum)

res1$BPP_or = res1$BPP_omission/res1$len_tau_true
res1$BPP_cr = res1$BPP_commission/res1$BPP_len_tau

res1[] <- lapply(res1, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}

#Main plot:
res11 = res1 #[res1$size_change>0.5, ] #res1$time_dist==1 & 
pal = brewer.pal(8, 'RdBu')##'PuOr'
pdf('main_study_gibbs_v3_nsamp500.pdf', width = 4, height = 4)
#layout(matrix(1, nrow = 1, ncol = 1,byrow=T))
#bottom,left,top,right
par(cex.lab = 1.25, cex.axis = 1., cex.main = 1.25,mgp=c(1.25,1.,0),mai = c(.75, 0.65, .5, 0.35))
x = res11$BPP_cr; y = 1-res11$BPP_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BPP: Gibbs sampler',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)
dev.off()



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# NORMAL DISCRETE BPPE plot for synthetic study
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################


#Discrete time, normal, and BPP_impliedpk
#Discrete time, normal, and BPP_impliedpk
res2 = read.csv('main_study_results_v1.csv')
res2[] <- lapply(res2, function(x) as.numeric(as.character(x)))
res2 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res2,FUN=sum)

res2$BPP_or = res2$BPP_omission/res2$len_tau_true
res2$BPP_cr = res2$BPP_commission/res2$BPP_len_tau

res2$BPP_impliedpk_or = res2$BPP_impliedpk_omission/res2$len_tau_true
res2$BPP_impliedpk_cr = res2$BPP_impliedpk_commission/res2$BPP_impliedpk_len_tau

res2$BinSeg_or = res2$BinSeg_omission/res2$len_tau_true
res2$BinSeg_cr = res2$BinSeg_commission/res2$BinSeg_len_tau

res2$PELT_or = res2$PELT_omission/res2$len_tau_true
res2$PELT_cr = res2$PELT_commission/res2$PELT_len_tau

res2[] <- lapply(res2, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})


res = read.csv('main_study_results_norm_disc_v2.csv')
res[] <- lapply(res, function(x) as.numeric(as.character(x)))
res1 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res,FUN=sum)

res1$Norm_or = res1$Norm_omission/res1$len_tau_true
res1$Norm_cr = res1$Norm_commission/res1$Norm_len_tau

res1$Disc_or = res1$Disc_omission/res1$len_tau_true
res1$Disc_cr = res1$Disc_commission/res1$Disc_len_tau

res1[] <- lapply(res1, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}

#Main plot:
pal = brewer.pal(8, 'RdBu')#'PuOr'
pdf('main_study_NDBPPE_v0.pdf', width = 12, height = 4)
layout(matrix(1:3, nrow = 1, ncol = 3,byrow=T))
#bottom,left,top,right
par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.675, 0.65, .6, 0.4))
x = res1$Norm_cr; y = 1-res1$Norm_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BPP: Normal Likelihood',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

x = res1$Disc_cr; y = 1-res1$Disc_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
	pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('Noninf Discrete',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

x = res2$BPP_impliedpk_cr; y = 1-res2$BPP_impliedpk_or;
plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
	pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
axis(1,c(0,1))
axis(2,c(0,1))
#axis(1, c(0.5) ,labels = FALSE)
axis(2, c(0.5) ,labels = FALSE)
title('BPPE',line=.5,adj = 0.09)#0.275)
abline(0,1)
z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

dev.off()



####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# NORMAL DISCRETE BPPE Subset plots for appendix
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#Discrete time, normal, and BPP_impliedpk
#Discrete time, normal, and BPP_impliedpk
res2 = read.csv('main_study_results_v1.csv')
res2[] <- lapply(res2, function(x) as.numeric(as.character(x)))
res2 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res2,FUN=sum)

res2$BPP_or = res2$BPP_omission/res2$len_tau_true
res2$BPP_cr = res2$BPP_commission/res2$BPP_len_tau

res2$BPP_impliedpk_or = res2$BPP_impliedpk_omission/res2$len_tau_true
res2$BPP_impliedpk_cr = res2$BPP_impliedpk_commission/res2$BPP_impliedpk_len_tau

res2$BinSeg_or = res2$BinSeg_omission/res2$len_tau_true
res2$BinSeg_cr = res2$BinSeg_commission/res2$BinSeg_len_tau

res2$PELT_or = res2$PELT_omission/res2$len_tau_true
res2$PELT_cr = res2$PELT_commission/res2$PELT_len_tau

res2[] <- lapply(res2, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})


res = read.csv('main_study_results_norm_disc_v2.csv')
res[] <- lapply(res, function(x) as.numeric(as.character(x)))
res1 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res,FUN=sum)

res1$Norm_or = res1$Norm_omission/res1$len_tau_true
res1$Norm_cr = res1$Norm_commission/res1$Norm_len_tau

res1$Disc_or = res1$Disc_omission/res1$len_tau_true
res1$Disc_cr = res1$Disc_commission/res1$Disc_len_tau

res1[] <- lapply(res1, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}

#Subplots
ev = c(0.1,0.2,0.3); rb = c(3,10,100)
res1_list = list();res2_list = list()
for(j in 1:3) res1_list[[j]] = res1[res1$time_dist==j, ]
for(j in 1:3) res1_list[[j+3]] = res1[res1$error_var==ev[j], ]
for(j in 1:3) res1_list[[j+6]] = res1[res1$robustness==rb[j], ]
#
for(j in 1:3) res2_list[[j]] = res2[res2$time_dist==j, ]
for(j in 1:3) res2_list[[j+3]] = res2[res2$error_var==ev[j], ]
for(j in 1:3) res2_list[[j+6]] = res2[res2$robustness==rb[j], ]

pal = brewer.pal(8, 'RdBu')#'PuOr'

for(ll in 1:3){
	pdf(paste0('main_study_NDBPPE_subset',ll,'_v0.pdf'), width = 12, height = 12)
	layout(matrix(1:9, nrow = 3, ncol = 3,byrow=T))
	#bottom,left,top,right
	par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.75, 0.65, .5, 0.35))
	for(lll in 1:3){
		res11 = res1_list[[lll + 3*(ll-1)]]; res12 = res2_list[[lll + 3*(ll-1)]]
		#	
		x = res11$Norm_cr; y = 1-res11$Norm_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
				pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('BPP: Normal Likelihood',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		print(z)
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

		x = res11$Disc_cr; y = 1-res11$Disc_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('Noninf Discrete',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		print(z)
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

		x = res12$BPP_impliedpk_cr; y = 1-res12$BPP_impliedpk_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('BPPE',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		print(z)
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")
	}
	dev.off()
}


# NUMBER OF SEGMENTS SUBPLOT
pal = brewer.pal(8, 'RdBu')#'PuOr'
pdf(paste0('main_study_NDBPPE_numSeg_v0.pdf'), width = 12, height = 12)
layout(matrix(1:12, nrow = 4, ncol = 3,byrow=T))
#bottom,left,top,right
for(ll in 1:4){
	par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.75, 0.65, .5, 0.35))
	res11 = res1[res1$num_regimes==ll, ]; res12 = res2[res2$num_regimes==ll, ]
	#	
	x = res11$Norm_cr; y = 1-res11$Norm_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('BPP: Normal Likelihood',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	print(z)
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

	x = res11$Disc_cr; y = 1-res11$Disc_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('Noninf Discrete',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	print(z)
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")

	x = res12$BPP_impliedpk_cr; y = 1-res12$BPP_impliedpk_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('BPPE',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	print(z)
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")
}
dev.off()

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Subset plots for Discrete, Normal, BPPE appendix
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Plot 1
# Main 3

# Plots 2:
# Time dist breakdown
# Robustness breakdown
# Error var breakdown
# Plot 2
# time dist 3

# Plot 3
# k=2,3,4

setwd('~/Documents/BPP/synthetic_study/')
#Discrete time, normal, and BPP_impliedpk
res2 = read.csv('main_study_results_v1.csv')
res2[] <- lapply(res2, function(x) as.numeric(as.character(x)))
res2 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res2,FUN=sum)

res2$BPP_or = res2$BPP_omission/res2$len_tau_true
res2$BPP_cr = res2$BPP_commission/res2$BPP_len_tau

res2$BPP_impliedpk_or = res2$BPP_impliedpk_omission/res2$len_tau_true
res2$BPP_impliedpk_cr = res2$BPP_impliedpk_commission/res2$BPP_impliedpk_len_tau

res2$BinSeg_or = res2$BinSeg_omission/res2$len_tau_true
res2$BinSeg_cr = res2$BinSeg_commission/res2$BinSeg_len_tau

res2$PELT_or = res2$PELT_omission/res2$len_tau_true
res2$PELT_cr = res2$PELT_commission/res2$PELT_len_tau

res2[] <- lapply(res2, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})


res = read.csv('main_study_results_norm_disc_v2.csv')
res[] <- lapply(res, function(x) as.numeric(as.character(x)))
res1 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res,FUN=sum)

res1$Norm_or = res1$Norm_omission/res1$len_tau_true
res1$Norm_cr = res1$Norm_commission/res1$Norm_len_tau

res1$Disc_or = res1$Disc_omission/res1$len_tau_true
res1$Disc_cr = res1$Disc_commission/res1$Disc_len_tau

res1[] <- lapply(res1, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}

#Subplots
ev = c(0.1,0.2,0.3); rb = c(3,10,100)
res1_list = list();res2_list = list()
#
for(j in 1:3){ res1_list[[j]] = res1[res1$time_dist==j, ];      res2_list[[j]] = res2[res2$time_dist==j,]       }
for(j in 1:3){ res1_list[[j+3]] = res1[res1$error_var==ev[j],]; res2_list[[j+3]] = res2[res2$error_var==ev[j],] }
for(j in 1:3){ res1_list[[j+6]] = res1[res1$robustness==rb[j],];res2_list[[j+6]] = res2[res2$robustness==rb[j],]}
pal = brewer.pal(8, 'RdBu')#'PuOr'

for(ll in 1:3){
	catn(ll)
	pdf(paste0('main_study_DiscNormBPPE_subset',ll,'_v1.pdf'), width = 12, height = 12)
	layout(matrix(1:9, nrow = 3, ncol = 3,byrow=T))
	#bottom,left,top,right
	par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.75, 0.65, .5, 0.35))
	for(lll in 1:3){
		catn(lll)
		res11 = res1_list[[lll + 3*(ll-1)]];
		res12 = res2_list[[lll + 3*(ll-1)]]
		#	
		x = res11$Norm_cr; y = 1-res11$Norm_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
				pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('BPP: nonrobust',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

		x = res11$Disc_cr; y = 1-res11$Disc_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
				pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('Noninf discrete time model',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)

		x = res12$BPP_impliedpk_cr; y = 1-res12$BPP_impliedpk_or;
		plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
		axis(1,c(0,1))
		axis(2,c(0,1))
		#axis(1, c(0.5) ,labels = FALSE)
		axis(2, c(0.5) ,labels = FALSE)
		title('BPPE',line=.5,adj = 0.09)#0.275)
		abline(0,1)
		z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
		contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")
	}
	dev.off()
}





####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
# Additional subset plots for Discrete, Normal, BPPE appendix
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

setwd('~/Documents/BPP/synthetic_study/')
#Discrete time, normal, and BPP_impliedpk
res2 = read.csv('main_study_results_v1.csv')
res2[] <- lapply(res2, function(x) as.numeric(as.character(x)))
res2 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res2,FUN=sum)

res2$BPP_or = res2$BPP_omission/res2$len_tau_true
res2$BPP_cr = res2$BPP_commission/res2$BPP_len_tau

res2$BPP_impliedpk_or = res2$BPP_impliedpk_omission/res2$len_tau_true
res2$BPP_impliedpk_cr = res2$BPP_impliedpk_commission/res2$BPP_impliedpk_len_tau

res2$BinSeg_or = res2$BinSeg_omission/res2$len_tau_true
res2$BinSeg_cr = res2$BinSeg_commission/res2$BinSeg_len_tau

res2$PELT_or = res2$PELT_omission/res2$len_tau_true
res2$PELT_cr = res2$PELT_commission/res2$PELT_len_tau

res2[] <- lapply(res2, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

res = read.csv('main_study_results_norm_disc_v2.csv')
res[] <- lapply(res, function(x) as.numeric(as.character(x)))
res1 = aggregate(.~time_dist+error_var+robustness+size_change+num_regimes,data=res,FUN=sum)

res1$Norm_or = res1$Norm_omission/res1$len_tau_true
res1$Norm_cr = res1$Norm_commission/res1$Norm_len_tau

res1$Disc_or = res1$Disc_omission/res1$len_tau_true
res1$Disc_cr = res1$Disc_commission/res1$Disc_len_tau

res1[] <- lapply(res1, function(x) {
  if (is.numeric(x)) x[is.nan(x)] <- 0
  x
})

jit<-function(x,sd = 0.01){
	x+rnorm(length(x),0,sd)
}

#Subplots
ns = 1:4; sc = c(0.1,0.3,0.5,0.7,0.9,1.1)
res1_list = list();res2_list = list()
#
for(j in 1:length(sc)){ res1_list[[j]] = res1[res1$size_change==sc[j], ];      res2_list[[j]] = res2[res2$size_change==sc[j],]}
for(j in 1:length(ns)){ res1_list[[j+length(sc)]] = res1[res1$num_regimes==ns[j],];res2_list[[j+length(sc)]] = res2[res2$num_regimes==ns[j],]}
pal = brewer.pal(8, 'RdBu')#'PuOr'

for(j in 1:length(res1_list)){
	pdf(paste0('main_study_DiscNormBPPE_ADDsubset',j,'_v1.pdf'), width = 12, height = 4)
	layout(matrix(1:3, nrow = 1, ncol = 3,byrow=T))
	#bottom,left,top,right
	par(cex.lab = 2., cex.axis = 2., cex.main = 2.25,mgp=c(1.5,1.,0),mai = c(.75, 0.65, .5, 0.35))
	#
	catn(ns[lll])
	res11 = res1_list[[j]];
	res12 = res2_list[[j]]
	#	
	x = res11$Norm_cr; y = 1-res11$Norm_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('BPP: nonrobust',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)
	catn('completed normal model')

	x = res11$Disc_cr; y = 1-res11$Disc_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
			pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('Noninf discrete time model',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal)
	catn('completed discrete model')

	x = res12$BPP_impliedpk_cr; y = 1-res12$BPP_impliedpk_or;
	plot(jit(x),jit(y),xlim=c(0,1),ylim=c(0,1),xlab='FPR',ylab='TPR',
		pch=3,col=adjustcolor('black',alpha.f=0.9),xaxt='n',yaxt='n',font.lab=2,cex=2)
	axis(1,c(0,1))
	axis(2,c(0,1))
	#axis(1, c(0.5) ,labels = FALSE)
	axis(2, c(0.5) ,labels = FALSE)
	title('BPPE',line=.5,adj = 0.09)#0.275)
	abline(0,1)
	z <- kde2d(x,y , n = 200,h=.5,lims=c(0,1,0,1))
	contour(z, lwd = 5, drawlabels = FALSE, add=T, col=pal) #, col=hcl.colors(10, "YlOrRd")
	catn('completed BPPE model')
	dev.off()
}


