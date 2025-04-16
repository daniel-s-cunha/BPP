library(parallel)
library(RColorBrewer)
library(changepoint)
setwd('~/Documents/BPP/')
source('./src/functions.R')
set.seed(17)

arg_list = list()

time_dist = 1:3 #c('Disc_Hypergeo','Cont_BPP_Beta(0.5,0.5)','Cont_BPP_Beta(2,2)')
error_var = c(0.1, 0.2, 0.3)
robustness = c(3,10,100)
size_change = c(0.1,0.3,0.5,0.7,0.9,1.1)
num_regimes = 1:4
repetition = 1:100

i = 1
for(a in time_dist){
	for(b in error_var){
		for(d in robustness){
			for(e in size_change){
				for(f in num_regimes){
					for(g in repetition){
						arg_list[[i]] = c(a,b,d,e,f,g)
						i = i+1
					}
				}
			}
		}
	}
}

#1 time dist
#2 error var
#3 robustness
#4 size_change
#5 num_regimes
#6 repetition
simulate_n_run <- function(l){
	source('./src/functions.R')
	n = 500
	k = l[5]
	chpts = c()
	if(l[1]==1){
		ti = (0:(n-1))/(n-1); 
		if(k>1){
			pos = 1
			for(j in 1:(k-1)){
				if(pos==(n-k+j-2)){
					cpt = pos
				}else{
					cpt = sample(pos:n-k+j-2,1) #leave at least 2 obs to estimate mean and variance
				}
				chpts = c(chpts,cpt)
				pos = cpt+2
			}			
		}
	}else{
		if(l[1]==2) ti = c(0,sort(rbeta(n-2,0.5,0.5)),1) else ti = c(0,sort(rbeta(n-2,2,2)),1)
		#ti = ti/max(ti)
		if(k>1){
			chpts = c(1)
			while(length(chpts)!=(k-1) | (1%in%chpts)) chpts = sample_BPP(ti,n,k)
		}	
	}
	#generate y
	if(length(chpts)==0){
		y = sqrt(l[2])*rt(n,df=l[3])
	}else{
		y = c(); prev=0; cnt = 1
		for(j in c(chpts,n)){
			y = c(y, sqrt(l[2])*rt(j-prev,df=l[3]) + cnt*l[4])
			prev = j
			cnt = cnt + 1
		}
	}
	tt = chpts/n
	w = 0.0225 #window 6 months total, 3 months before or after;
	#> 6/(500*16/30) #6 months times total_months^(-1)
	#total_months = total_days divided by 30 days/month
	#[1] 0.0225
	#
	#fit models
	h=2; psi=.1; lam=1; K=6;
	re_BPP = fit_EM2(ti,y,K,h,psi,lam,quants=c(0.5),intercept=T,geometric=F)
	taus_BPP = re_BPP$res_F[[1]]$taus/n
	oc_BPP = calc_omit_commit_window(taus_BPP,tt,window=w)
	#
	taus_nonInf = re_BPP$res_T[[1]]$taus/n
	oc_nonInf = calc_omit_commit_window(taus_nonInf,tt,window=w)
	#
	# re_GEO = fit_EM(ti,y,K,h,psi,lam,quants=c(0.5),intercept=T,geometric=T)
	# taus_GEO = re_GEO$res[[1]]$taus
	#
	res_binseg <- cpts(cpt.mean(y, method = "BinSeg",Q=(K-1) ))
    if(identical(res_binseg,numeric(0))) taus_binseg = NULL else taus_binseg = res_binseg/n
    oc_binseg = calc_omit_commit_window(taus_binseg,tt,window=w)
    #
    res_pelt <- cpts(cpt.mean(y, method = "PELT",Q=(K-1)))
    if(identical(res_pelt,numeric(0))) taus_pelt = NULL else taus_pelt = res_pelt/n
	oc_pelt = calc_omit_commit_window(taus_pelt,tt,window=w)
    #
    c(oc_BPP,oc_nonInf,oc_binseg,oc_pelt)
    return(c(oc_BPP,oc_nonInf,oc_binseg,oc_pelt,
    			length(tt),length(taus_BPP),length(taus_nonInf),length(taus_binseg),length(taus_pelt)))
}

safe_fnc <- function(l) {
  tryCatch({
    simulate_n_run(l)
  }, error = function(e) {
    paste("Error:", conditionMessage(e))
  })
}

ncores <- as.numeric(Sys.getenv("NSLOTS"))
cl <- makeCluster(ncores)
clusterExport(cl, c('df','cpts','cpt.mean','simulate_n_run'))
paste0('START TIME: ',format(Sys.time(), "%H:%M:%S"))
res_list <- clusterApply(cl, x = arg_list, fun = safe_fnc)
stopCluster(cl)
paste0('END TIME: ',format(Sys.time(), "%H:%M:%S"))
#
#
#
res = do.call(rbind, res_list)
colnames(res) = c('BPP_omission','BPP_commission',
					'BPP_impliedpk_omission','BPP_impliedpk_commission',
					'BinSeg_omission','BinSeg_commission',
					'PELT_omission','PELT_commission',
					'len_tau_true','BPP_len_tau','BPP_impliedpk_len_tau','BinSeg_len_tau','PELT_len_tau')

arg_df = do.call(rbind, arg_list)
colnames(arg_df) = c('time_dist','error_var','robustness','size_change','num_regimes','repetition')

res = cbind(arg_df,res)

write.csv(res,'main_study_results_v1.csv',row.names=F)




####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#This section runs the main study for the normal model and the discrete model
source('./src/functions.R')
set.seed(17)

arg_list = list()

time_dist = 1:3 #c('Disc_Hypergeo','Cont_BPP_Beta(0.5,0.5)','Cont_BPP_Beta(2,2)')
error_var = c(0.1, 0.2, 0.3)
robustness = c(3,10,100)
size_change = c(0.1,0.3,0.5,0.7,0.9,1.1)
num_regimes = 1:4
repetition = 1:100

i = 1
for(a in time_dist){
	for(b in error_var){
		for(d in robustness){
			for(e in size_change){
				for(f in num_regimes){
					for(g in repetition){
						arg_list[[i]] = c(a,b,d,e,f,g)
						i = i+1
					}
				}
			}
		}
	}
}


simulate_n_run2 <- function(l){
	source('./src/functions.R')
	n = 500
	k = l[5]
	chpts = c()
	if(l[1]==1){
		ti = (0:(n-1))/(n-1); 
		if(k>1){
			pos = 1
			for(j in 1:(k-1)){
				if(pos==(n-k+j-2)){
					cpt = pos
				}else{
					cpt = sample(pos:n-k+j-2,1) #leave at least 2 obs to estimate mean and variance
				}
				chpts = c(chpts,cpt)
				pos = cpt+2
			}			
		}
	}else{
		if(l[1]==2) ti = c(0,sort(rbeta(n-2,0.5,0.5)),1) else ti = c(0,sort(rbeta(n-2,2,2)),1)
		#ti = ti/max(ti)
		if(k>1){
			chpts = c(1)
			while(length(chpts)!=(k-1) | (1%in%chpts)) chpts = sample_BPP(ti,n,k)
		}	
	}
	#generate y
	if(length(chpts)==0){
		y = sqrt(l[2])*rt(n,df=l[3])
	}else{
		y = c(); prev=0; cnt = 1
		for(j in c(chpts,n)){
			y = c(y, sqrt(l[2])*rt(j-prev,df=l[3]) + cnt*l[4])
			prev = j
			cnt = cnt + 1
		}
	}
	tt = chpts/n
	w = 0.0225 #window 6 months total, 3 months before or after;
	#> 6/(500*16/30) #6 months times total_months^(-1)
	#total_months = total_days divided by 30 days/month
	#[1] 0.0225
	#
	#fit models
	h=2; psi=.1; lam=1; K=6;
	re_norm = fit_EM(ti,y,K,h,psi,lam,quants=c(0.5),intercept=T,normal=T,geometric=F)
	taus_norm = re_norm$res[[1]]$taus/n
	oc_norm = calc_omit_commit_window(taus_norm,tt,window=w)
	#
	re_disc = fit_EM(ti,y,K,h,psi,lam,quants=c(0.5),intercept=T,discrete=T,geometric=F)
	taus_disc = re_disc$res[[1]]$taus/n
	oc_disc = calc_omit_commit_window(taus_disc,tt,window=w)
    #
    c(oc_norm,oc_disc)
    return(c(oc_norm,oc_disc, length(tt),length(taus_norm),length(taus_disc)))
}

safe_fnc2 <- function(l) {
  tryCatch({
    simulate_n_run2(l)
  }, error = function(e) {
    paste("Error:", conditionMessage(e))
  })
}

ncores <- as.numeric(Sys.getenv("NSLOTS"))
cl <- makeCluster(ncores)
clusterExport(cl, c('df','cpts','cpt.mean','simulate_n_run2'))
source('./src/functions.R')
paste0('START TIME: ',format(Sys.time(), "%H:%M:%S"))
res_list <- clusterApply(cl, x = arg_list, fun = safe_fnc2)
stopCluster(cl)
paste0('END TIME: ',format(Sys.time(), "%H:%M:%S"))
#
#
#
res = do.call(rbind, res_list)
colnames(res) = c('Norm_omission','Norm_commission',
					'Disc_omission','Disc_commission',
					'len_tau_true','Norm_len_tau','Disc_len_tau')

arg_df = do.call(rbind, arg_list)
colnames(arg_df) = c('time_dist','error_var','robustness','size_change','num_regimes','repetition')

res = cbind(arg_df,res)

write.csv(res,'main_study_results_norm_disc_v2.csv',row.names=F)





####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################





#This section is for the gibbs sampler
source('./src/functions.R')
set.seed(17)

arg_list = list()

time_dist = 1:3 #c('Disc_Hypergeo','Cont_BPP_Beta(0.5,0.5)','Cont_BPP_Beta(2,2)')
error_var = c(0.1, 0.2, 0.3)
robustness = c(3,10,100)
size_change = c(0.1,0.3,0.5,0.7,0.9,1.1)
num_regimes = 1:4
repetition = 1:10

i = 1
for(a in time_dist){
	for(b in error_var){
		for(d in robustness){
			for(e in size_change){
				for(f in num_regimes){
					for(g in repetition){
						arg_list[[i]] = c(a,b,d,e,f,g)
						i = i+1
					}
				}
			}
		}
	}
}

#1 time dist
#2 error var
#3 robustness
#4 size_change
#5 num_regimes
#6 repetition
simulate_n_run <- function(l){
	source('./src/functions.R')
	n = 500
	k = l[5]
	chpts = c()
	if(l[1]==1){
		ti = (0:(n-1))/(n-1); 
		if(k>1){
			pos = 1
			for(j in 1:(k-1)){
				if(pos==(n-k+j-2)){
					cpt = pos
				}else{
					cpt = sample(pos:n-k+j-2,1) #leave at least 2 obs to estimate mean and variance
				}
				chpts = c(chpts,cpt)
				pos = cpt+2
			}			
		}
	}else{
		if(l[1]==2) ti = c(0,sort(rbeta(n-2,0.5,0.5)),1) else ti = c(0,sort(rbeta(n-2,2,2)),1)
		#ti = ti/max(ti)
		if(k>1){
			chpts = c(1)
			while(length(chpts)!=(k-1) | (1%in%chpts)) chpts = sample_BPP(ti,n,k)
		}	
	}
	#generate y
	if(length(chpts)==0){
		y = sqrt(l[2])*rt(n,df=l[3])
	}else{
		y = c(); prev=0; cnt = 1
		for(j in c(chpts,n)){
			y = c(y, sqrt(l[2])*rt(j-prev,df=l[3]) + cnt*l[4])
			prev = j
			cnt = cnt + 1
		}
	}
	tt = chpts/n
	w = 0.0225 #window 6 months total, 3 months before or after;
	#> 6/(500*16/30) #6 months times total_months^(-1)
	#total_months = total_days divided by 30 days/month
	#[1] 0.0225
	#
	#fit models
	h=2; psi=.1; lam=1; K=6;
	re_BPP = fit_gibbs(ti,y,K,h,psi,lam,discrete=F,geometric=F,normal=F,intercept=T,nsamp=500)
	taus_BPP = re_BPP$res$taus/n
	oc_BPP = calc_omit_commit_window(taus_BPP,tt,window=w)
  #
  c(oc_BPP)
  return(c(oc_BPP, length(tt),length(taus_BPP)))
}

safe_fnc <- function(l) {
  tryCatch({
    simulate_n_run(l)
  }, error = function(e) {
    paste("Error:", conditionMessage(e))
  })
}

ncores <- as.numeric(Sys.getenv("NSLOTS"))
cl <- makeCluster(ncores)
clusterExport(cl, c('df','cpts','cpt.mean','simulate_n_run'))
paste0('START TIME: ',format(Sys.time(), "%H:%M:%S"))
res_list <- clusterApply(cl, x = arg_list, fun = safe_fnc)
stopCluster(cl)
paste0('END TIME: ',format(Sys.time(), "%H:%M:%S"))
#
#
#
res = do.call(rbind, res_list)
colnames(res) = c('BPP_omission','BPP_commission',
					'len_tau_true','BPP_len_tau')

arg_df = do.call(rbind, arg_list)
colnames(arg_df) = c('time_dist','error_var','robustness','size_change','num_regimes','repetition')

res = cbind(arg_df,res)

write.csv(res,'main_study_results_gibbs_v2_nsamp500.csv',row.names=F)
