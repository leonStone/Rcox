# type.gaussian="naive" 20211213
#if((like_diff > 0) | (abs(like_diff) <= 1e-03) | (iter>100)) 20211213

# up.alpha, add ifelse, add alpha as a parameter. 20211216
# up.delta, add ifelse

# add updelta_optim & updelta_spg
# add parallel package
# release RAM



###################################################
#  inner  functions for PGKM-COX
datf <- function(datid,P,Q,PT,QT,N,Nvalid,Ntest){
  # N for training
  # Nvalid for validation
  # Ntest for prediction
  #X<-matrix(runif((N+Ntest+Nvalid)*P,-0.01,0.01),(N+Ntest+Nvalid),P);
  X<-matrix(rnorm((N+Ntest+Nvalid)*P,0,1),(N+Ntest+Nvalid),P);
  #X<-matrix(runif((N+Ntest+Nvalid)*P,-0.5,0.5),(N+Ntest+Nvalid),P);
  #Z<-matrix(runif((N+Ntest+Nvalid)*Q,0,0.5),(N+Ntest+Nvalid),Q);
  Z<-matrix(rnorm((N+Ntest+Nvalid)*Q,0,1),(N+Ntest+Nvalid),Q);
  true_beta<-c(rep(1,PT),rep(0,P-PT))
  #true_beta<-c(rep(0.01,PT),rep(0,P-PT))
  #lasso data
  XL<-cbind(X,Z)
  
  if(datid == 1){
    str <- 'H = 0.01 * (cos(Z[,1]) - 1.5 * (Z[,2])^2 + exp(-Z[,3]) * Z[,4] - 0.8 * sin(Z[,5]) * cos(Z[,3]) + 2 * Z[,1] * Z[,5] )'
  }else if(datid == 2){
    str <- 'H = cos(Z[,1]) - 1.5 * (Z[,2])^2 + exp(-Z[,3]) * Z[,4] - 0.8 * sin(Z[,5]) * cos(Z[,3]) + 2 * Z[,1] * Z[,5] + 0.9 * Z[,6] * sin(Z[,7]) - 0.8 * cos(Z[,6]) * Z[,7] + 2 * Z[,8] * sin(Z[,9]) * sin(Z[,10]) - 1.5 * Z[,8]^3 - Z[,8] * Z[,9] - 0.1 * exp(Z[,10]) * cos(Z[,10])'
  }else if(datid == 3){
    str <- 'H = 10*cos(Z[,1]) - 15*(Z[,2])^2 + 10*exp(-Z[,3])*Z[,4] - 8 * sin(Z[,5]) * cos(Z[,3]) + 20 * Z[,1] * Z[,5]'
  }else if(datid == 4){
    str <- 'H = - 0.01 * (cos(Z[,1]) - 1.5 * (Z[,2])^2 + exp(-Z[,3]) * Z[,4] - 0.8 * sin(Z[,5]) * cos(Z[,3]) + 2 * Z[,1] * Z[,5] ) '                
  }else if(datid == 599){
    str <- 'H = 0.06*(10*cos(Z[,1])* Z[,2] + 6* (Z[,1])^2 - 5*exp(Z[,1]) * Z[,2] - 6* sin(Z[,2]) * cos(Z[,3]) + 10* exp(Z[,3]) * sin(Z[,4])-8* Z[,2] * sin(Z[,4])- 2* cos(Z[,3])*Z[,4]^2 -2* exp(Z[,4]) * cos(Z[,5])-8*sin(Z[,4])* Z[,5]^2)' 
  }else if(datid == 699){
    str <- 'H = 5*cos(Z[,1])* Z[,2]^2 - 10* (Z[,1])^3 + 8*exp(Z[,1]) * cos(Z[,2]) - 8* sin(Z[,2]) * cos(Z[,3]) + 7* exp(Z[,3]) * sin(Z[,4])+6* Z[,2] ^2* sin(Z[,4])- 10* Z[,4]^3 -6* exp(Z[,4]) * cos(Z[,5])-6*sin(Z[,4])* Z[,5]^2' 
  }else if(datid == 899){
    str <- 'H = 0.06*(10*cos(Z[,1])* Z[,2] + 6* (Z[,1])^2 - 5*exp(Z[,1]) * Z[,2] - 6* sin(Z[,2]) * cos(Z[,3]) + 10* exp(Z[,3]) * sin(Z[,4])-8* Z[,2] * sin(Z[,4])- 2* cos(Z[,3])*Z[,4]^2 -2* exp(Z[,4]) * cos(Z[,5])-8*sin(Z[,4])* Z[,5]^2)' 
  }else if(datid==799){
    str <- 'H = 2*cos(Z[,1])* Z[,2] - 2* (Z[,1])^2 + 3.5*exp(Z[,1]) * Z[,2] - 3* sin(Z[,2]) * cos(Z[,3]) + 2* exp(Z[,3]) * sin(Z[,4])+5* Z[,2] * sin(Z[,4])- 1.5* cos(Z[,3])*Z[,4]^2 -3* exp(Z[,4]) * cos(Z[,5])-4*sin(Z[,4])* Z[,5]^2' 
  }else if(datid == 100){
    str <- 'H = 10*cos(Z[,1])* Z[,2] + 6* (Z[,1])^2 - 5*exp(Z[,1]) * Z[,2] - 6* sin(Z[,2]) * cos(Z[,3]) + 10* exp(Z[,3]) * sin(Z[,4])-8* Z[,2] * sin(Z[,4])- 2* cos(Z[,3])*Z[,4]^2 -2* exp(Z[,4]) * cos(Z[,5])-8*sin(Z[,4])* Z[,5]^2' 
  }else if(datid == 101){
    str <- 'H = 10 * Z[,1] + Z[,2]^2 + 20 * Z[,3] * Z[,4] + Z[,5]' 
  }
  
  eval(parse(text=str))
  #   print(paste(str,'P=',P,'Q=',Q,'PT=',PT,'QT=',QT))
  
  fx <- as.vector(X%*%true_beta) # X * beta
  y <- fx + H                    # X * beta + h(Z)
  g <- exp(y)                    # hazard function
#  D <- rexp(N+Ntest+Nvalid,g)    # death time
  D <- (unique(rexp((N+Ntest+Nvalid)*2,g)))[1:(N+Ntest+Nvalid)]
  U <- runif(N+Ntest+Nvalid,2000000000,3000000000)# censoring rate is about 0%
  #U <- runif(N+Ntest+Nvalid,5000,6000)# censoring rate is about 20%
  #U <- runif(N+Ntest+Nvalid,1,2)# censoring rate is about 50%
  C <- rexp(N+Ntest+Nvalid, 1/(U*g)) 
  # C <- runif(N+Ntest+Nvalid,0,100) # 0,100 needed to modified for censoring rate # censoring rate: 30%
  TIM <- ifelse(D <= C, D, C)      # survival time
  TAU <- ifelse(D <= C, 1, 0)    # death indicator 
  
  train <- 1:N                                                # taining data id
  valid <- (N+1):(N+Nvalid)                                   # validation data id
  test <- (N + Nvalid+1):(N + Nvalid + Ntest)                 # testing data id
  
  ZT <- as.matrix(Z[test,]);Zvalid <- as.matrix(Z[valid,]);Z <- as.matrix(Z[train,]) # train data must be listed at last
  XT <- as.matrix(X[test,]);Xvalid <- as.matrix(X[valid,]);X <- as.matrix(X[train,])
  XLT <- as.matrix(XL[test,]);XLvalid <- as.matrix(XL[valid,]);XL <- as.matrix(XL[train,])
  TIMvalid <- TIM[valid]; TIMT <- TIM[test]; TIM <- TIM[train]
  TAUvalid <- TAU[valid]; TAUT <- TAU[test]; TAU <- TAU[train]
  
  fxT <- as.matrix(fx[test]);fxvalid <- as.matrix(fx[valid]);fx <- as.matrix(fx[train]) 
  HT <- as.matrix(H[test]);Hvalid <- as.matrix(H[valid]);H <- as.matrix(H[train]) 
  
  
  orde_valid <- order(TIMvalid); orde_t <- order(TIMT); orde <- order(TIM) # order for Tvalid, TT, T
  TIMvalid <- sort(TIMvalid); TIMT <- sort(TIMT); TIM<- sort(TIM)                # sort Tvalid, TT, T in increasing order
  ZT <- as.matrix(ZT[orde_t,]); Zvalid <- as.matrix(Zvalid[orde_valid,]); Z <- as.matrix(Z[orde,])    # sort Z based on T increasing order
  XT <- as.matrix(XT[orde_t,]); Xvalid <- as.matrix(Xvalid[orde_valid,]); X <- as.matrix(X[orde,])
  XLT <- as.matrix(XLT[orde_t,]); XLvalid <- as.matrix(XLvalid[orde_valid,]); XL <- as.matrix(XL[orde,])
  TAUT <- TAUT[orde_t]; TAUvalid <- TAUvalid[orde_valid]; TAU <- TAU[orde]
  
  
  HT <- as.matrix(HT[orde_t]); Hvalid <- as.matrix(Hvalid[orde_valid]); H <- as.matrix(H[orde])
  fxT <- as.matrix(fxT[orde_t]); fxvalid <- as.matrix(fxvalid[orde_valid]); fx <- as.matrix(fx[orde])
  
  
  res <- list(X=X,Z=Z,XL=XL,TIM=TIM,TAU=TAU,
              XT=XT,ZT=ZT,XLT=XLT,TIMT=TIMT,TAUT=TAUT,
              Xvalid=Xvalid,Zvalid=Zvalid,XLvalid=XLvalid,TIMvalid=TIMvalid,TAUvalid=TAUvalid,
              HT=HT,Hvalid=Hvalid,H=H,
              fxT=fxT,fxvalid=fxvalid,fx=fx)
  res
}

# u Q*N*N
uf <- function(Z){  # uf, used in gradient of f function in updating delta, should be listed in data.R
  #   if(is.null(dim(Z))){Z <- matrix(Z,ncol=1)}
  #   N <- ifelse(is.null(nrow(Z)),length(Z),nrow(Z))
  #   Q<-ifelse(is.null(ncol(Z)),1,ncol(Z))#Z is N*Q
  N <- nrow(Z)
  Q <- ncol(Z)
  u <- array(0,c(Q,N,N))
  for(i in 1:N){
    u[,,i]<-(Z[i,]-t(Z))^2/Q
  }	
  u <- array(u,c(Q,N*N))
  gc()
   
	#cl <-makeCluster(4)
	#u <- parSapply(cl,1:N ,function(i,Z,Q){(Z[i,]-t(Z))^2/Q},Z=Z,Q=Q)
	#u <- array(u,c(Q,N*N))
	#stopCluster(cl)
	#gc()
	
	return(u)
}



# modi 20220808 short the computing time 10*  
# K N*N
Kf <- function(Z,delta){
  #N <- ifelse(is.null(nrow(Z)),length(Z),nrow(Z))
  N <- nrow(Z)
  Q <- ncol(Z)
  Zsqu <- Z^2 %*% delta %*% t(rep(1,N))
  #res <- Zsqu - 2 * Z %*% diag(delta,length(delta)) %*% t(Z) + t(Zsqu)  #  20220808 Q*Q too big to excute, cost so mucu time  Q = 20K, 100s
  res <- Zsqu -2 * (Z*matrix(delta,nrow=N,ncol=Q,byrow=TRUE) )%*% t(Z)  + t(Zsqu)  
  K <- exp(-res/Q) # modi Q = 1000
  return(K)  # N*N
}

# h N*1
hf <- function(Z,alpha,delta){
  h <- Kf(Z,delta) %*% alpha 
  return(h) # N *1
}

# output:eta=X %*% betaa + hf(Z,alpha,delta)
# eta N*1
etaf <- function(X,Z,alpha,betaa,delta){  # eta n*1 matrix!
  eta <- X %*% betaa + Kf(Z,delta) %*% alpha # modi 20220808
  return(eta)
}

#output:R=R, len=len, E=E, len_E=len_E
Rf <- function(TIM){   # Rf should be listed in data.R
  
  N <- length(TIM)
  R <- matrix(0,N,N)  # risk set
  len <- c()          # length of R[i,]
  
	for (i in 1:N){
		Ri <- which(TIM >= TIM[i]) 
		R[i,1:length(Ri)] <- Ri; 
	
		#R[i,-(1:len.r)] <- 0
		#len.r <- length(Ri) #
		#len <- c(len,len.r)
	} 
  
	#modi 20211216
	len <- apply(R,1,function(x) sum(x!=0))
  
  E <- matrix(0,N,N) # E_i: the risk sets containing individual i
  len_E <- c()
  
  for (i in 1:N){
    Ei <- c()
    for (m in 1:N){
      if( is.element(i,R[m,])){
        Ei=c(Ei,m)
      }else{Ei=Ei}	# needed to modified: ifelse
    }
  
	len.e <- length(Ei)
	E[i,1:len.e] <- Ei; 
	E[i,-(1:len.e)] <- 0
    len_E <- c(len_E,len.e)
  }
  #len_E <- apply(E,1,function(x),sum(x!=0))
  
  return(list( R=R, len=len, E=E, len_E=len_E ))
}

# output: sum_exp is a N*1 vector, where i-th elements is sum_{l in R_i} exp(X_l beta+K_l alpha) for j in R_i
# similar to Simon et al. (2011) Page 6.
sum.exp <- function(X,Z,R,len,alpha,betaa,delta){
	#N <- dim(X)[1]
	#sum_exp <- c()
	
	# modi 20220808 "+1e-5" to avoid zero
	exp.eta <- exp(etaf(X,Z,alpha,betaa,delta)) + 1e-5

	#for (i in 1:N) {
	#	l <- R[i, (1:len[i])]   # len[i] vector
	#	sum_expi <- sum(exp.eta[l])
	#	sum_exp <-c(sum_exp,sum_expi)
	#}
	
	# modi 20211216
	Ra <- ifelse(R>0,1,0)
	sum_exp <- as.vector(Ra %*% exp.eta)

	return(sum_exp)  # N*1 vector 
}

# output: B is a Q*N matrix, where qi-th elements is (sum_{l in R_i} {exp(X_l beta+K_l alpha) * [(K[l,]*(-U[q,l,]))%*%alpha]}) /(sum_{l in R_i} exp(X_l beta+K_l alpha)) 
# another: similar to Simon et al. (2011) Page 6.
up.B <- function(X,Z,u,R,len,alpha,betaa,delta){
	N <- dim(Z)[1]
	Q <- dim(Z)[2]
	B <- matrix(0,Q,N) 
  
	eta <- etaf(X,Z,alpha,betaa,delta)
	sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta)
	
	# 20211216   # to be modified
	# modi 20220808 "+1e-5" to avoid zero
	exp.eta <- exp(eta) + 1e-5 
	kfz <- Kf(Z,delta)
	for (q in 1:Q){
		for (i in 1:N){
			l <- R[i, (1:len[i])]   # len[i] vector
			#B[q,i] <- sum((exp.eta[l]) * ( (kfz[l,] * ((-1) * (u[q,l, ])) ) %*% alpha )) / sum_exp[i]
			B[q,i] <- sum((exp.eta[l]) * ( (kfz[l,] * (u[q,l, ]) ) %*% alpha ))
		}
	}
	B <- (-1) * t(t(B) / sum_exp)  #modi 20220808
  
#  for (q in 1:Q){
#    for (i in 1:N){
#      l <- R[i, (1:len[i])]   # len[i] vector
#      B[q,i] <- sum((exp(eta[l])) * ( ((Kf(Z,delta)[l,]) * ((-1) * (u[q,l, ])) ) %*% alpha )) / sum_exp[i]
#    }
#  }
  
  return(B)
  
}  

# type.gaussian="naive" 20211213
#output: updated beta for given lambda_1
up.beta<-function(Y,X,Z,alpha,delta,lambda_1,weig){
  #  require(glmnet);require(lattice);require(Matrix)
  yy <- Y - hf(Z,alpha,delta)
  if(ncol(X)==1){
    #betaa <- as.vector(lm(yy~X,weights=weig)$coef[2])
    #inter <- lm(yy~X)$coef[1]
	fit <- lm(yy~X,weights=weig)  # modi 20220808
	betaa <- as.vector(fit$coef[2])
	inter <- fit$coef[1]
  }else{
    glmfit <- glmnet(x=X,y=yy,family="gaussian",type.gaussian="naive",weights=weig,alpha=1,lambda=lambda_1,standardize=FALSE)
    betaa <- as.vector(glmfit$beta)
    inter <- glmfit$a0
  }
  #print(paste('beta zero is ',inter))
  return(betaa)
}

#output: updated alpha for given lambda_3

up.alpha <- function(Y,X,Z,alpha,betaa,delta,lambda_3,weig){
	N <- length(Y)
	K = Kf(Z,delta)
	W <- diag(weig,nrow=N,ncol=N) #N*N weight matrix
	A = (1/N) * K %*% W %*% K + lambda_3 * K + diag(1e-6,N)
	B = (1/N)* K %*% W %*% (Y - X %*% betaa)
  
	#alpha_n <- try(solve(A,B)) ## try to solve a singular A
	
	alpha_n <- try(solve(A,B),silent=TRUE)
 	
	if("try-error" %in% class(alpha_n)){
		alpha_n <- alpha
	}

	return(alpha_n)  # N*1 matrix
}

# #output: updated delta for given lambda_2,lambda_3 # nlm, nlminb, constrOptim, optimize
up.delta.spg <- function(X,Z,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3){
	N <- dim(X)[1]
	Q <- ncol(Z)
    
	# update 20220808 add
	TAU <- TAU; R <- R; len <- len; alpha <- alpha; betaa <- betaa; delta <- delta; lambda_2 <- lambda_2;lambda_3<- lambda_3

    updelta_time1 <- Sys.time()  
    
	# modi 20220808 change lambda_2*sum(delta) to lambda_2*sum(delta_n)
	f <- function(delta_n,X,Z,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){
	  sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
	  A <- log(sum_exp)
	  
	  res <- (1/N) * (t(TAU) %*% (X %*% betaa + hf(Z,alpha,delta_n) - A)) - lambda_2 * sum(delta_n) - (1/2) *lambda_3 * (t(alpha) %*% Kf(Z,delta_n) %*% alpha)
	  res
    }
	
	fg <- function(delta_n,X,Z,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){      		## gradient of function 'f'
      eta <- etaf(X,Z,alpha,betaa,delta_n)
      sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
  	  I = rep(1,Q)
      V <- matrix(1,N,N)
      V[lower.tri(V)]=0
      G <- t(TAU) * (1/t(sum_exp)) %*% V * exp(t(eta)) * t(alpha)
      IKF <-  I %*% t(as.vector(Kf(Z,delta_n))) * uf(Z)#array(u,c(Q,N*N))
	  
	  res <-  as.vector( -(1/N)* IKF %*% kronecker(TAU,alpha) + (1/N)*IKF %*% (rep(G,each=N) * rep(alpha,N)) - lambda_2 * I  + (1/2) * lambda_3 *IKF %*% kronecker(alpha,alpha) )	  
      res
    }
 
	fg_res <- fg(delta_n = delta,X,Z,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q)
	if ((max(fg_res) < 0) & (max(delta) ==0)){
		delta_n <- delta
	} else if((min(fg_res) > 0) & (min(delta)==1)) {
		delta_n <- delta
	} else {
		spg_res <- spg(par=delta,fn=f,gr=fg,X=X,Z=Z,TAU=TAU,R=R,len=len,alpha=alpha,betaa=betaa,lambda_2=lambda_2,lambda_3=lambda_3,N=N,Q=Q,lower=rep(0,Q),upper=rep(1,Q),control=list(maximize=TRUE,checkGrad=FALSE,eps=1e-4)) #modi 20220808 1e-4
		delta_n <- spg_res$par
		delta_n <- ifelse(delta_n < 1e-4,0,delta_n)  #1e-3 # modi 20220808 1e-4
		print(c(spg_res$message,spg_res$iter,spg_res$value,sum(delta_n != 0)))
	}
	
	gc()
	
	updelta_time2 <- Sys.time()
	updelta_time <- round(difftime(updelta_time2,updelta_time1,units='secs'))
	print(paste(min(fg_res),max(fg_res),'## time for updating delta=',updelta_time,'seconds##'))
 
	return(delta_n)
}



#input: X,Z,TAU,R,len,E,len_E,alpha,betaa,delta
#output: Y, weigths, A
#first and second derivatives of log-partial likelihood, working response, weight matrix
#first derivative: grad, n*1 vector; second derivative: hess, n*1 vector
upda.esti <- function(X,Z,TAU,R,len,E,len_E,alpha,betaa,delta){
	N <- dim(X)[1]
	#grad <- rep(1,N) # gradient of log-partial likelihood
	#hess <- rep(1,N) # diagonal elements of hessien matrix of log-partial likelihood
	#Y <- rep(0,N)
	#weig <- rep(1,N)
	
	#sum_E <- rep(0,N) # whose i-th element is sum_{m in E_i} TAU_m/(sum_{l in Rm} exp(eta_l))
	#sum_Esqua <- rep(0,N) # whose i-th element is sum_{m in E_i} TAU_m/(sum_{l in Rm} exp(eta_l))^2
   
	eta <- as.vector(etaf(X,Z,alpha,betaa,delta))
	expeta <- exp(eta) + 1e-5 # modi 20220809 
   
	#sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta)
	#A <- log(sum_exp)

#	for (i in 1:N){   # caculate sum_E
#	for (m in (intersect(E[i,(1:len_E[i])], which(TAU == 1)))){     
#	  sum_E[i] <-  sum_E[i] + (1/sum_exp[m])
#	  sum_Esqua[i] <-  sum_Esqua[i] + (1/(sum_exp[m])^2)	  
#	  }
#	} 
	
	## 20211216
	Ea <- ifelse(E!=0,1,0)
	Ra <- ifelse(R!=0,1,0)
	sum_E <- (Ea %*% TAU) / (Ra %*% expeta)
	sum_Esqua <- (Ea %*% TAU) / (Ra %*% expeta)^2
	
	grad <- TAU - expeta * sum_E
	hess <- expeta * (expeta * sum_Esqua - sum_E)   
	# 20211216
	hess <- ifelse(abs(hess) < 1e-5,(-1)*1e-5, hess)
#   for (i in which(abs(hess) < 1e-05)){
#      hess[i] <- (-1)*1e-05
#   }
   
	Y <- eta - as.vector(grad / hess)
	weig <- as.vector((-1) * hess)   # N*1 vector: whose elements could not be zero! 
	
	return(list(Y=Y, weig=weig))#, A=A))   
}

 #"No changes in the last 10 iterations, break iterations",see EPSGO step 6
Qf <- function(X,Z,TAU,R,len,E,len_E,alpha_nn,beta_nn,delta_nn,lambda_1,lambda_2,lambda_3){
	N <- length(TAU)
	P <- ncol(X)
	Q <- ncol(Z)
	#if(!is.matrix(Z)){Z <- matrix(Z,nrow=N)}
	#if(!is.matrix(X)){X <- matrix(X,nrow=N)}
	alpha <- alpha_nn  # initial value 
	betaa <- beta_nn
	delta <- delta_nn
	  
	### first update
	# update Y,weights
	upda_esti_r <- upda.esti(X,Z,TAU,R,len,E,len_E,alpha,betaa,delta)
	Y <- upda_esti_r$Y; 
	weig <- upda_esti_r$weig

	#upest_timeB1 <- Sys.time()
	beta_n <- up.beta(Y,X,Z,alpha,delta,lambda_1,weig)
	alpha_n <- up.alpha(Y,X,Z,alpha,betaa,delta,lambda_3,weig)
	delta_n <- up.delta(X,Z,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3)
	 
	A_n  <- log(sum.exp(X,Z,R,len,alpha_n,beta_n,delta_n))# old
	like_new <- (1/N) * (t(TAU)%*% (X %*% beta_n + hf(Z,alpha_n,delta_n) - A_n)) -  lambda_1  * sum(abs(beta_n)) - lambda_2 * sum(delta_n) - 1/2 *lambda_3 * t(alpha_n) %*% Kf(Z,delta_n) %*% alpha_n

	upest_time1 <- Sys.time()
	### iteration
	iter = 0 
	repeat{
	#	gc()
	#	print(paste( "Qf iter = ",iter, "like_new = ", like_new, "num of non-zero delta_n = ",sum(delta_n != 0)))
		
		iter = iter + 1
		alpha <- alpha_n
		betaa <- beta_n
		delta <- delta_n
		A <- A_n
		like_old <- like_new

		# update Y,weights
		upda_esti_r <- upda.esti(X,Z,TAU,R,len,E,len_E,alpha,betaa,delta)
		Y <- upda_esti_r$Y; 
		weig <- upda_esti_r$weig

		beta_n <- up.beta(Y,X,Z,alpha,delta,lambda_1,weig)
		alpha_n <- up.alpha(Y,X,Z,alpha,betaa,delta,lambda_3,weig)
		delta_n <- up.delta(X,Z,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3)
		
		A_n <- log(sum.exp(X,Z,R,len,alpha_n,beta_n,delta_n))# new
		like_new <- (1/N) * (t(TAU)%*% (X %*% beta_n + hf(Z,alpha_n,delta_n) - A_n)) -  lambda_1  * sum(abs(beta_n)) - lambda_2 * sum(delta_n) - 1/2 *lambda_3 * t(alpha_n) %*% Kf(Z,delta_n) %*% alpha_n
		#(print(c(like_new,(1/N) * (t(TAU)%*% (X %*% beta_n + hf(Z,alpha_n,delta_n) - A_n)),1/2 *lambda_3 * t(alpha_n) %*% Kf(Z,delta_n) %*% alpha_n)))

		like_diff <- like_new - like_old   
		par_diff <- sum(abs(beta_n - betaa))/P + sum(abs(delta_n - delta))/Q + sum(abs(alpha_n-alpha))/N
		if((abs(like_diff) <= 1e-03) | (iter>100) | (par_diff < 1e-4)){  
			break
		}
	}
	upest_time2 <- Sys.time()
	upest_time <- round(difftime(upest_time2,upest_time1,units='secs'))
	#print(paste('Qf like_new',like_new,'## parameters difference =',par_diff,',## partial log-likelihood difference =',like_diff,'; iter=',iter,'; lambda=c(',lambda_1,lambda_2,lambda_3,') has done ###'))
	print(paste('lambda=c(',lambda_1,lambda_2,lambda_3,') has done ###','iter=',iter,'##cost time=',upest_time,'seconds##','like_new = ',like_new,'num of non-zero delta',sum(delta_n != 0)))

	#result<-list(iter=iter,betaa=beta_n,alpha=alpha_n,delta=delta_n)
	result<-list(iter=iter,betaa=betaa,alpha=alpha,delta=delta)  #result of iter-1 iteration.
	return(result)
}
##############
###caculate the estimate of h(Z) and A on testing ZT or validation Zvalid data for gaussian kernel
##############
# R,len for validation or testing time(Tvalid,TT)
hathaf <- function(Xvalid,Z,Zvalid,Rvalid,lenvalid,alpha,betaa,delta){ 
  Nvalid <- nrow(Zvalid)
  N <- nrow(Z)
  zf <- t(Z^2 %*% delta %*% t(rep(1,Nvalid)))
  ztf <- Zvalid^2 %*% delta %*% t(rep(1,N))
  zz <- 2*Zvalid %*% diag(delta,length(delta)) %*% t(Z)
  hat_h <- exp(-(zf + ztf - zz)) %*% alpha #Nvalid*1 matrix= Nvalid*N matrix %*% N*1 vector
  
  sum_exp <- c()
    for (i in 1:Nvalid) {
    l <- Rvalid[i, (1:lenvalid[i])]   # lenT[i] vector
    sum_expi <- sum(exp( as.matrix(Xvalid[l,]) %*% betaa + hat_h[l]))
    sum_exp <-c(sum_exp,sum_expi)
  }
  hat_A <- as.matrix(log(sum_exp)) # Nvalid*1 matrix
  return(list( hat_h=hat_h, hat_A=hat_A ) ) 
}

# validated partial log-likelihood vpll=(1/Nvalid)*(t(TAU) %*% (Xvalid %*% betaa + hathf(Z,Zvalid,alpha,delta) - A_valid)) (in real data analysis, use cross-validated partial log-likelihood)
# Q.func: generalized error function to be minimized, -1* validated partial log-likelihood
#                  traning: X,Z,u,TAU,R,len,E,len_E,
#                  initial value: alpha_nn,beta_nn,delta_nn,
#                  validation: Xvalid,Zvalid,TAUvalid,Rvalid,lenvalid
# arguments for Q.func are lambda_1,lambda_2,lambda_3
# Q.func must satisfy:(1) containing arguments(bounds=bounds,parms.coding=parms.coding,seed=seed,maxevals=maxevals,verbose=verbose)
#                     (2) must use q.val  the generalized error and must be returned as "list" format
# Q.func <- function(lambda_1,lambda_2,lambda_3,
#                    X,Z,u,TAU,R,len,E,len_E,alpha_nn,beta_nn,delta_nn,
#                    Xvalid,Zvalid,TAUvalid,Rvalid,lenvalid,
#                    bounds=bounds,parms.coding=parms.coding,seed=seed,maxevals=maxevals,verbose=verbose)


lamf<- function(X,Z,TAU,TIM,Xvalid,Zvalid,TAUvalid,TIMvalid,#Hvalid, modi 20220809 Hvalid delete in real data
                 lambda_1_seq,lambda_2_seq,lambda_3_seq){
	N<-nrow(Z)# number of training observations
	Nvalid <- nrow(Zvalid)# number of testing observations
	if(!is.matrix(Z)){Z <- matrix(Z,nrow=N)}
	if(!is.matrix(X)){X <- matrix(X,nrow=N)}
	P <- ncol(X)
	Q <- ncol(Z)
	#u <- uf(Z)
	LH <- c()
	ITER<-c()
	#REG <- matrix(ncol=3)
	BETA <- matrix(ncol=P)
	DELTA <- matrix(ncol=Q)
	ALPHA <- matrix(ncol=N)
	LAMBDA <- matrix(ncol=3)
	  
	R<-Rf(TIM)$R
	len<-Rf(TIM)$len
	E<-Rf(TIM)$E
	len_E<-Rf(TIM)$len_E

	beta_nn<- runif(P,0,1)#rep(0.5,P)
	alpha_nn <- runif(N,0,1)#rep(0.5,N)
	delta_nn <- runif(Q,0,1) #rep(0.5,Q)

	for(lambda_1 in lambda_1_seq){
		for(lambda_2 in lambda_2_seq){
			for(lambda_3 in lambda_3_seq){
				a_b_d <- Qf(X,Z,TAU,R,len,E,len_E,alpha_nn,beta_nn,delta_nn,lambda_1,lambda_2,lambda_3)
				alpha <- as.vector(a_b_d$alpha) 
				betaa <- as.vector(a_b_d$betaa)
				delta <- as.vector(a_b_d$delta)
				iter  <- a_b_d$iter
				ITER<-c(ITER,iter)

				Rfvalid <- Rf(TIMvalid)
				Rvalid<-Rfvalid$R
				lenvalid<-Rfvalid$len

				hath <- hathaf(Xvalid,Z,Zvalid,Rvalid,lenvalid,alpha,betaa,delta)	
				hat_h <- hath$hat_h # hat_h is a estimate of h(Zvalid)
				hat_A <- hath$hat_A # hat_A is a estimate of A_valid
				#lsf <- lm(Hvalid~hat_h)

				# omitting intercept
				LH<-c(LH,c((-1/Nvalid) * (t(TAUvalid) %*% ( Xvalid%*%betaa+hat_h-hat_A ))))
				#REG <- rbind(REG,c(lsf$coef,summary(lsf)$adj.r.squared))
				BETA <- rbind(BETA,betaa)
				DELTA <- rbind(DELTA,delta)
				ALPHA <- rbind(ALPHA,alpha)
				LAMBDA <- rbind(LAMBDA,c(lambda_1,lambda_2,lambda_3))
				
				# warmstart
				alpha_nn <- alpha+(1e-5)
				beta_nn <- betaa+(1e-5)
				delta_nn <- delta+(1e-5)
			}
		}
	}
	
	id <- which.min(LH)
#	res <- list(LH = LH[id],REG = REG[id+1,],BETA = BETA[id+1,],DELTA = DELTA[id+1,],ALPHA=ALPHA[id+1,],LAMBDA = LAMBDA[id+1,],ITER = ITER[id]) 
	res <- list(LH = LH[id],BETA = BETA[id+1,],DELTA = DELTA[id+1,],ALPHA=ALPHA[id+1,],LAMBDA = LAMBDA[id+1,],ITER = ITER[id]) 
	res
}


cou <- function(x,PT){
  a <- all(x[1:PT]!=0)
  b <- all(x[-(1:PT)]==0)
  c <- any(x[1:PT]==0)
  c(a*b,a-a*b,c)
}

#lasso-COX
hata_lasso <- function(XLvalid,Rvalid,lenvalid,betaa){ 
  Nvalid <- nrow(XLvalid)
  sum_exp <- c()
  for (i in 1:Nvalid) {
    l <- Rvalid[i, (1:lenvalid[i])]   # lenT[i] vector
    sum_expi <- sum(exp(as.matrix(XLvalid [l,]) %*% betaa))
    sum_exp <-c(sum_exp,sum_expi)
  }
  hat_A <- as.matrix(log(sum_exp)) # Nvalid*1 matrix
  return(list(hat_A=hat_A)) 
}

lamf_lasso<- function(XL,TIM,TAU,XLvalid,TIMvalid,TAUvalid,lambda_lasso){
  Nvalid<-length(TIMvalid) # number of validation observations
  D <- ncol(XL)
  LH<-c()
  BETA <- matrix(ncol=D)
  LAMBDA <- matrix(ncol=1)
  S<-cbind(time = TIM, status = TAU)
  
  for(lambda_1 in lambda_lasso){
    a_b_c<-glmnet(x=XL,y=S,family = "cox",alpha=1,lambda = lambda_1,standardize=FALSE)
    betaa <- as.vector(a_b_c$beta)
    Rfvalid <- Rf(TIMvalid)
    Rvalid<-Rfvalid$R
    lenvalid<-Rfvalid$len
    hath <- hata_lasso(XLvalid,Rvalid,lenvalid,betaa)	
    hat_A <- hath$hat_A 
    LH<-c(LH,c((-1/Nvalid)* (t(TAUvalid) %*% (XLvalid%*%betaa-hat_A))))
    BETA <- rbind(BETA,betaa)
    LAMBDA <- rbind(LAMBDA,c(lambda_1))
  }
  id <- which.min(LH)
  res <- list(LH = LH[id],BETA = BETA[id+1,],LAMBDA = LAMBDA[id+1,]) 
  res
}

coul <- function(x,PT,P,QT,Q){
  a <- all(x[1:PT]!=0 & x[(P+1):(P+QT)]!=0)
  b <- all(x[(PT+1):P]==0 & x[(P+QT+1):(P+Q)]==0)
  c <- any(x[1:PT]==0 | x[(P+1):(P+QT)]==0)
  c(a*b,a-a*b,c)
}

# modi 20220808 add
ST<-function(dat){
  # for(i in 1:ncol(dat)){
  #   dat[,i]<-(dat[,i]-min(dat[,i]))/(max(dat[,i])-min(dat[,i]))
  # }
  # dat
  
  dat <- apply(dat,2,function(x) (x-min(x))/(max(x)-min(x)))
  dat
}

# optim
# optimize for one q step
up.delta.optim <- function(X,Z,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3){
	N <- dim(X)[1]
	#Q<-ifelse(is.null(ncol(Z)),1,ncol(Z))#Z is N*Q
	Q <- ncol(Z)
	delta_n <- delta
  
	updelta_time1 <- Sys.time()  
	for(q in which(delta_n !=0)){ # modi 20220808 kick out the zero delta to save time, greed
		gc()
		f <- function(delta_hat){
			delta_n[q] <- delta_hat
			sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
			A <- log(sum_exp)
		  
			res <- (1/N) * (t(TAU) %*% (X %*% betaa + hf(Z,alpha,delta_n) - A)) - lambda_2 * sum(delta_n)
			- (1/2) *lambda_3 * (t(alpha) %*% Kf(Z,delta_n) %*% alpha)
			-res
		}
		fg <- function(delta_hat){      		## gradient of function 'f'
			delta_n[q] <- delta_hat
			eta <- etaf(X,Z,alpha,betaa,delta_n)
			sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
		  
			V<-matrix(1,N,N)
			V[lower.tri(V)]=0
			G <- V%*%(exp(eta) * ( ((Kf(Z,delta_n)) * ((-1) * (u[q,,])) ) %*% alpha ))
		  
			res <- (1/N)*(t(TAU) %*% (Kf(Z,delta_n) * (-u[q,,])) %*% alpha) - (1/N)*(t(TAU) %*% (G/sum_exp))
			- lambda_2 - (1/2) * lambda_3 *(t(alpha) %*% (Kf(Z,delta_n) * (-u[q,,])) %*% alpha) 
			-res
		}
	
		#delta_n[q]<-optim(delta[q],f,fg,method="L-BFGS-B",lower=0,upper=1,control = list(factr=1e-5))$par	
		delta_n[q] <- optimize(f,c(0,1))$minimum
		if(q%%10 ==0){print(paste(q,Sys.time(),sep=','))}
	}
  
	updelta_time2 <- Sys.time()
	updelta_time <- round(difftime(updelta_time2,updelta_time1,units='secs'))
	print(paste('## time for updating delta=',updelta_time,'seconds##'))
  
	delta <- ifelse(delta_n<1e-4,0,delta_n)
	return(delta)
}

# optimParallel
up.delta.optimParallel <- function(X,Z,TAU,R,len,alpha,betaa,delta,lambda_2,lambda_3){
	N <- dim(X)[1]
	#Q<-ifelse(is.null(ncol(Z)),1,ncol(Z))#Z is N*Q
	Q <- ncol(Z)
	
	# maximation need checkGrad
	# without functions 
	f <- function(delta_n,X,Z,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){
		#sum_exp <- sum.exp(X,Z,R,len,alpha,betaa,delta_n)
		Zsqu <- Z^2 %*% delta_n %*% t(rep(1,N))
		res <- Zsqu -2 * (Z*matrix(delta_n,nrow=N,ncol=Q,byrow=TRUE) )%*% t(Z)  + t(Zsqu)  
		K <- exp(-res/Q) # modi Q = 1000
		h <- K %*% alpha
		eta <- X %*% betaa + h
		exp.eta <- exp(eta) + 1e-5
		sum_exp <- as.vector(ifelse(R>0,1,0) %*% exp.eta)
		A <- log(sum_exp)
	
		res <- (1/N) * (t(TAU) %*% (X %*% betaa + h - A)) - lambda_2 * sum(delta_n)
		- (1/2) *lambda_3 * (t(alpha) %*% K %*% alpha)
		-res
	}
		
	fg <- function(delta_n,X,Z,TAU,R,len,alpha,betaa,lambda_2,lambda_3,N,Q){      		## gradient of function 'f'
		Zsqu <- Z^2 %*% delta_n %*% t(rep(1,N))
		res <- Zsqu -2 * (Z*matrix(delta_n,nrow=N,ncol=Q,byrow=TRUE) )%*% t(Z)  + t(Zsqu)  
		K <- exp(-res/Q) # modi Q = 1000
		h <- K %*% alpha
		eta <- X %*% betaa + h
		exp.eta <- exp(eta) + 1e-5
		sum_exp <- as.vector(ifelse(R>0,1,0) %*% exp.eta)
		
		u <- array(0,c(Q,N,N))
		for(i in 1:N){
		u[,,i]<-(Z[i,]-t(Z))^2/Q
		}	
		u <- array(u,c(Q,N*N))
		gc()

		I = rep(1,Q)
		V <- matrix(1,N,N)
		V[lower.tri(V)]=0
		G <- t(TAU) * (1/t(sum_exp)) %*% V * exp(t(eta)) * t(alpha)
		IKF <-  I %*% t(as.vector(K)) * u#f(Z)#array(u,c(Q,N*N))
	  
		res <-  as.vector( -(1/N)* IKF %*% kronecker(TAU,alpha) + (1/N)*IKF %*% (rep(G,each=N) * rep(alpha,N)) - lambda_2 * I  + (1/2) * lambda_3 *IKF %*% kronecker(alpha,alpha) )	  
		-res
    }
	
	ca <- makeCluster(7)
	setDefaultCluster(cl=ca)
	opt.par <-optimParallel(par=delta,fn=f,gr=fg,X=X,Z=Z,TAU=TAU,R=R,len=len,alpha=alpha,betaa=betaa,lambda_2=lambda_2,lambda_3=lambda_3,N=N,Q=Q,lower=rep(0,Q),upper=rep(1,Q),control = list(factr=1e-5))	# maximization default
	delta_n <- opt.par$par
	setDefaultCluster(cl=NULL)
	stopCluster(ca)
	print(paste(opt.par$value,opt.par$message,opt.par$convergence,sep=','))

  
#	updelta_time2 <- Sys.time()
#	updelta_time <- round(difftime(updelta_time2,updelta_time1,units='secs'))
#	print(paste('## time for updating delta=',updelta_time,'seconds##'))
  
	delta <- ifelse(delta_n<1e-4,0,delta_n)
	
	gc(TRUE)
	return(delta)
}
