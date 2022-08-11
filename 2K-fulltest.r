

rm(list=ls())
gc(TRUE)
gc()


#setwd("D:/rong-paper/2022COXrealdata")

Hdir <- '../'

#require('mice')
require('glmnet')
require('BB') # spg
require('optimParallel')
require('survivalROC')
require('survC1')

lambda_1_seq <-seq(2^-8,2^-3,length.out = 5)
#lambda_1_seq <-seq(2^-7,2^-2,length.out = 10)
lambda_2_seq <-seq(2^-9,2^-3,length.out = 5)
#lambda_2_seq <-seq(2^-9,2^-4,length.out = 10)
lambda_3_seq <-seq(2^-6.5,2^-3,length.out = 5)
#lambda_3_seq <-seq(2^-6.5,2^-2,length.out = 10)
#lambda_lasso=0
#lambda_lasso <-seq(2^-8.5,2^-3.5,length.out = 10)
#lambda_1_seq <-seq(2^-5,2^-1.5,length.out = 10)
#lambda_2_seq <-seq(2^-7,2^-2.5,length.out = 10)
#lambda_3_seq <-seq(2^-6.5,2^-3,length.out = 10)
#lambda_cosso <-seq(2^0,2^3,length.out = 10)


source('funs-real.r')
load(paste(Hdir,'realdata.RData',sep=''))

up.delta <- up.delta.optimParallel
#up.delta <- up.delta.spg

# remove the imcomplete data
id0 <- unique(which(clin=='NA',arr.ind=TRUE)[,1])
dat <- dat[-id0,]
clin <- clin[-id0,]
rm(id0)

# unknow why
id0 <- which((clin$OS.TIME != 0) &(clin$Status=='Training'))
dat <- dat[id0,]
clin <- clin[id0,]
rm(id0)

## 2K  
subid <- 2000#ncol(dat)#1:2000 
dat <- dat[,1:subid]

## repetition
nre <- 1000
nam <- paste(subid,format(Sys.time(), "%y%m%d"),sep='-')
RES <- matrix(0,ncol=(3+5+20162))[-1,] 
#RES_lasso <- matrix(0,ncol=(1+5+20162))[-1,]
#RES_cosso <- matrix(0,ncol=(1))[-1,]

LH_PGKM <- AUC_PGKM <-  Cs_PGKM <- c()
FPS_PGKM<-c()
TPS_PGKM<-c()
#FPS_LASSO<-c()
#TPS_LASSO<-c()
#FPS_COSSO<-c()
#TPS_COSSO<-c()
pgckmTime <- matrix(0,ncol=1)[-1,] ## PGKM takes time

sink(file=paste(Hdir,nam,'.log',sep=''))
for(re in 1:nre){
	set.seed(re)
########### attention! original data need initial every repeat	  
    TIM<-clin$OS.TIME
	TAU<-clin$OS.CENSOR
	X<-apply(as.matrix(clin[,3:7]),2,as.numeric)  # 
	Z<-apply(dat,2,as.numeric)

#########train set
	train<- sample(1:nrow(Z),nrow(Z)*0.6)
	remain_X <- as.matrix(X[-train,])
	remain_Z <- as.matrix(Z[-train,])
	remain_TIM <- TIM[-train]
	remain_TAU <- TAU[-train]
  
	X<- as.matrix(X[train,])
	Z<- as.matrix(Z[train,])
	TIM <- TIM[train]
	TAU <- TAU[train]

#######vali and test set
	valid<- sample(1:nrow(remain_Z), nrow(remain_Z)*0.5)
	Xvalid <- as.matrix(remain_X[valid,])
	XT <- as.matrix(remain_X[-valid,])
	Zvalid <- as.matrix(remain_Z[valid,])
	ZT <- as.matrix(remain_Z[-valid,])

	TIMvalid <- remain_TIM[valid]
	TIMT <- remain_TIM[-valid]
	TAUvalid <- remain_TAU[valid]
	TAUT <- remain_TAU[-valid]

	Ntest<-length(TIMT)

#############
#####PGKM-cox method
	pgckmTime1 <- Sys.time()
	qt <- lamf(X,Z,TAU,TIM,Xvalid,Zvalid,TAUvalid,TIMvalid,lambda_1_seq,lambda_2_seq,lambda_3_seq)
	pgckmTime2 <- Sys.time()
	pgckmtime <- round(difftime(pgckmTime2,pgckmTime1,units='secs'))
	pgckmTime <- rbind(pgckmTime,as.vector(pgckmtime))
	save(qt,file=paste(Hdir,nam,'-qt.RData',sep=''))
	
	RT<-Rf(TIMT)$R;lenT<-Rf(TIMT)$len
	hath <- hathaf(XT,Z,ZT,RT,lenT,as.vector(qt$ALPHA),as.vector(qt$BETA),as.vector(qt$DELTA))
	#likelihood
	hat_h <- hath$hat_h 
	hat_A <- hath$hat_A 
	LH1 <-(1/Ntest)* (t(TAUT) %*% ( XT%*%as.vector(qt$BETA)+hat_h-hat_A ))
	LH_PGKM <- c(LH_PGKM,LH1)
	
	#AUC
	cutoff<-max(TIMT)*0.9
	mark<-XT%*%as.vector(qt$BETA)+hat_h
	ROC_PGKM<-survivalROC(Stime=TIMT,status = TAUT,marker = mark,predict.time = cutoff,span=0.25*100^(-0.2))
	AUC_PGKM <- c(AUC_PGKM,ROC_PGKM$AUC)
	
	FPS_PGKM<-rbind(FPS_PGKM,PGKM_ROC$FP)
	TPS_PGKM<-rbind(TPS_PGKM,PGKM_ROC$TP)
	
	#C-statistics
	ta<-quantile(TIMT,0.7)
	mydata<-cbind(TIMT,TAUT,mark)
	Cs_PGKM<-c(Cs_PGKM,Est.Cval(mydata, ta, nofit=TRUE)$Dhat)
	
	## PGKM-cox result	  
	RES <- rbind(RES, c(as.vector(qt$LAMBDA),as.vector(qt$BETA),as.vector(qt$DELTA)))

#############	
	print(paste('#### ',date(),'##re=',re,' has done ##PGKM takes',pgckmtime,'secs##'))
	rm(.Random.seed)

	save(FPS_PGKM,TPS_PGKM,LH_PGKM,AUC_PGKM,Cs_PGKM,RES,file=paste(Hdir,nam,'-res.RData',sep=''))
}
sink()  

print('finished'!)
