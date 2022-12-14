##############################################################
###   implementing PGKM procedure of Cox model in simulation #
##############################################################
source('funs.R')
library(BB)
library(glmnet)
library(survivalROC)
library(survC1)
datid = 11
P = 200
PT = 5
Q = 200  
QT = 3
N = 100
Nvalid = 300  
Ntest = 100

lambda_1_seq <-2^seq(-5,-1,length.out = 3)
lambda_2_seq <-2^seq(-10,-1,length.out = 10)
lambda_3_seq <-2^seq(-5,-1,length.out = 3)
lambda_lasso <-2^seq(-6,-1,length.out = 10)

up.delta <- up.delta.spg

## repetition
nre <- 200
nam <- paste(datid,format(Sys.time(), "%b%d"),N,P,Q,PT,QT,nre,sep=';')
RES <- matrix(0,ncol=(4+P+Q+3))[-1,] 
RES_lasso <- matrix(0,ncol=(1+P+Q))[-1,] 
COUBET <- matrix(0,ncol=3)[-1,]   ## COU for BETA 
COU <- matrix(0,ncol=3)[-1,]      ## COU for DELTA in PGKM 
COU_lasso <- matrix(0,ncol=3)[-1,] 
pgckmTime <- matrix(0,ncol=1)[-1,] ## PGKM takes time
pgckmTime_la <- matrix(0,ncol=1)[-1,]
LH<-matrix(0,ncol=2)[-1,]
AUC<-matrix(0,ncol=2)[-1,]
Cs<-matrix(0,ncol=2)[-1,]
FPS_PGKM<-c()
TPS_PGKM<-c()
FPS_LASSO<-c()
TPS_LASSO<-c()

for(re in 1:nre){
  set.seed(re)
  dat <- datf(datid,P,Q,PT,QT,N,Nvalid,Ntest)
  X <- dat$X; Z <- dat$Z; TAU<-dat$TAU; TIM<-dat$TIM; fx<-dat$fx; H<-dat$H;
  XT <- dat$XT; ZT<-dat$ZT; HT<-dat$HT; fxT<-dat$fxT;TAUT<-dat$TAUT; TIMT<-dat$TIMT;
  Xvalid<-dat$Xvalid; Zvalid<-dat$Zvalid; Hvalid <- dat$Hvalid;fxvalid<-dat$fxvalid;
  TAUvalid<-dat$TAUvalid; TIMvalid<-dat$TIMvalid;
  XL<-dat$XL;XLT <- dat$XLT;XLvalid<-dat$XLvalid;CR<-dat$CR;
  #####PGLKM-cox method
  pgckmTime1 <- Sys.time()
  qt <- lamf(X,Z,TAU,TIM,Xvalid,Zvalid,Hvalid,TAUvalid,TIMvalid,
             lambda_1_seq,lambda_2_seq,lambda_3_seq)
  pgckmTime2 <- Sys.time()
  pgckmtime <- round(difftime(pgckmTime2,pgckmTime1,units='secs'))
  pgckmTime <- rbind(pgckmTime,as.vector(pgckmtime))
  RT<-Rf(TIMT)$R;lenT<-Rf(TIMT)$len
  hath <- hathaf(XT,Z,ZT,RT,lenT,as.vector(qt$ALPHA),as.vector(qt$BETA),as.vector(qt$DELTA))
  #likelihood
  hat_h <- hath$hat_h 
  hat_A <- hath$hat_A 
  LH1 <-(1/Ntest)* (t(TAUT) %*% ( XT%*%as.vector(qt$BETA)+hat_h-hat_A ))
  #AUC
  cutoff<-max(TIMT)*0.9
  mark<-XT%*%as.vector(qt$BETA)+hat_h
  PGKM_ROC<-survivalROC(Stime=TIMT,status = TAUT,marker = mark,predict.time = cutoff,span=0.25*100^(-0.2))
  PGKM_AUC<-PGKM_ROC$AUC
  FPS_PGKM<-rbind(FPS_PGKM,PGKM_ROC$FP)
  TPS_PGKM<-rbind(TPS_PGKM,PGKM_ROC$TP)
  #C-statistics
  ta<-quantile(TIMT,0.7)
  mydata<-cbind(TIMT,TAUT,mark)
  PGKM_C<-Est.Cval(mydata, ta, nofit=TRUE)$Dhat
  
  COUBET <- rbind(COUBET, cou(as.vector(qt$BETA),PT))
  COU <- rbind(COU, cou(as.vector(qt$DELTA),QT))
  ## PGKM-cox result	  
  RES <- rbind(RES, c(qt$ITER,as.vector(qt$LAMBDA),as.vector(qt$BETA),as.vector(qt$DELTA),as.vector(qt$REG)))
  
  #####lasso-cox method
  pgckmTime3 <- Sys.time()
  qt_lasso <- lamf_lasso(XL,TIM,TAU,XLvalid,TIMvalid,TAUvalid,lambda_lasso)
  pgckmTime4 <- Sys.time()
  pgckmtime_la <- round(difftime(pgckmTime3,pgckmTime4,units='secs'))
  pgckmTime_la <- rbind(pgckmTime_la,as.vector(pgckmtime_la))
  RT_la<-Rf(TIMT)$R
  lenT_la<-Rf(TIMT)$len
  hath_la <- hata_lasso(XLT,RT_la,lenT_la,as.vector(qt_lasso$BETA))
  hat_A_la <- hath_la$hat_A 
  #likelihood
  LH_lasso <-(1/Ntest)* (t(TAUT) %*% (XLT%*%as.vector(qt_lasso$BETA)-hat_A_la))
  #AUC
  mark_la<-XLT%*%as.vector(qt_lasso$BETA)
  LASSO_ROC<-survivalROC(Stime=TIMT,status = TAUT,marker = mark_la,predict.time = cutoff,span=0.25*100^(-0.2))
  LASSO_AUC<-LASSO_ROC$AUC
  FPS_LASSO<-rbind(FPS_LASSO,LASSO_ROC$FP)
  TPS_LASSO<-rbind(TPS_LASSO,LASSO_ROC$TP)
  #C-statistics
  mydat<-cbind(TIMT,TAUT,mark_la)
  LASSO_C<-Est.Cval(mydat, ta, nofit=TRUE)$Dhat
  
  COU_lasso <- rbind(COU_lasso, coul(as.vector(qt_lasso$BETA),PT,P,QT,Q))
  ## lasso-cox result	  
  RES_lasso <- rbind(RES_lasso, c(as.vector(qt_lasso$LAMBDA),as.vector(qt_lasso$BETA)))
  
  LH <- rbind(LH,c(LH1,LH_lasso)) 
  AUC<- rbind(AUC,c(PGKM_AUC,LASSO_AUC))
  Cs<- rbind(Cs,c(PGKM_C,LASSO_C))
  COUU <- cbind(COU,COUBET)
  print(paste('#### ',date(),'##re=',re,' has done ##PGKM takes',pgckmtime,'secs##','lasso-cox takes',pgckmtime_la,'secs####'))
  rm(.Random.seed)
  
  print('++AUC++')
  print(AUC)
  print("++Cs++")
  print(Cs)
  
  save(FPS_PGKM,TPS_PGKM,FPS_LASSO,TPS_LASSO,LH,AUC,Cs,RES,RES_lasso,COUU,COU_lasso,pgckmTime,pgckmTime_la,file=paste(nam,'.RData',sep=''))
}