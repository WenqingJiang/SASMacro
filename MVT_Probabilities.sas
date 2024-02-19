submit / R;

library(mvtnorm);
corrt <- as.matrix(r);
vn <- as.numeric(v)
ti <- as.numeric(ti)
colnames(corrt) <- NULL


t2p <- pmvt(lower=-ti,upper=ti,df=vn, delta=0,corr=corrt);

t2pdata <- data.frame(t2p)
endsubmit;
