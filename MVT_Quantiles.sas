submit / R;

library(mvtnorm);
corrt <- as.matrix(r);
vn <- as.numeric(v)
prob <- as.numeric(prob)
colnames(corrt) <- NULL


t2 <- qmvt(p=prob,df=vn, delta=0,corr=corrt, tail = "both");

t2data <- data.frame(t2)
endsubmit;
