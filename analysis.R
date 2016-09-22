# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
#  ____                                 _                         #
# |  _ \ _ __ ___       _ __   ___  ___| |_    ___ ___  _ __      #
# | |_) | '__/ _ \_____| '_ \ / _ \/ __| __|  / __/ _ \| '_ \     #
# |  __/| | |  __/_____| |_) | (_) \__ \ |_  | (_| (_) | |_) |    #
# |_|   |_|  \___|     | .__/ \___/|___/\__|  \___\___/| .__/     #
#                      |_|                             |_|        #
#                                                                 #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

rm(list=ls()); 
setwd("~/Dropbox/DuthieManu/pre_post_cop");
gens      <- 39999;
beta_val <- 3;

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX      PLOT FIGURE 1  XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
setEPS(); # postscript below for final publication?
cairo_ps("figures/M00_P01_F00.eps",family="Arial",height=6,width=6);
par(oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig1/evo_M00_P01_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.01 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0,10), cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.00)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.01)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
rm(evo_cc_02);
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX      PLOT FIGURE 2  XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

setEPS(); # postscript below for final publication?
cairo_ps("figures/evo_combs.eps",family="Arial",height=8,width=5.5);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0,10), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.00)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-300, y=11.0, label="A", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,]; 
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.02)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-300, y=11.0, label="B", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0,10), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-18.0, label=expression(paste(c[M]==0.00)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-300, y=11.0, label="C", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val,]; 
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.00)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-300, y=11.0, label="D", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,]; 
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0,10), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.02)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-300, y=11.0, label="E", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val,]; 
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.02)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-300, y=11.0, label="F", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val,]; 
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0,10), cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-18.0, label=expression(paste(c[M]==0.00)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-300, y=11.0, label="G", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,14),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==gens));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==gens));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==gens));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==gens));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==gens));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==gens));
polygon(y=c(lPO,rev(uPO)),x=c(1:(gens+1),(gens+1):1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:(gens+1),(gens+1):1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-18.0, label=expression(paste(c[M]==0.02)),cex=1, pos=4);
text(x=0, y=-21.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-24.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-300, y=11.0, label="H", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
dev.off();
#==========================================================================

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #





















# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX SUPPORTING INFORMATION: POLYANDRY PHENOTYPES    XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #


rm(list=ls()); 
setwd("~/Dropbox/DuthieManu/pre_post_cop");
gens      <- 39999;
beta_val  <- 3;

evo_cc_02 <- read.table(file="results/fig2/evo_M00_P02_F00_b.txt", header=FALSE);
evgens <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evlast <- evgens[evgens[,6]==gens,]; 
above1 <- which(evlast[,8] == min(evlast[(evlast[,8]-1)>0,8]))[1];
below1 <- which(evlast[,8] == max(evlast[(evlast[,8]-0)<0,8]))[1];
simabo <- evlast[above1,1];
simbel <- evlast[below1,1];

ped_cc_02 <- read.table(file="results/fig2/ped_M00_P02_F00_b.txt", header=FALSE);
xxx <- 0:gens;

setEPS(); # postscript below for final publication?
cairo_ps("figures/Pp_vals.eps",family="Arial",height=11,width=9);
par(mfrow=c(2,1),mar=c(5,5,1,3));
evo_bel <- evo_cc_02[evo_cc_02[,1]==simbel,8];
plot(x=xxx,y=evo_bel,type="l",lwd=2,ylim=c(-4.2,4.2),col="red",xaxt="n",
     xlab=expression(paste("Generation")),
     ylab=expression(paste("Mean allele value")),
     cex.lab=1.75,cex.axis=1.5,yaxt="n",lty="solid");
abline(h=0,lty="dotted",lwd=0.8);
arrows(x0=gens+5000,x1=gens+250,y0=evo_bel[gens],y1=evo_bel[gens],length=0.1,lwd=2);
evo_bel_end <- as.character(round(x=evo_bel[gens],digits=3))
axis(side=4,at=evo_bel[gens],labels=evo_bel_end,lwd=2,cex.axis=1.25);
axis(side=2,at=c(-4,-2,0,2,4), cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-1600,y=4.1,pos=4,labels="A",cex=1.5,col="black");

pedbe     <- ped_cc_02[ped_cc_02[,1]==simbel,];
pedbe     <- pedbe[pedbe[,3]==0 & pedbe[,4]==1 & pedbe[,8] > 0,]; # Living females
hist(pedbe[,11],breaks=10,main="",cex.lab=1.75,yaxs="i",xaxs="i",ylim=c(0,22),
     xlab=expression(paste("Tendency for polyandry phenotype (",P[p],")")),
     col="red",cex.axis=1.5);
polygon(x=c(0:10,10:0),y=c(rep(0,11),rep(22,11)),col="grey30");
hist(pedbe[,11],breaks=10,main="",cex.lab=1.5,yaxs="i",xaxs="i",ylim=c(0,22),
     xlab=expression(paste("Tendency for polyandry phenotype (",P[p],")")),
     col="red",add=TRUE);
box();
abline(v=0,lwd=3);
text(x=0,y=21,pos=4,labels="Polyandrous",cex=1.5,col="white");
text(x=-16,y=21,pos=4,labels="B; Monandrous",cex=1.5,col="black");
arrows(x0=mean(pedbe[,11]),x1=mean(pedbe[,11]),y1=0,y0=2,length=0.1,lwd=2)
evo_bel_phn <- as.character(round(x=mean(pedbe[,11]),digits=3));
text(x=mean(pedbe[,11])+0.4,y=2.5,labels=evo_bel_phn,cex=0.8);
dev.off();

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

gens     <- 39999;
beta_val <- 3;
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX SUPPORTING INFORMATION FOR P00 M00 F00 SIMULATIONS  XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,]; 
c_WP_PO <- rep(x=0, times=gens);
c_WP_EP <- rep(x=0, times=gens);
c_PO_EP <- rep(x=0, times=gens);
for(i in 0:gens){ # XXX XXX  NOTE: THIS LOOP TAKES A LONG TIME TO RUN XXX #
    use          <- evo[evo[,6]==i,];
    c_WP_PO[i+1] <- cor(use[,7],use[,8]);
    c_WP_EP[i+1] <- cor(use[,7],use[,9]);
    c_PO_EP[i+1] <- cor(use[,8],use[,9]);
    if(i %% 1000 == 0){
        print(paste("Generation ",i));
    }
}

xxx  <- 0:gens;
sims <- unique(evo_cc_02[,1]);
# ------------------------------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("figures/SI_M00_P00_E00.eps",family="Arial",height=7,width=7.5);
par(mfrow=c(2,2),oma=c(5,5,1,4), mar=c(0.5,0.5,0.5,0.5),lwd=1);
# ----------------------------------------------------------------
mPo <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mPo,type="l",lwd=2,ylim=c(-50,30),col="red",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,8];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="thistle");
}
points(x=xxx,y=mPo,type="l",lwd=2,col="red");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
text(x=-900, y=27.0, label="A", cex=2, pos=4);
#text(x=0, y=-35.0, label=expression(paste(c[M]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-42.0, label=expression(paste(c[P]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-49.0, label=expression(paste(c[F]==0.00)),cex=1.2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[P]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-50,30),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,7];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="lightblue");
}
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=-900, y=27.0, label="B", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[M]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mEP,type="l",lwd=2,ylim=c(-50,30),col="black",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,9];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="grey70");
}
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-900, y=27.0, label="C", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[F]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
par(lwd=2);
plot(x=xxx,y=c_WP_PO,type="l",lwd=1,ylim=c(-1,1),col="purple",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
axis(side=4,at=c(-1,-0.5,0,0.5,1), cex.axis=1.5);
points(x=xxx,y=c_WP_EP,type="l",lwd=1,col="royalblue4",lty="solid");
points(x=xxx,y=c_PO_EP,type="l",lwd=1,col="red4",lty="solid");
abline(h=0,lwd="0.8",lty="dotted");
text(x=-300, y=0.88, label="D", cex=2, pos=4);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
mtext(expression(paste("Mean allele value correlation\t\t\t\t\t\t")),
	outer=TRUE,side=4,line=2.5,cex=1.5);
# ------------------------------------------------------------------------#
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX SUPPORTING INFORMATION FOR P00 M02 F00 SIMULATIONS  XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,]; 
c_WP_PO <- rep(x=0, times=gens);
c_WP_EP <- rep(x=0, times=gens);
c_PO_EP <- rep(x=0, times=gens);
for(i in 0:gens){ # XXX XXX  NOTE: THIS LOOP TAKES A LONG TIME TO RUN XXX #
    use          <- evo[evo[,6]==i,];
    c_WP_PO[i+1] <- cor(use[,7],use[,8]);
    c_WP_EP[i+1] <- cor(use[,7],use[,9]);
    c_PO_EP[i+1] <- cor(use[,8],use[,9]);
    if(i %% 1000 == 0){
        print(paste("Generation ",i));
    }
}

xxx  <- 0:gens;
sims <- unique(evo_cc_02[,1]);
# ------------------------------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("figures/SI_M02_P00_E00.eps",family="Arial",height=7,width=7.5);
par(mfrow=c(2,2),oma=c(5,5,1,4), mar=c(0.5,0.5,0.5,0.5),lwd=1);
# ----------------------------------------------------------------
mPo <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mPo,type="l",lwd=2,ylim=c(-50,30),col="red",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,8];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="thistle");
}
points(x=xxx,y=mPo,type="l",lwd=2,col="red");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
text(x=-900, y=27.0, label="A", cex=2, pos=4);
#text(x=0, y=-35.0, label=expression(paste(c[M]==0.02)),cex=1.2, pos=4);
#text(x=0, y=-42.0, label=expression(paste(c[P]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-49.0, label=expression(paste(c[F]==0.00)),cex=1.2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[P]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-50,30),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,7];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="lightblue");
}
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=-900, y=27.0, label="B", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[M]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mEP,type="l",lwd=2,ylim=c(-50,30),col="black",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,9];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="grey70");
}
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-900, y=27.0, label="C", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[F]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
par(lwd=2);
plot(x=xxx,y=c_WP_PO,type="l",lwd=1,ylim=c(-1,1),col="purple",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
axis(side=4,at=c(-1,-0.5,0,0.5,1), cex.axis=1.5);
points(x=xxx,y=c_WP_EP,type="l",lwd=1,col="royalblue4",lty="solid");
points(x=xxx,y=c_PO_EP,type="l",lwd=1,col="red4",lty="solid");
abline(h=0,lwd="0.8",lty="dotted");
text(x=-300, y=0.88, label="D", cex=2, pos=4);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
mtext(expression(paste("Mean allele value correlation\t\t\t\t\t\t")),
	outer=TRUE,side=4,line=2.5,cex=1.5);
# ------------------------------------------------------------------------#
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX SUPPORTING INFORMATION FOR P02 M00 F00 SIMULATIONS  XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.00 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val,]; 
c_WP_PO <- rep(x=0, times=gens);
c_WP_EP <- rep(x=0, times=gens);
c_PO_EP <- rep(x=0, times=gens);
for(i in 0:gens){ # XXX XXX  NOTE: THIS LOOP TAKES A LONG TIME TO RUN XXX #
    use          <- evo[evo[,6]==i,];
    c_WP_PO[i+1] <- cor(use[,7],use[,8]);
    c_WP_EP[i+1] <- cor(use[,7],use[,9]);
    c_PO_EP[i+1] <- cor(use[,8],use[,9]);
    if(i %% 1000 == 0){
        print(paste("Generation ",i));
    }
}

xxx  <- 0:gens;
sims <- unique(evo_cc_02[,1]);
# ------------------------------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("figures/SI_M00_P02_E00.eps",family="Arial",height=7,width=7.5);
par(mfrow=c(2,2),oma=c(5,5,1,4), mar=c(0.5,0.5,0.5,0.5),lwd=1);
# ----------------------------------------------------------------
mPo <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mPo,type="l",lwd=2,ylim=c(-50,30),col="red",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,8];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="thistle");
}
points(x=xxx,y=mPo,type="l",lwd=2,col="red");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
text(x=-900, y=27.0, label="A", cex=2, pos=4);
#text(x=0, y=-35.0, label=expression(paste(c[M]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-42.0, label=expression(paste(c[P]==0.02)),cex=1.2, pos=4);
#text(x=0, y=-49.0, label=expression(paste(c[F]==0.00)),cex=1.2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[P]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-50,30),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,7];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="lightblue");
}
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=-900, y=27.0, label="B", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[M]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mEP,type="l",lwd=2,ylim=c(-50,30),col="black",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,9];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="grey70");
}
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-900, y=27.0, label="C", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[F]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
par(lwd=2);
plot(x=xxx,y=c_WP_PO,type="l",lwd=1,ylim=c(-1,1),col="purple",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
axis(side=4,at=c(-1,-0.5,0,0.5,1), cex.axis=1.5);
points(x=xxx,y=c_WP_EP,type="l",lwd=1,col="royalblue4",lty="solid");
points(x=xxx,y=c_PO_EP,type="l",lwd=1,col="red4",lty="solid");
abline(h=0,lwd="0.8",lty="dotted");
text(x=-300, y=0.88, label="D", cex=2, pos=4);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
mtext(expression(paste("Mean allele value correlation\t\t\t\t\t\t")),
	outer=TRUE,side=4,line=2.5,cex=1.5);
# ------------------------------------------------------------------------#
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX SUPPORTING INFORMATION FOR P00 M00 F02 SIMULATIONS  XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.00 & evo_cc_02[,4]==0.00 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val,]; 
c_WP_PO <- rep(x=0, times=gens);
c_WP_EP <- rep(x=0, times=gens);
c_PO_EP <- rep(x=0, times=gens);
for(i in 0:gens){ # XXX XXX  NOTE: THIS LOOP TAKES A LONG TIME TO RUN XXX #
    use          <- evo[evo[,6]==i,];
    c_WP_PO[i+1] <- cor(use[,7],use[,8]);
    c_WP_EP[i+1] <- cor(use[,7],use[,9]);
    c_PO_EP[i+1] <- cor(use[,8],use[,9]);
    if(i %% 1000 == 0){
        print(paste("Generation ",i));
    }
}

xxx  <- 0:gens;
sims <- unique(evo_cc_02[,1]);
# ------------------------------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("figures/SI_M00_P00_E02.eps",family="Arial",height=7,width=7.5);
par(mfrow=c(2,2),oma=c(5,5,1,4), mar=c(0.5,0.5,0.5,0.5),lwd=1);
# ----------------------------------------------------------------
mPo <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mPo,type="l",lwd=2,ylim=c(-50,30),col="red",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,8];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="thistle");
}
points(x=xxx,y=mPo,type="l",lwd=2,col="red");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
text(x=-900, y=27.0, label="A", cex=2, pos=4);
#text(x=0, y=-35.0, label=expression(paste(c[M]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-42.0, label=expression(paste(c[P]==0.02)),cex=1.2, pos=4);
#text(x=0, y=-49.0, label=expression(paste(c[F]==0.00)),cex=1.2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[P]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-50,30),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,7];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="lightblue");
}
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=-900, y=27.0, label="B", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[M]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mEP,type="l",lwd=2,ylim=c(-50,30),col="black",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,9];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="grey70");
}
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-900, y=27.0, label="C", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[F]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
par(lwd=2);
plot(x=xxx,y=c_WP_PO,type="l",lwd=1,ylim=c(-1,1),col="purple",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
axis(side=4,at=c(-1,-0.5,0,0.5,1), cex.axis=1.5);
points(x=xxx,y=c_WP_EP,type="l",lwd=1,col="royalblue4",lty="solid");
points(x=xxx,y=c_PO_EP,type="l",lwd=1,col="red4",lty="solid");
abline(h=0,lwd="0.8",lty="dotted");
text(x=-300, y=0.88, label="D", cex=2, pos=4);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
mtext(expression(paste("Mean allele value correlation\t\t\t\t\t\t")),
	outer=TRUE,side=4,line=2.5,cex=1.5);
# ------------------------------------------------------------------------#
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX SUPPORTING INFORMATION FOR P02 M02 F00 SIMULATIONS  XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.00,];
evo <- evB[evB[,2]==beta_val,]; 
c_WP_PO <- rep(x=0, times=gens);
c_WP_EP <- rep(x=0, times=gens);
c_PO_EP <- rep(x=0, times=gens);
for(i in 0:gens){ # XXX XXX  NOTE: THIS LOOP TAKES A LONG TIME TO RUN XXX #
    use          <- evo[evo[,6]==i,];
    c_WP_PO[i+1] <- cor(use[,7],use[,8]);
    c_WP_EP[i+1] <- cor(use[,7],use[,9]);
    c_PO_EP[i+1] <- cor(use[,8],use[,9]);
    if(i %% 1000 == 0){
        print(paste("Generation ",i));
    }
}

xxx  <- 0:gens;
sims <- unique(evo_cc_02[,1]);
# ------------------------------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("figures/SI_M02_P02_E00.eps",family="Arial",height=7,width=7.5);
par(mfrow=c(2,2),oma=c(5,5,1,4), mar=c(0.5,0.5,0.5,0.5),lwd=1);
# ----------------------------------------------------------------
mPo <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mPo,type="l",lwd=2,ylim=c(-50,30),col="red",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,8];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="thistle");
}
points(x=xxx,y=mPo,type="l",lwd=2,col="red");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
text(x=-900, y=27.0, label="A", cex=2, pos=4);
#text(x=0, y=-35.0, label=expression(paste(c[M]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-42.0, label=expression(paste(c[P]==0.02)),cex=1.2, pos=4);
#text(x=0, y=-49.0, label=expression(paste(c[F]==0.00)),cex=1.2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[P]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-50,30),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,7];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="lightblue");
}
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=-900, y=27.0, label="B", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[M]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mEP,type="l",lwd=2,ylim=c(-50,30),col="black",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,9];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="grey70");
}
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-900, y=27.0, label="C", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[F]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
par(lwd=2);
plot(x=xxx,y=c_WP_PO,type="l",lwd=1,ylim=c(-1,1),col="purple",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
axis(side=4,at=c(-1,-0.5,0,0.5,1), cex.axis=1.5);
points(x=xxx,y=c_WP_EP,type="l",lwd=1,col="royalblue4",lty="solid");
points(x=xxx,y=c_PO_EP,type="l",lwd=1,col="red4",lty="solid");
abline(h=0,lwd="0.8",lty="dotted");
text(x=-300, y=0.88, label="D", cex=2, pos=4);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
mtext(expression(paste("Mean allele value correlation\t\t\t\t\t\t")),
	outer=TRUE,side=4,line=2.5,cex=1.5);
# ------------------------------------------------------------------------#
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #



# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX SUPPORTING INFORMATION FOR P02 M00 F02 SIMULATIONS  XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.00 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val,]; 
c_WP_PO <- rep(x=0, times=gens);
c_WP_EP <- rep(x=0, times=gens);
c_PO_EP <- rep(x=0, times=gens);
for(i in 0:gens){ # XXX XXX  NOTE: THIS LOOP TAKES A LONG TIME TO RUN XXX #
    use          <- evo[evo[,6]==i,];
    c_WP_PO[i+1] <- cor(use[,7],use[,8]);
    c_WP_EP[i+1] <- cor(use[,7],use[,9]);
    c_PO_EP[i+1] <- cor(use[,8],use[,9]);
    if(i %% 1000 == 0){
        print(paste("Generation ",i));
    }
}

xxx  <- 0:gens;
sims <- unique(evo_cc_02[,1]);
# ------------------------------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("figures/SI_M02_P00_E02.eps",family="Arial",height=7,width=7.5);
par(mfrow=c(2,2),oma=c(5,5,1,4), mar=c(0.5,0.5,0.5,0.5),lwd=1);
# ----------------------------------------------------------------
mPo <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mPo,type="l",lwd=2,ylim=c(-50,30),col="red",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,8];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="thistle");
}
points(x=xxx,y=mPo,type="l",lwd=2,col="red");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
text(x=-900, y=27.0, label="A", cex=2, pos=4);
#text(x=0, y=-35.0, label=expression(paste(c[M]==0.00)),cex=1.2, pos=4);
#text(x=0, y=-42.0, label=expression(paste(c[P]==0.02)),cex=1.2, pos=4);
#text(x=0, y=-49.0, label=expression(paste(c[F]==0.00)),cex=1.2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[P]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-50,30),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,7];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="lightblue");
}
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
abline(h=0,lwd=0.8,lty="dotted");
text(x=-900, y=27.0, label="B", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[M]==0.00)),cex=1.5, pos=4);
# ----------------------------------------------------------------
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
plot(x=xxx,y=mEP,type="l",lwd=2,ylim=c(-50,30),col="black",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
for(i in sims){
 iter <- evo[evo[,1]==i,9];
 points(x=xxx,y=iter,type="l",lwd=0.5,col="grey70");
}
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-40,-20,0,20),cex.axis=1.5);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
text(x=-900, y=27.0, label="C", cex=2, pos=4);
text(x=0, y=-47.0, label=expression(paste(c[F]==0.02)),cex=1.5, pos=4);
# ----------------------------------------------------------------
par(lwd=2);
plot(x=xxx,y=c_WP_PO,type="l",lwd=1,ylim=c(-1,1),col="purple",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
axis(side=4,at=c(-1,-0.5,0,0.5,1), cex.axis=1.5);
points(x=xxx,y=c_WP_EP,type="l",lwd=1,col="royalblue4",lty="solid");
points(x=xxx,y=c_PO_EP,type="l",lwd=1,col="red4",lty="solid");
abline(h=0,lwd="0.8",lty="dotted");
text(x=-300, y=0.88, label="D", cex=2, pos=4);
axis(side=1,at=c(0,10000,20000,30000), cex.axis=1.5);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
mtext(expression(paste("Mean allele value correlation\t\t\t\t\t\t")),
	outer=TRUE,side=4,line=2.5,cex=1.5);
# ------------------------------------------------------------------------#
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #


