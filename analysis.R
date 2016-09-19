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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0), cex.axis=1.5);
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.01)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
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
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P00_F00_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0), cex.axis=1.5);
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-100, y=3.0, label="A", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P00_F00_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",yaxt="n",
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-100, y=3.0, label="B", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P02_F00_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0), cex.axis=1.5);
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-100, y=3.0, label="C", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P00_F02_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",yaxt="n",
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-100, y=3.0, label="D", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P02_F00_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0), cex.axis=1.5);
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.00)),cex=1, pos=4);
text(x=-100, y=3.0, label="E", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P00_F02_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",yaxt="n",
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.00)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-100, y=3.0, label="F", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M00_P02_F02_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-20,-10,0), cex.axis=1.5);
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-100, y=3.0, label="G", cex=2, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="results/fig2/evo_M02_P02_F02_a.txt", header=FALSE);
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-24,4),col="blue",xaxt="n",yaxt="n",
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
text(x=0, y=-20.0, label=expression(paste(c[P]==0.02)),cex=1, pos=4);
text(x=0, y=-22.0, label=expression(paste(c[F]==0.02)),cex=1, pos=4);
text(x=-100, y=3.0, label="H", cex=2, pos=4);
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








