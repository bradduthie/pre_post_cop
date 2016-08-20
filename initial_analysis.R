
rm(list=ls()); 
library(scatterplot3d);
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX Multi-panel figure of traits over genreations   XXX XXX XXX #
# XXX XXX XXX with fixed polyandry costs, random other costs  XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

evo_cc_02 <- read.table(file="evo_cost_combs02.txt", header=FALSE);

# ------------------------------------------------------------------------#
# -----------  BETA = 0  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 0;
setEPS(); # postscript below for final publication?
cairo_ps("evo_combs02_B0.eps",family="Arial",height=8,width=5.5);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
#polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
#legend("topleft", inset=0.025, bty="n", cex=0.95, 
#       title=expression(paste("Evolving trait (",beta==0,")")),
#      c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
#       fill=c("blue", "red", "black"));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="A", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="B", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="C", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="D", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="E", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="F", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="G", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="H", cex=2, pos=4);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
dev.off();
#==========================================================================

# ------------------------------------------------------------------------#
# -----------  BETA = 1  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 1;
setEPS(); # postscript below for final publication?
cairo_ps("evo_combs02_B1.eps",family="Arial",height=8,width=5.5);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
#polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
#legend("topleft", inset=0.025, bty="n", cex=0.95, 
#       title=expression(paste("Evolving trait (",beta==0,")")),
#      c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
#       fill=c("blue", "red", "black"));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="A", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="B", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="C", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="D", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="E", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="F", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="G", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="H", cex=2, pos=4);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
dev.off();
#==========================================================================


# ------------------------------------------------------------------------#
# -----------  BETA = 2  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 2;
setEPS(); # postscript below for final publication?
cairo_ps("evo_combs02_B2.eps",family="Arial",height=8,width=5.5);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
#polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
#legend("topleft", inset=0.025, bty="n", cex=0.95, 
#       title=expression(paste("Evolving trait (",beta==0,")")),
#      c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
#       fill=c("blue", "red", "black"));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="A", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="B", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="C", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="D", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="E", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="F", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="G", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="H", cex=2, pos=4);
# ------------------------------------------------------------------------#
mtext(expression(paste("Generation")),
	outer=TRUE,side=1,line=3.0,cex=1.5);
mtext(expression(paste("Mean allele value")),
	outer=TRUE,side=2,line=2.5,cex=1.5);
dev.off();
#==========================================================================

# ------------------------------------------------------------------------#
# -----------  BETA = 3  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 3;
setEPS(); # postscript below for final publication?
cairo_ps("evo_combs02_B3.eps",family="Arial",height=8,width=5.5);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
#polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
#legend("topleft", inset=0.025, bty="n", cex=0.95, 
#       title=expression(paste("Evolving trait (",beta==0,")")),
#      c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
#       fill=c("blue", "red", "black"));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="A", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="B", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="C", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="D", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0",cex=1, pos=4);
text(x=-100, y=3.0, label="E", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="F", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",
     xlab="",ylab="Allele value",cex.lab=2,cex.axis=1.5,yaxt="n");
axis(side=2,at=c(-8,-4,0,4), cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
text(x=0, y=-6.0, label="Prec-cost = 0",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="G", cex=2, pos=4);
# ------------------------------------------------------------------------#
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
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5);
axis(side=1,at=c(0,2000,4000), cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
text(x=0, y=-6.0, label="Prec-cost = 0.02",cex=1, pos=4);
text(x=0, y=-7.25, label="Poly-cost = 0.02",cex=1, pos=4);
text(x=0, y=-8.5, label="Post-cost = 0.02",cex=1, pos=4);
text(x=-100, y=3.0, label="H", cex=2, pos=4);
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

rm(list=ls()); 
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX Multi-panel figure of randomised inbreeding     XXX XXX XXX #
# XXX XXX XXX costs with fixed polyandry costs                XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

# Random inbreeding traits with a fixed pcost of zero
evo_rand_pc_00 <- read.table(file="evo_rand_pc_00.txt", header=FALSE);
ela_rand_pc_00 <- evo_rand_pc_00[evo[,6]==4999,];
write.table(x=ela_rand_pc_00, file="ela_rand_pc_00.txt", col.names=FALSE, row.names=FALSE);

# -----------  BETA = 0  -------------------------------------------------#
ela_rand_pc_00 <- read.table(file="ela_rand_pc_00.txt", header=FALSE);

elaB0 <- ela_rand_pc_00[ela_rand_pc_00[,2]==0,];
resB0 <- rep(x=0, times=dim(elaB0)[1]);
for(i in 1:dim(elaB0)[1]){ # Doing this as a loop just to show the logic
   if(elaB0[i,7] >= 1 & elaB0[i,9] >= 1){
       resB0[i] <- "purple"; # Pre and post-copulatory evolved
   }
   if(elaB0[i,7] >= 1 & elaB0[i,9] < 1){
       resB0[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB0[i,7] < 1 & elaB0[i,9] >= 1){
       resB0[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB0[i,7] < 1 & elaB0[i,9] < 1){
       resB0[i] <- "grey80"; # Neither pre nor post evolved
   }
}
polB0 <- rep(x=0, times=dim(elaB0)[1]);
for(i in 1:dim(elaB0)[1]){ # Doing this as a loop just to show the logic
   if(elaB0[i,8] >= 1){
       polB0[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       polB0[i] <- "black"; # Pre but not post-copulatory evolved
   }
}
plot(x=elaB0[,3], y=elaB0[,5], cex=1.9, lwd=2, col=polB0, pch=21, bg=resB0,
     xlab="Cost of pre-cop inbreeding preference",
     ylab="Cost of post-cop inbreeding preference");



# -----------  BETA = 1  -------------------------------------------------#
elaB1 <- ela_rand_pc_00[ela_rand_pc_00[,2]==1,];
resB1 <- rep(x=0, times=dim(elaB1)[1]);
for(i in 1:dim(elaB1)[1]){ # Doing this as a loop just to show the logic
   if(elaB1[i,7] <= -1 & elaB1[i,9] <= -1){
       resB1[i] <- "purple"; # Pre and post-copulatory evolved
   }
   if(elaB1[i,7] <= -1 & elaB1[i,9] > -1){
       resB1[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB1[i,7] > -1 & elaB1[i,9] <= -1){
       resB1[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB1[i,7] > -1 & elaB1[i,9] > -1){
       resB1[i] <- "grey80"; # Neither pre nor post evolved
   }
}
polB1 <- rep(x=0, times=dim(elaB1)[1]);
for(i in 1:dim(elaB1)[1]){ # Doing this as a loop just to show the logic
   if(elaB1[i,8] >= 1){
       polB1[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       polB1[i] <- "black"; # Pre but not post-copulatory evolved
   }
}
plot(x=elaB1[,3], y=elaB1[,5], cex=1.9, lwd=2, col=polB1, pch=21, bg=resB1,
     xlab="Cost of pre-cop inbreeding preference",
     ylab="Cost of post-cop inbreeding preference");






# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX  OLDER ANALYSES, PROBABLY NOT TO BE USED IN THE MANUSCRIPT  XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #



#==========================================================================

rm(list=ls()); 
                                                                     
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

evo  <- read.table(file="evo_zallele_0.txt", header=FALSE);
ela  <- evo[evo[,6]==4999,];
write.table(x=ela, file="ela_zallele_0.txt", col.names=FALSE, row.names=FALSE);

#==========================================================================


rm(list=ls()); 
library(scatterplot3d);
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

ela <- read.table(file="ela_zallele_0.txt", header=FALSE);

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX Three dimensional scatter plot  XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# First, let's look at the results for Beta = 0 (no inbreeding depression)
# Note, we're using past modelling to assume that inbreeding preference is
# favoured in the absence of inbreeding depression.
# For now, we also have an arbitrary threshold of 1 for indicating selection
# for or against a trait -- might find a less arbitrary way to do this,
# but given the type of modelling, I don't think it's unreasonable.
# -----------  BETA = 0  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_0.eps",family="Arial",height=7,width=7);
elaB0 <- ela[ela[,2]==0,]; # Now we have juat the Beta = 0 results in ela
resB0 <- rep(0,dim(elaB0)[1]); # This vector will find one of four results
for(i in 1:dim(elaB0)[1]){ # Doing this as a loop just to show the logic
   if(elaB0[i,7] >= 1 & elaB0[i,9] >= 1){
       resB0[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB0[i,7] >= 1 & elaB0[i,9] < 1){
       resB0[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB0[i,7] < 1 & elaB0[i,9] >= 1){
       resB0[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB0[i,7] < 1 & elaB0[i,9] < 1){
       resB0[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB0[,3:5], color=resB0, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==0,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 1  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_1.eps",family="Arial",height=7,width=7);
elaB1 <- ela[ela[,2]==1,]; # Now we have juat the Beta = 0 results in ela
resB1 <- rep(0,dim(elaB1)[1]); # This vector will find one of four results
for(i in 1:dim(elaB1)[1]){ # Doing this as a loop just to show the logic
   if(elaB1[i,7] <= -1 & elaB1[i,9] <= -1){
       resB1[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB1[i,7] <= -1 & elaB1[i,9] > -1){
       resB1[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB1[i,7] > -1 & elaB1[i,9] <= -1){
       resB1[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB1[i,7] > -1 & elaB1[i,9] > -1){
       resB1[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB1[,3:5], color=resB1, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding avoidance",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding avoidance");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==1,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 2  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_2.eps",family="Arial",height=7,width=7);
elaB2 <- ela[ela[,2]==2,]; # Now we have juat the Beta = 0 results in ela
resB2 <- rep(0,dim(elaB2)[1]); # This vector will find one of four results
for(i in 1:dim(elaB2)[1]){ # Doing this as a loop just to show the logic
   if(elaB2[i,7] <= -1 & elaB2[i,9] <= -1){
       resB2[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB2[i,7] <= -1 & elaB2[i,9] > -1){
       resB2[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB2[i,7] > -1 & elaB2[i,9] <= -1){
       resB2[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB2[i,7] > -1 & elaB2[i,9] > -1){
       resB2[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB2[,3:5], color=resB2, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding avoidance",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding avoidance");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==2,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 3  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_3.eps",family="Arial",height=7,width=7);
elaB3 <- ela[ela[,2]==3,]; # Now we have juat the Beta = 0 results in ela
resB3 <- rep(0,dim(elaB3)[1]); # This vector will find one of four results
for(i in 1:dim(elaB3)[1]){ # Doing this as a loop just to show the logic
   if(elaB3[i,7] <= -1 & elaB3[i,9] <= -1){
       resB3[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB3[i,7] <= -1 & elaB3[i,9] > -1){
       resB3[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB3[i,7] > -1 & elaB3[i,9] <= -1){
       resB3[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB3[i,7] > -1 & elaB3[i,9] > -1){
       resB3[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB3[,3:5], color=resB3, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding avoidance",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding avoidance");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==3,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

evo  <- read.table(file="evo_zallele_1.txt", header=FALSE);
ela  <- evo[evo[,6]==4999,];
write.table(x=ela, file="ela_zallele_1.txt", col.names=FALSE, row.names=FALSE);

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# Now will do the same thing for random initial conditions

rm(list=ls()); 
library(scatterplot3d);
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

ela <- read.table(file="ela_zallele_1.txt", header=FALSE);

# -----------  BETA = 0  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_0_z1.eps",family="Arial",height=7,width=7);
elaB0 <- ela[ela[,2]==0,]; # Now we have juat the Beta = 0 results in ela
resB0 <- rep(0,dim(elaB0)[1]); # This vector will find one of four results
for(i in 1:dim(elaB0)[1]){ # Doing this as a loop just to show the logic
   if(elaB0[i,7] >= 1 & elaB0[i,9] >= 1){
       resB0[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB0[i,7] >= 1 & elaB0[i,9] < 1){
       resB0[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB0[i,7] < 1 & elaB0[i,9] >= 1){
       resB0[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB0[i,7] < 1 & elaB0[i,9] < 1){
       resB0[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB0[,3:5], color=resB0, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==0,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 1  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_1_z1.eps",family="Arial",height=7,width=7);
elaB1 <- ela[ela[,2]==1,]; # Now we have juat the Beta = 0 results in ela
resB1 <- rep(0,dim(elaB1)[1]); # This vector will find one of four results
for(i in 1:dim(elaB1)[1]){ # Doing this as a loop just to show the logic
   if(elaB1[i,7] <= -1 & elaB1[i,9] <= -1){
       resB1[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB1[i,7] <= -1 & elaB1[i,9] > -1){
       resB1[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB1[i,7] > -1 & elaB1[i,9] <= -1){
       resB1[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB1[i,7] > -1 & elaB1[i,9] > -1){
       resB1[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB1[,3:5], color=resB1, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding avoidance",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding avoidance");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==1,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 2  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_2_z1.eps",family="Arial",height=7,width=7);
elaB2 <- ela[ela[,2]==2,]; # Now we have juat the Beta = 0 results in ela
resB2 <- rep(0,dim(elaB2)[1]); # This vector will find one of four results
for(i in 1:dim(elaB2)[1]){ # Doing this as a loop just to show the logic
   if(elaB2[i,7] <= -1 & elaB2[i,9] <= -1){
       resB2[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB2[i,7] <= -1 & elaB2[i,9] > -1){
       resB2[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB2[i,7] > -1 & elaB2[i,9] <= -1){
       resB2[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB2[i,7] > -1 & elaB2[i,9] > -1){
       resB2[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB2[,3:5], color=resB2, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding avoidance",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding avoidance");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==2,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 3  -------------------------------------------------#
setEPS();
cairo_ps("EC_Beta_3_z1.eps",family="Arial",height=7,width=7);
elaB3 <- ela[ela[,2]==3,]; # Now we have juat the Beta = 0 results in ela
resB3 <- rep(0,dim(elaB3)[1]); # This vector will find one of four results
for(i in 1:dim(elaB3)[1]){ # Doing this as a loop just to show the logic
   if(elaB3[i,7] <= -1 & elaB3[i,9] <= -1){
       resB3[i] <- "red"; # Pre and post-copulatory evolved
   }
   if(elaB3[i,7] <= -1 & elaB3[i,9] > -1){
       resB3[i] <- "blue"; # Pre but not post-copulatory evolved
   }
   if(elaB3[i,7] > -1 & elaB3[i,9] <= -1){
       resB3[i] <- "orange"; # Post but not pre-copulatory evolved
   }
   if(elaB3[i,7] > -1 & elaB3[i,9] > -1){
       resB3[i] <- "grey80"; # Neither pre nor post evolved
   }
}
scatterplot3d(elaB3[,3:5], color=resB3, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding avoidance",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding avoidance");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==3,")")),
       c("Pre and post", "Pre not post", "Post not pre", "Neither"),
       fill=c("red", "blue", "orange", "grey80"));
dev.off();
# --------------------------------------------------------------------------#
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX    Now look plots over time with no costs   XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

rm(list=ls()); 
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

evo_nc <- read.table(file="evo_no_cost.txt", header=FALSE);
# -----------  BETA = 0  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B0.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==0,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==0,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();

# -----------  BETA = 1  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B1.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==1,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==1,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();

# -----------  BETA = 2  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B2.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==2,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==2,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();

# -----------  BETA = 3  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B3.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==3,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==3,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX  Plots over time -- no cost rand initial traits XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

rm(list=ls()); 
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

evo_nc <- read.table(file="evo_no_cost_z1.txt", header=FALSE);
# -----------  BETA = 0  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B0_z1.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==0,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==0,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();

# -----------  BETA = 1  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B1_z1.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==1,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==1,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();

# -----------  BETA = 2  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B2_z1.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==2,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==2,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();

# -----------  BETA = 3  -------------------------------------------------#
setEPS(); # postscript below for final publication?
cairo_ps("evo_no_cost_B3_z1.eps",family="Arial",height=7,width=7);
evo <- evo_nc[evo_nc[,2]==3,];
mWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=mean);
mPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=mean);
mEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=mean);
sWP <- tapply(X=evo[,7],INDEX=evo[,6],FUN=sd);
sPO <- tapply(X=evo[,8],INDEX=evo[,6],FUN=sd);
sEP <- tapply(X=evo[,9],INDEX=evo[,6],FUN=sd);
xxx <- tapply(X=evo[,6],INDEX=evo[,6],FUN=mean);
len <- dim(evo)[1];
par(mar=c(5,5,1,1));
plot(x=xxx,y=mWP,type="l",lwd=2,ylim=c(-9,4),col="blue",
     xlab="Generation",ylab="Allele value",cex.lab=2,cex.axis=1.5);
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
abline(h=0,lwd=0.8,lty="dotted");
uWP <- mWP + sWP / sqrt(sum(evo[,6]==4999));
lWP <- mWP - sWP / sqrt(sum(evo[,6]==4999));
uPO <- mPO + sPO / sqrt(sum(evo[,6]==4999));
lPO <- mPO - sPO / sqrt(sum(evo[,6]==4999));
uEP <- mEP + sEP / sqrt(sum(evo[,6]==4999));
lEP <- mEP - sEP / sqrt(sum(evo[,6]==4999));
polygon(y=c(lPO,rev(uPO)),x=c(1:5000,5000:1),border=NA,col="thistle");
points(x=xxx,y=mPO,type="l",lwd=2,col="red");
polygon(y=c(lEP,rev(uEP)),x=c(1:5000,5000:1),border=NA,col="grey70");
points(x=xxx,y=mEP,type="l",lwd=2,col="black");
polygon(y=c(lWP,rev(uWP)),x=c(1:5000,5000:1),border=NA,col="lightblue");
points(x=xxx,y=mWP,type="l",lwd=2,col="blue");
xleg <- seq(from=0,to=2250,by=1);
yleg <- rep(2,length(xleg));
polygon(y=c(yleg,rev(yleg+2.25)),x=c(xleg,rev(xleg)),border="black",col="grey80", lwd=2);
legend("topleft", inset=0.025, bty="n", cex=0.95, 
       title=expression(paste("Evolving trait (",beta==3,")")),
       c("Pre-copulatory inbreeding", "Tendency for polyandry", "Post-copulatory inbreeding"),
       fill=c("blue", "red", "black"));
dev.off();


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX      Now look at 3D scatters of polyandry   XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

rm(list=ls()); 
library(scatterplot3d);
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

ela <- read.table(file="ela_zallele_0.txt", header=FALSE);
# -----------  BETA = 0  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_0.eps",family="Arial",height=7,width=7);
elaB0 <- ela[ela[,2]==0,]; # Now we have juat the Beta = 0 results in ela
resB0 <- rep(0,dim(elaB0)[1]); # This vector will find one of four results
for(i in 1:dim(elaB0)[1]){ # Doing this as a loop just to show the logic
   if(elaB0[i,8] >= 1){
       resB0[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB0[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB0[,3:5], color=resB0, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==0,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 1  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_1.eps",family="Arial",height=7,width=7);
elaB1 <- ela[ela[,2]==1,]; # Now we have juat the Beta = 0 results in ela
resB1 <- rep(0,dim(elaB1)[1]); # This vector will find one of four results
for(i in 1:dim(elaB1)[1]){ # Doing this as a loop just to show the logic
   if(elaB1[i,8] >= 1){
       resB1[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB1[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB1[,3:5], color=resB1, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==1,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 2  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_2.eps",family="Arial",height=7,width=7);
elaB2 <- ela[ela[,2]==2,]; # Now we have juat the Beta = 0 results in ela
resB2 <- rep(0,dim(elaB2)[1]); # This vector will find one of four results
for(i in 1:dim(elaB2)[1]){ # Doing this as a loop just to show the logic
   if(elaB2[i,8] >= 1){
       resB2[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB2[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB2[,3:5], color=resB2, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==2,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 3  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_3.eps",family="Arial",height=7,width=7);
elaB3 <- ela[ela[,2]==3,]; # Now we have juat the Beta = 0 results in ela
resB3 <- rep(0,dim(elaB3)[1]); # This vector will find one of four results
for(i in 1:dim(elaB2)[1]){ # Doing this as a loop just to show the logic
   if(elaB3[i,8] >= 1){
       resB3[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB3[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB3[,3:5], color=resB3, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==3,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #

rm(list=ls()); 
library(scatterplot3d);
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

ela <- read.table(file="ela_zallele_1.txt", header=FALSE);

# -----------  BETA = 0  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_0_z1.eps",family="Arial",height=7,width=7);
elaB0 <- ela[ela[,2]==0,]; # Now we have juat the Beta = 0 results in ela
resB0 <- rep(0,dim(elaB0)[1]); # This vector will find one of four results
for(i in 1:dim(elaB0)[1]){ # Doing this as a loop just to show the logic
   if(elaB0[i,8] >= 1){
       resB0[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB0[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB0[,3:5], color=resB0, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==0,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 1  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_1_z1.eps",family="Arial",height=7,width=7);
elaB1 <- ela[ela[,2]==1,]; # Now we have juat the Beta = 0 results in ela
resB1 <- rep(0,dim(elaB1)[1]); # This vector will find one of four results
for(i in 1:dim(elaB1)[1]){ # Doing this as a loop just to show the logic
   if(elaB1[i,8] >= 1){
       resB1[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB1[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB1[,3:5], color=resB1, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==1,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 2  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_2_z1.eps",family="Arial",height=7,width=7);
elaB2 <- ela[ela[,2]==2,]; # Now we have juat the Beta = 0 results in ela
resB2 <- rep(0,dim(elaB2)[1]); # This vector will find one of four results
for(i in 1:dim(elaB2)[1]){ # Doing this as a loop just to show the logic
   if(elaB2[i,8] >= 1){
       resB2[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB2[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB2[,3:5], color=resB2, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==2,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# -----------  BETA = 3  -------------------------------------------------#
setEPS();
cairo_ps("PO_Beta_3_z1.eps",family="Arial",height=7,width=7);
elaB3 <- ela[ela[,2]==3,]; # Now we have juat the Beta = 0 results in ela
resB3 <- rep(0,dim(elaB3)[1]); # This vector will find one of four results
for(i in 1:dim(elaB2)[1]){ # Doing this as a loop just to show the logic
   if(elaB3[i,8] >= 1){
       resB3[i] <- "red"; # Pre and post-copulatory evolved
   }else{
       resB3[i] <- "grey80"; # Pre but not post-copulatory evolved
   }
}
scatterplot3d(elaB3[,3:5], color=resB3, pch=20, type="h", lty.hplot=3, lwd=0.9,
              cex.symbols=1.5, cex.axis=1.0, cex.lab=1.25, box=TRUE,
              xlab="Cost of pre-cop inbreeding preference",
              ylab="Cost of tendency for polyandry",
              zlab="Cost of post-cop inbreeding preference");
#box3d();
legend("topleft", inset=0, bty="n", cex=0.95, 
       title=expression(paste("Evolution (",beta==3,")")),
       c("Polyandry", "Monandry"),
       fill=c("red", "grey80"));
dev.off();
# --------------------------------------------------------------------------#

# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX  Large distribution and correlation plots   XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
#http://musicroamer.com/blog/2011/01/16/r-tips-and-tricks-modified-pairs-plot/
panel.cor.scale <- function(x, y, digits=2, prefix="", cex.cor){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r = (cor(x, y,use="pairwise"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex * abs(r))
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(0, 1, 0, 1))
	r = (cor(x, y,use="pairwise"))
	txt <- format(c(r, 0.123456789), digits=digits)[1]
	txt <- paste(prefix, txt, sep="")
	if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
	text(0.5, 0.5, txt, cex = cex )
}

panel.hist <- function(x, ...){
	usr <- par("usr"); on.exit(par(usr))
	par(usr = c(usr[1:2], 0, 1.5) )
	h <- hist(x, plot = FALSE)
	breaks <- h$breaks; nB <- length(breaks)
	y <- h$counts; y <- y/max(y)
	rect(breaks[-nB], 0, breaks[-1], y, col="grey70", ...)
}

pairs.panels <- function (x,y,smooth=TRUE,scale=FALSE){
	if(smooth){
		if (scale) {
			pairs(x,
			diag.panel  = panel.hist,
			upper.panel = panel.cor.scale,
			lower.panel = panel.smooth)
		}
		else{
			pairs(x,
				diag.panel  = panel.hist,
				upper.panel = panel.cor,
				lower.panel = panel.smooth)
		}
	}else{ 
		if(scale){
			pairs(x,
				diag.panel=panel.hist,
				upper.panel=panel.cor.scale)
			}else{pairs(x,
				diag.panel=panel.hist,
				upper.panel=panel.cor)
			}
	} #end of else (smooth)
} #end of function

dat <- cbind(ela[,2:5],ela[,7:9]);
#colnames(dat) <- c("Process","Generation","WP-Pref","Polyandry","EP-Pref","F-coef");
setEPS(); # postscript below for final publication?
#postscript("temp.eps",family="Arial",paper="special",onefile=FALSE);
postscript("tempGEN.eps",family="Arial");
pairs.panels(x=dat);
dev.off();
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #








