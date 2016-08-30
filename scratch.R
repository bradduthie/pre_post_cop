# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX CORRELATIONS BETWEEN M[a] and F[a]  XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# ------------------------------------------------------------------------#
# -----------  BETA = 0  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 0;
setEPS(); # postscript below for final publication?
cairo_ps("end_cor_WP_PO_B0.eps",family="Arial",height=10,width=8);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="A", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="B", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="C", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="D", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="E", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.00 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="F", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.00 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="G", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="H", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
mtext(expression(paste("End mean inbreeding mating strategy allele value (",M[a],")")),
	outer=TRUE,side=1,line=3.0,cex=1.2);
mtext(expression(paste("End mean tendency for polyandry allele value (",P[a],")")),
	outer=TRUE,side=2,line=2.5,cex=1.2);
dev.off();
# ------------------------------------------------------------------------#

# ------------------------------------------------------------------------#
# -----------  BETA = 1  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 1;
setEPS(); # postscript below for final publication?
cairo_ps("end_cor_WP_PO_B1.eps",family="Arial",height=10,width=8);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="A", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="B", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="C", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="D", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="E", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.00 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="F", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.00 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-10,0,10), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="G", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-20,-10,0,10), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="H", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
mtext(expression(paste("End mean inbreeding mating strategy allele value (",M[a],")")),
	outer=TRUE,side=1,line=3.0,cex=1.2);
mtext(expression(paste("End mean tendency for polyandry allele value (",P[a],")")),
	outer=TRUE,side=2,line=2.5,cex=1.2);
dev.off();

# ------------------------------------------------------------------------#
# -----------  BETA = 2  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 2;
setEPS(); # postscript below for final publication?
cairo_ps("end_cor_WP_PO_B2.eps",family="Arial",height=10,width=8);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="A", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="B", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="C", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="D", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="E", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.00 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="F", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.00 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-10,0,10), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="G", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-20,-10,0,10), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="H", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
mtext(expression(paste("End mean inbreeding mating strategy allele value (",M[a],")")),
	outer=TRUE,side=1,line=3.0,cex=1.2);
mtext(expression(paste("End mean tendency for polyandry allele value (",P[a],")")),
	outer=TRUE,side=2,line=2.5,cex=1.2);
dev.off();

# ------------------------------------------------------------------------#
# -----------  BETA = 3  -------------------------------------------------#
# ------------------------------------------------------------------------#
beta_val <- 3;
setEPS(); # postscript below for final publication?
cairo_ps("end_cor_WP_PO_B3.eps",family="Arial",height=10,width=8);
par(mfrow=c(4,2),oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5));
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="A", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="B", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="C", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0 & evo_cc_02[,4]==0 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="D", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F00.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.00)),cex=1.25, pos=4);
text(x=-25, y=18, label="E", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P00_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.00 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
#axis(side=1,at=c(-10,0,10,20), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="F", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M00_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.00 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-10,0,10), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.00)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="G", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
evo_cc_02 <- read.table(file="res_c02s/evo_M02_P02_F02.txt", header=FALSE);
evB <- evo_cc_02[evo_cc_02[,3]==0.02 & evo_cc_02[,4]==0.02 & evo_cc_02[,5]==0.02,];
evo <- evB[evB[,2]==beta_val & evB[,6]==gens,]; 
WP  <- evo[,7];
PO  <- evo[,8];
EP  <- evo[,9];
plot(x=WP,y=PO,lwd=2,ylim=c(-20,20),xlim=c(-25,10),col="black",xaxt="n",yaxt="n",
     xlab="",ylab="",cex.lab=2,cex.axis=1.5, lty="solid");
mod <- lm(PO~WP);
a   <- mod$coefficients[1];
b   <- mod$coefficients[2];
xxx <- seq(from=min(WP),to=max(WP), by=0.01);
yyy <- a + xxx * b;
points(x=xxx,y=yyy,type="l",lwd=3);
abline(h=0,lwd=0.8,lty="dotted");
#axis(side=2,at=c(-10,0,10,20), cex.axis=1.5);
axis(side=1,at=c(-20,-10,0,10), cex.axis=1.5);
text(x=-25, y=-14, label=expression(paste(c[M]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-17, label=expression(paste(c[P]==0.02)),cex=1.25, pos=4);
text(x=-25, y=-20, label=expression(paste(c[F]==0.02)),cex=1.25, pos=4);
text(x=-25, y=18, label="H", cex=2, pos=4);
rho1 <- as.character(round(cor(WP,EP),digits=3));
text(x=-1, y=18, label=expression(paste("r = ")), cex=1.6, pos=4);
text(x=2, y=18, label=rho1, cex=1.6, pos=4);
rm(evo_cc_02);
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
mtext(expression(paste("End mean inbreeding mating strategy allele value (",M[a],")")),
	outer=TRUE,side=1,line=3.0,cex=1.2);
mtext(expression(paste("End mean tendency for polyandry allele value (",P[a],")")),
	outer=TRUE,side=2,line=2.5,cex=1.2);
dev.off();

