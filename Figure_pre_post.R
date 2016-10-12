
rm(list=ls()); 

setwd("~/Dropbox/DuthieProg/c/PolyInbreed/Results");
old_evo <- read.table(file="evo.txt",header=FALSE);
old_evo <- old_evo[old_evo[,6]==4999,];
sd_evo  <- sd(old_evo[,9]);                  
                                                     
setwd("~/Dropbox/DuthieManu/pre_post_cop/B3Results");

nme <- "evo310389751.txt";

evo  <- read.table(nme, header=FALSE);
cpre <- evo[1,3];
cpol <- evo[1,4];
cpos <- evo[1,5];
gen  <- evo[,6];
pre  <- evo[,7];
pol  <- evo[,8];
pos  <- evo[,9];

par(mar=c(5,5,1,1));
plot(x=gen, y=pre, type="l", lwd=2, col="black", ylim=c(-8,10),
     xlab="Generation", ylab="Trait", cex.axis=1.5, cex.lab=1.5);
polygon(x=c(-500,5500,5500,-500), y=c(-1,-1,1,1),border=NA,col="grey90");
box();
points(x=gen, y=pre, type="l", lwd=2, col="black");
points(x=gen, y=pol, type="l", lwd=2, col="red");
points(x=gen, y=pos, type="l", lwd=2, col="blue");
abline(h=0, lty="dotted", lwd=0.8);
legend(x=3400,y=10.3,legend=c("Pre-cop","Polyandry","Post-cop"), cex=1.25,
       col=c("black","red","blue"), lty=c("solid","solid","solid"), lwd=3);


text(x=800, y=10.0, labels=paste("Pre_cost = ",cpre), cex=1.25);
text(x=800, y=9.00, labels=paste("Pol_cost = ",cpol), cex=1.25);
text(x=800, y=8.00, labels=paste("Pos_cost = ",cpos), cex=1.25);




rm(list=ls()); 

setwd("~/Dropbox/DuthieManu/pre_post_cop/presentation");

x  <- seq(from=0, to=1, by=0.0001);
b0 <- exp(-0 * x);
b1 <- exp(-1 * x);
b2 <- exp(-2 * x);
b3 <- exp(-3 * x);
setEPS();
cairo_ps("betas.eps",family="Arial",height=6,width=6);
par(mar=c(5,5,1,1), lwd=3);
plot(x=x, y=b0, cex.lab=1.5, cex.axis=1.5, type="l", lwd=3, ylim=c(0,1),
     ylab=expression(paste("Offspring viability (",Psi[off],")")),
     xlab=expression(paste("Kinship between parents (",k[ij],")")));
points(x=x, y=b1, cex.lab=1.5, cex.axis=1.5, type="l", lwd=3)
points(x=x, y=b2, cex.lab=1.5, cex.axis=1.5, type="l", lwd=3)
points(x=x, y=b3, cex.lab=1.5, cex.axis=1.5, type="l", lwd=3)
dev.off();

