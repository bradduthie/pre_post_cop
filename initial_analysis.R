#==========================================================================

rm(list=ls()); 
                                                                     
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

evo  <- read.table(file="evo.txt", header=FALSE);
ela  <- evo[evo[,6]==4999,];
write.table(x=ela, file="ela.txt", col.names=FALSE, row.names=FALSE);

#==========================================================================


rm(list=ls()); 
library(scatterplot3d);
setwd("~/Dropbox/DuthieManu/pre_post_cop/Initial_Results");

ela <- read.table(file="ela.txt", header=FALSE);

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


# XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX XXX #
# XXX XXX XXX  Large distribution and correlation plots  XXX XXX XXX XXX  #
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








