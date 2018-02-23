initialise_inds <- function(N = 200){
    ID      <- matrix(data = 0, nrow = N, ncol = 30);
    ID[,1]  <- 1:N;
    ID[,2]  <- rbinom(n = N, size = 1, prob = 0.5);
    ID[,3]  <- sample(x = 1:10, size = N, replace = TRUE);
    ID[,4]  <- sample(x = 1:10, size = N, replace = TRUE);
    ID[,5]  <- rep(x = 0, times = N);
    ID[,6]  <- rep(x = -1, times = N);
    ID[,7]  <- rep(x = -1, times = N);
    ID[,8]  <- rep(x = 0, times = N);
    ID[,9]  <- rep(x = -1, times = N);
    ID[,10] <- rep(x = 0, times = N);
    return(ID);
}

initialise_Rmat <- function(ID){
    rows <- dim(ID)[1];
    Rmat <- matrix(data = 0, nrow = rows, ncol = rows);
    diag(Rmat) <- 1;
    Rmat <- cbind(ID[,1], Rmat);
    return(Rmat);
}

Rbuild <- function(Rmat, ID, Rmat_new){
    inds <- dim(ID)[1];
    for(i in 1:inds){
        if(ID[i, 6] %in% Rmat[,1] == FALSE){
            MomP <- -1;
        }else{
            MomP <- which(Rmat[,1] == ID[i, 6]);
        }
        if(ID[i, 7] %in% Rmat[,1] == FALSE){
            DadP <- -1;
        }else{
            DadP <- which(Rmat[,1] == ID[i, 7]);
        }
        # Inbreeding coefficient
        if(MomP < 0 | DadP < 0){
            Rmat_new[i, i+1] <- 1; 
        }else{
            Rmat_new[i, i+1] <- 1 + Rmat[MomP, DadP+1]; 
        }
        # Kinship coefficient
        j <- 1;
        while(j < i){
            if(ID[j, 6] %in% Rmat[,1] == FALSE){
                Om <- -1;
            }else{
                Om <- which(Rmat[,1] == ID[j, 6]);
            }
            if(ID[j, 7] %in% Rmat[,1] == FALSE){
                Od <- -1;
            }else{
                Od <- which(Rmat[,1] == ID[j, 7]);
            }
            
            if(ID[i,5] == 0 & ID[j,5] == 0){
                if(Om < 0 | MomP < 0){
                    w <- 0;
                }else{
                    if(MomP == Om){
                        w <- 0.5 * Rmat[MomP, Om + 1];
                    }else{
                        w <- Rmat[MomP, Om + 1];
                    }
                }
                if(Od < 0 | MomP < 0){
                    y <- 0;
                }else{
                    y <- Rmat[MomP, Od + 1];
                }
                if(Om < 0 | DadP < 0){
                    x <- 0;
                }else{
                    x <- Rmat[DadP, Om + 1];
                }
                if(Od < 0 | DadP < 0){
                    z <- 0;
                }else{
                    if(DadP == Od){
                        z <- 0.5 * Rmat[DadP, Od + 1];
                    }else{
                        z <- Rmat[DadP, Od + 1];
                    }
                }
            }
            rMmuP <- 0.5 * (w + y);
            rDmuP <- 0.5 * (x + z);
            rval  <- 0.5 * (rMmuP + rDmuP);
            
            Rmat_new[i, j+1] <- rval;
            Rmat_new[j, i+1] <- rval;
            j <- j + 1;
        }
    }
    
    for(i in 1:inds){
        for(j in 1:inds){
            if(ID[i, 6] == ID[j, 1]){
                if(ID[j, 6] == -1){
                    Rmat_new[i, j+1] = 0.25;
                    Rmat_new[j, i+1] = 0.25;
                }    
            }
            if(ID[i, 7] == ID[j, 1]){
                if(ID[j, 7] == -1){
                    Rmat_new[i, j+1] = 0.25;
                    Rmat_new[j, i+1] = 0.25;
                }    
            }
        }
    }
    
    return(Rmat_new);
}

mateselect <- function(ID, Rmat_new){
    fem <- ID[ID[,2] == 0 & ID[,5] == 0,];
    mal <- ID[ID[,2] == 1 & ID[,5] == 0,];
    for(i in 1:dim(fem)[1]){
        strat <- sum(fem[i,11:20]);
        vals  <- rep(0, dim(mal)[1]);
        fpos  <- which(Rmat_new[,1] == i); 
        for(j in 1:length(vals)){
            mpos    <- which(Rmat_new[,1] == j);
            kval    <- Rmat_new[i, j+1];
            vals[j] <- kval * strat;
        }
        if(strat < 0){
            weight <- 1 / (1 + (-1 * vals));
            probs  <- weight / sum(weight);
        }else{
            weight <- 1 + vals;
            probs  <- weight / sum(weight);
        }
        mate    <- sample(x = 1:dim(mal)[1], size = 1, pr = probs);
        mateID  <- mal[mate,1];
        matepos <- which(ID[,1] == mateID);
        ownID   <- fem[i, 1];
        ownpos  <- which(ID[,1] == ownID);
        ID[ownpos, 9]  <- mateID;
        ID[matepos, 9] <- ownID;
    }
    return(ID);
}

make_offspring <- function(ID, Rmat_new, nn = 6){
    start_ID <- max(ID[,1]) + 1;
    for(i in 1:dim(ID)[1]){
        if(ID[i, 2] == 0 & ID[i, 9] > 0 & ID[i, 5] == 0){
            mum    <- ID[i,];
            dadpos <- which(ID[,1] == ID[i, 9]);
            dad    <- ID[dadpos,];
            offs   <- matrix(data = 0, nrow = nn, ncol = dim(ID)[2]);
            offs[,1] <- start_ID:(start_ID+nn-1);
            start_ID <- start_ID + nn;
            offs[,2] <- rbinom(n = nn, size = 1, prob = 0.5);
            offs[,3] <- mum[3];
            offs[,4] <- mum[4];
            offs[,6] <- mum[1];
            offs[,7] <- dad[1];
            mrpos    <- which(Rmat_new[,1] == mum[1]);
            drpos    <- which(Rmat_new[,1] == dad[1]);
            offs[,8] <- Rmat_new[mrpos, drpos+1];
            offs[,11:30] <- mum[11:30];
            dum_vec  <- rbinom(n = length(offs[,11:30]), size = 1, prob = 0.5);
            dummy    <- matrix(data = dum_vec, nrow = dim(offs)[1]);
            dummy[dummy > 0] <- NA;
            
        }
    }
}

produce_offs <- function(ID, Rmat_new, nn){
    # NOTE: This toy model has unrealistic recombination between parents!
    mums  <- sum(ID[,2] == 0 & ID[, 5] == 0 & ID[, 9] != 0);
    offs  <- matrix(data = 0, nrow = mums * nn, ncol = dim(ID)[2]);
    start <- 1;
    for(i in 1:dim(ID)[1]){  # Unrealistic recombination saves some time.
        if(ID[i,2] == 0 & ID[i, 5] == 0 & ID[i, 9] != 0){
            her_off     <- matrix(data = ID[i,], nrow = nn, ncol = dim(ID)[2], 
                                  byrow = TRUE);
            her_off[,2] <- rbinom(n = nn, size = 1, prob = 0.5);
            her_off[,6] <- ID[i,1];
            her_off[,7] <- ID[i,9];
            mum_p       <- which(Rmat_new[,1] == ID[i,1]);
            dad_p       <- which(Rmat_new[,1] == ID[i,9]);
            her_off[,8] <- Rmat_new[mum_p, dad_p + 1];
            dad_o       <- which(ID[,1] == ID[i,9]);
            dad_overlay <- matrix(data = ID[dad_o,], nrow = nn, 
                                  ncol = dim(ID)[2], byrow = TRUE);
            to_vec      <- rbinom(n = nn*dim(ID)[2], size = 1, prob = 0.5);
            to_overlay  <- matrix(data = to_vec, nrow = nn);
            to_overlay[,1:10] <- 1;
            her_off[to_overlay < 1] <- NA;
            her_off[is.na(her_off)] <- dad_overlay[is.na(her_off)];
            offs[start:(start+(nn-1)), 1:dim(ID)[2]] <- her_off;
            start <- start + nn;
        }
    }
    max_ID   <- max(ID[,1]);
    offs[,1] <- (max_ID + 1):(max_ID + dim(offs)[1]); 
    return(offs);
}

mortality <- function(offs, beta){
    depression        <- exp(-1 * beta * offs[,8]);
    alive             <- rbinom(n = dim(offs)[1], size = 1, prob = depression);
    offs[alive < 1,5] <- -1;
    offs              <- offs[offs[,5] > -1,];
    return(offs);
}

mutation  <- function(offs, mu){
    mutations <- rbinom(n = 20*dim(offs)[1], size = 1, prob = mu);
    tot_mut   <- sum(mutations);
    mut_mat   <- matrix(data = mutations, nrow = dim(offs)[1]);
    mu_eff    <- rnorm(n = tot_mut, mean = 0, sd = 2);
    mut_mat[mut_mat > 0] <- mu_eff;
    offs[,11:30] <- offs[,11:30] + mut_mat;
    return(offs);
}

retain <- function(ID, offs){
    retainers <- rep(x = 0, times = dim(ID)[1]);
    for(i in 1:dim(ID)[1]){
        ID_num <- ID[i, 1];
        ismum  <- sum(offs[,6] == ID_num);
        isdad  <- sum(offs[,7] == ID_num);
        ispar  <- ismum + isdad;
        if(ispar < 1){
            retainers[i] <- 1;
        }
    }
    to_remove <- which(retainers == 1);
    ID <- ID[-to_remove,];
    return(ID);
}

get_stats <- function(ID){
    mean_inbr <- mean(ID[,8]);
    mean_aval <- mean(ID[,11:20]);
    sd_aval   <- sd(ID[,11:20]);
    mean_nval <- mean(ID[,21:30]);
    sd_nval   <- sd(ID[,21:30]);
    stat_vec  <- c(mean_inbr, mean_aval, sd_aval, mean_nval, sd_nval);
    return(stat_vec);
}

immigration <- function(ID, imm, IDst){
    max_ID     <- max(ID[,1]);
    immigrants <- initialise_inds(N = imm);
    avalues    <- rnorm(n = 10 * imm, mean = IDst[2], sd = IDst[3]);
    nvalues    <- rnorm(n = 10 * imm, mean = IDst[4], sd = IDst[5]);
    a_mat      <- matrix(data = avalues, nrow = imm);
    n_mat      <- matrix(data = nvalues, nrow = imm);
    immigrants[,2]     <- 1;
    immigrants[,11:20] <- a_mat;
    immigrants[,21:30] <- n_mat;
    imm_IDs            <- (max_ID + 1):(max_ID + imm);
    immigrants[,1]     <- imm_IDs;
    ID <- rbind(ID, immigrants);
    return(ID);
}

new_gen <- function(ID, Rmat, cost = 0, imm = 5, beta = 1, 
                    mu = 0.01, Kf = 100, Km = 100, nn = 6){
    N       <- dim(ID)[1];
    death   <- rep(x = 0, times = N);
    females <- sum(ID[,2] == 0 & ID[,5] <= 1);
    males   <- sum(ID[,2] == 1 & ID[,5] <= 1);
    strat   <- apply(X = ID[,11:20], MARGIN = 1, FUN = sum);
    strat_c <- abs(strat) * cost;
    strat_d <- rbinom(n = N, size = 1, prob = strat_c);
    death   <- death + strat_d;
    death   <- death + (ID[,5] > 1);
    ID[,5]  <- -1 * death;
    ID      <- ID[order(ID[,2], decreasing = FALSE),];
    ID      <- ID[order(ID[,5], decreasing = FALSE),];
    livef   <- sum(ID[,2] == 0 & ID[,5] < 2 & ID[,5] > -1);
    livem   <- sum(ID[,2] == 1 & ID[,5] < 2 & ID[,5] > -1);
    if(livef > Kf){
        rm   <- livef - Kf;
        rows <- which(ID[,2] == 0 & ID[,5] < 2 & ID[,5] > -1);
        ID[rows[1:rm],5] <- -1;
    }
    if(livem > Km){
        rm   <- livem - Km;
        rows <- which(ID[,2] == 1 & ID[,5] < 2 & ID[,5] > -1);
        ID[rows[1:rm],5] <- -1;
    }
    ID <- rbind(ID[ID[,5] >= 0,], ID[ID[,5] < 0,]);
    
    Rmat_new <- initialise_Rmat(ID);
    Rmat_new <- Rbuild(Rmat, ID, Rmat_new);
    
    ID   <- mateselect(ID, Rmat_new);
    offs <- produce_offs(ID, Rmat_new, nn);
    offs <- mortality(offs, beta);
    offs <- mutation(offs, mu);
    IDst <- get_stats(ID);
    
    ID[ID[,5] >= 0, 5] <- ID[ID[,5] >= 0,5] + 1;
    ID[ID[,5] > 0, 5]  <- -1;
    
    ID <- retain(ID, offs);
    ID <- rbind(ID, offs);
    
    ID <- immigration(ID, imm, IDst);
    
    results <- list(ID = ID, Rmat = Rmat_new, stats = IDst);
}









females  <- sum(ID[,2] == 0 & ID[,5] == 0);
males    <- sum(ID[,2] == 1 & ID[,5] == 0);
distrb   <- seq(from = 0, to = 2*pi, length.out = females + males + 1);
pchs     <- c(rep(15, females), rep(16, males));
cols     <- c(rep("blue", females), rep("red", males));


plot(x = 0, y = 0, xlim = c(-2, 2), ylim = c(-2, 2), xaxt = "n",
     yaxt = "n", xlab = "", ylab = "", type = "n");

points(x = 2 * cos(distrb), y = 2 * sin(distrb), pch = pchs, col = cols, 
       cex = 3);



mate_connect <- function(ID){
    living     <- ID[ID[,5] == 0,];
    fem_count  <- sum(living[,2] == 0);
    livsor     <- living[order(living[,2]),];
    livsor[,3] <- 1:dim(livsor)[1];
    pointpairs <- matrix(data = 0, nrow = fem_count, ncol = 4);
    for(i in 1:dim(livsor)[1]){
        if(livsor[i, 2] == 0){
            femID <- livsor[i, 1];
            malID <- livsor[i, 9];
            femps <- i;
            malps <- which(livsor[,1] == malID);
            pointpairs[i, 1] <- femID;
            pointpairs[i, 2] <- malID;
            pointpairs[i, 3] <- femps;
        }
    }
    uniquem <- unique(pointpairs[,2]);
    for(i in 1:dim(pointpairs)[1]){
        pointpairs[i, 4] <- which(uniquem == pointpairs[i, 2]) + fem_count;
    }
    return(pointpairs);
}























iallofem <- c(rep(1, females), rep(2, males));
lallofem <- 2 * pi * (1:length(iallofem)) / length(iallofem);
xallofem <- cos(allofem);
yallofem <- sin(allofem);

points(x = xallofem, y = yallofem);






sc3 <- seq(from = 0, to = 2*pi, by = 0.001);
points(x = 2 * cos(sc3), y = 2 * sin(sc3), lwd = 1.5, type = "l");



ID   <- initialise_inds(N = 10);
Rmat <- initialise_Rmat(ID);



gen  <- new_gen(ID = ID, Rmat = Rmat, Kf = 10, Km = 10, beta = 0, imm = 1);
ID   <- gen$ID;
Rmat <- gen$Rmat;

tmat <- Rmat[,2:dim(Rmat)[1]]
diag(tmat) <- diag(tmat)/2;


par(mfrow = c(1,1));



plot(x = x0:x1, y = x0:x1, type = "n", ylim = c(0, 1), cex.axis = 1.5,
     cex.lab = 1.5, xlab = "Generation", ylab = "k");
points(x = ts$df[x0:x1,1], y = ts$df[x0:x1,2], pch = 20, cex = 1.5, type = "b",
       lwd = 2);

par(mar = c(4, 1, 1, 1));
hist(tmat, xlim = c(0, 1), ylim = c(0, 40), freq = FALSE, col = "black",
     yaxt = "n", ylab = "", main = "", xlab = "", 
     breaks = seq(from = 0, to = 1, by = 0.01));
title(xlab = "Kinship (k)", line = 2);






hist(ID[,11:20]);




which(Rmat[,2:dim(Rmat)[1]] > 1)












