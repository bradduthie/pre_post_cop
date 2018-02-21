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

new_gen <- function(ID, Rmat, cost, imm, beta, mu, Kf = 100, Km = 100, nn = 6){
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
    
}





