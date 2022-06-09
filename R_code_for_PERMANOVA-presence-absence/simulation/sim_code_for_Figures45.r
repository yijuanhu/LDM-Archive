library(LDM) # permanovaFL
library(dirmult) # simPop
library(MASS)
library(adephylo) # distTips
library(cluster) # pam
library(ips) # average.unifrac
library(castor) # average.unifrac
library(phangorn) # average.unifrac


#-----------------------------
# simulation parameters
#-----------------------------

options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

if(length(args)==0){
    print("No arguments supplied")
    
    data <- 1           # 1: throat data with 856 taxa; # 2:dune data (vegan) with 30 species
    causal_type <- 1    # 1: a random set of 100 or 5 otus; 4: a cluster of otus
    dist_type <- 5      # 5: jaccard; 4: unwt-unifrac

    lib_mu1 <- 10000     # mean of library size for controls
    lib_mu2 <- 5000      # mean of library size for cases
    lib_low <- 2500      # lowest depth
    lib_rff <- 2500      # rarefaction depth
    
    n_rarefy <- 100     # number of rarefaction
    beta <- 0.2         # effect size 
    disp <- 0.02        # overdispersion
    
    n_sim <- 1
    i_seed <- 1
    
    save_all_permanovaD = 1   # save results for different numbers of rarefactions? 
    
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

Y_type <- 1         # binary trait
n_sam <- 100        # total sample size
if (dist_type==5) {
    dist_method = "jaccard"
} else if (dist_type==4) {
    dist_method = "unwt-unifrac"
} 

filename.prefix <- paste("Y",Y_type,"_cau",causal_type,"_d", dist_type, "_data", data, 
                         "_libK-",lib_mu1/1000,"-",lib_mu2/1000,"-",lib_low/1000,"-",lib_rff/1000, 
                         "_nrff", n_rarefy ,"_b", beta, "_d", disp, "_seed", i_seed, sep="")


#-----------------------------
# read in data on otus and trees
# read in estimated pi (frequencies) 
#        and disp (dispersion) 
#-----------------------------
    
tree=NULL
if (dist_method=="unwt-unifrac") {
    data(throat.tree)
    tree=throat.tree
}

if (data == 1) { # throat microbiome data
    
    data(throat.otu.tab)
    label <- colnames(throat.otu.tab)
    otu.tab <- throat.otu.tab
    
    null_dirmult <- NULL
    null_dirmult$pi <- read.table("input_throat/fit_dirmult_pi.txt", header=FALSE, as.is=TRUE)[,1]
    null_dirmult$disp <- read.table("input_throat/fit_dirmult_theta.txt", header=FALSE, as.is=TRUE)[1,1]
    
} else if (data == 2) { # dune ecology data
    
    data(dune)
    data(dune.env)
    dune1 <- dune[which(dune.env$Management=="HF"),]
    dune1 <- dune1[,which(colSums(dune1)>0)]
    label <- colnames(dune1)
    otu.tab <- dune1
    
    ### estimating parameters and save
    # res <- dirmult(dune1)
    # write.table(res$pi, "input_dune/fit_dirmult_pi.txt", row.names=FALSE, col.names=FALSE)
    # write.table(res$theta, "input_dune/fit_dirmult_theta.txt", row.names=FALSE, col.names=FALSE)
    null_dirmult <- NULL
    null_dirmult$pi <- read.table("input_dune/fit_dirmult_pi.txt", header=FALSE, as.is=TRUE)[,1]
    null_dirmult$disp <- read.table("input_dune/fit_dirmult_theta.txt", header=FALSE, as.is=TRUE)[1,1]
}

pi = null_dirmult$pi # baseline relative abundances for simulation
n_otus = length( pi ) 


#-----------------------------
# causal mechanisms
#-----------------------------

if (causal_type==1) { # a random set of 100 (throat data) or 5 (dune data) otus
    
    set.seed(123)
    freq_table <- t( scale( t(otu.tab), center=FALSE, scale=rowSums(otu.tab) ) )
    otu_abundance <- colMeans(freq_table)
    n_common_otus <- ifelse(data==1, 20, 1)
    common_otus <- order(otu_abundance, decreasing=TRUE)[1:n_common_otus]
    uncommon_otus = setdiff(1:n_otus, common_otus)
    top1_otu <- common_otus[1]
    
    n_otus_causal = ifelse(data==1, 100, 5)
    n_otus_confdg = ifelse(data==1, 100, 5)
    causal_otus = sample(uncommon_otus, size=n_otus_causal, replace=FALSE)
    confdg_otus = sample(uncommon_otus, size=n_otus_confdg, replace=FALSE)
    
}

if (causal_type==4) { # tree-based mechanism
    set.seed(123)
    
    freq_table <- t( scale( t(otu.tab), center=FALSE, scale=rowSums(otu.tab) ) )
    otu_abundance <- colMeans(freq_table)
    n_common_otus <- ifelse(data==1, 20, 1)
    common_otus <- order(otu_abundance, decreasing=TRUE)[1:n_common_otus]
    uncommon_otus = setdiff(1:n_otus, common_otus)
    top1_otu <- common_otus[1]
    ### determining the cluster
    # D_patristic=distTips(throat.tree, tips=label, method="patristic", useC=TRUE)
    # cluster20=pam(D_patristic, 20, diss = TRUE, 
    #              cluster.only = TRUE,do.swap = TRUE,keep.diss = FALSE,keep.data = FALSE,pamonce = FALSE, trace.lev = 0)
    # cluster20_names=names(cluster20)
    # A_names=cluster20_names[which(cluster20==14)]
    # causal_otus=match(A_names, label)
    causal_otus = c(308,306,307,300,303,302,301,304,305,286,285,284,280,279,281,282,283,278,224,277,276,298,
                    299,336,335,339,337,338,440,402,401,327,325,326,323,324,333,334,426,427,428,431,430,429,
                    423,425,424,470,471,472,398,399,397,394,395,396,474,421,422,415,419,418,420,416,417,441,
                    439,433,434,435,438,437,332,445,444,446,448,449,447,443,442,467,463,468,469,464,466,465,
                    450,451,436,473,459,461,457,458,462,460,455,454,453,452,456,432,403,404,414,413,321,320,
                    318,317,316,319,311,310,313,314,312,315,309,322,350,349,348,347,352,351,410,407,408,409,
                    411,412,406,405,328,331,330,329,400,377)
    sum(otu_abundance[causal_otus]) # 25.86% abundance
    causal_otus = intersect(causal_otus, uncommon_otus)
    sum(otu_abundance[causal_otus]) # 6% abundance
    
    confdg_otus = sample(uncommon_otus, size=100, replace=FALSE)
    n_otus_causal = length(causal_otus)
    n_otus_confdg = length(confdg_otus)
    
}


##################################
# simulation
##################################

for (sim in c(1:n_sim)) {
  
    print(sim) 
    
    set.seed(i_seed*1000+sim)
        
    #--------------------------------------------------------
    # simulate Y (trait), X (confounder), relative abundances
    #--------------------------------------------------------
    
    if (causal_type %in% c(1,4)) {
        
        if (Y_type == 1) {
            n_sam_group0 = n_sam/2
            n_sam_group1 = n_sam/2
            w0 = 1:n_sam_group0
            w1 = (n_sam_group0+1):n_sam
            
            Y = c(rep(0, n_sam_group0), rep(1, n_sam_group1))
            X = c(rep(0, n_sam_group0*0.3), rep(1, n_sam_group0*0.7), rep(0, n_sam_group1*0.7), rep(1, n_sam_group1*0.3))
            wX = which(X==1)
            len_wX = length(wX)
            wY = which(Y==1)
            len_wY = length(wY)
            
            betaX = 0.5
            betaY = beta
            
            pi_mixture = rep(1,n_sam) %o% pi
            pi_mixture[wX, confdg_otus] = pi_mixture[wX, confdg_otus] * (matrix(runif(len_wX*n_otus_confdg),nrow=len_wX) > betaX)
            pi_mixture[wY, causal_otus] = pi_mixture[wY, causal_otus] * (matrix(runif(len_wY*n_otus_causal),nrow=len_wY) > betaY)
            pi_mixture[,top1_otu] = pi_mixture[,top1_otu] + 1-rowSums(pi_mixture)
            
        }
        
    } 
    
    
    #--------------------------------------------------------
    # simulate otu data
    #--------------------------------------------------------
    
    lib_raw = rep(NA, n_sam)
    if (lib_mu2==lib_low) {
        lib_raw[w0] = lib_mu1
        lib_raw[w1] = lib_mu2
    } else {
        lib_raw[w0] = round(rnorm(length(w0), mean=lib_mu1, sd=lib_mu1/3 ))
        lib_raw[w1] = round(rnorm(length(w1), mean=lib_mu2, sd=lib_mu2/3 ))
    }
    
    lib <- ifelse(lib_raw < lib_low, lib_low, lib_raw)
    
    otu_table_sim <- NULL
    if (disp > 0) {
        for (i in 1:n_sam) {
            otu_table_sim_i <- simPop( J=1, K=n_otus, n=lib[i], pi=pi_mixture[i,], theta=disp)$data # theta is rho^2 in wiki
            otu_table_sim <- rbind(otu_table_sim, otu_table_sim_i)
        }
    } else {
        for (i in 1:n_sam) {
            otu_table_sim_i <- t(rmultinom( n=1, size=lib[i], prob=pi_mixture[i,]))
            otu_table_sim <- rbind(otu_table_sim, otu_table_sim_i)
        }
    }
    colnames(otu_table_sim) = label
    
    # remove samples with library sizes lower than rarefaction depth
    
    w_rm = which(lib < lib_rff)
    if (n_rarefy > 0 & length(w_rm) > 0) {
        otu_table_sim = otu_table_sim[-w_rm,]
        X = X[-w_rm]
        Y = Y[-w_rm]
    }
    
    
    ###############################################
    # permanova-F
    ###############################################
    
    res_permanovaF <- permanovaFL(formula = otu_table_sim | X ~ Y, 
                       tree=tree, n.perm.max=500, scale.otu.table=FALSE, 
                       dist.method=dist_method, binary=TRUE, n.rarefy=n_rarefy)
    
    p_permanovaF <- res_permanovaF$p.permanova

    
    ###############################################
    # permanova-D (including permanova-D-A)
    # permanova-D2 (including permanova-D2-A)
    # adonis2
    ###############################################
    
    set.seed(res_permanovaF$seed)
    
    n_sam_new = nrow(otu_table_sim)
    if (save_all_permanovaD) {
        dist_all = array(0, dim=c(n_sam_new, n_sam_new, n_rarefy))
        dist_sq_all = array(0, dim=c(n_sam_new, n_sam_new, n_rarefy))
    } else {
        dist_all = matrix(0, n_sam_new, n_sam_new)
        dist_sq_all = matrix(0, n_sam_new, n_sam_new)
    }

    for (r in 1:n_rarefy) {
        otu_table_rarefy = Rarefy(otu_table_sim)$otu.tab.rff
        otu_table_rarefy1 = (otu_table_rarefy>0)*1
        colnames(otu_table_rarefy1) = label
        colnames(otu_table_rarefy) = label
        if (dist_type == 5 ) {
            dist <- vegdist(x=otu_table_rarefy1, method=dist_method)
        } else if (dist_type == 4) {
            dist <- GUniFrac(otu_table_rarefy1, tree, alpha=c(1))$unifrac[,,"d_UW"]
        }
        dist <- as.matrix(dist)
        
        if (save_all_permanovaD) {
            if (r==1) {
                dist_all[,,1] = dist
                dist_sq_all[,,1] = dist^2
            } else {
                dist_all[,,r] = dist_all[,,r-1] + dist
                dist_sq_all[,,r] = dist_sq_all[,,r-1] + dist^2
            }
        } else {
            dist_all = dist_all + dist
            dist_sq_all = dist_sq_all + dist^2
        }
    }
    
    if (dist_type == 5 ) {
        res <- jaccard.mean( otu_table_sim )
        dist_inf = res$jac.mean.o2
        dist_sq_inf = res$jac.mean.sq.o2
    } else if (dist_type == 4) {
        res <- unifrac.mean(otu_table_sim, tree )
        dist_inf = res$unifrac.mean.o2
        dist_sq_inf = res$unifrac.mean.sq.o2
    }

    p_permanovaD = NULL
    p_permanovaD2 = NULL
    p_adonis2D = NULL
    p_adonis2D2 = NULL
    if (save_all_permanovaD) {
        for (r in c(1:10, 50, 100, 100^2)) { # 100^2 corresponds to our proposed approach of using infinite number of rarefactions
            if (r == 100^2) {
                dist_r = dist_inf
            } else {
                dist_r = dist_all[,,r]/r
            }
            res_permanovaD = permanovaFL(formula = dist_r | X ~ Y, 
                                          tree=tree, dist.method=dist_method, n.perm.max=500, 
                                          n.rarefy=0, seed=res_permanovaF$seed)
            
            set.seed(res_permanovaF$seed)
            res_adonis2D = adonis2(formula = dist_r ~ X+Y)
            
            # recommended
            if (r == 100^2) {
                dist_sq_r = dist_sq_inf
            } else {
                dist_sq_r = dist_sq_all[,,r]/r
            }
            res_permanovaD2 = permanovaFL(formula = dist_sq_r | X ~ Y, 
                                          tree=tree, dist.method=dist_method, n.perm.max=500, 
                                          square.dist=FALSE, n.rarefy=0, seed=res_permanovaF$seed)
            
            dist_sq_r_sqrt = sqrt(dist_sq_r)
            set.seed(res_permanovaF$seed)
            res_adonis2D2 = adonis2(formula = dist_sq_r_sqrt ~ X+Y)
            
            p_permanovaD = c(p_permanovaD, res_permanovaD$p.permanova)
            p_permanovaD2 = c(p_permanovaD2, res_permanovaD2$p.permanova)
            p_adonis2D = c(p_adonis2D, res_adonis2D$`Pr(>F)`[2])
            p_adonis2D2 = c(p_adonis2D2, res_adonis2D2$`Pr(>F)`[2])
        }
    } else {
        res_permanovaD <- permanovaFL(formula = dist_all | X ~ Y, 
                                      tree=tree, dist.method=dist_method, n.perm.max=500, 
                                      n.rarefy=0, seed=res_permanovaF$seed)
        
        res_adonis2D = adonis2(formula = dist_all ~ X+Y)
        
        # recommended
        res_permanovaD2 <- permanovaFL(formula = dist_sq_all | X ~ Y, 
                                       tree=tree, dist.method=dist_method, n.perm.max=500, 
                                       square.dist=FALSE, n.rarefy=0, seed=res_permanovaF$seed)
        
        dist_sq_all_sqrt = sqrt(dist_sq_all)
        res_adonis2D2 = adonis2(formula = dist_sq_all_sqrt ~ X+Y)
        
        p_permanovaD = res_permanovaD$p.permanova
        p_permanovaD2 = res_permanovaD2$p.permanova
        p_adonis2D = res_adonis2D$`Pr(>F)`[2]
        p_adonis2D2 = res_adonis2D2$`Pr(>F)`[2]
    }
    
    
    ###########################
    # permanova-UR, permanova-L
    ###########################
    
    otu_table_sim1 = (otu_table_sim>0)*1
    colnames(otu_table_sim1) = label
    if (dist_type == 5 ) {
        dist_ur <- vegdist(x=otu_table_sim1, method=dist_method)
    } else if (dist_type == 4) {
        dist_ur <- GUniFrac(otu_table_sim1, tree, alpha=c(1))$unifrac[,,"d_UW"]
    }
    dist_ur <- as.matrix(dist_ur)
    
    # permanova-UR
    res_permanovaUR <- permanovaFL(formula = dist_ur | X ~ Y, 
                                   tree=tree, dist.method=dist_method, n.perm.max=500, 
                                   n.rarefy=0, seed=res_permanovaF$seed)
    p_permanovaUR = res_permanovaUR$p.permanova
    
    # permanova-L
    L <- rowSums(otu_table_sim)
    res_permanovaL <- permanovaFL(formula = dist_ur | X + L ~ Y, 
                                  tree=tree, dist.method=dist_method, n.perm.max=500, 
                                  n.rarefy=0, seed=res_permanovaF$seed)
    p_permanovaL = res_permanovaL$p.permanova
    
    
    #########################
    # Output
    #########################
    
    tab <- c(p_permanovaF, p_permanovaD, p_adonis2D, p_permanovaD2, p_adonis2D2, p_permanovaUR, p_permanovaL)
    
    write.table(t(round(tab, 3)), paste(filename.prefix, ".txt", sep=""), append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE)
    
} # simulation

