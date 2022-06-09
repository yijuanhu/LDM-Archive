library(LDM)
library(biomformat)
library(ape) #to read in tree
library(ggplot2)
library(ips) # average.unifrac
library(castor) # average.unifrac
library(phangorn) # average.unifrac


options(echo=TRUE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)

if(length(args)==0){
    print("No arguments supplied")
    
    dist_type <- 4      # 5: "jaccard"; 4: "unwt-unifrac"
    data_full <- 0      # 1: full sample; 0: a subsample of 50 subjects
    i_seed <- 1         # seed for generating different subsamples
    
} else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

if (dist_type==5) {
    dist_method = "jaccard"
} else if (dist_type==4) {
    dist_method = "unwt-unifrac"
} 


###############################################
# read in data
###############################################

# read in otu table, taxonomy info
dat <- read_biom("rawData/otu_table.biom")
otu_table <- as.data.frame(as.matrix(biom_data(dat)))
taxonomy <- observation_metadata(dat)
crohn.tree <- read.tree("rawData/insertion_tree.relabelled.tre")

otu_table_col <- as.data.frame(t(otu_table)) #make samples rows, OTUs columns

# read in metadata
metadata <- read.delim("rawData/1939_20180418-110402.txt", header = TRUE, sep = "\t", dec = ".")

# create dichotomous variable for diagnosis
metadata <- cbind(metadata,ifelse(metadata$diagnosis=="CD" | metadata$diagnosis=="IC" | metadata$diagnosis=="UC", 1, 0))
colnames(metadata)[61] <- "disease" 

# sort metadata to match otu table
meta.sorted <- metadata[order(match(metadata$sample_name, rownames(otu_table_col))),]

# remove duplicates at the end
meta.site.dup <- subset(meta.sorted, collection=="RISK")
meta.sub.dup <- subset(meta.site.dup, type_sample=="biopsy" & biopsy_location=="Rectum")

# want subset with only biopsy samples, no follow-up measurements, no missing values in key covariates
no.follow <- meta.sorted[!duplicated(meta.sorted$anonymized_name),]
meta.site <- subset(no.follow, collection=="RISK")
meta.sub <- subset(meta.site, type_sample=="biopsy" & biopsy_location=="Rectum")

# match otu table with cleaned meta data
otu.sub <- subset(otu_table_col, rownames(otu_table_col) %in% meta.sub$sample_name)

# filter otu table 
pa.table <- 1*(otu.sub > 0)
otu.sub.f <- otu.sub[, -which(colSums(pa.table)<5)]  
dim(otu.sub.f)  

# remove obs with < 5000 reads
otu.full <- otu.sub.f[rowSums(otu.sub.f)>5000,]
dim(otu.full) 
min(rowSums(otu.full))
meta.full <- subset(meta.sub, meta.sub$sample_name %in %rownames(otu.full))
p.t.full = t.test(rowSums(otu.full) ~ meta.full$disease, var.equal=FALSE)$p.value
p.wilcox.full = wilcox.test(rowSums(otu.full) ~ meta.full$disease)$p.value

# --------------------------
# a subsample of 50 subjects
# --------------------------

set.seed(i_seed)
otu.sam60 <- otu.full[sort(sample(nrow(otu.full), 50)), ]  
meta.sam60 <- subset(meta.sub,meta.sub$sample_name %in% rownames(otu.sam60))


########################
# prepare data for anal
########################

if (data_full==1) {
    otu.tab = otu.full
    meta.data = meta.full
} else {
    otu.tab = otu.sam60
    meta.data = meta.sam60
}
dim(otu.tab)
dim(meta.data)


########################
# permanova-F
########################

if (data_full==1) {
    formula = otu.tab | (sex + antibiotics) ~ disease
} else {
    formula = otu.tab | sex ~ disease
}

res.F <- permanovaFL(formula = formula,
                                data=meta.data, tree=crohn.tree, dist.method = dist_method,
                                scale.otu.table=FALSE, n.rarefy=100, n.rej.stop=100, seed=123)
p_permanovaF <- res.F$p.permanov


# -----------------------------
# dist for permanova-D2 (1-100)
# -----------------------------

n_rarefy = 100
label = colnames(otu.full)

n_sam_new = nrow(otu.tab)
dist_sq_all = array(0, dim=c(n_sam_new, n_sam_new, n_rarefy))

set.seed(123)

for (r in 1:n_rarefy) {
    otu_table_rarefy = Rarefy(otu.tab)$otu.tab.rff
    otu_table_rarefy1 = (otu_table_rarefy>0)*1
    colnames(otu_table_rarefy1) = label
    colnames(otu_table_rarefy) = label
    if (dist_type == 5 ) {
        dist <- vegdist(x=otu_table_rarefy1, method=dist_method)
    } else if (dist_type == 4) {
        dist <- GUniFrac(otu_table_rarefy1, crohn.tree, alpha=c(1))$unifrac[,,"d_UW"]
    }
    dist <- as.matrix(dist)
    
    if (r==1) {
        dist_sq_all[,,1] = dist^2
    } else {
        dist_sq_all[,,r] = dist_sq_all[,,r-1] + dist^2
    }
}

# --------------------------
# dist for permanova-D2-A
# --------------------------
    
if (dist_type == 5 ) {
    res <- jaccard.mean( otu.tab)
    dist_sq_inf_o1 = res$jac.mean.sq.o1
    dist_sq_inf_o2 = res$jac.mean.sq.o2
} else if (dist_type == 4) {
    res <- unifrac.mean(otu.tab, crohn.tree)
    dist_sq_inf_o1 = res$unifrac.mean.sq.o1
    dist_sq_inf_o2 = res$unifrac.mean.sq.o2
}


#########################
# permanova-D2, D2-A
#########################

if (data_full==1) {
    formula = dist_sq_r | (sex + antibiotics) ~ disease
} else {
    formula = dist_sq_r | sex ~ disease
}

p_permanovaD2 = NULL
p_permanovaD2A2 = NULL
p_permanovaD2A1 = NULL

for (r in c(1:10, 50, 100, 100^2, 100^3)) { # "100^2" corresponds to D2-A2, "100^3" corresponds to D2-A1
    cat("number of rarefaction:", r, "\n")

    if (r == 100^2) {
        dist_sq_r = dist_sq_inf_o2
    } else if (r == 100^3) {
        dist_sq_r = dist_sq_inf_o1
    } else {
        dist_sq_r = dist_sq_all[,,r]/r
    }
    res_permanovaD2 = permanovaFL(formula = formula, 
                                  tree=tree, data = meta.data, n.perm.max=5000, 
                                  square.dist=FALSE, n.rarefy=0, n.rej.stop=100, seed=123)
    
    if (r == 100^2) {
        p_permanovaD2A2 = res_permanovaD2$p.permanova
    } else if (r == 100^3) {
        p_permanovaD2A1 = res_permanovaD2$p.permanova
    } else {
        p_permanovaD2 = c(p_permanovaD2, res_permanovaD2$p.permanova)
    }
}

tab = c(p_permanovaF, p_permanovaD2, p_permanovaD2A2, p_permanovaD2A1, i_seed)
write.table(t(tab), file = paste("results/datafull", data_full, "_dist", dist_type, ".txt", sep=""), 
            quote=FALSE, append=TRUE, row.names=FALSE, col.names=FALSE)
