rm(list=ls())
source("simul_functions.R")
source("Library_CorMulti.R")

############################################################################
###################  codes for simulation, single node   ###################
############################################################################


#### Simulation comparison table for single level testing of DM ####

## under the null
rate1 <- rate2 <- c(0.15, 0.05, 0.22, 0.3, 0.03, 0.1, 0.07, 0.08)
theta1 <- 3
theta2 <- 5
ell <- c(0.12, 0.06, 0.08, 0.43, 0.02, 0.14, 0.1, 0.05)
iter <- 5000
rho <- c(0:6)/10
n_sample <- c(20, 50, 100)
res_table_null <- simu_table_DM(iter, n_sample, rate1, rate2, theta1, theta2, ell, rho, seed=0)

## under the alternative
rate1 <- c(0.15, 0.05, 0.22, 0.3, 0.03, 0.1, 0.07, 0.08)
rate2 <- c(0.1, 0.1, 0.22, 0.3, 0.03, 0.1, 0.07, 0.08)
theta1 <- 3
theta2 <- 5
ell <- c(0.12, 0.06, 0.08, 0.43, 0.02, 0.14, 0.1, 0.05)
iter <- 5000
rho <- c(0:6)/10
n_sample <- c(20, 50, 100)
res_table_alt <- simu_table_DM(iter, n_sample, rate1, rate2, theta1, theta2, ell, rho, seed=0)
save(res_table_null, res_table_alt, file="result/res_table_DM.Rdata")

#### Simulation comparison table for single level testing of NormM ####

## under the null
mu1 <- mu2 <- c(3,1,0.5,1,0,1,1,0)
sigma <- rep(1,8)
iter <- 5000
rho <- c(0:6)/10
n_sample <- c(20, 50, 100)
res_table_null <- simu_table_NormM(iter, n_sample, mu1, mu2, sigma, rho, seed=0)

## under the alternative
mu1 <- c(3,1,0.5,1,0,1,1,0)
mu2 <- c(3,1,1,0.5,0,1,1,0)
sigma <- rep(1,8)
iter <- 5000
rho <- c(0:6)/10
n_sample <- c(20, 50, 100)
res_table_alt <- simu_table_NormM(iter, n_sample, mu1, mu2, sigma, rho, seed=0)
save(res_table_null, res_table_alt, file="result/res_table_NormM.Rdata")

#### size and power plots for single level testing ####

# for Figure 2
ind <- c('_DM', '_NormM')
for (ii in 1:2){
    load(paste('result/res_table', ind[ii], '.Rdata', sep=''))
    combo <- res_table_alt$combo
    size_pval <- rowMeans(res_table_null$pval<0.05)
    size_pval.p <- rowMeans(res_table_null$pval.p<0.05)
    power_pval <- rowMeans(res_table_alt$pval<0.05)
    power_pval.p <- rowMeans(res_table_alt$pval.p<0.05)
    
    pdf(paste('single', ind[ii], '.pdf', sep=''), width=7, height=2.5)
    par(mfrow=c(1,3), mar=c(4,4,2,0.5))
    plot(combo[1,1:7],size_pval[1:7], type="b", lty=2, ylim=c(0,1), 
         xlab="", ylab="size/power")
    lines(combo[1,1:7],size_pval.p[1:7], type="b", col="red")
    lines(combo[1,1:7],power_pval[1:7], type="b", lty=2, pch=8)
    lines(combo[1,1:7],power_pval.p[1:7], type="b", pch=8, col="red")
    abline(h=0.05, lty=3, col="blue")
    title("sample size=20")
    legend(0, 1, col=c(2,1,2,1), lty=c(1,2,1,2), pch=c(8,8,1,1),
        c("Paired Power","Unpaired Power","Paired Size","Unpaired Size"))
    plot(combo[1,8:14],size_pval[8:14], type="b", lty=2, ylim=c(0,1), 
         xlab="Correlation Parameter", ylab="")
    lines(combo[1,8:14],size_pval.p[8:14], type="b", col="red")
    lines(combo[1,8:14],power_pval[8:14], type="b", lty=2, pch=8)
    lines(combo[1,8:14],power_pval.p[8:14], type="b", pch=8, col="red")
    abline(h=0.05, lty=3, col="blue")
    title("sample size=50")
    plot(combo[1,15:21],size_pval[15:21], type="b", lty=2, ylim=c(0,1),
         xlab="", ylab="")
    lines(combo[1,15:21],size_pval.p[15:21], type="b", col="red")
    lines(combo[1,15:21],power_pval[15:21], type="b", lty=2, pch=8)
    lines(combo[1,15:21],power_pval.p[15:21], type="b", pch=8, col="red")
    abline(h=0.05, lty=3, col="blue")
    title("sample size=100")
    dev.off()
}



############################################################################
#################  codes for simulation based on real data #################
############################################################################


#### load the data on which the simulation is based ####

extra_gt <- read.csv("dataset/mapsheet.csv",  row.names=1, check.names=F)
cts3 <- read.csv("dataset/tree_count.csv", row.names=1, check.names=F)
tree2 <- read.csv("dataset/tree_struc.csv")

# get rid of some bad samples
gone <- which(cts3[,'g__Lactobacillus']/cts3[,1]>0.25)
extra_gt <- extra_gt[-gone,]
cts3 <- cts3[-gone,]

# find the overlaping subjects
{
	weeks <- c(0,1,2,3,5,8)
	overlap <- extra_gt$PersonalID[extra_gt$Week==0&extra_gt$BodySite=="gut"]
	length(overlap)
	for (i in weeks){
		aim <- extra_gt$PersonalID[extra_gt$Week==i&extra_gt$BodySite=="gut"]
		temp <- match(aim, overlap)
		temp <- temp[!is.na(temp)]
		overlap <- overlap[temp]
		print(length(overlap))
	}

	temp <- match(extra_gt$PersonalID, overlap)
	extra_gt <- extra_gt[!is.na(temp),]
	temp <- match(extra_gt$Week, weeks)
	extra_gt <- extra_gt[!is.na(temp),]
	extra_gt <- extra_gt[order(extra_gt$Personal),]
	rownames(extra_gt) <- extra_gt$Sample
}

# get the proportion database
prp_base <- total2sub(cts3,tree2)
prp_base <- prp_base/cts3[,1]
ttl_ct <- cts3[,1]


#### simulate data from week 0 vs week 5 gut data, sparse setting ####

# get the indices of different nodes
ind <- c('g__Streptococcus')
find_all_parent(ind, tree2)
pct_change <- c(0:10)/500
n_sample <- c(20,50,100)
simu_table_from_real(100, n_sample, tree2, ttl_ct, prp_base, ind,
		pct_change, seed=0, filename="result/res_simu_from_real.Rdata")


#### simulate data from week 0 vs week 5 gut data, dense setting ####

# get the indices of nodes with differential subcompositions
ind1 <- c('g__Streptococcus', 'g__Eubacterium', 'g__Parabacteroides')
ind2 <- c('g__Porphyromonas', 'g__Moraxella', 'g__Ruminococcus')
find_all_parent(ind1[1], tree2)
find_all_parent(ind1[2], tree2)
find_all_parent(ind1[3], tree2)
find_all_parent(ind2[1], tree2)
find_all_parent(ind2[2], tree2)
find_all_parent(ind2[3], tree2)
pct_change <- c(0:10)/500
n_sample <- c(20,50,100)
simu_table_from_real2(100, n_sample, tree2, ttl_ct, prp_base, ind1, ind2,
		pct_change, seed=0, filename="result/res_simu_from_real2.Rdata")

        
#### plot the results for the paper ####

## sparse setting

load("result/res_simu_from_real.Rdata")
combo <- matrix(0, 2, length(n_sample)*length(pct_change))
combo[1,] <- rep(pct_change, length(n_sample))
combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(pct_change))))

loc <- res_table$location
loc.01 <- res_table$location.01
node_freq <- loc[,find_all_parent(ind,tree2)]
node_freq <- node_freq/100

rej <- apply(res_table$pval<=0.05,1,mean)
rej.p <- apply(res_table$pval.p<=0.05,1,mean)
rej.Fisher <- apply(res_table$pval.Fisher<=0.05,1,mean)
rej.uni.p <- apply(res_table$pval.uni.p<=0.05,1,mean)

fdr_mean <- function(x){ x[is.na(x)]<-0; return(apply(x,1,mean)) }
fdr <- fdr_mean(res_table$fdr)


# for Figure 4(a)
pdf("nodes_real.pdf", width=7, height=2.5)
par(mfrow=c(1,3), mar=c(4,4,2,0.5))
plot(c(0, max(pct_change)), c(0,1), type="n", xlab="", ylab="Percent of Discovery/FDR")
for (i in c(1,3,5)){
	lines(combo[1,1:11], node_freq[1:11,i], type="l", lty=i)
}
lines(combo[1,1:11],fdr[1:11], type="b", lty=2)
abline(h=0.05,lty=4,col='blue')
title("sample size=20")
legend(0, 1, pch=c(rep(NA,3),1), lty=c(1,3,5,2), c("g__Streptococcus", "o__Lactobacillales", "p__Firmicutes", "FDR"))
plot(c(0, max(pct_change)), c(0,1), type="n", xlab="Percent of Perturbation", ylab="")
for (i in c(1,3,5)){
	lines(combo[1,12:22], node_freq[12:22,i], type="l", lty=i)
}
lines(combo[1,12:22],fdr[12:22], type="b", lty=2)
abline(h=0.05,lty=4,col='blue')
title("sample size=50")
plot(c(0, max(pct_change)), c(0,1), type="n", xlab="", ylab="")
for (i in c(1,3,5)){
	lines(combo[1,23:33], node_freq[23:33,i], type="l", lty=i)
}
lines(combo[1,23:33],fdr[23:33], type="b", lty=2)
abline(h=0.05,lty=4,col='blue')
title("sample size=100")
dev.off()


# for Figure 3(a)
pdf("tree_real.pdf", width=7, height=2.5)
par(mfrow=c(1,3), mar=c(4,4,2,0.5))
plot(combo[1,1:11],rej.uni.p[1:11], type="b", ylim=c(0,1), pch=8, 
     xlab="", ylab="Rejection Rate")
lines(combo[1,1:11],rej.p[1:11], type="b", pch=16)
lines(combo[1,1:11],rej.Fisher[1:11], type="b", pch=17)
lines(combo[1,1:11],rej[1:11], type="b", pch=7)
abline(h=0.05, lty=4, col="blue")
title("sample size=20")
legend(0, 1, pch=c(8,7,16,17), lty=c(1,1,1,1), cex=0.9,
	c("Rejection-PERMANOVA","Rejection-DM-2nd","Rejection-PairMN-2nd","Rejection-PairMN-Fisher"))
plot(combo[1,12:22],rej.uni.p[12:22], type="b", ylim=c(0,1), pch=8, 
     xlab="Percent of Perturbation", ylab="Rejection Rate")
lines(combo[1,12:22],rej.p[12:22], type="b", pch=16)
lines(combo[1,12:22],rej.Fisher[12:22], type="b", pch=17)
lines(combo[1,12:22],rej[12:22], type="b", pch=7)
abline(h=0.05, lty=4, col="blue")
title("sample size=50")
plot(combo[1,23:33],rej.uni.p[23:33], type="b", ylim=c(0,1), pch=8, 
     xlab="", ylab="Rejection Rate")
lines(combo[1,23:33],rej.p[23:33], type="b", pch=16)
lines(combo[1,23:33],rej.Fisher[23:33], type="b", pch=17)
lines(combo[1,23:33],rej[23:33], type="b", pch=7)
abline(h=0.05, lty=4, col="blue")
title("sample size=100")
dev.off()


## dense setting

load("result/res_simu_from_real2.Rdata")
combo <- matrix(0, 2, length(n_sample)*length(pct_change))
combo[1,] <- rep(pct_change, length(n_sample))
combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(pct_change))))

loc <- res_table$location
loc.01 <- res_table$location.01
ind <- ind_all <- c(ind1,ind2)
for (i in 1:length(ind)){
	ind_all <- c(ind_all, find_all_parent(ind[i],tree2))
}
ind_all <- unique(ind_all)
ind_all <- colnames(loc)[match(ind_all,colnames(loc))]
ind_all <- ind_all[!is.na(ind_all)]
node_freq <- loc[,ind_all]
node_freq <- node_freq/100

rej <- apply(res_table$pval<=0.05,1,mean)
rej.p <- apply(res_table$pval.p<=0.05,1,mean)
rej.Fisher <- apply(res_table$pval.Fisher<=0.05,1,mean)
rej.uni.p <- apply(res_table$pval.uni.p<=0.05,1,mean)

fdr_mean <- function(x){ x[is.na(x)]<-0; return(apply(x,1,mean)) }
fdr <- fdr_mean(res_table$fdr)

ind_node <- list(c(1,3,5),c(11,14,22),c(7,12,9,21))
ltycode <- c(1,3,5,6)


# for Figure 4(b)
pdf("nodes_real22.pdf", width=8, height=7.5)
par(mfrow=c(3,3), mar=c(4,4,2,0.5))
for (j in 1:3){
	ind_j <- ind_node[[j]]
	n_j <- length(ind_j)

	plot(c(0, max(pct_change)), c(0,1), type="n", xlab="", ylab="Percent of Discovery/FDR")
	for (i in 1:n_j){
		lines(combo[1,1:11], node_freq[1:11,ind_j[i]], type="l", lty=ltycode[i])
	}
	lines(combo[1,1:11],fdr[1:11], type="b", lty=2)
	abline(h=0.05,lty=4,col='blue')
	title("sample size=20")
	legend(0, 1, pch=c(rep(NA,length(ind_j)),1), lty=c(ltycode[1:n_j],2), c(colnames(node_freq[,ind_j]), "FDR"))
	plot(c(0, max(pct_change)), c(0,1), type="n", xlab="Percent of Perturbation", ylab="")
	for (i in 1:n_j){
		lines(combo[1,12:22], node_freq[12:22,ind_j[i]], type="l", lty=ltycode[i])
	}
	lines(combo[1,12:22],fdr[12:22], type="b", lty=2)
	abline(h=0.05,lty=4,col='blue')
	title("sample size=50")
	plot(c(0, max(pct_change)), c(0,1), type="n", xlab="", ylab="")
	for (i in 1:n_j){
		lines(combo[1,23:33], node_freq[23:33,ind_j[i]], type="l", lty=ltycode[i])
	}
	lines(combo[1,23:33],fdr[23:33], type="b", lty=2)
	abline(h=0.05,lty=4,col='blue')
	title("sample size=100")
}
dev.off()


# for Figure 3(b)
pdf("tree_real2.pdf", width=7, height=2.5)
par(mfrow=c(1,3), mar=c(4,4,2,0.5))
plot(combo[1,1:11],rej.uni.p[1:11], type="b", ylim=c(0,1), pch=8, 
     xlab="", ylab="Rejection Rate")
lines(combo[1,1:11],rej.p[1:11], type="b", pch=16)
lines(combo[1,1:11],rej.Fisher[1:11], type="b", pch=17)
lines(combo[1,1:11],rej[1:11], type="b", pch=7)
abline(h=0.05, lty=4, col="blue")
title("sample size=20")
plot(combo[1,12:22],rej.uni.p[12:22], type="b", ylim=c(0,1), pch=8, 
     xlab="Percent of Perturbation", ylab="Rejection Rate")
lines(combo[1,12:22],rej.p[12:22], type="b", pch=16)
lines(combo[1,12:22],rej.Fisher[12:22], type="b", pch=17)
lines(combo[1,12:22],rej[12:22], type="b", pch=7)
abline(h=0.05, lty=4, col="blue")
title("sample size=50")
plot(combo[1,23:33],rej.uni.p[23:33], type="b", ylim=c(0,1), pch=8, 
     xlab="", ylab="Rejection Rate")
lines(combo[1,23:33],rej.p[23:33], type="b", pch=16)
lines(combo[1,23:33],rej.Fisher[23:33], type="b", pch=17)
lines(combo[1,23:33],rej[23:33], type="b", pch=7)
abline(h=0.05, lty=4, col="blue")
title("sample size=100")
dev.off()



############################################################################
######################  codes for real data analyses  ######################
############################################################################


snsdat <- read.csv("dataset/cts_smoker.csv", row.names=1, check.names=F)
p <- ncol(snsdat)-2
cts3 <-snsdat[,1:p]
nonsmoker <- as.numeric(snsdat[,p+1])
bodysite <- as.numeric(snsdat[,p+2])
tree <- read.csv("dataset/tree_smoker.csv", row.names=1, check.names=F)
tree <- as.matrix(tree[,-3])
rownames(tree) <- NULL
node_name <- colnames(cts3)
tree2 <- cbind(node_name, node_name)
colnames(tree2) <- colnames(tree)
tree2 <- as.data.frame(rbind(tree, tree2))

sampleID <- rownames(cts3)
temp <- strsplit(sampleID,".",fixed=T)
subjectID <- as.character(lapply(temp, function(x) x[[3]]))


#### paired test of nasopharynx vs oropharynx, left side, nonsmoker ####
temp1 <- bodysite==1 & nonsmoker==1
temp2 <- bodysite==3 & nonsmoker==1
set1 <- cts3[temp1,]
set2 <- cts3[temp2,]
ID1 <- subjectID[temp1]
ID2 <- subjectID[temp2]
IDshare <- ID1[match(ID2,ID1)]
IDshare <- IDshare[!is.na(IDshare)]
set1 <- set1[match(IDshare, ID1),]
set2 <- set2[match(IDshare, ID2),]

res <- CorM_tree_test(set1, set2, tree2, paired=F, elim=1)
res$parent[!is.na(res$qvalue)&res$qvalue<=0.05]
res$pvalue.combined
res$pvalue.Fisher
nodes <- res$parent[!is.na(res$qvalue)&res$qvalue<=0.05]

res3 <- test_unifrac(set1, set2, tree2, FALSE)
res3[1,6]

# for Figure 5(a)
pdf('barplot_left.pdf',7,4)
barplot_nodes(nodes, set1, set2, tree2, 
    grp.names=c('Nasopharynx', 'Oropharynx'), elim=1, common=0.15)
dev.off()


#### unpaired test for smoker vs nonsmoker, nasopharynx ####
set1 <- cts3[nonsmoker==1&bodysite<=2,]
set2 <- cts3[nonsmoker==2&bodysite<=2,]

res <- CorM_tree_test(set1, set2, tree2, F, 1)
res$pvalue.combined
res$pvalue.Fisher
res$parent[!is.na(res$qvalue)&res$qvalue<=0.05]
nodes <- res$parent[!is.na(res$qvalue)&res$qvalue<=0.05]

res2 <- test_unifrac(set1, set2, tree2, FALSE)
res2[1,6]

# for Figure 5(b)
pdf('barplot_naso2.pdf',7,4)
barplot_nodes(nodes, set1, set2, tree2, 
    grp.names=c('Nonsmoker', 'Smoker'), elim=1, common=0.06)
dev.off()


#### unpaired test for smoker vs nonsmoker, oropharynx ####
set1 <- cts3[nonsmoker==1&bodysite>=3,]
set2 <- cts3[nonsmoker==2&bodysite>=3,]
res <- CorM_tree_test(set1, set2, tree2, F, 1)
res$pvalue.combined
res$pvalue.Fisher
res$parent[!is.na(res$qvalue)&res$qvalue<=0.05]
nodes <- res$parent[!is.na(res$qvalue)&res$qvalue<=0.05]

res2 <- test_unifrac(set1, set2, tree2, FALSE)
res2[1,6]

# for Figure 5(c)
pdf('barplot_oro2.pdf',7,4)
barplot_nodes(nodes, set1, set2, tree2, 
    grp.names=c('Nonsmoker', 'Smoker'), elim=1, common=0.06)
dev.off()


