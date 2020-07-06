
#### Simulate single level counts data ####

rdirichlet <- function(n, alpha){
	k <- length(alpha)
	datmat <- matrix(0, n, k)
	for (i in 1:k){
		datmat[,i] <- rgamma(n, alpha[i])
	}
	datmat <- datmat/rowSums(datmat)
	return(datmat)
}


rDM <- function(n, size, alpha){
	# size need to be length n, rate sums to 1
	if (length(size)!=n) stop("Length of size is not n!")
	k <- length(alpha)
	xx <- matrix(0,n,k)
	pp <- rdirichlet(n, alpha)
	for (i in 1:n){
		xx[i,] <- rmultinom(1, size[i], pp[i,])
	}
	return(xx)
}


rDM_Paired <- function(n, size1, size2, alpha1, alpha2, ell, rho){
	# size need to be length n, rate sums to 1
	if (length(size1)!=n | length(size2)!=n) stop("Length of size is not n!")
	k1 <- length(alpha1); k2 <- length(alpha2)
	k <- length(ell)
	if(k1!=k | k2!=k) stop("Dimension of alpha and ell don't match!")

	p12 <- rdirichlet(n, ell)
	pp1 <- (1-rho)*rdirichlet(n, alpha1) + rho*p12
	pp2 <- (1-rho)*rdirichlet(n, alpha2) + rho*p12

	xx1 <- matrix(0,n,k)
	xx2 <- matrix(0,n,k)
	for (i in 1:n){
		xx1[i,] <- rmultinom(1, size1[i], pp1[i,])
		xx2[i,] <- rmultinom(1, size2[i], pp2[i,])
	}
	return(list(sample1 = xx1, sample2 = xx2))
}


rNormM <- function(n, size1, size2, mu1, mu2, sigma, rho){
	# size need to be length n
	if (length(size1)!=n | length(size2)!=n) stop("Length of size is not n!")
	k <- length(sigma)
	if(length(mu1)!=k | length(mu2)!=k) stop("mu1, mu2 and sigma do not have same length!")

	Sig <- matrix(c(1,rho,rho,1),2,2)
	pp1 <- matrix(0,n,k)
	pp2 <- matrix(0,n,k)
	for (j in 1:k){
		pp <- mvrnorm(n, c(mu1[j], mu2[j]), sigma[j]*Sig)
		pp1[,j] <- pp[,1]
		pp2[,j] <- pp[,2]
	}
	pp1 <- exp(pp1); pp1 <- pp1/rowSums(pp1)
	pp2 <- exp(pp2); pp2 <- pp2/rowSums(pp2)

	xx1 <- matrix(0,n,k)
	xx2 <- matrix(0,n,k)
	for (i in 1:n){
		xx1[i,] <- rmultinom(1, size1[i], pp1[i,])
		xx2[i,] <- rmultinom(1, size2[i], pp2[i,])
	}
	return(list(sample1 = xx1, sample2 = xx2))
}


#### Simulate tree data ####

simu_DMtree <- function(tree2, alpha, n_sample){
	# names of all the nodes
	node_name <- unique(as.vector(as.matrix(tree2[,1:2])))
	# names of all the non-leaf nodes
	parent <- unique(tree2$parent[tree2$parent!=tree2$child])
	# number of reads per sample
	N_reads <- c(rpois(n_sample/4, 10000), rpois(n_sample/4, 8000), 
				rpois(n_sample/4, 6000), rpois(n_sample/4, 4000))
	# initialize the data matrix
	dm_sub <- matrix(0, n_sample, length(node_name))
	colnames(dm_sub) <- node_name
	# assign reads to root
	dm_sub[,node_name==1] <- N_reads
	for (j in parent){
		alp <- alpha[tree2$parent==j]
		ind_child <- find_child(j,tree2)
		temp <- rDM(n_sample, dm_sub[,node_name==j], alp)
		dm_sub[,ind_child[ind_child!=j]] <- temp[,ind_child!=j]
	}
	return(dm_sub)
}


simu_DMtree_Paired <- function(tree2, alpha1, alpha2, ell, rho, n_sample){
	# names of all the nodes
	node_name <- unique(as.vector(as.matrix(tree2[,1:2])))
	# names of all the non-leaf nodes
	parent <- unique(tree2$parent[tree2$parent!=tree2$child])
	# number of reads per sample
	N_reads1 <- c(rpois(n_sample/5, 10000), rpois(n_sample/5, 8000), 
				rpois(n_sample/5, 6000), rpois(n_sample/5, 4000), 
				rpois(n_sample/5, 2000))
	N_reads2 <- c(rpois(n_sample/5, 10000), rpois(n_sample/5, 8000), 
				rpois(n_sample/5, 6000), rpois(n_sample/5, 4000),
				rpois(n_sample/5, 2000))
	# initialize the data matrix
	dm_sub1 <- matrix(0, n_sample, length(node_name))
	dm_sub2 <- matrix(0, n_sample, length(node_name))
	colnames(dm_sub1) <- node_name
	colnames(dm_sub2) <- node_name
	# assign reads to root
	dm_sub1[,node_name==1] <- N_reads1
	dm_sub2[,node_name==1] <- N_reads2
	for (j in parent){
		ind_child <- find_child(j,tree2)
        ind_j <- tree2$parent==j
		alp1 <- alpha1[ind_j]
		alp2 <- alpha2[ind_j]
		ell0 <- ell[ind_j]
		temp <- rDM_Paired(n_sample, dm_sub1[,node_name==j], dm_sub2[,node_name==j], 
						alp1, alp2, ell0, rho)
                        ind_c_j <- ind_child!=j
		dm_sub1[,ind_child[ind_c_j]] <- temp$sample1[,ind_c_j]
		dm_sub2[,ind_child[ind_c_j]] <- temp$sample2[,ind_c_j]
	}
	return(list(dataset1 = dm_sub1, dataset2 = dm_sub2))
}


simu_NormMtree <- function(tree2, mu1, mu2, sigma, rho, n_sample){
	# names of all the nodes
	node_name <- unique(as.vector(as.matrix(tree2[,1:2])))
	# names of all the non-leaf nodes
	parent <- unique(tree2$parent[tree2$parent!=tree2$child])
	# number of reads per sample
	N_reads1 <- c(rpois(n_sample/5, 10000), rpois(n_sample/5, 8000), 
				rpois(n_sample/5, 6000), rpois(n_sample/5, 4000), 
				rpois(n_sample/5, 2000))
	N_reads2 <- c(rpois(n_sample/5, 10000), rpois(n_sample/5, 8000), 
				rpois(n_sample/5, 6000), rpois(n_sample/5, 4000),
				rpois(n_sample/5, 2000))
	# initialize the data matrix
	dm_sub1 <- matrix(0, n_sample, length(node_name))
	dm_sub2 <- matrix(0, n_sample, length(node_name))
	colnames(dm_sub1) <- node_name
	colnames(dm_sub2) <- node_name
	# assign reads to root
	dm_sub1[,node_name==1] <- N_reads1
	dm_sub2[,node_name==1] <- N_reads2
	for (j in parent){
		ind_child <- find_child(j,tree2)
        ind_j <- tree2$parent==j
		m1 <- mu1[ind_j]
		m2 <- mu2[ind_j]
		sig <- sigma[ind_j]
		temp <- rNormM(n_sample, dm_sub1[,node_name==j], dm_sub2[,node_name==j], 
				m1, m2, sig, rho)
        ind_c_j <- ind_child!=j
		dm_sub1[,ind_child[ind_c_j]] <- temp$sample1[,ind_c_j]
		dm_sub2[,ind_child[ind_c_j]] <- temp$sample2[,ind_c_j]
	}
	return(list(dataset1 = dm_sub1, dataset2 = dm_sub2))
}


#### Simulate tables for comparison ####

simu_table_DM <- function(iter, n_sample, rate1, rate2, theta1, theta2, ell, rho, seed){
	set.seed(seed)
	combo <- matrix(0, 2, length(n_sample)*length(rho))
	combo[1,] <- rep(rho, length(n_sample))
	combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(rho))))
	n_combo <- ncol(combo)
	pval <- matrix(0, n_combo, iter)
	pval.p <- matrix(0, n_combo, iter)
	for (ii in 1:n_combo){
		rho_i <- combo[1,ii]
		n_i <- combo[2,ii]
		cat("rho=", rho_i, ", n=", n_i, "\n", sep="")
		alpha1 <- (rate1-rho_i*ell)/(1-rho_i)*theta1
		alpha2 <- (rate2-rho_i*ell)/(1-rho_i)*theta2
		for (i in 1:iter){
			size1 <- rpois(n_i,1000); size2 <- rpois(n_i,1000)
			sample <- rDM_Paired(n_i, size1, size2, alpha1, alpha2, ell, rho_i)
			res <- CorM_test(sample$sample1, sample$sample2, FALSE)
			res.p <- CorM_test(sample$sample1, sample$sample2, TRUE)
			pval[ii, i] <- res$pvalue
			pval.p[ii, i] <- res.p$pvalue
		}
	}
	return(list(pval=pval, pval.p=pval.p, combo=combo))
}


simu_table_tree_DM <- function(iter, n_sample, tree2, rate1, rate2, theta1, theta2, ell, rho, 
						seed, filename){
	set.seed(seed)

	combo <- matrix(0, 2, length(n_sample)*length(rho))
	combo[1,] <- rep(rho, length(n_sample))
	combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(rho))))
	n_combo <- ncol(combo)
	diff.node <- c(4,9)

	pval <- matrix(0, n_combo, iter)
	pval.p <- matrix(0, n_combo, iter)
	pval.Fisher <- matrix(0, n_combo, iter)
	pval.uni <- matrix(0, n_combo, iter)
	parent <- as.vector(unique(tree2$parent[tree2$parent!=tree2$child]))
	loc <- matrix(0, n_combo, length(parent))
	loc.01 <- matrix(0, n_combo, length(parent))
	fdr <- matrix(0, n_combo, iter)
	fdr.01 <- matrix(0, n_combo, iter)

	for (ii in 1:n_combo){
		rho_i <- combo[1,ii]
		n_i <- combo[2,ii]
		cat("rho=", rho_i, ", n=", n_i, "\n", sep="")
		alpha1 <- (rate1-rho_i*ell)/(1-rho_i)*theta1
		alpha2 <- (rate2-rho_i*ell)/(1-rho_i)*theta2
		for (i in 1:iter){
			sets <- simu_DMtree_Paired(tree2, alpha1, alpha2, ell, rho_i, n_i)
			res <- CorM_tree_test(sets$dataset1, sets$dataset2, tree2, 1, FALSE)
			pval[ii, i] <- res$pvalue.combined
			res.p <- CorM_tree_test(sets$dataset1, sets$dataset2, tree2, 1, TRUE)
			pval.p[ii, i] <- res.p$pvalue.combined
			qval.p <- res.p$qval
			pval.Fisher[ii, i] <- res.p$pvalue.Fisher
			# res.uni <- test_unifrac(sets$dataset1, sets$dataset2, tree2, TRUE)
			# pval.uni[ii, i] <- res.uni[1,6]
			sel <- as.numeric((qval.p<0.05)&(!is.na(qval.p)))
			loc[ii,] <- loc[ii,] + sel
			fdr[ii, i] <- sum(sel[-diff.node])/sum(sel)
			sel.01 <- as.numeric((qval.p<0.01)&(!is.na(qval.p)))
			loc.01[ii,] <- loc.01[ii,] + sel.01
			fdr.01[ii, i] <- sum(sel.01[-diff.node])/sum(sel.01)
		}
		res_table <- list(pval=pval, pval.p=pval.p, pval.Fisher=pval.Fisher, pval.uni=pval.uni, 
					location=loc, location.01=loc.01, fdr=fdr, fdr.01=fdr.01, combo=combo)
		save(res_table, file=filename)
	}
}



simu_table_NormM <- function(iter, n_sample, mu1, mu2, sigma, rho, seed){
	set.seed(seed)
	combo <- matrix(0, 2, length(n_sample)*length(rho))
	combo[1,] <- rep(rho, length(n_sample))
	combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(rho))))
	n_combo <- ncol(combo)
	pval <- matrix(0, n_combo, iter)
	pval.p <- matrix(0, n_combo, iter)
	for (ii in 1:n_combo){
		rho_i <- combo[1,ii]
		n_i <- combo[2,ii]
		cat("rho=", rho_i, ", n=", n_i, "\n", sep="")
		for (i in 1:iter){
			size1 <- rpois(n_i,1000); size2 <- rpois(n_i,1000)
			sample <- rNormM(n_i, size1, size2, mu1, mu2, sigma, rho_i)
			res <- CorM_test(sample$sample1, sample$sample2, FALSE)
			res.p <- CorM_test(sample$sample1, sample$sample2, TRUE)
			pval[ii, i] <- res$pvalue
			pval.p[ii, i] <- res.p$pvalue
		}
	}
	return(list(pval=pval, pval.p=pval.p, combo=combo))
}



simu_table_tree_NormM <- function(iter, n_sample, tree2, mu1, mu2, sigma, rho, 
							seed, filename){
	set.seed(seed)
 
	combo <- matrix(0, 2, length(n_sample)*length(rho))
	combo[1,] <- rep(rho, length(n_sample))
	combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(rho))))
	n_combo <- ncol(combo)
	diff.node <- c(4,9)

	pval <- matrix(0, n_combo, iter)
	pval.p <- matrix(0, n_combo, iter)
	pval.Fisher <- matrix(0, n_combo, iter)
	pval.uni <- matrix(0, n_combo, iter)
	parent <- as.vector(unique(tree2$parent[tree2$parent!=tree2$child]))
	loc <- matrix(0, n_combo, length(parent))
	loc.01 <- matrix(0, n_combo, length(parent))
	fdr <- matrix(0, n_combo, iter)
	fdr.01 <- matrix(0, n_combo, iter)

	for (ii in 1:n_combo){
		rho_i <- combo[1,ii]
		n_i <- combo[2,ii]
		cat("rho=", rho_i, ", n=", n_i, "\n", sep="")
		for (i in 1:iter){
			sets <- simu_NormMtree(tree2, mu1, mu2, sigma, rho_i, n_i)
			res <- CorM_tree_test(sets$dataset1, sets$dataset2, tree2, 1, FALSE)
			pval[ii, i] <- res$pvalue.combined
			res.p <- CorM_tree_test(sets$dataset1, sets$dataset2, tree2, 1, TRUE)
			pval.p[ii, i] <- res.p$pvalue.combined
			qval.p <- res.p$qval
			pval.Fisher[ii, i] <- res.p$pvalue.Fisher
			# res.uni <- test_unifrac(sets$dataset1, sets$dataset2, tree2, TRUE)
			# pval.uni[ii, i] <- res.uni[1,6]
			sel <- as.numeric((qval.p<0.05)&(!is.na(qval.p)))
			loc[ii,] <- loc[ii,] + sel
			fdr[ii, i] <- sum(sel[-diff.node])/sum(sel)
			sel.01 <- as.numeric((qval.p<0.01)&(!is.na(qval.p)))
			loc.01[ii,] <- loc.01[ii,] + sel.01
			fdr.01[ii, i] <- sum(sel.01[-diff.node])/sum(sel.01)
		}
		res_table <- list(pval=pval, pval.p=pval.p, pval.Fisher=pval.Fisher, pval.uni=pval.uni, 
					location=loc, location.01=loc.01, fdr=fdr, fdr.01=fdr.01, combo=combo)
		save(res_table, file=filename)
	}
}




#### simulate from real data ####

resample_count2 <- function(mat, prob){
	mat <- as.matrix(mat)
	resample_single <- function(x){ return(rbinom(1,x,prob)) }
	mat2 <- sapply(mat, resample_single)
	mat2 <- matrix(mat2, nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat)))
	return(mat2)
}

resample_count <- function(n_sample, prp_base, ttl_ct, overlap=0.5){
	n_node <- ncol(prp_base); n_base <- nrow(prp_base)
	ttl1 <- sample(ttl_ct, n_sample)
	ttl2 <- sample(ttl_ct, n_sample)
	rr <- sample.int(n_base, n_sample*3, replace=T)
	prp1 <- prp_base[rr[1:n_sample],] #+ 
			#prp_base[rr[(n_sample+1):(n_sample*2)],] + 
			#prp_base[rr[(n_sample*2+1):(n_sample*3)],]
	prp2 <- prp1
	rr <- sample.int(n_base, n_sample*3, replace=T)
	prp1 <- prp1*overlap + prp_base[rr[1:n_sample],]*(1-overlap) #+ 
			#prp_base[rr[(n_sample+1):(n_sample*2)],] + 
			#prp_base[rr[(n_sample*2+1):(n_sample*3)],]
	rr <- sample.int(n_base, n_sample*3, replace=T)
	prp2 <- prp2*overlap + prp_base[rr[1:n_sample],]*(1-overlap) #+ 
			#prp_base[rr[(n_sample+1):(n_sample*2)],] + 
			#prp_base[rr[(n_sample*2+1):(n_sample*3)],]
	
	temp <- rmultinom(1, sum(ttl1), as.vector(as.matrix(prp1)*ttl1))
	set1 <- matrix(temp, n_sample, n_node, dimnames = list(rownames(prp1), colnames(prp_base)))
	temp <- rmultinom(1, sum(ttl2), as.vector(as.matrix(prp2)*ttl2))
	set2 <- matrix(temp, n_sample, n_node, dimnames = list(rownames(prp2), colnames(prp_base)))	
	
	return(list(dataset1 = set1, dataset2 = set2))
}

simu_table_from_real <- function(iter, n_sample, tree2, ttl_ct, prp_base, ind, pct_change,
							seed, filename){
	set.seed(seed)

	combo <- matrix(0, 2, length(n_sample)*length(pct_change))
	combo[1,] <- rep(pct_change, length(n_sample))
	combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(pct_change))))
	n_combo <- ncol(combo)
	
	n_node <- ncol(prp_base)

	pval <- matrix(0, n_combo, iter)
	pval.p <- matrix(0, n_combo, iter)
	pval.Fisher <- matrix(0, n_combo, iter)
	pval.uni <- matrix(0, n_combo, iter)
	pval.uni.p <- matrix(0, n_combo, iter)
	parent <- as.vector(unique(tree2$parent[tree2$parent!=tree2$child]))
	loc <- matrix(0, n_combo, length(parent))
	loc.01 <- matrix(0, n_combo, length(parent))
	fdr <- matrix(0, n_combo, iter)
	fdr.01 <- matrix(0, n_combo, iter)
	colnames(loc) <- colnames(loc.01) <- parent
	diff.node <- match(find_all_parent(ind,tree2), parent)

	for (ii in 1:n_combo){
		pct_i <- combo[1,ii]
		n_i <- combo[2,ii]
		cat("pct_change=", pct_i, ", n=", n_i, "\n", sep="")
		for (i in 1:iter){
			print(i)
			sets <- resample_count(n_i, prp_base, ttl_ct)
			set1 <- sets$dataset1
			set2 <- sets$dataset2
			set2[,ind] <- set2[,ind] + resample_count2(rowSums(set2), pct_i)
			# set1 <- resample_count(ori_set1, divi)
			# set2[,-ind] <- resample_count(ori_set2[,-ind], divi)
			# set2[,ind] <- resample_count(ori_set2[,ind], pct_i*divi)
			set1 <- sub2total(set1, tree2)
			set2 <- sub2total(set2, tree2)
			res <- CorM_tree_test(set1, set2, tree2, 1, F)
			pval[ii, i] <- res$pvalue.combined
			res.p <- CorM_tree_test(set1, set2, tree2, 1, T)
			pval.p[ii, i] <- res.p$pvalue.combined
			qval.p <- res.p$qval
			pval.Fisher[ii, i] <- res.p$pvalue.Fisher
			res.uni <- test_unifrac(set1, set2, tree2, FALSE)
			pval.uni[ii, i] <- res.uni[1,6]
			res.uni.p <- test_unifrac(set1, set2, tree2, TRUE)
			pval.uni.p[ii, i] <- res.uni.p[1,6]
			sel <- as.numeric((qval.p<0.05)&(!is.na(qval.p)))
			loc[ii,] <- loc[ii,] + sel
			fdr[ii, i] <- sum(sel[-diff.node])/sum(sel)
			sel.01 <- as.numeric((qval.p<0.01)&(!is.na(qval.p)))
			loc.01[ii,] <- loc.01[ii,] + sel.01
			fdr.01[ii, i] <- sum(sel.01[-diff.node])/sum(sel.01)
			res_table <- list(pval=pval, pval.p=pval.p, pval.Fisher=pval.Fisher, pval.uni=pval.uni, pval.uni.p=pval.uni.p, 
							location=loc, location.01=loc.01, fdr=fdr, fdr.01=fdr.01, combo=pct_change)
			save(res_table, file=filename)
		}
	}
}



simu_table_from_real2 <- function(iter, n_sample, tree2, ttl_ct, prp_base, ind1, ind2, pct_change,
							seed, filename){
	combo <- matrix(0, 2, length(n_sample)*length(pct_change))
	combo[1,] <- rep(pct_change, length(n_sample))
	combo[2,] <- as.vector(t(matrix(n_sample, length(n_sample), length(pct_change))))
	n_combo <- ncol(combo)
	
	n_node <- ncol(prp_base)

	pval <- matrix(0, n_combo, iter)
	pval.p <- matrix(0, n_combo, iter)
	pval.Fisher <- matrix(0, n_combo, iter)
	pval.uni <- matrix(0, n_combo, iter)
	pval.uni.p <- matrix(0, n_combo, iter)
	parent <- as.vector(unique(tree2$parent[tree2$parent!=tree2$child]))
	loc <- matrix(0, n_combo, length(parent))
	loc.01 <- matrix(0, n_combo, length(parent))
	fdr <- matrix(0, n_combo, iter)
	fdr.01 <- matrix(0, n_combo, iter)
	colnames(loc) <- colnames(loc.01) <- parent
    ind <- c(ind1,ind2)
    ind_all <- NULL
    for (i in 1:length(ind)){
        ind_all <- c(ind_all,find_all_parent(ind[i],tree2))
    }
    ind_all <- unique(ind_all)
	diff.node <- c(match(ind_all, parent))
    diff.node <- diff.node[!is.na(diff.node)]
	for (ii in 1:n_combo){
        set.seed(seed+ii)
		pct_i <- combo[1,ii]
		n_i <- combo[2,ii]
		cat("pct_change=", pct_i, ", n=", n_i, "\n", sep="")
		for (i in 1:iter){
			print(i)
			sets <- resample_count(n_i, prp_base, ttl_ct)
			set1 <- sets$dataset1
			for (j in 1:length(ind1)){
				set1[,ind1[j]] <- set1[,ind1[j]] + resample_count2(rowSums(set1), pct_i)
			}
			set2 <- sets$dataset2
			for (j in 1:length(ind2)){
				set2[,ind2[j]] <- set2[,ind2[j]] + resample_count2(rowSums(set2), pct_i)
			}
			# set1 <- resample_count(ori_set1, divi)
			# set2[,-ind] <- resample_count(ori_set2[,-ind], divi)
			# set2[,ind] <- resample_count(ori_set2[,ind], pct_i*divi)
			set1 <- sub2total(set1, tree2)
			set2 <- sub2total(set2, tree2)
			res <- CorM_tree_test(set1, set2, tree2, 1, F)
			pval[ii, i] <- res$pvalue.combined
			res.p <- CorM_tree_test(set1, set2, tree2, 1, T)
			pval.p[ii, i] <- res.p$pvalue.combined
			qval.p <- res.p$qval
			pval.Fisher[ii, i] <- res.p$pvalue.Fisher
			res.uni <- test_unifrac(set1, set2, tree2, FALSE)
			pval.uni[ii, i] <- res.uni[1,6]
			res.uni.p <- test_unifrac(set1, set2, tree2, TRUE)
			pval.uni.p[ii, i] <- res.uni.p[1,6]
			sel <- as.numeric((qval.p<0.05)&(!is.na(qval.p)))
			loc[ii,] <- loc[ii,] + sel
			fdr[ii, i] <- sum(sel[-diff.node])/sum(sel)
			sel.01 <- as.numeric((qval.p<0.01)&(!is.na(qval.p)))
			loc.01[ii,] <- loc.01[ii,] + sel.01
			fdr.01[ii, i] <- sum(sel.01[-diff.node])/sum(sel.01)
			res_table <- list(pval=pval, pval.p=pval.p, pval.Fisher=pval.Fisher, pval.uni=pval.uni, pval.uni.p=pval.uni.p, 
							location=loc, location.01=loc.01, fdr=fdr, fdr.01=fdr.01, combo=pct_change)
			save(res_table, file=filename)
		}
	}
}
