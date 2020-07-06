library(qiimer) # For reading QIIME files
library(vegan) # Ecology
library(ape) # For PCoA
library(reshape2) # Like pivot table in MS Excel
library(ggplot2) # Plots

#### Testing ####

CorM_test <- function(sample1, sample2, paired=FALSE){
	if(paired){
		res <- Wald_paired(sample1, sample2)
	}else{
		res <- Wald_unpaired(sample1, sample2)
	}
	return(res)
}


Wald_unpaired <- function(sample1, sample2){
	k <- ncol(sample1)
	k2 <- ncol(sample2)
	if (k!=k2) stop("The two sets have different number of columns!")
	n1 <- nrow(sample1)
	n2 <- nrow(sample2)
	N1 <- rowSums(sample1)
	N2 <- rowSums(sample2)
	NN1 <- sum(N1)
	NN2 <- sum(N2)
	Nc1 <- (NN1-sum(N1^2)/NN1)/(n1-1)
	Nc2 <- (NN2-sum(N2^2)/NN2)/(n2-1)
	m1 <- colSums(sample1)/NN1
	m2 <- colSums(sample2)/NN2
	mi1 <- sample1/N1
	mi2 <- sample2/N2
	S1 <- colSums((t(t(mi1)-m1))^2*N1) / (n1-1)
	S2 <- colSums((t(t(mi2)-m2))^2*N2) / (n2-1)
	G1 <- colSums(mi1*(1-mi1)*N1) / (NN1-n1)
	G2 <- colSums(mi2*(1-mi2)*N2) / (NN2-n2)	
	theta1 <- sum(S1-G1)/sum(S1+(Nc1-1)*G1)
	theta2 <- sum(S2-G2)/sum(S2+(Nc2-1)*G2)
	theta1 <- max(theta1, 0); theta1 <- min(theta1, 1)
	theta2 <- max(theta2, 0); theta2 <- min(theta2, 1)
	con1 <- theta1*(sum(N1^2)/NN1^2) + (1-theta1)/NN1
	con2 <- theta2*(sum(N2^2)/NN2^2) + (1-theta2)/NN2
	#using new mean estimate
	#m1 <- colMeans(mi1)
	#m2 <- colMeans(mi2)
	#using handsolved general inverse
	stat.Wilson <- sum((m1-m2)^2/(con1*m1+con2*m2))
	stat.Wald <- stat.Wilson + stat.Wilson^2/(1/con1+1/con2-stat.Wilson)
	pval.Wilson <- 1-pchisq(stat.Wilson, k-1)
	pval.Wald <- 1-pchisq(stat.Wald, k-1)
	#using general inverse
	#Smat <- diag(con1*m1+con2*m2)-con1*m1%o%m1-con2*m2%o%m2
	#stat <- t(m1-m2)%*%ginv(Smat)%*%(m1-m2)

	#return(list(stat.Wald=stat.Wald, stat.Wilson=stat.Wilson, 
	#		pvalue.Wald=pval.Wald, pvalue.Wilson=pval.Wilson))
	return(list(stat=stat.Wilson, pvalue=pval.Wilson))
}


Wald_paired <- function(sample1, sample2){
	k <- ncol(sample1)
	k2 <- ncol(sample2)
	if (k!=k2) stop("The two sets have different number of columns!")
	n1 <- nrow(sample1)
	n2 <- nrow(sample2)
	if (n1!=n2) stop("The two sets are not paired!")
	N1 <- rowSums(sample1)
	N2 <- rowSums(sample2)
	NN1 <- sum(N1)
	NN2 <- sum(N2)
	eyes <- diag(rep(1,k))
	if(NN1==n1 | NN2==n2 | n1<k){
		stat <- NA
		pval <- NA
	}else{
		Nc1 <- (NN1-sum(N1^2)/NN1)/(n1-1)
		Nc2 <- (NN2-sum(N2^2)/NN2)/(n2-1)

		m1 <- colSums(sample1)/NN1
		m2 <- colSums(sample2)/NN2
		mi1 <- sample1/N1
		mi2 <- sample2/N2
		diff1 <- t(mi1)-m1
		diff2 <- t(mi2)-m2
		S1 <- diff1%*%diag(N1)%*%t(diff1)/(n1-1)
		S2 <- diff2%*%diag(N2)%*%t(diff2)/(n2-1)
		G1 <- (diag(as.vector(N1%*%mi1)) - t(mi1)%*%diag(N1)%*%mi1) / (NN1-n1)
		G2 <- (diag(as.vector(N2%*%mi2)) - t(mi2)%*%diag(N2)%*%mi2) / (NN2-n2)
		V <- diff1%*%diag(N1+N2)%*%t(diff2) / (n1-1)
		V <- (V + t(V))/(Nc1+Nc2)*sum(N1*N2)/NN1/NN2
		#Sig1 <- (S1+(Nc1-1)*G1)/Nc1/NN1 + (S1-G1) * (sum(N1^2)-NN1)/Nc1/NN1^2
		#Sig2 <- (S2+(Nc2-1)*G2)/Nc2/NN2 + (S2-G2) * (sum(N2^2)-NN2)/Nc2/NN2^2
		Sig1 <- S1*sum(N1^2)/Nc1/NN1^2 + G1*(Nc1-sum(N1^2)/NN1)/NN1/Nc1
		Sig2 <- S2*sum(N2^2)/Nc2/NN2^2 + G2*(Nc2-sum(N2^2)/NN2)/NN2/Nc2
		Sig1 <- trunc_matrix(Sig1)
		Sig2 <- trunc_matrix(Sig2)

		Smat <- Sig1 + Sig2 - V
		Smat <- trunc_matrix(Smat)

		stat <- t(m1-m2)%*%ginv(Smat)%*%(m1-m2)*(n1-k+1)/(n1-1)/(k-1)
		pval <- 1-pf(stat, k-1, n1-k+1)
	}
	return(list(stat=stat, pvalue=pval))
}


CorM_tree_test <- function(set1, set2, tree2, paired=FALSE, elim=1, comparable=FALSE){
	# names of all the non-leaf nodes
	parent <- as.vector(unique(tree2$parent[tree2$parent!=tree2$child]))
	# store the pvalues and test statistics
	pval <- rep(0, length(parent))
	stat <- rep(0, length(parent))
	na <- rep(0, length(parent))

	for (i in 1:length(parent)){
		j=parent[i]
		chd <- find_child(j,tree2)
		sample1 <- set1[,chd]
		sample2 <- set2[,chd]
        chd_j <- chd==j
		sample1[,chd_j] <- sample1[,chd_j] - rowSums(as.matrix(sample1[,!chd_j]))
		sample2[,chd_j] <- sample2[,chd_j] - rowSums(as.matrix(sample2[,!chd_j]))
		# vars with small total counts are elimniated
		zeros <- colMeans(sample1)>elim & colMeans(sample2)>elim
		# obs with small total counts are eliminated
		ind1 <- rowMeans(as.matrix(sample1[,zeros]))>elim
		ind2 <- rowMeans(as.matrix(sample2[,zeros]))>elim
		k <- sum(zeros)
		if(paired | comparable){
			ind12 <- ind1&ind2
			# check if there's enough degree of freedom
			if (sum(ind12)>=k & k>1){
				sample1 <- sample1[,zeros]
				sample2 <- sample2[,zeros]
				sample1 <- as.matrix(sample1[ind12,])
				sample2 <- as.matrix(sample2[ind12,])
				res <- CorM_test(sample1, sample2, paired)
				pval[i] <- res$pvalue
				stat[i] <- res$stat
			}else{
				pval[i] <- NA
				stat[i] <- NA
				na[i] <- ifelse(sum(ind12)<k, "obs<col", "one_col")
			}
		}else{
			if (sum(ind1)>=k & sum(ind2)>=k & k>1){
				sample1 <- sample1[,zeros]
				sample2 <- sample2[,zeros]
				sample1 <- as.matrix(sample1[ind1,])
				sample2 <- as.matrix(sample2[ind2,])
				res <- CorM_test(sample1, sample2, FALSE)
				pval[i] <- res$pvalue
				stat[i] <- res$stat
			}else{
				pval[i] <- NA
				stat[i] <- NA
				na[i] <- ifelse(sum(ind1)<k | sum(ind2)<k, "obs<col", "one_col")
			} 
		} 
	}
	pval.ranked <- sort(pval)
	pval.n <- length(pval.ranked)
	qval <- p.adjust(pval, 'BH') 
	pval.combined <- 1 - (1-pval.ranked[1])^pval.n
	if (pval.n>1){
		pval.combined <- 1 - ((1-pval.ranked[2])^(pval.n-1))*(1+(pval.n-1)*pval.ranked[2])
	}
	pval.Fisher <- 1 - pchisq(-2*sum(log(pval.ranked)), 2*pval.n) 
	return(list(qvalue=qval, pvalue=pval, stat=stat, pvalue.combined=pval.combined, pvalue.Fisher=pval.Fisher, na = na, parent=parent))
}



trunc_matrix <- function(Smatrix){
	eSmat <- eigen(Smatrix)
	eSmat$values[eSmat$values<0] <- 0
	Smatrix <- eSmat$vectors%*%diag(eSmat$values)%*%t(eSmat$vectors)
	return(Smatrix)
}


test_unifrac <- function(set1, set2, tree2, paired=FALSE){
	ind <- c(rep(1, nrow(set1)), rep(2, nrow(set2)))
	sets <- rbind(set1, set2)
	distmat <- get_dist(sets, sets, tree2)
	if (paired){
		stratum <- c(c(1:nrow(set1)),c(1:nrow(set1)))
		res <- adonis(distmat ~ ind, strata=stratum)$aov.tab
	}else{
		res <- adonis(distmat ~ ind)$aov.tab
	}
	return(res)
}


#### tree operation functions ####

get_dist <- function(set1, set2, tree2){
	distmat <- matrix(0, nrow(set1), nrow(set2))
	parent <- unique(tree2$parent[tree2$parent!=tree2$child])
	set1 <- as.matrix(set1/set1[,1])
	set2 <- as.matrix(set2/set2[,1])
	for (i in 1:nrow(set1)){
		d1 <- abs(t(set2) - as.vector(set1[i,]))
		distmat[i,] <- colSums(d1)
	}
	return(distmat)
}


find_all_child <- function(parent, tree2){
	all_child <- parent
	tree2_uni <- tree2[tree2$parent!=tree2$child,]
	# iterate the child nodes
	k_ind <- 1
	k_stop <- TRUE
	while (k_stop){
		ind3 <- as.vector(tree2_uni$child[tree2_uni$parent == all_child[k_ind]])
		if (length(ind3)==0 & k_ind == length(all_child)) k_stop=FALSE
		all_child <- c(all_child, ind3)
		k_ind <- k_ind + 1
	}
	return(all_child)
}

find_all_parent <- function(chd, tree2){
	all_parent <- chd
	tree2_uni <- tree2[tree2$parent!=tree2$child,]
	# iterate the parent nodes
	current <- as.character(tree2_uni[tree2_uni$child==chd,1])
	while (length(current)>0){
		all_parent <- c(all_parent, current)
		current <- as.character(tree2_uni[tree2_uni$child==current,1])
	}
	return(all_parent)
}

find_child <- function(parent, tree2){
	return(as.character(tree2$child[tree2$parent==parent]))
}

find_parent <- function(chd, tree2){
	return(as.character(tree2$parent[tree2$child==chd]))
}

total2sub <- function(set_total, tree2){
	parent <- as.vector(unique(tree2$parent[tree2$parent!=tree2$child]))
	set_sub <- set_total
	for (i in 1:length(parent)){
		j=parent[i]
		chd <- find_child(j,tree2)
		subsets <- set_total[,chd]
        chd_j <- chd==j
		subsets[,chd_j] <- subsets[,chd_j] - rowSums(as.matrix(subsets[,!chd_j]))
		set_sub[,chd] <- subsets
	}
	return(set_sub)
}

sub2total <- function(set_sub, tree2){
	set_total <- set_sub
	name_taxa <- colnames(set_sub)
	n_taxa <- length(name_taxa)
	for (i in n_taxa:1){
		chd <- find_child(name_taxa[i],tree2)
		if(length(chd)>1){
			set_total[,name_taxa[i]] <- rowSums(set_total[,chd])
		}
	}
	return(set_total)
}


#### plotting functions ####

heatmap_cor <- function(pnode, set1, set2, tree2, elim=1){
    chd <- find_child(pnode, tree2)
    sample1 <- set1[,chd]
    sample2 <- set2[,chd]
    chd_pnode <- chd==pnode
    sample1[,chd_pnode] <- sample1[,chd_pnode] - rowSums(as.matrix(sample1[,!chd_pnode]))
    sample2[,chd_pnode] <- sample2[,chd_pnode] - rowSums(as.matrix(sample2[,!chd_pnode]))
    sample10 <- sample1; sample20 <- sample2
    # vars with small total counts are elimniated
    zeros <- colMeans(sample1)>elim & colMeans(sample2)>elim
    # obs with small total counts are eliminated
    ind1 <- rowMeans(as.matrix(sample1[,zeros]))>elim
    ind2 <- rowMeans(as.matrix(sample2[,zeros]))>elim
    k <- sum(zeros)
    ind12 <- ind1&ind2
    # check if there's enough degree of freedom
    if (sum(ind12)>=k & k>1){
        sample1 <- sample1[,zeros]
        sample2 <- sample2[,zeros]
        sample1 <- as.matrix(sample1[ind12,])
        sample2 <- as.matrix(sample2[ind12,])
        
        k <- ncol(sample1)
        k2 <- ncol(sample2)
        if (k!=k2) stop("The two sets have different number of columns!")
        n1 <- nrow(sample1)
        n2 <- nrow(sample2)
        if (n1!=n2) stop("The two sets are not paired!")
        N1 <- rowSums(sample1)
        N2 <- rowSums(sample2)
        NN1 <- sum(N1)
        NN2 <- sum(N2)
        eyes <- diag(rep(1,k))
        if(NN1==n1 | NN2==n2 | n1<k){
            stop('Sample is too sparse to test!')
        }else{
            Nc1 <- (NN1-sum(N1^2)/NN1)/(n1-1)
            Nc2 <- (NN2-sum(N2^2)/NN2)/(n2-1)

            m1 <- colSums(sample1)/NN1
            m2 <- colSums(sample2)/NN2
            mi1 <- sample1/N1
            mi2 <- sample2/N2
            diff1 <- t(mi1)-m1
            diff2 <- t(mi2)-m2
            S1 <- diff1%*%diag(N1)%*%t(diff1)/(n1-1)
            S2 <- diff2%*%diag(N2)%*%t(diff2)/(n2-1)
            G1 <- (diag(as.vector(N1%*%mi1)) - t(mi1)%*%diag(N1)%*%mi1) / (NN1-n1)
            G2 <- (diag(as.vector(N2%*%mi2)) - t(mi2)%*%diag(N2)%*%mi2) / (NN2-n2)
            V <- diff1%*%diag(N1+N2)%*%t(diff2) / (n1-1)
            V <- V / (Nc1+Nc2)*sum(N1*N2)/NN1/NN2
            Sig1 <- S1*sum(N1^2)/Nc1/NN1^2 + G1*(Nc1-sum(N1^2)/NN1)/NN1/Nc1
            Sig2 <- S2*sum(N2^2)/Nc2/NN2^2 + G2*(Nc2-sum(N2^2)/NN2)/NN2/Nc2
            Sig1 <- trunc_matrix(Sig1)
            Sig2 <- trunc_matrix(Sig2)

            Smat1 <- trunc_matrix(Sig1+Sig2-V-t(V))
            Smat2 <- trunc_matrix(Sig1+Sig2)
        }
    }else{stop('Sample is too sparse to test!')}
    dSmat1 <- sqrt(diag(Sig1))
    dSmat2 <- sqrt(diag(Sig2))
    hmat <- t(V/dSmat1)/dSmat2
    hmat <- hmat+t(hmat)
    colnames(hmat) <- rownames(hmat) <- colnames(sample1)
    melt_hmat <- melt(hmat)
    colnames(melt_hmat)[c(1,2)] <- c('Var1', 'Var2')
    myheatmap <- ggplot(data = melt_hmat, aes(x=Var1, y=Var2, fill=value)) + 
    ggtitle(pnode) +
    geom_tile(color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
        midpoint = 0, limit = c(-1,1), name="Differential Correlation") +
        theme_minimal() + # minimal theme
        theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
        coord_fixed()
    print(myheatmap)
    sample10 <- sample10[,colnames(sample1)]
    sample20 <- sample20[,colnames(sample2)]
    p1 <- sample10/rowSums(sample10)
    p2 <- sample20/rowSums(sample20)
    return(list(hmat=hmat, p1=p1, p2=p2))
}


barplot_nodes <- function(nodes, set1, set2, tree2, grp.names=c('Group 1', 'Group 2'), elim=1, common=0.05){
    NAs <- list()
    hset1 <- matrix(, nrow(set1), 0)
    hset2 <- matrix(, nrow(set2), 0)

    for (i in 1:length(nodes)){
        chd <- find_child(nodes[i],tree2)
        sample1 <- set1[,chd]
        sample2 <- set2[,chd]
        chd_j <- chd==nodes[i]
        sample1[,chd_j] <- sample1[,chd_j] - rowSums(as.matrix(sample1[,!chd_j]))
        sample2[,chd_j] <- sample2[,chd_j] - rowSums(as.matrix(sample2[,!chd_j]))
        # vars with small total counts are elimniated
        zeros <- colMeans(sample1)>elim & colMeans(sample2)>elim
        # obs with small total counts are eliminated
        ind1 <- rowMeans(as.matrix(sample1[,zeros]))>elim
        ind2 <- rowMeans(as.matrix(sample2[,zeros]))>elim
        #ind12 <- ind1&ind2
        sample1 <- sample1[,zeros]
        sample2 <- sample2[,zeros]
        #sample1 <- as.matrix(sample1[ind12,])
        #sample2 <- as.matrix(sample2[ind12,])

        hs1 <- sample1/rowSums(sample1)
        hs2 <- sample2/rowSums(sample2)
        mymedian <- function(x){median(x[!is.na(x)])}
        pchange <- apply(hs1,2,mymedian) - apply(hs2,2,mymedian)
        if(nrow(hs1)==nrow(hs2)) pchange <- apply(hs1-hs2,2,mymedian)
        ind <- order(-abs(pchange))
        hs1 <- hs1[,ind]
        hs2 <- hs2[,ind]
        colnames(hs1) <- colnames(hs2) <- paste(colnames(hs1),i)
        hset1 <- cbind(hset1,hs1)
        hset2 <- cbind(hset2,hs2)
        NAs[[i]] <- colSums(is.na(rbind(hs1,hs2)))
    }

    noNAnames <- NULL
    emptynames <- paste('E',1:6,sep='')
    for (i in 1:length(NAs)){
        # insert empty columns so that missing values are displayed as color grey
        noNAnames <- c(noNAnames, names(NAs[[i]]), emptynames)
    }
    
    n1 <- nrow(hset1)
    n2 <- nrow(hset2)

    emptyblock <- as.data.frame(matrix(0, n1+n2, 6))
    colnames(emptyblock) <- emptynames
    subset12 <- rbind(hset1,hset2)
    subset12 <- cbind(subset12, emptyblock)
    subset12 <- subset12[,noNAnames]
        
    # insert empty columns so that missing values are displayed as color grey
    nc <- ncol(subset12); excol <- nc%/%5
    subset120 <- as.data.frame(matrix(0,nrow(subset12),ncol(subset12)+excol))
    for (i in 1:excol){
        colnames(subset120)[(i*6-5):(i*6-1)] <- colnames(subset12)[(i*5-4):(i*5)]
        subset120[,(i*6-5):(i*6-1)] <- subset12[,(i*5-4):(i*5)]
    }
    if(ncol(subset12)>=i*5+1){
        colnames(subset120)[(i*6+1):ncol(subset120)] <- colnames(subset12)[(i*5+1):ncol(subset12)]
        subset120[,(i*6+1):ncol(subset120)] <- subset12[,(i*5+1):ncol(subset12)]
    }
    for (i in 1:length(NAs)){
        namesel <- names(NAs[[i]])
        sel <- is.na(rowSums(subset120[,namesel]))
        ind <- match(namesel[length(namesel)],colnames(subset120))
        subset120[sel,ind+(-ind)%%6] <- 1
    }
    subset120[is.na(subset120)] <- 0
    # add empty column between the two groups and add a label column
    subset120 <- rbind(subset120[1:n1,],subset120[1:4,]*NA,subset120[(n1+1):(n1+n2),], subset120[1:8,]*NA, colMeans(subset120,na.rm=T))
    # position of right labels
    rlab <- colMeans(subset120,na.rm=T)
    common_ones <- rlab >= common
    rlab_y <- rlab/2
    for (i in 2:length(rlab_y)){
        rlab_y[i] <- sum(rlab[1:i]) - rlab_y[i]
    }
    rlab_name <- colnames(subset120)
    cut.tail <- function(x){substr(x,1,nchar(x)-2)}
    rlab_name <- unlist(lapply(rlab_name,cut.tail))

    mycolor <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#f0f0f0')
    par(mar=c(3,8,0.5,0.5), lwd=0.3)
    barplot(t(subset120), col=mycolor, xaxt='n',yaxt='n', xlim=c(0,(n1+n2+14)*1.2+45), 
        ylim=c(0, length(nodes)+0.1), space=0.2)
    axis(2, 0:length(nodes), rep('', length(nodes)+1), cex.axis=0.8)
    axis(2, 1:length(nodes)-0.5, nodes, las=1, cex.axis=0.8, tick=FALSE) # leftside label
    axis(1, c(0, n1)*1.2, c('',''), cex.axis=0.8)
    axis(1, c(n1+4, n1+n2+4)*1.2, c('',''), cex.axis=0.8)
    axis(1, c(n1/2, n1+n2/2+2)*1.2, grp.names, cex.axis=0.8, tick=FALSE) # bottom label
    text((n1+n2+14)*1.2, rlab_y[common_ones], rlab_name[common_ones], pos=4,cex=0.6)
}
