# The R code is applied to cluster differentially expressed genes(DEGs) by their expression in all cell clusters
# need 1. factoextra package, please install.packages("factoextra");2. fields package, please install.package("fields")
# Three input files are required: 1. the reads count matrix file (reads_count_matrix_demo.txt), 2. the cell cluster annotation file (cell_cluster_annotation_demo.txt), and 3. DEGs list (DEG_list_demo.txt)
# Please prepare your files in the format as the demo files of "reads_count_matrix_demo.txt", "cell_cluster_annotation_demo.txt" and "DEG_list_demo.txt".

# Output£º1. the heatmap of cluster DEGs across cell clusters (generated in the plot window); 2. the gene name of genes in each group (generated "GeneRes.txt")


# Load the input files

ReadsCount <- read.table ("reads_count_matrix_demo.txt",sep="\t",header=FALSE);
CellType <- read.table("cell_cluster_annotation_demo.txt",sep="\t");
DEGs <- read.table("DEG_list_demo.txt",sep="\t")

#############################################

# Get expression percentage matrix for all genes in the cell clusters
a <- dim(ReadsCount);
exp <- as.matrix(ReadsCount[2:a[1],2:a[2]]);
exp <- apply(exp,2,as.numeric);
colnames(exp) <- ReadsCount[1,2:a[2]];
rownames(exp) <- ReadsCount[2:a[1],1];

unique_cluster <- unique(CellType[,2]);
a <- length(unique_cluster);
Exp_per <- c();
cluster_id <- c();
for (i in 1:a)
{
	print (unique_cluster[i]);
	m=which(CellType[,2]==unique_cluster[i]);
	CellId<- CellType[m,1];
	aa <- length(CellId);
	if (aa>20)
	{	
		CellIndex <- c();
		cluster_id <- c(cluster_id,unique_cluster[i]);
		for (j in 1:aa)
		{
			n=which(colnames(exp)==CellId[j])
			CellIndex <- c(CellIndex,n);		
		}
		exp_type <- exp[,CellIndex];
		exp_type <- exp_type>0;
		exp_type_sum <- apply(exp_type,1,sum);		
		exp_type_per <- exp_type_sum/length(CellIndex);		
		Exp_per <- cbind(Exp_per,exp_type_per);
	}
}
colnames(Exp_per) <- cluster_id;

###############################################
	
# Get expression percentage matrix for DEGs in the cell clusters
DEGs <- as.matrix(DEGs);
Gene <- rownames(Exp_per);

Index <- c();
a <- length(DEGs);
for (i in 1:a)
{
	m <- which(Gene==DEGs[i]);
	if (length(m)>0) {Index <- c(Index, m); }
}
GeneList <- Gene[Index];
Exp_per1 <- Exp_per[Index,];

#################################################################	
	
# Cluster DEGs

Tisu_per <- Exp_per1;
library(factoextra)
df <- Tisu_per;
result <- dist(df, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D2")
fviz_dend(result_hc, cex = 0.6)

#################################################################

# Get expression divergence for gene groups under different k parameters


k_list <- c();
for (k in 2:40)
{
	print(k)
	clusterCut <- cutree(result_hc, k)
	a=sort(unique(clusterCut))
	zzz=length(a);
	GeneRes=c();
	PerRes=c();
	kx=0
	for (i in a)
	{
		aa=GeneList[clusterCut==i];
		bb=Tisu_per[clusterCut==i,];
		if (length(aa)>4) # Only analyze the gene group with more than 4 genes
		{
			if (kx%%2==1) {zx=1} else {zx=0}
			bbb=cbind(bb,zx);
			cc=cbind(aa,i);
			GeneRes=rbind(GeneRes,cc);
			PerRes=rbind(PerRes,bbb)
			kx=kx+1;
		}
	}

	# Get expression dissimilarity for each cell cluster 
	Cluster <- unique(GeneRes[,2]);
	a <- length(Cluster);
	CTyeScore <- c();
	absdif <- c();
	for (i in 1:a)
	{
		m <- which(GeneRes[,2]==Cluster[i]);
		Ge <- GeneRes[m,1];
		n <- length(Ge);
		Index <- c();
		for (j in 1:n)
		{
			x <- which(GeneList==Ge[j]);
			Index <- c(Index,x);
		}
		Exp_c <- Exp_per1[Index,];
		Median_c <- apply(Exp_c,2,mean);
		z <- c(Median_c);
		CTyeScore <- cbind(CTyeScore,z);
		Max_c <- apply(Exp_c,2,max);
		Min_c <- apply(Exp_c,2,min);
		absdif <- c(absdif, mean(abs(Max_c-Min_c)));
	}

	yyy <- max(absdif) # get expression divergence for gene groups under the k parameter
	k_list <- rbind(k_list,c(k,yyy));
}

######################################################################

# get kdeter parameter
m=which(k_list[,2]<0.5)
if (length(m)>0) {mm <- min(m); kdeter <- k_list[mm,1];} else {kdeter <- 40;}

########################################################################

# Get gene groups with more than 4 genes and their expression heatmap under kdeter parameter.
	
clusterCut <- cutree(result_hc,kdeter)
a=sort(unique(clusterCut))
zzz=length(a);
GeneRes=c();
PerRes=c();
kx=0
for (i in a)
{
	aa=GeneList[clusterCut==i];
	bb=Tisu_per[clusterCut==i,];
	if (length(aa)>4)
	{
	if (kx%%2==1) {Group=1} else {Group=0}
	bbb=cbind(bb,Group);
	cc=cbind(aa,i);
	GeneRes=rbind(GeneRes,cc);
	PerRes=rbind(PerRes,bbb)
	kx=kx+1;
	}
}


library(fields)
image.plot(t(PerRes),col = hcl.colors(12, "YlOrRd", rev = TRUE))
colnames(PerRes)
write.table(GeneRes,"GeneRes.txt",sep="\t", row.name=FALSE);

###################################################################