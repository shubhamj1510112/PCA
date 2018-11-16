# set working directory
setwd("C:/Users/user/Desktop/elucidata")
# check for correct working directory
getwd()

# load the data
genedata1 <- read.csv("Assignment-gene_data.csv", sep=",", header=T)
metadata1 <- read.csv("Assignment-Meta_data sheet.csv", sep=",", header=T)

####### data pre-processing ###########
# remove identifier column
genedata2 = genedata1[,-1]
# check for missing values; 1 means column with missing value
print(table(apply(genedata2, 2, function(x) {sum(is.na(x))}))) # for columns
print(table(apply(genedata2, 1, function(x) {sum(is.na(x))}))) # for rows
# drop observations with missing values
genedata2 = genedata2[complete.cases(genedata2), ]
# transpose data so that samples become observations and gene intensities become variables
genedata3 = t(genedata2)
genedata3 = as.data.frame(genedata3)
genedata3 <- genedata3[-1, ]
genedata3$sIdx = rownames(genedata3)
# define new column names from gene_1 to gene_22410
name_col =   as.character(genedata2[1:22410,1])
name_col = c(name_col, 'sIdx')
colnames(genedata3) = name_col
# metadata file is clean so no preprocessing required
# merge data from intensity and metadata file
df = merge(x = genedata3, y = metadata1, all.x=T, by = "sIdx")
# now sIdx (sample ID) is the first column followed by the gene intensity values
# and time data for the samples is also added to df
# df has a total of 22413 columns and 30 rows
# convert gene intensity values to the numeric data format
for (i in 2:22411){
  df[,i] = as.numeric(as.character(df[,i]))
}
# identify and remove constant or not varying variables
df1 = df[,2:22411]
if (!require("caret")) install.packages("caret")
require(caret)
if (!require("foreach")) install.packages("foreach")
require(foreach)
nzv <- nearZeroVar(df1, freqCut=95/5, uniqueCut=10, allowParallel=T)
df2 <- df1[,-nzv]
df3 = cbind(df[,1], df2, df[,22412:22413])
new_var_count = ncol(df2) + 1

######### PCA analysis ###################
# pca analysis using prcomp function
gene_pca <- prcomp(df3[,2:new_var_count], scale = TRUE)
# construct PCA plots
if (!require("factoextra")) install.packages("factoextra")
require(factoextra)
png(filename = "pca_plot_without_centroid.png", width = 7, height = 5, units='in', res=300, bg = "white")
print({p1 = fviz_pca_ind(gene_pca, habillage=df$Time,pointsize=3,pointshape = 19,
                         label = "none", legend.title = "Time", 
                         title="PCA Plot Without centroid")})
dev.off()
png(filename = "pca_plot_with_centroid.png", width = 7, height = 5, units='in', res=300, bg = "white")
print({p1 = fviz_pca_ind(gene_pca, habillage=df$Time,pointsize=3,pointshape = 19,
                         label = "none", legend.title = "Time", 
                         title="PCA Plot With centroid")})
dev.off()
# construct scree plot
png(filename = "scree_plot.png", width = 7, height = 5, units='in', res=300, bg = "white")
print({p1 = fviz_eig(gene_pca, addlabels = TRUE, ylim = c(0, 50))})
dev.off()
# construct cumulative scree plot
if (!require("ggfortify")) install.packages("ggfortify")
require(ggfortify)
std_dev <- gene_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
png(filename = "cumulative_scree_plot.png", width = 5, height = 5, units = "in", res=300, bg = "white")
plot(cumsum(prop_varex), xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", type = "b", col="brown4",cex = 1.5, pch=16)
dev.off()

# identify top contributing genes for PC1-5
if (!require("FactoMineR")) install.packages("FactoMineR")
require(FactoMineR)
gene_pca <- PCA(df3[,2:new_var_count], graph = FALSE)
cont = as.data.frame(gene_pca$var$contrib)
cont1 = as.data.frame(cont[order(cont[,1]),])
pc1 = rownames(cont1)
cont2 = as.data.frame(cont[order(cont[,2]),])
pc2 = rownames(cont2)
cont3 = as.data.frame(cont[order(cont[,3]),])
pc3 = rownames(cont3)
cont4 = as.data.frame(cont[order(cont[,4]),])
pc4 = rownames(cont4)
contrib = as.data.frame(cbind(PC1=pc1, PC2=pc2, PC3=pc3, PC4=pc4))
write.table(contrib, file="contribution_gene_to_PCs.csv", row.names=F, col.names=T, sep="\t", quote=F) # save as csv 

# identify quality of representation of samples in PC-1 and PC-2
# for PC1
png(filename = "representation_quality_PC1.png", width = 6, height = 4, units = "in", res=300, bg = "white")
fviz_cos2(gene_pca, choice = "ind", axes=1, title="Cos2 of Samples for PC-1")
dev.off()
# for PC2
png(filename = "representation_quality_PC2.png", width = 6, height = 4, units = "in", res=300, bg = "white")
fviz_cos2(gene_pca, choice = "ind", axes=2, title="Cos2 of Samples for PC-2")
dev.off()

