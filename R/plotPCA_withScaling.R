#' Construct an interactive PCA plot (without data scaling) using the gene intensity data from different samples or individuals and also save the plot as a static png file in the working directory
#'
#' This function makes an interactive PCA plot using the gene intensity data for different samples or individuals and also save the plot as a static png file in the working directory. The PCA analysis is done after scaling the data.
#' @param genedata Gene_intensity_data metadata the_sample_metadata.
#' @return The interactive PCA plot with scaling
#' @export
#' @examples
#' @POST
#' plotpca(genename="Assignment-gene_data.csv", metadata="Assignment-Meta_data sheet.csv")

plotpca_withscaling <- function(genedata, metadata){
  # data loading and error handling
  genedata1 <- if(is.character(genedata) && file.exists(genedata)){
    read.csv(genedata, sep=",", header=T)
  } else {
    stop ('intensity file does not exist')
  }
  metadata1 <- if(is.character(metadata) && file.exists(metadata)){
    read.csv(metadata, sep=",", header=T)
  } else {
    stop ('metadata file does not exist')
  }
  if (ncol(genedata1)==1) stop ('intensity file is not comma separated')
  if (ncol(metadata1)==1) stop ('metadata file is not comma separated')
  if (length(unique(metadata1['Unit']))!=1) stop ('time-point units are different in metadata file')
  if (tolower(colnames(genedata1[2])) != 'symbol') stop ('gene symbol column is missing in intensity file')
  # data preprocessing for intensity file
  # remove the identifier column
  genedata2 = genedata1[,-1]
  # check for missing values; 1 means column with missing value
  print(table(apply(genedata2, 2, function(x) {sum(is.na(x))})))
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
  # A total of 22413 columns and 30 rows
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
  # construct the pca plot
  gene_pca <- prcomp(df3[,2:new_var_count], scale = TRUE)
  if (!require("factoextra")) install.packages("factoextra")
  require(factoextra)
  png(filename = "pca_plot_withscaling.png", width = 10, height = 5, units='in', res=300, bg = "white")
  print({p1 = fviz_pca_biplot(gene_pca, habillage=df$Time,pointsize=3,pointshape = 19,
                         label = "var", legend.title = "Time", invisible="quali", 
                         title="PCA Plot With Scaling", col.var="black", select.var=list(contrib = 4), repel=T)})
  dev.off()
  if (!require("plotly")) install.packages("plotly")
  require(plotly)
  p1 = fviz_pca_biplot(gene_pca, habillage=df$Time,pointsize=3,pointshape = 19,
                         label = "var", legend.title = "Time", invisible="quali", 
                         title="PCA Plot With Scaling", col.var="black", select.var=list(contrib = 4), repel=T)
  return(ggplotly(p1, tooltip = c("x", "y", "colour")))
}
