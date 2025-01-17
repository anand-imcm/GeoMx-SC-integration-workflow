library("optparse")
library("AUCell")
library("NMF")
library("data.table")
library("Seurat")
library("dplyr")
library("patchwork")


option_list = list(
    make_option(c("-o", "--outdir"), type="character", default=".", help="Path to the output directory", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default=NULL, help="Sample ID", metavar="character"),
    make_option(c("-g", "--geneset"), type="character", default=NULL, help="Gene set file (tsv)", metavar="character"),
    make_option(c("-m", "--metadata"), type="character", default=NULL, help="Sample metadata file (tsv)", metavar="character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$outdir)){
    print_help(opt_parser)
    stop("Path to the output directory is mandatory", call.=FALSE)
}

if (is.null(opt$sample)){
    print_help(opt_parser)
    stop("Input sample ID is mandatory", call.=FALSE)
}

if (is.null(opt$geneset)){
    print_help(opt_parser)
    stop("Input gene set tsv file is mandatory", call.=FALSE)
}

if (is.null(opt$metadata)){
    print_help(opt_parser)
    stop("Input sample metadata tsv file is mandatory", call.=FALSE)
}

geneSets <- fread(opt$geneset) %>% as.list() %>% lapply(function(x) x[x != "" & !is.na(x)])

setwd(opt$outdir)

PI=paste0(opt$sample,'_CellType_Enrichment')
output.dir_perSample <- opt$outdir

MetaData = fread(opt$metadata, stringsAsFactors = F, header = T)
MetaData = as.data.frame(MetaData)

ListOfFolders_All = list.dirs(full.names = F, recursive = T)

toMatch <- 'QC' # MetaData[which(MetaData$Diagnosis == d & MetaData$Region == r),'Sample_ID']
ListOfFolders <- unique(grep(paste(toMatch,collapse="|"), ListOfFolders_All, value=TRUE))
i=1
for(i in 1:length(ListOfFolders))
{   
    Sample <- fread(paste0(ListOfFolders[i], '/', 'Gene_Expression_Matrix_passed_QC_noDoublet.txt'), stringsAsFactors = F, header = T)
    Sample = as.data.frame(Sample)
    row.names(Sample) = Sample$GeneSymbol
    Sample = Sample[,-which(colnames(Sample) == 'GeneSymbol')]

    # You should adjust this line per data set
    colnames(Sample) = gsub('-.*', paste0('-', gsub('\\/.*', '', ListOfFolders[i])), colnames(Sample))

    # ------------ filtering
    # Initialize the Seurat object with the raw (non-normalized data).
    pbmc <- CreateSeuratObject(counts = Sample, project = "pbmc3k", min.cells = 3, min.features = 200)

    # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc <- subset(pbmc, subset = percent.mt < 10)

    Sample = Sample[which(row.names(Sample) %in% row.names(pbmc)), which(colnames(Sample) %in% colnames(pbmc))]

    if(i==1){TP_profile = Sample}else{TP_profile = merge(TP_profile, Sample, by = "row.names"); row.names(TP_profile) = TP_profile$Row.names; TP_profile = TP_profile[,-which(colnames(TP_profile)[1] == "Row.names")]}
    rm(Sample)
    print(i)
    print(dim(TP_profile))
}

dim(TP_profile)
# TP_profile = TP_profile[apply(TP_profile, 1, function(x) !all(x==0)),]
# dim(TP_profile)

TP_profile[1:4,1:4]
fwrite(cbind(GeneSymbol = row.names(TP_profile), TP_profile), paste0(output.dir_perSample,'/', gsub('_CellType_Enrichment', '_Gene_Expression_All_Samples.txt', PI)), quote = F, row.names = F, sep = '\t')
print(paste0("Complete: ", opt$outdir,":", opt$sample))


set.seed(123)

TP_profile <- as.matrix(TP_profile)
cells_rankings <- AUCell_buildRankings(TP_profile)

# ---------------------------------- reports Genes in the gene sets NOT available in the dataset
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, nCores=1)

# ---------------------------------- cells assignment
selectedThresholds <- NULL

set.seed(123)
par(mfrow=c(3,3)) # PROVIDE ENOUGH SPACE IN PLOTS ENVIRONMENT OF THE RSTUDIO
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)

## ------ export cell Assignment as text
cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"

# ---------------------------------- save cell assignments as a heatmap
assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
#assignmentMat = assignmentMat[!duplicated(rowSums(assignmentMat)),]
dim(assignmentMat)

set.seed(123)
setwd(paste0(output.dir_perSample,'/'))
aheatmap(assignmentMat, scale="none", color="-RdYlBu2:100", legend=FALSE, filename = paste0('Celltype_assignment_HeatMap.pdf'), height = 10, width = 8)
getwd()

#----- make a summary

Freq_geneSet <- data.frame(table(assignmentTable[,"geneSet"]), stringsAsFactors = F)
Freq_geneSet$Var1 = as.character(Freq_geneSet$Var1)
colnames(Freq_geneSet)[1] = 'CellType'

Freq_geneSet$Porportion = Freq_geneSet$Freq/sum(Freq_geneSet$Freq)
sum(Freq_geneSet$Porportion)

options(digits=2)
colnames(Freq_geneSet)
Freq_geneSet$CellType = gsub('_',' ',Freq_geneSet$CellType)

#-----  Removing duplicate along with the original value
keep_singles <- function(v){ v[!(v %in% v[duplicated(v)])] }
assignmentTable_dedup = assignmentTable[which(assignmentTable$cell %in% keep_singles(assignmentTable$cell)),]

table(assignmentTable_dedup$geneSet)
table(assignmentTable_dedup$geneSet)/dim(assignmentTable_dedup)[1]

dim(assignmentTable_dedup)
dim(assignmentTable)

#----- save the results
fwrite(Freq_geneSet, paste0(output.dir_perSample, '/Table_Summary.txt'), quote = F, row.names = F, sep = '\t')

assignmentTable$Sample_ID = gsub('.*-', '', assignmentTable$cell)
assignmentTable = merge(assignmentTable, MetaData, by = 'Sample_ID')
assignmentTable$cell_orig = gsub('-.*', '-1', assignmentTable$cell)
fwrite(assignmentTable, paste0(output.dir_perSample, '/assignmentTable.txt'), quote = F, row.names = F, sep = '\t')

assignmentTable_dedup$Sample_ID = gsub('.*-', '', assignmentTable_dedup$cell)
assignmentTable_dedup = merge(assignmentTable_dedup, MetaData, by = 'Sample_ID')
assignmentTable_dedup$cell_orig = gsub('-.*', '-1', assignmentTable_dedup$cell)
fwrite(assignmentTable_dedup, paste0(output.dir_perSample, '/assignmentTable_dedup.txt'), quote = F, row.names = F, sep = '\t')

#====================================================================================
setwd(paste0(output.dir_perSample,'/'))

CellType = unique(assignmentTable_dedup$geneSet)
Regions = unique(assignmentTable_dedup$Region)
Sex = unique(assignmentTable_dedup$Sex)
Diagnosis = unique(assignmentTable_dedup$Diagnosis)


i=1
for(i in 1:length(CellType))
{
  for(j in 1:length(Regions))
  {
    for(d in 1:length(Diagnosis))
    {
      
      Temp <- assignmentTable_dedup[which(assignmentTable_dedup$geneSet == CellType[i] & assignmentTable_dedup$Region == Regions[j] & assignmentTable_dedup$Diagnosis == Diagnosis[d]),]
      
      print(paste(CellType[i], Regions[j], Diagnosis[d], sep = '-'))
      print(paste0('Dim: ', dim(Temp)[1]))
      
      if(dim(Temp)[1] > 10)
      {
        fwrite(Temp, paste0(output.dir_perSample,'/', CellType[i], '_', Regions[j], '_', Diagnosis[d], '_assignmentTable_dedup.txt'), quote = F, row.names = F, sep = '\t')
      }
      for(s in 1:length(Sex))
      {
        TempS <- Temp[which(Temp$Sex == Sex[s]),]
        print(paste0('Sex: ', dim(TempS)[1]))
        
        if(dim(TempS)[1] > 10)
        {
          fwrite(TempS, paste0(output.dir_perSample,'/', CellType[i], '_', Regions[j], '_', Diagnosis[d], '_', Sex[s], '_assignmentTable_dedup.txt'), quote = F, row.names = F, sep = '\t')
        }
      }
      
    }
  }
}


# ---------------------------------- visualization as tSNE plot
temp = TP_profile[,which(colnames(TP_profile) %in% assignmentTable_dedup$cell)]
temp = temp[apply(temp, 1, function(x) !all(x==0)),]
dim(temp)

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = temp, project = "pbmc3k", min.cells = 3, min.features = 10)
pbmc
dim(TP_profile)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = percent.mt < 5)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.0082)  # I changed the resolution to get 6 clusters
table(pbmc$seurat_clusters)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
#pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)


# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
# DimPlot(pbmc, reduction = "umap")
# DimPlot(pbmc, reduction = "tsne")

#===============================

cellsTsne = as.matrix(pbmc@reductions$tsne@cell.embeddings)
head(cellsTsne)

dim((cellsTsne))
table(is.na(match(row.names(cellsTsne), assignmentTable_dedup$cell)))

assignmentTable_dedup = assignmentTable_dedup[match(row.names(cellsTsne), assignmentTable_dedup$cell),]

identical(assignmentTable_dedup$cell, row.names(cellsTsne))

# ---------------------------------- select 6 top dominant enriched cell types
selectedThresholds <- getThresholdSelected(cells_assignment)
selectedCellTypes = unique(assignmentTable_dedup$geneSet)
selectedThresholds = selectedThresholds[which(names(selectedThresholds) %in% selectedCellTypes)]

# ---------------------------------- visualization of tSNE plot + functional annotation
setwd(paste0(output.dir_perSample,'/'))
set.seed(123)
pdf(paste0('tSNE_Plot.pdf'), height = 8, width = 10, useDingbats = F)
par(mfrow=c(2,3)) # Splits the plot into two rows and three columns

geneSetName="Oligodendrocytes"
for(geneSetName in selectedCellTypes)
{
  nBreaks <- 5 # Number of levels in the color palettes
  # Color palette for the cells that do not pass the threshold
  colorPal_Neg <- grDevices::colorRampPalette(c("black","blue", "skyblue"))(nBreaks)
  # Color palette for the cells that pass the threshold
  colorPal_Pos <- grDevices::colorRampPalette(c("pink", "magenta", "red"))(nBreaks)
  
  passThreshold <- getAUC(cells_AUC)[geneSetName,] >  selectedThresholds[geneSetName]
  
  if(sum(passThreshold) > 0 )
  {
aucSplit <- split(getAUC(cells_AUC)[geneSetName,], passThreshold)

# Assign cell color
cellColor <- c(setNames(colorPal_Neg[cut(aucSplit[[1]], breaks=nBreaks)], names(aucSplit[[1]])),
   setNames(colorPal_Pos[cut(aucSplit[[2]], breaks=nBreaks)], names(aucSplit[[2]])))

plot(cellsTsne, main=geneSetName,
 sub="Pink/red cells pass the threshold",
 col=cellColor[rownames(cellsTsne)], pch=16, cex = 0.5)
  }
}
dev.off()