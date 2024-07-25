library("optparse")
library("Seurat")
library("SingleCellExperiment")
library("scDblFinder")
library("data.table")

option_list = list(
    make_option(c("-d", "--data"), type="character", default=".", help="Path to the parent directory containing 'filtered_feature_bc_matrix' dir", metavar="character"),
    make_option(c("-s", "--sample"), type="character", default=NULL, help="Sample ID", metavar="character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$data)){
    print_help(opt_parser)
    stop("Input data directory is mandatory", call.=FALSE)
}

if (is.null(opt$sample)){
    print_help(opt_parser)
    stop("Input sample ID is mandatory", call.=FALSE)
}

eislet <- Read10X(paste(opt$data, '/filtered_feature_bc_matrix', sep=''))
eislet = as.data.frame(eislet)
dim(eislet)

eislet = eislet[apply(eislet, 1, function(x) !all(x==0)),]

counts.fname = as.matrix(eislet)

sce <- SingleCellExperiment(
list(counts=as(counts.fname, "dgCMatrix")), 
rowData=rownames(counts.fname), 
colData=DataFrame(Barcode=colnames(counts.fname))
)

rowData(sce)$Symbol = rownames(counts.fname)
colnames(sce) <- sce$Barcode
sce

#-------------------- Doublet detection by simulation
set.seed(100)

dbl.dens <- computeDoubletDensity(sce)#doubletCells
summary(dbl.dens)
sce$DoubletScore <- dbl.dens

table(dbl.dens < 2)
sce = sce[,dbl.dens < 2]
range(sce$DoubletScore)

#---------------------------------- Save results
eislet = eislet[,which(colnames(eislet) %in% colnames((assay(sce))))]

library(tibble)
eislet = add_column(eislet, GeneSymbol = rownames(counts.fname), .before = 1)

eislet[1:5,1:5]
dim(eislet)

fwrite(eislet, paste(opt$sample, '_Gene_Expression_Matrix_passed_QC_noDoublet.txt', sep = ''), quote = F, row.names = F, col.names = T, sep = '\t')
write.table(colnames((assay(sce))), paste(opt$sample, '_Cells_passed_QC_noDoublet.txt', sep = ''), quote = F, row.names = F, col.names = T, sep = '\t')
