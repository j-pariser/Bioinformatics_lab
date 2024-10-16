library(Seurat)
library(SeuratDisk)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)
library(CellChat)
library(patchwork)
library(graphics)
library(Matrix)
options(stringsAsFactors = FALSE)

library(circlize)
library(colorspace)
options(repr.plot.width = 12, repr.plot.height = 9, repr.plot.res = 300)

library(pracma)
library(glue)
library(anndata)


adata = read_h5ad('/extra/zhanglab0/CommonData/Multiome/AUD/AUD.h5ad')
filename = "./cellchat_AUD_computed.rds"


cellchat <- createCellChat(object = t(as(adata$X, 'CsparseMatrix')), meta = adata$obs, group.by = "NSForest")


condition <- "AUD"
type <- "triMean"        #  c("triMean", "truncatedMean", "thresholdedMean", "median")
trim <- 0.1              #c(0.05, 0.10, 0.15, 0.20, 0.25)

#choose database
CellChatDB <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
#set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

#do parallel
library(future)
options(future.globals.maxSize = +Inf)
future::plan("multicore", workers = 1)

#analysis
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


######################################
## KEY LINE ##########################
cellchat <- computeCommunProb(cellchat, type = type, trim = trim)
######################################


#Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = filename)


