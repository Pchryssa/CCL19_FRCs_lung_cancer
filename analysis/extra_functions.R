#Functions
library(tibble)
library(tidyr)
library(circlize)
library(magrittr)

Remove_ensebl_id <- function(data_conv){
  data_conv <- as.SingleCellExperiment(data_conv, assay = "RNA")
  
  gene_name_fixdf <- data.frame(Gene_Name = rownames(data_conv))
  gene_name_fixdf$Symbol <- gsub("^[^.]*.","", gene_name_fixdf$Gene_Name)
  
  # Check again for double duplicates
  duplidx <- rownames((gene_name_fixdf[duplicated(gene_name_fixdf$Symbol), ]))
  
  while (length(duplidx)>0){
    Fixed_dupl <- gene_name_fixdf %>% 
      filter(row_number() %in% duplidx)  %>% 
      mutate(Symbol = paste0(Symbol, ".4"))
    gene_name_fixdf[duplidx,"Symbol"] <- Fixed_dupl$Symbol
    rownames(data_conv) <- gene_name_fixdf$Symbol
    duplidx <- rownames((gene_name_fixdf[duplicated(gene_name_fixdf$Symbol), ]))
  }
  
  # Convert SingleCellexperiment object back to seurat object by retaining the reductions, normalization, etc
  data_conv <- as.Seurat(data_conv, assay = NULL)
  return(data_conv)
}

get_full_gene_name <- function(gene, obj){
  return(grep(paste0(".",gene,"$"), rownames(obj), value = TRUE))
}

preprocessing <- function(seurat_obj, resolution){
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj, npcs = 30, verbose = FALSE,seed.use = 8734)
  seurat_obj <- RunTSNE(object = seurat_obj, reduction = "pca", dims = 1:20, seed.use = 8734)
  seurat_obj <- RunUMAP(object = seurat_obj, reduction = "pca", dims = 1:20, seed.use = 8734)
  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "pca", dims = 1:20, seed.use = 8734)
  for(k in 1:length(resolution)){
    seurat_obj <- FindClusters(object = seurat_obj, resolution = resolution[k], random.seed = 8734)
  }
  return(seurat_obj)
}


Cellchat_Analysis <- function(data){
  data$annot<-data$cell_type
  data.input<-GetAssayData(data,assay = "RNA",slot = "data")
  names_init <- c(rownames(data))
  names_split <- sub("^[^.]*.", "", names_init)
  names_split <- unname(names_split)
  rownames(data.input) <- names_split
  labels <- Idents(data)
  meta <- data.frame(group = data$annot, row.names = names(labels))
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  return(cellchat)
}

CellChatDownstreamAnalysis <- function(cellchat,type){
  # Cell-Cell Communication Analysis
  cellchat <- CellToCellCommunicationAnalysis(cellchat,type)
  return(cellchat)
}

CellToCellCommunicationAnalysis <- function(cellchat,type){
  # Preprocessing of Expression data
  cellchat <- LoadCellChatDB(cellchat,type)
  
  # Set the ligand-receptor interaction database
  cellchat <- DataPreprocessing(cellchat)
  
  # Inference of cell-cell communication network
  cellchat <- CellToCellCommunicationInference(cellchat)
  
  return(cellchat)
}

LoadCellChatDB <- function(cellchat,type){
  if(type == "human"){
    CellChatDB <- CellChatDB.human
  }else{
    CellChatDB <- CellChatDB.mouse
  }
  showDatabaseCategory(CellChatDB)
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  return(cellchat)
}

DataPreprocessing <- function(cellchat){
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  return(cellchat)
}

CellToCellCommunicationInference <- function(cellchat){
  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat,population.size = FALSE)
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  # Compute network centrality scores - slot 'netP' : inferred intercellular communication network of signaling pathways
  cellchat <- netAnalysis_computeCentrality(cellchat)
  return(cellchat)
}

# dotPlot_construct is based on netAnalysis_dot function from Cellchat R package
dotPlot_construct <-function(object, slot.name = "netP", pattern = c("outgoing", 
                                                                     "incoming"), cutoff = NULL, color.use = palet, pathway.show = NULL, 
                             group.show = NULL ) {
  pattern <- match.arg(pattern)
  patternSignaling <- methods::slot(object, slot.name)$pattern[[pattern]]
  data1 = patternSignaling$pattern$cell
  data2 = patternSignaling$pattern$signaling
  data = patternSignaling$data
  if (is.null(color.use)) {
    color.use <- scPalette(nlevels(data1$CellGroup))
  }
  if (is.null(cutoff)) {
    cutoff <- 1/length(unique(data1$Pattern))
  }
  options(warn = -1)
  data1$Contribution[data1$Contribution < cutoff] <- 0
  data2$Contribution[data2$Contribution < cutoff] <- 0
  data3 = merge(data1, data2, by.x = "Pattern", by.y = "Pattern")
  data3$Contribution <- data3$Contribution.x * data3$Contribution.y
  data3 <- data3[, colnames(data3) %in% c("CellGroup", "Signaling", 
                                          "Contribution")]
  if (!is.null(pathway.show)) {
    data3 <- data3[data3$Signaling %in% pathway.show, ]
    pathway.add <- pathway.show[which(pathway.show %in% data3$Signaling == 
                                        0)]
    if (length(pathway.add) > 1) {
      data.add <- expand.grid(CellGroup = levels(data1$CellGroup), 
                              Signaling = pathway.add)
      data.add$Contribution <- 0
      data3 <- rbind(data3, data.add)
    }
    data3$Signaling <- factor(data3$Signaling, levels = pathway.show)
  }
  if (!is.null(group.show)) {
    data3$CellGroup <- as.character(data3$CellGroup)
    data3 <- data3[data3$CellGroup %in% group.show, ]
    data3$CellGroup <- factor(data3$CellGroup, levels = group.show)
  }
  data <- as.data.frame(as.table(data))
  data <- data[data[, 3] != 0, ]
  data12 <- paste0(data[, 1], data[, 2])
  data312 <- paste0(data3[, 1], data3[, 2])
  idx1 <- which(match(data312, data12, nomatch = 0) == 0)
  data3$Contribution[idx1] <- 0
  data3$id <- data312
  data3 <- data3 %>% group_by(id) %>% top_n(1, Contribution)
  data3$Contribution[which(data3$Contribution == 0)] <- NA
  if (pattern == "outgoing") {
    data3$cond <-rep("out",nrow(data3))
  }
  else if (pattern == "incoming") {
    data3$cond <-rep("in",nrow(data3))
  }
  data3$celltype_cond <-paste(data3$CellGroup,"_",data3$cond)
  return(data3)
}

netAnalysis_joint_dot <-function (object, slot.name = "netP", cutoff = NULL, color.use = palet, pathway.show = NULL, 
                                  group.show = NULL, shape = 21, dot.size = c(1, 3), dot.alpha = 1, 
                                  main.title = NULL, font.size = 10, font.size.title = 10, 
                                  pathways = pathways,order_list = order_list, exclude =NULL) 
{
  
  if (!is.null(exclude)) {
    order_list  <- order_list[!(order_list %in% exclude)]
  }
  palet <- palet[names(palet) %in% order_list]
  df1 <-dotPlot_construct(cellchat, pattern = "outgoing",color.use = palet)
  df2 <-dotPlot_construct(cellchat, pattern = "incoming",color.use = palet)
  df <-rbind(df2,df1)
  df <-df[df$Signaling %in% pathways,]
  if (!is.null(exclude)) {
    df <-df[df$CellGroup %in% order_list,]
  }

  y_values_order <-c(rbind(paste0(order_list, " _ ", "out"),paste0(order_list, " _ ", "in")))
  
  df <-df[order(match(df$CellGroup, order_list)), ]
  df$Signaling <- factor(df$Signaling, levels = pathways)
  df$celltype_cond <- factor(df$celltype_cond, levels = y_values_order)
  df$CellGroup <- droplevels(df$CellGroup)
  gg <- ggplot(data = df, aes(x = Signaling, y = celltype_cond)) + 
    geom_point(aes(size = Contribution, fill = CellGroup, 
                   colour = CellGroup), shape = shape) + scale_size_continuous(range = dot.size) + 
    theme_linedraw() + scale_x_discrete(position = "bottom") + 
    ggtitle("In- and outcoming communication patterns") + theme(plot.title = element_text(hjust = 0.5)) + 
    theme(text = element_text(size = font.size), plot.title = element_text(size = font.size.title, 
                                                                           face = "plain"), axis.text.x = element_text(angle = 45, 
                                                                                                                       hjust = 1), axis.text.y = element_text(angle = 0, 
                                                                                                                                                              hjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    theme(axis.line.x = element_line(size = 0.25), axis.line.y = element_line(size = 0.25)) + 
    theme(panel.grid.major = element_line(colour = "grey90", 
                                          size = (0.1)),
          axis.ticks.length=unit(.05, "cm")) + coord_fixed(ratio=0.8)
  
  gg <- gg + scale_fill_manual(values = ggplot2::alpha(color.use, 
                                                       alpha = dot.alpha), drop = FALSE, na.value = "white")
  gg <- gg + scale_colour_manual(values = color.use, drop = FALSE, 
                                 na.value = "white")
  gg <- gg + theme(legend.title = element_text(size = 8), 
                   legend.text = element_text(size = 8))
  gg <- gg + scale_shape_manual(values= c("25","50","70"))
  gg <- gg + scale_y_discrete(label = function(x) {
    x %>% gsub(".*_", " ", .) 
  })
  gg
  
  return(gg)
}

Heatmap<- function(data, Selgenes, order_list, palet = palet){
  genes <- data.frame(gene=rownames(data)) %>%
    mutate(geneID=gsub("^.*?\\.", "", gene))
  
  Selgenes_df <- filter(genes, geneID %in% Selgenes)
  # order by Selgenes
  colnames(Selgenes_df) <- c('geneID','gene')
  Selgenes_df <-Selgenes_df[match(Selgenes, Selgenes_df$gene),]
  
  Idents(data) <- data$cell_type
  
  palet <- palet[names(palet) %in% data$cell_type] 
  palet<-palet[order(match(names(palet),data$cell_type))]
  levels(data)<-levels(data)[order(match(levels(data),order_list))]
  pOut <- avgHeatmap2(seurat = data, selGenes = Selgenes_df,
                      colVecIdent = palet,
                      ordVec=levels(data),
                      gapVecR=NULL, gapVecC=NULL,cc=T,
                      cr=F, condCol=F,
                      order_list= order_list)
}

# Average gene expression heatmap 
avgHeatmap2 <- function(seurat, selGenes, colVecIdent, colVecCond=NULL,
                        ordVec=NULL, gapVecR=NULL, gapVecC=NULL,cc=FALSE,
                        cr=FALSE, condCol=FALSE,
                        order_list = order_list){
  
  selGenes <- selGenes$gene
  ## assay data
  clusterAssigned <- as.data.frame(Idents(seurat)) %>%
    dplyr::mutate(cell=rownames(.))
  colnames(clusterAssigned)[1] <- "ident"
  seuratDat <- GetAssayData(seurat)
  
  ## genes of interest
  genes <- data.frame(gene=rownames(seurat)) %>% 
    mutate(geneID=gsub("^.*\\.", "", gene)) %>% filter(geneID %in% selGenes)
 
  ## matrix with averaged cnts per ident
  logNormExpres <- as.data.frame(t(as.matrix(
    seuratDat[which(rownames(seuratDat) %in% genes$gene),])))
  
  logNormExpres <- logNormExpres %>% dplyr::mutate(cell=rownames(.)) %>%
    dplyr::left_join(.,clusterAssigned, by=c("cell")) %>%
    dplyr::select(-cell) %>% dplyr::group_by(ident) %>%
    dplyr::summarise_all(mean)
  
  
  logNormExpresMa <- logNormExpres %>% dplyr::select(-ident) %>% as.matrix()
  rownames(logNormExpresMa) <- logNormExpres$ident
  logNormExpresMa <- t(logNormExpresMa)
  rownames(logNormExpresMa) <- gsub("^.*?\\.","",rownames(logNormExpresMa))
  
  ## remove genes if they are all the same in all groups
  ind <- apply(logNormExpresMa, 1, sd) == 0
  logNormExpresMa <- logNormExpresMa[!ind,]
  
  genes <- genes[!ind,]
  
  ## color columns according to cluster
  annotation_col <- as.data.frame(gsub("(^.*?_)","",
                                       colnames(logNormExpresMa)))%>%
    dplyr::mutate(celltype=gsub("(_.*$)","",colnames(logNormExpresMa)))
  colnames(annotation_col)[1] <- "col1"
  annotation_col <- annotation_col %>%
    dplyr::mutate(cond = gsub("(^[0-9]_?)","",col1)) %>%
    dplyr::select(cond, celltype)
  
  annotation_col$final <- NULL
  for (x in 1:nrow(annotation_col)){
    if (annotation_col$cond[x] != annotation_col$celltype[x]){
      annotation_col$final[x] <- paste(annotation_col$celltype[x],annotation_col$cond[x],sep="_")
    }else{
      annotation_col$final[x]<-annotation_col$celltype[x]
    }
  }
  rownames(annotation_col) <- colnames(logNormExpresMa) 
  ann_colors = list(
    cond = colVecCond,
    celltype=colVecIdent)
  if(is.null(ann_colors$cond)){
    annotation_col$cond <- NULL
  }
  ## adjust order
  logNormExpresMa <- logNormExpresMa[selGenes,]
  if(is.null(ordVec)){
    ordVec <- levels(seurat)
  }
  logNormExpresMa <- logNormExpresMa[,ordVec]
  
  annotation_col$final <-NULL
  p <- pheatmap(logNormExpresMa, scale="row" ,treeheight_row = 0, cluster_rows = cr, treeheight_col = 0,
                cluster_cols = F,
                color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(50),
                annotation_col = annotation_col, cellwidth=15, cellheight=10,
                annotation_colors = ann_colors, gaps_row = gapVecR, gaps_col = gapVecC)
  return(p)
}

# Gaussian kermel
GK <- function(x){return(1/sqrt(2*pi)*(exp(-0.5*x*x)))}



Compute_average_mean_expression <- function(object,genes_list,slot_type,gene_signature){
  sig_expr<- FetchData(object, vars = genes_list, slot = slot_type)
  
  # Get mean expression of selected genes per cell
  mean_expr <- rowMeans(x = sig_expr, na.rm = TRUE)
  
  col_nam <- paste0("gene.set.score_",gene_signature,"_", slot_type)
  print(col_nam)
  
  # Add mean expression values in 'object@meta.data$gene.set.score'
  if (all(names(x = mean_expr) == rownames(x = object@meta.data))) {
    cat("Cell names order match in 'mean_expr' and 'object@meta.data':\n", 
        "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
    object@meta.data[col_nam] <- mean_expr
  }
  return(c(object, col_nam))
  
}

# This function is a modification on the plotSmoothers function from R package tradeSeq
plotSmoothers1<- function(models, counts, gene, nPoints = 100, lwd = 2,
                          size = 2/3, xlab = "Pseudotime",
                          ylab = "Log(expression + 1)", border = FALSE,
                          alpha = 2/3, sample = 1, pointCol = NULL,
                          curvesCols = NULL, plotLineages = TRUE, 
                          lineagesToPlot = NULL)
{
  
  #input is singleCellExperiment object.
  if (is.null(names(models))) {
    rownames(models) <- rownames(counts) <- seq_len(nrow(models))
    message(paste0(
      "The sce object has no rownames. Assuming that the counts and the sce ",
      "objects are ordered in the same way"))
  }
  #  if (length(gene) > 1) stop("Only provide a single gene's ID with the ",
  #              "gene argument.")
  # check if all gene IDs provided are present in the models object.
  if (is(gene, "character")) {
    if (!all(gene %in% names(models))) {
      stop("The gene ID is not present in the models object.")
    }
    id <- which(names(models) %in% gene)
  } else id <- gene
  
  dm <- colData(models)$tradeSeq$dm # design matrix
  y <- unname(counts[names(models),][id,])
  X <- colData(models)$tradeSeq$X # linear predictor
  slingshotColData <- colData(models)$crv
  pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
                                       pattern = "pseudotime")]
  if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
  nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
  betaMat <- rowData(models)$tradeSeq$beta[[1]]
  beta <- betaMat[id,]
  
  
  #construct time variable based on cell assignments.
  lcol <- timeAll <- rep(0, nrow(dm))
  for (jj in seq_len(nCurves)) {
    for (ii in seq_len(nrow(dm))) {
      if (dm[ii, paste0("l", jj)] == 1) {
        timeAll[ii] <- dm[ii, paste0("t", jj)]
        lcol[ii] <- jj
      } else {
        next
      }
    }
  }
  
  if (!is.null(pointCol)) {
    if (length(pointCol) == 1) {
      col <- colData(models)[,pointCol]
    } else if (length(pointCol) == ncol(models)) {
      col <- pointCol
    } else {
      col <- lcol
      message(paste("pointCol should have length of either 1 or the number of cells,",
                    "reverting to default color scheme."))
    }
  } else {
    col <- lcol
  }
  
  # plot raw data
  df <- data.frame("time" = timeAll,
                   "gene_count" = y,
                   "pCol" = as.character(col),
                   "lineage" = as.character(lcol))
  
  rows <- sample(seq_len(nrow(df)), nrow(df) * sample, replace = FALSE)
  df <- df[rows, ]
  if(!is.null(lineagesToPlot)){
    df <- df[df$lineage %in% lineagesToPlot,]
  }
  p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
    labs(x = xlab, y = ylab) +
    theme_classic()
  if(is.null(pointCol)){
    p <- p +
      geom_point(size = size, aes(col = lineage)) +
      scale_color_viridis_d(alpha = alpha)
  } else {
    p <- p +
      geom_point(size = size, alpha = alpha, aes(col = pCol)) +
      scale_color_discrete() +
      labs(col = "Cell labels")
  }
  
  if (plotLineages) {
    if (!is.null(curvesCols)) {
      if (length(curvesCols) != nCurves) {
        curvesCols <- viridis::viridis(nCurves)
        message("Incorrect number of lineage colors. Default to viridis")
      }
    } else {
      curvesCols <- viridis::viridis(nCurves)
    }
    if(is.null(lineagesToPlot)){
      lineagesToPlot <- seq_len(nCurves)
    }
    for (jj in lineagesToPlot) {
      df <- .getPredictRangeDf(dm, jj, nPoints = nPoints)
      Xdf <- predictGAM(lpmatrix = X,
                        df = df,
                        pseudotime = pseudotime)
      yhat <-  c(exp(t(Xdf %*% t(beta)) + df$offset))
      if (border) {
        p <- p +
          geom_line(data = data.frame("time" = df[, paste0("t", jj)],
                                      "gene_count" = yhat,
                                      "lineage" = as.character(jj),
                                      "pCol" = as.character(jj)),
                    lwd = lwd + 1, colour = "white")
      }
      data <-data.frame("time" = df[, paste0("t", jj)],
                        "gene_count" = yhat,
                        "lineage" = as.character(jj),
                        "pCol" = as.character(jj))
      p <- p +
        geom_line(data = data,
                  lwd = lwd, col = curvesCols[jj])
    }
  }
  
  return(data)
}

# internal function of TradeSeq R package
.getPredictRangeDf <- function(dm, lineageId, conditionId = NULL, nPoints = 100){
  vars <- dm[1, ]
  if ("y" %in% colnames(vars)) {
    vars <- vars[!colnames(vars) %in% "y"]
    off <- 1
  } else {
    off <- 0
  }
  offsetId <- grep(x = colnames(vars), pattern = "offset")
  offsetName <- colnames(vars)[offsetId]
  offsetName <- substr(offsetName, start = 8, stop = nchar(offsetName) - 1)
  names(vars)[offsetId] <- offsetName
  # set all times on 0
  vars[, grep(colnames(vars), pattern = "t[1-9]")] <- 0
  # set all lineages on 0
  vars[, grep(colnames(vars), pattern = "l[1-9]")] <- 0
  # duplicate to nPoints
  vars <- rbind(vars, vars[rep(1, nPoints - 1), ])
  # set range of pseudotime for lineage of interest
  if (is.null(conditionId)) {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId, "($|_)"))
  } else {
    lineageIds <- grep(colnames(vars), pattern = paste0("l", lineageId,
                                                        "_", conditionId, "$"))
  }
  if (length(lineageIds) == 1){
    lineageData <- dm[dm[, lineageIds + off] == 1,
                      paste0("t", lineageId)]
  } else {
    lineageData <- dm[rowSums(dm[, lineageIds + off]) == 1,
                      paste0("t", lineageId)]
  }
  # make sure lineage starts at zero
  if(min(lineageData) / max(lineageData) < .01) {
    lineageData[which.min(lineageData)] <- 0
  }
  vars[, lineageIds] <- 1 / length(lineageIds)
  # set lineage
  vars[, paste0("t", lineageId)] <- seq(min(lineageData),
                                        max(lineageData),
                                        length = nPoints)
  # set offset
  vars[, offsetName] <- mean(dm[, grep(x = colnames(dm),
                                       pattern = "offset")])
  return(vars)
}

#Function of TradeSeq R package
predictGAM <- function(lpmatrix, df, pseudotime, conditions = NULL){
  # this function is an alternative of predict.gam(model, newdata = df, type = "lpmatrix")
  # INPUT:
  # lpmatrix is the linear predictor matrix of the GAM model
  # df is a data frame of values for which we want the lpmatrix
  # pseudotime is the n x l matrix of pseudotimes
  # conditions is the vector of conditions, if present.
  
  # if pseudotime is vector, make it a matrix.
  if(is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime,ncol=1)
  
  condPresent <- !is.null(conditions)
  if(condPresent) nConditions <- nlevels(conditions)
  
  # for each curve, specify basis function IDs for lpmatrix
  allBs <- grep(x = colnames(lpmatrix), pattern = "[0-9]):l[1-9]")
  
  if(!condPresent){
    lineages <- sub(pattern = "s\\(", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
    lineages <- sub(pattern = "\\):.*", replacement = "",
                    x = lineages)
    nCurves <- length(unique(lineages))
    for (ii in seq_len(nCurves)) {
      assign(paste0("id",ii), allBs[which(lineages == paste0("t", ii))])
    }
  } else if(condPresent){
    lineages <- sub(pattern = "s\\(t", replacement = "",
                    x = colnames(lpmatrix[,allBs]))
    lineages <- sub(pattern = "\\):.*", replacement = "",
                    x = lineages)
    nLineages <- length(unique(lineages))
    curves <- sub(pattern = ".*:l", replacement = "",
                  x = colnames(lpmatrix[,allBs]))
    curves <- sub(pattern = "\\..*", replacement = "",
                  x = curves)
    nCurves <- length(unique(curves))
    for (ii in seq_len(nLineages)) {
      for(kk in seq_len(nConditions))
        assign(paste0("id", ii, "_", kk), allBs[which(curves == paste0(ii, "_", kk))])
    }
  }
  
  
  # specify lineage assignment for each cell (i.e., row of lpmatrix)
  if(!condPresent){
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nCurves)) {
        if (!all(x[get(paste0("id", ii))] == 0)) {
          return(ii)
        }
      }
    })
  } else if(condPresent){
    # first number is lineage, second number is condition.
    lineageID <- apply(lpmatrix, 1, function(x){
      for (ii in seq_len(nLineages)) {
        # loop over lineages
        for(kk in seq_len(nConditions)){
          # loop over conditions
          if (!all(x[get(paste0("id", ii, "_", kk))] == 0)) {
            return(as.numeric(paste0(ii, kk)))
          }
        }
      }
    })
  }
  
  
  # fit splinefun for each basis function based on assigned cells
  if(!condPresent) {
    for (ii in seq_len(nCurves)) { # loop over curves
      for (jj in seq_len(length(allBs) / nCurves)) { #within curve, loop over basis functions
        assign(paste0("l",ii,".",jj),
               stats::splinefun(x = pseudotime[lineageID == ii, ii],
                                y = lpmatrix[lineageID == ii, #only cells for lineage
                                             get(paste0("id", ii))[jj]],
                                ties = mean)) #basis function
      }
    }
  } else if(condPresent) {
    for (ii in  seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        for (jj in seq_len(length(allBs) / (nLineages * nConditions))) {
          #within curve, loop over basis functions
          assign(paste0("l",ii, "_", kk,".",jj),
                 stats::splinefun(
                   x = pseudotime[lineageID == as.numeric(paste0(ii, kk)), ii],
                   y = lpmatrix[lineageID == as.numeric(paste0(ii, kk)), #only cells for lineage
                                get(paste0("id", ii, "_", kk))[jj]],
                   ties = mean)) #basis function
        }
      }
    }
  }
  
  
  # use input to estimate X for each basis function
  Xout <- matrix(0, nrow = nrow(df), ncol = ncol(lpmatrix))
  if(!condPresent){
    for (ii in seq_len(nCurves)) { # loop over curves
      if (all(df[, paste0("l", ii)] == 1)) { # only predict if weight = 1
        for (jj in seq_len(length(allBs) / nCurves)) { # within curve, loop over basis functions
          f <- get(paste0("l", ii, ".", jj))
          Xout[, get(paste0("id", ii))[jj]] <- f(df[, paste0("t", ii)])
        }
      }
    }
  } else if(condPresent){
    # for (ii in (seq_len(nCurves)[seq(2, nCurves, by=2)])/2) {
    for (ii in seq_len(nLineages)) {
      # loop over curves
      for(kk in seq_len(nConditions)){
        # loop over conditions
        if (all(df[, paste0("l", ii, "_", kk)] != 0)) { # only predict if weight = 1
          for (jj in seq_len(length(allBs) / (nLineages * nConditions))) { 
            # within curve, loop over basis functions
            f <- get(paste0("l", ii, "_", kk, ".", jj))
            Xout[, get(paste0("id", ii, "_", kk))[jj]] <- f(df[, paste0("t", ii)])
          }
        }
      }
    }
  }
  
  
  # add fixed covariates as in df
  dfSmoothID <- grep(x = colnames(df), pattern = "[t|l][1-9]")
  dfOffsetID <- grep(x = colnames(df), pattern = "offset")
  Xout[, -allBs] <- df[, -c(dfSmoothID, dfOffsetID)]
  
  # return
  colnames(Xout) <- colnames(lpmatrix)
  return(Xout)
}

# Extract all mouse gene sets
ExtractMouseGeneSets<- function(output_path){
  all_gene_sets <- msigdbr(species = "Mus musculus")
  map_df <- all_gene_sets %>% dplyr::distinct(ensembl_gene, entrez_gene,gene_symbol) %>% as.data.frame()
  return(map_df)
}

customize_parameters <- function(Vec,DEmarkers,organism,datatype,disease_phase, output_path){
  
  # Initialize dictionary
  enrichcl_dict <- dict()
  
  # Contains aggregated enrichment analysis S4 objects after all subanalyses
  enrichcl_list <- list()
  
  de_list <- list()
  
  analysis_res = Enrichment_Analysis(Vec[1],DEmarkers, "GO",datatype,organism,enrichcl_dict,enrichcl_list, de_list, disease_phase,output_path)
  print(paste0('Finish Enrichment_Analysis for GO ', Vec[1]))
  
  analysis_res = Enrichment_Analysis(Vec[2],DEmarkers, "GO",datatype,organism,analysis_res$enrichcl_dict,analysis_res$enrichcl_list,analysis_res$de_list, disease_phase,output_path)
  print(paste0('Finish Enrichment_Analysis for GO ', Vec[2]))
  
  params <- list("qscore_min" = min(analysis_res$enrichcl_dict$qscore_min),
                 "qscore_max" = max(analysis_res$enrichcl_dict$qscore_max),
                 "GeneRatioNum_min" = min(analysis_res$enrichcl_dict$GeneRatioNum_min),
                 "GeneRatioNum_max" = max(analysis_res$enrichcl_dict$GeneRatioNum_max),
                 "Count_min" = min(analysis_res$enrichcl_dict$Count_min),
                 "Count_max" = max(analysis_res$enrichcl_dict$Count_max))
  
  return(list("params" = params, "enrichcl_list" = analysis_res$enrichcl_list, "de_list" = analysis_res$de_list))
  
}

Enrichment_Analysis <- function(anal_type,DEmarkers, analysisType,analysisData,organismDB, enrichcl_dict,enrichcl_list, de_list, disease_phase,output_path) {
  # Filter markers on anal_type
  DEmarkers <- subset(DEmarkers, cluster== anal_type) 
  Ensembl_id_presence <- DEmarkers%>%mutate(gene=map(., ~str_detect(.x, 'ENSMUSG'))%>%pmap_lgl(any))
  
  if (any(sum(Ensembl_id_presence$gene))){
    ENS <- vapply(strsplit(DEmarkers$gene,"\\."), `[`, 1, FUN.VALUE=character(1))
    SYMBOL <- vapply(strsplit(DEmarkers$gene,"\\."), `[`, 2, FUN.VALUE=character(1))
    DEmarkers <- cbind(DEmarkers, ENS, SYMBOL) 
  }
  else{
    DEmarkers <- MapGenetoENSEMBL(DEmarkers,organismDB,output_path) 
    
  }
  if (analysisType == "GO" & analysisData == 'ENSEMBL') {
    # Perform GO Enrichment Analysis
    enrichObj <- enrichGO(gene  = unique(DEmarkers$ENS),
                          OrgDb         = organismDB,
                          keyType       = analysisData,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
    
    
  }
  else if (analysisType == "GO" & analysisData == 'SYMBOL') {
    # Perform GO Enrichment Analysis
    enrichObj <- enrichGO(gene  = unique(DEmarkers$SYMBOL),
                          OrgDb         = organismDB,
                          keyType       = analysisData,
                          ont           = "BP",
                          pAdjustMethod = "fdr",
                          pvalueCutoff  = 0.05, 
                          qvalueCutoff  = 0.05)
    
    
  }
  else if (analysisType == "KEGG" & analysisData == 'ncbi-geneid') {
    DEmarkers <- MapENSEMBLtoENTREZID(DEmarkers,organismDB)
    # Perform KEGG pathways Enrichwment Analysis
    enrichObj <- enrichKEGG(gene= unique(DEmarkers$ENTREZID),
                            organism= "human",
                            keyType       = analysisData,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)
  }
  # Add qscore, GeneRatio, analysis, color as S4 Class elements
  sc1 <- setClass('qscore', contains = 'enrichResult', slots = c(qscore = 'numeric'))
  enrichPlusqscore <- as(enrichObj, 'qscore')
  enrichPlusqscore@result$qscore <- unlist(lapply(enrichPlusqscore@result$p.adjust,function(x) -log(x,base=10)))
  
  sc2 <- setClass('analysis', contains = 'qscore', slots = c(grp = 'numeric'))
  enrichclPlusanal <- as(enrichPlusqscore, 'analysis')
  enrichclPlusanal@result$analysis <- c(replicate(length(enrichclPlusanal@result$Count), anal_type))
  
  sc3 <- setClass('GeneRatioNum', contains = 'analysis', slots = c(GeneRatioNum = 'numeric'))
  enrichclPlusgration <- as(enrichclPlusanal, 'GeneRatioNum')
  enrichclPlusgration@result$GeneRatioNum <- parse_ratio(c(enrichclPlusgration@result$GeneRatio))
  
  sc4 <- setClass('Color', contains = 'GeneRatioNum', slots = c(Color = 'numeric'))
  enrichcl <- as(enrichclPlusgration, 'Color')
  
  
  
  enrichcl@result$Color <- c(replicate(length(enrichcl@result$Count), "gray"))

  enrichcl_dict[["qscore_min"]] <- c(enrichcl_dict[["qscore_min"]],min(enrichcl@result$qscore))
  enrichcl_dict[["qscore_max"]] <- c(enrichcl_dict[["qscore_max"]],max(enrichcl@result$qscore))
  
  enrichcl_dict[["GeneRatioNum_min"]] <- c(enrichcl_dict[["GeneRatioNum_min"]],min(enrichcl@result$GeneRatioNum))
  enrichcl_dict[["GeneRatioNum_max"]] <- c(enrichcl_dict[["GeneRatioNum_max"]],max(enrichcl@result$GeneRatioNum))
  
  enrichcl_dict[["Count_min"]] <- c(enrichcl_dict[["Count_min"]],min(enrichcl@result$Count))
  enrichcl_dict[["Count_max"]] <- c(enrichcl_dict[["Count_max"]],max(enrichcl@result$Count))
  
  enrichcl_list <- c(enrichcl_list,enrichcl)
  DE_genes_set <- unlist(list(DEmarkers$avg_log2FC))
  names(DE_genes_set) <- unlist(list(DEmarkers$SYMBOL))
  
  de_list <- c(de_list,list(DE_genes_set))
  return(list("enrichcl_dict" = enrichcl_dict, "enrichcl_list" = enrichcl_list, "de_list" = de_list))
}

CreateOutputDirectory <- function(output_path,subDir){
  saving_path <- paste0(output_path,"/",subDir)
  
  # Create directory for storing the results
  dir.create(file.path(output_path, subDir), showWarnings = FALSE,recursive = TRUE)
  
  return(saving_path)
}

MapGenetoENSEMBL<- function(DEGmarkers, organismDB,output_path){
  
  mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  
  G_list<- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                 filters = 'mgi_symbol',
                 values=DEGmarkers$gene,
                 mart=mart, uniqueRows=T)
  
  DEGmarkers <- left_join(DEGmarkers, G_list,
                          by= c("gene"="mgi_symbol"))
  
  
  # Check number of unmapped gene symbols
  print(sum(is.na(DEGmarkers$ensembl_gene_id)))
  
  
  umappped_genes <- filter(DEGmarkers, is.na(ensembl_gene_id))$gene
  print("umappped_genes :")
  print(umappped_genes)
  
  print(sum(is.na(DEGmarkers$ensembl_gene_id)))
  
  DEGmarkers <- rename(DEGmarkers,c('ensembl_gene_id'='ENS'))
  
  return(DEGmarkers)
}

MapENSEMBLtoENTREZID <- function(DEgenes, organismDB){
  entrezids<-bitr(DEgenes$ENS, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organismDB)
  mapped_ids = entrezids[!duplicated(entrezids[c("ENSEMBL")]),]
  names(mapped_ids)[names(mapped_ids) == 'ENSEMBL'] <- 'ENS'
  DEgenes <-merge(DEgenes,mapped_ids, by = "ENS")
  
  return(DEgenes)
}

BarPlot <- function(enrichobject,GO_terms,palet){
  enrichobject <-enrichobject[enrichobject$Description %in% GO_terms,]
  enrichobject <- enrichobject %>% group_by(Description) %>% slice(which.min(p.adjust))
  
  enrichobject$Description <- factor(enrichobject$Description, levels = rev(GO_terms))
  enrichobject <- enrichobject[order(enrichobject$Description), ]
  
  enrichobject$qscore_min<- min(enrichobject$qscore)
  enrichobject$qscore_max <- max(enrichobject$qscore)
  
  ggplot(data=enrichobject, aes(x=Description, y=qscore, fill = CellType)) + xlab(NULL) +
    geom_bar(stat="identity",position="dodge",width= 0.6, colour = "black",show.legend = TRUE,aes(fill = CellType)) + coord_flip() +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 80)) +
    scale_y_continuous(expand = expansion(c(0,0)), limits = c(0.0,5 ),breaks = c(0,1,2,3,4,5)) +
    theme(legend.justification = "top",
          plot.title = element_text(hjust = 0.5,size = 15, face = "bold"),axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 0.5,colour = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.8,colour = "black"),
          axis.title.y = element_text(size = rel(2), angle = 0),
          axis.title.x = element_text(size = rel(1.5), angle = 0),
          axis.text = element_text(size = 8),
          panel.background = element_blank(), legend.title=element_blank()) +
     scale_fill_manual(values = palet) + ggtitle(paste0("Enriched in Ccl19-EYFP (DTR","\U2212)"))

}

Visualize_GeneSignatures_sc <- function(object, genes_list, slot_type, score_type,gene_signature){
  
  ## Values for slot_type
  # "data" -> raw or normalized expression values
  # "scale.data" -> scaled expression values
  # "counts" -> raw counts
  
  ## Values for score_type
  # "average.mean" ->  average mean gene expression
  # "module.score" -> module score
  
  
  # Define gradient pallet colors (min,max)
  pal <-RColorBrewer::brewer.pal(n = 9, name = "RdYlBu") |> rev()
  
  col_nam <- ""
  if (score_type == 'average.mean'){
    res<- Compute_average_mean_expression(object,genes_list, slot_type,gene_signature)
  }
  else{
    print('Unknown score type, please try again!')
    return(-1)
  }
  p <-FeaturePlot(object = res[[1]], features = res[[2]], cols = pal, raster = FALSE) 
  return(p)
  
}

