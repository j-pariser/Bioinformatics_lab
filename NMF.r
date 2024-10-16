library(NMF)

#' Identification of major signals for specific cell groups and general communication patterns
#'
#' Identification of major signals for specific cell groups and general communication patterns; adapted from CellChat https://github.com/sqjin/CellChat
#' @param object NeuronChat object
#' @param slot.name the slot name of object that is used to store communication strength matrices, i.e.,  'net'
#' @param pattern "outgoing" or "incoming"
#' @param k the number of patterns
#' @param k.range a range of the number of patterns
#' @param heatmap.show whether showing heatmap
#' @param color.use the character vector defining the color of each cell group
#' @param color.heatmap a color name in brewer.pal
#' @param title.legend the title of legend in heatmap
#' @param width width of heatmap
#' @param height height of heatmap
#' @param font.size fontsize in heatmap
#' @importFrom methods slot
#' @importFrom NMF nmfEstimateRank nmf
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom stats setNames
#' @importFrom grid grid.grabExpr grid.newpage pushViewport grid.draw unit gpar viewport popViewport
#'
#' @return
#' @export
#'
#' @examples

identifyCommunicationPatterns_Neuron <- function(object, slot.name = "net", pattern = c("outgoing","incoming"), k = NULL, k.range = seq(2,10), heatmap.show = TRUE,
                                                 color.use = NULL, color.heatmap = "Spectral", title.legend = "Contributions",
                                                 width = 4, height = 6, font.size = 8,thresh_quantile=0, 
                                                 chosen_method = "brunet", seed = "random", nrun = 200, parallel = 'p16') {
  # heatmap.show = TRUE;color.use = NULL;color.heatmap = "Spectral";title.legend = "Contributions";
  # width = 4; height = 6; font.size = 8
  pattern <- match.arg(pattern)
  prob <- simplify2array(methods::slot(object, slot.name))
  if (pattern == "outgoing") {
    data_sender <- apply(prob, c(1,3), sum)
    data_sender = sweep(data_sender, 2L, apply(data_sender, 2, function(x) max(c(x,1e-6), na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_sender)
  } else if (pattern == "incoming") {
    data_receiver <- apply(prob, c(2,3), sum)
    data_receiver = sweep(data_receiver, 2L, apply(data_receiver, 2, function(x) max(c(x,1e-6), na.rm = TRUE)), '/', check.margin = FALSE)
    data0 = as.matrix(data_receiver)
  }
  options(warn = -1)
  data <- data0
  data <- data[rowSums(data)!=0,colSums(data)!=0]
  data <- data[,colSums(data) >= quantile(colSums(data),thresh_quantile)]
  data <- data[rowSums(data)!=0,]
  if (is.null(k)) {
    stop("Please run the function `selectK` for selecting a suitable k!")
  }
                                                   
  outs_NMF <- NMF::nmf(data, rank = k, method = chosen_method, seed = seed, nrun = nrun, .opt=parallel) 
  W <- scaleMat(outs_NMF@fit@W, 'r1')
  H <- scaleMat(outs_NMF@fit@H, 'c1')
  colnames(W) <- paste0("Pattern ", seq(1,ncol(W))); rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));
  if (heatmap.show) {
    net <- W
    if (is.null(color.use)) {
      color.use <- CellChat::scPalette(length(rownames(net)))
    }
    color.heatmap = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(255)

    df<- data.frame(group = rownames(net)); rownames(df) <- rownames(net)
    cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
    row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))

    ht1 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                  left_annotation = row_annotation,
                  cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  show_heatmap_legend = F,
                  column_title = "Cell patterns",column_title_gp = gpar(fontsize = 10)
    )


    net <- t(H)

    ht2 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                  cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = "Communication patterns",column_title_gp = gpar(fontsize = 10),
                  heatmap_legend_param = list(title = title.legend, title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, at = c(round(min(net, na.rm = T), digits = 1), round(max(net, na.rm = T), digits = 1)),
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 6),grid_width = unit(2, "mm"))
    )

    gb_ht1 = grid.grabExpr(draw(ht1))
    gb_ht2 = grid.grabExpr(draw(ht2))
    #grid.newpage()
    pushViewport(viewport(x = 0.1, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
    grid.draw(gb_ht1)
    popViewport()

    pushViewport(viewport(x = 0.6, y = 0.1, width = 0.2, height = 0.5, just = c("left", "bottom")))
    grid.draw(gb_ht2)
    popViewport()

  }

  data_W <- as.data.frame(as.table(W)); colnames(data_W) <- c("CellGroup","Pattern","Contribution")
  data_H <- as.data.frame(as.table(H)); colnames(data_H) <- c("Pattern","Signaling","Contribution")

  res.pattern = list("cell" = data_W, "signaling" = data_H)
  methods::slot(object, 'net_analysis')$pattern[[pattern]] <- list(data = data0, pattern = res.pattern)
  return(object)
}