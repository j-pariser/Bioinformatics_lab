library(NeuronChat)
library(CellChat)
library(ggalluvial)
library(glue)
library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)

MY_netAnalysis_river_Neuron <- function(object, slot.name = "net", pattern = c("outgoing","incoming"), cutoff.1 = 0.25, cutoff.2 = 0.5,top.n1=1e4,top.n2=1e4,
                                 sources.use = NULL, targets.use = NULL, signaling = NULL,
                                 color.use = NULL, color.use.pattern = NULL, color.use.signaling = "grey50",
                                 do.order = FALSE, main.title = NULL,
                                 font.size = 2.5, font.size.title = 12){
    
    requireNamespace("ggalluvial")
    res.pattern <- methods::slot(object, 'net_analysis')$pattern[[pattern]]
    data1 = res.pattern$pattern$cell
    data2 = res.pattern$pattern$signaling


    data1$Contribution[data1$Contribution <= sort(data1$Contribution,decreasing = T)[min(dim(data1)[1],top.n1)]] <- 0
    data1$Contribution[data1$Contribution < cutoff.1] <- 0
    plot.data <- data1
    
    nPatterns<-length(unique(plot.data$Pattern))
    nCellGroup<-length(unique(plot.data$CellGroup))
    cells.level = levels(object@idents)
    
    if (is.null(color.use)) {
        color.use <- CellChat::scPalette(length(cells.level))[cells.level %in% unique(plot.data$CellGroup)]
    }
    if (is.null(color.use.pattern)){
        color.use.pattern <- ggPalette(nPatterns)
    }
    if (!is.null(sources.use)) {
        if (is.numeric(sources.use)) {
            sources.use <- cells.level[sources.use]
        }
        plot.data <- subset(plot.data, CellGroup %in% sources.use)
    }
    if (!is.null(targets.use)) {
        if (is.numeric(targets.use)) {
            targets.use <- cells.level[targets.use]
        }
        plot.data <- subset(plot.data, CellGroup %in% targets.use)
    }
    
    ## connect cell groups with patterns--------------------------------------------------------------------------------------------------
    
    plot.data.long <- to_lodes_form(plot.data, axes = 1:2, id = "connection")

    if (do.order) {
        mat = tapply(plot.data[["Contribution"]], list(plot.data[["CellGroup"]], plot.data[["Pattern"]]), sum)
        d <- dist(as.matrix(mat))
        hc <- hclust(d, "ave")
        k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
        cluster <- hc %>% cutree(k)
        order.name <- order(cluster, decreasing = TRUE)
        plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(names(cluster)[order.name], colnames(mat)))
        color.use <- color.use[order.name]
    }

    
    color.use.all <- c(color.use, color.use.pattern)
    StatStratum <- ggalluvial::StatStratum
    gg1 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("CellGroup", "Pattern")),y=Contribution,
                                        stratum = stratum, alluvium = connection,
                                        fill = stratum, label = stratum)) +
        geom_flow(width = 1/3,aes.flow = "backward") +
        geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +
        geom_text(stat = "stratum", size = font.size) +
        scale_x_discrete(limits = c(),  labels=c("Cell groups", "Patterns")) +
        scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
        theme_bw()+
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text.y= element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor  = element_blank(),
              panel.border = element_blank(),
              axis.ticks = element_blank(),axis.text=element_text(size=10)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    ## connect patterns with signaling----------------------------------------------------------------------------------------------------
    data2$Contribution[data2$Contribution <= sort(data2$Contribution,decreasing = T)[min(dim(data2)[1],top.n2)]] <- 0
    data2$Contribution[data2$Contribution < cutoff.2] <- 0
    plot.data <- data2
    nPatterns<-length(unique(plot.data$Pattern))
    nSignaling<-length(unique(plot.data$Signaling))
    if (length(color.use.signaling) == 1) {
    color.use.all <- c(color.use.pattern, rep(color.use.signaling, nSignaling))
    } else {
    color.use.all <- c(color.use.pattern, color.use.signaling)
    }
    if (!is.null(signaling)) {
    plot.data <- plot.data[plot.data$Signaling %in% signaling, ]
    }
    
    plot.data.long <- ggalluvial::to_lodes_form(plot.data, axes = 1:2, id = "connection")

    if (do.order) {
        mat = tapply(plot.data[["Contribution"]], list(plot.data[["Signaling"]], plot.data[["Pattern"]]), sum)
        mat[is.na(mat)] <- 0; mat <- mat[-which(rowSums(mat) == 0), ]
        d <- dist(as.matrix(mat))
        hc <- hclust(d, "ave")
        k <- length(unique(grep("Pattern", plot.data.long$stratum[plot.data.long$Contribution != 0], value = T)))
        cluster <- hc %>% cutree(k)
        order.name <- order(cluster)
        plot.data.long$stratum <- factor(plot.data.long$stratum, levels = c(colnames(mat),names(cluster)[order.name]))

    }
    
    gg2 <- ggplot(plot.data.long,aes(x = factor(x, levels = c("Pattern", "Signaling")), y= Contribution,
                                    stratum = stratum, alluvium = connection,
                                    fill = stratum, label = stratum)) +
    geom_flow(width = 1/3,aes.flow = "forward") +
    geom_stratum(width=1/3,size=0.1,color="black", alpha = 0.8, linetype = 1) +

    geom_text(stat = "stratum", size = font.size) + # 2.5
    scale_x_discrete(limits = c(),  labels=c("Patterns", "Signaling")) +
    scale_fill_manual(values = alpha(color.use.all, alpha = 0.8), drop = FALSE) +
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor  = element_blank(),
          panel.border = element_blank(),
          axis.ticks = element_blank(),axis.text=element_text(size= 10))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
    
    gg <- cowplot::plot_grid(gg1, gg2, align = "h", nrow = 1)
    title <- cowplot::ggdraw() + cowplot::draw_label(main.title,size = font.size.title)
    gg <- cowplot::plot_grid(title, gg, ncol=1, rel_heights=c(0.1, 1))


    return(gg)
}