my_merged_river <- function(object, pattern, disease, inhibitory_celltypes, excitatory_celltypes, support_celltypes){
    if (pattern == "outgoing") {
      main.title = "Outgoing communication patterns of secreting cells"
      # The color palette is defined in the function ggPalette
      # color.use just needs three colors you can define them in a list with hex codes if you want 
      color.use <- ggPalette(3*2)[seq(1,3*2, by = 2)]
      color.use <- c("#F0E68C", "#66CDAA", "#D8BFD8")
      
    } else if (pattern == "incoming") {
      main.title = "Incoming communication patterns of target cells"
      # The color palette is defined in the function ggPalette
      # color.use just needs three colors you can define them in a list with hex codes if you want 
      color.use <- ggPalette(3*2)[seq(2,3*2, by = 2)]
      color.use <- c("#6A0DAD", "#FFD700", "#FF8C00")
    }

    slot.name = "net"
    res.pattern <- methods::slot(object, 'net_analysis')$pattern[[pattern]]
    
    data1 = res.pattern$pattern$cell
    data2 = res.pattern$pattern$signaling
    
    sorted_group <- data2[order(data2$Signaling, -data2$Contribution),]
    data2 <- sorted_group[!duplicated(sorted_group$Signaling),]
    
    sorted_group <- data1[order(data1$CellGroup, -data1$Contribution),]
    data1 <- sorted_group[!duplicated(sorted_group$CellGroup),]

    ##############################################################################################################################
    ### Throws out extra ligand-receptor pathways
    test_signal <- as.character(unique(data2$Signaling))
    test_signal <- sapply(strsplit(test_signal, "_"), `[`, 1)
    test_signal <- unique(test_signal)
    
    find_first_match <- function(small_elem, long_vec) {
      for (long_elem in long_vec) {
        if (startsWith(long_elem, paste0(small_elem, "_"))) {
          return(long_elem)
        }
      }
      return(NA)  # return NA if no match found
    }
    
    # Apply the function to each element in smaller_vector
    result <- sapply(test_signal, find_first_match, unique(data2$Signaling))
    data2 <- subset(data2, Signaling %in% result)
    data2$Signaling <- sapply(strsplit(as.character(data2$Signaling), split = "_"), "[", 1)
    ##############################################################################################################################

    # data2$Signaling <- gsub("_", "-", data2$Signaling)
    plot.data <- merge(data1,data2,by="Pattern")
    
    plot.data$Signaling <- as.character(plot.data$Signaling)
    print(unique(plot.data$Signaling))
    print(length(unique(plot.data$Signaling)))
    plot.data <- plot.data[order(plot.data$Pattern, plot.data$Signaling),]

    # order ligands by pattern
    if(pattern == "outgoing"){
      if(length(outgoing_ligand_order) > 0){
          plot.data$Signaling <- factor(plot.data$Signaling, levels = outgoing_ligand_order) 
      } 
      else {
        if(disease == "control"){
            outgoing_ligand_order <<- c()
            for (pattern in unique(plot.data$Pattern)){
              outgoing_ligand_order <<- c(outgoing_ligand_order, unique(plot.data$Signaling[plot.data$Pattern == pattern]))
            }
            plot.data$Signaling <- factor(plot.data$Signaling, levels = outgoing_ligand_order) 
        }
      }
    } else if(pattern == "incoming"){
      if(length(incoming_ligand_order) > 0){
          plot.data$Signaling <- factor(plot.data$Signaling, levels = incoming_ligand_order) 
      } 
      else {
        if(disease == "control"){
            incoming_ligand_order <<- c()
            for (pattern in unique(plot.data$Pattern)){
              incoming_ligand_order <<- c(incoming_ligand_order, unique(plot.data$Signaling[plot.data$Pattern == pattern]))
            }
            plot.data$Signaling <- factor(plot.data$Signaling, levels = incoming_ligand_order) 
        }
      }
    }

    #order cell types by cell type groups
    inhibitory_celltypes <- inhibitory_celltypes[order(inhibitory_celltypes)]
    excitatory_celltypes <- excitatory_celltypes[order(excitatory_celltypes)]
    celltypesOrder <- c(support_celltypes, inhibitory_celltypes, excitatory_celltypes)
    plot.data$CellGroup <- factor(plot.data$CellGroup, levels = celltypesOrder)
                                  
    #scale by cellgroup
    goal_sum <- 1
    for(celltype in unique(plot.data$CellGroup)){
        goal_sum <- pracma::Lcm(
            goal_sum,
            sum(plot.data$CellGroup == celltype)
        )
    }
    print(goal_sum)

    for(celltype in unique(plot.data$CellGroup)){
        index_to_copy <- which(match(plot.data$CellGroup, celltype) == 1)
        while(sum(plot.data$CellGroup == celltype) != goal_sum){
            for (index in index_to_copy){
                    plot.data <- rbind(plot.data, plot.data[rep(index, 1), ])
    
            }
        }
    }

    # Some signaling is NA, this will cause angles to be different length.
    plot.data <- subset(plot.data, !is.na(plot.data$Signaling))

    #make Pattern text Vertical
    angles <- rep(0, length(unique(plot.data$Pattern)) + length(unique(plot.data$CellGroup)) + length(unique(plot.data$Signaling)))
    for (i in seq(length(unique(plot.data$CellGroup))+1,length(unique(plot.data$CellGroup))+3)){
        angles[i] <- 90
    }

    #make our plot
    gg <- ggplot(plot.data,
        aes(axis1 = CellGroup, axis2 = Pattern, axis3 = Signaling)) +
    geom_flow(width = 1/3, aes.flow = "forward", aes(fill = Pattern)) + 
    scale_x_discrete(limits = c("Cell groups", "Pattern", "Signaling" )) +
    geom_stratum(alpha = 0.8, aes(fill = Pattern), width = 1/3, size=0.1) + 
    scale_fill_manual(values = color.use) +
    geom_text(angle=angles, size = 6, stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
        theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text.y= element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor  = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),axis.text=element_text(size=10)) + 
            theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(main.title)

    return(gg)
}