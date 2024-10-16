rankNet_Neuron <- function(object, slot.name = "net", measure = c("weight","count"), mode = c("comparison", "single"), comparison = c(1,2), color.use = NULL, stacked = FALSE, sources.use = NULL, targets.use = NULL,  signaling = NULL, pairLR = NULL, do.stat = FALSE, cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, show.raw = FALSE, return.data = FALSE, x.rotation = 90, title = NULL, bar.w = 0.75, font.size = 8,
                       do.flip = TRUE, x.angle = NULL, y.angle = 0, x.hjust = 1,y.hjust = 1,
                       axis.gap = FALSE, ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9,0,0.1)) {
  # slot.name = "net";measure = c("weight");mode = c("comparison");comparison = c(3,15);color.use = NULL;stacked = FALSE;sources.use = NULL; targets.use = NULL;  signaling = NULL;pairLR = NULL;do.stat = FALSE; cutoff.pvalue = 0.05;tol = 0.05; thresh = 0.05; show.raw = FALSE;return.data = FALSE; x.rotation = 90;title = NULL; bar.w = 0.75; font.size = 12;
  # do.flip = TRUE; x.angle = NULL; y.angle = 0; x.hjust = 1; y.hjust = 1;
  # axis.gap = FALSE; ylim = NULL; segments = NULL; tick_width = NULL; rel_heights = c(0.9,0,0.1)
  # ylim=c(0,90); axis.gap = TRUE; segments=c(10, 50);tick_width = c(1,10);rel_heights = c(0.9,0,0.1)
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (measure == "weight") {
    ylabel = "Information flow"
  } else if (measure == "count") {
    ylabel = "Number of links"
  }
  if (mode == "single") {
    object1 <- simplify2array(methods::slot(object, slot.name))
    object.names <- 'Single'
    prob = object1
    #prob[object1 > thresh] <- 0
    if (measure == "count") {
      prob <- 1*(prob > 0)
    }
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        } else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[1]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        } else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[ ,idx.t, ] <- 0
    }
    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }

    pSum <- apply(prob, 3, sum)
    pSum.original <- pSum
    if (measure == "weight") {
      # pSum <- -1/log(pSum)
      # pSum[is.na(pSum)] <- 0
      # idx1 <- which(is.infinite(pSum) | pSum < 0)
      # values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
      # position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      # pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    } else if (measure == "count") {
      pSum <- pSum.original
    }

    pair.name <- names(pSum)

    df<- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }

    if ((slot.name == "net") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "net") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }

    gg <- ggplot(df, aes(x=name, y=contribution.scaled)) + geom_bar(stat="identity",width = bar.w) +
      theme_classic() + theme(axis.text=element_text(size=10),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab(ylabel) + coord_flip()#+
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }

  } else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- simplify2array(object.list)
      #prob[object.list$pval > thresh] <- 0
      if (measure == "count") {
        prob <- 1*(prob > 0)
      }

      ##############################################
      ## This performs our normalization for DISEASED. 
      ## Control is always 1st. Diseased is always 2nd.
      ##############################################
      if (i == 2) {
          normalization = sum(simplify2array(methods::slot(object, slot.name)[[comparison[1]]])) / 
                          sum(simplify2array(methods::slot(object, slot.name)[[comparison[2]]]))
          prob <- prob * normalization
      }
      ##############################################
        
      prob.list[[i]] <- prob
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          } else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          } else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[ ,idx.t, ] <- 0
      }
      if (sum(prob) == 0) {
        stop("No inferred communications for the input!")
      }
      pSum.original[[i]] <- apply(prob, 3, sum)
      if (measure == "weight") {
        pSum[[i]] <- pSum.original[[i]]
        # pSum[[i]] <- -1/log(pSum.original[[i]])
        # pSum[[i]][is.na(pSum[[i]])] <- 0
        # idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
        # pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      } else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]]
      }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, object.names[comparison[i]])
    }
    # if (measure == "weight") {
    #   values.assign <- seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))
    #   position <- sort(pSum.original.all, index.return = TRUE)$ix
    #   for (i in 1:length(comparison)) {
    #     if (i == 1) {
    #       pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
    #     } else {
    #       pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1]))+1:length(unlist(idx[1:i])), position)]
    #     }
    #   }
    # }



    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0, group = object.names[comparison[i]], row.names = pair.name.all)
      df[[i]][pair.name[[i]],3] <- pSum[[i]]
      df[[i]][pair.name[[i]],2] <- pSum.original[[i]]
    }


    # contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution/abs(df[[1]]$contribution), digits=2))
    # #  contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution.scaled/abs(df[[1]]$contribution.scaled), digits=1))
    # contribution.relative2 <- as.numeric(format(df[[length(comparison)-1]]$contribution/abs(df[[1]]$contribution), digits=2))
    # contribution.relative[is.na(contribution.relative)] <- 0
    # for (i in 1:length(comparison)) {
    #   df[[i]]$contribution.relative <- contribution.relative
    #   df[[i]]$contribution.relative2 <- contribution.relative2
    # }
    # df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    # idx <- with(df[[1]], order(-contribution.relative,  -contribution.relative2, contribution, -contribution.data2))

    contribution.relative <- list()
    for (i in 1:(length(comparison)-1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison)-i+1]]$contribution/df[[1]]$contribution, digits=2))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, contribution, -contribution.data2))
    } else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2,contribution, -contribution.data2))
    } else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, contribution, -contribution.data2))
    } else {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, -contribution.relative.4, contribution, -contribution.data2))
    }



    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL

    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)

    if (is.null(color.use)) {
      color.use =  ggPalette(length(comparison))
    }

    # https://stackoverflow.com/questions/49448497/coord-flip-changes-ordering-of-bars-within-groups-in-grouped-bar-plot
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)

    # perform statistical analysis
    # if (do.stat) {
    #   pvalues <- c()
    #   for (i in 1:length(pair.name.all)) {
    #     df.prob <- data.frame()
    #     for (j in 1:length(comparison)) {
    #       if (pair.name.all[i] %in% pair.name[[j]]) {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(prob.list[[j]][ , , pair.name.all[i]]), group = comparison[j]))
    #       } else {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(matrix(0, nrow = nrow(prob.list[[j]]), ncol = nrow(prob.list[[j]]))), group = comparison[j]))
    #       }
    #
    #     }
    #     df.prob$group <- factor(df.prob$group, levels = comparison)
    #     if (length(comparison) == 2) {
    #       pvalues[i] <- wilcox.test(prob ~ group, data = df.prob)$p.value
    #     } else {
    #       pvalues[i] <- kruskal.test(prob ~ group, data = df.prob)$p.value
    #     }
    #   }
    #   df$pvalues <- pvalues
    # }
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          stop("Statistical test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE`")
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][ , , pair.name.all[i]])
          } else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) >3 & sum(is.na(prob.values)) == 0) {
          pvalues <- wilcox.test(prob.values[ ,1], prob.values[ ,2], paired = TRUE)$p.value
        } else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }


    if (length(comparison) == 2) {
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2], ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
      } else {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[2], ifelse(df$contribution.relative > 1+tol, color.use[1], "black"))
      }
    } else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }

    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), ]
      }
    }


    ################################################################
    # Take out Black *NS* pathways
    condition <- ((df$contribution.relative < 1-tol) | (df$contribution.relative > 1+tol))
    df <- subset(df, condition)
    colors.text <- colors.text[condition]
    ################################################################


    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }

    if (stacked) {
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position ="fill") # +
      # xlab("") + ylab("Relative information flow") #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      #  scale_y_discrete(breaks=c("0","0.5","1")) +
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      } else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }

      gg <- gg + geom_hline(yintercept = 0.5, linetype="dashed", color = "black", size=0.5)
    } else {
      if (show.raw) {
        gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      } else {
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      }

      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg,
                       ylim = ylim,
                       segments = segments,
                       tick_width = tick_width,
                       rel_heights = rel_heights)
      }
    }
    gg <- gg +  CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text))
      if (is.null(x.angle)) {
        x.angle = 0
      }

    } else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      gg <- gg + scale_x_discrete(limits = rev) + theme(axis.text.x = element_text(colour = rev(colors.text)))

    }

    gg <- gg + theme(axis.text.x=element_text(size=font.size), axis.title.y = element_text(size=font.size))
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE))
    gg <- gg + theme(axis.text.x = element_text(size=font.size*2/3, angle = x.angle, hjust=x.hjust),
                     axis.text.y = element_text(size=font.size, angle = y.angle, hjust=y.hjust))
    gg <- gg + theme(panel.border = element_rect(colour = "black", fill=NA, size=1), 
                     axis.line = element_blank())
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.6))
    }
  }

  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}