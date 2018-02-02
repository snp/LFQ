#' Title
#'
#' @param data
#' @param plot
#' @param groups
#'
#' @return
#' @export
#'
#' @examples
lfq_pca <- function(data, plot = TRUE, groups = NA) {
    tdata <- data %>%
      dplyr::select(id, Sample, Value) %>%
      tidyr::spread(Sample, Value)
    ttdata <- as.data.frame(t(dplyr::select(tdata, -id)))
    # names(ttdata) <- rnames pca.res <- PCA(t(na.exclude(tdata)), graph=FALSE)
    pca.res <- FactoMineR::PCA(ttdata, graph = FALSE)
    pca.coord <- as.data.frame(pca.res$ind$coord)
    pca.loadings <- as.data.frame(pca.res$var$coord) %>%
      dplyr::mutate(id = tdata$id) %>%
      dplyr::select(id, contains("Dim"))
    if ((class(groups) != "list") & !is.na(groups))
        groups <- list(groups)

    if (!is.na(groups)) {
      pca.coord$color_ <- groups[[1]]
      color_guide = names(groups)[1]
    } else {
      pca.coord$color_ <- unique(data$Sample)
      shape_guide = "Sample"
    }
    if (length(groups) > 1) {
        pca.coord$shape_ <- groups[[2]]
        shape_guide = names(groups)[2]
    } else {
        pca.coord$shape_ <- ""
        shape_guide = ""
    }

    pca.plot.12 <- ggplot(pca.coord, aes(Dim.1, Dim.2, color = color_, shape = shape_)) +
      geom_point(size = 5) +
      theme_minimal() +
      guides(color = guide_legend(title = color_guide),
             shape = guide_legend(title = shape_guide)) +
      labs(title = "PCA plot",
           x = sprintf("PC1, %.2f%%", pca.res$eig[1, 2]),
           y = sprintf("PC2, %.2f%%", pca.res$eig[2, 2]))

    pca.plot.23 <- ggplot(pca.coord, aes(Dim.2, Dim.3, color = color_, shape = shape_)) +
      geom_point(size = 5) +
      theme_minimal() +
      guides(color = guide_legend(title = color_guide),
             shape = guide_legend(title = shape_guide)) +
      labs(title = "PCA plot",
           x = sprintf("PC2, %.2f%%", pca.res$eig[2, 2]),
           y = sprintf("PC3, %.2f%%", pca.res$eig[3, 2]))

    lfq_pca_plot <- ggpubr::ggarrange(pca.plot.12, pca.plot.23, common.legend = T, ncol = 1, nrow=2, legend = 'right')
    if(plot)
      print(lfq_pca_plot)
    return(list(scores = pca.coord, loadings = pca.loadings, plot = lfq_pca_plot))
}


#' Title
#'
#' @param data
#' @param plot
#' @param groups
#'
#' @return
#' @export
#'
#' @examples
lfq_heatmap <- function(data, plot = TRUE, groups = NA) {
    tdata <- data %>% select(id, Sample, Value) %>% spread(Sample, Value) %>% select(-id)
    made4::heatplot(cor(tdata, use = "pair"), symm = TRUE, dualScale = FALSE, cexRow = 0.6, cexCol = 0.6, classvec = groups, key = T)
}
