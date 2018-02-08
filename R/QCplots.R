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
lfq_pca <- function(data, plot = TRUE, color = NA, shape=NA, plotly = TRUE) {
  color <- ifelse(is.na(color),"", color)
  shape <- ifelse(is.na(shape),"", shape)
  tdata <- data %>%
      dplyr::select(id, Sample, Value, one_of(color, shape)) %>%
      tidyr::spread(id, Value)
  # names(ttdata) <- rnames pca.res <- PCA(t(na.exclude(tdata)), graph=FALSE)
  pca.res <- FactoMineR::PCA(tdata %>% select(-one_of("Sample", color, shape)), graph = FALSE)
  pca.coord <- as.data.frame(pca.res$ind$coord)
  pca.loadings <- as.data.frame(pca.res$var$coord) %>%
    dplyr::mutate(id = rownames(pca.res$var$coord)) %>%
    dplyr::select(id, contains("Dim"))

  if (!is.na(color) & color != "") {
    pca.coord$color_ <- tdata[,color] %>% unlist(use.names = FALSE)
    color_guide <- color
  } else {
    pca.coord$color_ <- tdata$Sample
    color_guide <- "Sample"
  }
  if (!is.na(shape) & shape !="") {
      pca.coord$shape_ <- as.factor(tdata[,shape] %>% unlist(use.names = FALSE))

      shape_guide <- shape
  } else {
      pca.coord$shape_ <- ""
      shape_guide = ""
  }

  pca.plot.12 <- ggplot2::ggplot(pca.coord, ggplot2::aes(Dim.1, Dim.2, color = color_, shape = shape_, label=tdata$Sample)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::theme_minimal() +
    ggplot2::guides(color = ggplot2::guide_legend(title = color_guide),
                    shape = ggplot2::guide_legend(title = shape_guide)) +
    ggplot2::labs(title = "PCA plot",
         x = sprintf("PC1, %.2f%%", pca.res$eig[1, 2]),
         y = sprintf("PC2, %.2f%%", pca.res$eig[2, 2]))

  pca.plot.23 <- ggplot2::ggplot(pca.coord, ggplot2::aes(Dim.2, Dim.3, color = color_, shape = shape_, label=tdata$Sample)) +
    ggplot2::geom_point(size = 5) +
    ggplot2::guides(color = ggplot2::guide_legend(title = color_guide),
                    shape = ggplot2::guide_legend(title = shape_guide)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PCA plot",
         x = sprintf("PC2, %.2f%%", pca.res$eig[2, 2]),
         y = sprintf("PC3, %.2f%%", pca.res$eig[3, 2]))
  pca.plot.components <- as.data.frame(pca.res$eig) %>%
    dplyr::mutate(comp = as.factor(as.numeric(sub("comp ", "", rownames(pca.res$eig))))) %>%
    ggplot2::ggplot(ggplot2::aes(x=comp, y=`percentage of variance`)) +
    ggplot2::geom_col() + theme_bw()
  lfq_pca_plot <- ggpubr::ggarrange(pca.plot.12, pca.plot.components, common.legend = T, ncol = 1, nrow=2, legend = 'right')
  if(plot){
    if(!plotly) {
      print(lfq_pca_plot)
    }else{
      lfq_pca_plot <- plotly::subplot(plotly::ggplotly(pca.plot.12 + ggplot2::guides(color="none",shape="none")),
                                      plotly::ggplotly(pca.plot.components + ggplot2::guides(color="none",shape="none")),
                                      shareY = FALSE, nrows = 2)
      lfq_pca_plot
    }
  }
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
lfq_heatmap <- function(data, plot = TRUE, groups = NA, plotly=TRUE) {
  if(is.na(groups))
    groups <- ""

  tdata <- data %>%
    dplyr::select(id, Sample, Value, dplyr::one_of(groups)) %>%
    tidyr::spread(id, Value)

  if(groups %in% names(tdata)) {
    classvec_ <- unlist(tdata[,groups], use.names = F)
  }else{
    classvec_ <- unlist(tdata[,"Sample"], use.names = F)
    classvec_ <- NA
  }

  ttdata <- tdata %>% select(-Sample, -one_of(groups)) %>% t()
  colnames(ttdata) <- tdata$Sample

  if(plotly) {
    heatmaply::heatmaply(cor(ttdata, use = "pair"), symm = TRUE, dualScale = FALSE, cexRow = 0.6, cexCol = 0.6, classvec = classvec_, key = T)
  } else {
    made4::heatplot(cor(ttdata, use = "pair"), symm = TRUE, dualScale = FALSE, cexRow = 0.6, cexCol = 0.6, classvec = classvec_, key = T)
  }
}
