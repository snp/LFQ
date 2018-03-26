#' Title
#'
#' @param data
#' @param control
#' @param sample
#' @param var.equal
#' @param paired
#' @param alternative
#'
#' @return
#' @export
#'
#' @examples
lfq_ttest <- function(data, control, sample, var.equal = TRUE, paired = FALSE, alternative = "two.sided") {
    # compare_ <-
    data %>%
      dplyr::mutate(compare = ifelse(Sample %in% control,
                                            "Control",
                                            ifelse(Sample %in% sample,
                                                   "Sample",
                                                   NA))) %>%
      dplyr::filter(!is.na(compare) & !is.na(Value)) %>%
      dplyr::group_by(id) %>% do({
        data_ <- .
        if (all(table(data_$compare) > 1) & length(unique(data_$compare)) == 2) {
            tested <- t.test(Value ~ compare,
                             data = data_,
                             paired = paired,
                             var.equal = var.equal,
                             alternative = alternative)
            broom::tidy(tested) %>% dplyr::mutate(id = data_$id[[1]])
        } else {
            data.frame()
        }
    }) -> data_tested
    if(nrow(data_tested)>1)
      data_tested$`p.adj` <- p.adjust(data_tested$p.value, method='fdr')
}
