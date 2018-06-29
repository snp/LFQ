#' Title
#'
#' @param data
#' @param method
#'
#' @return
#' @export
#'
#' @examples
normalizeLF <- function(data, method = "lm") {
    numSamples <- length(unique(data$Sample))
    good_proteins <- data %>% filter(!is.na(Value)) %>% group_by(id) %>% summarize(score = n_distinct(Sample)) %>% filter(score == numSamples)

    message("Proteins for normalization: ", nrow(good_proteins))

    mean.data <- data %>% filter(id %in% good_proteins$id) %>% group_by(id) %>% summarize(total = mean(Value, na.rm=T))
    norm.data <- data %>% full_join(mean.data, by = "id") %>% group_by(Sample) %>% do({
        data_ <- .
        model <- lm(total ~ Value, data = data_)
        print(.$Sample[[1]])
        print(model)
        all_id <- data_$id
        all_Value <- data_$Value
        all_normValue <- all_Value * coef(model)[[2]] + coef(model)[[1]]
        data_frame(Sample = .$Sample[[1]], id = all_id, intercept = coef(model)[[1]], slope = coef(model)[[2]], Value = all_Value, normValue = all_normValue)
    })


    data_norm <- data %>% left_join(norm.data %>% select(id, Sample, normValue), by = c("id", "Sample")) %>% mutate(Value_ = Value) %>% mutate(Value = normValue) %>% select(-normValue)
    if("File" %in% names(data_norm) & nrow(data_norm %>% distinct(File))>1){
        irs <- data_norm %>%
            group_by(id, File ) %>%
            summarize(irs = mean(Value)) %>%
            mutate(average=mean(irs)) %>% ungroup() %>%
            mutate(fac = average - irs)

        data_norm <- data_norm %>%
            left_join(irs, by=c('id','File')) %>%
            mutate(Value = Value + fac) %>% select(-irs, -average)

    }
    data_norm
}
