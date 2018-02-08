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
    good_proteins <- data %>% group_by(id) %>% summarize(score = sum(!is.na(Value))) %>% filter(score == numSamples)

    message("Proteins for normalization: ", nrow(good_proteins))

    mean.data <- data %>% filter(id %in% good_proteins$id) %>% group_by(id) %>% summarize(total = mean(Value))
    norm.data <- data %>% full_join(mean.data, by = "id") %>% group_by(Sample) %>% do({
        data_ <- .
        model <- lm(total ~ Value, data = data_)
        all_id <- data_$id
        all_Value <- data_$Value
        all_normValue <- all_Value * coef(model)[[2]] + coef(model)[[1]]
        data_frame(Sample = .$Sample[[1]], id = all_id, intercept = coef(model)[[1]], slope = coef(model)[[2]], Value = all_Value, normValue = all_normValue)
    })


    data %>% left_join(norm.data %>% select(id, Sample, normValue), by = c("id", "Sample")) %>% mutate(Value_ = Value) %>% mutate(Value = normValue) %>% select(-normValue)
}
