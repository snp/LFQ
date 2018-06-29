#' Title
#'
#' @param proteinGroups
#' @param idVar
#' @param qPrefix
#' @param removeReverse
#' @param removeConaminants
#' @param cleanup
#'
#' @return
#' @export
#'
#' @examples
importMQ_LFQ <- function(proteinGroups = "proteinGroups.txt", idVar = "Majority protein IDs", qPrefix = "LFQ intensity", removeReverse = TRUE, removeConaminants = TRUE,
    cleanup = TRUE) {

    # Read datafile
    if (!file.exists(proteinGroups)) {
        message("Can't find file ", proteinGroups, ". Check the location")
        return(data_frame())
    }
    data <- suppressMessages(readr::read_tsv(proteinGroups)) %>% mutate_(id = sprintf("`%s`", idVar)) %>% mutate(id = sub(";.*", "", id)) %>% mutate(id = gsub("[^a-zA-Z0-9-]+",
        "_", id)) %>% gather(Sample, Value, matches(sprintf("^%s\\s*(.+)$", qPrefix))) %>% mutate(Sample = sub(sprintf("^%s\\s*(.+?)$", qPrefix),
        "\\1", Sample), Value = ifelse(Value > 10, log2(Value), NA))

    message("Data loaded. ", length(unique(data$id)), " proteins, ", length(unique(data$Sample)), " samples.")
    if (cleanup) {
        data <- data %>% select(-contains("Identification"), -contains("Sequence"), -contains("Peptide"), -contains("Intensity"), -contains("MS/MS"),
            -contains("site"), -contains("IDs"), -contains("Mol. weight"), -contains("Fasta"), matches("^Unique peptides$"))
        message("Removing outliers...")
        # N <- length(unique(data$Sample))

        na.proteins <- data %>% group_by(id) %>% summarise(num.na = sum(is.na(Value))) %>% ungroup()
        na.proteins <- na.proteins %>% mutate(selected = num.na < (mean(na.proteins$num.na) + sd(na.proteins$num.na))) %>% select(id, selected)
        message("    ", sum(!na.proteins$selected), " proteins have too many missing values and will be removed.")
        data <- data %>% full_join(na.proteins, by = c("id")) %>% filter(selected == TRUE) %>% select(-selected)

        na.samples <- data %>% group_by(Sample) %>% summarise(num.na = sum(is.na(Value))) %>% ungroup()
        na.samples <- na.samples %>% mutate(selected = num.na < (mean(na.samples$num.na) + 1.5 * sd(na.samples$num.na))) %>% select(Sample, selected)
        message("    ", sum(!na.samples$selected), " samples have too many missing values and will be removed.")
        data <- data %>% full_join(na.samples, by = "Sample") %>% filter(selected) %>% select(-selected)
    }
    return(data)
}

importMQ_TMT <- function(proteinGroups = "proteinGroups.txt", idVar = "Majority protein IDs", qPrefix = "Reporter intensity corrected", removeReverse = TRUE, removeConaminants = TRUE,
                         cleanup = TRUE) {

  # Read datafile
  if (!file.exists(proteinGroups)) {
    message("Can't find file ", proteinGroups, ". Check the location")
    return(data_frame())
  }
  data <- suppressMessages(readr::read_tsv(proteinGroups)) %>%
    dplyr::mutate_(id = sprintf("`%s`", idVar)) %>%
    dplyr::mutate(id = sub(";.*", "", id)) %>%
    dplyr::mutate(id = gsub("[^a-zA-Z0-9-]+", "_", id)) %>%
    tidyr::gather(Sample, Value, matches(sprintf("^%s\\s*(\\d+)\\s+(.+)$", qPrefix))) %>%
    dplyr::mutate(Channel = as.numeric(sub(sprintf("^%s\\s*(\\d+)\\s+(.+)$", qPrefix), "\\1", Sample)),
                  Sample = sub(sprintf("^%s\\s*(\\d+)\\s+(.+)$", qPrefix), "\\2", Sample),
                  Value = ifelse(Value > 10, log2(Value), NA))

  message("Data loaded. ", length(unique(data$id)), " proteins, ", length(unique(data$Sample)), " samples.")
  if (cleanup) {
    data <- data %>% select(-contains("Identification"), -contains("Sequence"), -contains("Peptide"), -contains("Intensity"), -contains("MS/MS"),
                            -contains("site"), -contains("IDs"), -contains("Mol. weight"), -contains("Fasta"), -contains("Fraction"), matches("^Unique peptides$"))
    message("Removing outliers...")
    # N <- length(unique(data$Sample))

    na.proteins <- data %>% group_by(id) %>% summarise(num.na = sum(is.na(Value))) %>% ungroup()
    na.proteins <- na.proteins %>% mutate(selected = num.na <= (mean(na.proteins$num.na) + sd(na.proteins$num.na))) %>% select(id, selected)
    message("    ", sum(!na.proteins$selected), " proteins have too many missing values and will be removed.")
    data <- data %>% full_join(na.proteins, by = c("id")) %>% filter(selected == TRUE) %>% select(-selected)

    na.samples <- data %>% group_by(Sample, Channel) %>% summarise(num.na = sum(is.na(Value))) %>% ungroup()
    na.samples <- na.samples %>% mutate(selected = num.na <= (mean(na.samples$num.na) + 2.5 * sd(na.samples$num.na))) %>% select(Sample,Channel, selected)
    message("    ", sum(!na.samples$selected), " samples have too many missing values and will be removed.")
    data <- data %>% full_join(na.samples, by = c("Sample", "Channel")) %>% filter(selected) %>% select(-selected)
  }
  return(data)
}

importPD_TMT <- function(proteinGroups = "Proteins.txt", idVar = "Accession", qPrefix = "Abundance", removeReverse = TRUE, removeConaminants = TRUE,
                         cleanup = TRUE) {

  # Read datafile
  if (!file.exists(proteinGroups)) {
    message("Can't find file ", proteinGroups, ". Check the location")
    return(data_frame())
  }
  data <- suppressMessages(readr::read_tsv(proteinGroups)) %>%
    dplyr::mutate_(id = sprintf("`%s`", idVar)) %>%
    dplyr::mutate(id = sub(";.*", "", id)) %>%
    dplyr::mutate(id = gsub("[^a-zA-Z0-9-]+", "_", id)) %>%
    tidyr::gather(Sample, Value, matches(sprintf("^%s:? (F\\d+):?\\s*(\\d+[NCnc]?),?\\s+(.+)$", qPrefix))) %>%
    dplyr::mutate(Channel = sub(sprintf("^%s:? (F\\d+):?\\s*(\\d+[NCnc]?),?\\s+(.+)$", qPrefix), "\\2", Sample),
                  Label = sub(sprintf("^%s:? (F\\d+):?\\s*(\\d+[NCnc]?),?\\s+(.+)$", qPrefix), "\\3", Sample),
                  File = sub(sprintf("^%s:? (F\\d+):?\\s*(\\d+[NCnc]?),?\\s+(.+)$", qPrefix), "\\1", Sample),
                  Sample = sprintf("%s_%s", File, Channel),
                  Value = ifelse(Value > 10, log2(Value), NA))

  message("Data loaded. ", length(unique(data$id)), " proteins, ", length(unique(data$Sample)), " samples.")
  if (cleanup) {
    data <- data %>% select(-contains("Identification"), -contains("Sequence"), -contains("Peptide"), -contains("Intensity"), -contains("MS/MS"),
                            -contains("site"), -contains("IDs"), -contains("Mol. weight"), -contains("Fasta"), -contains("Fraction"), matches("^Unique peptides$"))
    message("Removing outliers...")
    # N <- length(unique(data$Sample))

    na.proteins <- data %>% group_by(id) %>% summarise(num.na = sum(is.na(Value))) %>% ungroup()
    na.proteins <- na.proteins %>% mutate(selected = num.na < (mean(na.proteins$num.na) + 3*sd(na.proteins$num.na))) %>% select(id, selected)
    message("    ", sum(!na.proteins$selected), " proteins have too many missing values and will be removed.")
    data <- data %>% full_join(na.proteins, by = c("id")) %>% filter(selected == TRUE) %>% select(-selected)

    na.samples <- data %>% group_by(Sample, Channel) %>% summarise(num.na = sum(is.na(Value))) %>% ungroup()
    na.samples <- na.samples %>% mutate(selected = num.na <= (mean(na.samples$num.na) + 3.5 * sd(na.samples$num.na))) %>% select(Sample,Channel, selected)
    message("    ", sum(!na.samples$selected), " samples have too many missing values and will be removed.")
    data <- data %>% full_join(na.samples, by = c("Sample", "Channel")) %>% filter(selected) %>% select(-selected)
  }
  return(data)
}
