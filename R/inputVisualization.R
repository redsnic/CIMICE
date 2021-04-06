utils::globalVariables(c("sums"))

#' Histogram of genes' frequencies
#'
#' Create the histogram of the genes' mutational frequencies
#'
#' @param mutmatrix input dataset (mutational matrix)
#' @param binwidth binwidth parameter for the histogram (as in ggplot)
#'
#' @examples
#' gene_mutations_hist(example_dataset(), binwidth = 10)
#'
#' @return the newly created histogram
#'
#' @export gene_mutations_hist
gene_mutations_hist <- function(mutmatrix, binwidth = 1){
    df <- mutmatrix %>% as.matrix %>% as.data.frame 
    # transpose dataframe
    temp.df <- as.data.frame(t(df[,seq(1,nrow(t(df)))]))
    # compute number of alteration per gene
    temp.df$sums <- rowSums(temp.df[,seq(1,ncol(temp.df))], na.rm = TRUE)
    # count number of genes
    nGenes <- nrow(temp.df)
    # create the histogram
    ggplot(temp.df, aes(x=sums)) +
        geom_histogram(color="black", fill="white", binwidth = binwidth) +
        labs(title = glue("Number of mutations per gene ({nGenes} genes)"),
        x = "Number of mutations", y = "Count")
}

#' Histogram of samples' frequencies
#'
#' Create the histogram of the samples' mutational frequencies
#'
#' @param mutmatrix input dataset (mutational matrix)
#' @param binwidth binwidth parameter for the histogram (as in ggplot)
#'
#' @examples
#' sample_mutations_hist(example_dataset(), binwidth = 10)
#'
#' @return the newly created histogram
#'
#' @export sample_mutations_hist
sample_mutations_hist <- function(mutmatrix, binwidth = 1){
    df <- mutmatrix %>% as.matrix %>% as.data.frame
    temp.df <- df
    # compute number of alteration per sample
    temp.df$sums <- rowSums(temp.df[,seq(1,ncol(df))], na.rm = TRUE)
    # count number of samples
    nSample <- nrow(temp.df)
    # create the histogram
    ggplot(temp.df, aes(x=sums)) +
        geom_histogram(color="black", fill="white", binwidth = binwidth) +
        labs(title = glue("Number of mutations per sample ({nSample} samples)"),
            x = "Number of mutations",
            y = "Counts")
}

#' Correlation plot from mutational matrix
#'
#' Prepare correlation plot based on a mutational matrix
#'
#' @param mutmatrix input dataset
#'
#' @return the computed correlation plot
#'
#' @examples
#' corrplot_from_mutational_matrix(example_dataset())
#'
#' @export corrplot_from_mutational_matrix
corrplot_from_mutational_matrix <- function(mutmatrix){
    df <- mutmatrix %>% as.matrix %>% as.data.frame
    # prepare correlation matrix
    corr <- round(cor(df), 1)
    corr[is.nan(corr)] <- 0
    corr[is.na(corr)] <- 0
    corr[corr == Inf] <- 0
    # prepare corrplot
    ggcorrplot(corr, hc.order = TRUE)
}

#' Gene based correlation plot
#'
#' Prepare a correlation plot computed from genes' perspective
#' using a mutational matrix
#'
#' @param mutmatrix input dataset (mutational matrix)
#'
#' @return the computed correlation plot
#'
#' @examples
#' corrplot_genes(example_dataset())
#'
#' @export corrplot_genes
corrplot_genes <- corrplot_from_mutational_matrix

#' Sample based correlation plot
#'
#' Prepare a correlation plot computed from samples' perspective
#' using a mutational matrix
#'
#' @param mutmatrix input dataset (mutational matrix)
#'
#' @return the computed correlation plot
#'
#' @examples
#' corrplot_samples(example_dataset())
#'
#' @export corrplot_samples
corrplot_samples <- function (mutmatrix)
    corrplot_from_mutational_matrix(t(mutmatrix))






