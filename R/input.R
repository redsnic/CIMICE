#' Create mutational matrix from MAF file
#' 
#' Read a MAF (Mutation Annotation Format) file and create
#' a Mutational Matrix combining gene and sample IDs.
#' 
#' @param path path to MAF file
#' @param ... other maftools::mutCountMatrix arguments
#' 
#' @return the mutational (sparse) matrix associated to the MAF file
#' 
#' @examples 
#' read_MAF(system.file("extdata", "paac_jhu_2014_500.maf", package = "CIMICE", mustWork = TRUE))
#' 
#' @export read_MAF
read_MAF <- function(path, ...){
    maf <- read.maf(maf = path)
    # extract mutational matrix
    mat <- maf %>% mutCountMatrix() %>% 
        apply(c(1,2), function(x) ifelse(x>0,1,0)) %>%
        t 
    # make boolean
    mat <- Matrix(mat, sparse = TRUE)
    mat
} 


#' Gradually read a file from disk
#'
#' This function creates a reader to read a text file in batches (or chunks). 
#' It can be used for very large files that cannot fit in RAM.
#'
#' @param file_path path to large file
#'
#' @return a list-object containing the function `read` to 
#' read lines from the given file, and `close` to close the
#' connection to the file stream.
#' 
#' @examples 
#' # open connection to file
#' reader <- chunk_reader(
#'     system.file("extdata", "paac_jhu_2014_500.maf", package = "CIMICE", mustWork = TRUE)
#' )
#' 
#' while(TRUE){
#'     # read a chunk
#'     chunk <- reader$read(10)
#'     if(length(chunk) == 0){
#'         break
#'     }    
#'     # --- process chunk ---
#' }
#' # close connection
#' reader$close()
#' 
#' @export chunk_reader
chunk_reader <- function(file_path){
    f <- file(file_path)
    open(f)
    
    list(
        # Returns a character vector of the read lines.
        # The returned vector has length 0 if there is nothing to read.
        # This function keeps track of the current position in the file
        read = function(max_chunk = 10000){
            readLines(f, n = max_chunk)
        },
        # close connection to the file stream
        close = function(){
            close(f)
        }
    )
}

#' Add samples and genes names to a mutational matrix 
#'
#' Given M mutational matrix, add samples as row names, and genes as 
#' column names. If there are repetitions in row names, these are
#' solved by adding a sequential identifier to the names.
#' 
#' @param M mutational matrix
#' @param samples list of sample names
#' @param genes list of gene names
#'
#' @return N with the set row and column names
#' 
#' @examples 
#' require(Matrix)
#' genes <- c("A", "B", "C")
#' samples <- c("S1", "S2", "S2")
#' M <- Matrix(c(0,0,1,0,0,1,0,1,1), ncol=3, sparse=TRUE, byrow = TRUE)
#' 
#' annotate_mutational_matrix(M, samples, genes)
#' 
#' @export annotate_mutational_matrix
annotate_mutational_matrix <- function(M, samples, genes){
    assertthat::are_equal(ncol(M), length(genes))
    assertthat::are_equal(nrow(M), length(samples))
    # manage repeated rownames
    tryCatch({
        rownames(M) <- samples
    },
    error = function(e){
        # in case of repetitions, add a sequential number to the sample labels
        rownames(M) <- map2_chr(samples, seq(1,nrow(M)), ~ paste(.x, .y, sep="_"))
    })
    colnames(M) <- genes
    M
}


#' Read "CAPRI" file from a string
#'
#' Read a "CAPRI" formatted file, from a text string
#'
#' @param txt string in valid "CAPRI" format
#'
#' @return the described mutational matrix as a (sparse) matrix
#'
#' @examples
#' read_CAPRI_string("
#' s\\g A B C D
#' S1 0 0 0 1
#' S2 1 0 0 0
#' S3 1 0 0 0
#' S4 1 0 0 1
#' S5 1 1 0 1
#' S6 1 1 0 1
#' S7 1 0 1 1
#' S8 1 1 0 1
#' ")
#'
#' @export read_CAPRI_string
read_CAPRI_string <- function(txt){
    # read string as csv
    df <- read.csv(text=txt, sep="", strip.white = TRUE, blank.lines.skip = TRUE, row.names = NULL)
    # separate genes, samples and mutational matrix
    samples <- (df[1,])[-1]
    genes <- (df[,1])[-1]
    mutmatrix <- df[-1,-1] %>% as.matrix %>% Matrix(sparse = TRUE) 
    # glue components together
    annotate_mutational_matrix(mutmatrix, genes, samples)
}

#' Read a "CAPRI" file
#'
#' Read a "CAPRI" formatted file from the file system
#'
#' @param filepath path to file
#'
#' @return the described mutational matrix as a (sparse) matrix
#'
#' @examples
#' #          "pathToDataset/myDataset.CAPRI"
#' read_CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))
#'
#' @export read_CAPRI
read_CAPRI <- function(filepath){
    reader <- chunk_reader(filepath)
    
    # read header, remove first column name
    genes <- (reader$read(1) %>% strsplit("\\s+"))[[1]][-1] 
    mutmatrix <- Matrix(, ncol = length(genes), nrow=0)
    samples <- c()
    
    while(TRUE){
        chunk <- reader$read()
        if(length(chunk) == 0){
            break
        }
        
        # convert to list of char array
        chunk <- map(chunk, ~ strsplit(., "\\s+")[[1]])
        
        samples <- c(samples, map_chr(chunk, ~ .[1]))
        
        new_rows <- map(chunk, ~ .[-1]) %>% unlist %>% map_int(~strtoi(.))
        mutmatrix <- rbind(mutmatrix, Matrix(new_rows, sparse = TRUE, ncol = length(genes))) 
    }
    
    # close connection
    reader$close()
    annotate_mutational_matrix(mutmatrix, samples, genes)
}

#' Read a "CAPRI" file
#'
#' Read a "CAPRI" formatted file, as read_CAPRI
#'
#' @param filepath path to file
#'
#' @return the described mutational matrix as a (sparse) matrix
#'
#' @examples
#' read(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))
#'
#' @export read
read <- read_CAPRI

#' Read a "CAPRIpop" file
#'
#' Read a "CAPRIpop" formatted file from the file system
#'
#' @param filepath path to file
#'
#' @return a list containing the described mutational matrix as a (sparse) matrix
#'         and a list of the frequency of the genotypes
#'
#' @examples
#' #          "pathToDataset/myDataset.CAPRI"
#' read_CAPRI(system.file("extdata", "example.CAPRIpop", package = "CIMICE", mustWork = TRUE))
#'
#' @export read_CAPRIpop
read_CAPRIpop <- function(filepath){
    reader <- chunk_reader(filepath)
    
    # read header, remove first column name
    genes <- (reader$read(1) %>% strsplit("\\s+"))[[1]][-1] 
    mutmatrix <- Matrix(, ncol = length(genes), nrow=0)
    rown <- c()
    counts <- c()
    
    while(TRUE){
        chunk <- reader$read()
        if(length(chunk) == 0){
            break
        }
        
        # convert to list of char array
        chunk <- map(chunk, ~ strsplit(., "\\s+")[[1]])
        
        samples <- c(samples, map_chr(chunk, ~ .[1]))
        counts <- c(counts, map_chr(chunk, ~ .[length(counts)]))
        
        new_rows <- map(chunk, ~ .[-1, -length(counts)]) %>% unlist %>% map_int(~strtoi(.))
        mutmatrix <- rbind(mutmatrix, Matrix(new_rows, sparse = TRUE, ncol = length(genes))) 
    }
    
    # close connection
    reader$close()
    # list [mutational matrix, frequencies]
    list(matrix = annotate_mutational_matrix(mutmatrix, samples, genes), counts = counts)
}


#' Dataset line by line construction: initialization
#'
#' Initialize a dataset for "line by line" creation
#'
#' @param ... gene names (do not use '"', the input
#' is automatically converted to strings)
#'
#' @return a mutational matrix without samples structured as (sparse) matrix
#'
#' @examples
#' make_dataset(APC,P53,KRAS)
#'
#' @export make_dataset
make_dataset <- function(...){
    mutmatrix <- Matrix(, ncol = length(enexprs(...)), nrow=0)
    colnames(mutmatrix) <- enexprs(...)
    mutmatrix
}

#' Dataset line by line construction: add a sample
#'
#' Add a sample (a row) to an existing dataset.
#' This procedure is meant to be used with the "%>%" operator
#'
#' @param mutmatrix an existing (sparse) matrix (mutational matrix)
#' @param sampleName the row (sample) name
#' @param ... sample's genotype (0/1 numbers)
#'
#' @return the modified (sparse) matrix (mutational matrix)
#'
#' @examples
#'
#' require(dplyr)
#' make_dataset(APC,P53,KRAS)   %>%
#'     update_df("S1", 1, 0, 1) %>%
#'     update_df("S2", 1, 1, 1)
#'
#'
#' @export update_df
update_df <- function(mutmatrix, sampleName, ...){
    samples <- rownames(mutmatrix)
    mutmatrix <- rbind(mutmatrix, unlist(enexprs(...)))
    rownames(mutmatrix) <- c(samples, sampleName)
    mutmatrix
}



