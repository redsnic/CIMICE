
library(dplyr)
library(glue)
library(ggcorrplot)

read.CAPRI.string <- function(txt){
    df <- read.csv(text=txt, sep = " ", row.names=NULL)
    colnames(df)[1] = "Samples"
    tryCatch({
        rownames(df) = df[,1]
    },
    error = function(e){
        rownames(df) = 1:length(df[,1])
    })
    df <- df %>% select(-Samples)
    df
}

# Lettura dell'input in formato "CAPRI"
read.CAPRI <- function(filename){
    df <- read.csv(filename, sep = " ", row.names=NULL)
    colnames(df)[1] = "Samples"
    tryCatch({
        rownames(df) = df[,1]
    },
    error = function(e){
        rownames(df) = 1:length(df[,1])
    })
    df <- df %>% select(-Samples)
    df
}
# CAPRI è il formato di default
read <- read.CAPRI

# crea un dataset linea per linea
# richiede in input i nomi dei geni
make.dataset <- function(...){
    df <- data.frame(matrix(ncol = length(enexprs(...)), nrow = 0))
    colnames(df) <- enexprs(...)
    df
}

# add a line to the dataframe
update.df <- function(d, sampleName, ...){
    rnames <- rownames(d)
    cnames <- colnames(d)
    df <- rbind(d, unlist(enexprs(...)))
    rownames(df) <- c(rnames, sampleName)
    colnames(df) <- cnames
    df
}

# istogramma delle frequenze di mutazione dei singoli geni
gene.mutations.hist <- function(df, binwidth = 1){
    # trasposizione del dataframe
    temp.df <- as.data.frame(t(df[,1:nrow(t(df))]))
    # calcolo del numero di mutazioni per i singoli geni nei campioni
    temp.df$sums <- rowSums(temp.df[,1:ncol(temp.df)], na.rm = TRUE)
    # numero di geni considerati
    nGenes <- nrow(temp.df)
    # l'istogramma
    ggplot(temp.df, aes(x=sums)) +
        geom_histogram(color="black", fill="white", binwidth = binwidth) +
        labs(title = glue("Number of mutations per gene ({nGenes} genes)"),
             x = "Number of mutations", y = "Count")
}

# istogramma del numero di mutazioni per singolo campione
sample.mutations.hist <- function(df, binwidth = 1){
    temp.df <- df
    # conteggio delle mutazioni
    temp.df$sums <- rowSums(temp.df[,1:ncol(df)], na.rm = TRUE)
    # numero di campioni considerati
    nSample <- nrow(temp.df)
    # l'istogramma
    ggplot(temp.df, aes(x=sums)) +
        geom_histogram(color="black", fill="white", binwidth = binwidth) +
        labs(title = glue("Number of mutations per sample ({nSample} samples)"),
             x = "Number of mutations",
             y = "Counts")
}

# selezione sui geni (n geni più o meno mutati)
select.genes.on.mutations <- function(df, n, desc = T){
    # traspone il dataset per operare sui geni
    temp.df <- as.data.frame(t(df))
    # ordinamento per somma
    temp.df$sums <- rowSums(temp.df[,1:ncol(temp.df)], na.rm = TRUE)
    if(desc){
        # ordinamento decrescente
        temp.df <- temp.df[order(-temp.df$sums),][1:n,] %>% select(-sums)
    }else{
        # ordinamento crescente
        temp.df <- temp.df[order(temp.df$sums),][1:n,] %>% select(-sums)
    }
    # ripristino struttura
    temp.df <- as.data.frame(t(temp.df))
    # rimozione degli na se n è troppo grande
    temp.df %>% select_if(function(x) any(!is.na(x))) %>% drop_na()
}

# selezione sui campioni (n campioni più o meno mutati)
select.samples.on.mutations <- function(df, n, desc = T){
    temp.df <- df
    # ordinamento per somma
    temp.df$sums <- rowSums(temp.df[,2:ncol(temp.df)], na.rm = TRUE)
    if(desc){
        temp.df <- temp.df[order(-temp.df$sums),][1:n,] %>% select(-sums)
    }else{
        temp.df <- temp.df[order(temp.df$sums),][1:n,] %>% select(-sums)
    }
    # rimozione degli na se n è troppo grande
    temp.df %>% select_if(function(x) any(!is.na(x))) %>% drop_na()
}

# preparazione del plot della correlazione dei geni
corrplot.from.df <- function(df){
    # prepare correlation matrix
    corr <- round(cor(df), 1)
    corr[is.nan(corr)] <- 0
    corr[is.na(corr)] <- 0
    corr[corr == Inf] <- 0
    # prepare corrplot
    ggcorrplot(corr, hc.order = TRUE)
}

corrplot.genes <- corrplot.from.df
corrplot.samples <- function (df)
    corrplot.from.df(as.data.frame(t(df)))



# compatta un dataframe unificando le righe identiche
# pone i conteggi nella colonna freq
compact.dataset.easy <- function(dataset){
    dataset %>% plyr::count( colnames(dataset) )
}

build.topology.subset <- function(samples){
    # relazione di "contenuto"
    edges = list()
    index = 1
    # un semplice ciclo che valuta la relazione contenuto per tutte le coppie di nodi (diversi)
    for(i in 1:nrow(samples)){
        for(j in 1:nrow(samples)){
            if(i!=j){
                r1 = samples[i,]
                r2 = samples[j,]
                if( ! ((-1) %in% (r1 - r2)) ) {
                    edges[[index]] <- c(j,i)
                    index <- index + 1
                }
            }
        }
    }
    edges
}

# prepara le etichette dei nodi tenendo conto dei soli geni mutati
# aggiunge anche il nodo clonale
prepare.labels <- function(samples, genes){
    # preparo le etichette con i geni relativi ai nodi (seleziono i geni mutati)
    labels = apply(samples, MARGIN = 1, FUN = function(x) genes[x==1] )
    # preparazione etichette (da vettori di stringhe)
    labels = lapply(labels, function (x) paste(x, collapse=", "))
    labels
}


# compute number of children for each node given an adj matrix
get.no.of.children <- function(A,g){
    no.of.children <- numeric(length(V(g)))
    for( v in V(g) ){
        no.of.children[v] = sum(A[v,])
    }
    no.of.children
}


