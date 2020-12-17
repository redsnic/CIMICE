library(rly)
library(purrr)

# parser and lexer to read therapies effect information
# implementation is based on the rly package

TOKENS <- c('NAME', 'NUMBER', 'KILL', 'MUTATE')
LITERALS <- c('&', '-', '|', '(', ')', '[', ']', ',', '{', '}', ';')

Lexer <- R6::R6Class(
    "Lexer",
    public = list(
        tokens = TOKENS,
        literals = LITERALS,
        t_KILL = "kill",
        t_MUTATE = "mutate",
        t_NAME = '[a-zA-Z_][a-zA-Z0-9]*',
        t_NUMBER = function(re='[0-9]+\\.[0-9]*|[0-9]+', t){
            t$value = as.numeric(t$value)
            # check if the number actually represents a probability
            if(t$value < 0 || t$value >1){
                stop(paste("ERROR at line", t$lexer$lineno, ":", t$value, "is not a probability"))
            }
            return(t)
        },
        # comments start with '#'
        t_comment = function(re='#[^\\n]*', t) {
            return(NULL)
        },
        t_ignore = " \t\n",
        # keep track of line number for errors
        t_newline = function(re='\\n+', t) {
            t$lexer$lineno <- t$lexer$lineno + nchar(t$value)
            return(NULL)
        },
        # lexical error
        t_error = function(t){
            stop(paste("ERROR at line", t$lexer$lineno, ": illegal character", t$value[1]))
        }
    )
)

Parser <- R6::R6Class(
    "Parser",
    public = list(
        tokens = TOKENS,
        literals = LITERALS,
        # Parsing rules
        precedence = list(c('left', ',', ';'),
                          c('left','|'),
                          c('left','&'),
                          c('right','UMINUS')),
        # dictionary of names
        p_init = function(
            doc='init : drugs', p){
            # return list of therapies
            p$set(1, p$get(2))
        },
        p_drug_general = function(
            doc='drugs : NAME \'{\' rules \'}\' drugs', p){
            # construct list of therapies
            # each element is a named list with therapy name and associated rules
            p$set(1, c(list(list(name=p$get(2), rules=p$get(4))), p$get(6)))
        },
        p_drug_single = function(
            doc='drugs : NAME \'{\' rules \'}\'', p){
            # construct list of therapies (single element)
            # each element is a named list with therapy name and associated rules
            p$set(1, list(list(name=p$get(2), rules=p$get(4))))
        },
        p_rules_multiple = function(
            doc = 'rules : rule \';\' rules', p){
            # construct rule (effect) list for a specific therapy
            p$set(1,c(p$get(2), p$get(4)))
        },
        p_rules_last = function(
            doc = 'rules : rule \';\'', p){
            # construct rule (effect) list for a specific therapy (single element)
            p$set(1,p$get(2))
        },
        p_rule_kill = function(
            doc = 'rule : guard KILL NUMBER', p){
            # construct rule: kill
            # kill rule is structured as follow
            # guard : function that given a label (a list of genes) states if
            #         this effect applies to the described genotype
            # action : the string "kill"
            # eff : probability that this
            p$set(1, list(list(guard = p$get(2), action="kill", eff=p$get(4))))
        },
        p_rule_mutate = function(
            doc = 'rule : guard MUTATE op NUMBER', p){
            # construct rule: mutate
            # mutate rule is structured as follow
            # guard : function that given a label (a list of genes) states if
            #         this effect applies to the described genotype
            # action : the string "mutate"
            # changes : list of mutations to add and remove (list(add=...,remove=...))
            # eff : probability that this
            p$set(1, list(list(guard = p$get(2), action="mutate", changes=p$get(4), eff=p$get(5))))
        },
        p_guard = function(
            doc = 'guard : \'[\' expr \']\'', p){
            # returns guard function, a function that given a label (a list of genes) states if
            #         this effect applies to the described genotype
            p$set(1, p$get(3))
        },
        p_expr_last = function(
            doc = "expr : NAME", p){
            # base case for guard expression
            # returns a function that  given a label (a list of genes) states if
            # it contains the give gene name
            p$set(1, find.gene(toupper(p$get(2))))
        },
        p_expr_brackets = function(
            doc = 'expr : \'(\' expr \')\'', p){
            # manage brackets
            p$set(1, p$get(2))
        },
        p_expr_not = function(
            doc = 'expr : \'-\' expr %prec UMINUS', p){
            # add a negation to guard function
            f <- function(op){
                force(op)
                function(x){
                    !op(x)
                }
            }
            p$set(1, f(p$get(2)))
        },
        p_expr = function(
            doc = 'expr : expr \'&\' expr
                        | expr \'|\' expr', p){
            # add logical AND to guard function
            if(p$get(3) == "&"){
                f <- function(op1,op2){
                    force(op1)
                    force(op2)
                    function(x){
                        op1(x) && op2(x)
                    }
                }
                p$set(1,
                      f(p$get(2), p$get(4))
                )
            # add logical OR to guard function
            }else if(p$get(3) == "|"){
                f <- function(op1,op2){
                    force(op1)
                    force(op2)
                    function(x){
                        op1(x) || op2(x)
                    }
                }
                p$set(1,
                      f(p$get(2), p$get(4))
                )
            # stub to add additional operator, is never reached
            }else{
                stop(paste("ERROR: unknown operator", p$get(3)))
            }

        },
        p_op = function(doc='op : \'[\' ops \']\'', p){
            # returns the list of genes on which to operate for "mutate" effects
            p$set(1,p$get(3))
        },
        p_ops_single = function(doc='ops : NAME', p){
            # base case of mutate effect list (add mutation on target gene)
            # the list is divided in genes to mutate (add) and to correct (remove)
            p$set(1,list(add=list(toupper(p$get(2))), remove=list()))
        },
        p_ops_uminus = function(doc='ops : \'-\' NAME %prec UMINUS', p){
            # base case of mutate effect list (remove mutation on target gene)
            # the list is divided in genes to mutate (add) and to correct (remove)
            p$set(1,list(add=list(), remove=list(toupper(p$get(3)))))
        },
        p_ops_multiple = function(doc='ops : ops \',\' ops', p){
            # extend mutate effect list
            # the list is divided in genes to mutate (add) and to correct (remove)
            p$set(1,
                  list(add=c(p$get(2)$add, p$get(4)$add),
                       remove=c(p$get(2)$remove, p$get(4)$remove)
                  ))
        },
        p_error = function(p) {
            # manage syntax errors
            if(is.null(p)) stop("ERROR: Syntax error at EOF")
            else           stop(paste("ERROR: Syntax error at", p$value))
        }
    )
)

lexer <- rly::lex(Lexer)
parser <- rly::yacc(Parser)

#' Prepare "find a gene" functions
#'
#' Creates a function to look for a specific gene in a list
#'
#' @param gene gene to be searched for
#'
#' @return a function that looks for a specific gene (a string) in a list
#'         that takes that list in input. The output of that function is Boolean.
#'
#' @examples
#' require(purrr)
#' gene <- "A"
#' label <- c("B","C","A")
#' (find.gene(gene))(label)
#'
#' @export find.gene
find.gene <- function(gene){
    force(gene)
    function(label){
        any(map_lgl(label, ~ . == gene))
    }
}

#' Lex and parse treatment information
#'
#' Creates a data structure as defined by 'parser' object that represents
#' the information included in a treatment description file.
#' Check vignettes for format.
#'
#' @param path path to treatment description file
#'
#' @return a data structure with treatment information
#'
#' @examples
#' load.treatments('./inst/extdata/test.treatments')
#'
#' @export load.treatments
load.treatments <- function(path){
    txt <- paste(readLines(path), collapse="\n")
    parser$parse(txt, lexer)
}



