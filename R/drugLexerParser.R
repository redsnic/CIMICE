library(rly)
library(purrr)

TOKENS <- c('NAME', 'NUMBER', 'KILL', 'MUTATE')
LITERALS <- c('&', '-', '|', '(', ')', '[', ']', ',', '{', '}', ';')

find.gene <- function(gene){
    force(gene)
    function(label){
        any(map_lgl(label, ~ . == gene))
    }
}

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
            if(t$value < 0 || t$value >1){
                stop(paste("ERROR at line", t$lexer$lineno, ":", t$value, "is not a probability"))
            }
            return(t)
        },
        t_comment = function(re='#[^\\n]*', t) {
            return(NULL)
        },
        t_ignore = " \t\n",
        t_newline = function(re='\\n+', t) {
            t$lexer$lineno <- t$lexer$lineno + nchar(t$value)
            return(NULL)
        },
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
            p$set(1, p$get(2))
        },
        p_drug_general = function(
            doc='drugs : NAME \'{\' rules \'}\' drugs', p){
            p$set(1, c(list(list(name=p$get(2), rules=p$get(4))), p$get(6)))
        },
        p_drug_single = function(
            doc='drugs : NAME \'{\' rules \'}\'', p){
            p$set(1, list(list(name=p$get(2), rules=p$get(4))))
        },
        p_rules_multiple = function(
            doc = 'rules : rule \';\' rules', p){
            p$set(1,c(p$get(2), p$get(4)))
        },
        p_rules_last = function(
            doc = 'rules : rule \';\'', p){
            p$set(1,p$get(2))
        },
        p_rule_kill = function(
            doc = 'rule : guard KILL NUMBER', p){
            p$set(1, list(list(guard = p$get(2), action="kill", eff=p$get(4))))
        },
        p_rule_mutate = function(
            doc = 'rule : guard MUTATE op NUMBER', p){
            p$set(1, list(list(guard = p$get(2), action="mutate", changes=p$get(4), eff=p$get(5))))
        },
        p_guard = function(
            doc = 'guard : \'[\' expr \']\'', p){
            p$set(1, p$get(3))
        },
        p_expr_last = function(
            doc = "expr : NAME", p){
            p$set(1, find.gene(toupper(p$get(2))))
        },
        p_expr_brackets = function(
            doc = 'expr : \'(\' expr \')\'', p){
            p$set(1, p$get(2))
        },
        p_expr_not = function(
            doc = 'expr : \'-\' expr %prec UMINUS', p){
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
            }else{
                stop(paste("ERROR: unknown operator", p$get(3)))
            }

        },
        p_op = function(doc='op : \'[\' ops \']\'', p){
            p$set(1,p$get(3))
        },
        p_ops_single = function(doc='ops : NAME', p){
            p$set(1,list(add=list(toupper(p$get(2))), remove=list()))
        },
        p_ops_uminus = function(doc='ops : \'-\' NAME %prec UMINUS', p){
            p$set(1,list(add=list(), remove=list(toupper(p$get(3)))))
        },
        p_ops_multiple = function(doc='ops : ops \',\' ops', p){
            p$set(1,
                  list(add=c(p$get(2)$add, p$get(4)$add),
                       remove=c(p$get(2)$remove, p$get(4)$remove)
                  ))
        },
        p_error = function(p) {
            if(is.null(p)) stop("ERROR: Syntax error at EOF")
            else           stop(paste("ERROR: Syntax error at", p$value))
        }
    )
)

lexer <- rly::lex(Lexer)
parser <- rly::yacc(Parser)

load.treatments <- function(path){
    txt <- paste(readLines(path), collapse="\n")
    parser$parse(txt, lexer)
}



