## Function takes a formula, Formula, or a list of formulas
## and extracts only the necessary parts to create a model.frame,
## as well as the parts that are needed to setup smooth term
## specification lists, or any type of special model terms
## that should be used for fitting.
fake_formula <- function(formula, specials = NULL, nospecials = FALSE, onlyspecials = FALSE)
{
  if(is.list(formula))
    formula <- do.call("as.Formula", formula)

  if(any(grepl("|", as.character(formula), fixed = TRUE)))
    formula <- as.Formula(formula)

  n <- length(formula)

  outer <- FALSE
  if(length(n) > 1) {
    if(n[2L] > 1)
      outer <- TRUE
  }

  if(outer) {
    fl <- list()
    lhs <- formula(formula, rhs = 0)
    for(i in 1:n[2L]) {
      fl[[i]] <- fake_formula(formula(formula, rhs = i, lhs = 0), specials = specials,
        nospecials = nospecials, onlyspecials = onlyspecials)
    }
    if(length(fl)) {
      if(!is.character(fl[[1L]])) {
        formula <- . ~ .
        fl <- do.call("as.Formula", fl)
        fl <- formula(fl, lhs = 0, drop = FALSE)
        formula <- as.Formula(fl)
        formula[[3L]] <- fl
        attr(formula, "lhs") <- attr(as.Formula(lhs), "lhs")
      } else {
        formula <- fl
      }
    }
  } else {
    stn <- c("s", "te", "t2", "sx", "s2", "rs", "ti",
      "tx", "tx2", "tx3", "tx4", "la", "n", "lin",
      "pb", "pbc", "nn", "fk", "re", "ps", "pbz", "ga",
      "random", "ra", "lo", "tr", "tree", "cf", "NN", "pb2", "ct", "st")
    stn <- unique(c(stn, specials))

    mt <- terms(formula, specials = stn)

    os <- NULL
    st <- unlist(attr(mt, "specials"))
    if(!is.null(st)) {
      tls <- rownames(attr(mt, "factors"))[st]
      tls <- unique(tls)
      ff <- NULL
      for(j in tls) {
        p <- parse(text = j)
        v <- all.vars(p)
        e <- eval(parse(text = paste0("quote(", j, ")")))
        for(i in seq_along(e)) {
          av <- all.vars(e[[i]])
          av <- av[av != "|"]
          if(length(av)) {
            if(all(av %in% v)) {
              ef <- try(eval(e[[i]]), silent = TRUE)
              if(!inherits(ef, "try-error")) {
                if(inherits(ef, "formula")) {
                  vf <- attr(terms(eval(ef)), "variables")
                  for(k in 2:length(vf))
                    ff <- c(ff, vf[[k]])
                  next
                }
              }
              if(as.character(e[[i]])[1L] == "~") {
                vf <- attr(terms(eval(e[[i]])), "variables")
                for(k in 2:length(vf))
                  ff <- c(ff, vf[[k]])
              } else {
                ff <- c(ff, e[[i]])
              }
            }
          }
        }
        os <- c(os, j)
        eval(parse(text = paste0("formula <- update(formula, . ~ . -", j,")")))
        if(!nospecials) {
          for(i in ff) {
            eval(parse(text = paste0("formula <- update(formula, . ~ . +", deparse(i),")")))
          }
        }
      }
    }

    if(onlyspecials) {
      if(!is.character(os))
        os <- character(0)
      os <- gsub(" ", "", os)
      formula <- os
    } else {
      tf <- terms(formula, specials = stn)
      sj <- unlist(attr(tf, "specials"))
      if(!is.null(sj)) {
        formula <- fake_formula(formula)
      }
      tf <- terms(formula, specials = stn)
      if(length(j <- grep("list(", attr(tf, "term.labels"), fixed = TRUE, value = TRUE))) {
        fc <- paste("formula <- update(formula, . ~ . -", j, ")")
        eval(parse(text = fc))
      }
    }
  }

  return(formula)
}

