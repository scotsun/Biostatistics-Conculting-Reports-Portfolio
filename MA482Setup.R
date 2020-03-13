################################################################################
# File: MA482Setup
# Course: MA482 Biostatistics
# Description: Additional R functions for creating MA482 styles.
#
# Author: Eric Reyes
# Date: Spring 2019

## ---- Load Packages ----
pkgs <- c("boot", 
          "car",
          "skimr",
          "splines",
          "geepack",
          "rms",
          "lme4",
          "tidyverse",
          "broom",
          "broom.mixed",
          "knitr")

for(pkg in pkgs) library(pkg, character.only = TRUE, quietly = TRUE)


## ---- Change Options ----
# Suppress status bar in dplyr.
options(dplyr.show_progress = FALSE,
        contrasts = rep("contr.treatment", 2))
skim <- skim_with(numeric = sfl(hist = NULL), 
                  ts = sfl(line_graph = NULL),
                  character = get_skimmers()$factor)

# Change theme for plots
theme_set(theme_bw(12))
theme_update(legend.position = "bottom",
             legend.box = "vertical")


# Specify chunk options


knitr::opts_chunk$set(
  prompt = FALSE,
  comment = "",
  linewidth = 80)

if (is_latex_output()){
  knitr::opts_chunk$set(out.width = "0.8\\textwidth",
                        fig.align = "center")
} else if (is_html_output()){
  knitr::opts_chunk$set(out.width = "80%",
                        fig.align = "center")
}


## ---- Ensure Source Code Wraps ----
.hook_source = knitr::knit_hooks$get("source")

knitr::knit_hooks$set(
  source = function(x, options){
    # this hook is used only when linewidth option is not NULL
    if (!is.null(n <- options$linewidth)){
      # # only split on real line breaks
      # breakinstr <- '(\"([^\"]*\n)+[^\"]+\")|(\'([^\']*\n)+[^\']+\')'
      # if (stringr::str_detect(x, breakinstr)){
      #   t1 <- unlist(stringr::str_extract_all(x, breakinstr))
      #   t2 <- stringr::str_replace_all(t1, "\n", " ")
      #   
      #   for (i in 1:length(t1)){
      #     x <- stringr::str_replace(x, t1[i], t2[i])
      #   }
      # }
      # 
      # x = knitr:::split_lines(x)
      # x.length <- nchar(x)
      # spaces <- stringr::str_extract(x, "^[[:space:]]+")
      # spaces[is.na(spaces)] <- ""
      # spaces.length <- nchar(spaces)
      # 
      # # if all lines fine, skip
      # if (any(x.length > n)){
      #   # wrap lines which are longer than n characters
      #   x <- mapply(function(u, v, w, y){
      #     if (v > n){
      #       u = paste(w, stringr::str_wrap(u, width = n - y))
      #     } else {
      #       u
      #     }}, x, x.length, spaces, spaces.length)
      #   
      #   x <- paste(x, collapse = "\n")
      # }
      
      x = knitr:::split_lines(x)
      
      x = ifelse(nchar(x) > n, stringr::str_wrap(x, width = n, exdent = 2), x)
    }
    
    .hook_source(x, options)
  })



## ---- Create Special Blocks ----
eng_instructor <- function(options) {
  if (identical(options$echo, FALSE)) return()
  
  # Steal some ideas from block definition
  to = opts_knit$get("rmarkdown.pandoc.to")
  is_pandoc = !is.null(to)
  if(!is_pandoc){
    to = opts_knit$get("out.format")
    if(!(to %in% c("latex", "html", "markdown"))) to = NULL
  }
  
  if(is.null(to)) return(code)
  if(to=="beamer") to = "latex"
  if(grepl("(markdown)|(epub)|(html)|(revealjs)|(s5)|(slideous)|(slidy)",
           to)) to = "html"
  
  
  code = paste(options$code, collapse = '\n'); type = options$type
  if (is.null(type)) return(code)
  
  if(!is.null(type) && type=="solution"){
    code = paste("__SOLUTION__:", code)
  }
  
  if (is.null(opts_knit$get("rmarkdown.pandoc.to"))) stop('The engine "block2" is for R Markdown only')
  
  l1 = options$latex.options
  if (is.null(l1)) l1 = ''
  # protect environment options because Pandoc may escape the characters like
  # {}; when encoded in integers, they won't be escaped, but will need to
  # restore them later; see bookdown:::restore_block2
  if (l1 != '') l1 = paste(
    c('\\iffalse{', utf8ToInt(enc2utf8(l1)), '}\\fi{}'), collapse = '-'
  )
  h2 = ifelse(is.null(options$html.tag), 'div', options$html.tag)
  h3 = ifelse(is.null(options$html.before), '', options$html.before)
  h4 = ifelse(is.null(options$html.after), '', options$html.after)
  h5 = ifelse(is.null(options$html.before2), '', options$html.before2)
  h6 = ifelse(is.null(options$html.after2), '', options$html.after2)
  
  if(to=="latex"){
    sprintf('\\BeginKnitrBlock{%s}%s\n%s%s%s\n\\EndKnitrBlock{%s}',
            type, l1, h5, code, h6, type)
  } else {
    sprintf(
      '\\BeginKnitrBlock{%s}%s%s<%s class="%s" custom-style="%s">%s%s%s</%s>%s\\EndKnitrBlock{%s}',
      type, l1, h3, h2, type, type, h5, code, h6, h2, h4, type
    )
  }
}


knit_engines$set(c(knit_engines$get(), 
                   "instructor" = eng_instructor))


## ---- Special Functions ----

################################################################################
# function: confint.matrix
# description: method for confint that takes in a K-matrix to produce 
#              confidence intervals for linear combinations of the parameters.
confint.matrix <- function(object, parm, level = 0.95, coef., vcov., ...){
  Kb <- drop(object %*% coef.)
  KSK <- object %*% vcov. %*% t(object)
  
  parm <- seq_along(Kb)
  
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3))
  
  fac <- qnorm(a)
  
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(KSK))[parm]
  ci[] <- Kb[parm] + ses %o% fac
  
  ci 
}



################################################################################
# function: CoxSnellPlot
# description: Create a plot of the Cox-Snell residuals for examining the
#              goodness of fit for a cox PH model.
# 
# parameters
#   fit        coxph Model fit.  Object fit using coxph() in survival package.
CoxSnellPlot <- function(fit){
  # Update model, if necessary
  fit <- update(fit, x=TRUE, y=TRUE)
  
  # Obtain Cox-Snell residuals from fit
  resi <- fit$y[,2] - residuals(fit,type="martingale") 
  
  # Compute the cumulative hazard of these residuals
  Sfit <- survfit(Surv(time=resi, event=fit$y[,2])~1)
  CumHaz <- cumsum(Sfit$n.event/Sfit$n.risk)
  
  plot.dat <- data.frame(resi=Sfit$time,
                         cumhaz=CumHaz)
  
  # Construct cumulative hazard vs. residuals
  p1 <- ggplot(data=plot.dat,mapping=aes(y=cumhaz,x=resi)) +
    geom_abline(slope=1,intercept=0,colour="grey") +
    geom_point() +
    labs(x="Cox-Snell Residuals",y="Estimated Cumulative Hazard Rates") +
    theme_bw(16)
  
  return(p1)
}



################################################################################
# glance.geeglm
# description: Method for glance for the geeglm class.
glance.geeglm <- function(x, ...){
  s <- summary(x)
  data.frame(n=sum(s$clusz),
             sigma=sqrt(s$geese$scale$estimate),
             df.residual=s$df[2],
             n.clusters=length(s$clusz),
             max.cluster.size=max(s$clusz))
}



################################################################################
# KMest
# description: Compute kaplan-meier estimates for survival data. This is a front
#              end to survfit in order to return a nicer set of data.
#
# parameters:
#   form       Formula for describing the survival data.
#   data       Optional dataset where variables appearing in form should be
#                located
#   conflevel  Scalar. Confidence level for confidence intervals (default=0.95).
KMest <- function(form, data, conflevel=0.95){
  # Use survfit to get the estimates
  if(missing(data)) fit <- survfit(formula=form, conf.int=conflevel)
  if(!missing(data)){
    fit <- survfit(formula=form, data=data, conf.int=conflevel)
  } 
  
  # Obtain name of time in dataset
  survobj <- strsplit(as.character(form[2]), ",")[[1]][1]
  survobj <- strsplit(survobj,"\\(")[[1]][2]
  survobj <- strsplit(survobj,"=")[[1]][1]
  
  # Construct data
  dat <- data.frame(Survival=fit$surv,
                    UCL=fit$upper,
                    LCL=fit$lower)
  
  dat[,survobj] <- fit$time
  dat <- dat[,c(4,1,2,3)]
  
  # Add grouping if available
  if(!is.null(fit$strata)){
    names(fit$strata) <- sapply(strsplit(names(fit$strata),"="),"[",2)
    dat[,as.character(form[3])] <- rep(names(fit$strata), times=fit$strata)
  }
  
  return(dat)
}




################################################################################
# function: LifeTable
# description: Construct life table estimates for survival data.
#
# parameters:
#   nrisk      Number at risk at start of interval.
#   nevent     Number of events occurring during interval.
#   ncensor    Number of subjects censored during interval.
#   data       Dataframe. Optional dataset containing the variables above.
#   conflevel  Scalar. Number between 0 and 1 defining the confidence level.
LifeTable <- function(nrisk, nevent, ncensor, data, conflevel=0.95){
  data <- data.frame(data)
  
  # Form Dataset
  if(!missing(data)){
    dat <- data.frame(nrisk=data[,as.character(substitute(nrisk))],
                      nevent=data[,as.character(substitute(nevent))],
                      ncensor=data[,as.character(substitute(ncensor))])
  }
  
  if(missing(data)){
    dat <- data.frame(nrisk=nrisk,
                      nevent=nevent,
                      ncensor=ncensor)
  }
  
  
  # Compute Quantities of Interest
  dat <- within(dat, {
    # Mortality Rate
    mortality <- nevent/(nrisk - 0.5*ncensor)
    
    # Survival
    survival <- cumprod(1-mortality)
    
    # SE
    se <- survival*sqrt(cumsum(mortality/(nevent*(1-mortality))))
    
    # Confidence Limits
    LCL <- pmax(0, survival - qnorm(0.5*(1+conflevel))*se)
    UCL <- pmin(1, survival + qnorm(0.5*(1+conflevel))*se)
  })
  
  dat <- dat[,c(1,2,3,8,7,6,5,4)]
  names(dat) <- c("Number at Risk","Number of Events","Number Censored",
                  "Mortality", "Survival", "Std. Error",
                  "Lower CL", "Upper CL")
  
  class(dat) <- c("LifeTable","data.frame")
  return(dat)
}



################################################################################
# function: print.LifeTable
# description: Method for printing results of a LifeTable object.  This will not
#              be used by students directly.
#
# parameters:
#  x           LifeTable object.  The result of a LifeTable call.
#  ...         Other parameters to be passed to function. Currently ignored.
print.LifeTable <- function(x, ...){
  x[,-c(1:3)] <- round(x[,-c(1:3)], 3)
  print.data.frame(x, row.names=FALSE)
}



################################################################################
# fct: stat_halfnorm
# description: Make a stat for ggplot2 that makes the halfnormal plot.
StatHalfnorm <- ggproto("StatHalfnorm", Stat,
                        compute_group = function(data, scales, na.rm=FALSE){
                          # sort data
                          sample <- sort(data$sample)
                          
                          # obtain length
                          n <- length(sample)
                          
                          # determine quanitles
                          halfnorm <- qnorm((n + seq(n))/(2*n + 1))
                          
                          data.frame(sample, halfnorm)
                        },
                        required_aes = "sample",
                        default_aes = aes(x=..halfnorm.., y=..sample..))

stat_halfnorm <- function(mapping=NULL, data=NULL, geom="point",
                          position="identity", na.rm=FALSE, show.legend=NA,
                          inherit.aes=TRUE, ...){
  ggplot2::layer(data=data, mapping=mapping, stat=StatHalfnorm, geom=geom,
                 position=position, show.legend=show.legend,
                 inherit.aes=inherit.aes, params=list(na.rm=na.rm, ...))
}



################################################################################
# fct: stat_qqline
# description: Make a stat for ggplot2 that adds a best-fit line to the 
#              quantile plot to better see the linearity.
StatQqline <- ggproto("StatQqline", Stat,
                      compute_group = function(data, scales, quantiles=NULL, 
                                               distribution=stats::qnorm, dparams=list(), 
                                               na.rm=FALSE){
                        # replicate that from StatQq
                        sample <- sort(data$sample)
                        n <- length(sample)
                        if(is.null(quantiles)){
                          quantiles <- stats::ppoints(n)
                        } else {
                          stopifnot(length(quantiles)==n)
                        }
                        
                        theoretical <- do.call(distribution, c(list(p=quote(quantiles)), dparams))
                        
                        # compute best fit
                        fit <- lm(sample ~ theoretical)
                        intercept <- coef(fit)[1]
                        slope <- coef(fit)[2]
                        
                        data.frame(intercept, slope)
                      },
                      required_aes = "sample")

stat_qqline <- function(mapping=NULL, data=NULL, geom="abline",
                        position="identity", distribution=stats::qnorm,
                        dparams=list(), na.rm=FALSE, show.legend=NA,
                        inherit.aes=TRUE, ...){
  ggplot2::layer(data=data, mapping=mapping, stat=StatQqline, geom=geom,
                 position=position, show.legend=show.legend, 
                 inherit.aes=inherit.aes, 
                 params=list(distribution=distribution, dparams=dparams, 
                             na.rm=na.rm, ...))
}



################################################################################
# function: vcov.geeglm
# description: Method for obtaining the variance covariance matrix from a
#              geeglm object.
#
# parameters:
#   x          geeglm Object. The result of a gee fit to a glm model using the
#               geepack algorithms.
#   ...
vcov.geeglm <- function(x,...){
  return(x$geese$vbeta)
}




################################################################################
# function: WildBoot
# description: Methods for implementing a wild bootstrap for lm and nls objects.
#              This is a very rough approach based on the Boot() function.

WildBoot <- function(object, f = coef, labels = names(f(object)), R = 999,
                     ncores = 1, ...){
  UseMethod("WildBoot")
}

WildBoot.nls <- function(object, f = coef, labels = names(f(object)), R = 999,
                         ncores = 1, ...){
  f0 <- f(object)
  all.names <- all.vars(object$m$formula())
  param.names <- names(object$m$getPars())
  vars <- all.names[!(all.names %in% param.names)]
  obj <- update(object, data = na.omit(eval(object$data)[, 
                                                         vars]), start = coef(object))
  
  if (!(requireNamespace("car")))
    stop("The 'car' package is missing")
  if (!(requireNamespace("boot"))) 
    stop("The 'boot' package is missing")
  if (length(labels) != length(f0)) 
    labels <- paste("V", seq(length(f0)), sep = "")
  opt <- options(show.error.messages = FALSE)
  
  # adjustment for wild boostrap patterned after residual bootstrap
  boot.f <- function(data, indices, .fn){
    res <- residuals(object)
    
    .n <- length(res)
    .delta1 <- sqrt((3/4) + (sqrt(17)/12))
    .delta2 <- sqrt((3/4) - (sqrt(17)/12))
    .v <- matrix(rnorm(2*.n), nrow=.n, ncol=2)
    
    .mammen <- (.delta1 + (.v[,1]/sqrt(2)))*(.delta2 + (.v[,2]/sqrt(2))) -
      (.delta1*.delta2)
    
    val <- fitted(object) + res*.mammen
    if (!is.null(object$na.action)){
      pad <- object$na.action
      attr(pad, "class") <- "exclude"
      val <- naresid(pad, val)
    }
    
    assign(".y.boot", val, envir = .carEnv)
    mod <- try(update(object, get(".y.boot", envir = .carEnv) ~ 
                        ., start = coef(object)))
    if (inherits(mod, "try-error")) {
      out <- .fn(object)
      out <- rep(NA, length(out))
    }
    else {
      out <- .fn(mod)
    }
    out
  }
  
  if (ncores <= 1) {
    parallel_env = "no"
    ncores = getOption("boot.ncpus", 1L)
  }
  else {
    if (.Platform$OS.type == "unix") {
      parallel_env = "multicore"
    }
    else {
      parallel_env = "snow"
    }
  }
  
  b <- boot::boot(data.frame(update(object, model = TRUE)$model), 
                  boot.f, R, .fn = f, parallel = parallel_env, ncpus = ncores, 
                  ...)
  colnames(b$t) <- labels
  if (exists(".y.boot", envir = .carEnv)) 
    remove(".y.boot", envir = .carEnv)
  if (exists(".boot.indices", envir = .carEnv)) 
    remove(".boot.indices", envir = .carEnv)
  options(opt)
  d <- dim(na.omit(b$t))[1]
  if (d != R) 
    cat(paste("\n", "Number of bootstraps was", d, "out of", 
              R, "attempted", "\n"))
  b
}


WildBoot.lm <- function (object, f = coef, labels = names(f(object)), R = 999, 
                         ncores = 1, ...){
  if (!is.null(object$na.action)) 
    stop("The WildBoot function does not currently allow\n  missing values for lm or glm models.  Refit your model with rows \n  with missing values removed.  If you have a data frame called 'd', \n  then the argument data=na.omit(d) is likely to work.")
  
  if (!(requireNamespace("car")))
    stop("The 'car' package is missing")
  if (!(requireNamespace("boot"))) 
    stop("The 'boot' package is missing")
  f0 <- f(object)
  if (length(labels) != length(f0)) 
    labels <- paste0("V", seq_along(f0))
  
  # make alterations for wild bootstrap based on residual bootstrap in Boot.default
  boot.f <- function(data, indices, .fn){
    res <- residuals(object, type = "pearson")
    
    .n <- length(res)
    .delta1 <- sqrt((3/4) + (sqrt(17)/12))
    .delta2 <- sqrt((3/4) - (sqrt(17)/12))
    .v <- matrix(rnorm(2*.n), nrow=.n, ncol=2)
    
    .mammen <- (.delta1 + (.v[,1]/sqrt(2)))*(.delta2 + (.v[,2]/sqrt(2))) -
      (.delta1*.delta2)
    
    val <- fitted(object) + res*.mammen
    
    if (!is.null(object$na.action)) {
      pad <- object$na.action
      attr(pad, "class") <- "exclude"
      val <- naresid(pad, val)
    }
    assign(".y.boot", val, envir = .carEnv)
    mod <- update(object, get(".y.boot", envir = .carEnv) ~ .)
    out <- if (!is.null(object$qr) && (mod$qr$rank != object$qr$rank))
      f0 * NA
    else .fn(mod)
    out
  }
  
  nobs0 <- function(x, ...) {
    rval <- try(stats::nobs(x, ...), silent = TRUE)
    if (inherits(rval, "try-error") | is.null(rval)) 
      rval <- NROW(residuals(x, ...))
    return(rval)
  }
  n <- nobs0(object)
  dd <- data.frame(.zero = rep.int(0L, n))
  if (ncores <= 1) {
    parallel_env = "no"
    ncores = getOption("boot.ncpus", 1L)
  }
  else {
    if (.Platform$OS.type == "unix") {
      parallel_env = "multicore"
    }
    else {
      parallel_env = "snow"
    }
  }
  b <- boot::boot(dd, boot.f, R, .fn = f, parallel = parallel_env, 
                  ncpus = ncores, ...)
  colnames(b$t) <- labels
  if (exists(".y.boot", envir = .carEnv)) 
    remove(".y.boot", envir = .carEnv)
  if (exists(".boot.indices", envir = .carEnv)) 
    remove(".boot.indices", envir = .carEnv)
  b
}
