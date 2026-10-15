## InvaPact country-level analyses

## libraries
if (!requireNamespace("invapact", quietly = TRUE)) {
  invapact_package <- path.expand(
    "~/Documents/Papers/Invasive species/InvaPact/invapact package/invapact_0.2.0.tar.gz"
  )
  if (!file.exists(invapact_package)) {
    stop("InvaPact package archive not found: ", invapact_package, call. = FALSE)
  }
  install.packages(invapact_package, repos = NULL, type = "source")
}
library(invapact)
library(boot)
library(ggplot2)
library(dismo)
library(gbm)
library(readxl)
library(dplyr)
library(countrycode)


## functions
modifyVecFunc <- function(obj_name, index, new_value) {
  tryCatch({
    if (!object_exists(obj_name)) {
      stop("object does not exist: ", obj_name)
    }
    
    obj <- get(obj_name)
    
    # check object type and handle accordingly
    if (is.vector(obj) && !is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[index] <- new_value
    } else if (is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[[index]] <- new_value
    } else if (is.data.frame(obj)) {
      if (!all(index <= dim(obj))) stop("index out of bounds")
      obj[index[1], index[2]] <- new_value
    } else {
      stop("unsupported object type")
    }
    
    # Save modified object back to its original name
    assign(obj_name, obj, envir = .GlobalEnv)
    
  }, error = function(e) {
    message("error: ", e$message)
    return(FALSE)
  })
}

object_exists <- function(obj_name) {
  exists(obj_name, envir = .GlobalEnv)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

gbm.step.base <- get("gbm.step", envir = asNamespace("dismo"))

gbm.step.patched <- local({
  gbm.step.fn <- gbm.step.base
  body_lines <- deparse(body(gbm.step.fn), width.cutoff = 500)
  body_lines <- gsub("verbose = FALSE\\)",
                     "verbose = FALSE, n.minobsinnode = n.minobsinnode)",
                     body_lines)
  formals(gbm.step.fn)$n.minobsinnode <- 10L
  body(gbm.step.fn) <- parse(text = paste(body_lines, collapse = "\n"))
  gbm.step.fn
})

fit_cv_correlation <- function(fit) {
  if (is.null(fit) || is.null(fit$cv.statistics)) {
    return(-Inf)
  }
  cor.value <- fit$cv.statistics$correlation.mean
  if (!is.finite(cor.value)) {
    return(-Inf)
  }
  cor.value
}

fit_cv_correlation_se <- function(fit) {
  if (is.null(fit) || is.null(fit$cv.statistics)) {
    return(Inf)
  }
  cor.se <- fit$cv.statistics$correlation.se
  if (!is.finite(cor.se)) {
    return(Inf)
  }
  cor.se
}

fit_cv_jitter <- function(fit) {
  if (is.null(fit) || is.null(fit$cv.values)) {
    return(Inf)
  }
  cv.values <- fit$cv.values
  cv.values <- cv.values[is.finite(cv.values)]
  if (length(cv.values) < 7) {
    return(Inf)
  }
  
  global.k <- min(7L, length(cv.values) - (1L - length(cv.values) %% 2L))
  if (global.k < 3) {
    return(Inf)
  }
  global.smooth <- stats::runmed(cv.values, k = global.k, endrule = "median")
  scale.value <- diff(range(cv.values, na.rm = TRUE))
  if (!is.finite(scale.value) || scale.value == 0) {
    scale.value <- stats::sd(cv.values, na.rm = TRUE)
  }
  if (!is.finite(scale.value) || scale.value == 0) {
    scale.value <- 1
  }
  
  global.mad <- stats::mad(cv.values - global.smooth, center = 0, constant = 1,
                           na.rm = TRUE) / scale.value
  min.idx <- which.min(cv.values)[1]
  local.idx <- seq.int(max(1L, min.idx - 4L), min(length(cv.values), min.idx + 4L))
  local.values <- cv.values[local.idx]
  local.k <- min(5L, length(local.values) - (1L - length(local.values) %% 2L))
  local.mad <- 0
  if (local.k >= 3) {
    local.smooth <- stats::runmed(local.values, k = local.k, endrule = "median")
    local.mad <- stats::mad(local.values - local.smooth, center = 0, constant = 1,
                            na.rm = TRUE) / scale.value
  }
  
  cv.diffs <- diff(cv.values)
  cv.diffs <- cv.diffs[is.finite(cv.diffs) & abs(cv.diffs) > 0]
  sign.change.rate <- 0
  if (length(cv.diffs) > 1) {
    sign.change.rate <- sum(diff(sign(cv.diffs)) != 0, na.rm = TRUE) / length(cv.diffs)
  }
  
  max(global.mad, local.mad, sign.change.rate)
}

fit_has_jittery_convergence <- function(fit, jitter.threshold = 0.12) {
  isTRUE(resolve_fit_jitter(fit) > jitter.threshold)
}

resolve_fit_jitter <- function(fit, params = NULL) {
  jitter <- NULL
  if (!is.null(params)) {
    jitter <- params$jitter
  }
  jitter <- jitter[is.finite(jitter)]
  if (length(jitter) == 0L) {
    jitter <- fit_cv_jitter(fit)
  }
  jitter <- jitter[is.finite(jitter)]
  if (length(jitter) == 0L) {
    return(Inf)
  }
  jitter[[1]]
}

get_gbm_response <- function(dat, gbm.y) {
  if (length(gbm.y) == 1 && is.numeric(gbm.y)) {
    return(dat[[gbm.y]])
  }
  dat[[gbm.y]]
}

get_gbm_label <- function(dat, gbm.y) {
  if (length(gbm.y) == 1 && is.numeric(gbm.y)) {
    return(colnames(dat)[gbm.y])
  }
  paste(gbm.y, collapse = "_")
}

validate_gbm_columns <- function(dat, columns, arg.name) {
  if (is.null(columns)) {
    stop(arg.name, " must not be NULL", call. = FALSE)
  }
  if (is.numeric(columns)) {
    invalid <- columns[!is.finite(columns) | columns < 1 | columns > ncol(dat)]
    if (length(invalid) > 0) {
      stop(arg.name, " contains invalid column indices: ",
           paste(unique(invalid), collapse = ", "),
           call. = FALSE)
    }
    return(invisible(colnames(dat)[columns]))
  }
  missing.cols <- setdiff(columns, colnames(dat))
  if (length(missing.cols) > 0) {
    stop(arg.name, " contains columns not present in data: ",
         paste(missing.cols, collapse = ", "),
         call. = FALSE)
  }
  invisible(columns)
}

choose_gbm_n_folds <- function(dat, gbm.y, family = "gaussian", requested.n.folds = 10) {
  safe.n.folds <- min(as.integer(requested.n.folds), max(2L, nrow(dat) - 1L))
  if (identical(family, "bernoulli")) {
    response.values <- get_gbm_response(dat, gbm.y)
    class.counts <- table(response.values)
    if (length(class.counts) > 1) {
      safe.n.folds <- min(safe.n.folds, as.integer(min(class.counts)))
    }
  }
  max(2L, safe.n.folds)
}

fit_is_weak <- function(fit, weak.cor.threshold = 0.15, max.trees = 100000,
                        jitter.threshold = 0.12) {
  if (is.null(fit)) {
    return(TRUE)
  }
  if (fit_has_jittery_convergence(fit, jitter.threshold = jitter.threshold)) {
    return(TRUE)
  }
  if (fit_cv_correlation(fit) < weak.cor.threshold) {
    return(TRUE)
  }
  fit$n.trees >= floor(max.trees * 0.95)
}

fit_is_better <- function(candidate.fit, candidate.params, current.fit, current.params) {
  if (is.null(candidate.fit)) {
    return(FALSE)
  }
  if (is.null(current.fit)) {
    return(TRUE)
  }
  
  candidate.acceptable <- isTRUE(candidate.params$acceptable)
  current.acceptable <- isTRUE(current.params$acceptable)
  if (candidate.acceptable != current.acceptable) {
    return(candidate.acceptable)
  }
  
  candidate.cor <- fit_cv_correlation(candidate.fit)
  current.cor <- fit_cv_correlation(current.fit)
  if (abs(candidate.cor - current.cor) > 0.01) {
    return(candidate.cor > current.cor)
  }
  
  candidate.jitter <- resolve_fit_jitter(candidate.fit, candidate.params)
  current.jitter <- resolve_fit_jitter(current.fit, current.params)
  jitter.diff <- candidate.jitter - current.jitter
  if (is.finite(jitter.diff) && abs(jitter.diff) > 0.02) {
    return(candidate.jitter < current.jitter)
  }
  
  candidate.se <- fit_cv_correlation_se(candidate.fit)
  current.se <- fit_cv_correlation_se(current.fit)
  if (!isTRUE(all.equal(candidate.se, current.se))) {
    return(candidate.se < current.se)
  }
  
  candidate.distance <- abs(candidate.fit$n.trees - candidate.params$target.min)
  current.distance <- abs(current.fit$n.trees - current.params$target.min)
  candidate.distance < current.distance
}

build_bag_candidates <- function(requested.bag.fraction = NULL, n) {
  base.bag <- if (n < 75) {
    0.85
  } else if (n < 150) {
    0.75
  } else if (n < 300) {
    0.65
  } else {
    0.55
  }
  bag.candidates <- unique(round(c(
    requested.bag.fraction,
    base.bag,
    base.bag + 0.1,
    base.bag - 0.1
  ), 2))
  bag.candidates <- bag.candidates[is.finite(bag.candidates)]
  bag.candidates[bag.candidates >= 0.5 & bag.candidates <= 0.9]
}

build_lr_candidates <- function(requested.learning.rate = NULL, include.slow = TRUE) {
  lr.candidates <- unique(c(
    requested.learning.rate,
    0.01, 0.005, 0.001, 0.0005, 0.0002,
    if (isTRUE(include.slow)) 0.0001
  ))
  lr.candidates <- lr.candidates[is.finite(lr.candidates)]
  lr.candidates[lr.candidates > 0]
}

make_fit_params <- function(fit, strategy, learning.rate, bag.fraction, n.folds,
                            tree.complexity, target.min, max.trees, step.size = NA_integer_,
                            tolerance = NA_real_, n.minobsinnode = NA_integer_) {
  if (is.null(fit)) {
    return(list(
      strategy = strategy,
      learning.rate = learning.rate,
      bag.fraction = bag.fraction,
      n.folds = n.folds,
      tree.complexity = tree.complexity,
      target.min = target.min,
      max.trees = max.trees,
      step.size = step.size,
      tolerance = tolerance,
      n.minobsinnode = n.minobsinnode,
      acceptable = FALSE
    ))
  }
  list(
    strategy = strategy,
    learning.rate = learning.rate,
    bag.fraction = bag.fraction,
    n.folds = n.folds,
    tree.complexity = tree.complexity,
    target.min = target.min,
    max.trees = max.trees,
    step.size = step.size,
    tolerance = tolerance,
    n.minobsinnode = n.minobsinnode,
    n.trees = fit$n.trees,
    acceptable = fit$n.trees >= target.min && fit$n.trees < floor(max.trees * 0.95)
  )
}

call_gbm_step_safe <- function(gbm.args) {
  last.error <- NULL
  fit.try <- tryCatch(
    suppressMessages(suppressWarnings(do.call(gbm.step.patched, gbm.args))),
    error = function(e) {
      last.error <<- conditionMessage(e)
      NULL
    }
  )
  list(fit = fit.try, error = last.error)
}

adaptive_gbm_step <- function(dat, gbm.x, gbm.y, weights = NULL, family = "gaussian",
                              max.trees = 100000, tree.complexity = 2, n.folds = 10,
                              target.min = 1000, silent = FALSE, plot.main = FALSE,
                              plot.folds = FALSE, bag.candidates = NULL,
                              strategy = "default", stop.on.acceptable = TRUE,
                              n.minobs.multiplier = 1,
                              jitter.threshold = 0.12,
                              max.attempts = 9L,
                              extra.args = list(),
                              lr.candidates = NULL) {
  n.folds <- choose_gbm_n_folds(dat, gbm.y, family = family, requested.n.folds = n.folds)
  n <- nrow(dat)
  n.train <- floor(n * (n.folds - 1L) / n.folds)
  if (is.null(bag.candidates)) {
    bag.candidates <- build_bag_candidates(n = n)
  }
  if (is.null(lr.candidates)) {
    lr.candidates <- build_lr_candidates()
  }
  max.target.trees <- floor(max.trees * 0.95)
  best.fit <- NULL
  best.params <- NULL
  last.error <- NULL
  attempt.count <- 0L
  
  for (bag.fraction in bag.candidates) {
    step.size <- max(5L, round(n * bag.fraction / 2))
    tolerance <- max(1e-06, 0.001 * sqrt(30 / max(2, n * bag.fraction)))
    base.n.minobs <- max(2L, floor((n.train * bag.fraction - 2) / 2) - 1L)
    n.minobs <- max(2L, round(base.n.minobs * n.minobs.multiplier))
    gbm.args.base <- modifyList(extra.args, list(
      data = dat,
      gbm.x = gbm.x,
      gbm.y = gbm.y,
      family = family,
      max.trees = max.trees,
      tolerance = tolerance,
      bag.fraction = bag.fraction,
      tree.complexity = tree.complexity,
      step.size = step.size,
      n.folds = n.folds,
      n.minobsinnode = n.minobs,
      silent = silent,
      plot.main = plot.main,
      plot.folds = plot.folds
    ))
    if (!is.null(weights)) {
      gbm.args.base$weights <- weights
    }
    
    for (learning.rate in lr.candidates) {
      if (attempt.count >= max.attempts &&
          !is.null(best.fit) &&
          !is.null(best.params) &&
          isTRUE(best.params$acceptable)) {
        return(list(fit = best.fit, params = best.params, error = last.error,
                    attempts = attempt.count))
      }
      attempt.count <- attempt.count + 1L
      gbm.args.try <- gbm.args.base
      gbm.args.try$learning.rate <- learning.rate
      fit.try <- call_gbm_step_safe(gbm.args.try)
      if (!is.null(fit.try$fit)) {
        params.try <- make_fit_params(
          fit = fit.try$fit,
          strategy = strategy,
          learning.rate = learning.rate,
          bag.fraction = bag.fraction,
          n.folds = n.folds,
          tree.complexity = tree.complexity,
          target.min = target.min,
          max.trees = max.trees,
          step.size = step.size,
          tolerance = tolerance,
          n.minobsinnode = n.minobs
        )
        params.try$jitter <- fit_cv_jitter(fit.try$fit)
        params.try$jittery <- isTRUE(params.try$jitter > jitter.threshold)
        params.try$acceptable <- fit.try$fit$n.trees >= target.min &&
          fit.try$fit$n.trees < max.target.trees &&
          !params.try$jittery
        if (fit_is_better(fit.try$fit, params.try, best.fit, best.params)) {
          best.fit <- fit.try$fit
          best.params <- params.try
        }
        if (isTRUE(stop.on.acceptable) && isTRUE(params.try$acceptable)) {
          return(list(fit = best.fit, params = best.params, error = NULL,
                      attempts = attempt.count))
        }
      } else {
        last.error <- fit.try$error
      }
    }
  }
  
  list(fit = best.fit, params = best.params, error = last.error,
       attempts = attempt.count)
}

refine_weak_gbm_step <- function(dat, gbm.x, gbm.y, weights = NULL, family = "gaussian",
                                 max.trees = 100000, tree.complexity = 2, n.folds = 10,
                                 target.min = 1000, silent = FALSE, plot.main = FALSE,
                                 plot.folds = FALSE, weak.cor.threshold = 0.15,
                                 n.minobs.multiplier = 1, extra.args = list(),
                                 jitter.threshold = 0.12, bag.candidates = NULL,
                                 lr.candidates = NULL) {
  if (is.null(bag.candidates)) {
    bag.candidates <- build_bag_candidates(n = nrow(dat))
  }
  if (is.null(lr.candidates)) {
    lr.candidates <- build_lr_candidates()
  }
  primary.max.attempts <- max(9L, length(bag.candidates) * length(lr.candidates))
  primary.fit <- adaptive_gbm_step(
    dat,
    gbm.x = gbm.x,
    gbm.y = gbm.y,
    weights = weights,
    family = family,
    max.trees = max.trees,
    tree.complexity = tree.complexity,
    n.folds = n.folds,
    target.min = target.min,
    silent = silent,
    plot.main = plot.main,
    plot.folds = plot.folds,
    bag.candidates = bag.candidates,
    strategy = "default",
    stop.on.acceptable = FALSE,
    n.minobs.multiplier = n.minobs.multiplier,
    jitter.threshold = jitter.threshold,
    max.attempts = primary.max.attempts,
    extra.args = extra.args,
    lr.candidates = lr.candidates
  )
  
  best.fit <- primary.fit$fit
  best.params <- primary.fit$params
  best.error <- primary.fit$error
  total.attempts <- primary.fit$attempts %||% 0L
  weak.initial <- fit_is_weak(primary.fit$fit, weak.cor.threshold = weak.cor.threshold,
                              max.trees = max.trees, jitter.threshold = jitter.threshold)
  weak.refined <- FALSE
  
  if (weak.initial) {
    reduced.folds <- max(2L, min(5L, n.folds - 2L))
    alt.specs <- list(
      list(strategy = "simpler_tree",
           tree.complexity = max(1L, tree.complexity - 1L),
           n.folds = n.folds,
           target.min = max(300L, round(target.min * 0.75)),
           bag.candidates = unique(c(0.85, bag.candidates)),
           lr.candidates = c(0.005, 0.001, 0.0005),
           max.attempts = 4L),
      list(strategy = "smoother_curve",
           tree.complexity = tree.complexity,
           n.folds = n.folds,
           target.min = max(400L, round(target.min * 0.7)),
           bag.candidates = unique(c(0.9, 0.8)),
           lr.candidates = c(0.001, 0.0005),
           max.attempts = 4L),
      list(strategy = "slow_learning",
           tree.complexity = tree.complexity,
           n.folds = n.folds,
           target.min = max(400L, round(target.min * 0.7)),
           bag.candidates = unique(c(0.75, 0.7, 0.5, bag.candidates)),
           lr.candidates = build_lr_candidates(include.slow = TRUE)[
             build_lr_candidates(include.slow = TRUE) <= 0.0005
           ],
           stop.on.acceptable = FALSE,
           max.attempts = 6L),
      list(strategy = "reduced_folds",
           tree.complexity = 1L,
           n.folds = reduced.folds,
           target.min = max(200L, round(target.min * 0.5)),
           bag.candidates = unique(c(0.9, 0.8)),
           lr.candidates = c(0.01, 0.005, 0.001),
           max.attempts = 4L)
    )
    
    for (spec in alt.specs) {
      fit.try <- adaptive_gbm_step(
        dat,
        gbm.x = gbm.x,
        gbm.y = gbm.y,
        weights = weights,
        family = family,
        max.trees = max.trees,
        tree.complexity = spec$tree.complexity,
        n.folds = spec$n.folds,
        target.min = spec$target.min,
        silent = silent,
        plot.main = plot.main,
        plot.folds = plot.folds,
        bag.candidates = spec$bag.candidates,
        strategy = spec$strategy,
        stop.on.acceptable = spec$stop.on.acceptable %||% TRUE,
        n.minobs.multiplier = n.minobs.multiplier,
        jitter.threshold = jitter.threshold,
        max.attempts = spec$max.attempts,
        extra.args = extra.args,
        lr.candidates = spec$lr.candidates
      )
      total.attempts <- total.attempts + (fit.try$attempts %||% 0L)
      if (fit_is_better(fit.try$fit, fit.try$params, best.fit, best.params)) {
        best.fit <- fit.try$fit
        best.params <- fit.try$params
        best.error <- fit.try$error
        weak.refined <- TRUE
      }
      if (!fit_is_weak(best.fit, weak.cor.threshold = weak.cor.threshold,
                       max.trees = max.trees, jitter.threshold = jitter.threshold)) {
        break
      }
    }
  }
  
  list(
    fit = best.fit,
    params = best.params,
    error = best.error,
    attempts = total.attempts,
    weak.initial = weak.initial,
    weak.refined = weak.refined,
    initial.cv.cor = fit_cv_correlation(primary.fit$fit),
    final.cv.cor = fit_cv_correlation(best.fit)
  )
}

robust_gbm_step <- function(data, ...) {
  gbm.args <- list(...)
  gbm.args$data <- data
  if (is.null(gbm.args$gbm.x) || is.null(gbm.args$gbm.y)) {
    return(do.call(gbm.step.patched, gbm.args))
  }
  validate_gbm_columns(data, gbm.args$gbm.x, "gbm.x")
  validate_gbm_columns(data, gbm.args$gbm.y, "gbm.y")
  
  family <- gbm.args$family %||% "gaussian"
  max.trees <- gbm.args$max.trees %||% 100000
  tree.complexity <- gbm.args$tree.complexity %||% 1L
  n.folds <- choose_gbm_n_folds(
    data,
    gbm.args$gbm.y,
    family = family,
    requested.n.folds = gbm.args$n.folds %||% 10L
  )
  silent <- isTRUE(gbm.args$silent)
  plot.main <- if (is.null(gbm.args$plot.main)) TRUE else isTRUE(gbm.args$plot.main)
  plot.folds <- if (is.null(gbm.args$plot.folds)) FALSE else isTRUE(gbm.args$plot.folds)
  target.min <- if (nrow(data) < 120) 500L else 1000L
  model.label <- get_gbm_label(data, gbm.args$gbm.y)
  bag.candidates <- build_bag_candidates(gbm.args$bag.fraction %||% NULL, n = nrow(data))
  lr.candidates <- build_lr_candidates(gbm.args$learning.rate %||% NULL)
  
  extra.args <- gbm.args[setdiff(names(gbm.args), c(
    "data", "gbm.x", "gbm.y", "weights", "family", "max.trees", "tolerance",
    "learning.rate", "bag.fraction", "tree.complexity", "step.size", "n.folds",
    "silent", "plot.main", "plot.folds", "n.minobsinnode"
  ))]
  
  adaptive.fit <- refine_weak_gbm_step(
    data,
    gbm.x = gbm.args$gbm.x,
    gbm.y = gbm.args$gbm.y,
    weights = gbm.args$weights %||% NULL,
    family = family,
    max.trees = max.trees,
    tree.complexity = tree.complexity,
    n.folds = n.folds,
    target.min = target.min,
    silent = silent,
    plot.main = plot.main,
    plot.folds = plot.folds,
    bag.candidates = bag.candidates,
    lr.candidates = lr.candidates,
    extra.args = extra.args
  )
  
  final.fit <- adaptive.fit$fit
  final.params <- adaptive.fit$params
  if (is.null(final.fit)) {
    stop(
      "adaptive gbm.step failed for ", model.label,
      if (!is.null(adaptive.fit$error)) paste0(": ", adaptive.fit$error) else "",
      call. = FALSE
    )
  }
  
  final.fit$adaptive.meta <- list(
    model = model.label,
    strategy = final.params$strategy %||% "default",
    attempts = adaptive.fit$attempts %||% NA_integer_,
    jitter = fit_cv_jitter(final.fit),
    weak.initial = isTRUE(adaptive.fit$weak.initial),
    weak.refined = isTRUE(adaptive.fit$weak.refined),
    initial.cv.cor = if (is.finite(adaptive.fit$initial.cv.cor)) adaptive.fit$initial.cv.cor else NA_real_,
    final.cv.cor = if (is.finite(fit_cv_correlation(final.fit))) fit_cv_correlation(final.fit) else NA_real_,
    params = final.params,
    error = adaptive.fit$error
  )
  
  if (!silent) {
    message(
      "adaptive gbm.step [", model.label, "] strategy=",
      final.fit$adaptive.meta$strategy,
      ", attempts=", final.fit$adaptive.meta$attempts,
      ", jitter=", round(final.fit$adaptive.meta$jitter, 4),
      ", weak.initial=", final.fit$adaptive.meta$weak.initial,
      ", weak.refined=", final.fit$adaptive.meta$weak.refined,
      ", initial.cv.cor=", round(100 * final.fit$adaptive.meta$initial.cv.cor, 2),
      ", final.cv.cor=", round(100 * final.fit$adaptive.meta$final.cv.cor, 2)
    )
  }
  
  final.fit
}

gbm.step.adaptive <- robust_gbm_step
gbm.step <- robust_gbm_step


mean_fun <- function(data, indices) {
  mean(data[indices])
}
  
  AICc <- function(...) {
    models <- list(...)
    num.mod <- length(models)
    AICcs <- numeric(num.mod)
    ns <- numeric(num.mod)
    ks <- numeric(num.mod)
    AICc.vec <- rep(0,num.mod)
    for (i in 1:num.mod) {
      if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
      if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
      AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
      ns[i] <- n
      ks[i] <- k
      AICc.vec[i] <- AICcs[i]
    }
    return(AICc.vec)
  }
  
  delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
  weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
  ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object
  
  linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
    fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
    AIC.vec <- c(AICc(fit.full),AICc(fit.null))
    dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
    ER <- wAIC.vec[1]/wAIC.vec[2]
    r.sq.adj <- as.numeric(summary(fit.full)[9])
    return(c(ER,r.sq.adj))
  }
  


## upload InvaPact database
setwd("~/Documents/Papers/Invasive species/InvaPact/invapact_v2.0")
ip <- get_invapact('invapact_v2.0.rds')

## filter options
ip_filt <- filter_invapact(ip,
                           duplicates = "strict",
                           geography = "remove_uncertain",
                           evidence = "demonstrated",
                           spatial_aggregate = 'aggregate_only')

## expand by site
ip_by_site   <- expand_by_site(ip_filt, "invapact_v2.0_sites.rds")

## colnames
colnames(ip_by_site)

## record-level country ISO3 assignment
cntry_vec <- names(table(ip_by_site$countryCode3))
lcntry <- length(table(ip_by_site$countryCode3))

boot_med <- boot_lo <- boot_up <- rep(NA, lcntry)
for (c in 1:lcntry) {
  cntry_dat <- subset(ip_by_site, countryCode3 == cntry_vec[c])
  if (dim(cntry_dat)[1] >= 10) {
    IPscores <- na.omit(as.vector(cntry_dat$invapactScore))
    bootIP <- boot(IPscores, statistic = mean_fun, R = 1000)
    boot_med[c] <- bootIP$t0
    boot_ci <- boot.ci(bootIP, type = "perc")
    boot_lo[c] <- boot_ci$percent[4]
    boot_up[c] <- boot_ci$percent[5]
  }
  if (dim(cntry_dat)[1] < 10) {
    boot_med[c] <- boot_lo[c] <- boot_up[c] <- NA
  }
}
cntry.dat <- data.frame('ISO3'=cntry_vec, 'IPscore_med'=boot_med, 'IPscore_up'=boot_up, 'IPscore_lo'=boot_lo)
head(cntry.dat)

##
## import socio-economic variables

## population
# population data
setwd("~/Documents/Papers/Health/pop trend & wealth/data/pop")
popdat <- read.csv("popXage.csv", header = TRUE, stringsAsFactors = FALSE)
head(popdat)

## pop total
popdat$Ntot <- round(1000 * apply(popdat[,4:dim(popdat)[2]], MARGIN=1, sum, na.rm=T), 0)
head(popdat)
dim(popdat)

## remove duplicates
popdat <- popdat[!duplicated(popdat[,c("cntry.code", "year")]),]
dim(popdat)

## Ntot
Ntot2020 <- popdat$Ntot[popdat$year == 2020]
pop_dat20 <- data.frame('ISO3'=popdat$cntry.code[popdat$year == 2020], 'Ntot'=Ntot2020)
head(pop_dat20)

# wealth data (domestic comprehensive wealth index)
# https://data360.worldbank.org/en/indicator/WB_CWON_NW_DOW?compBreak1=WB_CWON_PC
setwd("~/Documents/Papers/Health/pop trend & wealth/data/wealth")
wealthdat <- read.csv("DCWI.csv", header = TRUE, stringsAsFactors = FALSE)
head(wealthdat)

## compare to PPP per-capita GDP
setwd("~/Documents/Papers/Health/pop trend & wealth/data/wealth")
gdppcPPP <- read.csv("gdppcPPP.csv", header = TRUE, stringsAsFactors = FALSE)
head(gdppcPPP)
gdppcPPP2020 <- data.frame(cntry.code=gdppcPPP$cntry.cod, gdppcPPP2020=gdppcPPP$a2020)
head(gdppcPPP2020)

## take latest year only
cntry.vec2 <- attr(table(wealthdat$cntry.code), "names")
max.year <- rep(NA, length(cntry.vec2))
for (c in 1:length(cntry.vec2)) {
  cntry <- cntry.vec2[c]
  max.year[c] <- max(wealthdat$year[wealthdat$cntry.code == cntry])
} # end c loop
max.year

## latest year only (2020)
wealthdat2020 <- wealthdat[wealthdat$year == 2020,]
head(wealthdat2020)
dim(wealthdat2020)

## chained only
wealthdat2020chained <- wealthdat2020[wealthdat2020$UNIT_MEASURE == "USD_REAL_CHAINED_2019",]

## per capita only
wealthdat2020chainedpc <- wealthdat2020chained[wealthdat2020chained$COMP_BREAKDOWN_1_LABEL == "Aggregation: per capita",]
dim(wealthdat2020chainedpc)
head(wealthdat2020chainedpc)
wealthDCWI <- wealthdat2020chainedpc[,c("cntry.code", "DCWI")]
head(wealthDCWI)
dim(wealthDCWI)

## remove duplicates
wealthDCWI <- wealthDCWI[!duplicated(wealthDCWI$cntry.code),]
dim(wealthDCWI)

## compare to GDP pc PPP
head(wealthDCWI)
DCWIgdp <- merge(wealthDCWI, gdppcPPP2020, by = "cntry.code", all.x = TRUE)
head(DCWIgdp)

ggplot(DCWIgdp, aes(x = gdppcPPP2020, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "per-capita GDP PPP (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()
fitDCWIgdp <- lm(log10(DCWI) ~ log10(gdppcPPP2020), data = DCWIgdp)
summary(fitDCWIgdp)

## merge DCWIgdp with population data
DCWIgdppop <- merge(DCWIgdp, pop_dat20, by.x = "cntry.code", by.y = "ISO3", all.x = TRUE)
head(DCWIgdppop)

## merge with InvaPact country-level data
IPcntry_wlth_gdp_pop <- merge(cntry.dat, DCWIgdppop, by.x = "ISO3", by.y = "cntry.code", all.x = TRUE)
head(IPcntry_wlth_gdp_pop)

IPcntry_wlth_gdp_pop_noZero <- subset(IPcntry_wlth_gdp_pop, !is.na(IPscore_med) & IPscore_med > 0)


##  bubble ggplot of InvaPact score vs wealth index, with bubble size proportional to populationp
ggplot(IPcntry_wlth_gdp_pop_noZero, aes(x = DCWI, y = (IPscore_med), size=Ntot)) +
  geom_point(alpha = 0.5) +
  #geom_errorbar(aes(ymin = IPscore_lo, ymax = IPscore_up), width = 0.01, alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  scale_size_continuous(range = c(1, 10), name = "population (2020)") +
  scale_x_log10() +
  labs(x = "per-capita domestic comprehensive wealth index (2020)", y = "InvaPact score (median)") +
  theme_minimal()

ggplot(IPcntry_wlth_gdp_pop_noZero, aes(x = DCWI, y = (IPscore_med))) +
  geom_point(alpha = 0.5) +
  geom_errorbar(aes(ymin = IPscore_lo, ymax = IPscore_up), width = 0.01, alpha = 0.5, show.legend = FALSE) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  #scale_size_continuous(range = c(1, 10), name = "population (2020)") +
  scale_x_log10() +
  labs(x = "per-capita domestic comprehensive wealth index (2020)", y = "InvaPact score (median)") +
  theme_minimal()

linreg.ER(log10(IPcntry_wlth_gdp_pop_noZero$DCWI), (IPcntry_wlth_gdp_pop_noZero$IPscore_med))

## Gini index
## https://databank.worldbank.org/source/jobs/Series/SI.POV.GINI
setwd("~/Documents/Papers/Health/pop trend & wealth/data/Gini/")
gini.dat <- read.csv("WB_WDI_SI_POV_GINI.csv", header=T, stringsAsFactors = FALSE)
head(gini.dat)
gini.names <- colnames(gini.dat)
table(gini.dat[,gini.names[6]]) # cntry.code
table(gini.dat[,gini.names[24]])

gini.clean <- gini.dat[, c(gini.names[6], gini.names[24], gini.names[25])]
colnames(gini.clean) <- c("cntry.code","year","gini")
head(gini.clean)

## take average over last y years to maximise sample size
gini.cntry.vec <- attr(table(gini.dat[,gini.names[6]]), "names")

## earliest year back in time to consider for means
earliest.year <- 2010
gini.dat.mn <- data.frame(cntry.code=NA, giniMn=NA, yrMax=NA, yrMin=NA) # empty data frame
for (c in 1:length(gini.cntry.vec)) {
  gini.cntry <- subset(gini.clean, cntry.code == gini.cntry.vec[c])
  gini.cntry.update <- gini.cntry[which(gini.cntry$year >= earliest.year),]
  gini.mn <- mean(gini.cntry.update$gini, na.rm = TRUE)
  gini.yr.min <- min(gini.cntry.update$year, na.rm=T)
  gini.yr.max <- max(gini.cntry.update$year, na.rm=T)
  gini.cntry.mn <- data.frame(cntry.code=gini.cntry.vec[c],
                              giniMn=gini.mn, yrMax=gini.yr.max, yrMin=gini.yr.min)
  gini.dat.mn <- rbind(gini.dat.mn, gini.cntry.mn)
}
gini.dat.mn <- gini.dat.mn[-1,] # remove first row
gini.dat.mn <- na.omit(gini.dat.mn)
head(gini.dat.mn)
dim(gini.dat.mn)

## remove duplicates
gini.dat.mn <- gini.dat.mn[!duplicated(gini.dat.mn$cntry.code),]
head(gini.dat.mn)
dim(gini.dat.mn)
gini.use <- data.frame(cntry.code=gini.dat.mn$cntry.code, gini=gini.dat.mn$giniMn)
head(gini.use)


######################################
## Human Development Index          ##
## adjusted for planetary pressures ##
######################################
## https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
setwd("~/Documents/Papers/Health/pop trend & wealth/data/hdi")
hdipp <- read.csv("HDIPP.csv", header=T, stringsAsFactors=F)
head(hdipp)
hist(logit(hdipp$HDIPP23), main="")
hdipp.use <- data.frame(cntry.code=hdipp$cntry.code, HDIPP=hdipp$HDIPP23)
head(hdipp.use)

## freedom (Freedom House total [aggregate score for all categories], 2013-2025)									)
## https://freedomhouse.org/report/freedom-world#Data
## For each country and territory, Freedom in the World analyses the electoral process, political pluralism 
## and participation, the functioning of the government, freedom of expression and of belief, associational 
## and organizational rights, the rule of law, and personal autonomy and individual rights.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/Freedom House/")
freed2025 <- read.csv("freedom2025.csv", header=T, stringsAsFactors = F)
head(freed2025)
freed <- read.csv("freedom.csv", header=T, stringsAsFactors = F)
head(freed)

## fill in cntry.code column from associated cntry in freed2025 to freed
freedmrg <- merge(freed, freed2025[,c("cntry","cntry.code")], by="cntry", all.x=T, all.y=F)
head(freedmrg)

# remove cntry.code.x
freedfix <- freedmrg[,c("cntry", "cntry.code.y", "year", "freedom")]
head(freedfix)
colnames(freedfix)[2] <- "cntry.code"
head(freedfix)
freedclean <- freedfix[freedfix$cntry.code != "", ]
head(freedclean)

# sort by cntry.code & year
freedsort <- freedclean[order(freedclean$cntry.code, freedclean$year), ]
head(freedsort)

## return most recent year's values for each cntry.code
table(freedsort$year)
freed2025sort <- freedsort[freedsort$year == 2025, c("cntry.code", "freedom")]
head(freed2025sort)
hist(freed2025sort$freedom)
hist(logit(freed2025sort$freedom/100))
freed.use <- data.frame(cntry.code=freed2025sort$cntry.code, freedom=freed2025sort$freedom)
head(freed.use)

## corruption
## set working directory
setwd("~/Documents/Papers/Climate change/corruption")

## import corruption perception index 2024
cpi <- read.csv("cpi24.csv", header = TRUE)
head(cpi)

## corruption perceptions index
## https://www.transparency.org/
setwd("/Users/brad0317/Documents/GitHub/wealthwellbeingageingpop/data/")
cpi <- read.csv("cpi24.csv", header = TRUE)
head(cpi)
hist(cpi$cpi24)
hist(logit(cpi$cpi24/100))
range(cpi$cpi24)

## imports of goods and services (IGS)
## https://data.worldbank.org/indicator/NE.IMP.GNFS.ZS
setwd("~/Documents/GitHub/InvaPactCountry/data/IGS")
igs_dat <- read.csv("IGS.csv", header = TRUE, stringsAsFactors = FALSE)
head(igs_dat)

# average of latest n years to maximise number of countries included
earliest.year <- "y2016"
igs_cntry.vec <- attr(table(igs_dat$cntry.code), "names")
igs.mn <- igs.sd <- rep(NA, length(igs_cntry.vec))
for (c in 1:length(igs_cntry.vec)) {
  igs.cntry <- subset(igs_dat, cntry.code == igs_cntry.vec[c])
  igs.cntry.igs <- as.numeric(igs.cntry[,which(colnames(igs.cntry) >= earliest.year)])
  igs.mn[c] <- mean(igs.cntry.igs, na.rm = TRUE)
  igs.sd[c] <- sd(igs.cntry.igs, na.rm = TRUE)
}
igs.dat.mn <- data.frame(cntry.code=igs_cntry.vec, igsMn=igs.mn, igsSD=igs.sd)
head(igs.dat.mn)


## ecological footprint
## footprintnetwork.org
## https://footprint.info.yorku.ca/data/
setwd("~/Documents/GitHub/InvaPactCountry/data/footprint")
ef_cntry <- read_excel(
  "NEFBA_Data_global_national_ecological_footprint_biocapacity.xlsx",
  sheet = "national_data"
)
head(ef_cntry)
colnames(ef_cntry)

ef_world <- read_excel(
  "NEFBA_Data_global_national_ecological_footprint_biocapacity.xlsx",
  sheet = "world_data"
)
head(ef_world)
colnames(ef_world)

bc_world <- ef_world %>%
  select(
    year,
    bc_pc_world = biocap_total_percapita
  )

ef_earths <- ef_cntry %>%
  left_join(bc_world, by = "year") %>%
  mutate(
    Earths = efc_total_percapita / bc_pc_world
  )
head(ef_earths)

# take average of last 10 years to maximise number of countries with data
ef_earths_recent <- ef_earths %>%
  filter(year >= 2016) %>%
  group_by(country_ISO3) %>%
  summarise(EarthsMn = mean(Earths, na.rm = TRUE))
head(ef_earths_recent)


## GBIF records (diversity index)
setwd("~/Documents/GitHub/InvaPactCountry/data/GBIF")
div <- read.table("gbif_records_country_year.csv", header=T, sep = ",", stringsAsFactors = FALSE)
head(div)

## most recent year and iso3c only
div_recent <- div %>%
  filter(year == max(year)) %>%
  select(cntry.code=iso3c, gbif=gbifRecords)
head(div_recent)

##########
## merge
gini_hdipp <- merge(gini.use, hdipp.use, by = "cntry.code", all.x = TRUE)
gini_hdipp_freed <- merge(gini_hdipp, freed.use, by = "cntry.code", all.x = TRUE)
gini_hdipp_freed_cpi <- merge(gini_hdipp_freed, cpi, by = "cntry.code", all.x = TRUE)
gini_hdpipp_freed_cpi_igs <- merge(gini_hdipp_freed_cpi, igs.dat.mn, by = "cntry.code", all.x = TRUE)
gini_hdipp_freed_cpi_igs_ef <- merge(gini_hdpipp_freed_cpi_igs, ef_earths_recent, by.x = "cntry.code", by.y = "country_ISO3", all.x = TRUE)
gini_hdipp_freed_cpi_igs_ef_gbif <- merge(gini_hdipp_freed_cpi_igs_ef, div_recent, by = "cntry.code", all.x = TRUE)
head(gini_hdipp_freed_cpi_igs_ef_gbif)

## merge with InvaPact country-level data
IPdat_merge <- merge(IPcntry_wlth_gdp_pop_noZero, gini_hdipp_freed_cpi_igs_ef_gbif, by.x = "ISO3", by.y = "cntry.code", all.x = TRUE) 
IPdat_merge$lNtot <- log10(IPdat_merge$Ntot)
IPdat_merge$lDCWI <- log10(IPdat_merge$DCWI)
IPdat_merge$lIGS <- log10(IPdat_merge$igsMn)
IPdat_merge$lEarthsMn <- log10(IPdat_merge$EarthsMn)
IPdat_merge$lGBIF <- log10(IPdat_merge$gbif)
head(IPdat_merge)

## missingness
colnames(IPdat_merge)
length(which(is.na(IPdat_merge$lDCWI) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$lNtot) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$gini) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$lIGS) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$HDIPP) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$freedom) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$cpi24) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$igsMn) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$EarthsMn) == T))/dim(IPdat_merge)[1]
length(which(is.na(IPdat_merge$lGBIF) == T))/dim(IPdat_merge)[1]

## multiple imputation of socio-economic predictors
if (!requireNamespace("mice", quietly = TRUE)) {
  stop(
    "mice package is required for multiple imputation; install with install.packages('mice').",
    call. = FALSE
  )
}

IPdat_merge_mice_predictors <- c(
  "DCWI", "Ntot", "gini", "HDIPP", "freedom", "cpi24",
  "igsMn", "EarthsMn", "gbif"
)
IPdat_merge_mice_derived_predictors <- c("lDCWI", "lNtot", "lIGS", "lEarthsMn", "lGBIF")
IPdat_merge_mice_variables <- c(
  "IPscore_med",
  IPdat_merge_mice_predictors,
  IPdat_merge_mice_derived_predictors
)
missing_mice_variables <- setdiff(IPdat_merge_mice_variables, names(IPdat_merge))
if (length(missing_mice_variables)) {
  stop(
    "variables required for MICE absent: ",
    paste(missing_mice_variables, collapse = ", "),
    call. = FALSE
  )
}
non_numeric_mice_variables <- IPdat_merge_mice_variables[
  !vapply(IPdat_merge[IPdat_merge_mice_variables], is.numeric, logical(1))
]
if (length(non_numeric_mice_variables)) {
  stop(
    "MICE variables must be numeric: ",
    paste(non_numeric_mice_variables, collapse = ", "),
    call. = FALSE
  )
}
positive_mice_variables <- c("DCWI", "Ntot", "igsMn", "EarthsMn", "gbif")
nonpositive_mice_values <- positive_mice_variables[
  vapply(
    IPdat_merge[positive_mice_variables],
    function(x) any(!is.na(x) & x <= 0),
    logical(1)
  )
]
if (length(nonpositive_mice_values)) {
  stop(
    "variables must be positive to derive log10 transforms: ",
    paste(nonpositive_mice_values, collapse = ", "),
    call. = FALSE
  )
}

# log-transformed predictors  recalculated after each imputation rather
# than independently imputed, preserving their deterministic relationship to
# DCWI, Ntot, and igsMn.
IPdat_merge_mice_data <- IPdat_merge[, IPdat_merge_mice_variables, drop = FALSE]
non_finite_mice_values <- vapply(
  IPdat_merge_mice_data[IPdat_merge_mice_predictors],
  function(x) sum(!is.finite(x) & !is.na(x)),
  integer(1)
)
if (any(non_finite_mice_values > 0L)) {
  message(
    "treating non-finite predictor values as missing for MICE: ",
    paste(
      names(non_finite_mice_values)[non_finite_mice_values > 0L],
      collapse = ", "
    )
  )
  for (variable in names(non_finite_mice_values)[non_finite_mice_values > 0L]) {
    IPdat_merge_mice_data[[variable]][
      !is.finite(IPdat_merge_mice_data[[variable]])
    ] <- NA_real_
  }
}

imputed_mice_variables <- setdiff(IPdat_merge_mice_predictors, "lIGS")
insufficient_observed_values <- imputed_mice_variables[
  vapply(
    IPdat_merge_mice_data[imputed_mice_variables],
    function(x) sum(is.finite(x)) < 2L,
    logical(1)
  )
]
if (length(insufficient_observed_values)) {
  stop(
    "MICE requires at least two observed values for: ",
    paste(insufficient_observed_values, collapse = ", "),
    call. = FALSE
  )
}

IPdat_merge_mice_method <- rep("", length(IPdat_merge_mice_variables))
names(IPdat_merge_mice_method) <- IPdat_merge_mice_variables
# CART conditional models avoid singular donor regressions in this small,
# highly collinear country-level dataset while retaining chained imputation.
IPdat_merge_mice_method[imputed_mice_variables] <- "cart"

IPdat_merge_mice_predictor_matrix <- matrix(
  0L,
  nrow = length(IPdat_merge_mice_variables),
  ncol = length(IPdat_merge_mice_variables),
  dimnames = list(IPdat_merge_mice_variables, IPdat_merge_mice_variables)
)
IPdat_merge_mice_predictor_matrix[imputed_mice_variables, "IPscore_med"] <- 1L
IPdat_merge_mice_predictor_matrix[
  imputed_mice_variables,
  imputed_mice_variables
] <- 1L
diag(IPdat_merge_mice_predictor_matrix) <- 0L

set.seed(20261015)
IPdat_merge_mice <- mice::mice(
  IPdat_merge_mice_data,
  m = 20L,
  maxit = 20L,
  method = IPdat_merge_mice_method,
  predictorMatrix = IPdat_merge_mice_predictor_matrix,
  printFlag = FALSE,
  seed = 20261015
)

# keep all 20 completed datasets so downstream models can pool uncertainty
# rather than treating a single completed dataset as observed data.
IPdat_merge_imputed <- lapply(seq_len(IPdat_merge_mice$m), function(imputation) {
  completed_predictors <- mice::complete(IPdat_merge_mice, action = imputation)
  completed_data <- IPdat_merge
  completed_data[IPdat_merge_mice_variables] <- completed_predictors
  completed_data$lDCWI <- log10(completed_data$DCWI)
  completed_data$lNtot <- log10(completed_data$Ntot)
  completed_data$lIGS <- log10(completed_data$igsMn)
  completed_data$lEarthsMn <- log10(completed_data$EarthsMn)
  completed_data$lGBIF <- log10(completed_data$gbif)
  completed_data
})
names(IPdat_merge_imputed) <- paste0("imputation_", seq_along(IPdat_merge_imputed))

# Mean-completed values provide one deterministic dataset for downstream
# analyses, while the individual completed datasets retain imputation variance.
IPdat_merge_imputed_mean <- IPdat_merge
for (variable in imputed_mice_variables) {
  IPdat_merge_imputed_mean[[variable]] <- rowMeans(vapply(
    IPdat_merge_imputed,
    `[[`,
    numeric(nrow(IPdat_merge)),
    variable
  ))
}
IPdat_merge_imputed_mean$lDCWI <- log10(IPdat_merge_imputed_mean$DCWI)
IPdat_merge_imputed_mean$lNtot <- log10(IPdat_merge_imputed_mean$Ntot)
IPdat_merge_imputed_mean$lIGS <- log10(IPdat_merge_imputed_mean$igsMn)
IPdat_merge_imputed_mean$lEarthsMn <- log10(IPdat_merge_imputed_mean$EarthsMn)
IPdat_merge_imputed_mean$lGBIF <- log10(IPdat_merge_imputed_mean$gbif)
head(IPdat_merge_imputed_mean)

# this long table records value used per country and whether it was
# originally missing, allowing downstream analyses to filter imputed entries.
IPdat_merge_imputation_means <- do.call(
  rbind,
  lapply(imputed_mice_variables, function(variable) {
    data.frame(
      ISO3 = IPdat_merge$ISO3,
      variable = variable,
      was_imputed = is.na(IPdat_merge[[variable]]),
      mean_imputed_value = IPdat_merge_imputed_mean[[variable]],
      stringsAsFactors = FALSE
    )
  })
)
rownames(IPdat_merge_imputation_means) <- NULL


## individual variable ggplots
ggplot(IPdat_merge_imputed_mean, aes(x = Ntot, y = (IPscore_med))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  scale_x_log10() +
  labs(x = "total population (2020)", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(log10(IPdat_merge_imputed_mean$Ntot),
          IPdat_merge_imputed_mean$IPscore_med)

ggplot(IPdat_merge_imputed_mean, aes(x = lDCWI, y = (IPscore_med))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  scale_x_log10() +
  labs(x = "per-capita domestic comprehensive wealth index (2020)", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(log10(IPdat_merge_imputed_mean$DCWI), (IPdat_merge_imputed_mean$IPscore_med))

ggplot(IPdat_merge_imputed_mean, aes(x = gini, y = IPscore_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  labs(x = "Gini index (2020)", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(IPdat_merge_imputed_mean$gini, (IPdat_merge_imputed_mean$IPscore_med))

ggplot(IPdat_merge_imputed_mean, aes(x = HDIPP, y = IPscore_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  labs(x = "Human Development Index adjusted for planetary pressures (2023)", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(IPdat_merge_imputed_mean$HDIPP, (IPdat_merge_imputed_mean$IPscore_med))

ggplot(IPdat_merge_imputed_mean, aes(x = freedom, y = IPscore_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  labs(x = "Freedom House score (2025)", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(IPdat_merge_imputed_mean$freedom, (IPdat_merge_imputed_mean$IPscore_med))

ggplot(IPdat_merge_imputed_mean, aes(x = cpi24, y = IPscore_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  labs(x = "Corruption Perceptions Index (2024)", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(IPdat_merge_imputed_mean$cpi24, (IPdat_merge_imputed_mean$IPscore_med))

ggplot(IPdat_merge_imputed_mean, aes(x = lIGS, y = IPscore_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  labs(x = "mean intl trade of goods & services", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(IPdat_merge_imputed_mean$lIGS, (IPdat_merge_imputed_mean$IPscore_med))

ggplot(IPdat_merge_imputed_mean, aes(x = lEarthsMn, y = IPscore_med)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", linetype="dashed", se = T, color = "blue", show.legend = FALSE) +
  labs(x = "biocapacity deficit", y = "InvaPact score (median)") +
  theme_minimal()
linreg.ER(IPdat_merge_imputed_mean$lEarthsMn, (IPdat_merge_imputed_mean$IPscore_med))

## create Spearman's correlation matrix among the socio-economic variables
cor_mat <- cor(IPdat_merge_imputed_mean[,c("lNtot", "lDCWI", "gini", "HDIPP", "freedom", "cpi24", "lIGS", "lEarthsMn", "lGBIF")], 
               use = "pairwise.complete.obs", method = "spearman")
## add NA to upper triangle and diagonal
diag(cor_mat) <- NA
cor_mat[upper.tri(cor_mat)] <- NA
cor_mat <- cor_mat[,-dim(cor_mat)[2]]
cor_mat <- cor_mat[-1,]
cor_mat

## ggplot histograms of each socio-economic variable
ggplot(IPdat_merge_imputed_mean, aes(x = log10(Ntot))) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(x = "population size (log10)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = log10(DCWI))) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(x = "per-capita domestic comprehensive wealth index (2020)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = gini)) +
  geom_histogram(bins = 30, fill = "lightgreen", color = "black", alpha = 0.7) +
  labs(x = "Gini index (2020)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = (HDIPP))) +
  geom_histogram(bins = 30, fill = "lightyellow", color = "black", alpha = 0.7) +
  labs(x = "Human Development Index adjusted for planetary pressures (2023)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = freedom)) +
  geom_histogram(bins = 30, fill = "lightcoral", color = "black", alpha = 0.7) +
  labs(x = "Freedom House score (2025)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = cpi24)) +
  geom_histogram(bins = 30, fill = "lightgray", color = "black", alpha = 0.7) +
  labs(x = "Corruption Perceptions Index (2024)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = lIGS)) +
  geom_histogram(bins = 30, fill = "lightgray", color = "black", alpha = 0.7) +
  labs(x = "mean intl trade of goods & services (log10)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = log10(EarthsMn))) +
  geom_histogram(bins = 30, fill = "lightgray", color = "black", alpha = 0.7) +
  labs(x = "biocapacity deficit (log10)", y = "Frequency") +
  theme_minimal()

ggplot(IPdat_merge_imputed_mean, aes(x = log10(gbif))) +
  geom_histogram(bins = 30, fill = "lightgray", color = "black", alpha = 0.7) +
  labs(x = "GBIF records (log10)", y = "Frequency") +
  theme_minimal()


### boosted regression tree (BRT) model of InvaPact score vs socio-economic variables
IPdat_merge.predictors <- c("lNtot","lDCWI", "gini", "HDIPP", "freedom", "cpi24", "lIGS", "lEarthsMn", "lGBIF")
IPdat_merge.response <- "IPscore_med"

# This is a country-level model. Restrict it to one fully observed row per
# country so missing socio-economic data cannot alter the cross-validation folds.
IPdat_merge_brt_data <- IPdat_merge_imputed_mean[, c(
  "ISO3", IPdat_merge.response, IPdat_merge.predictors
), drop = FALSE]

model_variables <- c(IPdat_merge.response, IPdat_merge.predictors)
non_numeric_variables <- model_variables[
  !vapply(IPdat_merge_brt_data[model_variables], is.numeric, logical(1))
]
if (length(non_numeric_variables)) {
  stop(
    "BRT variables must be numeric: ",
    paste(non_numeric_variables, collapse = ", "),
    call. = FALSE
  )
}

complete_brt_rows <- complete.cases(IPdat_merge_brt_data[model_variables]) &
  apply(
    IPdat_merge_brt_data[model_variables],
    1,
    function(x) all(is.finite(x))
  )
IPdat_merge_brt_data <- IPdat_merge_brt_data[complete_brt_rows, , drop = FALSE]
rownames(IPdat_merge_brt_data) <- NULL

if (anyDuplicated(IPdat_merge_brt_data$ISO3)) {
  stop("BRT data must contain one row per ISO3 country code.", call. = FALSE)
}
if (nrow(IPdat_merge_brt_data) < 10L * length(IPdat_merge.predictors)) {
  stop(
    "BRT requires at least 10 complete country records per predictor; found ",
    nrow(IPdat_merge_brt_data), ".",
    call. = FALSE
  )
}
constant_variables <- model_variables[
  vapply(
    IPdat_merge_brt_data[model_variables],
    function(x) diff(range(x)) == 0,
    logical(1)
  )
]
if (length(constant_variables)) {
  stop(
    "BRT variables have no variation: ",
    paste(constant_variables, collapse = ", "),
    call. = FALSE
  )
}

# 5 folds provide sufficiently large held-out country sets for this small
# country-level dataset. The seed makes fold allocation and bagging reproducible.
set.seed(20261013)
IPdat_merge_brt <- gbm.step.adaptive(
  IPdat_merge_brt_data,
  gbm.x = IPdat_merge.predictors,
  gbm.y = IPdat_merge.response,
  family = "gaussian",
  max.trees = 100000,
  tree.complexity = 2,
  n.folds = 5,
  tolerance.method = "auto"
)
IPdat_merge.summ <- summary(IPdat_merge_brt)
IPdat_merge.summ

IPdat_merge.plot <- ggplot(IPdat_merge.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
IPdat_merge.plot.flip <- IPdat_merge.plot + coord_flip()
IPdat_merge.plot.flip

IPdat_merge.CV.cor <- 100 * IPdat_merge_brt$cv.statistics$correlation.mean
IPdat_merge.CV.cor.se <- 100 * IPdat_merge_brt$cv.statistics$correlation.se
print(c(IPdat_merge.CV.cor, IPdat_merge.CV.cor.se))

predictor_labels <- c(
  lNtot = "log10(total population)",
  lDCWI = "log10(per-capita domestic comprehensive wealth index)",
  gini = "Gini coefficient",
  HDIPP = "Human Development Index adjusted for planetary pressures",
  freedom = "Freedom House score",
  lIGS = "log10(mean intl trade of goods & services)",
  cpi24 = "Corruption Perceptions Index",
  lEarthsMn = "log10(mean earth system metrics)",
  lGBIF = "log10(mean GBIF records)"
)
for (predictor in IPdat_merge.predictors) {
  gbm.plot(
    IPdat_merge_brt,
    variable.no = match(predictor, IPdat_merge.predictors),
    smooth = TRUE,
    rug = TRUE,
    common.scale = TRUE,
    write.title = FALSE,
    show.contrib = TRUE,
    y.label = "median InvaPact Score",
    x.label = predictor_labels[[predictor]],
    plot.layout = c(1, 1)
  )
}

## Resampled BRT: response uncertainty and continent-stratified country bootstrap
cont_cntry <- read.csv(
  "~/Documents/GitHub/InvaPactCountry/data/cont_cntry.csv",
  stringsAsFactors = FALSE
)
cont_lookup <- cont_cntry[, c("cntry.code", "cont")]
continents_by_code <- split(cont_lookup$cont, cont_lookup$cntry.code)
conflicting_continents <- names(continents_by_code)[
  vapply(
    continents_by_code,
    function(x) length(unique(stats::na.omit(x))) > 1L,
    logical(1)
  )
]
if (length(conflicting_continents)) {
  stop(
    "conflicting continent assignments for: ",
    paste(conflicting_continents, collapse = ", "),
    call. = FALSE
  )
}
cont_lookup <- cont_lookup[!duplicated(cont_lookup$cntry.code), ]

IPdat_merge_brt_bounds <- IPdat_merge_imputed_mean[, c("ISO3", "IPscore_lo", "IPscore_up")]
IPdat_merge_brt_bounds <- IPdat_merge_brt_bounds[
  !duplicated(IPdat_merge_brt_bounds$ISO3),
]
IPdat_merge_brt_resample_data <- merge(
  IPdat_merge_brt_data,
  IPdat_merge_brt_bounds,
  by = "ISO3",
  all.x = TRUE,
  sort = FALSE
)
IPdat_merge_brt_resample_data <- merge(
  IPdat_merge_brt_resample_data,
  cont_lookup,
  by.x = "ISO3",
  by.y = "cntry.code",
  all.x = TRUE,
  sort = FALSE
)

missing_continents <- IPdat_merge_brt_resample_data$ISO3[
  is.na(IPdat_merge_brt_resample_data$cont)
]
if (length(missing_continents)) {
  stop(
    "no continent assignment for: ",
    paste(missing_continents, collapse = ", "),
    call. = FALSE
  )
}

invalid_score_intervals <- with(
  IPdat_merge_brt_resample_data,
  !is.finite(IPscore_lo) | !is.finite(IPscore_up) |
    IPscore_lo > IPscore_med | IPscore_up < IPscore_med | IPscore_lo > IPscore_up
)
if (any(invalid_score_intervals)) {
  stop(
    "invalid InvaPact score intervals for: ",
    paste(IPdat_merge_brt_resample_data$ISO3[invalid_score_intervals],
          collapse = ", "),
    call. = FALSE
  )
}

# each iteration resamples countries with replacement within cont; retains
# every continent and its observed sample size, unlike equal-size regional draws
# that would discard most countries in sparsely represented continents.
IPdat_merge_brt_resamples_n <- 250L ## number of resampled BRT fits
IPdat_merge_brt_resample_seed <- 20261014L
continent_rows <- split(
  seq_len(nrow(IPdat_merge_brt_resample_data)),
  IPdat_merge_brt_resample_data$cont
)
if (length(continent_rows) < 2L) {
  stop("spatial resampling requires at least two continents.", call. = FALSE)
}

IPdat_merge_brt_resamples <- data.frame(
  iteration = seq_len(IPdat_merge_brt_resamples_n),
  cv_correlation = NA_real_,
  cv_correlation_se = NA_real_,
  model_error = NA_character_,
  stringsAsFactors = FALSE
)
IPdat_merge_brt_resample_influence <- matrix(
  NA_real_,
  nrow = IPdat_merge_brt_resamples_n,
  ncol = length(IPdat_merge.predictors),
  dimnames = list(NULL, IPdat_merge.predictors)
)
IPdat_merge_brt_resample_grid_points <- 100L
IPdat_merge_brt_resample_partial_x <- array(
  NA_real_,
  dim = c(
    IPdat_merge_brt_resample_grid_points,
    length(IPdat_merge.predictors),
    IPdat_merge_brt_resamples_n
  ),
  dimnames = list(
    NULL,
    IPdat_merge.predictors,
    NULL
  )
)
IPdat_merge_brt_resample_partial_prediction <- IPdat_merge_brt_resample_partial_x

for (iteration in seq_len(IPdat_merge_brt_resamples_n)) {
  set.seed(IPdat_merge_brt_resample_seed + iteration)
  resampled_rows <- unlist(
    lapply(
      continent_rows,
      function(rows) sample(rows, size = length(rows), replace = TRUE)
    ),
    use.names = FALSE
  )
  brt_resample <- IPdat_merge_brt_resample_data[resampled_rows, , drop = FALSE]
  brt_resample$IPscore_med <- brt_resample$IPscore_lo +
    runif(nrow(brt_resample)) *
    (brt_resample$IPscore_up - brt_resample$IPscore_lo)

  brt_fit <- tryCatch(
    gbm.step.adaptive(
      brt_resample,
      gbm.x = IPdat_merge.predictors,
      gbm.y = IPdat_merge.response,
      family = "gaussian",
      max.trees = 100000,
      tree.complexity = 2,
      n.folds = 5,
      tolerance.method = "auto",
      silent = TRUE,
      plot.main = FALSE,
      plot.folds = FALSE
    ),
    error = function(error) error
  )

  if (inherits(brt_fit, "error")) {
    IPdat_merge_brt_resamples$model_error[iteration] <- conditionMessage(brt_fit)
  } else {
    brt_summary <- summary(brt_fit, plotit = FALSE)
    IPdat_merge_brt_resamples$cv_correlation[iteration] <-
      brt_fit$cv.statistics$correlation.mean
    IPdat_merge_brt_resamples$cv_correlation_se[iteration] <-
      brt_fit$cv.statistics$correlation.se
    IPdat_merge_brt_resample_influence[iteration, brt_summary$var] <-
      brt_summary$rel.inf

    for (predictor_index in seq_along(IPdat_merge.predictors)) {
      partial_grid <- gbm::plot.gbm(
        brt_fit,
        i.var = predictor_index,
        continuous.resolution = IPdat_merge_brt_resample_grid_points,
        return.grid = TRUE
      )
      IPdat_merge_brt_resample_partial_x[, predictor_index, iteration] <-
        partial_grid[, 1]
      IPdat_merge_brt_resample_partial_prediction[, predictor_index, iteration] <-
        partial_grid[, 2]
    }
  }

  if (iteration %% 25L == 0L) {
    message("resampled BRT iteration ", iteration, " of ",
            IPdat_merge_brt_resamples_n)
  }
}

successful_brt_resamples <- is.finite(IPdat_merge_brt_resamples$cv_correlation)
if (!any(successful_brt_resamples)) {
  stop(
    "all resampled BRT fits failed: ",
    paste(unique(stats::na.omit(IPdat_merge_brt_resamples$model_error)),
          collapse = "; "),
    call. = FALSE
  )
}

resample_interval <- function(x) {
  x <- x[is.finite(x)]
  c(
    lower = unname(stats::quantile(x, 0.025)),
    median = stats::median(x),
    upper = unname(stats::quantile(x, 0.975))
  )
}

IPdat_merge_brt_resample_cv <- resample_interval(
  100 * IPdat_merge_brt_resamples$cv_correlation[successful_brt_resamples]
)
IPdat_merge_brt_resample_importance <- t(vapply(
  IPdat_merge.predictors,
  function(predictor) {
    resample_interval(IPdat_merge_brt_resample_influence[
      successful_brt_resamples, predictor
    ])
  },
  numeric(3)
))
IPdat_merge_brt_resample_importance <- data.frame(
  variable = rownames(IPdat_merge_brt_resample_importance),
  IPdat_merge_brt_resample_importance,
  row.names = NULL
)
IPdat_merge_brt_resample_importance <- IPdat_merge_brt_resample_importance[
  order(IPdat_merge_brt_resample_importance$median, decreasing = TRUE),
]

print(IPdat_merge_brt_resample_cv)
print(IPdat_merge_brt_resample_importance)

IPdat_merge_brt_resample_plot <- ggplot(
  IPdat_merge_brt_resample_importance,
  aes(x = reorder(variable, median), y = median)
) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  labs(x = "", y = "relative influence (median and 95% interval)") +
  theme_minimal()
IPdat_merge_brt_resample_plot

IPdat_merge_brt_resample_output_dir <- file.path(
  path.expand("~/Documents/GitHub/InvaPactCountry"),
  "output"
)
if (!dir.exists(IPdat_merge_brt_resample_output_dir) &&
    !dir.create(IPdat_merge_brt_resample_output_dir, recursive = TRUE)) {
  stop(
    "could not create resampled BRT output directory: ",
    IPdat_merge_brt_resample_output_dir,
    call. = FALSE
  )
}

IPdat_merge_brt_resample_partial_curves <- list()
IPdat_merge_brt_resample_partial_plots <- list()
for (predictor_index in seq_along(IPdat_merge.predictors)) {
  predictor <- IPdat_merge.predictors[predictor_index]
  partial_x <- IPdat_merge_brt_resample_partial_x[, predictor_index, , drop = FALSE]
  partial_prediction <- IPdat_merge_brt_resample_partial_prediction[
    , predictor_index, , drop = FALSE
  ]
  partial_x <- matrix(
    partial_x,
    nrow = IPdat_merge_brt_resample_grid_points,
    ncol = IPdat_merge_brt_resamples_n
  )
  partial_prediction <- matrix(
    partial_prediction,
    nrow = IPdat_merge_brt_resample_grid_points,
    ncol = IPdat_merge_brt_resamples_n
  )

  usable_grid_rows <- rowSums(is.finite(partial_prediction)) > 0L
  if (!any(usable_grid_rows)) {
    warning("no successful partial-dependence curves for ", predictor)
    next
  }

  partial_curve <- data.frame(
    x = apply(partial_x[usable_grid_rows, , drop = FALSE], 1, median, na.rm = TRUE),
    median = apply(partial_prediction[usable_grid_rows, , drop = FALSE], 1,
                   median, na.rm = TRUE),
    lower = apply(partial_prediction[usable_grid_rows, , drop = FALSE], 1,
                  quantile, probs = 0.025, na.rm = TRUE),
    upper = apply(partial_prediction[usable_grid_rows, , drop = FALSE], 1,
                  quantile, probs = 0.975, na.rm = TRUE)
  )
  partial_curve <- partial_curve[order(partial_curve$x), ]
  rownames(partial_curve) <- NULL
  IPdat_merge_brt_resample_partial_curves[[predictor]] <- partial_curve

  write.csv(
    partial_curve,
    file = file.path(
      IPdat_merge_brt_resample_output_dir,
      paste0("InvaPact_BRT_partial_", predictor, ".csv")
    ),
    row.names = FALSE
  )

  partial_plot <- ggplot(partial_curve, aes(x = x, y = median)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.25) +
    geom_line(colour = "steelblue", linewidth = 0.8) +
    labs(
      x = predictor_labels[[predictor]],
      y = "partialled median InvaPact score"
    ) +
    theme_minimal()
  IPdat_merge_brt_resample_partial_plots[[predictor]] <- partial_plot
  print(partial_plot)
}

## print the partial plots for the top 3 predictors by influence
## summary outputs
print(IPdat_merge_brt_resample_cv)
print(IPdat_merge_brt_resample_importance)

IPdat_merge_brt_resample_plot

predictors_by_influence <- IPdat_merge_brt_resample_importance$variable
ggplot_list <- IPdat_merge_brt_resample_partial_plots[predictors_by_influence[1:3]]
partial_plot_grid <- cowplot::plot_grid(plotlist = ggplot_list, ncol = 1)
print(partial_plot_grid)


