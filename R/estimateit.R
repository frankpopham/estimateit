#' Estimate marginal effects for binary exposure and outcome
#'
#' @param weightitobj A WeightIt object
#' @param data A data frame containing the outcome
#' @param outcome A binary outcome variable
#' @return Returns a summary table as a tibble, the model, and the individual effects
#'    (Control, Treatment, D=their difference, logRR= their log relative "risk",
#'    logOR=their log odds ratio). The individual effects are svrepstat objects that
#'    can be further analysed using relevant survey functions. The table displays the exponential
#'    of the log relative effects. In other words the relative risk and the odds ratio.
#' @importFrom magrittr %>%
#' @export
#' @examples
#' library(WeightIt)
#' library(tibble)
#' library(tidyr)
#' dfvi <- tibble(
#' C = rep(0:1, each = 4),
#' X = rep(0:1, times = 4),
#' Y = rep(0:1, times = 2, each = 2),
#' N = c(96, 36, 64, 54, 120, 120, 30, 480)
#' ) %>%
#'  uncount(N)
#' W1 <- weightit(X ~ C, data = dfvi,
#'               method = "ps", estimand = "ATE")
#' summary(W1)
#' E1 <- estimateit(weightitobj=W1, outcome=Y, data=dfvi)
#' E1


estimateit <- function(weightitobj, outcome, data) {
  df <- dplyr::select(data, Y={{outcome}}) %>%
    dplyr::bind_cols(treat=as.factor(weightitobj$treat)) %>%
    dplyr::bind_cols(w=weightitobj$weights) %>%
    dplyr::bind_cols(scale(tibble::
              as_tibble(stats::model.matrix(stats::as.formula(weightitobj),
              data=data)), scale=F))

forma <- df %>%
  dplyr::select(-c("(Intercept)", Y, treat, w)) %>%
  names()

out_model_ds <- survey::svydesign(id = ~1, data = df, weights = df$w)

out_model <- survey::svyglm(Y ~ treat, family = binomial, design = out_model_ds)

names(out_model$coefficients)[2]  <- "treat1"

outta <- function(model) {
out <- survey::svycontrast(model, list(X0l=c("(Intercept)"=1),
                                           X1l=c("(Intercept)"=1, treat1=1)))
logOR <- survey::svycontrast(out, quote(`X1l` - `X0l`))
Control <- survey::svycontrast(out, quote(exp(`X0l`) / (1 + (exp(`X0l`)))))
Treated <- survey::svycontrast(out, quote(exp(`X1l`) / (1 + (exp(`X1l`)))))
logRR <- survey::svycontrast(out, quote(log(exp(`X1l`) / (1 + (exp(`X1l`)))) -
                                     log(exp(`X0l`) / (1 + (exp(`X0l`))))))
D <- survey::svycontrast(out, quote(exp(`X1l`) / (1 + (exp(`X1l`))) -
                                     exp(`X0l`) / (1 + (exp(`X0l`)))))

list(Control=Control, Treated=Treated, D=D, RR = logRR, OR=logOR)
}

marginal <- outta(out_model)


#no tidy for svyrep - this based on https://www.tidymodels.org/learn/develop/broom/
tidyup <- function(x, conf.int = TRUE, conf.level = 0.95) {
   result <- tibble::as_tibble(x) %>%
     dplyr::select(-SE)
     colnames(result) <- c("estimate")
  if (conf.int) {
    ci <- tibble::as_tibble(confint(x, level = conf.level))
    colnames(ci) <- c("conf.low", "conf.high")
    result <- dplyr::bind_cols(result, ci)
  }
result
}

tabletime <- function(alist) {
  purrr::map_dfr(alist, ~ tidyup(.x), .id="effect") %>%
  dplyr::mutate(dplyr::across(c(estimate, conf.low, conf.high),
                              ~ dplyr::if_else(effect == "RR" | effect =="OR", exp(.x), .x))) %>%
  dplyr::mutate(estimand = weightitobj$estimand) %>%
  dplyr::select(estimand, tidyselect::everything())
}

marginal2 <- tabletime(marginal)


return(list(marginal=list("Outcome model"=out_model, effects=marginal,
                            "Effects table"=marginal2)))
}

