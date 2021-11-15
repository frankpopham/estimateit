#' Estimate marginal effects for binary exposure and outcome
#'
#' @param weightitobj A Weightit object
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
#' library(cobalt)
#' library(WeightIt)
#' data("lalonde", package = "cobalt")
#' W1 <- weightit(treat ~ age + educ + married +
#'                 re74, data = lalonde,
#'               method = "ps", estimand = "ATE")
#' summary(W1)
#' bal.tab(W1)
#' E1 <- estimateit(W1, nodegree, lalonde)
#' E1



estimateit <- function(weightitobj, outcome, data) {
  df <- dplyr::select(data, Y={{outcome}}) %>%
    dplyr::bind_cols(treat=weightitobj$treat) %>%
    dplyr::bind_cols(w=weightitobj$weights)

out_model_ds <- survey::svydesign(id = ~1, data = df, weights = df$w)

out_model <- survey::svyglm(Y ~ treat, family = binomial, design = out_model_ds)

out <- survey::svycontrast(out_model, list(X0l=c(1,0), X1l=c(1,1)))


logOR <- survey::svycontrast(out, quote(`X1l` - `X0l`))
Control <- survey::svycontrast(out, quote(exp(`X0l`) / (1 + (exp(`X0l`)))))
Treated <- survey::svycontrast(out, quote(exp(`X1l`) / (1 + (exp(`X1l`)))))
logRR <- survey::svycontrast(out, quote(log(exp(`X1l`) / (1 + (exp(`X1l`)))) -
                                     log(exp(`X0l`) / (1 + (exp(`X0l`))))))
D <- survey::svycontrast(out, quote(exp(`X1l`) / (1 + (exp(`X1l`))) -
                                     exp(`X0l`) / (1 + (exp(`X0l`)))))




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

summarylist <- list(Control=Control, Treated=Treated, D=D, RR = logRR, OR=logOR)

contrasts <- purrr::map_dfr(summarylist, ~ tidyup(.x), .id="effect") %>%
  dplyr::mutate(dplyr::across(c(estimate, conf.low, conf.high),
                              ~ dplyr::if_else(effect == "RR" | effect =="OR", exp(.x), .x))) %>%
  dplyr::mutate(estimand = weightitobj$estimand) %>%
  dplyr::select(estimand, tidyselect::everything())

return(list(Outcome_model=out_model, Control=Control, Treated=Treated, D=D,
            logRR=logRR, logOR=logOR, Table=contrasts))
}
