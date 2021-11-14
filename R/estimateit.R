#' Estimate marginal effects for binary exposure and outcome
#'
#' @param weightitobj A Weightit object
#' @param df A data frame containing the variables
#' @param outcome A binary outcome variable coded 0/1
#' @return Returns a Summary tibble, and the individual effects
#' @importFrom magrittr %>%
#' @examples
#' library(cobalt)
#' library(WeightIt)
#'
#' W1 <- weightit(treat ~ age + educ + married +
#'                 re74, data = lalonde,
#'               method = "ps", estimand = "ATT")
#' summary(W1)
#' bal.tab(W1)
#' E1 <- estimateit(W1, nodegree, lalonde)
#' E1

estimateit <- function(weightitobj, outcome, data) {
  df <- select(data, Y={{outcome}}) %>%
    dplyr::bind_cols(X=weightitobj$treat) %>%
    dplyr::bind_cols(w=weightitobj$weights)

out_model_ds <- survey::svydesign(id = ~1, data = df, weights = df$w)

out_model <- survey::svyglm(Y ~ factor(X), family = binomial, design = out_model_ds)

sum_model <- broom::tidy(out_model, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::filter(term == "factor(X)1") %>%
    dplyr::mutate(Effect = "OR") %>%
    dplyr::select(Effect, Estimate = estimate, Conf.low=conf.low, Conf.high=conf.high)



X0 <- survey::svycontrast(out_model, quote(exp(`(Intercept)`) / (1 + exp(`(Intercept)`))))
X1 <- survey::svycontrast(out_model, quote((exp(`(Intercept)` + `factor(X)1`))
    / (1 + (exp(`(Intercept)` + `factor(X)1`)))))
rr <- survey::svycontrast(out_model, quote(log((exp(`(Intercept)` + `factor(X)1`))
    / (1 + (exp(`(Intercept)` + `factor(X)1`)))) -
      log((exp(`(Intercept)`) / (1 + exp(`(Intercept)`))))))
rd <- survey::svycontrast(out_model, quote(((exp(`(Intercept)` + `factor(X)1`))
    / (1 + (exp(`(Intercept)` + `factor(X)1`)))) -
      (exp(`(Intercept)`) / (1 + exp(`(Intercept)`)))))


summary <- list(X0 = X0, X1 = X1, D = rd, RR = rr)

contrasts_ci <- purrr::map(summary, confint) %>%
    purrr::map_dfr(~ tibble::as_tibble(.x), .id = "Effect")

contrasts <- purrr::map_dfr(summary, ~ tibble::as_tibble(.x), .id="Effect") %>%
  dplyr::bind_cols(select(contrasts_ci, -Effect)) %>%
  dplyr::select(Effect, Estimate = nlcon, Conf.low = "2.5 %", Conf.high = "97.5 %") %>%
  dplyr::mutate(dplyr::across(c(Estimate, Conf.low, Conf.high),
                              ~ dplyr::if_else(Effect == "RR", exp(.x), .x))) %>%
  dplyr::bind_rows(sum_model) %>%
  dplyr::mutate(Estimand = weightitobj$estimand) %>%
  dplyr::select(Estimand, everything())
return(list(Outcome_model=out_model, X0 = X0, X1 = X1, D = rd, RR = rr, Table=contrasts))
}
