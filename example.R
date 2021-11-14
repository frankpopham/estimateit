wrt <-map(list(ATE="ATE", ATT="ATT", ATC="ATC", AT0="ATO", ATM="ATM", ATOS="ATOS"),
    ~weightit(formula = X ~ C, data = df, method = "ps", estimand = .x))

test45 <- map(wrt, ~estimateit(.x, Y, df))

test46 <- map(test45, ~.x$Summary) %>%
  map_dfr(~.x, .id="Estimand")
