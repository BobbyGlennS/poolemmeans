pool_emmeans <- function (fit, formula, data) {
  library(dplyr)
  # set up basic statistics ====================================================
  # n imputations
  m <- length(fit$analyses) 
  
  # n of each cell over which marginal means are estimated
  # !!NEED TO MAKE THIS WORK FOR MULTIPLE FACTORS
  # # extract all predictors from formula
  # factors <- str_extract_all(formula, ".*\\|")
  # factors <- str_split(factors, "[^\\w]", simplify = TRUE)
  # factors <- glue::glue_collapse(factors[str_detect(factors, "\\w")], sep = ", ")
  # get complete data 
  data_complete <- mice::complete(data, action='long', include=TRUE)
  # count n cells
  n_cell <- 
    filter(data_complete, .imp == 1) %>% 
    group_by(condition) %>% 
    summarise(n = n())
  
  # find marginal means for each imputed dataset ===============================
  fit_em <- purrr::map(fit$analyses, 
                       function (x) emmeans::emmeans(x, as.formula(formula)))
  
  # get paramater estimate
  Q_hat <- 
    bind_cols(
      purrr::map(fit_em, 
          function (x) {
            as.data.frame(x)$emmean
          }
      )
    )
  Q_mu <- rowMeans(Q_hat)
  names(Q_mu) <- as.data.frame(fit_em[[1]])[, 1]
   
  # get variance
  # within imputation
  U_hat <- 
    bind_cols(
      map(fit_em, 
          function (x) {
            (as.data.frame(x)$SE * sqrt(n_cell$n))^2
          }
      )
    )
  
  U_mu <- rowMeans(U_hat)
  
  # between imputation
  B <- rowSums((Q_hat - Q_mu)^2) / (m - 1)
  
  # variance
  Q_var <- U_mu + (1 + 1 / m) * B
  Q_SE <- sqrt(Q_var) / sqrt(n_cell$n)
  
  # return output ============================================================
  tibble(Condition = as.data.frame(fit_em[[1]])[, 1], M = Q_mu, SE = Q_SE)
}
