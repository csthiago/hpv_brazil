set.seed(seed = 1234)
library(tidyverse)
library(brms)
# read data
dados <- read_rds("dados.rds")

# create variables
process_data <- function(data, cancer) {
    data |> 
    filter(DIAG_DETH == {{cancer}}) |> 
    mutate(
      cohort = case_when(
        YEAR_NASC >= 2001 ~ "2001-2003",
        YEAR_NASC >= 2000 ~ "2000",
        YEAR_NASC >= 1999 ~ "1999",
        YEAR_NASC >= 1994 ~ "1994-1998"
      ),
      cohort = forcats::fct_relevel(
        cohort,
        "1994-1998"
      ),
      age_group = case_when(
        IDADE < 23 ~ "20-22",
        IDADE < 25 ~ "23-24",
      ),
      age_group = forcats::fct_relevel(age_group, "20-22"),
      IDADE = as.factor(IDADE),
      IDADE = fct_relevel(IDADE, "22"),
      ANO_DIAGN = as.factor(ANO_DIAGN)
    ) 
}


dados_c53_mod <- process_data(dados, "C53")
dados_d06_mod <- process_data(dados, "D06")
dados_c50_mod <- process_data(dados, "C50")

# Fix cmd_stan for using "update"
fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
options(cmdstanr_write_stan_file_dir = getwd())

# priors
prior1 <- c(
  set_prior("normal(0,5)", class = "b"),
  set_prior("normal(-10,5)", class = "Intercept")
)
mod <- brm(
  n ~ 1+ 
    month_cat+
    ANO_DIAGN+
    age_group+
    cohort + offset(log(value)), #value = population at single age/calendar year
  data = dados_c53_mod, family = "negbinomial",
  prior = prior1,
  chains = 0,
  cores = 4,
  warmup = 2000,
  threads = threading(4),
  control = list(adapt_delta=0.99,
                 max_treedepth = 15),
  iter = 4000,
  seed = 1234
)
fit_model_results <- function(df, cancer) {
  upd_mod <- update(mod,
                    chains = 4,
                    iter = 4000,
                    warmup = 2000,
                    threads = threading(4),
                    seed = 1234,
                    newdata=df)
  parameters <- parameters::model_parameters(upd_mod, exponentiate=T) |> 
    select(Parameter, Component, Median,CI_low, CI_high, pd, Rhat, ESS)
  parameters |> mutate(cancer = cancer)
  
}
# D06 #####
d06_coef <- fit_model_results(dados_d06_mod, "D06")
d06_coef |> clipr::write_clip()
# C53 ####
c53_coef <- fit_model_results(dados_c53_mod, "C53")
c53_coef |> clipr::write_clip()
# C50 ####
c50_coef <- fit_model_results(dados_c50_mod, "C50")
c50_coef |> clipr::write_clip()


# Flat priors ####

mod <- brm(
  n ~ 1+ month_cat+
    ANO_DIAGN+
    age_group+
    cohort + offset(log(value)),
  data = dados_d06_mod, family = "negbinomial",
  chains = 0,
  cores = 4,
  seed = 1234,
  threads = threading(4),
  control = list(adapt_delta=0.99,
                 max_treedepth = 15),
  iter = 4000,
)
fit_model_results <- function(df, cancer) {
  upd_mod <- update(mod,
                    chains = 4,
                    iter = 4000,
                    warmup = 2000,
                    threads = threading(4),
                    seed = 1234,
                    newdata=df)
  parameters <- parameters::model_parameters(upd_mod, exponentiate=T) |> 
    select(Parameter, Component, Median,CI_low, CI_high, pd, Rhat, ESS)
  parameters |> mutate(cancer = cancer)
  
}
# D06 #####
d06_coef <- fit_model_results(dados_d06_mod, "D06")
d06_coef |> clipr::write_clip()
# C53 ####
c53_coef <- fit_model_results(dados_c53_mod, "C53")
c53_coef |> clipr::write_clip()
# C50 ####
c50_coef <- fit_model_results(dados_c50_mod, "C50")
c50_coef |> clipr::write_clip()
