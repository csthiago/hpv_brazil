library(tidyverse)
# https://www.gov.br/inca/pt-br/assuntos/gestor-e-profissional-de-saude/controle-do-cancer-do-colo-do-utero/acoes/deteccao-precoce
# "dados.rds" contains all combinations found in the raw data, with the female population by age/ diagnosis year joined
dados_diag <- read_rds("dados.rds") |> 
  mutate(
    cohort = case_when(
      YEAR_NASC >= 2001 ~ "2001-2003",
      YEAR_NASC >= 2000 ~ "2000",
      YEAR_NASC >= 1999 ~ "1999",
      YEAR_NASC >= 1994 ~ "1994-1998",
      YEAR_NASC >= 1980 ~ "1980-1989",
      YEAR_NASC >= 1970 ~ "1970-1979",
      YEAR_NASC >= 1945 ~ "1945-1959"
    ),
    cohort = forcats::fct_relevel(
      cohort,
      "1994-1998"
    ),
    age_group = case_when(
      IDADE < 23 ~ "20-22",
      IDADE < 25 ~ "23-24",
    ),
    age_group = forcats::fct_relevel(age_group, "20-22"))
dados_diag <- dados_diag|> 
  mutate(YEAR_NASC=as.factor(YEAR_NASC),
         IDADE=as.factor(IDADE),
         ANO_DIAGN=as.factor(ANO_DIAGN),
         month_cat=as.factor(month_cat))
dados_c53_mod <- dados_diag |> filter(DIAG_DETH=="C53")
dados_c50_mod <- dados_diag |> filter(DIAG_DETH=="C50")
dados_d06_mod <- dados_diag |> filter(DIAG_DETH=="D06")

# Fix combinations with no cases
## get the birth year range
cohort_range <- range(as.numeric(as.character(dados_c53_mod$YEAR_NASC)))

valid <- expand_grid(
  IDADE     = sort(unique(as.numeric(as.character(dados_c53_mod$IDADE)))),
  ANO_DIAGN = sort(unique(as.numeric(as.character(dados_c53_mod$ANO_DIAGN)))),
  tri       = 0:1                     # older / younger in that month/year;
) |>
  mutate(YEAR_NASC = ANO_DIAGN - IDADE - tri) |>
  filter(between(YEAR_NASC, cohort_range[1], cohort_range[2])) |>
  select(-tri) |>
  expand_grid(month_cat = 1:4)

# get all valid combinations for the data 50 combinations of age/year birth/diagnosis. e.g. in 2019, those born in 1994 can only enter with diagnosis month before birthday month - similar to 2003 and 2023 (only diagnosis month after birthday month), other combinations allows both, in those 50 comb there is 4 trimester each = 200 possible cells
valid <- valid|> 
  mutate(YEAR_NASC=as.factor(YEAR_NASC),
         IDADE=as.factor(IDADE),
         ANO_DIAGN=as.factor(ANO_DIAGN),
         month_cat=as.factor(month_cat))
# get those without any case in the combination age/diagnosis year/trimester
missing <- left_join(valid, dados_c53_mod,
                     by = c("YEAR_NASC", "IDADE", "ANO_DIAGN", "month_cat"))

# replace with 0 the combination with missing cases and add the female population in that cell
dados_c53_mod <- missing |> group_by(YEAR_NASC, IDADE, ANO_DIAGN) |> 
  fill(c(value,cohort,age_group), .direction = "downup") |> 
  mutate(n=replace_na(n,0))

missing <- left_join(valid, dados_d06_mod,
                     by = c("YEAR_NASC", "IDADE", "ANO_DIAGN", "month_cat"))

dados_d06_mod <- missing |> group_by(YEAR_NASC, IDADE, ANO_DIAGN) |> 
  fill(c(value,cohort,age_group), .direction = "downup") |> 
  mutate(n=replace_na(n,0))

missing <- left_join(valid, dados_c50_mod,
                     by = c("YEAR_NASC", "IDADE", "ANO_DIAGN", "month_cat"))

dados_c50_mod <- missing |> group_by(YEAR_NASC, IDADE, ANO_DIAGN) |> 
  fill(c(value,cohort,age_group), .direction = "downup") |> 
  mutate(n=replace_na(n,0))


# Fix cmd_stan
fname <- paste0("fit_cmdstanr_", sample.int(.Machine$integer.max, 1))
options(cmdstanr_write_stan_file_dir = getwd())
library(brms)
prior1 <- c(
  set_prior("normal(0,5)", class = "b"),
  set_prior("normal(-10,5)", class = "Intercept")
)
set.seed(seed = 123)
mod <- brm(
  n ~ 1+ 
    month_cat+
    ANO_DIAGN+
    age_group+
    cohort + offset(log(value)),
  data = dados_c53_mod, family = "negbinomial",
  prior = prior1,
  chains = 0,
  cores = 2,
  control = list(adapt_delta=0.99,
                 max_treedepth = 15),
  warmup = 2000,
  iter = 6000,
)
mod_flat <- brm(
  n ~ 1+ 
    month_cat+
    ANO_DIAGN+
    age_group+
    cohort + offset(log(value)),
  data = dados_c53_mod, family = "negbinomial",
  chains = 0,
  cores = 2,
  control = list(adapt_delta=0.99,
                 max_treedepth = 15),
  warmup = 2000,
  iter = 6000,
)
fit_model_results <- function(df, cancer) {
  upd_mod <- update(mod,
                    chains = 4,
                    newdata=df,
                    threads = threading(4),
                    seed = 1234,
                    iter=6000,
                    warmup=2000,
                    control = list(adapt_delta=0.99,
                                   max_treedepth = 15))
  parameters <- parameters::model_parameters(upd_mod, exponentiate=T) |> 
    select(Parameter, Median,CI_low, CI_high, pd, Rhat, ESS)
  parameters |> mutate(cancer = cancer)
  
}
fit_model_flat_results <- function(df, cancer) {
  upd_mod <- update(mod_flat,
                    chains = 4,
                    newdata=df,
                    threads = threading(4),
                    seed = 1234,
                    iter=6000,
                    warmup=2000,
                    control = list(adapt_delta=0.99,
                                   max_treedepth = 15))
  parameters <- parameters::model_parameters(upd_mod, exponentiate=T) |> 
    select(Parameter, Median,CI_low, CI_high, pd, Rhat, ESS)
  parameters |> mutate(cancer = cancer)
  
}
format_results <- function(mod){
  mod |> mutate(
    across(Median:CI_high,\(x)round(x,2)),
    fmt = paste0(Median," (",CI_low," - ",CI_high,")")
  )
}
# D06 #####

d06_coef <- fit_model_results(dados_d06_mod, "D06")
format_results(d06_coef) |> clipr::write_clip()
# C53 ####
c53_coef <- fit_model_results(dados_c53_mod, "C53")
format_results(c53_coef) |> clipr::write_clip()
# C50 ####
c50_coef <- fit_model_results(dados_c50_mod, "C50")
format_results(c50_coef) |> clipr::write_clip()


# Flat priors
d06_coef <- fit_model_flat_results(dados_d06_mod, "D06")
format_results(d06_coef) |> clipr::write_clip()
# C53 ####
c53_coef <- fit_model_flat_results(dados_c53_mod, "C53")
format_results(c53_coef) |> clipr::write_clip()
# C50 ####
c50_coef <- fit_model_flat_results(dados_c50_mod, "C50")
format_results(c50_coef) |> clipr::write_clip()

# Incidence ####

denom <- dados_c53_mod |> ungroup() |> 
  distinct(ANO_DIAGN, IDADE, value) |>
  group_by(ANO_DIAGN) |>
  summarise(wyears = sum(value), .groups = "drop")
denom|> clipr::write_clip()

denom_ag <- dados_c53_mod |> ungroup() |> 
  distinct(age_group, ANO_DIAGN, IDADE, value) |>
  summarise(.by = age_group, wyears = sum(value))
denom_ag|> clipr::write_clip()


denom_coh <- dados_c53_mod |> ungroup() |> 
  distinct(YEAR_NASC, cohort, IDADE, ANO_DIAGN, value) |>
  mutate(pyears = value / 2) |>                 # one Lexis triangle per YEAR_NASC
  summarise(.by = cohort, wyears = sum(pyears))
denom_coh|> clipr::write_clip()
# c53 ####
inc_by_year <- dados_c53_mod |>
  group_by(ANO_DIAGN) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(denom, by = "ANO_DIAGN")
inc_by_age <- dados_c53_mod |>
  group_by(age_group) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(denom_ag, by = "age_group")
inc_by_cohort <- dados_c53_mod |>
  group_by(cohort) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(denom_coh, by = "cohort")

mod_inc <- brm(n ~ 1 + offset(log(wyears / 100000)),
               data = inc_by_age, family = "poisson",
               warmup = 2000,
               iter = 6000
)
format_results(parameters::model_parameters(mod_inc, exponentiate = T)) |> clipr::write_clip()



mod_inc <- brm(n ~ -1 + ANO_DIAGN + offset(log(wyears / 100000)), 
               data = inc_by_year, family = "poisson",
               warmup = 2000,
               iter = 6000)
format_results(parameters::model_parameters(mod_inc, exponentiate = T))  |> clipr::write_clip()

mod_inc <- brm(n ~ -1 + age_group + offset(log(wyears / 100000)), 
               data = inc_by_age, family = "poisson",
               warmup = 2000,
               iter = 6000)
format_results(parameters::model_parameters(mod_inc, exponentiate = T)) |> clipr::write_clip()

mod_inc <- brm(n ~ -1 + cohort + offset(log(wyears / 100000)), 
               data = inc_by_cohort, family = "poisson",
               warmup = 2000,
               iter = 6000)
format_results(parameters::model_parameters(mod_inc, exponentiate = T))  |> clipr::write_clip()


# d06 ####
inc_by_year <- dados_d06_mod |>
  group_by(ANO_DIAGN) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(denom, by = "ANO_DIAGN")
inc_by_age <- dados_d06_mod |>
  group_by(age_group) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(denom_ag, by = "age_group")
inc_by_cohort <- dados_d06_mod |>
  group_by(cohort) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(denom_coh, by = "cohort")

mod_inc <- brm(n ~ 1 + offset(log(wyears / 100000)),
               data = inc_by_age, family = "poisson",
               warmup = 2000,
               iter = 6000
)
format_results(parameters::model_parameters(mod_inc, exponentiate = T)) |> clipr::write_clip()

mod_inc <- brm(n ~ -1 + ANO_DIAGN + offset(log(wyears / 100000)), 
               data = inc_by_year, family = "poisson",
               warmup = 2000,
               iter = 6000)
format_results(parameters::model_parameters(mod_inc, exponentiate = T))  |> clipr::write_clip()

mod_inc <- brm(n ~ -1 + age_group + offset(log(wyears / 100000)), 
               data = inc_by_age, family = "poisson",
               warmup = 2000,
               iter = 6000)
format_results(parameters::model_parameters(mod_inc, exponentiate = T)) |> clipr::write_clip()

mod_inc <- brm(n ~ -1 + cohort + offset(log(wyears / 100000)), 
               data = inc_by_cohort, family = "poisson",
               warmup = 2000,
               iter = 6000)
format_results(parameters::model_parameters(mod_inc, exponentiate = T))  |> clipr::write_clip()



# Figure 2 ####

bind_rows(d06_coef,
          c53_coef,
          c50_coef) |> 
  filter(Parameter != "shape",
         Parameter != "b_Intercept",
         str_detect(Parameter, "b_cohort")) |> 
  mutate(Parameter = fct_recode(Parameter,
                                "1999" = "b_cohort1999",
                                "2000" = "b_cohort2000",
                                "2001-2003"="b_cohort2001M2003")) |> 
  add_row(Parameter = "1994-1998\nReference", Median = 1, cancer = "ref") |> 
  ggplot(aes(x=Parameter,y=Median))+
  geom_hline(aes(yintercept = 1), linetype=2)+
  scale_y_continuous(breaks = c(0.2,0.4,0.6,.8,1,1.2,1.4),
                     limits = c(.1,1.5))+
  ggstats::geom_stripped_cols()+
  labs(y="Incidence Rate Ratio (95% CredI)", x="Birth Cohort", color="Cancer")+
  scale_x_discrete(limits=rev)+
  geom_vline(aes(xintercept=c(2.5)), linetype=3)+
  geom_vline(aes(xintercept=c(3.5)), linetype=3)+
  geom_pointrange(aes(color = cancer,
                      ymin = CI_low, ymax=CI_high),
                  position = position_dodge(width=1),
                  size=1,
                  linewidth=1)+
  coord_flip()+
  scale_color_manual(values = c("C53"="#DD5129FF",
                                "C50"="#0F7BA2FF",
                                "D06"="#FAB255FF"),
                     breaks = c("C53","D06","C50"),
                     labels = c("C53"="Cervical",
                                "C50"="Breast",
                                "D06"="CIN3"))+
  firatheme::theme_fira()+
  theme(legend.position = "bottom")+
  guides(color=guide_legend(override.aes = list(linetype=0),
                            title.position = "left",
                            label.position = "bottom"))

# Flat priors ####



# D06 #####
d06_coef <- fit_model_results(dados_d06_mod, "D06")
d06_coef
# C53 ####
c53_coef <- fit_model_results(dados_c53_mod, "C53")
c53_coef
# C50 ####
c50_coef <- fit_model_results(dados_c50_mod, "C50")
c50_coef


# Table ####
library(gtsummary)

dados_c53_mod |> 
  uncount(weights = n) |> 
  select(age_group,
         ANO_DIAGN,
         cohort) |> 
  tbl_summary(digits = all_categorical()~1)
process_data_tbl(dados_d06) |> 
  select(age_group,
         ANO_DIAGN,
         cohort) |> 
  tbl_summary(digits = all_categorical()~1)

process_data_tbl(dados_c50) |> 
  select(age_group,
         ANO_DIAGN,
         cohort) |> 
  tbl_summary(digits = all_categorical()~1)

# AIH #####

aih_data <- read_rds("aih_data.rds")
aih_data <- aih_data|> 
  mutate(YEAR_NASC=as.factor(YEAR_NASC),
         IDADE=as.factor(IDADE),
         ANO_DIAGN=as.factor(ANO_DIAGN),
         month_cat=as.factor(month_cat))
dados_c53_mod <- aih_data |> filter(cancer=="c53")
dados_c50_mod <- aih_data |> filter(cancer=="c50")
dados_d06_mod <- aih_data |> filter(cancer=="d06")

# structurally possible, but absent from the data = incidental zeros, dropped
missing <- left_join(valid, dados_c53_mod,
                     by = c("YEAR_NASC", "IDADE", "ANO_DIAGN", "month_cat"))

# replace with 0 the combination with missing cases
dados_c53_mod <- missing |> group_by(YEAR_NASC, IDADE, ANO_DIAGN) |> 
  fill(c(value,cohort,age_group), .direction = "downup") |> 
  mutate(n=replace_na(n,0))

missing <- left_join(valid, dados_d06_mod,
                     by = c("YEAR_NASC", "IDADE", "ANO_DIAGN", "month_cat"))

dados_d06_mod <- missing |> group_by(YEAR_NASC, IDADE, ANO_DIAGN) |> 
  fill(c(value,cohort,age_group), .direction = "downup") |> 
  mutate(n=replace_na(n,0))

missing <- left_join(valid, dados_c50_mod,
                     by = c("YEAR_NASC", "IDADE", "ANO_DIAGN", "month_cat"))

dados_c50_mod <- missing |> group_by(YEAR_NASC, IDADE, ANO_DIAGN) |> 
  fill(c(value,cohort,age_group), .direction = "downup") |> 
  mutate(n=replace_na(n,0))


d06_coef <- fit_model_results(dados_d06_mod, "D06")
format_results(d06_coef) |> clipr::write_clip()
# C53 ####
c53_coef <- fit_model_results(dados_c53_mod, "C53")
format_results(c53_coef) |> clipr::write_clip()
# C50 ####
c50_coef <- fit_model_results(dados_c50_mod, "C50")
format_results(c50_coef) |> clipr::write_clip()


# Flat priors AIH
d06_coef <- fit_model_flat_results(dados_d06_mod, "D06")
format_results(d06_coef) |> clipr::write_clip()
# C53 ####
c53_coef <- fit_model_flat_results(dados_c53_mod, "C53")
format_results(c53_coef) |> clipr::write_clip()
# C50 ####
c50_coef <- fit_model_flat_results(dados_c50_mod, "C50")
format_results(c50_coef) |> clipr::write_clip()

