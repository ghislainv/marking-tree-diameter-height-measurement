# Loading libraries
library(readr)
library(here)
library(ggplot2)
library(lme4)
library(dplyr)
library(nlme)

# Read data-sets
df <- read_csv(here("dbh-measurement-error.csv"))

# Function to get tree id
get_tree_id <- function(x) {
  id <- substr(x, 3, 5)
  id <- paste0("t", id)
  return(id)
}

# Manipulate variables
df <- df |> mutate(tree_id=get_tree_id(tree))

# Model
mod_1 <- lmer(dbh ~ -1 + tree_id, data=df)
ranef(mod_1)

# Residual variance
var(residuals(mod_1))

# Plot obs vs. mean
df <- df |> mutate(dbh_mean=rep(fixef(mod_1), 5))
p <- ggplot(data=df, aes(dbh_mean, dbh, color=pom)) +
  geom_point() +
  theme_bw(base_size=14)
ggsave("obs_mean_first.pdf")

# Strange point
df |> dplyr::filter(dbh_mean < 15 & dbh > 15)

# Replace strange point
df2 <- df |>
  mutate(dbh=ifelse((group == "gp3") & (tree_id == "t235"), 10.4, dbh))

# Model
mod_2 <- lmer(dbh ~ -1 + tree_id + (1 | grp_name), data=df2)
df2 <- df2 |> mutate(dbh_mean=rep(fixef(mod_2), 5))
# Plot
p <- ggplot(data=df2, aes(dbh_mean, dbh, color=pom)) +
  geom_point() +
  theme_bw(base_size=15)

# lme model with different variances
# https://fukamilab.github.io/BIO202/03-C-heterogeneity.html
mod_3 <- lme(dbh ~ -1 + tree_id,
             random=~1|grp_name,
             weights=varIdent(form=~1|pom),
             na.action=na.exclude, 
             data=df2)

mod_4 <- lme(dbh ~ -1 + tree_id,
             random=~1|grp_name,
             weights=varFixed(~dbh),
             na.action=na.exclude, 
             data=df2)

# ==================
# Residuals analysis
# ==================

# Function to compute the standard error
std_err <- function(x) {
  sd(x)/sqrt(length(x))
}

# Model
mod_res <- lm(dbh ~ -1 + tree_id, data=df2)
# Pred and residuals
df2 <- df2 |>
  mutate(res=residuals(mod_res)) |>
  mutate(dbh_mean=rep(coefficients(mod_res), 5))
# Plot
p <- ggplot(data=df2, aes(dbh_mean, dbh, color=pom)) +
  geom_point() +
  theme_bw(base_size=15)
ggsave("obs_mean.pdf")

# Histogram plot
p <- ggplot(data=df2, aes(res, fill=pom)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position="identity") +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  theme_bw(base_size=15) + labs(fill="")
ggsave("histogram.pdf")
  
# Residual variance per pom
# pom_yes
res_pom_yes <- df2 |>
  dplyr::filter(pom == "pom_yes") |>
  select(res) |>
  pull()
var_pom_yes <- var(res_pom_yes)
se_pom_yes <- std_err(res_pom_yes) # 0.01326107, 0.13 mm error on average
# pom_no
res_pom_no <- df2 |>
  dplyr::filter(pom == "pom_no") |>
  select(res) |>
  pull()
var_pom_no <- var(res_pom_no)
se_pom_no <- std_err(res_pom_no) # 0.0620643, 0.62 mm error on average

# Residual variance for each tree
df3 <- df2 |>
  group_by(tree_id) |>
  summarise(std_err=std_err(res), dbh_mean=unique(dbh_mean)) |>
  ungroup()
mod_str_err_dbh_mean <- lm(std_err~dbh_mean, data=df3)
summary(mod_str_err_dbh_mean)
p <- ggplot(data=df3, aes(dbh_mean, std_err)) +
  geom_point() +
  geom_smooth(method="lm", alpha=.15) +
  theme_bw(base_size=15)
ggsave("std_err_dbh_mean.pdf")
  
# End of file
