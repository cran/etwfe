## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)
fixest::setFixest_notes(FALSE)

## -----------------------------------------------------------------------------
# install.packages("did")
data("mpdta", package = "did")
head(mpdta)

## -----------------------------------------------------------------------------
library(etwfe)

mod =
  etwfe(
    fml  = lemp ~ lpop, # outcome ~ (time-invariant) controls
    tvar = year,        # time variable
    gvar = first.treat, # group variable
    data = mpdta,       # dataset
    vcov = ~countyreal  # vcov adjustment (here: clustered)
    )

## -----------------------------------------------------------------------------
mod

## -----------------------------------------------------------------------------
emfx(mod)

## -----------------------------------------------------------------------------
mod_es = emfx(mod, type = "event")

# Summary print looks a bit nicer, so let's use that here.
summary(mod_es) 

## -----------------------------------------------------------------------------
library(modelsummary)

# Quick renaming function to replace ".Dtreat" with something more meaningful
rename_fn = function(old_names) {
  new_names = gsub(".Dtreat", "Years post treatment =", old_names)
  setNames(new_names, old_names)
}

modelsummary(
  mod_es,
  shape       = term:event:statistic ~ model,
  coef_rename = rename_fn,
  gof_omit    = "Adj|Within|IC|RMSE",
  title       = "Event study",
  notes       = "Std. errors are clustered at the county level"
)

## -----------------------------------------------------------------------------
library(ggplot2)
theme_set(
  theme_minimal() + theme(panel.grid.minor = element_blank())
)

ggplot(mod_es, aes(x = event, y = dydx, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(col = "darkcyan") +
  labs(x = "Years post treatment", y = "Effect on log employment")

## -----------------------------------------------------------------------------
# Use post_only = FALSE to get the "zero" pre-treatment effects
mod_es2 = emfx(mod, type = "event", post_only = FALSE)

ggplot(mod_es2, aes(x = event, y = dydx, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = -1, lty = 2) +
  geom_pointrange(col = "darkcyan") +
  labs(
    x = "Years post treatment", y = "Effect on log employment",
    caption = "Note: Zero pre-treatment effects for illustrative purposes only."
  )

## ---- warning=FALSE, message=FALSE--------------------------------------------
mpdta$emp = exp(mpdta$lemp)

etwfe(
  emp ~ lpop, tvar = year, gvar = first.treat, data = mpdta, vcov = ~countyreal,
  family = "poisson"
  ) |>
  emfx("event") |>
  summary()

## -----------------------------------------------------------------------------
mod$fml_all

## -----------------------------------------------------------------------------
# First construct the dataset
mpdta2 = mpdta |>
  transform(
    .Dtreat = year >= first.treat & first.treat != 0,
    lpop_dm = ave(lpop, first.treat, year, FUN = \(x) x - mean(x, na.rm = TRUE))
  )

# Then estimate the manual version of etwfe
mod2 = fixest::feols(
  lemp ~ .Dtreat:i(first.treat, i.year, ref = 0, ref2 = 2003) / lpop_dm |
    first.treat[lpop] + year[lpop],
  data = mpdta2,
  vcov = ~countyreal
)

## -----------------------------------------------------------------------------
modelsummary(
  list("etwfe" = mod, "manual" = mod2),
  gof_map = NA # drop all goodness-of-fit info for brevity
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# fe = "feo" (fixed effects only)
mod_feo = etwfe(
  lemp ~ lpop, tvar = year, gvar = first.treat, data = mpdta, vcov = ~countyreal,
  fe = "feo"
)
# ... which is equivalent to the manual regression
mod_feo2 = fixest::feols(
  lemp ~ .Dtreat:i(first.treat, i.year, ref = 0, ref2 = 2003) / lpop_dm +
    lpop + i(first.treat, lpop, ref = 0) + i(year, lpop, ref = 2003) |
    first.treat + year,
  data = mpdta2, vcov = ~countyreal
)

# fe = "none"
mod_none = etwfe(
  lemp ~ lpop, tvar = year, gvar = first.treat, data = mpdta, vcov = ~countyreal,
  fe = "none"
)
# ... which is equivalent to the manual regression
mod_none2 = fixest::feols(
  lemp ~ .Dtreat:i(first.treat, i.year, ref = 0, ref2 = 2003) / lpop_dm +
    lpop + i(first.treat, lpop, ref = 0) + i(year, lpop, ref = 2003) +
    i(first.treat, ref = 0) + i(year, ref = 2003),
  data = mpdta2, vcov = ~countyreal
)

## -----------------------------------------------------------------------------
mods = list(
  "etwfe"         = mod,
  "manual"        = mod2,
  "etwfe (feo)"   = mod_feo,
  "manual (feo)"  = mod_feo2,
  "etwfe (none)"  = mod_none,
  "manual (none)" = mod_none2
)

modelsummary(mods, gof_map = NA)

## -----------------------------------------------------------------------------
mod_es_i = etwfe(
  lemp ~ lpop, tvar = year, gvar = first.treat, data = mpdta,
  ivar = countyreal  # NEW: Use unit-level (county) FEs
  ) |>
  emfx("event")

modelsummary(
  list("Group-level FEs (default)" = mod_es, "Unit-level FEs" = mod_es_i),
  shape       = term:event:statistic ~ model,
  coef_rename = rename_fn,
  gof_omit    = "Adj|Within|IC|RMSE"
)

