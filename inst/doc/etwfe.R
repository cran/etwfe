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
# mod_es
summary(mod_es) # Summary print looks a bit nicer

## -----------------------------------------------------------------------------
library(modelsummary)

rename_fn = function(old_names) {
  new_names = gsub(".Dtreat", "Years post treatment =", old_names)
  setNames(new_names, old_names)
}

modelsummary(
  mod_es,
  shape = term:event:statistic ~ model,
  coef_rename =  rename_fn,
  gof_omit = "Adj|Within|RMSE",
  title = "Event study",
  notes = "Std. errors are clustered at the county level"
)

## -----------------------------------------------------------------------------
library(ggplot2)

ggplot(mod_es, aes(x = event, y = dydx, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, col = "grey50") +
  geom_pointrange() +
  labs(x = "Years post treatment", y = "Effect on log employment") +
  theme_minimal()

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
    .Dtreat = as.integer(year >= first.treat & first.treat != 0),
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

