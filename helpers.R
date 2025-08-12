library(ggplot2)
library(marquee)

# colors
int_col <- "#B02A3B"
lightgray <- "gray60"
darkgray <- "gray30"
white <- "white"
creme <- "#f5ebe0"
lightblack <- "#212529"

title_style <- modify_style(
  classic_style(),
  "int_highlight",
  background = int_col,
  color = white,
  border_size = trbl(em(0.1)),
  border_radius = 3,
  padding = trbl(em(0.08))
)

title_style <- modify_style(
  title_style,
  "cont_highlight",
  color = white,
  background = darkgray,
  border_size = trbl(em(0.1)),
  border_radius = 3,
  padding = trbl(em(0.15))
)

title_style <- modify_style(
  title_style,
  "int_col",
  color = int_col,
)

title_style <- modify_style(
  title_style,
  "cont_col",
  color = darkgray,
)

custom_theme <- function(gridlines = c("y", "x", "both", "scatter"), base_size = 12,
                         family = "Roboto Condensed", title_family = "Roboto", margins = TRUE,
                         plot.title.position = "plot", axis_titles = TRUE, multiplot = FALSE,
                         ...) {
  grd <- match.arg(gridlines)
  element_gridline <- element_line(colour = creme, size = 0.3)
  base_theme <- theme_minimal(base_size = base_size, base_family = family)
  custom_theme <- theme(
    plot.title.position = plot.title.position,
    text = element_text(colour = lightblack),
    plot.title = element_marquee(
      size = base_size * 1.2, family = title_family, lineheight = 1.1,
      width = 1,
      margin = margin(b = 6),
      style = title_style
    ), 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = if (grd != "y") {
      element_gridline
    } else {
      element_blank()
    }, 
    panel.grid.major.y = if (grd != "x") {
      element_gridline
    } else {
      element_blank()
    }, 
    panel.background = element_rect(
      fill = white,
      colour = NA
    ),
    plot.background = element_rect(
      fill = white,
      colour = NA
    ),
    axis.title = if (axis_titles) {
      element_text(colour = darkgray)
    } else {
      element_blank()
    }, 
    strip.text = element_text(hjust = 0, colour = lightblack),
    plot.margin = if (margins) {
      margin(l = 1, b = 1, t = 4)
    } else {
      base_theme$plot.margin
    }, 
    strip.background = if (multiplot) {
      element_rect(fill = white, colour = NA)
    } else {
      element_blank()
    }, 
    plot.subtitle = element_marquee(
      family = title_family,
      lineheight = 1.1, margin = margin(b = 6),
      width = 1,
      colour = darkgray,
      style = title_style
    ), 
    plot.caption = element_marquee(
      colour = lightgray,
      lineheight = 1, margin = margin(t = 5)
    ),
    axis.text = element_text(colour = darkgray),
    ...
  )
  base_theme + custom_theme
}


# IRT ---------------------------------------------------------------------

vic_codebook <- tibble::tribble(
  ~item, ~dim, ~rev,
  "vic_1", "vic", FALSE, 
  "vic_2", "vic", FALSE,
  "vic_3", "vic", FALSE,
  "vic_4", "vic", FALSE,
  "vic_5", "vic", FALSE,
  "vic_6", "vic", FALSE,
  "vic_7", "vic", FALSE,
  "vic_8", "vic", FALSE,
  "vic_9", "vic", FALSE,
  "vic_10", "vic", FALSE
)

make_mirt_data <- function(.data, codebook, na_level = NULL,
                           ignore_case = TRUE, negate = FALSE, plot = TRUE,
                           check_lvls = FALSE) {
  # make a "backup" of original factor variables for barplots etc.
  .data <- .data |>
    mutate(across(all_of(codebook$item), .names = ".{.col}"))

  common_lvls <- .data %>%
    select(all_of(codebook$item)) |>
    map(levels) |>
    unique()

  if (length(common_lvls) != 1L) rlang::abort("levels of items differ")

  common_lvls <- common_lvls[[1L]]


  if (!is.null(na_level)) {
    match <- str_subset(common_lvls, regex(na_level, ignore_case = ignore_case),
                        negate = negate
    )

    .data <- .data |>
      mutate(
        across(
          all_of(codebook$item),
          \(item) fct_nanify(item,
                             level = na_level,
                             ignore_case = ignore_case, negate = negate
          ) |>
            fct_drop(match) # remove the level after data are removed
        )
      )
  }
  new_common_lvls <- .data %>%
    select(all_of(codebook$item)) |>
    map(levels) |>
    unique()

  if (length(new_common_lvls) != 1L) rlang::abort("new levels of items differ")

  new_common_lvls <- new_common_lvls[[1L]]

  if(check_lvls) {
    message("check new integers (for nonreversed items only):")
    names(new_common_lvls) <- seq_along(new_common_lvls)
    
    print(new_common_lvls)
}

  # now reversed items
  rev_items <- codebook |>
    filter(rev) |>
    pull(item)
  max_val <- length(new_common_lvls)

  .data <- .data |>
    mutate(
      across(all_of(codebook$item), \(item) as.integer(item)),
      across(any_of(rev_items), \(item) (max_val + 1L) - item)
    )

  if (plot) {
    print(.data |>
      select(any_of(codebook$item)) |>
      psych::polychoric() |>
      magrittr::extract2("rho") |>
      as.data.frame() |> 
      rownames_to_column("var1") |> 
      pivot_longer(-var1, names_to = "var2", values_to = "corr") |>
      mutate(across(starts_with("var"), 
                    \(x) factor(x, levels = vic_codebook$item))) |> 
      ggplot(aes(x = var1, y = var2, fill = corr)) +
      geom_tile() +
      scale_fill_gradient2(
        low = "#E53935",
        mid = "#F5F5F5",
        high = "#1E88E5",
        midpoint = 0,
        limits = c(-1, 1)
      ) +
      scale_y_discrete(name = NULL, limits = rev) +
      scale_x_discrete(name = NULL, position = "top") +
      custom_theme() +
      theme(axis.text.x.top = element_marquee(angle = 90, vjust = 0.5)))
  }

  .data
}


# fun for building mirt model out of codebook
mod_from_codebook <- function(codebook, cov_terms = TRUE) {
  out <- codebook %>%
    mutate(id = row_number(), dim = fct_inorder(dim)) %>%
    group_by(dim) %>%
    summarise(item = str_flatten(id, collapse = ",")) %>%
    glue::glue_data("{dim} = {item}") %>%
    str_flatten(collapse = "\n")

  if (cov_terms) {
    cov_terms <- codebook$dim %>%
      unique() %>%
      str_flatten(collapse = "*")

    out <- str_c(out, "\nCOV = ", cov_terms)
  }

  out
}


# fun for fitting GRM on W1 data accoring to codebook and model specified above
fit_grm <- function(.data, codebook, drop_empty = TRUE, method = "MHRM", ...) {
  mod <- codebook %>% mod_from_codebook()

  # ensure the order is the same as in the codebook which defines
  # item indices mirt acts upon
  .data <- .data %>% select(all_of(codebook$item))

  if (drop_empty) {
    # refactor: .data |> drop_na()
    .data <- .data |>
      filter(if_any(everything(), ~ !is.na(.x)))
  }

  .data |>
    mirt(mod, itemtype = "graded", method = method, ...)
}



bind_thetas <- function(.data, fit_object, codebook, restore_items = TRUE, plausible_draw = FALSE, method = NULL, ...) {
  resp_patt <- .data %>% select(all_of(codebook$item))

  # if plausible vals, mirt return raw matrix without names
  factor_names <- fit_object |> extract.mirt("factorNames")

  if (is.null(method)) {
    method <- if (length(factor_names) > 1L) "MAP" else "EAP"
    method <- if (plausible_draw) "plausible" else method
  }

  qmc <- if (length(factor_names) > 1L) TRUE else FALSE

  set.seed(123)
  res <- fscores(fit_object,
                 response.pattern = resp_patt, method = method, QMC = qmc, ...
  )


  # we have to treat the output diferently if using PVs (its a mistake of mirt to some extent)
  if (plausible_draw) {
    # get indices of empty rows
    empty_case <- rowMeans(is.na(resp_patt)) == 1

    # make a form to cast the results - the same dims as the resp_patt
    out <- matrix(nrow = nrow(resp_patt), ncol = length(factor_names))

    # fill non-empty rows with the results - other rows are NAs, as they would be in normal fscores output
    out[!empty_case, ] <- res

    # name cols
    colnames(out) <- str_c(factor_names, "_pv")
  }

  # restore original variables, as the mirt ones are not needed anymore
  # (we have the thetas)
  if (restore_items) {
    .data <- .data |>
      select(-any_of(codebook$item)) |> # drop mirt vars
      rename_with(
        \(item) str_remove(item, "\\."),
        any_of(str_c(".", codebook$item))
      )
  }

  if (plausible_draw) {
    bind_cols(.data, out)
  } else {
    bind_cols(.data, res)
  }
}


# hurdle-gaussian ---------------------------------------------------------

# the custom hurdle_gaussian fam is copied from Andrew Heiss's blog
# https://www.andrewheiss.com/blog/2022/05/09/hurdle-lognormal-gaussian-brms/
hurdle_gaussian <-
  # Create a custom family that is logit if y = 0, normal/gaussian if not
  custom_family(
    "hurdle_gaussian",
    dpars = c("mu", "sigma", "hu"),
    links = c("identity", "log", "logit"),
    lb = c(NA, 0, NA),
    type = "real"
  )

# Stan code
stan_funs <- "
  real hurdle_gaussian_lpdf(real y, real mu, real sigma, real hu) {
    if (y == 0) {
      return bernoulli_lpmf(1 | hu);
    } else {
      return bernoulli_lpmf(0 | hu) +
             normal_lpdf(y | mu, sigma);
    }
  }
"

posterior_predict_hurdle_gaussian <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  theta <- brms::get_dpar(prep, "hu", i = i)

  hu <- runif(prep$ndraws, 0, 1)
  ifelse(hu < theta, 0, rnorm(prep$ndraws, mu, sigma))
}

posterior_epred_hurdle_gaussian <- function(prep) {
  with(prep$dpars, mu * (1 - hu))
}

# prepare Stan code for use in brm()
stanvars <- stanvar(scode = stan_funs, block = "functions")
