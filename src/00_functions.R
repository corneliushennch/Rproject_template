# Header start ==============================================================================
# 00_survival_analysis_functions.R
#
# Author: Hennch Cornelius (cornelius.hennch@charite.de)
#
# Description: Custom functions for the survival analysis of the PMBCL dataset
#
# Code written according to Hadley Wickhams "tidyverse style guide"
# Header end ==============================================================================

# 1. Data wrangling functions --------------------------------------------------
## 1.1 fct_count_all  ----------------------------------------------------------
# This function counts all n of all factors in a dataset
# Input: dataframe/tibble containing at least one factor variable
fct_count_all <- function(data) {
  factor_list <- data %>%
    #    filter(!is.na(.data[[event_variables[i]]])) %>% filtering to be done outside the function
    select_if(is.factor) %>% as.list()

  # counting
  factor_counts <- purrr::map(factor_list, fct_count)

  # binding back into a table
  factor_table <- Map(cbind, factor_counts, variable = names(factor_counts)) %>% bind_rows()

  # polishing
  factor_table <- factor_table %>%
    select(variable, everything()) %>%
    mutate(percentage = (n / length(factor_list[[1]])) * 100) %>%
    dplyr::rename(level = "f")

  return(factor_table)
}

## 1.2 as.numeric.factor   -----------------------------------------------------
# convert factor back to numeric
as.numeric.factor <- function(x) {
  as.numeric(levels(x))[x]
}

## 1.3 tbl_summary_list -----------------------------------------------------

# wrapper for table summary with styling

tbl_summary_list <- function(data, dependent, explanatory) {

  # remove missing data -> doesn't make it faster...
  data <- filter(data, !is.na(.data[[dependent]]))

  # TODO add total of examined vars

  var_n <- length(explanatory)

  # construct summary
  tbl_summary(
    data = data,
    by = dependent,
    include = all_of(explanatory),
    missing = "no"
  ) %>%
    # add_difference(pvalue_fun = ~style_pvalue(.x, digits = 3)) %>%
    # add_q(pvalue_fun = ~style_pvalue(.x, digits = 3), quiet = TRUE) %>%
    # sort_p(q = F) %>%
    modify_header(label ~ paste0("Variable (Total = ", var_n, ")")) %>%
    modify_spanning_header(c("stat_1", "stat_2") ~ dependent)
}

# 2. Regression model functions --------------------------------------------
## 2.1 coxfit -----------------------------------------------------------
# convenient wrapper for building cox models
# 1) data -> df containg outcome vars (time, event) and variables
# 2) time, event -> to construct "Surv" object
# 3) variabels -> character string of variables to be included in the model
# 4) ... -> additional args to coxph()

coxfit <- function(data, surv, variables, ...) {
  vars <- paste(variables, collapse = " + ")

  cox_formula_string <- paste0(surv, " ~ ", vars)
  cox_formula <- as.formula(cox_formula_string)

  data_name <- deparse(substitute(data))
  cox_call <- paste0(
    "coxph(formula = ", cox_formula_string,
    ", data = ", data_name, ")"
  ) %>% str2lang()

  cox_model <- coxph(cox_formula, data = data, ...)
  cox_model$call <- cox_call
  return(cox_model)
}


## 2.2 multi_coxph -------------------------------------------------------------------
# computes coxph models for list of variables supplied as a character vector
# 1) data -> df containing event vars and variables of interest as factors
# 2) formula -> character vector with Surv formula to which variables are being mapped to
# 3) variables -> one or more variables as character vector (or formula?)
# return: list of coxph models
multi_coxph <- function(data, formula, variables) {
  # purrr style -> might be slower
  # cox_models_molecular <- purrr::map(variables, ~coxph(as.formula(paste(formula, .x)), data = survival_dataset))

  # lapply style -> uses super fast multi-core computing mclapply
  formula_list <- map(variables, ~ as.formula(paste(formula, .x)))
  cox_model_list <- mclapply(formula_list, coxph, data = data)
  names(cox_model_list) <- variables
  return(cox_model_list)
}

## 2.3 tidy_coxph ---------------------------------------------------------------------
# takes a list of coxph model objects and returns a tidy dataframe with n and percentage added
tidy_coxph <- function(cox_models, factor_table) {
  # extracting tidy information of the cox model (HR, conf.int, p-value setc.)
  tidy_table <- lapply(cox_models, broom::tidy, exp = TRUE, conf.int = TRUE)
  names(tidy_table) <- names(cox_models)

  # adding a column for the variable name
  tidy_tabl_list <- Map(cbind, tidy_table, variable = as.list(names(cox_models)))

  # binding it back into a column
  tidy_table <- bind_rows(tidy_tabl_list) %>%
    # filter(!grepl(pattern = c("ipi|ldh"), term)) %>% apparently not necessary
    mutate(level = str_remove(term, variable)) %>%
    select(variable, level, everything())

  # return(tidy_table)

  # join summary  -> needs filtering
  tidy_table <- factor_table %>%
    filter(variable %in% tidy_table$variable) %>%
    left_join(tidy_table, by = c("variable", "level"))

  return(tidy_table)
}


## 2.4 filter_coxph ------------------------------------------------------------------
# takes the tidy list of coxph output (tidy_coxph) and returns a filtered list
# filtering options: prevalence, alpha
# if alpha == 0.05, confidence intervals are taken for filtering
# returns
filter_coxph <- function(tidy_coxph_table, prevalence = 25, alpha = 0.05) {

  # looking for potential effects
  if (alpha == 0.05) {
    potential_effects <- tidy_coxph_table %>%
      filter(conf.low > 1 & conf.high > 1 | conf.low < 1 & conf.high < 1) %>%
      filter(n >= prevalence) %>%
      # filtering for prevalence of respective lesion
      arrange(std.error)
  } else {
    potential_effects <- tidy_coxph_table %>%
      filter(p.value <= alpha) %>%
      filter(n >= prevalence) %>%
      # filtering for prevalence of respective lesion
      arrange(std.error)
  }

  # filtering for potential effects
  filtered_coxph <- tidy_coxph_table %>%
    filter(variable %in% potential_effects$variable)

  return(filtered_coxph)
}

## 2.5 extract_exp -----------------------------------------------------
# helper to extract explanatory variables of a model (e.g. after var selection)
extract_exp <- function(model) {
  explanatory <- paste(model$call[2]) %>%
    str_split(" ~ ") %>%
    flatten() %>%
    pluck(2) %>%
    str_split(" [+] ") %>%
    Reduce(c, .)
  return(explanatory)
}

## 2.6 compare cox models ------------------------------------------------------
# first model will be taken as reference
compare_models <- function(models, model_names = NULL, compact = TRUE) {

  # error message
  stopifnot("Models need to be stored in a list" = class(models) == "list")

  # create model names if not supplied
  if (is.null(model_names)) {
    if (is.null(names(models))) {
      model_names <- paste("Model", seq_along(models))
    } else {
      model_names <- names(models)
    }
  }
  # extract model formula
  model_formulas <- map(models, ~ pluck(.x, "call")) %>%
    map(as.character) %>%
    lapply("[", 2) %>%
    Reduce(c, .)

  # first part, basic info of the models
  model_info <- data.frame(name = model_names, formula = model_formulas) %>%
    separate(col = "formula", into = c("dependent", "explanatory"), sep = " ~ ") %>%
    as_tibble()

  # metrics of the model
  gmodel <- map(models, ~ broom::glance(.x)) %>%
    bind_rows() %>%
    janitor::remove_empty("cols") %>%
    dplyr::rename(loglik = "logLik")

  # binding it together
  model_df <- gmodel %>%
    select(AIC, BIC, concordance, everything()) %>%
    bind_cols(model_info, .)

  list_elements <- paste0("models[[", seq_along(models), "]]", collapse = ", ")

  # add degrees of freedom and compare with anova
  aic_info <- eval(parse(text = paste0("AIC(", list_elements, ")")))

  # choose model with lowest AIC for Anova reference
  anova_order <- aic_info %>%
    arrange(AIC) %>%
    rownames() %>%
    paste(collapse = ", ")

  # compare models
  anova_info <- eval(parse(text = paste0("anova(", anova_order, ")")))

  # final comparison df
  model_df <- left_join(model_df, anova_info, by = "loglik") %>%
    left_join(aic_info, by = "AIC") %>%
    dplyr::mutate(name = if_else(is.na(Df), paste(name, "(reference)"), name)) %>%
    dplyr::mutate(seq_p.value = if_else(round(`P(>|Chi|)`, 3) < 0.001, "<0.001",
      as.character(round(`P(>|Chi|)`, 3))
    )) %>%
    select(
      name, df, AIC, BIC, concordance, n, nevent,
      seq_p.value, loglik, everything()
    ) %>%
    relocate(dependent, explanatory, .after = last_col()) %>%
    arrange(AIC)

  # drop some info
  if (compact) {
    drop_cols <- c(
      "statistic.log", "statistic.sc", "p.value.sc", "r.squared.max",
      "std.error.concordance", "nobs", "Df"
    )
    model_df <- model_df %>%
      dplyr::select(!(all_of(drop_cols)))
  }
  return(model_df)
}

## 2.7 pairwise_anova ----------------------------------------------------------
# function for pairwise model comparison

pairwise_anova <- function(models) {
  stopifnot("Provide at least 3 models." = length(models) >= 3)

  # names -> will be row/column names
  name_pairs <- combn(names(models),
    m = 2,
    FUN = function(x) {
      paste(x, collapse = ", ")
    },
    simplify = FALSE
  )
  # pairs to test
  model_pairs <-
    combn(length(models),
      m = 2,
      FUN = function(x) {
        paste0("models[[", x, "]]", collapse = ", ")
      },
      simplify = TRUE
    )

  anova_list <- map(model_pairs, ~ eval(parse(text = paste0("anova(", .x, ")")))) %>%
    set_names(nm = model_pairs)

  anova_df <- anova_list %>%
    map2(., seq_along(.), ~ add_column(.x, pair = name_pairs[[.y]], .after = "Df")) %>%
    map(~ select(.x, pair, `P(>|Chi|)`)) %>%
    map(na.omit) %>%
    bind_rows() %>%
    as_tibble() %>%
    tidyr::separate(pair, into = c("Var1", "Var2"), sep = ", ") %>%
    dplyr::mutate(p.value = if_else(`P(>|Chi|)` < 0.001, 0.001, round(`P(>|Chi|)`, 3)))


  anova_matrix <- xtabs(p.value ~ Var1 + Var2, data = anova_df) %>%
    as.data.frame.matrix() %>%
    dplyr::mutate(across(everything(), ~ if_else(. == 0, true = NaN, .))) # %>%

  return(anova_matrix)
}

## 2.8 list wrapper for compare_performance ------------------------------------

compare_performance_list <- function(model_list, rank = FALSE) {

  # get model names
  model_names <- data.frame(
    original = names(model_list),
    Name = janitor::make_clean_names(names(model_list))
  )

  names(model_list) <- model_names$Name

  # pass list objects to the environment
  list2env(model_list, envir = environment())

  # collapse model names in one string
  model_str <- paste0(model_names$Name, collapse = ", ")

  # function call
  fn_call <- paste0("performance::compare_performance(", model_str, ", rank =", rank, ")")

  # call compare_performance
  compare_df <- eval(parse(text = fn_call))

  # retrieve original names
  compare_df$Name <- model_names$original[match(compare_df$Name, model_names$Name)]

  return(compare_df)
}

# the functions of the performance package perform the same tasks as the two
# functions above. I just wrote a convenient wrapper for
# performance::compare_performance which takes list as an input
## 2.9 logreg ------------------------------------------------------------------
# convenient wrapper for building logistic regression models with glm
# 1) data -> df containg outcome vars (time, event) and variables
# 2) dependent -> dependent variables
# 3) explanatory -> character string of variables to be included in the model
# 4) ... -> additional args to coxph()

logreg <- function(data, data_name = NULL, dependent, explanatory, ...) {
  vars <- paste(explanatory, collapse = " + ")

  formula_string <- paste0(dependent, " ~ ", vars)
  formula <- as.formula(formula_string)

  # name of the data needs to be supplied for vectorized handling
  if (is.null(data_name)) {
    data_name <- deparse(substitute(data))
  }

  call <- paste0(
    "glm(formula = ", formula_string, ", family = binomial",
    ", data = ", data_name, ")"
  ) %>% str2lang()

  model <- glm(formula, data = data, family = "binomial", ...)
  model$call <- call
  return(model)
}


# 3. plotting functions --------------------------------------------------------
## 3.1 custom ggplot forest plot ---------------------------------------------

forest_plot <- function(coxph_table, plot_title = "Hazard ratios",
                        plot_subtitle = "",
                        n_caption = FALSE) {
  if (!("variable_labels" %in% colnames(coxph_table))) {
    coxph_table$variable_labels <- coxph_table$variable
  }

  plot_data <- coxph_table %>%
    # add reference label variable
    mutate(level = if_else(is.na(term) & !is.na(level),
      true = paste(level, "(reference)", sep = " "),
      false = level
    )) %>%
    unite(col = term, variable_labels, level, sep = ": ", remove = FALSE) %>%
    unite(col = term, term, n, sep = "\n n = ", remove = FALSE) %>%
    # add dummy data for the reference levels
    mutate(estimate = if_else(grepl("reference", level), true = 1, false = estimate)) %>%
    mutate(std.error = if_else(grepl("reference", level), true = 0, false = std.error)) %>%
    mutate(prognostic_factor = case_when(
      conf.low > 1 & conf.high > 1 ~ "harmful",
      conf.low < 1 & conf.high < 1 ~ "favorable",
      conf.low < 1 & conf.high > 1 ~ "ns.",
      is.na(statistic) ~ "ns."
    )) %>%
    mutate(prognostic_factor = factor(prognostic_factor, levels = c("harmful", "favorable", "ns."))) %>%
    mutate(across(contains("conf"), ~ if_else(grepl("reference", level), true = 1, false = .))) %>%
    mutate(term = fct_inorder(term)) %>%
    filter(!conf.high == Inf)

  # calculating n samples for plot caption
  if (n_caption) {
    n_samples <- sum(coxph_table$n) / length(unique(coxph_table$variable))
    plot_caption <- paste("n = ", n_samples)
  } else {
    plot_caption <- ""
  }


  # creating the ggplot
  plot <- plot_data %>%
    ggplot(aes(y = fct_rev(term), x = estimate, xmin = conf.low, xmax = conf.high, color = prognostic_factor)) +
    # this adds the effect sizes to the plot
    geom_point(aes(size = n), shape = "diamond") + # aes(shape = variable_type) , color = p.value

    # adds the CIs
    geom_errorbarh(height = .2) +

    # sets the scales
    # scale_y_discrete(labels = toupper(rev(potential_effect_list[[i]]$variable)))+
    scale_x_log10(
      name = "Hazard Ratio (log10)", n.breaks = 6,
      labels = function(n) {
        format(n, scientific = FALSE, drop0trailing = TRUE)
      }
    ) +
    scale_color_manual(values = c(harmful = "red", favorable = "darkgreen", ns. = "black")) +

    # adding a vertical line at the effect = 0 mark
    geom_vline(xintercept = 1, color = "black", linetype = "dashed", alpha = .5) +

    # faceting based on my subgroups
    # facet_grid(Comparison~., scales= "free", space="free")+

    # thematic stuff
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = plot_caption,
      y = "Variables",
      shape = "variable type"
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 12, color = "black"),
      legend.position = "right"
    ) +
    theme(panel.spacing = unit(1, "lines"))

  return(plot)
}

## 3.2 customized hr_plot from the finalfit package -------------------
# TODO color direction of prognostic factor
# add metric
hr_plot_custom <- function(.data, dependent, explanatory, factorlist = NULL,
                           coxfit = NULL, remove_ref = FALSE, breaks = NULL,
                           n_breaks = 5, metrics = TRUE, print = FALSE,
                           add_n = TRUE, linebreak_n = FALSE,
                           prog_color = TRUE, return_df = FALSE,
                           column_space = c(-0.3, 0, 0.5),
                           dependent_label = "Survival", prefix = "",
                           suffix = ": HR (95% CI, p-value)", table_text_size = 4,
                           title_text_size = 16, plot_opts = NULL, table_opts = NULL,
                           ...) {
  # required packages for the function
  requireNamespace("ggplot2")
  requireNamespace("stringr")
  requireNamespace("magrittr")
  requireNamespace("dplyr")

  # chcek if correct factorlist is supplied
  if (!is.null(factorlist)) {
    if (is.null(factorlist$fit_id)) {
      stop("summary_factorlist function must include fit_id=TRUE")
    }
  }
  if (is.null(factorlist)) {
    factorlist <- summary_factorlist(.data, dependent, explanatory,
      fit_id = TRUE
    )
  }
  if (remove_ref) {
    factorlist <- factorlist %>%
      dplyr::mutate(label = ifelse(label ==
        "", NA, label)) %>%
      tidyr::fill(label) %>%
      dplyr::group_by(label) %>%
      dplyr::filter(dplyr::row_number() != 1 | dplyr::n() >
        2 | levels %in% c("Mean (SD)", "Median (IQR)")) %>%
      rm_duplicate_labels()
  }
  if (is.null(breaks)) {
    breaks <- scales::breaks_pretty(n_breaks)
  }
  factorlist$Total <- as.numeric(stringr::str_extract(
    as.character(factorlist$all),
    "^[:digit:]*"
  ))
  factorlist$Total[factorlist$levels == "Mean (SD)" | factorlist$levels ==
    "Median (IQR)"] <- dim(.data)[1]
  drop <- grepl("Mean \\(SD\\)|Median \\(IQR\\)", factorlist$levels)
  factorlist$levels[drop] <- "-"
  factorlist$all <- NULL

  # check if coxfit is supplied, create if not
  if (is.null(coxfit)) {
    coxfit <- coxphmulti(.data, dependent, explanatory)
  }
  coxfit_df_c <- fit2df(coxfit,
    condense = TRUE, estimate_suffix = " (multivariable)",
    estimate_name = "HR", exp = TRUE, ...
  )
  coxfit_df <- fit2df(coxfit,
    condense = FALSE, estimate_name = "HR",
    exp = TRUE, ...
  )
  df.out <- finalfit_merge(factorlist, coxfit_df_c)
  df.out <- finalfit_merge(df.out, coxfit_df, ref_symbol = "1.0")

  # add shortened metrics
  if (metrics) {
    gmodel <- broom::glance(coxfit)
    AIC <- paste(" \nAIC: ", round(gmodel$AIC, 2))

    metrics <- ff_metrics(coxfit) %>%
      as.character() %>%
      str_remove("\\).*") %>%
      paste0(")") %>%
      str_replace(", Co", replacement = paste0(AIC, ", Co"))
  } else {
    metrics <- ""
  }

  # specify labels
  if (any(is.na(df.out$label))) {
    remove_rows <- which(is.na(df.out$label))
    df.out <- df.out[-remove_rows, ]
  }
  else {
    df.out
  }
  df.out$levels <- as.character(df.out$levels)
  df.out$fit_id <- factor(df.out$fit_id, levels = df.out$fit_id[order(-df.out$index)])
  # add number of observations
  if (add_n) {
    if (linebreak_n) {
      df.out$levels <- paste0(df.out$levels, "\n(n = ", df.out$Total, ")")
    }
    else {
      df.out$levels <- paste0(df.out$levels, " (n = ", df.out$Total, ")")
    }
  }

  # add color according to direction of the prognostic factor
  if (prog_color) {
    df.out <- df.out %>%
      mutate(prognostic_factor = case_when(
        L95 > 1 & U95 > 1 ~ "harmful",
        L95 < 1 & U95 < 1 ~ "favorable",
        L95 < 1 & U95 > 1 ~ "ns.",
        is.na(p) ~ "ns."
      ) %>%
        factor(levels = c("harmful", "favorable", "ns.")))
  } else {
    df.out$prognostic_factor <- "ns."
  }
  # create forest plot
  g1 <- ggplot(df.out, aes(
    x = as.numeric(HR), y = fit_id,
    xmin = as.numeric(L95),
    xmax = as.numeric(U95),
    color = prognostic_factor
  )) +
    geom_point(aes(size = Total), shape = "diamond") + # , color = "black"
    geom_errorbarh(height = 0.2) +
    # annotate("text", x = HR_pos, y = df.out$fit_id,
    #          label = df.out[, 6], hjust = 1, size = table_text_size) +
    geom_vline(xintercept = 1, linetype = "longdash", colour = "black") +
    scale_x_continuous(
      trans = "log10", breaks = breaks,
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_color_manual(values = c(harmful = "red", favorable = "darkgreen", ns. = "black")) +
    xlab("Hazard ratio (95% CI, log scale)") +
    theme_classic(11) +
    theme(
      axis.title.x = element_text(), axis.title.y = element_blank(),
      axis.text.y = element_blank(), axis.line.y = element_blank(),
      axis.ticks.y = element_blank(), legend.position = "none",
      plot.margin = unit(c(5.5, 5.5, 5.5, -11), "points")
    )

  # create table with metrics
  t1 <- ggplot(df.out, aes(x = as.numeric(HR), y = fit_id)) +
    annotate("text",
      x = column_space[1], y = df.out$fit_id,
      label = df.out[, 2], hjust = 0, size = table_text_size,
      fontface = 2
    ) +
    annotate("text",
      x = column_space[2], y = df.out$fit_id,
      label = df.out[, 3], hjust = 0, size = table_text_size
    ) +
    annotate("text",
      x = column_space[3], y = df.out$fit_id,
      label = df.out[, 6], hjust = 1, size = table_text_size
    ) +
    theme_classic(12) +
    theme(
      axis.title.x = element_text(colour = "white"),
      axis.text.x = element_text(colour = "white"), axis.title.y = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      line = element_blank(), plot.margin = unit(c(5.5, -5.5, 5.5, -5.5), "points")
    )

  # t2 <- ggplot(df.out, aes(x = as.numeric(HR), y = fit_id)) +
  #   annotate("text", x = 0, y = df.out$fit_id,
  #            label = df.out[, 6], hjust = 1, size = table_text_size) +
  #   theme_classic(12) + theme(axis.title.x = element_text(colour = "white"),
  #                             axis.text.x = element_text(colour = "white"), axis.title.y = element_blank(),
  #                             axis.text.y = element_blank(), axis.ticks.y = element_blank(),
  #                             line = element_blank())
  # add options, if supplied
  g1 <- g1 + plot_opts
  t1 <- t1 + table_opts
  title <- plot_title(
    .data = .data, dependent = dependent,
    dependent_label = dependent_label, prefix = prefix,
    suffix = suffix
  )

  # save plot objects in list
  hr_plot_list <- lst(g1, t1, title, metrics)
  class(hr_plot_list) <- "ff_plot"

  # print or return plot objects
  if (print) {
    print(hr_plot_list, title_text_size = title_text_size)
  } else {
    if (return_df) {
      return(df.out)
    } else {
      return(hr_plot_list)
    }
  }
}


## 3.3 customized or_plot from the finalfit package ---------------------------
or_plot_custom <- function(.data, dependent, explanatory,
                           random_effect = NULL, factorlist = NULL, glmfit = NULL,
                           confint_type = NULL, remove_ref = FALSE, breaks = NULL,
                           n_breaks = 5, metrics = TRUE, print = FALSE, add_n = TRUE,
                           linebreak_n = FALSE, prog_color = TRUE, return_df = FALSE,
                           column_space = c(-0.3, 0, 0.5), dependent_label = NULL,
                           prefix = "", suffix = ": OR (95% CI, p-value)",
                           table_text_size = 4,
                           title_text_size = 16,
                           plot_opts = NULL, table_opts = NULL, ...) {
  requireNamespace("ggplot2")

  # Generate or format factorlist object
  if (!is.null(factorlist)) {
    if (is.null(factorlist$Total)) stop("summary_factorlist function must include total_col=TRUE")
    if (is.null(factorlist$fit_id)) stop("summary_factorlist function must include fit_id=TRUE")
  }

  if (is.null(factorlist)) {
    factorlist <- summary_factorlist(.data, dependent, explanatory, total_col = TRUE, fit_id = TRUE)
  }

  if (remove_ref) {
    factorlist <- factorlist %>%
      dplyr::mutate(label = ifelse(label == "", NA, label)) %>%
      tidyr::fill(label) %>%
      dplyr::group_by(label) %>%
      dplyr::filter(dplyr::row_number() != 1 |
        dplyr::n() > 2 |
        levels %in% c("Mean (SD)", "Median (IQR)")) %>%
      rm_duplicate_labels()
  }

  if (is.null(breaks)) {
    breaks <- scales::pretty_breaks(n_breaks)
  }

  # Confidence intervals, default to "profile" for glm and "Wald" for glmer
  if (is.null(confint_type) && is.null(random_effect)) {
    confint_type <- "profile"
  } else if (is.null(confint_type) && (!is.null(random_effect) | class(glmfit) == "glmerMod")) {
    confint_type <- "default"
  }

  # Generate or format glm
  if (is.null(glmfit) && is.null(random_effect)) {
    glmfit <- glmmulti(.data, dependent, explanatory)
    glmfit_df_c <- fit2df(glmfit,
      condense = TRUE, estimate_suffix = " (multivariable)",
      confint_type = confint_type, ...
    )
  } else if (is.null(glmfit) && !is.null(random_effect)) {
    glmfit <- glmmixed(.data, dependent, explanatory, random_effect)
    glmfit_df_c <- fit2df(glmfit,
      condense = TRUE, estimate_suffix = " (multilevel)",
      confint_type = confint_type, ...
    )
  }
  if (!is.null(glmfit) && is.null(random_effect)) {
    glmfit_df_c <- fit2df(glmfit,
      condense = TRUE, estimate_suffix = " (multivariable)",
      confint_type = confint_type, estimate_name = "OR", exp = TRUE, ...
    )
  } else if (!is.null(glmfit) && !is.null(random_effect)) {
    glmfit_df_c <- fit2df(glmfit,
      condense = TRUE, estimate_suffix = " (multilevel)",
      confint_type = confint_type, estimate_name = "OR", exp = TRUE, ...
    )
  }

  # add shortened metrics
  if (metrics) {
    metrics <- ff_metrics(glmfit) %>%
      as.character() %>%
      str_replace(", C-", replacement = paste0(", \nC-"))
  } else {
    metrics <- ""
  }

  glmfit_df <- fit2df(glmfit, condense = FALSE, confint_type = confint_type, estimate_name = "OR", exp = TRUE, ...)

  # Merge
  df.out <- finalfit_merge(factorlist, glmfit_df_c)
  df.out <- finalfit_merge(df.out, glmfit_df, ref_symbol = "1.0")

  # Remove proportions from total column and make continuous explanatory reflect dataset
  df.out$Total <- stringr::str_remove(df.out$Total, " \\(.*\\)") %>%
    as.numeric()
  df.out$Total[which(df.out$levels %in% c("Mean (SD)", "Median (IQR)"))] <- dim(.data)[1]

  # For continuous variables, remove level label
  df.out$levels[which(df.out$levels %in% c("Mean (SD)", "Median (IQR)"))] <- "-"

  # Remove unwanted lines, where there are more variables in model than wish to display.
  # These not named in factorlist, creating this problem. Interactions don't show on plot.
  if (any(is.na(df.out$label))
  ) {
    remove_rows <- which(is.na(df.out$label)) # This row doesn't work when is.na == FALSE, hence if()
    df.out <- df.out[-remove_rows, ]
  } else {
    df.out
  }

  # Fix order
  df.out$levels <- as.character(df.out$levels)
  df.out$fit_id <- factor(df.out$fit_id, levels = df.out$fit_id[order(-df.out$index)])

  # add number of observations
  if (add_n) {
    if (linebreak_n) {
      df.out$levels <- paste0(df.out$levels, "\n(n = ", df.out$Total, ")")
    }
    else {
      df.out$levels <- paste0(df.out$levels, " (n = ", df.out$Total, ")")
    }
  }

  # add color according to direction of the prognostic factor
  if (prog_color) {
    df.out <- df.out %>%
      mutate(prognostic_factor = case_when(
        L95 > 1 & U95 > 1 ~ "harmful",
        L95 < 1 & U95 < 1 ~ "favorable",
        L95 < 1 & U95 > 1 ~ "ns.",
        is.na(p) ~ "ns."
      ) %>%
        factor(levels = c("harmful", "favorable", "ns.")))
  } else {
    df.out$prognostic_factor <- "ns."
  }

  # Plot
  g1 <- ggplot(df.out, aes(
    x = as.numeric(OR),
    xmin = as.numeric(L95),
    xmax = as.numeric(U95),
    y = fit_id, color = prognostic_factor
  )) +
    geom_point(aes(size = Total), shape = "diamond") +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "longdash", colour = "black") +
    scale_x_continuous(
      trans = "log10", breaks = breaks,
      labels = scales::label_number(accuracy = 0.01)
    ) +
    scale_color_manual(values = c(harmful = "red", favorable = "darkgreen", ns. = "black")) +
    xlab("Odds ratio (95% CI, log scale)") +
    theme_classic(14) +
    theme(
      axis.title.x = element_text(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(5.5, 5.5, 5.5, -11), "points")
    )

  t1 <- ggplot(df.out, aes(x = as.numeric(OR), y = fit_id)) +
    annotate("text", x = column_space[1], y = df.out$fit_id, label = df.out[, 2], hjust = 0, size = table_text_size) +
    annotate("text", x = column_space[2], y = df.out$fit_id, label = df.out[, 3], hjust = 0, size = table_text_size) +
    annotate("text", x = column_space[3], y = df.out$fit_id, label = df.out[, 8], hjust = 1, size = table_text_size) +
    theme_classic(14) +
    theme(
      axis.title.x = element_text(colour = "white"),
      axis.text.x = element_text(colour = "white"),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      line = element_blank()
    )

  # Add optional arguments
  g1 <- g1 + plot_opts
  t1 <- t1 + table_opts

  # Add dependent name label
  title <- plot_title(.data, dependent,
    dependent_label = dependent_label,
    prefix = prefix, suffix = suffix
  )

  # save plot objects in list
  plot_list <- lst(g1, t1, title, metrics)
  class(plot_list) <- "ff_plot"

  # print or return plot objects
  if (print) {
    print(plot_list, title_text_size = title_text_size)
  } else {
    if (return_df) {
      return(df.out)
    } else {
      return(plot_list)
    }
  }
}


## 3.4 print.ff_plot -----------------------------------------------------
# common printing method for the ff plots (or and hr custom),
# if components of the plot were generated
print.ff_plot <- function(hr_plot_list, title_text_size = 16,
                          plot_widths = c(3.4, 1.6)) {

  # export list objects to function envir.
  list2env(hr_plot_list, envir = environment())

  # arrange plot
  gridExtra::grid.arrange(t1, g1,
    ncol = 2, widths = plot_widths,
    top = grid::textGrob(title,
      x = 0.02, y = 0.2,
      gp = grid::gpar(fontsize = title_text_size),
      just = "left"
    ),
    bottom = grid::textGrob(metrics,
      x = 0.02, y = 0.5,
      gp = grid::gpar(fontsize = title_text_size * 0.6),
      just = "left"
    )
  )
}


## 3.3 custom ggsurvplot wrapper ----------------------------------------------
# adds pairwise p-values (taken from the movics package)
km_plot_custom <-
  function(fit,
           data = NULL,
           pval = FALSE,
           pval.size = 4,
           pval.table = FALSE,
           pval.table.just = "left", # specifies alignment of the pval.table
           p.adjust.method = "BH",
           risk.table = TRUE,
           conf.int = FALSE,
           xlim = c(0, 240),
           break.time.by = 24,
           risk.table.y.text.col = TRUE,
           risk.table.y.text = FALSE,
           ggtheme =
             theme_survminer(
               font.main = c(18, "plain", "black"),
               font.submain = c(12, "plain", "black"),
               font.caption = c(11, "plain", "black")
             ),
           ...) {
    km_plot <- ggsurvplot(
      fit = fit,
      data = data,
      pval = pval,
      pval.size = pval.size,
      risk.table = risk.table,
      conf.int = conf.int,
      xlim = xlim,
      break.time.by = break.time.by,
      risk.table.y.text.col = risk.table.y.text.col,
      risk.table.y.text = risk.table.y.text,
      ggtheme = ggtheme,
      ...
    )

    # extract survfit call for pairwise logrank test
    call <- as.character(fit$call)
    surv_formula <- as.formula(call[2])

    if (pval.table) {
      # calculate pair-wise survival comparison
      ps <- pairwise_survdiff(
        formula = surv_formula,
        data = data,
        p.adjust.method = p.adjust.method
      )

      # add pair-wise comparison table
      # options(stringsAsFactors = FALSE)
      addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
        round(ps$p.value, 3)
      )))
      addTab[is.na(addTab)] <- "-"
      # options(stringsAsFactors = TRUE)
      if (pval.table.just == "left") {
        df <- tibble(x = 0, y = 0, tb = list(addTab))
      }
      if (pval.table.just == "right") {
        df <- tibble(x = xlim[2], y = 0, tb = list(addTab))
      }
      # add table to plot
      km_plot$plot <- km_plot$plot +
        ggpmisc::geom_table(
          data = df,
          aes(x = x, y = y, label = tb),
          hjust = pval.table.just,
          vjust = "bottom",
          table.rownames = TRUE
        )

      return(km_plot)
    } else {
      #   #calculate overall p.value
      #   p.val <- survdiff(surv_formula, data = data)
      #
      #   # add nominal pvalue for log-rank test
      #   p.lab <- paste0("Overall P",
      #                   ifelse(p.val < 0.001, " < 0.001",
      #                          paste0(" = ",round(p.val, 3))))
      #
      #   p$plot <- p$plot + annotate("text",
      #                               x = 0, y = 0.55,
      #                               hjust = 0,
      #                               fontface = 4,
      #                               label = p.lab)
      return(km_plot)
    }
  }

### 2.3.1 add p.value table ----------------------------------------------------
# geom like wrapper to supply pairwise logrank p-values to KM plot
# needs formula of the fitted survival curve and the same data
geom_ptable <- function(surv_formula, data,
                        p.adjust.method = "BH") {
  ps <- pairwise_survdiff(
    formula = surv_formula,
    data = data,
    p.adjust.method = p.adjust.method
  )

  # add pair-wise comparison table
  # options(stringsAsFactors = FALSE)
  addTab <- as.data.frame(as.matrix(ifelse(round(ps$p.value, 3) < 0.001, "<0.001",
    round(ps$p.value, 3)
  )))
  addTab[is.na(addTab)] <- "-"
  # options(stringsAsFactors = TRUE)
  if (pval.table.just == "left") {
    df <- tibble(x = 0, y = 0, tb = list(addTab))
  }
  if (pval.table.just == "right") {
    df <- tibble(x = xlim[2], y = 0, tb = list(addTab))
  }
  # add table to plot
  p.table <- ggpmisc::geom_table(
    data = df,
    aes(x = x, y = y, label = tb),
    hjust = pval.table.just,
    vjust = "bottom",
    table.rownames = TRUE
  )
  return(p.table)
}

# # Applied statistics functions
# ## RMSE
# calc_loocv_rmse = function(model) {
#   sqrt(mean((resid(model) / (1 - hatvalues(model))) ^ 2))
# }
