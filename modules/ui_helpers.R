# =============================================================================
# ui_helpers.R
# Minor shared UI utility functions
# =============================================================================

# Status badge helpers
status_ok  <- function(msg) tags$span(class = "text-success small", icon("check"), " ", msg)
status_err <- function(msg) tags$span(class = "text-danger  small", icon("xmark"), " ", msg)
status_run <- function(msg) tags$span(class = "text-warning small",
                                      icon("spinner", class = "fa-spin"), " ", msg)
