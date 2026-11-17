# =============================================================================
# user_paths.R (module)
#
# Per-user default file paths.  Each person keeps their own paths in
# `user_paths.R` at the project root, which is gitignored, so nobody's paths
# end up in the repository and nobody opens the app to someone else's.
# Start from `user_paths_template.R` (tracked).
#
# The app pre-fills every path box from that file at start-up; the
# diagnostic and analysis scripts read the same file.  Missing file = all
# boxes blank.  Paths can be written with forward slashes on any OS.
# =============================================================================

# Read <app_dir>/user_paths.R and return a named list of non-empty paths
# (backslashes normalised to "/").  Never errors: a broken file gives a
# warning and an empty list.
load_user_paths <- function(app_dir = ".") {
  f <- file.path(app_dir, "user_paths.R")
  if (!file.exists(f)) return(list())
  e <- new.env()
  ok <- tryCatch({
    source(f, local = e, encoding = "UTF-8")
    TRUE
  }, error = function(err) {
    warning("Could not read ", f, ": ", conditionMessage(err),
            "\n  (On Windows, write paths with forward slashes: C:/data/fwi.nc)")
    FALSE
  })
  if (!ok || !is.list(e$user_paths)) return(list())
  up   <- e$user_paths
  keep <- vapply(up, function(v) is.character(v) && length(v) == 1 &&
                   !is.na(v) && nzchar(trimws(v)), logical(1))
  lapply(up[keep], function(v) gsub("\\\\", "/", trimws(v)))
}

# One path from the list, or `default` when it is not set.
user_path <- function(paths, id, default = "") {
  v <- paths[[id]]
  if (is.null(v)) default else v
}

# Start-up report: entries the app does not know (typos) and paths that do
# not exist on this machine.  Returns the messages (character), invisibly.
check_user_paths <- function(paths, known_ids) {
  msgs <- character(0)
  unknown <- setdiff(names(paths), known_ids)
  if (length(unknown))
    msgs <- c(msgs, paste0("user_paths.R: unknown entries (typo?): ",
                           paste(unknown, collapse = ", ")))
  missing <- names(paths)[!vapply(paths, function(p) file.exists(p) || dir.exists(p),
                                   logical(1))]
  missing <- intersect(missing, known_ids)
  if (length(missing))
    msgs <- c(msgs, paste0("user_paths.R: not found on this machine: ",
                           paste(sprintf("%s = %s", missing, unlist(paths[missing])),
                                 collapse = "; ")))
  if (length(msgs)) message(paste(msgs, collapse = "\n"))
  invisible(msgs)
}
