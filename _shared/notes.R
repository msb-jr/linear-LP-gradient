# Figure/table notes for outputs.
#
# Each script calls build_notes() with a list of legends, producing
# a plain text file with captions and footnotes for each output.

wrap_text <- function(text, width = 78) {
  if (is.null(text) || nchar(text) == 0) return("")
  words <- strsplit(text, " ")[[1]]
  lines <- character()
  current_line <- ""
  for (word in words) {
    test_line <- if (nchar(current_line) == 0) word else paste(current_line, word)
    if (nchar(test_line) <= width) {
      current_line <- test_line
    } else {
      if (nchar(current_line) > 0) lines <- c(lines, current_line)
      current_line <- word
    }
  }
  if (nchar(current_line) > 0) lines <- c(lines, current_line)
  lines
}

#' Build a notes file from figure/table legends
#'
#' @param legends List of lists, each with:
#'   - target: str (e.g., "Figure 1", "Table 1", or filename)
#'   - caption: str (the legend/caption text)
#'   - footnotes: character vector, optional (table footnotes)
#'   - abbreviations: str, optional (abbreviation definitions)
#' @param title Optional header title
#' @return Formatted notes text as a single string
build_notes <- function(legends, title = NULL) {
  lines <- character()

  if (!is.null(title) && nchar(title) > 0) {
    lines <- c(lines, toupper(title))
    lines <- c(lines, paste(rep("=", 78), collapse = ""))
    lines <- c(lines, paste0("Generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
    lines <- c(lines, "")
  }

  for (i in seq_along(legends)) {
    legend <- legends[[i]]
    target <- legend[["target"]] %||% paste0("Output ", i)
    caption <- legend[["caption"]] %||% ""
    footnotes <- legend[["footnotes"]] %||% character()
    abbreviations <- legend[["abbreviations"]] %||% ""

    # Target header
    lines <- c(lines, paste0("[", target, "]"))
    lines <- c(lines, "")

    # Caption (wrapped)
    if (nchar(caption) > 0) {
      wrapped <- wrap_text(caption, width = 78)
      lines <- c(lines, wrapped)
      lines <- c(lines, "")
    }

    # Footnotes
    if (length(footnotes) > 0) {
      for (fn in footnotes) {
        fn_wrapped <- wrap_text(fn, width = 74)
        lines <- c(lines, fn_wrapped[1])
        if (length(fn_wrapped) > 1) {
          for (cont in fn_wrapped[-1]) {
            lines <- c(lines, paste0("    ", cont))
          }
        }
      }
      lines <- c(lines, "")
    }

    # Abbreviations
    if (nchar(abbreviations) > 0) {
      abbr_line <- paste0("Abbreviations: ", abbreviations)
      abbr_wrapped <- wrap_text(abbr_line, width = 78)
      lines <- c(lines, abbr_wrapped)
      lines <- c(lines, "")
    }

    # Separator
    lines <- c(lines, paste(rep("-", 78), collapse = ""))
    lines <- c(lines, "")
  }

  paste(trimws(lines, which = "right"), collapse = "\n")
}

# Null-coalescing operator
`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
}

# Helper formatters for dynamic values
format_p <- function(value, digits = 3) {
  if (is.null(value) || !is.finite(value)) return("p = NA")
  if (value < 0.001) return("p < 0.001")
  sprintf(paste0("p = %.", digits, "f"), value)
}

format_ci <- function(estimate, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f (95%% CI: %.", digits, "f–%.", digits, "f)"),
          estimate, lower, upper)
}

format_median_iqr <- function(median, q1, q3, digits = 1) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
          median, q1, q3)
}
