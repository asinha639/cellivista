#' Delete a Cellivista Session Folder Safely
#'
#' Removes the per-session temporary folder created by the Shiny app when a
#' session ends. The folder is deleted only if all checks pass: it exists,
#' its normalized path sits directly inside `root` (normally `tempdir()`), it
#' is not `root` itself, and its name matches `pattern`. Symbolic links are
#' resolved first, so a link pointing outside `root` is refused. Any failed
#' check or error skips deletion without raising an error.
#'
#' @param path Path of the session folder.
#' @param root Directory the folder must sit directly inside. Default is
#'   `tempdir()`.
#' @param pattern Regular expression the folder name must match. Default
#'   matches names from `tempfile(pattern = "cellivista_session_")`.
#'
#' @return TRUE if the folder was deleted, FALSE otherwise (invisibly).
#' @export
#'
#' @examples
#' session_dir <- tempfile(pattern = "cellivista_session_")
#' dir.create(session_dir)
#' cleanup_session_dir(session_dir)
cleanup_session_dir <- function(path,
                                root = tempdir(),
                                pattern = "^cellivista_session_[0-9a-f]+$") {
  tryCatch({
    if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
      return(invisible(FALSE))
    }
    if (!dir.exists(path)) {
      return(invisible(FALSE))
    }
    target <- normalizePath(path, winslash = "/", mustWork = TRUE)
    base <- normalizePath(root, winslash = "/", mustWork = TRUE)
    if (identical(target, base)) {
      return(invisible(FALSE))
    }
    if (!identical(normalizePath(dirname(target), winslash = "/", mustWork = TRUE), base)) {
      return(invisible(FALSE))
    }
    if (!grepl(pattern, basename(target))) {
      return(invisible(FALSE))
    }
    unlink(target, recursive = TRUE)
    invisible(!dir.exists(target))
  }, error = function(e) invisible(FALSE))
}

#' Zip a Folder While Keeping Its Relative Paths
#'
#' Creates `zipfile` containing `folder` (a subfolder of `parent`) with paths
#' stored relative to `parent`, e.g. `GSE132044/Raw Data/...`. On Linux and
#' macOS the `zip` command runs in a child shell that changes into `parent`,
#' so the R process's working directory never changes. On Windows, where `sh`
#' is usually unavailable, the working directory is changed temporarily and
#' restored on exit.
#'
#' @param parent Directory that contains `folder`.
#' @param folder Name of the folder to zip, relative to `parent`.
#' @param zipfile Path of the zip file to create.
#' @param use_shell Logical; use the child-shell method. Default is TRUE on
#'   all platforms except Windows.
#'
#' @return The path of the zip file, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' zip_folder_relative(staging_dir, "GSE132044", "GSE132044_example_dataset.zip")
#' }
zip_folder_relative <- function(parent, folder, zipfile,
                                use_shell = .Platform$OS.type != "windows") {
  if (!dir.exists(file.path(parent, folder))) {
    stop("Folder to zip not found: ", file.path(parent, folder))
  }
  parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  zipfile <- file.path(normalizePath(dirname(zipfile), winslash = "/", mustWork = TRUE), basename(zipfile))
  files <- file.path(folder, list.files(file.path(parent, folder), recursive = TRUE, all.files = FALSE))

  if (isTRUE(use_shell)) {
    zip_cmd <- Sys.getenv("R_ZIPCMD", "zip")
    if (!nzchar(zip_cmd)) zip_cmd <- "zip"
    cmd <- paste(
      "cd", shQuote(parent), "&&",
      shQuote(zip_cmd), "-r9Xq", shQuote(zipfile), paste(shQuote(files), collapse = " ")
    )
    status <- system2("sh", c("-c", shQuote(cmd)))
  } else {
    old_wd <- setwd(parent)
    on.exit(setwd(old_wd), add = TRUE)
    status <- utils::zip(zipfile = zipfile, files = files)
  }

  if (!identical(as.integer(status), 0L) || !file.exists(zipfile)) {
    stop("Could not create the zip file (zip exit status ", status, "). Check that the 'zip' command is installed.")
  }
  invisible(zipfile)
}
