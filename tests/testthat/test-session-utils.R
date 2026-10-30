test_that("cleanup_session_dir refuses unsafe paths and deletes a valid session folder", {
  root <- tempfile(pattern = "cleanup-root-")
  dir.create(root)
  outside <- tempfile(pattern = "cellivista_session_", tmpdir = tempfile(pattern = "other-"))
  dir.create(outside, recursive = TRUE)
  wrong_name <- file.path(root, "not_a_session")
  dir.create(wrong_name)
  valid <- tempfile(pattern = "cellivista_session_", tmpdir = root)
  dir.create(file.path(valid, "qc"), recursive = TRUE)
  writeLines("x", file.path(valid, "qc", "plot.txt"))

  expect_false(cleanup_session_dir(outside, root = root))
  expect_true(dir.exists(outside))
  expect_false(cleanup_session_dir(root, root = root))
  expect_true(dir.exists(root))
  expect_false(cleanup_session_dir(wrong_name, root = root))
  expect_true(dir.exists(wrong_name))
  expect_false(cleanup_session_dir(NULL, root = root))
  expect_false(cleanup_session_dir(file.path(root, "cellivista_session_abc123"), root = root))

  expect_true(cleanup_session_dir(valid, root = root))
  expect_false(dir.exists(valid))
  expect_true(dir.exists(root))

  skip_on_os("windows")
  link <- file.path(root, "cellivista_session_feed01")
  file.symlink(outside, link)
  expect_false(cleanup_session_dir(link, root = root))
  expect_true(dir.exists(outside))
})

test_that("zip_folder_relative keeps the example dataset folder structure", {
  skip_if(!nzchar(Sys.which(Sys.getenv("R_ZIPCMD", "zip"))), "zip command not available")
  parent <- tempfile(pattern = "zip-staging-")
  rel_files <- c(
    file.path("Raw Data", "GSE132044_pbmc_hg38_count_matrix.mtx.gz"),
    file.path("Raw Data", "GSE132044_pbmc_hg38_cell.tsv.gz"),
    file.path("Raw Data", "GSE132044_pbmc_hg38_gene.tsv.gz"),
    "rds_convert.R",
    "GSE132044_pbmc_hg38.rds"
  )
  for (rel in rel_files) {
    dest <- file.path(parent, "GSE132044", rel)
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    writeLines(rel, dest)
  }
  wd_before <- getwd()
  methods_to_test <- if (.Platform$OS.type == "windows") FALSE else c(TRUE, FALSE)
  for (use_shell in methods_to_test) {
    zipfile <- tempfile(fileext = ".zip")
    zip_folder_relative(parent, "GSE132044", zipfile, use_shell = use_shell)
    entries <- utils::unzip(zipfile, list = TRUE)$Name
    expect_setequal(entries[!grepl("/$", entries)], file.path("GSE132044", rel_files))
    expect_identical(getwd(), wd_before)
  }
  expect_error(zip_folder_relative(parent, "missing", tempfile(fileext = ".zip")), "Folder to zip not found")
})
