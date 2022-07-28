function(input, output, session) {

  # Load Server files ####
  files <- c(
    "neutralsEdit.R" ,
    "neutralsParse.R" ,
    "neutralsSample.R",
    "ionsEdit.R",
    "ionsSample.R",
    "photo.R"
  )

  for (f in files)
    source(
      file.path("server_files", f),
      local = TRUE
    )
}
