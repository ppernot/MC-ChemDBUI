function(input, output, session) {

  # Load Server files ####
  files <- c(
    "neutrals.R" ,
    "ions.R",
    "photo.R"
  )

  for (f in files)
    source(
      file.path("server_files", f),
      local = TRUE
    )
}
