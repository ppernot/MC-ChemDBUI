function(input, output, session) {

  # Load Server files ####
  files <- c(
    "neutralsLoad.R" ,
    "neutralsEdit.R" ,
    "neutralsSample.R",
    "ionsLoad.R",
    "ionsEdit.R",
    "ionsSample.R",
    "photoLoad.R",
    "photoEdit.R",
    "photoSample.R"
  )
  
  for (f in files)
    source( file.path("server_files", f), local = TRUE )
}
