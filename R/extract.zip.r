extract.zip = function (file, unzip = getOption("unzip")) {
  if (!is.character(unzip) || length(unzip) != 1) 
    stop("'unzip' must be a single character string")
  if (!nzchar(unzip)) 
    unzip <- "internal"
  path <- dirname(file)
  topic <- basename(file)
  if (file.exists(file.path(path, topic))) {
    if (unzip != "internal") {
      cmd <- paste(unzip, "-oq", file, " -d ", path)
      res <- if (.Platform$OS.type == "windows") 
        system(cmd, invisible = TRUE)
      else system(paste(cmd, "> /dev/null"))
      if (!res) 
        file <- file.path(path, topic)
    }
    else {
      rc <- .Internal(int.unzip(file.path(path, topic), 
                                NULL, path))
      if (rc == 0) 
        file <- file.path(path, topic)
    }
  }
  path
}
