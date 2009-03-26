extract.zip = function (file, extract_path = dirname(file)) {
  topic = basename(file)
  path = dirname(file)
  if (file.exists(file.path(path, topic))) {
    unzip(zipfile = file.path(path, topic), exdir = extract <- path)
    rc = unzip(zipfile = file.path(path, topic), list = TRUE, exdir = extract <- path)
    rc = as.character(rc[,1])
    if(inherits(rc,"try-error"))
      stop("Cannot extract the files the downloaded archive. Please install unzip on your machine")   
  } else stop(sprintf("%s does not exist",file.path(path, topic)))  
  return(rc)
}
