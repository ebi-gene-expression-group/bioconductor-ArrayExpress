extract.zip = function (file, extractpath = dirname(file)[1]) {
  topic = basename(file)
  path = dirname(file)
  if (file.exists(file.path(path, topic))) {
    lapply(1:length(file), function(i) unzip(zipfile = file.path(path[i], topic[i]), exdir = extractpath))
    rc = lapply(1:length(file), function(i) unzip(zipfile = file.path(path[i], topic[i]), list = TRUE, exdir = extractpath))
    if(length(rc) > 1)
      rc = do.call(c, sapply(1:length(file), function(i) as.character(rc[[i]][,1]))) else rc = as.character(rc[[1]][,1])    
    if(inherits(rc,"try-error"))
      stop("Cannot extract the files the downloaded archive. Please install unzip on your machine")   
  } else stop(sprintf("%s does not exist",file.path(path, topic)))  
  return(rc)
}
