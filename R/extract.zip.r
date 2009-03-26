extract.zip = function (file, unzip = getOption("unzip"), extract_path = dirname(file)) {
  if (!is.character(unzip) || length(unzip) != 1) 
    stop("'unzip' must be a single character string")
  if (!nzchar(unzip)) 
    unzip = "internal"
  topic = basename(file)
  path = dirname(file)
  if (file.exists(file.path(path, topic))) {
    if (unzip != "internal") {
      list = paste(unzip, "-l", file)
      cmd = paste(unzip, "-oq", file, " -d ", extract_path)
      if (.Platform$OS.type == "windows")
	{
          rc = try(system(list, intern=T))
          if(inherits(rc,"try-error"))
            stop("Cannot extract the files the downloaded archive. Please install unzip on your machine")   
          rc = paste("E-",na.omit(sapply(1:length(rc), function(i) strsplit(rc," E-")[[i]][2])),sep="")
          if(length(grep("zip",rc)) != 0)
            rc = rc[-grep("zip",rc)]
          res = system(cmd, invisible = TRUE)
	} else {
          rc = system(list, intern=T)##works
          rc = paste("E-",na.omit(sapply(1:length(rc), function(i) strsplit(rc," E-")[[i]][2])),sep="")
          if(length(grep("zip",rc)) != 0)
            rc = rc[-grep("zip",rc)]
          res = system(paste(cmd, "> /dev/null"))
        }
    }
    else {
      rc = unzip(zipfile = file.path(path, topic), files = NULL, list = TRUE, exdir = path) ##works
      rc = as.character(rc[,1])
    }
  } else stop(sprintf("%s does not exist",file.path(path, topic)))
  
  return(rc)
}
