extract.zip = function (file, unzip = getOption("unzip")) {
  if (!is.character(unzip) || length(unzip) != 1) 
    stop("'unzip' must be a single character string")
  if (!nzchar(unzip)) 
    unzip = "internal"
  path = dirname(file)
  topic = basename(file)
  if (file.exists(file.path(path, topic))) {
    if (unzip != "internal") {
      list = paste(unzip, "-l", file)
      cmd = paste(unzip, "-oq", file, " -d ", path)
      if (.Platform$OS.type == "windows")
        {
          rc = system(list, intern=T, invisible = TRUE)##does it work?
          rc = sapply(4:(length(rc)-2), function(i) strsplit(rc,"   ")[[i]][2])
          res = system(cmd, invisible = TRUE)
        } else {
          rc = system(list, intern=T)##works
          rc = paste("E-",sapply(4:(length(rc)-2), function(i) strsplit(rc," E-")[[i]][2]),sep="")
          res = system(paste(cmd, "> /dev/null"))
        }
    }
    else {
      rc = .Internal(int.unzip(file.path(path, topic),NULL, path)) ##works
      rc = attr(rc,"extracted")
    }
  } else stop(sprintf("%s does not exist",file.path(path, topic)))

  return(rc)
}
