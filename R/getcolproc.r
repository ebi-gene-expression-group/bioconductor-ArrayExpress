getcolproc = function(files)
  {
    path = files$path
    procfile = files$procfile[1]
    coln = scan(file.path(path,procfile),what = "",nlines = 1, skip = 1, sep = "\t")
    return(coln)
  }

getcolraw = function(path, rawfiles)
  {
    rawfile = rawfiles[1]
    coln = scan(file.path(path,rawfile),what = "",nlines = 1, sep = "\t")
    return(coln)
  }
