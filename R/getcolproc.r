getcolproc = function(files){
    path = files$path
    procfile = files$processedFiles[1]
    coln = scan(file.path(path,procfile),what = "",nlines = 1, skip = 1, sep = "\t")
    return(unique(coln))
}

getcolraw = function(path, rawfiles){
    rawfile = rawfiles[1]
    coln = scan(file.path(path,rawfile),what = "",nlines = 1, sep = "\t")
    return(coln)
}
