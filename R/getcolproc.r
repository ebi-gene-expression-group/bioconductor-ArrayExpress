getcolproc = function(files){
  path = paste(files$path, "/", sep="")
  procfile = files$processedArchive[1]
  ph = try(read.AnnotatedDataFrame(basename(files$sdrf), path = path, row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
  derivedMatrixCol = getSDRFcolumn("DerivedArrayMatrix",varLabels(ph))
  derivedFileCol = getSDRFcolumn("DerivedArrayFile",varLabels(ph))
  
  if(length(derivedMatrixCol) != 0)
    skiplines = 1
  else if(length(derivedFileCol) != 0)
    skiplines = 0
  else
    warning("Unable to find Derived Data column in SDRF")
  
  coln = scan(file.path(procfile),what = "",nlines = 1, skip = skiplines, sep = "\t")
  return(unique(coln))
}

getcolraw = function(rawfiles){
  rawfile = rawfiles[1]
  coln = scan(file.path(rawfile),what = "",nlines = 1, sep = "\t")
  return(coln)
}