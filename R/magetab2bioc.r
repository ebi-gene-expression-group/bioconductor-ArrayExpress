magetab2bioc = function(files, rawcol = NULL, save = TRUE)  
{
  if(!save) on.exit(cleanupAE(as.list(files)))
  
  rawfiles = files$rawfiles
  sdrf = files$sdrf
  idf = files$idf
  adf = files$adf
  path = files$path

  if (!is.null(rawcol)) {
    if (length(rawcol) > 1 && !is(rawcol, "list")) 
      stop("The argument 'rawcol' must be a list if multiple column name are given.")
    if (length(rawcol) == 1 && !is(rawcol, "character")) 
      stop("The argument 'rawcol' must be a character if one column name is given.")
    if (is(rawcol, "list") && !("R" %in% names(rawcol) && 
                                "G" %in% names(rawcol))) 
      stop("The names of the columns must contain R and G.")
  }

  notuse = grep(rawfiles, pattern = "info.txt$|idf.txt$|processed|sdrf.txt$|.log$|RData|class|log")
    
  if (length(notuse) != 0) 
    files = rawfiles[-notuse]
  if (length(notuse) == 0) 
    files = rawfiles

  basefilenames = basename(files)

  ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$"))

  if(inherits(ph, 'try-error') && length(grep(".cel",files, ignore.case = TRUE)) == 0)
  {
    warning(sprintf("No sdrf file available. The object may not be built."))
    adr = "Empty"
    raweset = try(nonAB(i = 1, files, path, ph, rawcol, adr, adf, idf))
  }
  if(inherits(ph, 'try-error') && length(grep(".cel",files, ignore.case = TRUE)) != 0)
    {
      warning(sprintf("No sdrf file available. The object may not be built."))
      adr = "Empty"
      raweset = try(AB(i = 1, files, path, ph, adr, adf, idf))

    }
  if(!inherits(ph, 'try-error'))
    {
      emptylines = which(sapply(seq_len(nrow(pData(ph))), function(i)
        all(pData(ph)[i,] == "",na.rm = TRUE)))
      if(length(emptylines) != 0)
        pData(ph) = pData(ph)[-emptylines,]
      
      if(!all(basefilenames %in% pData(ph)$Array.Data))
        warning("Some files from the zip archive do not have annotation in the sdrf file. The object may not be built.")

      adr = unique(pData(ph)$Array.Design.REF)
      adr = adr[adr != ""]
      if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",files, ignore.case = TRUE)) != 0)
        warning("Cannot find the array design reference in the sdrf file. The object may not be built.")
      if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",files, ignore.case = TRUE)) == 0)
        stop("Cannot find the array design reference in the sdrf file. The object cannot be built.")
    }
    
  ## Building the S4 class object
  if(length(grep(".cel",files, ignore.case = TRUE)) == length(files))
    {
      if(length(adr) == 1)
        raweset = try(AB(i = 1, files, path, ph, adr, adf, idf))

      if(length(adr) > 1)
        raweset = try(lapply(seq_len(length(adr)), function(i) try(AB(i, files, path, ph, adr, adf, idf))))
    }
   
  ## Non Affymetrix data
  if(length(grep(".cel",files, ignore.case = TRUE)) == 0)
    {
      if(length(adr) == 1)
        raweset = try(nonAB(i = 1, files, path, ph, rawcol, adr, adf, idf))
      if(length(adr) > 1)
        raweset = try(lapply(seq_len(length(adr)), function(i) try(nonAB(i, files, path, ph, rawcol, adr, adf, idf))))
    }
  return(raweset)
}#end of magetab2bioc
