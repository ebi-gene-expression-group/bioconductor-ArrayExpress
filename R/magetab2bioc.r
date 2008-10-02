magetab2bioc = function(rawfiles, sdrf, idf, path = ".", columns = NULL, save = TRUE) {
  if(!is.null(columns))
    {
      if(length(columns) > 1 && !is(columns,"list"))
        stop("The argument 'columns' must be a list if multiple column name are given.")
      if(length(columns) == 1 && !is(columns,"character"))
        stop("The argument 'columns' must be a character if one column name is given.")
      if(is(columns,"list") && !("R" %in% names(columns) && "G" %in% names(columns)))
        stop("The names of the columns must contain R and G.")
    }

  if(!save) on.exit(try(file.remove(file.path(path, rawfiles))))

  notuse = grep(rawfiles, pattern = "info.txt$|idf.txt$|processed|sdrf.txt$|.log$|RData|class|log")
  if(length(notuse) != 0)
    files = rawfiles[-notuse]
  if(length(notuse) == 0)
    files = rawfiles
  
  ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))

  if(inherits(ph, 'try-error') && length(grep(".cel",files)) == 0)
    stop(sprintf("No sdrf file available. The object cannot be built."))
  if(inherits(ph, 'try-error') && length(grep(".cel",files)) != 0)
    {
      warning(sprintf("No sdrf file available. The object may not be built."))
      adr = "Empty"
      raweset = try(AB(i=1, files, path, ph, adr))

    }
  if(!inherits(ph, 'try-error'))
    {
      emptylines = which(sapply(seq_len(nrow(pData(ph))), function(i)
        all(pData(ph)[i,]=="",na.rm=T)))
      if(length(emptylines) != 0)
        pData(ph) = pData(ph)[-emptylines,]
      
      if(!all(files %in% pData(ph)$Array.Data))
        warning("Some files from the zip archive do not have annotation in the sdrf file. The object may not be built.")

      adr = unique(pData(ph)$Array.Design.REF)
      adr = adr[adr!=""]
      if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",files)) != 0)
        warning("Cannot find the array design reference in the sdrf file. The object may not be built.")
      if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",files)) == 0)
        stop("Cannot find the array design reference in the sdrf file. The object cannot be built.")
    }
    
  ## Building the S4 class object
  if(length(grep(".cel",files)) == length(files))
    {
      if(length(adr) == 1)
        raweset = try(AB(i=1, files, path, ph, adr))

      if(length(adr) > 1)
        raweset = try(lapply(seq_len(length(adr)), function(i) try(AB(i, files, path, ph, adr))))
    }
   
  ## Non Affymetrix data
  if(length(grep(".cel",files)) == 0)
    {
      if(length(adr) == 1)
        raweset = try(nonAB(i=1, files, path, ph, columns, adr))
      if(length(adr) > 1)
        raweset = try(lapply(seq_len(length(adr)), function(i) try(nonAB(i, files, path, ph, columns, adr))))
    }

 rawesetex = try(creating_experiment(idf = idf, eset = raweset, path = path))
  if(inherits(raweset, 'try-error'))
    warning(sprintf("Cannot create the experimentData."))
  if(!inherits(ph, 'try-error'))
    raweset = rawesetex
  
  return(raweset)
}#end of magetab2bioc
