ArrayExpress = function(input, path = ".", save = TRUE, rawcol = NULL)
  {
    if(!save) on.exit({try(file.remove(file.path(path,extract$sdrf))); try(file.remove(file.path(path,extract$idf)))})
   
    ## Extracting the data if the checking was fine
    extract = getAE(input, path, save, type = "raw")
  
    raweset = magetab2bioc(rawfiles = extract$rawfiles,
      sdrf = extract$sdrf,
      idf = extract$idf,
      adf = extract$adf,
      path = path,
      rawcol = rawcol,
      save = save)
    
    if(!inherits(raweset, 'try-error'))
      message(paste("\n The object containing experiment ", input," has been built.\n"))
    
    return(raweset)
  }#end of ArrayExpress
