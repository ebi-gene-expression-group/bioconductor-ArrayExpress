ArrayExpress = function(input, path = ".", save = TRUE, rawcol = NULL)
  {
    if(!save) on.exit({try(file.remove(file.path(path,extract$drf))); try(file.remove(file.path(path,extract$idf))); try(file.remove(file.path(path,extract$adf)))})
   
    ## Extracting the data if the checking was fine
    extract = getAE(input, path, save, type = "raw")
  
    raweset = magetab2bioc(files = extract,
      rawcol = rawcol,
      save = save)
    
    if(!inherits(raweset, 'try-error'))
      message(paste("\n The object containing experiment ", input," has been built.\n"))
    
    return(raweset)
  }#end of ArrayExpress
