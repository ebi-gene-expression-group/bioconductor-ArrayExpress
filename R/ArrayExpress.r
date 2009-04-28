ArrayExpress = function(input, path = getwd(), save = FALSE, rawcol = NULL)
  {
    if(!save) on.exit(cleanupAE(extract)) 

    ## Extracting the data if the checking was fine
    extract = getAE(input, path, type = "raw")
  
    raweset = magetab2bioc(files = extract,
      rawcol = rawcol,
      save = save)
    
    if(!inherits(raweset, 'try-error'))
      message(paste("\n The object containing experiment ", input," has been built.\n"))
    
    return(raweset)
  }#end of ArrayExpress
