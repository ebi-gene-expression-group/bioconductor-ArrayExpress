ArrayExpress = function(input, path = getwd(), save = FALSE, rawcol = NULL)
  {
    ## Extracting the data if the checking was fine
    extract = getAE(input, path, type = "raw")
  
    raweset = magetab2bioc(files = extract,
      rawcol = rawcol,
      save = save)
    
    if(inherits(raweset, 'try-error') && length(grep("following", raweset))!=0)
    save = TRUE

    if(!save) on.exit(cleanupAE(extract))

    if(!inherits(raweset, 'try-error'))
      message(paste("\n The object containing experiment ", input," has been built.\n"))
    
    return(raweset)
  }#end of ArrayExpress


