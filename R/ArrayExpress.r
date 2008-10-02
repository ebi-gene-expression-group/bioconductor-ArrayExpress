ArrayExpress = function(input, path = ".", save = TRUE, columns = NULL, type = "full")
  {
    if(!save) on.exit({try(file.remove(file.path(path,extract$sdrf))); try(file.remove(file.path(path,extract$idf)))})
   
    ## Extracting the data if the checking was fine
    extract = getAE(input, path, save, type)
  
    if(type == "full" || type == "raw")
      raweset = magetab2bioc(rawfiles = extract$rawfiles,
        sdrf = extract$sdrf,
        idf = extract$idf,
        path = path,
        columns = columns,
        save = save)
    
    if(type == "full" || type == "processed")
      proceset = AEproc2bioc(procfiles = extract$procfiles,
        sdrf = extract$sdrf,
        idf = extract$idf,
        path = path,
        save = save)

    if(type == "full")
      return(list(raw = raweset, processed = proceset))
    
    if(type == "raw")
      return(raweset)
    
    if(type == "processed")
      return(proceset)

    
  }#end of ArrayExpress
