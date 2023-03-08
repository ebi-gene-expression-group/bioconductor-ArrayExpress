ArrayExpress = function(accession, path = tempdir(), save = FALSE, dataCols = NULL, drop = TRUE){
  
  expFiles = getAE(accession, path = path, type = "full")
  
  raweset = try(ae2bioc(mageFiles = expFiles, dataCols = dataCols, drop=drop))
  
  if(inherits(raweset, 'try-error')){
    save = TRUE
    raweset = NULL
  }
  else{
    if(!is.list(raweset)){
      message(paste("\n",accession," was successfully loaded into ",class(raweset),"\n"))
    }
    else {
      message(paste("\n",accession," was successfully loaded into ",length(raweset)," ",unlist(lapply(raweset,function(x){class(x)})),"\n"))
    }
  }
  
  
  if(!save) 
    on.exit(cleanupAE(expFiles))

  return(raweset)
}