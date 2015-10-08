ArrayExpress = function(accession, path = tempdir(), save = FALSE, dataCols = NULL, drop = TRUE){
	
	expFiles = getAE(accession, path, type = "raw")
	
#	if(!save) 
#		on.exit(cleanupAE(expFiles))
	
	raweset = try(ae2bioc(mageFiles = expFiles, dataCols = dataCols))
	
	if(inherits(raweset, 'try-error')){
		save = TRUE
		raweset = NULL
	}
	else{
		if(!is.list(raweset)){
#			raweset = raweset[[1]];
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
