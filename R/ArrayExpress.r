ArrayExpress = function(input, path = getwd(), save = FALSE, dataCols = NULL){
	
	expFiles = getAE(input, path, type = "raw")
	
#	if(!save) 
#		on.exit(cleanupAE(expFiles))
	
	raweset = try(ae2bioc(mageFiles = expFiles, dataCols = dataCols))
	
	if(inherits(raweset, 'try-error')){
		save = TRUE
		raweset = expFiles
	}
	else{
		if(length(raweset)==1){
			raweset = raweset[[1]];
			message(paste("\n",input," was successfully loaded into ",class(raweset),"\n"))
		}
		else if(length(raweset)>1){
			message(paste("\n",input," was successfully loaded into ",length(raweset)," ",class(raweset[[1]]),"(s)\n"))
		}
	}
		
	
	if(!save) 
		on.exit(cleanupAE(expFiles))
	
	return(raweset)
}