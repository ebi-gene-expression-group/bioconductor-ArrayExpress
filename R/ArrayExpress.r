ArrayExpress = function(input, path = getwd(), save = FALSE, dataCols = NULL){
	
	## Extracting the data if the checking was fine
	extract = getAE(input, path, type = "raw")
	
	raweset = ae2bioc(mageFiles = extract,
			dataCols = dataCols,
			save = save)
	
	if(inherits(raweset, 'try-error') && length(grep("following", raweset))!=0)
		save = TRUE
	
	if(!save) 
		on.exit(cleanupAE(extract))
	
	if(!inherits(raweset, 'try-error'))
		message(paste("\n The object containing experiment ", input," has been built.\n"))
	
	return(raweset)
}#end of ArrayExpress