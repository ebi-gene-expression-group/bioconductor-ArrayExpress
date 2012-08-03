procset = function(files, procol){
	
  stopifnot(length(procol)==1)
  with(files, 
   {
	
	if(length(processedFiles) > 1)
	   stop("The processed files contain different numbers/subsets of reporters and cannot be automatically assembled.")
	
   #Read Processed data file
	if(length(processedFiles) == 1){
        proctot = try(read.table(file.path(path,processedFiles),header=TRUE,sep="\t", quote="", row.names=1, stringsAsFactors=FALSE))
        if(inherits(proctot, 'try-error'))
          if(length(grep("duplicate",proctot)) != 0) 
			  stop("The probe identifiers are not unique. The processed file cannot automatically be treated.")
	  	  else stop("Cannot read the processed file automatically.") 
    }
      

    #Create matrix of data and construct an expressionSet
    procsel = matrix(as.numeric(as.matrix(proctot[-1, procol == proctot[1,]])), nrow=nrow(proctot)-1)
    colnames(procsel) = gsub(".[a-z][a-z][a-z]$","",colnames(proctot[, procol == proctot[1,]]),ignore.case=T)
    rownames(procsel) = rownames(proctot[-1,])
    proceset = new("ExpressionSet", exprs = procsel)

	
	#read and attach sample data
	samples = readPhenoData(sdrf,path);
	if(inherits(samples, 'try-error'))
		warning("ArrayExpress: Parsing SDRF failed. Please make sure SDRF file ",sdrf," exists in ",path, " and is not corrupt.")
	else{
		if(all(sampleNames(proceset) %in% rownames(pData(samples)))){
			samples = samples[sampleNames(proceset),]
			phenoData(proceset) = samples
		} else warning("Cannot attach phenoData")
	}

	
	#Choose ADF ref column in processed data matrix to link to from ADF
	ADFrefCol = rownames(proctot)[1]
	ADFrefCol = gsub(" ",".",gsub(" REF"," Name",ADFrefCol))
	
	#read and attach feature data
	features= try(readFeatures(adf=adf,path=path, procADFref=ADFrefCol))
	if(inherits(features, "try-error")){
		warning("ArrayExpress: Unable to attach feature data")
	}else{
		if(all(featureNames(proceset) %in% featureNames(features))){
			features = features[featureNames(proceset),]
			featureData(proceset) = features	
		} else warning("Cannot attach feature data")
	}
		
	#read and attach experiment meta data
	experimentData = try(readExperimentData(idf=idf,path=path))
	if(inherits(experimentData, "try-error")){
		warning("ArrayExpress: Unable to read experiment data");
		experimentData = new("MIAME");
	}
	experimentData(proceset) = experimentData;
	
	if(!validObject(proceset)) 
		warning(validObject(proceset))
	else
		message(paste("\n",gsub(".sdrf.txt","",sdrf)," processed data was successfully loaded into ",class(proceset),"\n"))
	 
	return(proceset)
  }) ## with
}
