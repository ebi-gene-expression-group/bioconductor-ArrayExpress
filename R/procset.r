procset = function(files, procol){
	
  stopifnot(length(procol)==1)
  with(files, 
		  {
			  if(length(processedFiles) > 1)
		   		stop("The processed files contain different numbers/subsets of reporters and cannot be automatically assembled.")
			
			#Read Processed data file
			if(length(processedFiles) == 1){
				message("ArrayExpress: Reading data file, ",file.path(path,processedFiles))
		
				#READ HEADER
				proctot.header <- try(read.table(file.path(path, processedFiles), sep = "\t",
						header = TRUE, # read in and interprete header
						row.names = 1, # header located in row 1
						nrows = 1,     # num rows after header
						stringsAsFactors = FALSE));
		
				if (inherits(proctot.header, "try-error")) {
					stop("Cannot read header ", file.path(path, processedFiles));
				}
		
				#READ a few rows to determine column classes for faster read.table
				proctot.class <- try(read.table(file.path(path, processedFiles), sep = "\t",
								header = TRUE, # interprete header
								skip = 1,      # starting from second line
								row.names = 1, # header located in line 1
								nrows = 5,     # read in 5 lines
								stringsAsFactors = FALSE));        
				
				if (inherits(proctot.class, "try-error")) {
					stop("Cannot read classes ", file.path(path, processedFiles));
				}
		
				classes <- sapply(proctot.class, class) # figure out classes

				#now read all table
				proctot <- try(read.table(file.path(path, processedFiles), sep = "\t", 
											header = TRUE, # interpret header
											skip = 1,      # starting from second line
											row.names = 1, # header located in line 1
											stringsAsFactors = FALSE, colClasses = classes));
		
		        if(inherits(proctot, 'try-error'))
		          if(length(grep("duplicate",proctot)) != 0) 
					  stop("The probe identifiers are not unique. The processed file cannot automatically be treated.")
			  	  else stop("Cannot read the processed file automatically.")
			  
			  	#Attach column header to data table
	  			colnames(proctot) <- colnames(proctot.header)
			}
      

		    #Create matrix of data and construct an expressionSet
			procsel = matrix(as.numeric(as.matrix(proctot[-1, procol == proctot.header[1, ]])), nrow = nrow(proctot) - 1)
			colnames(procsel) = colnames(proctot.header[, procol == proctot.header[1, ]])
			rownames(procsel) = rownames(proctot[-1, ])
		    proceset = new("ExpressionSet", exprs = procsel)

	
			#read and attach sample data
			samples = readPhenoData(sdrf,path);
			if(inherits(samples, 'try-error'))
				warning("ArrayExpress: Parsing SDRF failed. Please make sure SDRF file ",sdrf," exists in ",path, " and is not corrupt.")
			else{
				dataSamples = gsub(".[a-z][a-z][a-z]$","",sampleNames(proceset),ignore.case=T)
				if(all(dataSamples %in% rownames(pData(samples)))){
					samples = samples[dataSamples,] #Reorder samples in phenoData to match those in data
					sampleNames(proceset) = dataSamples #rename sample names to match phenoData samples
					phenoData(proceset) = samples
				} else warning("Cannot attach phenoData")
			}

	
			#Choose ADF ref column in processed data matrix to link to from ADF
			ADFrefCol = rownames(proctot.header)
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
