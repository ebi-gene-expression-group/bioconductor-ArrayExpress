procset = function(files, procol){
	
	stopifnot(length(procol)==1)
	with(files,
			{
				#ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
				#Read SDRF
				samples = readPhenoData(sdrf,path);
				if(!inherits(samples, 'try-error')){
					derivedMatrixCol = getSDRFcolumn("DerivedArrayMatrix",varLabels(samples))
					derivedFileCol = getSDRFcolumn("DerivedArrayFile",varLabels(samples))	
				}
				else
					warning("Unable to read SDRF")
				
				if(length(processedFiles) > 1 && length(derivedFileCol) != 0){
					#stop("The processed files contain different numbers/subsets of reporters and cannot be automatically assembled.")
					# READ individual files for each assay
					proceset = readDerivedDataFiles(processedFiles,procol,path)
					assayCol = derivedFileCol
					ADFrefCol = "Reporter.Name"
				}
				else if(length(processedFiles) == 1 && length(derivedMatrixCol) != 0){
					#Read Processed data file
					res = readDerivedDataMatrixFile(processedFiles,procol,path)
					proceset = res$proceset
					assayREF = res$assayREF
					assayCol = grep(assayREF,varLabels(samples))
					ADFrefCol = res$ADFrefCol
				}		
				else
					stop("Unable to read processed data")

				
				#attach sample data
				dataSamples = gsub("\\.[a-z][a-z][a-z]$","",sampleNames(proceset),ignore.case=T)
				SDRFsamples = gsub("\\.[a-z][a-z][a-z]$","",pData(samples)[[assayCol]])
				
				if(all(dataSamples %in% SDRFsamples)){
					#if(all(dataSamples %in% rownames(pData(samples)))){
					rownames(pData(samples)) = SDRFsamples
					samples = samples[dataSamples,] #Reorder samples in phenoData to match those in data
					phenoData(proceset) = samples
				} else warning("Cannot attach phenoData")

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

readDerivedDataFiles = function(processedFiles,procol,path){

	for(procFile in processedFiles){
		cat("Reading processed File: ",procFile,"\n")
		procData <- try(read.table(file.path(path, procFile), sep = "\t", 
						header = TRUE, # interpret header
						row.names = 1, # header located in line 1
						stringsAsFactors = FALSE, 
						na.strings = c('NA','NULL','null')));
		
		procsel = matrix(as.numeric(as.matrix(procData[, procol == colnames(procData)])), nrow = nrow(procData))
		
		if(grep(procFile,processedFiles) == 1){
			E = matrix(procsel)
			rownames(E) = rownames(procData)
			
		}else
			E = cbind(E,procsel)
	}
	
	###TEMP###
	colnames(E) = gsub("\\.txt","",processedFiles)
	##########
	
	proceset = new("ExpressionSet", exprs = E)
	return(proceset)	
}

readDerivedDataMatrixFile = function(processedFile,procol,path){
	#Convention for ArrayDataMatrix file
	# First row contains assay refs (eg. SCAN REF)
	# Second row contains gene annotation values (eg. GEO:AFFYMETRIX_VALUE)
	# avoided using header=T because of automatic name 
	
	
	cat("Reading processed data matrix file, ",file.path(path,processedFile), "\n")
	#READ HEADER
	matrix.header <- try(read.table(file.path(path, processedFile), sep = "\t",
					header = FALSE, # read in and interprete header
					row.names = 1, # header located in row 1
					nrows = 2,     # num rows after header
					stringsAsFactors = FALSE));
	
	if (inherits(matrix.header, "try-error")) {
		stop("Cannot read header ", file.path(path, processedFile));
	}
	
	#READ a few rows to determine column classes for faster read.table
	proctot.class <- try(read.table(file.path(path, processedFile), sep = "\t",
					header = TRUE, # interprete header
					skip = 1,      # starting from second line
					row.names = 1, # header located in line 1
					nrows = 5,     # read in 5 lines
					stringsAsFactors = FALSE));        
	
	if (inherits(proctot.class, "try-error")) {
		stop("Cannot read classes ", file.path(path, processedFile));
	}
	
	classes <- sapply(proctot.class, class) # figure out classes
	
	#now read all table
	data.matrix <- try(read.table(file.path(path, processedFile), sep = "\t", 
					header = TRUE, # interpret header
					skip = 1,      # starting from second line
					row.names = 1, # header located in line 1
					stringsAsFactors = FALSE, 
					na.strings = c('NA','NULL','null'),
					colClasses = classes));
	
	if(inherits(data.matrix, 'try-error'))
		if(length(grep("duplicate",data.matrix)) != 0) 
			stop("The probe identifiers are not unique. The processed file cannot automatically be treated.")
		else stop("Cannot read the processed file automatically.")

	#Create matrix of data and construct an expressionSet
	datasel = matrix(as.numeric(as.matrix(data.matrix[, procol == matrix.header[2, ]])), nrow = nrow(data.matrix))
	#colnames(datasel) = colnames(proctot.header[, procol == proctot.header[1, ]])
	colnames(datasel) = matrix.header[, procol == matrix.header[2, ]][1,]
	rownames(datasel) = rownames(data.matrix)
	
	proceset = new("ExpressionSet", exprs = datasel)
	
	
	assayREF = rownames(matrix.header)[1]
	assayREF = gsub('\\sREF','',assayREF)
	
	#Choose ADF ref column in processed data matrix to link to from ADF
	ADFrefCol = rownames(matrix.header)[2]
	ADFrefCol = gsub(" ",".",gsub(" REF"," Name",ADFrefCol))
	
	res = list()
	res$proceset = proceset
	res$assayREF = assayREF
	res$ADFrefCol = ADFrefCol
	
	return(res)	
}
