# TODO: Add comment
# 
# Author: iemam
###############################################################################

headers<-list(
		ae1=c("metaColumn","metaRow","row","column"),
		genepix=c("Block","Column","Row","X","Y"),
		arrayvision=c("Primary","Secondary"),
		agilent=c("Row","Col","PositionX","PositionY"),
		scanalyze=c("GRID","COL","ROW","LEFT","TOP","RIGHT","BOT"),
		scanarray=c('Array Column','Array Row','Spot Column','Spot Row','X','Y'),
		quantarray=c('Array Column', 'Array Row', 'Column', 'Row'),
		spotfinder=c("MC","MR","SC","SR","C","R"),
		mev=c("MC","MR","C","R","UID"),
		codelink=c("Logical_row","Logical_col","Center_X","Center_Y"),
		bluefuse=c("COL","ROW","SUBGRIDCOL","SUBGRIDROW"),
		UCSFSpot=c("Arr-colx","Arr-rowy","Spot-colx","Spot-rowy"),
		NimbleScanFeature=c("X","Y","PROBE_ID","X_PIXEL","Y_PIXEL"),
		NimblegenNASA=c("X_BC","Y_BC","Feature_ID","ProbID_BC"),
		imagene=c('Meta Column', 'Meta Row', 'Column', 'Row', 'Field', 'Gene ID'),
		ImaGene3=c("Meta_col","Meta_row","Sub_col","Sub_row","Name","Selected"),
		ImaGene7=c("Block","Column","Row","Ch1 XCoord","Ch1 YCoord", "Ch2 XCoord", "Ch2 YCoord"),
		ImaGeneFields=c("Field","Column","Row","XCoord","YCoord"),
		CSIRO_Spot=c("grid_c","grid_r","spot_c","spot_r","indexs")

#add illumina
)

isOneChannel = function(sdrf,path){
	ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
	labelCol = getSDRFcolumn("label",varLabels(ph))
	return(length(unique(tolower(ph[[labelCol]])))==1)
}

readPhenoData = function(sdrf,path){
	
	message("ArrayExpress: Reading pheno data from SDRF")
	ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
		
	#Remove empty rows from pheno data
	emptylines = which(sapply(seq_len(nrow(pData(ph))), function(i) all(pData(ph)[i,] == "",na.rm = TRUE)))
	if(length(emptylines) != 0)
		pData(ph) = pData(ph)[-emptylines,]
		
	phenoData = pData(ph)
	arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))
	labelCol = getSDRFcolumn("label",varLabels(ph))
	
	if(length(arrayDataCol)==0 || length(labelCol)==0)
		if(length(arrayDataCol)==0)
			warning("ArrayExpress: Cannot find 'Array Data File' column in SDRF. Object might not be created correctly.")
		if(length(labelCol)==0)
			warning("ArrayExpress: Cannot find 'Label' column in SDRF. Object might not be created correctly.")
		
	if(length(arrayDataCol)!=0 && length(labelCol)!=0){
		#filter out duplicated rows where multiple derived data files are available per one array data file
		ph=ph[!duplicated(phenoData[,c(arrayDataCol,labelCol)])]
	}
	
	
	if(length(arrayDataCol)!=0 & length(unique(tolower(ph[[labelCol]])))==1)
		#set rownames of phenoData annotated data frame to array data files			
		rownames(pData(ph)) = gsub(".[a-z][a-z][a-z]$","",ph[[arrayDataCol]],ignore.case=T)
	
	
	#treat SDRF for two channel experiments
	if(length(unique(tolower(ph[[labelCol]])))==2){
		
		arrayFilesNum = length(unique(ph[[arrayDataCol]]))
		
		si = pData(ph)[1:(arrayFilesNum*2),]
		lab = split(si,si[,"Label"])
		
		if(nrow(lab[[1]]) != nrow(lab[[2]])){
			stop("Number of CY3/CY5 is not equal")
		}
		
		#Reorder rows in each group (Cy3,Cy5) to the same order
		lab[[1]] = lab[[1]][order(lab[[1]][,arrayDataCol]),]
		lab[[2]] = lab[[2]][order(lab[[2]][,arrayDataCol]),]
		
		same = which(lapply(1:ncol(lab[[1]]), function(i) all(lab[[1]][i] == lab[[2]][i])) == TRUE)
		all = lab[[1]][same]
		gspe = lab[[1]][-same]
		colnames(gspe) = paste(colnames(gspe),names(lab)[1],sep = ".")
		rspe = lab[[2]][-same]
		colnames(rspe) = paste(colnames(rspe),names(lab)[2],sep = ".")
		
		metaData = data.frame(labelDescription = c(rep("_ALL_",ncol(all)),rep("G",ncol(gspe)),rep("R",ncol(rspe))))
		ph = new("AnnotatedDataFrame", data = cbind(all,gspe,rspe), varMetadata = metaData)
		
		arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))	
		rownames(pData(ph)) = gsub(".[a-z][a-z][a-z]$","",ph[[arrayDataCol]],ignore.case=T)
	}
		
		
#		if(!all(basename(files) %in% pData(ph)[,getSDRFcolumn("ArrayDataFile",varLabels(ph))]))
#			warning("Some data files in the zip archive are missing from the SDRF. The object may not be built.")
	return(ph)
}

readAEdata = function(path,files,dataCols,green.only){
	
	message("ArrayExpress: Reading data files")
	source = getDataFormat(path,files)

#process extra arguments passed to read.maimages
#		if (!is.null(dataCols)) {
#		if (length(rawcol) > 1 && !is(rawcol, "list")) 
#			stop("The argument 'rawcol' must be a list if multiple column name are given.")
#		if (length(rawcol) == 1 && !is(rawcol, "character")) 
#			stop("The argument 'rawcol' must be a character if one column name is given.")
#		if (is(rawcol, "list") && !("R" %in% names(rawcol) && "G" %in% names(rawcol))) 
#			stop("The names of the columns must contain R and G.")
#	}
	
	if(source == "affy"){
		rawdata = try(ReadAffy(filenames = file.path(path,unique(files))))
		if(inherits(rawdata, 'try-error')){
			stop("Unable to read cel files in",path)
		}
		return(rawdata)
		
	}else if(source=="ae1"){
		#Old AE1 formatted data files
		if (is.null(dataCols)){
			dataCols= try(getDataColsForAE1(path,files))
			if(inherits(dataCols,'try-error')) return()
		}
		rawdata = try(read.maimages(files=files,path=path,source="generic",columns=dataCols,annotation=headers$ae1))
		if(!inherits(rawdata, 'try-error')) rawdata$source="ae1"

	}else if(source %in% c("agilent","arrayvision","bluefuse","genepix","bluefuse","imagene","quantarray")){
		#data format can be directly read by limma
		rawdata = try(read.maimages(files=files,path=path,source=source,columns=dataCols,green.only=green.only))
		
	}else if(!is.null(dataCols) && source=="unknown"){
		#read generic source given columns specified by user
		rawdata = try(read.maimages(files=files,path=path,source="generic",columns=dataCols,green.only=green.only))
		
	}else
		# TODO: add more formats (illumina)
		stop("Unable to recognize data file format")
		
	if(inherits(rawdata, 'try-error')){
		stop("Unable to read data files in",path)
	}
	
	#reorder rows of RGList to reflect the same order in ADF as Block Row/Block Column/Row/Column
	if(is.null(rawdata$genes))
		stop("Unable to read probe annotation from RGList")
	if(is.null(rawdata$source))
		stop("Unable to read source from RGList")
	
	rawdata<-switch(source,
			agilent = rawdata<-rawdata[with(rawdata$genes,order(Row,Col)),],
			genepix = rawdata<-rawdata[with(rawdata$genes,order(Block,Row,Column)),],
			ae1 = rawdata<-rawdata[with(rawdata$genes,order(metaRow,metaColumn,row,column)),],
			rawdata)
	#if generic user has to specify feature columns in data
	
	
	return(rawdata)
}

readFeatures<-function(adf,path,procADFref=NULL){
	
	message("ArrayExpress: Reading feature metadata from ADF")
		
	lines2skip = skipADFheader(adf,path,!is.null(procADFref))
	features = try(read.table(file.path(path, adf), row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, sep="\t", na.strings=c('?','NA'), skip = lines2skip, header=TRUE, quote=""))
	
	if('Block.Column' %in% colnames(features) & 'Reporter.Name' %in% colnames(features)){
		ommittedRows = which(is.na(features[,'Block.Column']) | is.na(features[,'Reporter.Name']))
		if(length(ommittedRows)!=0){
			message("ArrayExpress: Ommitting NA rows from ADF")
			features = features[-ommittedRows,]
		}
	}
	

	
	if(!is.null(procADFref)){
		if(procADFref %in% colnames(features))
			rownames(features) = features[,procADFref]
		else{
			repCol = getSDRFcolumn("reporter",colnames(features))
			if(length(repCol) != 0)
				repCol= repCol[1]
			else
				repCol = getSDRFcolumn("composite",colnames(features))
			if(length(repCol) != 0)
				repCol = repCol[1]
			else
				repCol = NULL
			if(!is.null(repCol))
				rownames(features) = features[[repCol]]
		}
	}
	

	
	#Sort ADF features by columns Block row/Block column/Row/Column (only applicable for raw data exps, processed data is ordered by reporter/composite name)
	if("Block.Row" %in% colnames(features))
		features = features[with(features,order(Block.Row,Block.Column,Row,Column)),]
	
	
	#Row names of featureData must match row names of the matrix / matricies in assayData
#	ri1 = grep("reporter.identifier|reporter.name", colnames(adff), ignore.case=TRUE)
#	ri2 = grep("reporter.identifier|reporter.name", colnames(fn), ignore.case=TRUE)	
#	if(all(adff2[,ri1] == fn[,ri2])) 
#		featureData(eset) = new("AnnotatedDataFrame",adff2) 
#	else stop("Do not manage to map the reporter identifier between the annotation and the data files.\n")

	return(new("AnnotatedDataFrame",features))
}

## Assign experiment Data
## By Juok Cho
readExperimentData = function(idf, path){
	idffile = scan(file.path(path,idf),character(),sep = "\n")
	idf.data = list()
	for(g in idffile) { 
		e = unlist(strsplit(g,"\t"))
		key = e[1] 
		if(length(e)>1)
			values = e[2:length(e)]
		else
			values = NA
		idf.data[[key]] = values
	}
	
	## making key matches #
	## (Person Last Name, Person First Name, Person Mid Initials), Person Email, Person Phone, Person Address, Person Affiliation, 
	Person_Name = c(idf.data$"Person First Name", idf.data$"Person Last Name",idf.data$"Person Mid Initials") 
	Personal_contact = c(idf.data$"Person Email", idf.data$"Person Phone", idf.data$"Person Address")
	
	## making experimentData object #		
	SubmitterIndex = which(idf.data$"Person Roles"=="submitter")
	experimentData = new("MIAME", 
			name = as.character(paste(idf.data$"Person Last Name"[SubmitterIndex],", ",idf.data$"Person First Name"[SubmitterIndex], sep = "")), #performer
			lab = as.character(idf.data$"Person Affiliation"[SubmitterIndex]) , #Person Affiliation 
			contact = as.character(idf.data$"Person Email"[SubmitterIndex]), # Person Email(Person Phone, Person Address)
			title = as.character(idf.data$"Investigation Title") , #description #Investigation Title
			##abstract= "",	#not provided in the idf.data
			##url	= "",
			other = list(
					accession = gsub(".sdrf.txt","",idf.data$"SDRF File"), #Experiment Name
					identifier = gsub(".sdrf.txt","",idf.data$"SDRF File"), #Experiment Name
					##Experimental Factor Type
					experimentalFactor = c(idf.data$"Experimental Factor Type"), 
					##Experimental Design
					type = c(idf.data$"Experimental Design")
					#measurementType = experimentData(eset)@other$measurementType #from processed data.zip depending on user answer about QT type
			)
	)
	#experimentData(eset) = experimentData
	return(experimentData)	  
}

skipADFheader<-function(adf,path,proc=F){
	if(!proc)
		columns = list('Block Column','Block Row','Column','Row')
	else
		columns = list('Composite Element Name')
	
	con = file(file.path(path, adf), "r")	
	on.exit(close(con))
	
	Found = FALSE
	i = 0
	repeat {
		i = i+1
		txt <- readLines(con,n=1)
		if(!length(txt))
			stop("Failed to recognize ADF file format")
		Found = TRUE
		for(a in columns) 
			Found = Found && length(grep(a,txt))
		Found2 = length(grep("^Reporter[[:punct:]|[:blank:]]*Name",txt,ignore.case=TRUE))
		if(Found || Found2)
			break
	}
	return(i-1)
}

getPhenoDataPerAD<-function(ad,ph,dataFiles){
	phenoData = pData(ph)
	arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))
	arrayDesignRefCol = getSDRFcolumn("ArrayDesignREF",varLabels(ph))
	
	if(length(arrayDataCol)==0)
		stop("Cannot find 'Array.Data.File' column in SDRF. Experiment uses multiple array designs. Cannot distinguish arrays with similar array design.")
	
	if(length(arrayDesignRefCol)==0)
		stop("Cannot find 'Array.Design.REF' column in SDRF. Experiment uses multiple array designs. Cannot distinguish arrays with similar array design.")
	
	if(length(phenoData[phenoData[arrayDesignRefCol]==ad,arrayDataCol]) == 0)
		stop("Cannot find array data file names in the sdrf file. Experiment uses multiple array designs. Cannot distinguish arrays with similar array design.") 
	
	#Subselect data files for current ArrayDesign REF
	selectFiles = phenoData[phenoData[arrayDesignRefCol]==ad,arrayDataCol]
	
	
	if(!all(selectFiles %in% dataFiles)){
		stop("Some or all data files for ",ad," array are missing.")
	}
	
	#subselect phenoData Frame for files 
	ph = ph[phenoData[arrayDesignRefCol]==ad]
	phenoData = pData(ph)
	
	return(list(pheno=ph,dataFiles=selectFiles))
}

getDataFormat=function(path,files){
	
	if(length(grep(".cel",files, ignore.case = TRUE)) == length(files)){
		return("affy")
	}
	else{
		#Retrieve Column names from first data file
		allcnames = scan(file.path(path,files[1]),what = "",nlines = 200, sep = "\t",quiet=TRUE)
		allcnamesL = try(tolower(gsub("^\\s","",allcnames)),silent=TRUE)
		if(inherits(allcnamesL, 'try-error')){
			allcnames = scan(file.path(path,files[1]),what = "",nlines = 200, sep = "\t",quiet=TRUE,encoding="latin1")
			allcnamesL = try(tolower(gsub("^\\s","",allcnames)))
			if(inherits(allcnamesL, 'try-error')){
				allcnamesL = allcnames
			}
		}
		
		#Find source of data
		for(source in names(headers)){
			allthere<-tolower(headers[[source]]) %in% allcnamesL
			if(all(allthere))
				return(source)
		}
	}

	#Unable to detect source
	allcnames = scan(file.path(path,files[1]),what = "",nlines = 1, sep = "\n",quiet=TRUE)
	return("unknown")
	
	#stop("Unable to recognize data file format. First line:\n",allcnames)
	
}

getDataColsForAE1 = function(path,files){
	url2 = "http://tab2mage.svn.sourceforge.net/viewvc/tab2mage/trunk/Tab2MAGE/lib/ArrayExpress/Datafile/QT_list.txt" 
	
	qt = try(read.table(url2, sep = "\t", quote = "",
					check.names = FALSE, fill = TRUE,
					comment.char = "#",               
					stringsAsFactors =  FALSE)) ##read the QT file from the web
	if(inherits(qt, "try-error"))
		qt = try(read.table(
						file.path(system.file("doc", package = "ArrayExpress"),"QT_list.txt"),
						sep = "\t", quote = "",
						check.names = FALSE, fill = TRUE,
						comment.char = "#",               
						stringsAsFactors =  FALSE)) ##read the QT file from the package
	
	
	scanners = grep(">>>",qt[,1],value = TRUE) ## list all the scanner names
	sl = grep(">>>",qt[,1]) ## list all the line numbers wherea scanner type starts
	scanners = gsub(">","",scanners)
	
	## Parsing the first line of the expression file
	allcnames = scan(file.path(path,files[1]),what = "",nlines = 1, sep = "\t",quiet=TRUE)
	
	## Looking for the right column to use
	scanname = allcnames[grep(":",allcnames)]
	if(length(grep("Database|Reporter",scanname)) != 0)
		scanname = scanname[-grep("Database|Reporter",scanname)] 
	
	##Feature is a problem because of Feature Extraction
	feature = grep("^Feature^",scanname)
	fe = grep("Feature Extraction",scanname)
	if(length(feature) != 0)
		scanname = scanname[-feature[!feature %in% fe]]
	
	##Ready to read data
	if(length(scanname) == 0) 
		stop(sprintf("No scanner name is given. It is not possible to handle such a case. Try to set the argument 'dataCols' by choosing among the following columns names: \n") ,
				sprintf("\"%s\" \n",allcnames))
	
	##Image Analysis Program
	software = unique(sapply(seq_len(length(scanname)), function(i) strsplit(scanname,":")[[i]][1]))     
	st = NULL

    for(x in software){
      if(length(grep(x, scanners)) != 0){
        st = x
        break
      } 
    }

	#if(length(grep(st, scanners)) == 0)
	if(is.null(st))
    	stop(sprintf("Scanner name is, ",software,". This scanner type is not valid. \nTry to set the argument 'dataCols' by choosing among the following columns names: \n", st),sprintf("\"%s\" \n",scanname))
	
	if(length(st) != 1)
		stop(sprintf("%s scanner names are given ( ",length(st)), sprintf("\"%s\" ",st), sprintf("). It is not possible to handle such a case. Try to set the argument 'dataCols' by choosing among the following columns names: \n") ,sprintf("\"%s\" \n",scanname))
	
	if(length(grep(st, scanners)) > 1)
		stop(sprintf("Scanner name can be '%s'. \nTry to set the argument 'dataCols' by choosing among the following columns names: \n", scanners[grep(st, scanners)]),sprintf("\"%s\" \n",scanname))
	
	gs = qt[((sl[grep(st,scanners)]+1):(sl[grep(st,scanners)+1]-1)),] ## extract the QTs of the specific scanner type
	foreground = gs[(gs[,4] == "MeasuredSignal" & (is.na(gs[,5]) | gs[,5] == 0)),c(1,7)] ## the colnames to use in the read.column
	background = gs[(gs[,4] == "MeasuredSignal" & (gs[,5] == 1)),c(1,7)] ## the colnames to use in the read.column
	
	foreground[,1] = paste(st,":",foreground[,1],sep = "")
	colnamesf = foreground[which(foreground[,1] %in% allcnames),]
	df = dim(colnamesf)        
	
	if(dim(background)[1] != 0){
		background[,1] = paste(st,":",background[,1],sep = "")
		colnamesb = background[which(background[,1] %in% allcnames),]
		db = dim(colnamesb)
	}else 
		db = 0
	
	if(length(files) != 1){
		if(!all(sapply(2:length(files), function(i) readLines(file.path(path,files[1]),1) == readLines(file.path(path,files[i]),1))))
			warning(sprintf("The files do not all have the same headings whereas the array design is the same. It may cause the object not being created."))
	}
	
	#Two channel data
	if(df[1] == 2){
		rawcol = if(db[1] == 2) 
					list(R = colnamesf[colnamesf[,2] == "Cy5",1], 
						 G = colnamesf[colnamesf[,2] == "Cy3",1],
						 Rb = colnamesb[colnamesb[,2] == "Cy5",1],
						 Gb = colnamesb[colnamesb[,2] == "Cy3",1])
			 	 else 
					list(R = colnamesf[colnamesf[,2] == "Cy5",1],
						 G = colnamesf[colnamesf[,2] == "Cy3",1])
		
		if(length(rawcol) == 0 || (0 %in% sapply(seq_len(length(rawcol)), function(i) length(rawcol[[i]]))))
			stop(sprintf("The known column names for this scanner are not in the heading of the files.\nTry to set the argument 'rawcol' by choosing among the following columns names: \n"),
					sprintf("\"%s\" \n",scanname))
	}
	## one channel data
	if(df[1] == 1){
		rawcol = if(db != 0)
					list(G = colnamesf[,1],
						 Gb = colnamesb[,1])
			 	 else
					 list(G = colnamesf[,1])
	}
	
	if(df[1] == 0)
		stop(sprintf("None of the columns names of the expression files is matching a valid known quantitation type.\nTry to set the argument 'rawcol' by choosing among the following columns names: \n"), sprintf("\"%s\" \n",scanname))
	if(df[1] > 2)
		stop(sprintf("There are too many columns that could be read in the files.\n Try to set the argument 'rawcol' by choosing among the following columns names: \n"),sprintf("\"%s\" \n",scanname))
	
	return(rawcol)		
}

## assign phenoData to Nchannelset
preparePhenoDataFor2channel = function(ph,files){
	#if(length(unique(tolower(ph[[labelCol]])))==2){
		
		arrayFilesNum = length(unique(ph[[arrayDataCol]]))
		
		si = pData(ph)[1:(arrayFilesNum*2),]
		lab = split(si,si[,"Label"])
		
		if(nrow(lab[[1]]) != nrow(lab[[2]])){
			stop("Number of CY3/CY5 is not equal")
		}
		
		#Reorder rows in each group (Cy3,Cy5) to the same order
		lab[[1]] = lab[[1]][order(lab[[1]][,arrayDataCol]),]
		lab[[2]] = lab[[2]][order(lab[[2]][,arrayDataCol]),]
		
		same = which(lapply(1:ncol(lab[[1]]), function(i) all(lab[[1]][i] == lab[[2]][i])) == TRUE)
		all = lab[[1]][same]
		gspe = lab[[1]][-same]
		colnames(gspe) = paste(colnames(gspe),names(lab)[1],sep = ".")
		rspe = lab[[2]][-same]
		colnames(rspe) = paste(colnames(rspe),names(lab)[2],sep = ".")
		
		metaData = data.frame(labelDescription = c(rep("_ALL_",ncol(all)),rep("G",ncol(gspe)),rep("R",ncol(rspe))))
		ph = new("AnnotatedDataFrame", data = cbind(all,gspe,rspe), varMetadata = metaData)
		
		arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))	
		rownames(pData(ph)) = gsub(".[a-z][a-z][a-z]$","",ph[[arrayDataCol]],ignore.case=T)
	#}
	return(ph)
}

getSDRFcolumn = function(col,headers){
	pattern<-switch(col,
			ArrayDataFile = "^Array[[:punct:]|[:blank:]]*Data[[:punct:]|[:blank:]]*File",
			ArrayDesignREF = "^Array[[:punct:]|[:blank:]]*Design[[:punct:]|[:blank:]]*REF",
			ArrayDataMatrixFile = "^Array[[:punct:]|[:blank:]]*Data[[:punct:]|[:blank:]]*Matrix[[:punct:]|[:blank:]]*File",
			label = "^Label$",
			factorValues = "^Factor[[:punct:]|[:blank:]]*Value",
			DerivedArrayMatrix = "^Derived[[:punct:]|[:blank:]]*Array[[:punct:]|[:blank:]]*Data[[:punct:]|[:blank:]]*Matrix[[:punct:]|[:blank:]]*File",
			DerivedArrayFile = "^Derived[[:punct:]|[:blank:]]*Array[[:punct:]|[:blank:]]*Data[[:punct:]|[:blank:]]*File",
			reporter = "^Reporter[[:punct:]|[:blank:]]*[Name | Identifier]",
			composite = "^Composite[[:punct:]|[:blank:]]*Element[[:punct:]|[:blank:]]*[Name | Identifier]")
	colIndex = grep(pattern,headers,ignore.case = TRUE)
	return(colIndex)
}

## remove all downloaded files
cleanupAE = function(mageFiles){
	path = mageFiles$path
	
	try(file.remove(file.path(path, mageFiles$rawFiles)))
	try(file.remove(file.path(path, mageFiles$processedFiles)))
	
	try(file.remove(file.path(path, mageFiles$sdrf)))
	try(file.remove(file.path(path, mageFiles$idf)))
	try(file.remove(file.path(path, mageFiles$adf)))
	try(file.remove(file.path(path, mageFiles$rawArchive)))
	try(file.remove(file.path(path, mageFiles$processedArchive)))
}


