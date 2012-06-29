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

readPhenoData = function(sdrf,path){
	
	ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, varMetadata.char = "$", quote="\""))
		
	if(!inherits(ph, 'try-error')){
		#Remove empty rows from pheno data
		emptylines = which(sapply(seq_len(nrow(pData(ph))), function(i) all(pData(ph)[i,] == "",na.rm = TRUE)))
		if(length(emptylines) != 0)
			pData(ph) = pData(ph)[-emptylines,]
		
		phenoData = pData(ph)
		arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))
		labelCol = getSDRFcolumn("label",varLabels(ph))
		if(length(arrayDataCol)!=0){
			#filter out duplicated rows where multiple derived data files are available per one array data file
			ph=ph[!duplicated(phenoData[,c(arrayDataCol,labelCol)])]
			
			phenoData = pData(ph)
			#set rownames of phenoData annotated data frame to array data files
			#rownames(pData(ph)) = phenoData[,arrayDataCol]
			#rownames(pData(ph)) = gsub(".[a-z][a-z][a-z]$","",ph$Array.Data)
		}else
			warning("Cannot find 'Array Data File' column in SDRF. Object might not be created correctly.")
		
#		if(!all(basename(files) %in% pData(ph)[,getSDRFcolumn("ArrayDataFile",varLabels(ph))]))
#			warning("Some data files in the zip archive are missing from the SDRF. The object may not be built.")
		return(ph)
	}else
		return(NULL)
}

readAEdata = function(path,files,dataCols,...){
	
	message("ArrayExpress: Reading data files")
	if(length(grep(".cel",files, ignore.case = TRUE)) == length(files)){
		#Affy experiments
		rawdata = try(ReadAffy(filenames = file.path(path,unique(files)))) 
		if(inherits(rawdata, 'try-error')){
			stop("Unable to read Affy CEL files in ",path)
		}
	}else{
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
		
		#Old AE1 formatted data files
		if(source=="ae1"){
			if (is.null(dataCols))
				dataCols= getDataColsForAE1(path,files)
			rawdata = read.maimages(files=files,path=path,source="generic",columns=dataCols,annotation=headers$ae1)
			rawdata$source="ae1"
		}
		
		else if(source %in% c("agilent","arrayvision","bluefuse","genepix","bluefuse","imagene","quantarray")){
			#data format can be directly read by limma
			rawdata = read.maimages(files=files,path=path,source=source)
		}else
			# TODO: add more formats
			stop("Unable to read data file from ",source)
		
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
	}
	
	return(rawdata)
}

readFeatures<-function(rawdata,adf,path){
	
	if(class(rawdata) == "AffyBatch")
		return(rawdata)
	
	message("ArrayExpress: Reading feature metadata from ADF")
#	if(rawdata$source=="agilent"){
#		rawdata<-rawdata[with(rawdata$genes,order(Row,Col)),]		
#	}
#	if(rawdata$source=="genepix"){	
#		rawdata<-rawdata[with(rawdata$genes,order(Block,Row,Column)),]
#	}
#	if(rawdata$source=="ae1"){
#		rawdata<-rawdata[with(rawdata$genes,order(metaRow,metaColumn,row,column)),]
#	}
		
	#Sort ADF features by columns Block row/Block column/Row/Column
	lines2skip = skipADFheader(adf,path)
	features = try(read.table(file.path(path, adf), row.names = NULL, blank.lines.skip = TRUE, fill = TRUE, sep="\t", skip = lines2skip, header=TRUE, quote=""))
	features<-features[with(features,order(Block.Row,Block.Column,Row,Column)),]
	ommittedRows<-which(is.na(features[,'Reporter.Name']))
	
	if(length(ommittedRows)!=0){
		message("ArrayExpress: Ommitting NA rows from ADF")
		features<- features[-ommittedRows,]
	}
	
	#Row names of featureData must match row names of the matrix / matricies in assayData
	#rownames(features) = 
	
#	ri1 = grep("reporter.identifier|reporter.name", colnames(adff), ignore.case=TRUE)
#	ri2 = grep("reporter.identifier|reporter.name", colnames(fn), ignore.case=TRUE)
#	
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

skipADFheader<-function(adf,path){
	columns<-list('Block Column','Block Row','Column','Row')
	con <- file(file.path(path, adf), "r")	
	on.exit(close(con))
	
	Found <- FALSE
	i <- 0
	repeat {
		i <- i+1
		txt <- readLines(con,n=1)
		if(!length(txt)) stop("ADF file empty")
		Found <- TRUE
		for(a in columns) Found <- Found && length(grep(a,txt))
		if(Found) break
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
	
	if(!all(selectFiles) %in% dataFiles){
		stop("Some or all data files for ",ad," array are missing.")
	}
	
#	if("Array.Data.Matrix.File" %in% colnames(pht))
#		files = pht[pht$Array.Design.REF == adr[i],"Array.Data.Matrix.File"]
#	if("Array.Data.File" %in% colnames(pht))
#		files = pht[pht$Array.Design.REF == adr[i],"Array.Data.File"]
#	
	#subselect phenoData Frame for files 
	ph = ph[phenoData[arrayDesignRefCol]==ad]
	phenoData = pData(ph)
	
	return(list(pheno=ph,dataFiles=selectFiles))
}

getDataFormat=function(path,files){
	allcnames = scan(file.path(path,files[1]),what = "",nlines = 200, sep = "\t")
	for(sw in names(headers)){
		allthere<-tolower(headers[[sw]]) %in% tolower(gsub("\\s","",allcnames))
		if(all(allthere))
			return(sw)
	}
	return(FALSE)
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
	allcnames = scan(file.path(path,files[1]),what = "",nlines = 1, sep = "\t")
	
	## Looking for the right column to use
	scanname = allcnames[grep(":",allcnames)]
	if(length(grep("Database|Reporter",scanname)) != 0)
		scanname = scanname[-grep("Database|Reporter",scanname)] 
	
	##Feature is a problem because of Feature Extraction
	feature = grep("Feature",scanname)
	fe = grep("Feature Extraction",scanname)
	if(length(feature) != 0)
		scanname = scanname[-feature[!feature %in% fe]]
	
	##Ready to read data
	if(length(scanname) == 0) 
		stop(sprintf("No scanner name is given. It is not possible to handle such a case. Try to set the argument 'rawcol' by choosing among the following columns names: \n") ,
				sprintf("\"%s\" \n",allcnames))
	
	##Image Analysis Program
	st = unique(sapply(seq_len(length(scanname)), function(i) strsplit(scanname,":")[[i]][1]))     
	
	if(length(grep(st, scanners)) == 0)
		stop(sprintf("Scanner name is '%s'. This scanner type is not valid. \nTry to set the argument 'rawcol' by choosing among the following columns names: \n", st),sprintf("\"%s\" \n",scanname))
	
	if(length(st) != 1)
		stop(sprintf("%s scanner names are given ( ",length(st)), sprintf("\"%s\" ",st), sprintf("). It is not possible to handle such a case. Try to set the argument 'rawcol' by choosing among the following columns names: \n") ,sprintf("\"%s\" \n",scanname))
	
	if(length(grep(st, scanners)) > 1)
		stop(sprintf("Scanner name can be '%s'. \nTry to set the argument 'rawcol' by choosing among the following columns names: \n", scanners[grep(st, scanners)]),sprintf("\"%s\" \n",scanname))
	
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
	
	## Building NChannelSet when two colours
	if(df[1] == 2){
		rawcol = if(db[1] == 2) 
					list(R = colnamesf[colnamesf[,2] == "Cy5",1], 
							G = colnamesf[colnamesf[,2] == "Cy3",1],
							Rb = colnamesb[colnamesb[,2] == "Cy5",1],
							Gb = colnamesb[colnamesb[,2] == "Cy3",1]) 
				else 
					list(R = colnamesf[colnamesf[,2] == "Cy5",1], G = colnamesf[colnamesf[,2] == "Cy3",1])
		
		if(length(rawcol) == 0 || (0 %in% sapply(seq_len(length(rawcol)), function(i) length(rawcol[[i]]))))
			stop(sprintf("The known column names for this scanner are not in the heading of the files.\nTry to set the argument 'rawcol' by choosing among the following columns names: \n"),
					sprintf("\"%s\" \n",scanname))
	}
	## Building ExpressionSet when one colour
	if(df[1] == 1){
		rawcol = list(R = colnamesf[,1], G = colnamesf[,1])
	}
	
	if(df[1] == 0)
		stop(sprintf("None of the columns names of the expression files is matching a valid known quantitation type.\nTry to set the argument 'rawcol' by choosing among the following columns names: \n"), sprintf("\"%s\" \n",scanname))
	if(df[1] > 2)
		stop(sprintf("There are too many columns that could be read in the files.\n Try to set the argument 'rawcol' by choosing among the following columns names: \n"),sprintf("\"%s\" \n",scanname))
	
	return(rawcol)		
}

## assign phenoData to Nchannelset
assign.pheno.ncs = function(ph,files,raweset){
	si = pData(ph)[1:(length(files)*2),]
	lab = split(si,si[,"Label"])
	arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))
	
	#Reorder rows in each group (Cy3,Cy5) to the same order
#	seq<-with(lab[[1]],order(ArrayData.File))
	lab[[1]] = lab[[1]][order(lab[[1]][,arrayDataCol]),]
	lab[[2]] = lab[[2]][order(lab[[2]][,arrayDataCol]),]
	
	
	
	same = which(lapply(1:ncol(lab[[1]]), function(i) all(lab[[1]][i] == lab[[2]][i])) == TRUE)
	all = lab[[1]][same]
	gspe = lab[[1]][-same]
	colnames(gspe) = paste(colnames(gspe),names(lab)[1],sep = ".")
	rspe = lab[[2]][-same]
	colnames(rspe) = paste(colnames(rspe),names(lab)[2],sep = ".")
	
	metaData = data.frame(labelDescription = c(rep("_ALL_",ncol(all)),rep("G",ncol(gspe)),rep("R",ncol(rspe))))
	samples = new("AnnotatedDataFrame", data = cbind(all,gspe,rspe), varMetadata = metaData)
	
	arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(samples))	
	rownames(pData(samples)) = gsub(".[a-z][a-z][a-z]$","",samples[[arrayDataCol]])
	
	if(nrow(samples) != length(sampleNames(raweset)))
		warning("Number of data files read does not match available sample annotations")
	
	pData(samples) = pData(samples)[sampleNames(raweset),]
	phenoData(raweset) = samples
	return(raweset)
}

prepPhenoDataFor2channel = function(ph,files,raweset){
	si = pData(ph)[1:(length(files)*2),] # This is a cheat
	lab = split(si,si[,"Label"])
	same = which(lapply(1:ncol(lab[[1]]), function(i) all(lab[[1]][i] == lab[[2]][i])) == TRUE)
	all = lab[[1]][same]
	gspe = lab[[1]][-same]
	colnames(gspe) = paste(colnames(gspe),names(lab)[1],sep = ".")
	rspe = lab[[2]][-same]
	colnames(rspe) = paste(colnames(rspe),names(lab)[2],sep = ".")
	
	metaData = data.frame(labelDescription = c(rep("_ALL_",ncol(all)),rep("G",ncol(gspe)),rep("R",ncol(rspe))))
	
	samples = new("AnnotatedDataFrame", data = cbind(all,gspe,rspe), varMetadata = metaData)
	
	rownames(pData(samples)) = gsub(".[a-z][a-z][a-z]$","",samples$Array.Data)
	pData(samples) = pData(samples)[sampleNames(raweset),]
	return(samples)
}

getSDRFcolumn = function(col,headers){
	pattern<-switch(col,
			ArrayDataFile = "^Array[[:punct:]|[:blank:]]*Data[[:punct:]|[:blank:]]*File",
			ArrayDesignREF = "^Array[[:punct:]|[:blank:]]*Design[[:punct:]|[:blank:]]*REF",
			label = "^Label$")
	colIndex = grep(pattern,headers,ignore.case = TRUE)
	return(colIndex)
}


