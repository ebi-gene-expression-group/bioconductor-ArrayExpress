# TODO: Add comment
# 
# Author: iemam
###############################################################################


ae2bioc = function(mageFiles, dataCols, save = TRUE){
	
#	if(!save) 
#		on.exit(cleanupAE(as.list(mageFiles)))
	
	dataFiles = mageFiles$rawFiles
	sdrf = mageFiles$sdrf
	idf = mageFiles$idf
	adf = mageFiles$adf
	path = mageFiles$path
	
	notuse = grep(dataFiles, pattern = "info.txt$|idf.txt$|processed|sdrf.txt$|.log$|RData|class")
	if (length(notuse) != 0) 
		dataFiles = dataFiles[-notuse]
	if (length(notuse) == 0) 
		dataFiles = dataFiles
	dataFiles = dataFiles[dataFiles != ""]
	
	#read sample annotations
	message("ArrayExpress: Reading pheno data from SDRF")
	ph = readPhenoData(sdrf,path);
	
	if(is.null(ph)){
		stop("Parsing SDRF failed. Please make sure SDRF file ",sdrf," exists in ",path, " and is not corrupt.")
	}
	
	#adf ref in SDRF
	adr = unique(pData(ph)[,getSDRFcolumn("ArrayDesignREF",varLabels(ph))])
	adr = adr[adr != ""]
	if(length(adr)>1)
		message("ArrayExpress: Experiment uses multiple Array Designs. A separate expressionSet will be created for each")
	
	if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",dataFiles, ignore.case = TRUE)) != 0)
		warning("Cannot find the array design reference in the sdrf file. The object may not be built.")
	if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",dataFiles, ignore.case = TRUE)) == 0)
		stop("Cannot find the array design reference in the sdrf file. The object cannot be built.")
	
	for (ad in adr){
		#Subselect SDRF data for current ArrayDesign REF
		if(length(adr)>1){
			res=getPhenoDataPerAD(ad,ph,dataFiles)
			dataFiles = res["files"]
			ph = res["ph"]	
		}	
		
		#read data files
		message("ArrayExpress: Reading data files")
		rawdata= try(readAEdata(path = path,files = dataFiles,dataCols=dataCols))
		if(inherits(rawdata, "try-error"))
			stop("Cannot read assay Data")
		
		#read and match array feature metadata to raw data
		adfFile = adf[grep(ad,adf)]
		features= try(readFeatures(rawdata=rawdata,adf=adfFile,path=path))
		if(inherits(features, "try-error"))
			warning("Cannot read feature Data")
		
		#read experiment meta data
		experimentData = try(readExperimentData(idf=idf,path=path))
		if(inherits(experimentData, "try-error"))
			warning("Cannot read experiment Data")
		
		
		
		#Finally build ExpressionSet
		
		if(class(rawdata) == "AffyBatch"){
			raweset=rawdata
			experimentData(raweset) = experimentData
			phenoData(raweset) = ph
		}
		
		if(class(rawdata) == "RGList"){
			#construct nchannelset
			assayData = if("Rb" %in% names(rawdata))
						with(rawdata, assayDataNew(R = R, G = G, Rb = Rb, Gb = Gb)) 
					else 
						with(rawdata, assayDataNew(R = R, G = G))
			
			raweset = new("NChannelSet",
					assayData = assayData,
					featureData = features,
					experimentData = experimentData)
			
			raweset = try(assign.pheno.ncs(files = dataFiles,ph = ph,raweset = raweset))
		}
		
		if(class(rawdata) == "EListRaw"){
			#construct expressionSet
			assayData = rawdata$G
			raweset = new("ExpressionSet",
					expr = rawdata$G,
					phenoData = ph,
					featureData = features,
					experimentData = experimentData)
		}
		
		
		
#		if(class(rawdata) == "RGList"){
#			assayData = if("Rb" %in% names(rawdata))
#							with(rawdata, assayDataNew(R = R, G = G, Rb = Rb, Gb = Gb)) 
#						else 
#							with(rawdata, assayDataNew(R = R, G = G))}
#		
#		if(class(rawdata) == "EListRaw"){
#			raweset = try(new("ExpressionSet",exprs = rawdata$G))
#			if(inherits(raweset, "try-error"))
#				stop(sprintf(raweset[1]))}
		
		
		#attach phenoData		
#		if(class(rawset)=="NChannelSet"){
#			rawesetph = try(assign.pheno.ncs(files = dataFiles,ph = ph,raweset = rawset))
#			if(!inherits(rawesetph, "try-error"))
#				rawset = rawesetph
#		}else{
#			##pData(ph) = pData(ph)[sampleNames(raweset),]
#			#check number of assays in rawset is equal to those in ph
#			phenoData(rawset) = ph
#		}
		
		#attach experiment meta data
#		rawesetex = try(creating_experiment(idf = idf, eset = rawset, path = path))
#		if(!inherits(rawesetex, "try-error"))
#			rawset = rawesetex
#		else warning("Cannot attach experimentData")
		
		#attach array meta data
#		adfFile = adf[grep(ad,adf)]
#		rawesetex = try(addADF(adf = adfFile, eset = rawset, path = path, files=dataFiles))
#		if(!inherits(rawesetex, "try-error"))
#			rawset = rawesetex 
#		else warning("Cannot attach featureData")
		
	}
	
	return(raweset)
}
