# A function that takes ArrayExpress MAGETAB files for a specific experiment and returns an equivalent R object representation (NChannelSet, ExpressionSet or AffyBatch) 
# 
# Author: iemam
###############################################################################


ae2bioc = function(mageFiles, dataCols=NULL, drop=TRUE){
	
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
	dataFiles = dataFiles[dataFiles != ""]
	allDataFiles = dataFiles;
	
	if(length(dataFiles)==0)
		stop("ArrayExpress: Experiment has no raw files available. Consider using processed data instead by following procedure in the vignette")
	
	#check for duplicates in sdrf
	sdrfData<-read.delim(paste(path,sdrf,sep="/"), check.names=FALSE)
	if(any(duplicated(sdrfData[, 'Array Data File']))){
	  message("Duplicates found in SDRF file")
	  #remove duplicates based on Array Data Files and update SDRF
	  sdrfData<-sdrfData[!duplicated(sdrfData[, 'Array Data File']), ]
	  write.table(sdrfData, file=sdrf, row.names=FALSE, quote=FALSE, sep='\t')
	  message("Removed duplicates in SDRF file")
	}

	#read sample annotations
	ph = readPhenoData(sdrf,path)
	if(inherits(ph, 'try-error')){
		ph=NULL
		stop("ArrayExpress: Parsing SDRF failed. Please make sure SDRF file ",sdrf," exists in ",path, " and is not corrupt.")
	}
	fullPhenoData = ph;
		
	#Checks
	arrayDataCol = getSDRFcolumn("ArrayDataFile",varLabels(ph))
#	labelCol = getSDRFcolumn("label",varLabels(ph))
	
	#ArrayDesign REF in SDRF
	adr = unique(pData(ph)[,getSDRFcolumn("ArrayDesignREF",varLabels(ph))])
	adr = adr[adr != ""]
	
	if(!all(dataFiles %in% ph[[arrayDataCol]]))
		warning("Some data files in the zip archive are missing from the SDRF. The object may not be built.")
	if(length(adr)>1)
		message("ArrayExpress: Experiment uses multiple Array Designs. A separate expressionSet will be created for each")	
	if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",dataFiles, ignore.case = TRUE)) != 0)
		warning("ArrayExpress: Cannot find the array design reference in the sdrf file. The object may not be built.")
	if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",dataFiles, ignore.case = TRUE)) == 0)
		stop("ArrayExpress: Cannot find the array design reference in the sdrf file. The object cannot be built.")
	
	#list of return R objects
	robjs=list();
	
	for (ad in adr){
		#Subselect SDRF data for current ArrayDesign REF
		if(length(adr)>1){
			res=getPhenoDataPerAD(ad,fullPhenoData,allDataFiles)
			dataFiles = unique(res$dataFiles)
			ph = res$pheno
		}	
		
		#read data files
		green.only = isOneChannel(sdrf,path)
		rawdata= try(readAEdata(path = path,files = dataFiles,dataCols=dataCols,green.only=green.only))
		if(inherits(rawdata, "try-error"))
			stop("ArrayExpress: Unable to read assay data")
		
		#read and match array feature metadata to raw data
		if(!inherits(rawdata,"FeatureSet")){
			adfFile = adf[grep(ad,adf)]
			features= try(readFeatures(adf=adfFile,path=path))
			if(inherits(features, "try-error")){
				warning("ArrayExpress: Unable to read feature data")
				features = NULL;
			}
		}
			
		#read experiment meta data
		experimentData = try(readExperimentData(idf=idf,path=path))
		if(inherits(experimentData, "try-error")){
			warning("ArrayExpress: Unable to read experiment data");
			experimentData = new("MIAME");
		}
			

		#Finally build ExpressionSet
		
		#Attach pheno and feature data to oligo::FeatureSet
		if(inherits(rawdata,"FeatureSet")){
			raweset=rawdata
			phenoData(raweset) = ph[sampleNames(rawdata)]
		}
		
		
		if(class(rawdata) == "RGList" | class(rawdata) == "EListRaw"){
			#construct nchannelset
			if(class(rawdata) == "RGList"){
				assayData = if("Rb" %in% names(rawdata)) #FIXME: keep all
								with(rawdata, assayDataNew(R = R, G = G, Rb = Rb, Gb = Gb)) #will not work if datacolumns where user specified
							else 
							with(rawdata, assayDataNew(G = G, R = R))
				
				raweset = new("NChannelSet",
								assayData = assayData,
								experimentData = experimentData)
			}
			
			#construct expressionSet
			if(class(rawdata) == "EListRaw"){
				assayData = with(rawdata, assayDataNew(E = E, Eb = Eb))
				raweset = new("NChannelSet",
								assayData = assayData,
								experimentData = experimentData)
			}
			
			#Attach pheno data
			if(!is.null(ph)){ #FIXME:?
				#imagene doesnt have targets slot
				ph = ph[rawdata$targets$FileName,]
				phenoData(raweset) = ph
			}
			
			#Attach features
			if(!is.null(rawdata$genes) && rawdata$source != "ae1"){
				features2 = new("AnnotatedDataFrame",rawdata$genes)
				featureData(raweset) = features2;
			}
			else
				if(!is.null(features))
					featureData(raweset) = features;
			#consistency of order between features in featureData and assayData is established via previous sorting of feature columns
			
			if(length(featureNames(assayData(raweset))) != length(featureNames(featureData(raweset))))
				warning("Number of features in assayData and featureData are not equal. Check control features (NA) that might have been removed from either assayData or featureData.");
			
		}
		
		
		robjs[[ad]]=raweset
			
	}

    if(drop && length(robjs) == 1)
        robjs = robjs[[1]]
	
	return(robjs)
}
