# TODO: Add comment
# 
# Author: iemam
###############################################################################


getAE = function (input, path = getwd(), type = "raw", extract = TRUE) {
	baseURL = "http://www.ebi.ac.uk/arrayexpress/xml/v2/files"
	#input="E-TABM-25"
	xmlURL = paste(baseURL,input,sep="/")
	xml = xmlTreeParse(xmlURL,useInternalNodes=TRUE)
	
	#node = getNodeSet(xml,"/files/experiment/file[kind='sdrf' and extension='txt']/url")
	#sdrfURL = xmlSApply(node[[1]], xmlValue)
	#sdrfName = getNodeSet(xml,"/files/experiment/file[kind='sdrf' and extension='txt']/name")
	
	sdrfURL = xpathSApply(xml,"/files/experiment/file[kind='sdrf' and extension='txt']/url", xmlValue)
	sdrfName = xpathSApply(xml,"/files/experiment/file[kind='sdrf' and extension='txt']/name", xmlValue)
	
	idfURL = xpathSApply(xml,"/files/experiment/file[kind='idf' and extension='txt']/url", xmlValue)
	idfName = xpathSApply(xml,"/files/experiment/file[kind='idf' and extension='txt']/name", xmlValue)
	
	adfURL_list = xpathApply(xml,"/files/experiment/file[kind='adf' and extension='txt']/url", xmlValue)
	adfName_list = xpathApply(xml,"/files/experiment/file[kind='adf' and extension='txt']/name", xmlValue)
	
	rawArchiveURL = xpathApply(xml,"/files/experiment/file[kind='raw' and extension='zip']/url", xmlValue)
	procArchiveURL = xpathApply(xml,"/files/experiment/file[kind='fgem' and extension='zip']/url", xmlValue)
	
	#Download files
	sdrffile = paste(path,sdrfName,sep="/")
	sdrf = try(download.file(sdrfURL, sdrffile, mode="wb"))
	
	## Download sdrf checking
	if(inherits(sdrf, 'try-error') || file.info(sdrffile)$size == 0) {
		warning(paste(sdrf, " does not exist or is empty. The object will not have featureData or phenoData. \n"),sep="")
		sdrffile = NULL
		adffile = NULL 
	} else sdrffile = basename(sdrffile)
	
	idffile = paste(path,idfName,sep="/")
	idf = try(download.file(idfURL, idffile, mode="wb"))
	
	## Download idf checking
	if(inherits(idf, 'try-error') || file.info(idffile)$size == 0) {
		warning(paste(idf, " does not exist or is empty. \n"),sep="")
		idffile = NULL 
	} else idffile = basename(idffile)
	
	
	adffiles<-lapply(adfURL_list, function(adfURL){
				filedest = paste(path,basename(adfURL),sep="/")
				dnld = try(download.file(adfURL, filedest, mode="wb"))
				
				## Download adf checking
				if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
					warning(paste(adfURL, " does not exist or is empty. \n"),sep="")
					adffile = NULL } 
				else {	
					adffile = basename(filedest)
				}
				return(adffile);
			})
	
	if(!is.null(adffiles))
		adffiles = unlist(adffiles)
	
	## Download data files
	rawArchive = NULL
	processedArchive = NULL
	
	if(!is.null(rawArchiveURL) && (type == "full" || type == "raw")){
		##RAW DATA######################
		## Saving temporarily the raw data
		message("Downloading raw data files\n")
		rawArchive<-lapply(rawArchiveURL, function(url){
					filedest = paste(path,basename(url),sep="/")
					dnld = try(download.file(url, filedest, mode="wb"))
					
					## Download raw.x.zip checking
					if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
						warning(paste(url, " does not exist or is empty. \n"),sep="")
						#rawfiles = NULL
					} else  {
						return (filedest)
						#rawfiles = rbind(rawdata)
						#uz = try(rawfiles = unzip(filedest)) 
					}
					
				})
		if(!is.null(rawArchive)){
			rawArchive = unlist(rawArchive)
			rawArchive = basename(rawArchive)
		}
		
	}
	
	if((type == "full" || type == "processed") && !is.null(procArchiveURL)){
		##PROCESSED DATA######################
		## Saving temporarily the processed data
		message("Downloading processed data files\n")
		processedArchive<-lapply(procArchiveURL, function(url){
					filedest = paste(path,basename(url),sep="/")
					dnld = try(download.file(url, filedest, mode="wb"))
					
					## Download processed.x.zip checking
					if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
						warning(paste(url, " does not exist or is empty. \n"),sep="")
						#procfiles = NULL
						#return(NULL)
					} else  {
						return(filedest)
#						cat("Extracting zip file\n")
#						procdata = extract.zip(file = filedest)
						#procfiles = rbind(procdata)
						#uz = try(rawfiles = unzip(filedest)) 
					}
					
				})
		if(!is.null(processedArchive)){
			processedArchive = unlist(processedArchive)
			processedArchive = basename(processedArchive)
		}
		
	}
	
	rawFiles = NULL
	processedFiles = NULL
	
	if(extract){
		message("Unzipping files")
		if(!is.null(rawArchive))
			rawFiles<-lapply(rawArchive, function(zipfile){
						rawfiles = extract.zip(file = zipfile)
						return(rawfiles)
					})
		if(!is.null(processedArchive))
			processedFiles<-lapply(processedArchive, function(zipfile){
						procfiles = extract.zip(file = zipfile)
						return(procfiles)
					})
		
		if(!is.null(rawFiles))
			rawFiles = unlist(rawFiles)
		if(!is.null(processedFiles))
			processedFiles = unlist(processedFiles)
	}
	
	res = list(path = path,  
			rawFiles = rawFiles,
			rawArchive = rawArchive,
			processedFiles = processedFiles,
			processedArchive = processedArchive,
			sdrf = sdrffile, 
			idf = idffile, 
			adf = adffiles)
	return(res)
}

