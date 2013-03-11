# TODO: Add comment
# 
# Author: iemam
###############################################################################


getAE = function (accession, path = getwd(), type = "full", extract = TRUE, local = FALSE, sourcedir=path) {
	
	if(!local){
		baseURL = "http://www.ebi.ac.uk/arrayexpress/xml/v2/files"
		xmlURL = paste(baseURL,accession,sep="/")
		xml = xmlTreeParse(xmlURL,useInternalNodes=TRUE)
		
		sdrfURL = xpathSApply(xml,"/files/experiment/file[kind='sdrf' and extension='txt']/url", xmlValue)
		sdrfFile = xpathSApply(xml,"/files/experiment/file[kind='sdrf' and extension='txt']/name", xmlValue)
		
		idfURL = xpathSApply(xml,"/files/experiment/file[kind='idf' and extension='txt']/url", xmlValue)
		idfFile = xpathSApply(xml,"/files/experiment/file[kind='idf' and extension='txt']/name", xmlValue)
		
		adfURL = xpathApply(xml,"/files/experiment/file[kind='adf' and extension='txt']/url", xmlValue)
		adfFiles = xpathApply(xml,"/files/experiment/file[kind='adf' and extension='txt']/name", xmlValue)
		
		rawArchiveURL = xpathApply(xml,"/files/experiment/file[kind='raw' and extension='zip']/url", xmlValue)
		procArchiveURL = xpathApply(xml,"/files/experiment/file[kind='processed' and extension='zip']/url", xmlValue)
		
	}else{
		
		allfiles = list.files(sourcedir)
		
		#SDRF
		sdrfFile = allfiles[grep(paste(accession,".sdrf.txt$",sep=""),allfiles)]
		if(length(sdrfFile)==0)
			stop("SDRF file not found in directory ",sourcedir)
		sdrfURL=paste("file:/",sourcedir,sdrfFile,sep="/")
		
		#IDF
		idfFile = allfiles[grep(paste(accession,".idf.txt$",sep=""),allfiles)]
		if(length(idfFile)==0)
			warning("IDF file not found in directory ",sourcedir)
		idfURL=paste("file:/",sourcedir,idfFile,sep="/")
		
		#ADF
		ph = try(read.AnnotatedDataFrame(sdrfFile, path = sourcedir, row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))
		if(inherits(ph,'try-error')){
			warning("Unable to retrieve ADF reference from SDRF. Reading any ADF in directory.")
			adfFiles = allfiles[grep(".adf.txt$",allfiles)]
		}
		else{
			adr = unique(pData(ph)[,getSDRFcolumn("ArrayDesignREF",varLabels(ph))])
			adfFiles = paste(adr,".adf.txt",sep="");
		}
			
		if(all(file.exists(file.path(sourcedir,adfFiles)))){
			adfURL = paste("file:/",sourcedir,adfFiles,sep="/")
			downloadADF = FALSE
		}
			
		else{
			#ADF not found in local directory. Attempt to retrieve it from FTP
			filesURL="http://www.ebi.ac.uk/arrayexpress/files";
			adfURL = paste(filesURL,adr,adfFiles,sep="/");
			downloadADF = TRUE
		}
			
			
		
		rawArchiveURL = NULL
		procArchiveURL = NULL
		#RAW files
		rawArchive = allfiles[grep(paste(accession,".raw.[0-9]{1,}.zip",sep=""),allfiles)]
		if(length(rawArchive)!=0)
			rawArchiveURL = paste("file:/",sourcedir,rawArchive,sep="/")
		else
			warning("No raw files found in directory ",sourcedir)
		
		
		#Processed files
		processedArchive = allfiles[grep(paste(accession,".processed.[0-9]{1,}.zip",sep=""),allfiles)]
		if(length(processedArchive)!=0)
			procArchiveURL = paste("file:/",sourcedir,processedArchive,sep="/")
		else
			warning("No processed data files found in directory ", sourcedir)
	}
	
	#a temporary solution for old GEO imports with seq and array files
	if(length(sdrfURL) > 1){
		warning("Found two SDRF files: \n",paste(sdrfURL,"\n"))
		hybSDRF = grep("hyb.sdrf",sdrfURL)
		if(length(hybSDRF)>0){
			message("Choosing ",sdrfURL[hybSDRF])
			sdrfURL=sdrfURL[hybSDRF];
			sdrfFile=sdrfFile[hybSDRF];
		}
		else{
			warning("Unable to choose SDRF file. Please report experiment to miamexpress@ebi.ac.uk")
		}
	}
	
	
	#Download/copy files to working directory specified by 'path'
	if(!local || path!=sourcedir || downloadADF){
		#Download/copy ADF
		adfFiles<-lapply(adfURL, function(url){
					filedest = paste(path,basename(url),sep="/")
					dnld = try(download.file(url, filedest, mode="wb"))
					if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
						warning(paste(url, " does not exist or is empty. \n"),sep="")
						adffile = NULL} 
					else{	
						adffile = basename(filedest)}
					return(adffile);
				})
		if(!is.null(adfFiles))
			adfFiles = unlist(adfFiles)
	}
		
	if(!local || path!=sourcedir){
		#Download/copy SDRF
		sdrfFileDest = paste(path,sdrfFile,sep="/")
		dnld = try(download.file(sdrfURL, sdrfFileDest, mode="wb"))
		if(inherits(dnld, 'try-error') || file.info(sdrfFileDest)$size == 0) {
			warning(paste(sdrfFile, " does not exist or is empty. The object will not have featureData or phenoData. \n"),sep="")
			sdrfFile = NULL
			adffile = NULL 
		}
		
		#Download/copy IDF
		idfFileDest = paste(path,idfFile,sep="/")
		dnld = try(download.file(idfURL, idfFileDest, mode="wb"))
		if(inherits(dnld, 'try-error') || file.info(idfFileDest)$size == 0) {
			warning(paste(idfFile, " does not exist or is empty. \n"),sep="")
			idfFile = NULL 
		}
		
		
		## Download data files
		rawArchive = NULL
		processedArchive = NULL
		
		##RAW DATA
		if(type!="mageFilesOnly" && !is.null(rawArchiveURL) && (type == "full" || type == "raw")){
			message("Copying raw data files\n")
			rawArchive<-lapply(rawArchiveURL, function(url){
						filedest = paste(path,basename(url),sep="/")
						dnld = try(download.file(url, filedest, mode="wb"))
						
						## Download raw.x.zip checking
						if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
							warning(paste(url, " does not exist or is empty. \n"),sep="")
						} else  {
							return (filedest)
						}
						
					})
			if(!is.null(rawArchive)){
				rawArchive = unlist(rawArchive)
				rawArchive = basename(rawArchive)
			}
			
		}
		
		##PROCESSED DATA
		if((type!="mageFilesOnly" && type == "full" || type == "processed") && !is.null(procArchiveURL)){
			message("Copying processed data files\n")
			processedArchive<-lapply(procArchiveURL, function(url){
						filedest = paste(path,basename(url),sep="/")
						dnld = try(download.file(url, filedest, mode="wb"))
						
						## Download processed.x.zip checking
						if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
							warning(paste(url, " does not exist or is empty. \n"),sep="")
						} else  {
							return(filedest)
						}
						
					})
			if(!is.null(processedArchive)){
				processedArchive = unlist(processedArchive)
				processedArchive = basename(processedArchive)
			}
		}	
	}
		
	
	rawFiles = NULL
	processedFiles = NULL
	
	if(extract){
		message("Unpacking data files")
		if(!is.null(rawArchive))
			rawFiles<-lapply(rawArchive, function(zipfile){
						rawfiles = extract.zip(file = paste(path, zipfile, sep="/"))
						return(rawfiles)
					})
		if(!is.null(processedArchive))
			processedFiles<-lapply(processedArchive, function(zipfile){
					procfiles = extract.zip(file = paste(path, zipfile, sep="/"))
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
			sdrf = sdrfFile, 
			idf = idfFile, 
			adf = adfFiles)
	return(res)
}
