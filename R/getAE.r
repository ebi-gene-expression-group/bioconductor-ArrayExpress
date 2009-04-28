getAE = function (input, path = getwd(), type = "full", extract = TRUE)
{
  oldopt = options()$HTTPUserAgent
  on.exit(options(HTTPUserAgent = oldopt))
  
  ## Building the link with the input name
  url = "http://www.ebi.ac.uk/microarray-as/ae/files"
  exp = paste(url,input,input,sep="/")

  if(type == "full" || type == "raw")
    {
      ##RAW DATA######################
      ## Saving temporarily the raw data
      raw = paste(exp,".raw.zip",sep="")
      rawdata = file.path(path,paste(input,".raw.zip",sep=""))
      options(HTTPUserAgent = paste(oldopt, "ArrayExpress", sessionInfo("ArrayExpress")$otherPkgs[[1]]$Version))
      dnld = try(download.file(raw, rawdata, mode="wb"))
      
      ## Download raw.zip checking
      if(inherits(dnld, 'try-error') || file.info(rawdata)$size == 0) {
        warning(paste(raw, " does not exist or is empty. \n"),sep="")
        rawdata = NULL
        rawfiles = NULL
      } else  {
	if(extract==TRUE)
          rawfiles = extract.zip(file = rawdata)
	else
          rawfiles = rawdata
        rawdata = basename(rawdata) }
    }
  
  if(type == "full" || type == "processed")
    {
      ##PROCESSED DATA######################
      ## Saving temporarily the processed data
      proc = paste(exp,".processed.zip",sep="")
      procdata = file.path(path,paste(input,".processed.zip",sep=""))
      dnldp = try(download.file(proc, procdata, mode="wb"))

      ## Download processed.zip checking
      if(inherits(dnldp, 'try-error') || file.info(procdata)$size == 0) {
        warning(paste(proc, " does not exist or is empty. \n"),sep="")
        procdata = NULL
        procfile = NULL
      } else  {
        procfile = extract.zip(file = procdata)
        procdata = basename(procdata)
      }
    }

  ##SDRF DATA######################
  ## Download sdrf file
  samples = paste(exp,".sdrf.txt",sep="")
  sdrffile = paste(path,basename(samples),sep="/")
  sdrf = try(download.file(samples, sdrffile, mode="wb"))
  
  ## Download sdrf checking
  if(inherits(sdrf, 'try-error') || file.info(sdrffile)$size == 0) {
    warning(paste(sdrf, " does not exist or is empty. The object will not have featureData or phenoData. \n"),sep="")
    sdrffile = NULL
    adffile = NULL } else sdrffile = basename(sdrffile)

  ##ADF DATA######################
  ## Download adf file
  if(!inherits(sdrf, 'try-error') && file.info(file.path(path,sdrffile))$size != 0)
    {
      ph = try(read.AnnotatedDataFrame(sdrffile, path = path, row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))
      adr = try(unique(pData(ph)$Array.Design.REF))
      adr = adr[adr != ""]
      if(!inherits(adr, 'try-error'))
        {
          url3 = "http://www.ebi.ac.uk/microarray-as/ae/files"
          features = paste(adr,".adf.txt",sep="")
          featureannot = paste(url3,adr,features,sep="/")
          
          adffile = paste(path, basename(featureannot), sep="/")
          adf = try(lapply(1:length(featureannot), function(i) download.file(featureannot[[i]], adffile[i], mode="wb")))
        } else {
          adffile = NULL
          warning("Cannot retrieve the array design information from the sdrf file, the object will not have featureData attached. \n",sep="")
        }
      ## Download adf checking
      if(inherits(adf, 'try-error') || file.info(adffile)$size == 0) {
        warning(paste(adf, " does not exist or is empty. \n"),sep="")
        adffile = NULL } else adffile = basename(adffile)
    }
  
  
  ##IDF DATA######################
  ## Download idf file
  annot = paste(exp,".idf.txt",sep="")
  idffile = paste(path,basename(annot),sep="/")
  idf = try(download.file(annot, idffile, mode="wb"))

  ## Download idf checking
  if(inherits(idf, 'try-error') || file.info(idffile)$size == 0) {
    warning(paste(idf, " does not exist or is empty. \n"),sep="")
    idffile = NULL } else idffile = basename(idffile)
  
  ##EXPORTING RESULTS######################
  if(type == "raw")
    {
      procdata = NULL
      procfile = NULL
    }
  if(type == "processed")
    {
      rawdata = NULL
      rawfiles = NULL
    }

  res = list(path = path, rawdata = rawdata, rawfiles = rawfiles, procdata = procdata, procfile = procfile, sdrf = sdrffile, idf = idffile, adf = adffile)
  return(res)
  
}#end of getAE
