getAE = function (input, path = ".", save = TRUE, type = "full", extract = TRUE) {
  if(!save) on.exit({try(file.remove(file.path(path, rawdata)));try(file.remove(file.path(path, procdata)))})
  
  ## Building the link with the input name
  dir = gsub("^E-|-[0-9]{1,10}","",input)
  if(nchar(dir) == 5)
    dir = gsub("[a-z]$","",dir)
  url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment"
  exp = paste(url,dir,input,input,sep="/")

  if(type == "full" || type == "raw")
    {
      ##RAW DATA######################
      ## Saving temporarily the raw data
      raw = paste(exp,".raw.zip",sep="")
      rawdata = file.path(path,paste(input,".raw.zip",sep=""))
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
        procfiles = NULL
      } else  {
        procfiles = extract.zip(file = procdata)
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
          dira = gsub("^A-|-[0-9]{1,10}","",adr)
          url3 = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/array"
          features = paste(adr,".adf.txt",sep="")
          featureannot = paste(url3,dira,adr,features,sep="/")
          
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
      procfiles = NULL
    }
  if(type == "processed")
    {
      rawdata = NULL
      rawfiles = NULL
    }
  res = list(path = path, rawdata = rawdata, rawfiles = rawfiles, procdata = procdata, procfiles = procfiles, sdrf = sdrffile, idf = idffile, adf = adffile)
  return(res)
  
}#end of getAE
