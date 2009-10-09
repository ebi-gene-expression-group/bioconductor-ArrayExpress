getAE = function (input, path = getwd(), type = "full", extract = TRUE)
{
  oldopt = options()$HTTPUserAgent
  on.exit(options(HTTPUserAgent = oldopt))
  
  ## Building the link with the input name
  url = "http://www.ebi.ac.uk/microarray-as/ae/files"
  exp = paste(url,input,input,sep="/")
  indexf = file.path(path,paste("index",input,".html",sep=""))
  d = try(download.file(paste(url, input, "index.html",sep="/"), indexf , mode="wb"))
  ind = read.table(indexf, sep="\n", colClasses="character")
  
  if(type == "full" || type == "raw")
    {
      nraw = length(gregexpr("raw.[0-9]{1,}.zip",ind[1,], fixed=FALSE)[[1]])/2

      ##RAW DATA######################
      ## Saving temporarily the raw data
      options(HTTPUserAgent = paste(oldopt, "ArrayExpress", sessionInfo("ArrayExpress")$otherPkgs[[1]]$Version))

      rawdata = lapply(1:nraw, function(i) file.path(path,paste(input,".raw.",i,".zip",sep="")))

      dnld =  try(sapply(1:nraw, function(i) download.file(paste(exp,".raw.",i,".zip",sep=""), rawdata[[i]], mode="wb")))

      rawdata = unlist(rawdata)
       
      ## Download raw.zip checking
      if(inherits(dnld, 'try-error') || file.info(rawdata[[1]])$size == 0) {
        warning(paste(rawdata, " does not exist or is empty. \n"),sep="")
        rawdata = NULL
        rawfiles = NULL
	    rawcol = NULL
      } else  {
	if(extract==TRUE)
	{
	 rawfiles = extract.zip(file = rawdata)
	 rawcol = getcolraw(path, rawfiles)
	 }
	else  {
          rawfiles = rawdata
	  rawcol = NULL
	  }
        rawdata = basename(rawdata)
	}
    }
  
  if(type == "full" || type == "processed")
    {
      ##PROCESSED DATA######################
      ## Saving temporarily the processed data
      npro = length(gregexpr("processed.[0-9]{1,}.zip",ind[1,], fixed=FALSE)[[1]])/2

      procdata = lapply(1:npro, function(i) file.path(path,paste(input,".processed.",i,".zip",sep="")))

      dnldp = try(sapply(1:npro, function(i) download.file(paste(exp,".processed.",i,".zip",sep=""), procdata[[i]], mode="wb")))

      procdata = unlist(procdata)

      ## Download processed.zip checking
      if(inherits(dnldp, 'try-error') || file.info(procdata[[1]])$size == 0) {
        warning(paste(procdata, " does not exist or is empty. \n"),sep="")
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
  sdrffile = file.path(path,basename(samples))
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
          
          adffile = file.path(path, basename(featureannot))
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
  idffile = file.path(path,basename(annot))
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
      rawcol = NULL
    }

  res = list(path = path, rawdata = rawdata, rawfiles = rawfiles, rawcol = rawcol, procdata = procdata, procfile = procfile, sdrf = sdrffile, idf = idffile, adf = adffile)
  return(res)
  
}#end of getAE
