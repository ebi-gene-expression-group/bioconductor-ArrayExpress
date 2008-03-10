ArrayExpress = function(input, tempoutdir = ".")
  {
    dir = gsub("^E-|-[0-9]{1,10}","",input)
    url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment"
    exp = paste(url,dir,input,input,sep="/")
    raw = paste(exp,".raw.zip",sep="")
    ## Saving temporarily the raw data
    rawdata = tempfile(tmpdir = tempoutdir)
    dnld = try(download.file(raw, rawdata, mode="wb"))

    if(inherits(dnld, 'try-error'))
      {
        file.remove(rawdata)
        stop(paste(raw, " does not exist. \n"),sep="")
      }

    if(file.info(rawdata)$size == 0)
      {
        file.remove(rawdata)
        stop(paste(raw, " is empty. \n"),sep="")
      }

    extract.zip(file = rawdata)
   
    ## creating the expression set of raw data
    allfiles = dir(tempoutdir, pattern = input)
    notuse = grep(allfiles, pattern = "info.txt$|idf.txt$|processed|sdrf.txt$|.log$")
    if(length(notuse)!=0)
      files = allfiles[-notuse]
    if(length(notuse)==0)
      files = allfiles

    if(length(grep(".cel",files)) == 0)
      {
        file.remove(allfiles)
        file.remove(rawdata)
        stop(sprintf("'%s' non Affymetrix experiments are not handled yet.",input))
      }

    if(length(grep(".cel",files)) == length(files))
      {
        raweset = try(ReadAffy(filenames = paste(tempoutdir,files,sep="/")))
        samples = paste(exp,".sdrf.txt",sep="")
        ph = try(read.AnnotatedDataFrame(samples,row.names=NULL, blank.lines.skip = TRUE, fill=TRUE))       
       
       if(inherits(ph, 'try-error'))
          warning(sprintf("Cannot attach phenoData to the object as the sdrf file cannot be read or is empty."))            
          
        if(!inherits(ph, 'try-error'))
          {
            pData(ph) = pData(ph)[1:length(files),]
            if(length(ph$Array.Data.File)==length(unique(ph$Array.Data.File)))
              {
                rownames(pData(ph)) = ph$Array.Data.File
                pData(ph) = pData(ph)[sampleNames(raweset),]
                phenoData(raweset) = ph
              } else {
                warning(sprintf("Cannot attach phenoData to the object as the Array Data File column of the sdrf file for %s contains duplicated elements.", input)) 
              }
          }
      }
    
    if(inherits(raweset,'try-error'))
      {
        file.remove(allfiles)
        file.remove(rawdata)
        if(length(grep("dimensions",raweset[1])) != 0)
          stop(sprintf("%s \n This file may not be a binary cel file. \n", raweset[1]))
        if(length(grep("type",raweset[1])) != 0)
          stop(sprintf("%s \n ArrayExpress package does not handle multiple chip types data sets yet. \n", raweset[1]))
        if(length(grep("type",raweset[1])) == 0 && length(grep("dimensions",raweset[1])) == 0)         
        stop(sprintf(raweset[1]))
      }

    file.remove(allfiles)
    file.remove(rawdata)
    return(raweset)

  }
