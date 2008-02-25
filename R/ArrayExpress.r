ArrayExpress = function(input, tempoutdir = ".")
  {
    dir = gsub("^E-|-[0-9]{1,10}","",input)
    url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment"
    exp = paste(url,dir,input,input,sep="/")

    ##saving temporarily the raw data
    raw = paste(exp,".raw.zip",sep="")
    if(!exists(tempoutdir))
     tempoutdir =  tempdir()
    rawdata = tempfile(tmpdir = tempoutdir)
    dnld = try(download.file(raw, rawdata))

    if(class(dnld) == 'try-error')
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
        stop(sprintf("'%s' is not an Affymetrix experiment.",input))
      }

    if(length(grep(".cel",files)) == length(files))
      raweset = try(ReadAffy(filenames = paste(tempoutdir,files,sep="/")))
    
    if(class(raweset) == 'try-error')
      {
        file.remove(allfiles)
        file.remove(rawdata)
        stop(sprintf(raweset[1]))
      }

    file.remove(allfiles)
    file.remove(rawdata)
    return(raweset)

  }
