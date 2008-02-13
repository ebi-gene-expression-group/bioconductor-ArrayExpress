extract.zip.from.AE = function (input, unzip = getOption("unzip"), tempoutdir = tempdir()) {
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
    return(sprintf("The %s data have been extracted in %s", input, tempoutdir))
    
}
