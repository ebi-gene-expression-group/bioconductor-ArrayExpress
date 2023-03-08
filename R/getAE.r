get_request <-function(url){
  qr = URLencode(url)
  response <- GET( qr )
  if(status_code(response) != 200) {
    stop(
      paste(
        "Error running query. Received HTTP error code",
        status_code(response),
        "from server. Please try again later. If you continue to experience 
        problems please contact us 
        at https://www.ebi.ac.uk/about/contact/support/arrayexpress"
      )
    )
  }
  
  json_parsed <- fromJSON(txt = qr) 
  return(json_parsed)
}

download_file <- function(url, filedest, overwrite){
  if (file.exists(filedest) && overwrite == FALSE){
    message(paste(filedest, ' already exists, skipping.', sep =""))
    retVal <- filedest
  }
  else{
    dnld = try(download.file(url, filedest))
    if(inherits(dnld, 'try-error') || file.info(filedest)$size == 0) {
      retVal <- NULL
    } else  {
      retVal <- filedest
    }    
  }
  return(retVal)
}


get_filenames <-function(jsonData, ftpUrl, type){
  subsections = jsonData$section$subsections
  fileNames = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(fileNames) <- c("type", "file", "url")
  for (s in subsections) {
    if (s$type[1] == "Assays and Data"){
      data_files = s$subsections
      break
    }
  }
  
  for (row in 1:nrow(data_files)){
    file_type = data_files[row, "type"]
    files=data_files[row, "files"]
    for (i in files){
      if (is.data.frame(i)){
        df=i
      }
      else{
        for (n in i){
          df=n
        }
      }    
    }
    
    for(f in df$path){
      if (file_type %in% c("Processed Data","MAGE-TAB Files", "Raw Data")) {
        fileNames[nrow(fileNames) + 1,] = c(file_type, f, paste(ftpUrl,"Files",f,sep="/"))
      }
    }
    if (file_type == "Array Designs"){
      for (i in data_files[row, "links"]){
        for (n in i){
          adf_accession = n$url
          if(!is_empty(adf_accession)){
            adf_name = paste(adf_accession, "adf","txt", sep=".")
            adf_link = paste("https://www.ebi.ac.uk/biostudies/files", adf_accession, adf_name, sep="/")
            fileNames[nrow(fileNames) + 1,] = c(file_type, adf_name, adf_link)
          }
        }
      }
    }
  }
  
  if (tolower(type)=='raw' || tolower(type) == 'full'){
    enaAccession <- NULL
    jsonLinks = unlist(jsonData$section$links)
    if(!is_empty(jsonLinks)){
      if (toupper(jsonLinks['attributes.value']) == 'ENA'){
        enaAccession = jsonLinks['url']
      }     
    }
    if (!is_empty(enaAccession)){
      enaUrl=paste("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=",enaAccession,"&download=false&fields=fastq_ftp&format=json&limit=0&result=read_run", sep="")
      enaData= get_request(enaUrl)
      for(i in 1:nrow(enaData)) { 
        ftpLikns = strsplit(enaData[i, "fastq_ftp"], split = ";")
        for (enaFileLinkList in ftpLikns){
          for (enaFileLink in enaFileLinkList){
            enaFileName = sub(".*/", "", enaFileLink)
            fileNames[nrow(fileNames) + 1,] = c("Raw Data", enaFileName, paste("ftp:", enaFileLink, sep="//"))
          }
        }
      }  
    }
    
    if (tolower(type)=='raw'){
      fileNames= fileNames[fileNames$type == "Raw Data", ]
    }
  }
  else if (tolower(type)=='processed'){
    fileNames= fileNames[fileNames$type == "Processed Data", ]
  }
  else if (tolower(type)=='mage'){
    fileNames= fileNames[fileNames$type == "MAGE-TAB Files", ]
  }
  
  
  return(fileNames)
}

getAE = function (accession, path = getwd(), type = "full", extract = TRUE, sourcedir=path, overwrite = FALSE) {
  
  typeOptions <- c("raw", "full", "processed", "mage")
  if (!tolower(type) %in% typeOptions){
    print(paste("Option'",type,"is not a valid option for parameter 'type'", sep=" "))
    return()
  }
  
  baseURL = "https://www.ebi.ac.uk/biostudies/api/v1/studies"
  httpURL = paste(baseURL, accession, sep="/")
  submissionMetaData = get_request(httpURL)
  
  infoURL = paste(baseURL, accession, "info", sep="/")
  infoData = get_request(infoURL)
  ftpLink = infoData$ftpLink
  allFileNames = get_filenames(submissionMetaData, ftpLink, type)
  
  #Download files
  rawArchive <- NULL
  processedArchive <- NULL
  mageTabFiles <- NULL
  adfFiles <- NULL
  
  for (row_file in 1:nrow(allFileNames)){
    
    filedest = paste(path,allFileNames[row_file,"file"],sep="/")
    url = allFileNames[row_file,"url"]
    message(paste("Downloading file: ", allFileNames[row_file,"file"], "\n", sep=""))
    file_ok <- download_file(url, filedest, overwrite)
    
    if (is.null(file_ok)){
      warning(paste(url, " does not exist or is empty. \n"),sep="")
    }
    else{
      if (allFileNames[row_file,"type"] == 'Raw Data'){
        rawArchive <- c(rawArchive, file_ok)
      }
      else if(allFileNames[row_file,"type"] == 'MAGE-TAB Files'){
        
        mageTabFiles <- c(mageTabFiles, file_ok)         
      }
      else if(allFileNames[row_file,"type"] == 'Processed Data'){
        processedArchive <- c(processedArchive, file_ok)
      }
      else if(allFileNames[row_file,"type"] == 'Array Designs'){
        adfFiles <- c(adfFiles, file_ok)          
      }
    }
  }
  
  rawFiles = NULL
  processedFiles = NULL
  
  if(extract){
    message("Unpacking data files")
    if(!is.null(rawArchive))
      rawFiles<-lapply(rawArchive, function(zipfile){
        rFiles = zipfile
        if (file_ext(zipfile) == "zip"){
          rFiles = extract.zip(file = zipfile)
        }
        return(rFiles)
      })
    if(!is.null(processedArchive))
      processedFiles<-lapply(processedArchive, function(zipfile){
        procFiles = zipfile
        if (file_ext(zipfile) == "zip"){
          procFiles = extract.zip(file = zipfile)
        }
        return(procFiles)
      })
    
    if(!is.null(rawFiles))
      rawFiles = unlist(rawFiles)
    if(!is.null(processedFiles))
      processedFiles = unlist(processedFiles)
  }
  
  sdrfFile <- NULL
  idfFile <- NULL
  
  if (!is_empty(mageTabFiles)){
    sdrfFile = mageTabFiles[grep(mageTabFiles, pattern = "sdrf.txt$")]
    idfFile = mageTabFiles[grep(mageTabFiles, pattern = "idf.txt$")]
  }
  
  res = list(path = path, 
             rawFiles = rawFiles,
             rawArchive = rawArchive,
             processedFiles = processedFiles,
             processedArchive = processedArchive,
             mageTabFiles = mageTabFiles, 
             sdrf = sdrfFile,
             idf = idfFile,
             adf = adfFiles,
             dataFiles = allFileNames)
  return(res)
}