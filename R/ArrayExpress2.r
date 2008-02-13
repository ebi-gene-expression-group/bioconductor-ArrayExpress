ArrayExpress2 = function(input, tempoutdir = tempdir())
  {

    ## Building the link with the input name
    dir = gsub("^E-|-[0-9]{1,10}","",input)
    url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment"
    exp = paste(url,dir,input,input,sep="/")
    raw = paste(exp,".raw.zip",sep="")
    ## Saving temporarily the raw data
    rawdata = tempfile(tmpdir = tempoutdir)
    dnld = try(download.file(raw, rawdata))

    ## Download checking
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
    
    ## Extracting the data if the checking was fine
    extract.zip(file = rawdata)
   
    ## Listing the expression files names
    allfiles = dir(tempoutdir, pattern = input)
    notuse = grep(allfiles, pattern = "info.txt$|idf.txt$|processed|sdrf.txt$|.log$")
    if(length(notuse) != 0)
      files = allfiles[-notuse]
    if(length(notuse) == 0)
      files = allfiles

    ## Building the S4 class object
    
    ## Create AffyBatch for Affymetrix data sets
    if(length(grep(".cel",files)) == length(files))
      raweset = try(ReadAffy(filenames = paste(tempoutdir,files,sep="/")))    
   
    ## Non Affymetrix data
    if(length(grep(".cel",files)) == 0)
      {
        ## Downloading and reading the sample information file .sdrf
        samples = paste(exp,".sdrf.txt",sep="")
        ph = try(read.AnnotatedDataFrame(samples,row.names=NULL, blank.lines.skip = T, fill=TRUE))
        ## Checking the number of dyes in the sdrf file
        if(!is.null(class(ph$Label)))
          dyes = levels(as.factor(ph$Label))           
        if(is.null(class(ph$Label)))
          stop("Cannot determine the number of dyes from the .sdrf file (usually in Label column).")
        
        ndyes = length(dyes[dyes!=""])
        
        ## Building ExpressionSet when one colour
        if(ndyes == 1)
          {
            columns = 0
            ## Parsing the first line of the expression file
            allcnames = scan(paste(tempoutdir,files[1],sep="/"),what="",nlines=1, sep="\t")
            
            ## Looking for the right column to use
            if(length(grep("raw_intensity",allcnames, ignore.case=TRUE)) != 0)  columns = "raw_intensity"         
            if(length(grep("raw intensity",allcnames, ignore.case=TRUE)) != 0)  columns = "raw intensity"
            if(length(grep("mean intensity",allcnames, ignore.case=TRUE)) != 0)  columns = "mean intensity"
            if(length(grep("signal mean",allcnames, ignore.case=TRUE)) != 0)  columns = "signal mean"
            if(length(grep("AVG_Signal",allcnames, ignore.case=TRUE)) != 0)  columns = "AVG_Signal"
            if(length(grep("gMeanSignal",allcnames, ignore.case=TRUE)) != 0)  columns = "gMeanSignal"
                        
            ## Check that the column name is of known format 
            if(columns == 0)
              {
                file.remove(allfiles)
                file.remove(rawdata)
                stop(sprintf("Data are not from 1 colour standard source type."))
              }
            
            ## Creating the S4 object
            raweset = try(read.AE_1col(path=tempoutdir,files=files, columns=columns))
          }#end ndyes ==1
                
        ## Building NChannelSet when two colour
        if(ndyes == 2)
          {
            idsource = 0
            ## Parsing the first line of the expression file
            allcnames = scan(paste(tempoutdir,files[1],sep="/"),what="",nlines=1, sep="\t")
            
            ## Looking for the right source to use
            if(length(grep("F635 Mean$",allcnames)) != 0)  idsource = "genepix"
            if(length(grep("rMeanSignal$",allcnames)) != 0)  idsource = "agilent"
            if(length(grep("MTM Dens - Levels$",allcnames)) != 0)  idsource = "arrayvision"
            if(length(grep("AMPCH1$",allcnames)) != 0)  idsource = "bluefuse"
            if(length(grep("ch2 Intensity$",allcnames)) != 0)  idsource = "quantarray"
            if(length(grep("Ch1 Mean$",allcnames)) != 0)  idsource = "scanarrayexpress"
            if(length(grep("CH1I_MEAN$",allcnames)) != 0)  idsource = "smd.old"
            if(length(grep("Ch1 Intensity (Mean)$",allcnames)) != 0)  idsource = "smd"
            if(length(grep("morphR.close.open$",allcnames)) != 0)  idsource = "spot.close.open"
            if(idsource != "spot" && length(grep("Rmean$",allcnames)) != 0)  idsource = "spot"
            
            ## Check that the source type is of known format 
            if(idsource == 0)
              {
                if(length(grep("Mean [QL] Cy3$",allcnames)) != 0)  idsource = "generic1"
                if(length(grep("Cy3_Volume$",allcnames)) != 0)  idsource = "generic2"
                if(length(grep("Signal med Cy3$",allcnames)) != 0)  idsource = "generic3"
                if(length(grep("Med.Sig.Cy3$",allcnames)) != 0)  idsource = "generic4"
                if(length(grep("Signal Mean_Cy3$",allcnames)) != 0)  idsource = "generic5"
                if(length(grep("Signal_Mean_Cy3$",allcnames)) != 0)  idsource = "generic6" ## E-MEXP-181 works
                if(length(grep("F633 Mean$",allcnames)) != 0)  idsource = "generic7" ## other channel is F543 Mean E-MEXP-406 works
                if(length(grep("Cy3 Intensity$",allcnames)) != 0)  idsource = "generic8" ## E-TABM-106 works
                if(length(grep("Cy3 Mean$",allcnames)) != 0)  idsource = "generic9" ## E-TABM-138 works
                if(length(grep("Cy3_Intensity$",allcnames)) != 0)  idsource = "generic10" ## E-TIGR-200 works
                #if(idsource != 0) stop(sprintf("New 2 colour standard source type."))
              }
            if(idsource == 0)
              {
                file.remove(allfiles)
                file.remove(rawdata)                
                stop(sprintf("Data are not from 2 colour standard source type."))
              }
              

            ## Creating the RGList
            eset_RGList = try(read.AE_2col(files=files,
              path = tempoutdir,
              source = idsource))
    
            ## Checking that the object has been successfully created
            if(class(eset_RGList) == 'try-error')
              {
                file.remove(allfiles)
                file.remove(rawdata)
                stop(sprintf("Error in read.AE_2col: '%s'.", eset_RGList[1]))

              }
            
            ## Creating the S4 object
            assayData = with(eset_RGList, assayDataNew(R=R, G=G, Rb=Rb, Gb=Gb))
            raweset = try(new("NChannelSet",
              assayData = assayData))
          }#end ndyes ==2


        if(ndyes != 1 && ndyes != 2)
          {
            file.remove(allfiles)
            file.remove(rawdata)
            stop(sprintf("Label column of the .sdrf file is not useable"))
          }
      }#end of non Affymetrix objects

    ## Removing the files downloaded and extracted
    file.remove(allfiles)
    file.remove(rawdata)
    
    ## Checking that the object has been successfully created
    if(class(raweset) == 'try-error')      
      stop(sprintf(raweset[1]))      

    return(raweset)
    
  }#end of ArrayExpress2


# Read 1 color files from AE to ExpressionSet
# Modificated version of read.maimages from limma

read.AE_1col <- function(files=NULL,path=NULL,columns=NULL,verbose=TRUE,sep="\t",quote=NULL,...)
{
  ## Begin checking input arguments
  slides <- as.vector(as.character(files))
  nslides <- length(slides)
  names <- removeExt(files)        
  ## End checking input arguments

  ## Read first file to get nspots
  fullname <- slides[1]
  if(!is.null(path)) fullname <- file.path(path,fullname)
  id.col = c("Name","Identifier")
  required.col <- unique(c(id.col,columns))
  obj <- read.columnsAE(fullname,required.col,sep="\t",as.is=TRUE,fill=TRUE,flush=TRUE,...)
  nspots = nrow(obj)
        
  ## Initialize RG list object (object.size for matrix of NAs is smaller)
  Y <- matrix(NA,nspots,nslides)
  colnames(Y) = names
  if(length(grep(columns,colnames(obj)[1], ignore.case=TRUE)) == 0)
    {
      if(length(obj[,1][duplicated(obj[,1])])!=0)               
        {
          temp = split(1:length(obj[,1]), obj[,1])
          for(n in names(temp))
            obj[temp[[n]],1] = paste(n, 1:length(temp[[n]]),sep=".")
          rownames(Y) = obj[,1]
        }  else rownames(Y) = obj[,1]

    }

  ## Read remainder of files
  for (i in 1:nslides) {
    fullname <- slides[i]
    if(!is.null(path)) fullname <- file.path(path,fullname)
    obj <- read.columnsAE(fullname,required.col,as.is=TRUE,nrows=nspots,flush=TRUE,...)
       
    Y[,i] <- as.numeric(obj[,grep(columns,colnames(obj), ignore.case=TRUE)])
    if(verbose) cat(paste("Read",fullname,"\n"))          
  }
  new("ExpressionSet",exprs=Y)
}#end of read.AE_1col

# Read 2 colors files from AE into RGList
# Modificated version of read.maimages from limma

read.AE_2col <- function(files=NULL,source="generic",path=NULL,columns=NULL,verbose=TRUE,sep="\t",quote=NULL,...)
# Extracts an RG list from a series of image analysis output files
# Gordon Smyth, modified by Audrey Kauffmann for AE 
{
  ## Begin checking input arguments
  slides <- as.vector(as.character(files))
  nslides <- length(slides)
  names <- removeExt(files)        

  source <- match.arg(source,c("generic","agilent","arrayvision","bluefuse","genepix","genepix.median","genepix.custom","imagene","quantarray","scanarrayexpress","smd.old","smd","spot","spot.close.open","generic1","generic2","generic3","generic4","generic5","generic6","generic7","generic8","generic9","generic10"))
  ## source2 is the source type with qualifications removed
  source2 <- strsplit(source,split=".",fixed=TRUE)[[1]][1]
  if(source=="agilent") quote <- "" else quote <- "\""
  if(source2=="imagene") return(read.imagene(files=files,path=path,ext=ext,names=names,columns=columns,other.columns=other.columns,verbose=verbose,sep=sep,quote=quote,...))

  if(is.null(columns)) {
    if(source2=="generic") stop("must specify columns for generic input")
    columns <- switch(source,
                      agilent = list(G="gMeanSignal$",Gb="gBGMedianSignal$",R="rMeanSignal$",Rb="rBGMedianSignal$"),
                      arrayvision = list(G="MTM Dens - Levels$",Gb="Bkgd$",R="MTM Dens - Levels$",Rb="Bkgd$"),
                      bluefuse = list(G="AMPCH1$",R="AMPCH2$"),
                      genepix = list(R="F635 Mean$",G="F532 Mean$",Rb="B635 Mean$",Gb="B532 Mean$"),
                      quantarray = list(R="ch2 Intensity$",G="ch1 Intensity$",Rb="ch2 Background$",Gb="ch1 Background$"),
                      scanarrayexpress = list(G="Ch1 Mean$",Gb="Ch1 B Median$",R="Ch2 Mean$",Rb="Ch2 B Median$"),
                      smd.old = list(G="CH1I_MEAN$",Gb="CH1B_MEDIAN$",R="CH2I_MEAN$",Rb="CH2B_MEDIAN$"),
                      smd = list(G="Ch1 Intensity (Mean)$",Gb="Ch1 Background (Median)$",R="Ch2 Intensity (Mean)$",Rb="Ch2 Background (Median)$"),
                      spot = list(R="Rmean$",G="Gmean$",Rb="morphR$",Gb="morphG$"),
                      spot.close.open = list(R="Rmean$",G="Gmean$",Rb="morphR.close.open$",Gb="morphG.close.open$"),
                      generic1 = list(R="Mean [QL] Cy5$",G="Mean [QL] Cy3$",Rb="Bkg [QL] Cy5$",Gb="Bkg [QL] Cy3$"),
                      generic2 = list(R="Cy5_Volume$",G="Cy3_Volume$",Rb="Cy5_BG_Tot$",Gb="Cy3_BG_Tot$"),                      
                      generic3 = list(R="Signal med Cy5$",G="Signal med Cy3$",Rb="Background med Cy5$",Gb="Background med Cy3$"),                     
                      generic4 = list(R="Med.Sig.Cy5$",G="Med.Sig.Cy3$",Rb="BG.Cy5$",Gb="BG.Cy3$"),                      
                      generic5 = list(R="Signal Mean_Cy5$",G="Signal Mean_Cy3$",Rb="Background Mean_Cy5$",Gb="Background Mean_Cy3$"),                      
                      generic6 = list(R="Signal_Mean_Cy5$",G="Signal_Mean_Cy3$",Rb="Background_Mean_Cy5$",Gb="Background_Mean_Cy3$"),                     
                      generic7 = list(R="F633 Mean$",G="F543 Mean$",Rb="B633 Mean$",Gb="B543 Mean$"),                      
                      generic8 = list(R="Cy5 Intensity$",G="Cy3 Intensity$",Rb="Cy5 Background$",Gb="Cy3 Background$"),                      
                      generic9 = list(R="Cy5 Mean$",G="Cy3 Mean$",Rb="Cy5 Mean - B$",Gb="Cy3 Mean - B$"),                      
                      generic10 = list(R="Cy5_Intensity$",G="Cy3_Intensity$",Rb="Cy5_Background$",Gb="Cy3_Background$"),
                    NULL
                      )
  } else {
    if(!is.list(columns)) stop("columns must be a list")
    names(columns)[names(columns)=="Gf"] <- "G"
    names(columns)[names(columns)=="Rf"] <- "R"
    if(is.null(columns$G) || is.null(columns$R)) stop("columns must specify foreground G and R")
    if(!all(names(columns) %in% c("G","R","Gb","Rb"))) warning("non-standard columns specified")
  }
  cnames <- names(columns)
	
  annotation <- switch(source,
                       agilent = c("Row","Col","Start","Sequence","SwissProt","GenBank","Primate","GenPept","ProbeUID","ControlType","ProbeName","GeneName","SystematicName","Description"),
                       arrayvision = "Spot labels",
                       bluefuse = c("ROW","COL","SUBGRIDROW","SUBGRIDCOL","BLOCK","NAME","ID"),   
                       genepix=,genepix.median=,genepix.custom = c("Block","Row","Column","ID","Name"),
                       quantarray= c("Array Row","Array Column","Row","Column","Name"),
                       scanarrayexpress = c("Array Row","Array Column","Spot Row","Spot Column"), 	
                       smd = c("Spot","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","Locuslink ID","Name","Sequence Type","X Grid Coordinate (within sector)","Y Grid Coordinate (within sector)","Sector","Failed","Plate Number","Plate Row","Plate Column","Clone Source","Is Verified","Is Contaminated","Luid"),
                       smd.old = c("SPOT","NAME","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","SUID"),
                       NULL
                       )

  ## End checking input arguments

  ## Read first file to get nspots
  fullname <- slides[1]
  if(!is.null(path)) fullname <- file.path(path,fullname)
  required.col <- unique(c(annotation,unlist(columns)))
  obj <- read.columnsAE(fullname,required.col,sep="\t",as.is=TRUE,fill=TRUE,flush=TRUE,...)
  nspots = nrow(obj)

  ## Initialize RG list object (object.size for matrix of NAs is smaller)
  Y <- matrix(NA,nspots,nslides)
  colnames(Y) <- names
  RG <- columns
  for (a in cnames) RG[[a]] <- Y
  RG$targets <- data.frame(FileName=I(files),row.names=names)

  ## Set annotation columns
  if(!is.null(annotation)) {
    j <- match(annotation,colnames(obj),0)
    if(any(j>0)) RG$genes <- data.frame(obj[,j,drop=FALSE],check.names=FALSE)
  }

  RG$source <- source

  ## Set printer layout, if possible
  if(source2=="agilent") {
    if(!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
      nr <- length(unique(RG$genes$Row))
      nc <- length(unique(RG$genes$Col))
      if(nspots==nr*nc) RG$printer <- list(ngrid.r=1,ngrid.c=1,nspot.r=nr,nspot.c=nc)
    }
  }
  if(source2=="genepix") {
    if(!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && !is.null(RG$genes$Column)) {
      RG$printer <- getLayout(RG$genes,guessdups=FALSE)
      nblocks <- RG$printer$ngrid.r*RG$printer$ngrid.c
      if(!is.na(nblocks) && (nblocks>1) && !is.null(obj$X)) {
        blocksize <- RG$printer$nspot.r*RG$printer$nspot.c
        i <- (1:(nblocks-1))*blocksize
        ngrid.r <- sum(obj$X[i] > obj$X[i+1]) + 1
        if(!is.na(ngrid.r) && nblocks%%ngrid.r==0) {
          RG$printer$ngrid.r <- ngrid.r
          RG$printer$ngrid.c <- nblocks/ngrid.r
        } else {
          warning("Can't determine number of grid rows")
          RG$printer$ngrid.r <- RG$printer$ngrid.c <- NA
        }
      }
    }
  }

  ## Read remainder of files
  for (i in 1:nslides) {
    if(i > 1) {
      fullname <- slides[i]
      if(!is.null(path)) fullname <- file.path(path,fullname)
      if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
      obj <- read.columnsAE(fullname,required.col,sep="\t",as.is=TRUE,fill=TRUE,flush=TRUE)      
    }
    for (a in cnames) RG[[a]][,i] <- obj[,grep(unlist(columns)[a],colnames(obj))]
    if(verbose) cat(paste("Read",fullname,"\n"))
  }
  new("RGList",RG)
}#end of read.AE_2col


read.columnsAE = function(file,required.col=NULL,sep="\t",quote="\"",skip=0,fill=TRUE,blank.lines.skip=TRUE,allowEscapes=FALSE,...)
# Read specified columns from a delimited text file with header line
# Gordon Smyth, modified by Audrey Kauffmann for AE
{
  ## Read header to get column names
  allcnames <- scan(file,what="",sep=sep,quote=quote,nlines=1,quiet=TRUE,skip=skip,strip.white=TRUE,blank.lines.skip=blank.lines.skip,allowEscapes=allowEscapes)
  ncn <- length(allcnames)
  colClasses <- rep("NULL",ncn)

  ## Are required columns in the header?
  if(is.numeric(required.col)) {
    colClasses[required.col] <- NA
  } else {
    colClasses[unlist(lapply(1:length(required.col),function(i) grep(required.col[i],allcnames, ignore.case=TRUE)))] <- NA
  }

  ## Read specified columns
  read.table(file=file,header=TRUE,col.names=allcnames,check.names=FALSE,colClasses=colClasses,sep=sep,quote=quote,skip=skip,fill=fill,blank.lines.skip=blank.lines.skip,allowEscapes=allowEscapes,...)
}#end of read.columnsAE



