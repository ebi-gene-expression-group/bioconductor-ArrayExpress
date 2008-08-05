ArrayExpress = function(input, tempoutdir = ".", save = FALSE, columns = NULL)
  {
    if(!is.null(columns))
      {
        if(length(columns) > 1 && !is(columns,"list"))
          stop("The argument 'columns' must be a list if multiple column name are given.")
        if(length(columns) == 1 && !is(columns,"character"))
          stop("The argument 'columns' must be a character if one column name is given.")
        if(is(columns,"list") && !("R" %in% names(columns) && "G" %in% names(columns)))
          stop("The names of the columns must contain R and G.")
      }

    if(save) on.exit(file.remove(rawdata)) else on.exit({file.remove(rawdata);try(file.remove(allfiles))})
    
    ## Building the link with the input name
    ##dircontentold = dir()[grep(input,dir())]
    dir = gsub("^E-|-[0-9]{1,10}","",input)
    if(nchar(dir) == 5)
      dir = gsub("[a-z]$","",dir)
    url = "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment"
    exp = paste(url,dir,input,input,sep="/")
    raw = paste(exp,".raw.zip",sep="")
    ## Saving temporarily the raw data
    rawdata = tempfile(tmpdir = tempoutdir)
    dnld = try(download.file(raw, rawdata, mode="wb"))

    ## Download checking
    if(inherits(dnld, 'try-error'))
      stop(paste(raw, " does not exist. \n"),sep="")
    if(file.info(rawdata)$size == 0)
      stop(paste(raw, " is empty. \n"),sep="")
    
    ## Extracting the data if the checking was fine
    allfiles = extract.zip(file = rawdata)
   
    ## Listing the expression files names
    ##allfiles = dir(tempoutdir, pattern = input)
    ##allfiles= allfiles[!allfiles %in% dircontentold]

    ##dircontent = dir(tempoutdir, pattern = input) #list all the files from the current directory
    ##directories = which(file.info(dircontent)$isdir ==TRUE) #list the directories    
    notuse = grep(allfiles, pattern = "info.txt$|idf.txt$|processed|sdrf.txt$|.log$|RData|class|log")
    if(length(notuse) != 0)
      files = allfiles[-notuse]
    if(length(notuse) == 0)
      files = allfiles
    
    samples = paste(exp,".sdrf.txt",sep="")
    ph = try(read.AnnotatedDataFrame(samples,row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))
        
    if(inherits(ph, 'try-error') && length(grep(".cel",files)) == 0)
      stop(sprintf("No sdrf file available. The object cannot be built."))
    if(inherits(ph, 'try-error') && length(grep(".cel",files)) != 0)
      {
        warning(sprintf("No sdrf file available. The object may not be built."))
        adr = "Empty"
        raweset = try(AB(i=1,input, files, tempoutdir, ph, adr))

      }
    if(!inherits(ph, 'try-error'))
      {
        emptylines = which(sapply(seq_len(nrow(pData(ph))), function(i)
          all(pData(ph)[i,]=="",na.rm=T)))
        if(length(emptylines) != 0)
          pData(ph) = pData(ph)[-emptylines,]
        ##pData(ph) = pData(ph)[!pData(ph)$Array.Data=="",]
    
        if(!all(files %in% pData(ph)$Array.Data))
          warning("Some files from the zip archive do not have annotation in the sdrf file. The object may not be built.")

        adr = unique(pData(ph)$Array.Design.REF)
        adr = adr[adr!=""]
        if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",files)) != 0)
          warning("Cannot find the array design reference in the sdrf file. The object may not be built.")
         if((length(adr) == 0 || is.na(adr)) && length(grep(".cel",files)) == 0)
          stop("Cannot find the array design reference in the sdrf file. The object cannot be built.")
     }
    
    ## Building the S4 class object
    if(length(grep(".cel",files)) == length(files))
      {
        if(length(adr) == 1)
          raweset = try(AB(i=1,input, files, tempoutdir, ph, adr))

        if(length(adr) > 1)
          raweset = try(lapply(seq_len(length(adr)), function(i) try(AB(i,input, files, tempoutdir, ph, adr))))
      }
   
    ## Non Affymetrix data
    if(length(grep(".cel",files)) == 0)
      {
        if(length(adr) == 1)
          raweset = try(nonAB(i=1,input, files, tempoutdir, ph, columns, adr))
        if(length(adr) > 1)
          raweset = try(lapply(seq_len(length(adr)), function(i) try(nonAB(i,input, files, tempoutdir, ph, columns, adr))))
      }    
    return(raweset)
  }#end of ArrayExpress


## Build NChannelSet
build.ncs = function(ev,ph,files,raweset)
  {
    assayData = if("Rb" %in% names(ev))with(ev, assayDataNew(R=R, G=G, Rb=Rb, Gb=Gb)) else with(ev, assayDataNew(R=R, G=G))

    raweset = try(new("NChannelSet", assayData = assayData))

    if(inherits(raweset, 'try-error'))
      stop(sprintf(raweset[1]))
    
    if(!inherits(ph, 'try-error'))
      {
        rawesetph = try(assign.pheno.ncs(files=files,ph=ph,raweset=raweset))
        if(!inherits(rawesetph, 'try-error'))
          raweset = rawesetph
        
        if(inherits(rawesetph, 'try-error'))
          warning(sprintf("Cannot attach phenoData to the object as the Array Data File column of the sdrf file contains duplicated elements."))
      }
    return(raweset)       
  }

## Build ExpressionSet
build.es = function(ev,ph,files,raweset)
  {
    raweset = try(new("ExpressionSet",
      exprs = ev$G))

    if(inherits(raweset, 'try-error'))
      stop(sprintf(raweset[1]))
            
    if(!inherits(ph, 'try-error'))
      {
        pData(ph) = pData(ph)[1:length(files),]
        if(length(ph$Array.Data)==length(unique(ph$Array.Data)))
          {
            rownames(pData(ph)) = gsub(".[a-z][a-z][a-z]$","",ph$Array.Data)
            pData(ph) = pData(ph)[sampleNames(raweset),]
            phenoData(raweset) = ph
          } else {
            warning(sprintf("Cannot attach phenoData to the object as the Array Data File column of the sdrf file contains duplicated elements.")) 
          }
      }
    return(raweset)       
  }

## assign phenoData to Nchannelset
assign.pheno.ncs = function(ph,files,raweset)
  {
    si = pData(ph)[1:(length(files)*2),]
    lab = split(si,si[,'Label'])
    same = which(lapply(1:ncol(lab[[1]]), function(i) all(lab[[1]][i]==lab[[2]][i]))==TRUE)
    all = lab[[1]][same]
    gspe = lab[[1]][-same]
    colnames(gspe) = paste(colnames(gspe),names(lab)[1],sep=".")
    rspe = lab[[2]][-same]
    colnames(rspe) = paste(colnames(rspe),names(lab)[2],sep=".")
                
    metaData = data.frame(labelDescription=c(rep("_ALL_",ncol(all)),rep("G",ncol(gspe)),rep("R",ncol(rspe))))
                            
    samples = new("AnnotatedDataFrame", data=cbind(all,gspe,rspe), varMetadata=metaData)

    rownames(pData(samples)) = gsub(".[a-z][a-z][a-z]$","",samples$Array.Data)
    pData(samples) = pData(samples)[sampleNames(raweset)$R,]
    phenoData(raweset) = samples
    return(raweset)
  }


## Create AffyBatch for Affymetrix data sets
AB = function(i,input, files, tempoutdir, ph, adr)
  {
    if(adr == "Empty" || is.na(adr))
      {
        files = files[files!=""]
        raweset = try(ReadAffy(filenames = paste(tempoutdir,unique(files),sep="/")))
      } else {  
        pht = pData(ph)
        if(!"Array.Data.File" %in% colnames(pht))
          warning("Cannot find array data file names in the sdrf file. The object may not be built.")
    
        if(length(pht[pht$Array.Design.REF==adr[i],"Array.Data.File"]!="")==0)
          warning("Cannot find array data file names in the sdrf file. The object may not be built.") else files = pht[pht$Array.Design.REF==adr[i],"Array.Data.File"]
     
        files = files[files!=""]

        pht = ph[pht[,"Array.Design.REF"]==adr[i],]
        raweset = try(ReadAffy(filenames = paste(tempoutdir,unique(files),sep="/")))
      }
    if(!inherits(raweset, 'try-error'))
      {
        if(!inherits(ph, 'try-error'))
          {
            pData(ph) = pData(ph)[1:length(files),]
            if(length(ph$Array.Data)==length(unique(ph$Array.Data)))
              {
                rownames(pData(ph)) = ph$Array.Data
                pData(ph) = pData(ph)[sampleNames(raweset),]
                phenoData(raweset) = ph
              } else {
                warning(sprintf("Cannot attach phenoData to the object as the Array Data File column of the sdrf file contains duplicated elements.")) 
              }              
          }
      }
    return(raweset)
  }

## Create NCS or ES for non Affymetrix data sets
nonAB = function(i, input, files, tempoutdir, ph, columns, adr)
  {
    pht = pData(ph)
    if(!"Array.Data.Matrix.File" %in% colnames(pht))
      warning("Cannot find array data file names in the sdrf file. The object may not be built.")


    if(length(pht[pht$Array.Design.REF==adr[i],"Array.Data.Matrix.File"]!="")==0)
      warning("Cannot find array data file names in the sdrf file. The object may not be built.") else files = pht[pht$Array.Design.REF==adr[i],"Array.Data.Matrix.File"]
     
    files = files[files!=""]
    pht = ph[pht$Array.Design.REF==adr[i],]
  
    url2 = "http://tab2mage.svn.sourceforge.net/viewvc/*checkout*/tab2mage/trunk/Tab2MAGE/lib/ArrayExpress/Datafile/QT_list.txt" 
        
    qt = try(read.table(url2, sep = "\t", quote = "",
      check.names = FALSE, fill = TRUE,
      comment.char = "#",               
      stringsAsFactors =  FALSE)) ##read the QT file from the web
    if(inherits(qt, 'try-error'))
      qt = try(read.table(
        file.path(system.file("doc", package="ArrayExpress"),"QT_list.txt"),
        sep = "\t", quote = "",
        check.names = FALSE, fill = TRUE,
        comment.char = "#",               
        stringsAsFactors =  FALSE)) ##read the QT file from the package

        
    scanners = grep(">>>",qt[,1],value=T) ## list all the scanner names
    sl = grep(">>>",qt[,1]) ## list all the line numbers wherea scanner type starts
    scanners = gsub(">","",scanners)
    
    ## Parsing the first line of the expression file
    allcnames = scan(paste(tempoutdir,files[1],sep="/"),what="",nlines=1, sep="\t")
        
    ## Looking for the right column to use
    scanname = allcnames[grep(":",allcnames)]
    if(length(grep("Database|Reporter",scanname)) != 0)
      scanname = scanname[-grep("Database|Reporter",scanname)] ##Feature is a problem because of Feature Extraction
    feature = grep("Feature",scanname)
    fe = grep("Feature Extraction",scanname)
    if(length(feature) != 0)
      scanname = scanname[-feature[!feature %in% fe]]
    if(is.null(columns) && length(scanname) == 0)
      stop(sprintf("No scanner name is given. It is not possible to handle such a case. Try to set the argument 'columns' by choosing among the following columns names: \n") ,sprintf("\"%s\" \n",allcnames))

    st = unique(sapply(seq_len(length(scanname)), function(i) strsplit(scanname,":")[[i]][1]))
    
    #tw =  read.table(paste("tempwcl",input,sep=""),nrow=length(files),sep="E")
    #nlines = unique(tw[,1])
    #file.remove(paste("tempwcl",input,sep=""))    
    #if(length(nlines)>1)
      #stop(sprintf("The files have different number of lines."))         

    
    if(!is.null(columns))
      {
        if(length(columns) == 1)
          columnsn = list(G=columns,R=columns)
        if(length(columns) > 1)
          columnsn = columns
            
        ev = try(read.maimages(files=unique(files), path=tempoutdir,columns=columnsn))
            
        if(inherits(ev, 'try-error'))
          stop(sprintf("Error in read.maimages: %s The files have probably different columns header.", ev[1]))
            
        if(length(columns) > 1)
          raweset = build.ncs(ev,pht,files,raweset)
            
        if(length(columns) == 1)
          raweset = build.es(ev,pht,files,raweset)
      }
        
    if(is.null(columns))
      {
        if(length(grep(st, scanners)) == 0)
          stop(sprintf("Scanner name is '%s'. This scanner type is not valid. \nTry to set the argument 'columns' by choosing among the following columns names: \n", st),sprintf("\"%s\" \n",scanname))
            
        if(length(st) != 1)
          stop(sprintf("%s scanner names are given ( ",length(st)), sprintf("\"%s\" ",st), sprintf("). It is not possible to handle such a case. Try to set the argument 'columns' by choosing among the following columns names: \n") ,sprintf("\"%s\" \n",scanname))
        
        if(length(grep(st, scanners)) > 1)
          stop(sprintf("Scanner name can be '%s'. \nTry to set the argument 'columns' by choosing among the following columns names: \n", scanners[grep(st, scanners)]),sprintf("\"%s\" \n",scanname))
            
        gs = qt[((sl[grep(st,scanners)]+1):(sl[grep(st,scanners)+1]-1)),] ## extract the QTs of the specific scanner type
        foreground = gs[(gs[,4]=="MeasuredSignal" & (is.na(gs[,5]) | gs[,5]==0)),c(1,7)] ## the colnames to use in the read.column
        background = gs[(gs[,4]=="MeasuredSignal" & (gs[,5]==1)),c(1,7)] ## the colnames to use in the read.column

        foreground[,1] = paste(st,":",foreground[,1],sep="")
        colnamesf = foreground[which(foreground[,1] %in% allcnames),]
        df = dim(colnamesf)        

        if(dim(background)[1] != 0)
          {
            background[,1] = paste(st,":",background[,1],sep="")
            colnamesb = background[which(background[,1] %in% allcnames),]
            db = dim(colnamesb)
          } else db = 0

        if(length(files) != 1)
          {
            if(!all(sapply(2:length(files), function(i) readLines(files[1],1)==readLines(files[i],1))))
              warning(sprintf("The files do not all have the same headings whereas the array design is the same. It may cause the object not being created."))
          }
        
        ## Building NChannelSet when two colours
        if(df[1] == 2)
          {
            columns = if(db[1] == 2) list(R=colnamesf[colnamesf[,2]=="Cy5",1], G=colnamesf[colnamesf[,2]=="Cy3",1],Rb=colnamesb[colnamesb[,2]=="Cy5",1], Gb=colnamesb[colnamesb[,2]=="Cy3",1]) else list(R=colnamesf[colnamesf[,2]=="Cy5",1], G=colnamesf[colnamesf[,2]=="Cy3",1])
            
            if(length(columns) == 0 || (0 %in% sapply(seq_len(length(columns)), function(i) length(columns[[i]]))))
              stop(sprintf("The known column names for this scanner are not in the heading of the files.\nTry to set the argument 'columns' by choosing among the following columns names: \n"),sprintf("\"%s\" \n",scanname))

            ev = try(read.maimages(files=unique(files), path=tempoutdir,columns=columns))
            if(inherits(ev, 'try-error'))
              stop(sprintf("Error in read.maimages: %s.", ev[1]))
            raweset = build.ncs(ev,pht,files,raweset)
          }
        
        ## Building ExpressionSet when one colour
        if(df[1] == 1)
          {
            columns = list(R=colnamesf[,1], G=colnamesf[,1])
            ev = try(read.maimages(files=unique(files), path=tempoutdir,columns=columns))
            if(inherits(ev, 'try-error'))
              stop(sprintf("Error in read.maimages: %s", ev[1]))
            raweset = build.es(ev,pht,files,raweset)

          }

        if(df[1] == 0)
          stop(sprintf("None of the columns names of the expression files is matching a valid known quantitation type.\nTry to set the argument 'columns' by choosing among the following columns names: \n"), sprintf("\"%s\" \n",scanname))
        if(df[1] > 2)
          stop(sprintf("There are too many columns that could be read in the files.\n Try to set the argument 'columns' by choosing among the following columns names: \n"),sprintf("\"%s\" \n",scanname))
     }
    return(raweset)
  }#end of non Affymetrix objects

