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
AB = function(i, files, path, ph, adr)
  {
    if(adr == "Empty" || is.na(adr))
      {
        files = files[files!=""]
        raweset = try(ReadAffy(filenames = paste(path,unique(files),sep="/")))
      } else {  
        pht = pData(ph)
        if(!"Array.Data.File" %in% colnames(pht))
          warning("Cannot find array data file names in the sdrf file. The object may not be built.")
    
        if(length(pht[pht$Array.Design.REF==adr[i],"Array.Data.File"]!="")==0)
          warning("Cannot find array data file names in the sdrf file. The object may not be built.") else files = pht[pht$Array.Design.REF==adr[i],"Array.Data.File"]
     
        files = files[files!=""]

        pht = ph[pht[,"Array.Design.REF"]==adr[i],]
        raweset = try(ReadAffy(filenames = paste(path,unique(files),sep="/")))
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
  }#end of AffyBatch

## Create NCS or ES for non Affymetrix data sets
nonAB = function(i, files, path, ph, rawcol, adr)
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
    allcnames = scan(paste(path,files[1],sep="/"),what="",nlines=1, sep="\t")
        
    ## Looking for the right column to use
    scanname = allcnames[grep(":",allcnames)]
    if(length(grep("Database|Reporter",scanname)) != 0)
      scanname = scanname[-grep("Database|Reporter",scanname)] ##Feature is a problem because of Feature Extraction
    feature = grep("Feature",scanname)
    fe = grep("Feature Extraction",scanname)
    if(length(feature) != 0)
      scanname = scanname[-feature[!feature %in% fe]]
    if(is.null(rawcol) && length(scanname) == 0)
      stop(sprintf("No scanner name is given. It is not possible to handle such a case. Try to set the argument 'rawcol' by choosing among the following columns names: \n") ,sprintf("\"%s\" \n",allcnames))

    st = unique(sapply(seq_len(length(scanname)), function(i) strsplit(scanname,":")[[i]][1]))
       
    if(!is.null(rawcol))
      {
        if(length(rawcol) == 1)
          rawcoln = list(G=rawcol,R=rawcol)
        if(length(rawcol) > 1)
          rawcoln = rawcol
            
        ev = try(read.maimages(files=unique(files), path=path,columns=rawcoln))
            
        if(inherits(ev, 'try-error'))
          stop(sprintf("Error in read.maimages: %s The files have probably different columns header.", ev[1]))
            
        if(length(rawcol) > 1)
          raweset = build.ncs(ev,pht,files,raweset)
            
        if(length(rawcol) == 1)
          raweset = build.es(ev,pht,files,raweset)
      }
        
    if(is.null(rawcol))
      {
        if(length(grep(st, scanners)) == 0)
          stop(sprintf("Scanner name is '%s'. This scanner type is not valid. \nTry to set the argument 'rawcol' by choosing among the following columns names: \n", st),sprintf("\"%s\" \n",scanname))
            
        if(length(st) != 1)
          stop(sprintf("%s scanner names are given ( ",length(st)), sprintf("\"%s\" ",st), sprintf("). It is not possible to handle such a case. Try to set the argument 'rawcol' by choosing among the following columns names: \n") ,sprintf("\"%s\" \n",scanname))
        
        if(length(grep(st, scanners)) > 1)
          stop(sprintf("Scanner name can be '%s'. \nTry to set the argument 'rawcol' by choosing among the following columns names: \n", scanners[grep(st, scanners)]),sprintf("\"%s\" \n",scanname))
            
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
            rawcol = if(db[1] == 2) list(R=colnamesf[colnamesf[,2]=="Cy5",1], G=colnamesf[colnamesf[,2]=="Cy3",1],Rb=colnamesb[colnamesb[,2]=="Cy5",1], Gb=colnamesb[colnamesb[,2]=="Cy3",1]) else list(R=colnamesf[colnamesf[,2]=="Cy5",1], G=colnamesf[colnamesf[,2]=="Cy3",1])
            
            if(length(rawcol) == 0 || (0 %in% sapply(seq_len(length(rawcol)), function(i) length(rawcol[[i]]))))
              stop(sprintf("The known column names for this scanner are not in the heading of the files.\nTry to set the argument 'rawcol' by choosing among the following columns names: \n"),sprintf("\"%s\" \n",scanname))

            ev = try(read.maimages(files=unique(files), path=path,rawcol=rawcol))
            if(inherits(ev, 'try-error'))
              stop(sprintf("Error in read.maimages: %s.", ev[1]))
            raweset = build.ncs(ev,pht,files,raweset)
          }
        
        ## Building ExpressionSet when one colour
        if(df[1] == 1)
          {
            rawcol = list(R=colnamesf[,1], G=colnamesf[,1])
            ev = try(read.maimages(files=unique(files), path=path,rawcol=rawcol))
            if(inherits(ev, 'try-error'))
              stop(sprintf("Error in read.maimages: %s", ev[1]))
            raweset = build.es(ev,pht,files,raweset)

          }

        if(df[1] == 0)
          stop(sprintf("None of the columns names of the expression files is matching a valid known quantitation type.\nTry to set the argument 'rawcol' by choosing among the following columns names: \n"), sprintf("\"%s\" \n",scanname))
        if(df[1] > 2)
          stop(sprintf("There are too many columns that could be read in the files.\n Try to set the argument 'rawcol' by choosing among the following columns names: \n"),sprintf("\"%s\" \n",scanname))
     }
    return(raweset)
  }#end of non Affymetrix objects

## Assign experiment Data
## By Juok Cho
creating_experiment = function(idf, eset, path)
  {
    idffile = scan(file.path(path,idf),character(),sep="\n")
    idf.data = list()
    for(g in idffile) { 
      e = unlist(strsplit(g,"\t"))
      key=e[1] 
      if(length(e)>1)
        values=e[2:length(e)]
      else
        values=NA
      idf.data[[key]]=values
    }
    
    ## making key matches #
    ## (Person Last Name, Person First Name, Person Mid Initials), Person Email, Person Phone, Person Address, Person Affiliation, 
    Person_Name = c(idf.data$"Person First Name", idf.data$"Person Last Name",idf.data$"Person Mid Initials") 
    Personal_contact=c(idf.data$"Person Email", idf.data$"Person Phone", idf.data$"Person Address")
	  
    ## making experimentData object #		
    SubmitterIndex=which(idf.data$"Person Roles"=="submitter")
    experimentData = new("MIAME", 
      name=paste(idf.data$"Person Last Name"[SubmitterIndex],", ",idf.data$"Person First Name"[SubmitterIndex], sep=""), #performer
      lab = idf.data$"Person Affiliation"[SubmitterIndex] , #Person Affiliation 
      contact = idf.data$"Person Email"[SubmitterIndex], # Person Email(Person Phone, Person Address)
      title = idf.data$"Investigation Title" , #description #Investigation Title
      ##abstract= "",	#not provided in the idf.data
      ##url	= "",
      other = list(
        accession = gsub(".sdrf.txt","",idf.data$"SDRF File"), #Experiment Name
        identifier = gsub(".sdrf.txt","",idf.data$"SDRF File"), #Experiment Name
        ##Experimental Factor Type
        experimentalFactor = c(idf.data$"Experimental Factor Type"), 
        ##Experimental Design
        type = c(idf.data$"Experimental Design"),
        measurementType = experimentData(eset)@other$measurementType
        )
      )
    experimentData(eset) = experimentData
    return(eset)	  
  }

