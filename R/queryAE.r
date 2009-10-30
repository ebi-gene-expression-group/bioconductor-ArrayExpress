getelt = function(x, node, element)
  {
    elt = sapply(1:length(xmlRoot(x)), function(i){
      if(length(grep(node,names(xmlRoot(x)[[i]])))  != 0)
        unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node))[names(unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node)))== element]  else "NA" })

    elt2 = lapply(elt, function(i) if(length(i) == 0) "NA" else i)

    elt3 = unlist(lapply(elt2, function(i) do.call("paste",c(as.list(i),sep=" | "))))
    
    names(elt3) = NULL
    return(elt3)
  }

geteltmulti = function(x, node, element1, element2)
 {
elt = sapply(1:length(xmlRoot(x)), function(i){
     if(length(grep(node,names(xmlRoot(x)[[i]])))  != 0){
          lapply(1:length(xmlElementsByTagName(xmlRoot(x)[[i]], node)), function(j) {
        extr = unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node)[[j]])

        e1 = extr[names(extr)== element1]
        e2 = extr[names(extr)== element2]

        paste(e1, e2, sep="=")
       })} else "NA" })

elt2 = lapply(elt, function(i) if(length(i) == 0) "NA" else i)

elt3 = lapply(elt2, function(i) unlist(i))

elt4 = unlist(lapply(elt3, function(i) do.call("paste",c(as.list(i),sep=" | "))))

   names(elt4) = NULL
   return(elt4)
 }

queryAE = function(keywords = NULL, species = NULL)
  {
    if(!is.null(keywords))
      {
        qr = paste("http://www.ebi.ac.uk/microarray-as/ae/xml/experiments?keywords=",keywords,sep="")
        if(!is.null(species))
          qr = paste(qr,"&species=",species,sep="")
      }
    if(is.null(keywords) && !is.null(species))
      qr = paste("http://www.ebi.ac.uk/microarray-as/ae/xml/experiments?species=",species,sep="")
    
    queryfilename = paste("query",keywords,species,".xml",sep="")
    query = try(download.file(qr, queryfilename, mode="wb"))
    
    x = xmlTreeParse(queryfilename)
    
    ID = sapply(1:length(xmlRoot(x)), function(i) unlist(xmlElementsByTagName(xmlRoot(x)[[i]], "accession"))["accession.children.text.value"])
    names(ID) = NULL
    x2 = xmlTreeParse(queryfilename, useInternalNodes = TRUE)
    ra = getNodeSet(x2,"/experiments//raw[@count]")
    Raw = sapply(ra, function(r) xmlGetAttr(r, "count"))
    
    pr = getNodeSet(x2,"/experiments//fgem[@count]")
    Processed = sapply(pr, function(p) xmlGetAttr(p, "count"))
    
    date = getelt(x, node = "releasedate",
      element = "releasedate.children.text.value")
    
    pmid = getelt(x, node = "bibliography",
      element = "bibliography.children.accession.children.text.value")

    spec = getelt(x, node = "species",
      element = "species.children.text.value")
    
    experimentdesign = getelt(x, node = "experimentdesign",
      element = "experimentdesign.children.text.value")   
    
    ID = do.call("paste",c(as.list(ID),sep=" | "))
    Raw = do.call("paste",c(as.list(Raw),sep=" | "))
    Processed = do.call("paste",c(as.list(Processed),sep=" | "))
    date = do.call("paste",c(as.list(date),sep=" | "))
    pmid = do.call("paste",c(as.list(pmid),sep=" | "))
    spec = do.call("paste",c(as.list(spec),sep=" | "))
    experimentdesign = do.call("paste",c(as.list(experimentdesign),sep=" | "))

    experimentalfactor = geteltmulti(x, node = "experimentalfactor",
      element1 = "children.name.children.text.value",
      element2 = "children.value.children.text.value")

    experimentalfactor = do.call("paste",c(as.list(experimentalfactor),sep=" || "))

    xmlparsed = data.frame(ID = ID, Raw = Raw, Processed = Processed, ReleaseDate = date, PubmedID = pmid, Species = spec, ExperimentDesign = experimentdesign, ExperimentFactors = experimentalfactor)
    return(xmlparsed)
  }
