getelt = function(x, node, element)
  {
    elt = sapply(1:length(xmlRoot(x)), function(i){
      if(length(grep(node,names(xmlRoot(x)[[i]])))  != 0)
        unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node))[names(unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node)))== element]  else "NA" })

    sizeelt = sapply(1:length(xmlRoot(x)), function(i){
      if(length(grep(node,names(xmlRoot(x)[[i]])))  != 0)
        length(unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node))[names(unlist(xmlElementsByTagName(xmlRoot(x)[[i]], node)))== element])  else "NA" })

	if(inherits(elt, "list") || (inherits(elt, "character") && max(sizeelt[sizeelt!="NA"]) == 1))
	{
	    elt2 = lapply(elt, function(i) if(length(i) == 0) "NA" else i)

	    elt3 = unlist(lapply(elt2, function(i) do.call("paste",c(as.list(i),sep=" | "))))
    	} 
	if(inherits(elt, "matrix") && max(sizeelt) > 1)
	   elt3 = do.call("paste",c(as.list(elt),sep=" | "))

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

	if(inherits(elt, "list"))
	{
		elt2 = lapply(elt, function(i) if(length(i) == 0) "NA" else i)
		elt3 = lapply(elt2, function(i) unlist(i))
		elt4 = unlist(lapply(elt3, function(i) do.call("paste",c(as.list(i),sep=" | "))))
	} else elt4 = do.call("paste",c(as.list(do.call("paste",c(as.list(elt),sep=" | ")), sep="||")))

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
    Rawids = sapply(ra, function(r) xmlGetAttr(r, "name"))
    Rawids = gsub(".raw.*.zip","",Rawids)
    names(Raw) = Rawids
   
    pr = getNodeSet(x2,"/experiments//fgem[@count]")
    Processed = sapply(pr, function(p) xmlGetAttr(p, "count"))
    Procids = sapply(pr, function(r) xmlGetAttr(r, "name"))
    Procids = gsub(".processed.*.zip","",Procids)
    names(Processed) = Procids

    date = getelt(x, node = "releasedate",  element = "releasedate.children.text.value")
    
    pmid = getelt(x, node = "bibliography",
      element = "bibliography.children.accession.children.text.value")

    spec = getelt(x, node = "species",
      element = "species.children.text.value")

    experimentdesign = getelt(x, node = "experimentdesign",
      element = "experimentdesign.children.text.value")   

    experimentalfactor = geteltmulti(x, node = "experimentalfactor",
      element1 = "children.name.children.text.value",
      element2 = "children.value.children.text.value")

    xmlparsed = data.frame(ID = ID, Raw = Raw[ID], Processed = Processed[ID], ReleaseDate = date, PubmedID = pmid, Species = spec, ExperimentDesign = experimentdesign, ExperimentFactors = experimentalfactor)
    return(xmlparsed)
  }
