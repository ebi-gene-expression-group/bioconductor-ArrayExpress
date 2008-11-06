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
    
    xml = xmlTreeParse(queryfilename)
    ID = sapply(1:length(xmlRoot(xml)), function(i) unlist(xmlElementsByTagName(xmlRoot(xml)[[i]], "accession"))["accession.children.text.value"])
    names(ID) = NULL
    xml2 = xmlTreeParse(queryfilename, useInternalNodes = TRUE)
    ra = getNodeSet(xml2,"/experiments//raw[@count]")
    Raw = sapply(ra, function(r) xmlGetAttr(r, "count"))
    
    pr = getNodeSet(xml2,"/experiments//fgem[@count]")
    Processed = sapply(pr, function(p) xmlGetAttr(p, "count"))
    
    Raw[Raw != "0"] = "Yes"
    Raw[Raw == "0"] = "No"
    Processed[Processed != "0"] = "Yes"
    Processed[Processed == "0"] = "No"
    
    xmlparsed = data.frame(ID = ID, Raw = Raw, Processed = Processed)
    return(xmlparsed)
  }
