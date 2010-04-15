procset = function(files, procol)
{
  stopifnot(length(procol)==1)
  with(files, {
    
    proctot = read.table(file.path(path,procfile),header=TRUE,sep="\t", quote="", row.names=1, stringsAsFactors=FALSE)
    procsel = matrix(as.numeric(as.matrix(proctot[-1, procol == proctot[1,]])), nrow=nrow(proctot)-1)
    colnames(procsel) = colnames(proctot[, procol == proctot[1,]])
    rownames(procsel) = rownames(proctot[-1,])
   
    proceset = new("ExpressionSet",
      exprs = procsel)

    ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))
    ph = ph[which(pData(ph)$Array.Data.File !=""),]
    sampleNames(ph) = ph$Array.Data.File
    if(all(sampleNames(proceset) %in% ph$Array.Data.File))
    {
	ph = ph[sampleNames(proceset),]
	phenoData(proceset) = ph
    } else warning("Cannot attach phenoData")
    procesetex = try(creating_experiment(idf = idf, eset = proceset, path = path)) 
    if(!inherits(procesetex, 'try-error'))
      proceset = procesetex else warning("Cannot attach experimentData")
    
    procesetex = try(addADFproc(adf = adf, eset = proceset, path = path, procfile=procfile))
    if(!inherits(procesetex, 'try-error'))
      proceset = procesetex else warning("Cannot attach featureData")
    
    if(!validObject(proceset)) warning(validObject(proceset))
    return(proceset)
  }) ## with
}
