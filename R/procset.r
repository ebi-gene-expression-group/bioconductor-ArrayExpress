procset = function(files, procol)
{
  stopifnot(length(procol)==1)
  with(files, {

    if(length(processedFiles) == 1)
      {
        proctot = try(read.table(file.path(path,processedFiles),header=TRUE,sep="\t", quote="", row.names=1, stringsAsFactors=FALSE))
        if(inherits(proctot, 'try-error'))
          if(length(grep("duplicate",proctot)) != 0) stop("The probe identifiers are not unique. The processed file cannot automatically be treated.") else stop("Cannot read the processed file automatically.") 
      }
      
     if(length(processedFiles) > 1)
       stop("The processed files contain different numbers/subsets of reporters and cannot be automatically assembled.")
   
    procsel = matrix(as.numeric(as.matrix(proctot[-1, procol == proctot[1,]])), nrow=nrow(proctot)-1)
    colnames(procsel) = colnames(proctot[, procol == proctot[1,]])
    rownames(procsel) = rownames(proctot[-1,])
   
    proceset = new("ExpressionSet", exprs = procsel)

    ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))
    if("Array.Data.Matrix.File" %in% colnames(pData(ph)) || "Array.Data.File" %in% colnames(pData(ph)))
{
    if("Array.Data.File" %in% colnames(pData(ph)))
    {
      ph = ph[which(pData(ph)$Array.Data.File !=""),]
      sampleNames(ph) = ph$Array.Data.File
      }
    if("Array.Data.Matrix.File" %in% colnames(pData(ph)))
    {
       ph = ph[which(pData(ph)$Array.Data.Matrix.File !=""),]
      sampleNames(ph) = ph$Array.Data.Matrix.File
      }
    if(all(sampleNames(proceset) %in% ph$Array.Data.File) || all(sampleNames(proceset) %in% ph$Array.Data.Matrix.File))
    {
	ph = ph[sampleNames(proceset),]
	phenoData(proceset) = ph
    } else warning("Cannot attach phenoData")
} else warning("Cannot attach phenoData")
    procesetex = try(creating_experiment(idf = idf, eset = proceset, path = path)) 
    if(!inherits(procesetex, 'try-error'))
      proceset = procesetex else warning("Cannot attach experimentData")
    
    procesetex = try(addADFproc(adf = adf, eset = proceset, path = path, procfile=processedFiles))
    if(!inherits(procesetex, 'try-error'))
      proceset = procesetex else warning("Cannot attach featureData")
    
    if(!validObject(proceset)) warning(validObject(proceset))
    return(proceset)
  }) ## with
}
