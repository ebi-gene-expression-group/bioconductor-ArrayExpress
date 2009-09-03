procset = function(files, procol)
  {
    path = files$path
    procfile = files$procfile
    sdrf = files$sdrf
    idf = files$idf
    adf = files$adf

    proctot = read.table(file.path(path,procfile),header=TRUE,sep="\t")
    procsel = proctot[-1, procol == proctot[1,]]
    colnames(procsel) = colnames(proctot[, procol == proctot[1,]])
    rownames(procsel) = proctot[-1,1]

    proceset = new("ExpressionSet",
      exprs = as.matrix(procsel))

    ph = try(read.AnnotatedDataFrame(sdrf, path = path, row.names=NULL, blank.lines.skip = TRUE, fill=TRUE, varMetadata.char="$"))

    phenoData(proceset) = ph

    procesetex = try(creating_experiment(idf = idf, eset = proceset, path = path))
    if(!inherits(procesetex, 'try-error'))
      proceset = procesetex else warning("Cannot attach experimentData")
    
    procesetex = try(addADF(adf = adf, eset = proceset, path = path))
    if(!inherits(procesetex, 'try-error'))
      proceset = procesetex else warning("Cannot attach featureData")   
  }
