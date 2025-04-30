plot.estimateEffect <- function(x, covariate, model=NULL,
                                topics=x$topics,
                                method=c("pointestimate", "difference","continuous"),
                                ref=NULL, alt=NULL,
                                moderator=NULL, moderator.value=NULL,
                                npoints=100, nsims=500, ci.level=.95,
                                xlim=NULL, ylim=NULL, xlab="",ylab=NULL,
                                main="", printlegend=T,
                                labeltype="numbers", n=7, frexw=.5,
                                add=F, linecol="black", width=25,
                                verbose.labels=T, family=NULL,
                                custom.labels=NULL, omit.plot=FALSE,...){
  
  method <- match.arg(method)
  if(method=="difference" && (is.null(ref) | is.null(alt))) {
    # stop("For method='difference' both ref and alt must be specified.")
    if(! covariate %in% names(x$ref.vec)) stop("For method='difference', either the covariate needs to be an covariate specified in the `estimateEffect`, or both ref and alt must be specified.")
    ref.vec <- x$ref.vec
    ref <- ref.vec[match(covariate, names(ref.vec))]
    metadata <- colData(model$settings$sce)
    alt <- setdiff(metadata[[covariate]], ref)
    if(length(alt) != 1) stop("Both ref and alt must be specified, since more than 2 levels are found in metadata")
  }

  #Produce cdata (data in original form) and
  #cmatrix (data in design matrix form)
  cthis <- produce_cmatrix(prep=x, covariate=covariate, method=method,
                           ref=ref,
                           alt=alt, npoints=npoints,
                           moderator=moderator, moderator.value=moderator.value)
  cdata <- cthis$cdata
  cmat <- cthis$cmatrix
  #Simulate betas
  simbetas <- simBetas(x$parameters, nsims=nsims)

  #Find offset for confidence level
  offset <- (1-ci.level)/2

  
  #Plot for each method
  if(method=="continuous"){
    toreturn <- plotContinuous(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,custom.labels=custom.labels,model=model,frexw=frexw,printlegend=printlegend,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
  if(method=="pointestimate"){
    toreturn <- plotPointEstimate(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                      offset=offset,xlab=xlab, ylab=ylab, main=main,
                      xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                      labeltype=labeltype,n=n,
                                  custom.labels=custom.labels,model=model,frexw=frexw,width=width,
                                  verbose.labels=verbose.labels,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
  if(method=="difference"){
    if(missing(ref)) stop("Missing a value for ref. See documentation.")
    if(missing(alt)) stop("Missing a value for alt. See documentation.")
    toreturn <- plotDifference(prep=x,covariate=covariate,topics=topics, cdata=cdata, cmat=cmat, simbetas=simbetas,
                   offset=offset,xlab=xlab, ylab=ylab, main=main,
                   xlim=xlim, ylim=ylim, linecol=linecol, add=add,
                   labeltype=labeltype,n=n,
                   custom.labels=custom.labels, printlegend=printlegend,
                   model=model,frexw=frexw,width=width,
                   ref=ref,
                               alt=alt,verbose.labels=verbose.labels,omit.plot=omit.plot,...)
    return(invisible(toreturn))
  }
}