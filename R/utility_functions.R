# Not very general function that is made to plot a histogram given discrete (integer)
# data.
discrete_hist <- function(x, xlim, col="skyblue", lwd=3, x_marked = c(), marked_col = "red", yaxt="n",...) {
  hist_data <- hist(x, (xlim[1] - 1):(xlim[2]) + 0.5 , plot=FALSE)
  cols <- ifelse(hist_data$mids %in% x_marked, marked_col, col )
  plot(hist_data$mids, hist_data$density, type="h", col=cols, lwd=lwd, bty = "L",...)
  invisible(hist_data)
}

# Takes coda samples generates a data frame with different statistics for the
# parameters
mcmc_stats <- function(samples, cred_mass = 0.95, comp_val = 0) {
  samples_mat <- as.matrix(samples)
  stats <- data.frame(mean = colMeans(samples_mat))
  stats$sd <- apply(samples_mat, 2, sd)
  hdi_lim <- apply(samples_mat, 2, HDIofMCMC, credMass=cred_mass)
  stats$"HDI%" <- round(cred_mass * 100)
  stats$HDIlo <- hdi_lim[1,]
  stats$HDIup <- hdi_lim[2,]
  stats$comp <- comp_val
  stats$"%>comp" <- apply(samples_mat, 2, function(x) { mean(x > comp_val) })
  stats$"%<comp" <- apply(samples_mat, 2, function(x) { mean(x < comp_val) })
  stats$"q2.5%" <- apply(samples_mat, 2, quantile,  probs= 0.025)
  stats$"q25%" <- apply(samples_mat, 2, quantile,  probs= 0.25)
  stats$median <- apply(samples_mat, 2, median)
  stats$"q75%" <- apply(samples_mat, 2, quantile,  probs= 0.75)
  stats$"q97.5%" <- apply(samples_mat, 2, quantile,  probs= 0.975)
  stats$mcmc.se <- summary(samples)$statistics[,"Time-series SE"]
  stats$Rhat <- gelman.diag(samples)$psrf[, 1]
  stats$n.eff <- effectiveSize(samples)
  as.matrix(stats) # 'cause it's easier to index
}


#' Kruschke
HDIofICDF = function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
  # Arguments:
  #   ICDFname is R's name for the inverse cumulative density function
  #     of the distribution.
  #   credMass is the desired mass of the HDI region.
  #   tol is passed to R's optimize function.
  # Return value:
  #   Highest density iterval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30,12) distribution, type
  #   HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
  #   Notice that the parameters of the ICDFname must be explicitly named;
  #   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
  # Adapted and corrected from Greg Snow's TeachingDemos package.
  incredMass =  1.0 - credMass
  intervalWidth = function( lowTailPr , ICDFname , credMass , ... ) {
    ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
  }
  optInfo = optimize( intervalWidth , c( 0 , incredMass ) , ICDFname=ICDFname ,
                      credMass=credMass , tol=tol , ... )
  HDIlowTailPr = optInfo$minimum
  return( c( ICDFname( HDIlowTailPr , ... ) ,
             ICDFname( credMass + HDIlowTailPr , ... ) ) )
}                  # Kruschke, J. K. (2011). Doing Bayesian data analysis: A
# Tutorial with R and BUGS. Elsevier Science/Academic Press.


#' Kruschke
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = floor( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

#' Kruschke
mcmcSummary = function( paramSampleVec , compVal=NULL ) {
  meanParam = mean( paramSampleVec )
  medianParam = median( paramSampleVec )
  dres = density( paramSampleVec )
  modeParam = dres$x[which.max(dres$y)]
  hdiLim = HDIofMCMC( paramSampleVec )
  if ( !is.null(compVal) ) {
    pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                    / length( paramSampleVec ) )
  } else {
    pcgtCompVal=NA
  }
  return( c( meanParam , medianParam , modeParam , hdiLim , pcgtCompVal ) )
}



#' Author John Kruschke
plotPost = function( param_sample_vec , cred_mass=0.95 , comp_val=NULL ,
                     HDI_text_place=0.7 , ROPE=NULL , yaxt=NULL , ylab=NULL ,
                     xlab=NULL , cex.lab=NULL , cex=NULL , xlim=NULL , main=NULL ,
                     col=NULL , border=NULL , show_mode=F , show_curve=F , breaks=NULL , 
                     ... ) {
  # Override defaults of hist function, if not specified by user:
  # (additional arguments "..." are passed to the hist function)
  if ( is.null(xlab) ) xlab="Parameter"
  if ( is.null(cex.lab) ) cex.lab=1.5
  if ( is.null(cex) ) cex=1.4
  if ( is.null(xlim) ) xlim=range( c( comp_val , param_sample_vec ) )
  if ( is.null(main) ) main=""
  if ( is.null(yaxt) ) yaxt="n"
  if ( is.null(ylab) ) ylab=""
  if ( is.null(col) ) col="skyblue"
  if ( is.null(border) ) border="white"
  
  postSummary = matrix( NA , nrow=1 , ncol=11 , 
                        dimnames=list( c( xlab ) , 
                                       c("mean","median","mode",
                                         "hdiMass","hdiLow","hdiHigh",
                                         "comp_val","pcGTcomp_val",
                                         "ROPElow","ROPEhigh","pcInROPE")))              
  postSummary[,"mean"] = mean(param_sample_vec)
  postSummary[,"median"] = median(param_sample_vec)
  mcmcDensity = density(param_sample_vec)
  postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]
  HDI = HDIofMCMC( param_sample_vec , cred_mass )
  postSummary[,"hdiMass"]=cred_mass
  postSummary[,"hdiLow"]=HDI[1]
  postSummary[,"hdiHigh"]=HDI[2]
  
  # Plot histogram.
  if ( is.null(breaks) ) {
    breaks = c( seq( from=min(param_sample_vec) , to=max(param_sample_vec) ,
                     by=(HDI[2]-HDI[1])/18 ) , max(param_sample_vec) )
  }
  if ( !show_curve ) {
    par(xpd=NA)
    histinfo = hist( param_sample_vec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                     freq=F , border=border , col=col ,
                     xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                     breaks=breaks , ... )
  }
  if ( show_curve ) {
    par(xpd=NA)
    histinfo = hist( param_sample_vec , plot=F )
    densCurve = density( param_sample_vec , adjust=2 )
    plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
          xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
          main=main , cex=cex , cex.lab=cex.lab , ... )
  }
  cenTendHt = 0.9*max(histinfo$density)
  cvHt = 0.7*max(histinfo$density)
  ROPEtextHt = 0.55*max(histinfo$density)
  # Display mean or mode:
  if ( show_mode==F ) {
    meanParam = mean( param_sample_vec )
    text( meanParam , cenTendHt ,
          bquote(mean==.(signif(meanParam,3))) , adj=c(.5,0) , cex=cex )
  } else {
    dres = density( param_sample_vec )
    modeParam = dres$x[which.max(dres$y)]
    text( modeParam , cenTendHt ,
          bquote(mode==.(signif(modeParam,3))) , adj=c(.5,0) , cex=cex )
  }
  # Display the comparison value.
  if ( !is.null( comp_val ) ) {
    cvCol = "darkgreen"
    pcgtcomp_val = round( 100 * sum( param_sample_vec > comp_val )
                          / length( param_sample_vec )  , 1 )
    pcltcomp_val = 100 - pcgtcomp_val
    lines( c(comp_val,comp_val) , c(0.96*cvHt,0) ,
           lty="dotted" , lwd=1 , col=cvCol )
    text( comp_val , cvHt ,
          bquote( .(pcltcomp_val)*"% < " *
                   .(signif(comp_val,3)) * " < "*.(pcgtcomp_val)*"%" ) ,
          adj=c(pcltcomp_val/100,0) , cex=0.8*cex , col=cvCol )
    postSummary[,"comp_val"] = comp_val
    postSummary[,"pcGTcomp_val"] = ( sum( param_sample_vec > comp_val ) 
                                     / length( param_sample_vec ) )
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    ropeCol = "darkred"
    pcInROPE = ( sum( param_sample_vec > ROPE[1] & param_sample_vec < ROPE[2] )
                 / length( param_sample_vec ) )
    lines( c(ROPE[1],ROPE[1]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol )
    lines( c(ROPE[2],ROPE[2]) , c(0.96*ROPEtextHt,0) , lty="dotted" , lwd=2 ,
           col=ropeCol)
    text( mean(ROPE) , ROPEtextHt ,
          bquote( .(round(100*pcInROPE))*"% in ROPE" ) ,
          adj=c(.5,0) , cex=1 , col=ropeCol )
    
    postSummary[,"ROPElow"]=ROPE[1] 
    postSummary[,"ROPEhigh"]=ROPE[2] 
    postSummary[,"pcInROPE"]=pcInROPE
  }
  # Display the HDI.
  lines( HDI , c(0,0) , lwd=4 )
  text( mean(HDI) , 0 , bquote(.(100*cred_mass) * "% HDI" ) ,
        adj=c(.5,-1.7) , cex=cex )
  text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
        adj=c(HDI_text_place,-0.5) , cex=cex )
  text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
        adj=c(1.0-HDI_text_place,-0.5) , cex=cex )
  par(xpd=F)
  #
  #return( postSummary )
  return(invisible())
}



