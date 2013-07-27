plotPost = function( paramSampleVec , credMass=0.95 , compVal=NULL ,
           HDItextPlace=0.7 , ROPE=NULL , yaxt=NULL , ylab=NULL ,
           xlab=NULL , cex.lab=NULL , cex=NULL , xlim=NULL , main=NULL ,
           col=NULL , border=NULL , showMode=F , showCurve=F , breaks=NULL , 
                     ... ) {
    # Override defaults of hist function, if not specified by user:
    # (additional arguments "..." are passed to the hist function)
    if ( is.null(xlab) ) xlab="Parameter"
    if ( is.null(cex.lab) ) cex.lab=1.5
    if ( is.null(cex) ) cex=1.4
    if ( is.null(xlim) ) xlim=range( c( compVal , paramSampleVec ) )
    if ( is.null(main) ) main=""
    if ( is.null(yaxt) ) yaxt="n"
    if ( is.null(ylab) ) ylab=""
    if ( is.null(col) ) col="skyblue"
    if ( is.null(border) ) border="white"

    postSummary = matrix( NA , nrow=1 , ncol=11 , 
                          dimnames=list( c( xlab ) , 
                                         c("mean","median","mode",
                                           "hdiMass","hdiLow","hdiHigh",
                                           "compVal","pcGTcompVal",
                                           "ROPElow","ROPEhigh","pcInROPE")))              
    postSummary[,"mean"] = mean(paramSampleVec)
    postSummary[,"median"] = median(paramSampleVec)
    mcmcDensity = density(paramSampleVec)
    postSummary[,"mode"] = mcmcDensity$x[which.max(mcmcDensity$y)]

    HDI = hdi_of_mcmc( paramSampleVec , credMass )
    postSummary[,"hdiMass"]=credMass
    postSummary[,"hdiLow"]=HDI[1]
    postSummary[,"hdiHigh"]=HDI[2]

    # Plot histogram.
    if ( is.null(breaks) ) {
    breaks = c( seq( from=min(paramSampleVec) , to=max(paramSampleVec) ,
                   by=(HDI[2]-HDI[1])/18 ) , max(paramSampleVec) )
    }
    if ( !showCurve ) {
      par(xpd=NA)
      histinfo = hist( paramSampleVec , xlab=xlab , yaxt=yaxt , ylab=ylab ,
                       freq=F , border=border , col=col ,
                       xlim=xlim , main=main , cex=cex , cex.lab=cex.lab ,
                       breaks=breaks , ... )
    }
    if ( showCurve ) {
      par(xpd=NA)
      histinfo = hist( paramSampleVec , plot=F )
      densCurve = density( paramSampleVec , adjust=2 )
      plot( densCurve$x , densCurve$y , type="l" , lwd=5 , col=col , bty="n" ,
            xlim=xlim , xlab=xlab , yaxt=yaxt , ylab=ylab ,
            main=main , cex=cex , cex.lab=cex.lab , ... )
    }
    cenTendHt = 0.9*max(histinfo$density)
    cvHt = 0.7*max(histinfo$density)
    ROPEtextHt = 0.55*max(histinfo$density)
    # Display mean or mode:
    if ( showMode==F ) {
        meanParam = mean( paramSampleVec )
        text( meanParam , cenTendHt ,
              bquote(mean==.(signif(meanParam,3))) , adj=c(.5,0) , cex=cex )
    } else {
        dres = density( paramSampleVec )
        modeParam = dres$x[which.max(dres$y)]
        text( modeParam , cenTendHt ,
              bquote(mode==.(signif(modeParam,3))) , adj=c(.5,0) , cex=cex )
    }
    # Display the comparison value.
    if ( !is.null( compVal ) ) {
      cvCol = "darkgreen"
      pcgtCompVal = round( 100 * sum( paramSampleVec > compVal )
                            / length( paramSampleVec )  , 1 )
       pcltCompVal = 100 - pcgtCompVal
       lines( c(compVal,compVal) , c(0.96*cvHt,0) ,
              lty="dotted" , lwd=1 , col=cvCol )
       text( compVal , cvHt ,
             bquote( .(pcltCompVal)*"% < " *
                     .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ) ,
             adj=c(pcltCompVal/100,0) , cex=0.8*cex , col=cvCol )
      postSummary[,"compVal"] = compVal
      postSummary[,"pcGTcompVal"] = ( sum( paramSampleVec > compVal ) 
                                  / length( paramSampleVec ) )
    }
    # Display the ROPE.
    if ( !is.null( ROPE ) ) {
      ropeCol = "darkred"
       pcInROPE = ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                            / length( paramSampleVec ) )
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
    text( mean(HDI) , 0 , bquote(.(100*credMass) * "% HDI" ) ,
          adj=c(.5,-1.7) , cex=cex )
    text( HDI[1] , 0 , bquote(.(signif(HDI[1],3))) ,
          adj=c(HDItextPlace,-0.5) , cex=cex )
    text( HDI[2] , 0 , bquote(.(signif(HDI[2],3))) ,
          adj=c(1.0-HDItextPlace,-0.5) , cex=cex )
    par(xpd=F)
    #
    return( postSummary )
}
