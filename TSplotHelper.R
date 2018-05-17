# Simple time series plots of data and model-based forecasts
# (uses R base graphics)
#
# Allows for up to 2 forecasts, with CIs or 2*se intervals
# and allows plotting of "factual" data over the forecast period
#
# Christopher Adolph
# 28 July 2014
#
ctrfactTS <- function(observed,
                      predicted,
                      lower=NULL,
                      upper=NULL,
                      se=NULL,
                      predicted0=NULL,
                      lower0=NULL,
                      upper0=NULL,
                      se0=NULL,
                      factual=NULL,
                      at.xaxis=NULL,
                      lab.xaxis=NULL,
                      ylim=NULL,
                      main="Counterfactual time series",
                      xlab="x",
                      ylab="y",
                      col="red",
                      col0="blue"
                      ) {


    ndata <- length(observed)
    ncft <- max(c(length(predicted), length(predicted0), length(factual)))

    if (!is.null(se)) {
        lower <- predicted - 2*se
        upper <- predicted + 2*se
    }
    if (!is.null(se0)) {
        lower0 <- predicted0 - 2*se0
        upper0 <- predicted0 + 2*se0
    }
    
    if (is.null(ylim)) {
        ydata <- c(observed, predicted, lower, upper, predicted0, lower0, upper0)
        ylim <- c(min(ydata), max(ydata))
    }

    plot.new()
    par(usr = c(0, ndata+ncft, ylim[1], ylim[2]) )

    # make the x-axis
    if (!is.null(at.xaxis)) {
        if (!is.null(lab.xaxis)) {
            axis(1, at=at.xaxis, labels=lab.xaxis)
        } else {
            axis(1, at=at.xaxis)
        }
    } else {
        axis(1)
    }

    # Make the y-axis
    axis(2)
                      
    # Make titles
    title(main=main, xlab=xlab, ylab=ylab)

    # Polygon of confidence interval
    if (!is.null(lower)) {
        x0 <- (ndata+1):(ndata + length(lower))
        polygon(x = c(ndata, x0, rev(x0), ndata),
                y = c(observed[ndata], lower, rev(upper), observed[ndata]),
                border=NA,
                col=lighten(col)
                )
    }

    # Plot the actual data (including factual data if present)
    if (!is.null(factual)) {
        lines(x = 1:(ndata+length(factual)),
              y = c(observed, factual))
    } else {
        lines(x = 1:ndata,
              y = observed
              )
    }

    # Add the predictions
    lines(x = ndata:(ndata+length(predicted)),
          y = c(observed[length(observed)],predicted),  # link up the actual data to the prediction
          col = col
          )

    # Add the lower predictive interval for alt scenario
    if (!is.null(lower0)) {
        lines(x = ndata:(ndata+length(predicted0)),
              y = c(observed[length(observed)],lower0),
              col = col0,
              lty = "dashed"
              )
    }
    
    # Add the upper predictive interval for no law
    if (!is.null(upper0)) {
        lines(x = ndata:(ndata+length(predicted0)),
              y = c(observed[length(observed)],upper0),
              col = col0,
              lty = "dashed"
              )
    }
    
    # Add the predictions for keeping law
    if (!is.null(predicted0)) {
        lines(x = ndata:(ndata+length(predicted0)),
              y = c(observed[length(observed)],predicted0), 
              col = col0
              )
    }

}




lighten <- function (col, pct = 0.75, alpha = 1) 
{
    if (abs(pct) > 1) {
        print("Warning:  Error in Lighten; invalid pct")
        pcol <- col2rgb(col)/255
    }
    else {
        col <- col2rgb(col)/255
        if (pct > 0) {
            pcol <- col + pct * (1 - col)
        }
        else {
            pcol <- col * pct
        }
    }
    pcol <- rgb(pcol[1], pcol[2], pcol[3], alpha)
    pcol
}
