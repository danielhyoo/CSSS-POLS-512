# Dynamic Panel Data Models (Arellano-Bond, etc.)
# for fixed effects in short-T panels
# Model Fitting and Interpretation via Simulation
#
# Helper functions 
#
# Christopher Adolph   faculty.washington.edu/cadolph
# 13 May 2015
             
# Graph expected values, baseline, and first diffs for all 8 models using tile...
cigLineplots <- function(EV, EVbase, FD, RR, limitsEV, limitsFD, limitsRR,
                         initialY, main, file) {

# Get 3 nice colors for traces
require(RColorBrewer)
col <- brewer.pal(4,"Dark2")
periods <- 1:3

# Set up lineplot trace of EV packpc
EVtrace <- lineplot(x=c(0, periods),
                  y=c(initialY, EV$pe),
                  lower=c(initialY, EV$lower),
                  upper=c(initialY, EV$upper),
                  col=col[1],
                  plot=1)

EVptrace <- pointsTile(x=c(0, periods),
                     y=c(initialY, EV$pe),
                     pch=16,
                     col=col[1],
                     plot=1)

BASEtrace <- lineplot(x=c(0, periods),
                    y=c(initialY, EVbase$pe),
                    lower=c(initialY, EVbase$lower),
                    upper=c(initialY, EVbase$upper),
                    ci=list(mark="dashed"),
                    col=col[2],
                    plot=1)

BASEptrace <- pointsTile(x=c(0, periods),
                       y=c(initialY, EVbase$pe),
                       pch=16,
                       col=col[2],
                       plot=1)

FDtrace <- lineplot(x=c(0, periods),
                  y=c(0, FD$pe),
                  lower=c(0, FD$lower),
                  upper=c(0, FD$upper),
                  col=col[3],
                  plot=2)

FDptrace <- pointsTile(x=c(0, periods),
                     y=c(initialY, FD$pe),
                     pch=16,
                     col=col[3],
                     plot=2)

RRtrace <- lineplot(x=c(0, periods),
                  y=c(0, 100*(RR$pe-1)),
                  lower=c(0, 100*(RR$lower-1)),
                  upper=c(0, 100*(RR$upper-1)),
                  col=col[4],
                  plot=3)

RRptrace <- pointsTile(x=c(0, periods),
                     y=c(initialY, 100*(RR$pe-1)),
                     pch=16,
                     col=col[4],
                     plot=3)

# Set up text labels of lines
text1 <- textTile(x=3.9,
                  y=EV$pe[length(EV$pe)],
                  labels="+$.60/pack",
                  col=col[1],
                  cex=0.75,
                  plot=1)

text2 <- textTile(x=3.85,
                  y=EVbase$pe[length(EVbase$pe)],
                  labels="baseline",
                  col=col[2],
                  cex=0.75,
                  plot=1)

text3 <- textTile(x=3.85,
                  y=FD$pe[length(FD$pe)],
                  labels="change\nin packs",
                  col=col[3],
                  cex=0.75,
                  plot=2)

text4 <- textTile(x=3.85,
                  y=100*(RR$pe[length(RR$pe)] -1), 
                  labels="percent\nchange",
                  col=col[4],
                  cex=0.75,
                  plot=3)

# Set up baseline: for first difference, this is 0
baselineFD <- linesTile(x=c(-1,5),
                      y=c(0,0),
                      plot=2)

# Set up baseline: for percent change, this is 0
baselineRR <- linesTile(x=c(-1,5),
                      y=c(0, 0),
                      plot=3)

# Plot all traces using tile
tile(EVtrace, BASEtrace, FDtrace, RRtrace,
     EVptrace, BASEptrace, FDptrace, RRptrace,
     text1, text2, text3, text4,
     baselineFD, baselineRR,
     limits=rbind(c(-0.25,4.8,limitsEV),
                  c(-0.25,4.8,limitsFD),
                  c(-0.25,4.8,limitsRR)),
     xaxis=list(at=c(0,periods), labels=c("t", "t+1", "t+2", "t+3")),
     yaxis=list(label.loc=-0.5, major=FALSE),
     xaxistitle=list(labels="Forecast Year"),
     yaxistitle=list(labels1="Expected Packs per capita",
                     labels2="Change, Packs per capita",
                     labels3="%Change, Packs per capita"),
     maintitle=list(labels=main),
     gridlines=list(type="y"),
     height=list(maintitle=4),
     width=list(null=5,yaxistitle=4,yaxis.labelspace=-0.5),
     output=list(file=file,width=10)
     )

}


