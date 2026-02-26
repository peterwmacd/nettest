#' Plot Time Series with Shewhart-style Control Bands
#'
#' Computes summary statistics and produces plots for Shewhart-style monitoring of time series derived from sequences
#' of networks. Plots for each time series visualize a sample mean plus/minus 3 sample standard deviations, calibrated
#' from a baseline window. Additional options to add known changepoints, and save plot output to PDF.
#'
#' @param series A matrix, each column contains a univariate time series of length \code{m}.
#' @param baseline An integer, sample size of the leading portion used to estimate control limits. Defaults to the minimum of \code{ceiling(m/2)} and \code{25}.
#' @param makeplot A Boolean, should a plot be produced. Defaults to \code{TRUE}.
#' @param makeCP A numeric vector of changepoints to mark on plots. Defaults to \code{NULL}.
#' @param pdfname A string, file location and name to save as a pdf plot. Defaults to \code{NULL}.
#'
#' @return A dataframe containing baseline sample means and standard deviations for each time series.
#' If \code{makeplot=TRUE}, also produces a plot of Shewhart control charts.
#'
#' @export
CP_shewhart_plots <- function(series,
                              baseline=NULL,
                              makeplot=TRUE,
                              makeCP=NULL,
                              pdfname=NULL) {
  # dimensions
  m <- nrow(series)
  nplot <- ncol(series)
  # update colnames if not provided
  if(is.null(colnames(series))){
    colnames(series) <- paste0('series',seq_len(nplot))
  }
  # update baseline sample size if not provided
  if(is.null(baseline)){
    baseline <- min(ceiling(m/2),25)
  }
  # additional plot setup
  if(makeplot & !is.null(pdfname)){
    pdf(file=pdfname,width=11,height=8)
  }
  nr <- ceiling(sqrt(nplot)); nc <- ceiling(nplot/nr)
  par(mfrow = c(nr,nc),mar = c(4, 4, 2, 1))
  pwmcol <- c('dodgerblue3','sienna2','plum4','darkseagreen3','lightsteelblue2','lightgoldenrod2')
  # allocate space for output
  shew_matrix <- matrix(NA,nplot,2); colnames(shew_matrix) <- c('mean','sd')

  for (jj in seq_len(nplot)) {
    y <- series[,jj]
    yname <- colnames(series)[jj]

    y0 <- y[1:baseline]
    m0 <- shew_matrix[jj,1] <- mean(y0)
    s0 <- shew_matrix[jj,2] <- sd(y0)

    UCL <- m0 + 3 * s0
    LCL <- m0 - 3 * s0
    yrange <- range(c(y, UCL, LCL), na.rm = TRUE)
    yexpand <- yrange + c(-1,1)*(diff(yrange)/10) # default expansion

    if(makeplot){
      plot(1:m,1:m,type='n',ylim = yexpand,
           ylab = yname, xlab = "Index",
           main = '')
      usr <- par("usr")
      rect(0, usr[3], baseline, usr[4],
           col = adjustcolor('gray',0.2),
           border = NA)
      lines(1:m, y, type = "b",col=pwmcol[2])
      abline(h = m0, col = pwmcol[1], lty=2)
      abline(h = c(UCL, LCL), col = pwmcol[3], lty =3)
      if(!is.null(makeCP)){
        abline(v=makeCP,lty=1,col=pwmcol[4])
      }
    }
  }
  # revert plot options
  par(mfrow=c(1,1))
  if(makeplot & !is.null(pdfname)){
    dev.off()
  }
  # return Shewhart statistics
  return(shew_matrix)
}
