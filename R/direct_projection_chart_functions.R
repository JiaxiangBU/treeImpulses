# File: direct_projection_chart_functions.R
# Author: Tyler Pike
# Date: 2/5/2020
# Note(s): Houses functions to chart different types of Jorda (2005) style direct projections

#------------------------------------------------------------------
# function to plot hooks at top of charts
#------------------------------------------------------------------
#' plot.irf
#'
#' Helper function to draw border onto IRF chart
#'
#' @export

plotHookBox = function(){
    box(bty = "u")
    lines(grconvertX(c(0.018, 0, 0, 1, 1, 0.982), "npc", "user"),
          grconvertY(c(1, 1, 0, 0, 1, 1), "npc", "user"), col = 'black',
          lwd = 1.6)
}

#------------------------------------------------------------------
# function to create simple plot of localProjectionIRF output
#------------------------------------------------------------------
#' plot.irf
#'
#' Plot a single impulse response
#'
#' @param Data An object of class IRF to be charted
#'
#' @export

plot.irf = function(Data){

  # set charting option
  Title = paste0('Response of ', Data$target, ' to ', Data$shock)
  Ylim = c( min(Data$impulseResponse$lowerBound) - abs(0.1*min(Data$impulseResponse$lowerBound)),
            max(Data$impulseResponse$upperBound) + abs(0.1*max(Data$impulseResponse$upperBound)))

  # the underlying base plot
  plot(x = Data$impulseResponse$Horizon,
       y = Data$impulseResponse$Coef,
       type = 'l',
       main = Title,
       xlab = "",
       ylab = "",
       col  = "black",
       bty  = 'l',
       yaxs = "i",
       xaxs = "i",
       xlim = c(.8, (nrow(Data$impulseResponse) + .2)),
       ylim = Ylim,
       pch  = 16
  )

  # confidence interval bands
  polygon(c(Data$impulseResponse$Horizon,rev(Data$impulseResponse$Horizon)),
          c(Data$impulseResponse$upperBound,rev(Data$impulseResponse$lowerBound)),
          col="lightblue", border=NA)

  # the overlying base plot
  par(new=TRUE)
  plot(x = Data$impulseResponse$Horizon,
       y = Data$impulseResponse$Coef,
       type = 'l',
       main = Title,
       xlab = "",
       ylab = "",
       col  = "black",
       bty  = 'l',
       yaxs = "i",
       xaxs = "i",
       lwd = 2,
       xlim = c(.8, (nrow(Data$impulseResponse) + .2)),
       ylim = Ylim,
       pch  = 16
  )

  # plot the nice looking hook box
  plotHookBox()

  # labels
  abline(h=0, lwd = 1, lty = 2)

}

#------------------------------------------------------------------
# function to plot localProjectionIRF_var output
#   a wrapper for plot.irf
#------------------------------------------------------------------
#' plot.irf_var
#'
#' Plot a collection of impulse responses
#'
#' @param Data An object of class IRF_VAR to be charted
#' @param state State label used to define groups in state-dependent IRF_var
#'
#' @export

plot.irf_var = function(Data, state = NULL){

  # number of shocks
  shocks = Data %>%
    purrr::map(
      function(X){
        return(X$shock)
      }
    ) %>%
    purrr::reduce(rbind) %>%
    unique() %>%
    length()

  # number of targets
  targets = Data %>%
    purrr::map(
      function(X){
        return(X$target)
      }
    ) %>%
    purrr::reduce(rbind) %>%
    unique() %>%
    length()

  # set charting layout
  par(mfrow = c(targets, shocks))

  # chart
  for(data in Data){

    # chart specific state
    if(!is.null(state)){
      data$impulseResponse = as.data.frame(data$impulseResponse[paste0('state_',state)])
      colnames(data$impulseResponse) = gsub("^.*\\.", '', colnames(data$impulseResponse))
    }

    plot(data)
  }

}

