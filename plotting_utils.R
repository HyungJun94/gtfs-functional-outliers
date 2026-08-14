#' Plot Generated Curves and Spot Outliers
#'
#' @param Y Matrix. Functional data matrix where rows represent curves.
#' @param main Character. Title of the plot.
#' @param n Integer. Number of observations.
#' @param c Numeric. Contamination ratio.
#' @export
plot_mod <- function(Y, main='', n=200, c=0.025) {
  t_seq <- seq(0, 1, length.out = ncol(Y))
  
  # Spot outliers in red
  plot(t_seq, Y[1,], type='l', ylim=c(min(Y), max(Y)), main=main, ylab='f(t)', xlab='t')
  for(i in 1:n) {
    lines(t_seq, Y[i,], col='grey')
  }
  for(i in (n*(1-c)+1):n) {
    lines(t_seq, Y[i,], col='red')
  }
  lines(t_seq, Y[1,], col='blue')
}

#' Plot 3D Surface using Plotly
#'
#' @param x Matrix. Data matrix to visualize as a 3D surface.
#' @export
splot <- function(x) {
  plotly::plot_ly(z = ~x) %>% plotly::add_surface()
}

#' Plot Heatmap using Plotly
#'
#' @param x Matrix. Data matrix to visualize as a heatmap.
#' @export
hplot <- function(x) {
  plotly::plot_ly(z = ~x) %>% plotly::add_heatmap()
}