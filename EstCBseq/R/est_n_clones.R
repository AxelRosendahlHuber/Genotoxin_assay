#' est_n_clones
#' @name est_n_clones
#' @description estimate the number of control and exposed clones needed for a significant result
#' based on the number of mutations in a few inital cells
#' @param ctrl_muts a numeric vector containing the number of mutations found in individual control clones.
#' At least 3 values are needed.
#' @param cond_muts similar to \code{ctrl_muts}, but for exposed clones
#' @param range_ctrl a numeric vector containing the numbers of potential control clones to be tested.
#' Default: \code{2:8}
#' @param range_cond similar to \code{range_ctrl}, but for exposed clones
#' @param n_repeat the number of iterations of tests on which the final P-value is based.
#' At least 100 is recommended. A higher number gives a more reliable result.
#' @return a heatmap plot indicating the number of control and exposed clones to sequence to determine a statistical significant increase.
#' @import ggplot2 reshape2
#' @export

est_n_clones <- function(ctrl_muts, cond_muts, range_ctrl = 2:8, range_cond = 2:8, n_repeat = 100) {
    stopifnot(is.numeric(range_ctrl),
              is.numeric(range_cond),
              is.numeric(ctrl_muts),
              is.numeric(n_repeat),
              is.numeric(cond_muts),
              all(ctrl_muts > 0),
              all(cond_muts > 0),
              length(ctrl_muts) > 2,
              length(cond_muts) > 2,
              all(range_cond > 1),
              all(range_ctrl > 1))

  set.seed(12345)
  # calculate sd and mean for control and condition values
  sd_ctrl <- sd(ctrl_muts)
  mean_ctrl <- mean(ctrl_muts)
  sd_cond <- sd(cond_muts)
  mean_cond <- mean(cond_muts)
  # make a matrix for
  avgP_value <- matrix(nrow = length(range_ctrl), ncol = length(range_cond))
  colnames(avgP_value) <- as.character(range_cond)
  rownames(avgP_value) <- as.character(range_ctrl)
  for (n_control in range_ctrl) {
    for (n_exposed in range_cond) {
      ptest <- numeric()
      for (i in seq(n_repeat)) {
        ctrl_muts_generated <- rnorm(n_control, mean = mean_ctrl, sd = sd_ctrl)
        cond_muts_generated <- rnorm(n_exposed, mean = mean_cond, sd = sd_cond)
        ptest[i] <- t.test(ctrl_muts_generated, cond_muts_generated)$p.value
      }
      avgP_value[as.character(n_control), as.character(n_exposed)] <- mean(ptest)
    }
  }

  df <- reshape2::melt(avgP_value,value_name <- "p_value")
  df$Var1 <- factor(df$Var1, levels = sort(as.numeric(unique(df$Var1))))
  df$Var2 <- factor(df$Var2, levels = sort(as.numeric(unique(df$Var2))))

  df$signif <- ifelse(df$value < 0.05, "*", "N/S")
  df$signif[df$value < 0.01] <- "**"
  df$signif[df$value < 0.01] <- "***"

  plot <- ggplot2::ggplot(df, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + geom_text(aes(label = signif), color = "white") +
    theme_minimal() + ylab("Exposed clones") + xlab("Control clones")
  return(plot)
}
