## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  warning = FALSE,
  message = FALSE,
  fig.align = "center",
  out.width = "90%"
)

performance <- data.frame(
  scenario = c(
    "Moderate tree",
    "More histories",
    "1,000-tip tree",
    "16 states / 74 parameters",
    "500 tips / 8 states"
  ),
  tips = c(96, 96, 1000, 250, 500),
  histories = c(25, 100, 25, 25, 50),
  states = c(4, 4, 4, 16, 8),
  parameters = c(8, 8, 8, 74, 22),
  previous_seconds = c(0.335, 1.283, 4.713, 1.908, 5.513),
  updated_seconds = c(0.131, 0.463, 1.372, 0.365, 1.267),
  stringsAsFactors = FALSE
)
performance$speedup <- performance$previous_seconds /
  performance$updated_seconds

performance_table <- data.frame(
  Scenario = performance$scenario,
  Tips = performance$tips,
  Histories = performance$histories,
  States = performance$states,
  Parameters = performance$parameters,
  `Previous version (s)` = sprintf("%.3f", performance$previous_seconds),
  `Updated version (s)` = sprintf("%.3f", performance$updated_seconds),
  Speedup = sprintf("%.2fx", performance$speedup),
  check.names = FALSE
)

## ----performance-table--------------------------------------------------------
knitr::kable(
  performance_table,
  align = c("l", rep("r", 7)),
  caption = "Median elapsed time across repeated runs. Lower times are better."
)

## ----performance-plot, fig.width=8, fig.height=4.8, fig.cap="Speedup of the updated implementation over the previous version."----
bar_colors <- ifelse(performance$speedup >= 4, "#2b8cbe", "#7bccc4")
old_par <- par(mar = c(8, 4.2, 1, 0.5))
bars <- barplot(
  performance$speedup,
  names.arg = performance$scenario,
  las = 2,
  ylim = c(0, max(performance$speedup) * 1.18),
  ylab = "Speedup (times faster)",
  col = bar_colors,
  border = NA
)
abline(h = 1, lty = 2, col = "grey50")
text(
  bars,
  performance$speedup,
  labels = sprintf("%.2fx", performance$speedup),
  pos = 3,
  cex = 0.9
)
par(old_par)

