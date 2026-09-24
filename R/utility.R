# ---- Themes ----
#' @title Set modern default ggplot theme.
#' @description
#' Applies a custom modern ggplot theme. This will unify themes across track plots, but will also apply to all ggplots in the R session. To reset the theme, run `set_theme()`.
#' @import ggplot2
#' @export
set_theme_modern <- function() {
  set_theme(
    theme_classic()+
      theme(
        panel.grid = element_blank(),
        axis.text = element_text(size = 12, color = "black"), axis.title = element_text(size = 12),
        legend.text = element_text(size = 12, color = "black"), legend.title = element_text(size = 12),
        axis.line = element_line(linewidth = 0.4), panel.border = element_rect(linewidth = 0.4, color = "black"),
        strip.background = element_blank(), strip.text = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")
      )
  )
}
