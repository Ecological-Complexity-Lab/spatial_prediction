# ---- useful for plotting ----
## ---- combine plots ----
combine_plots <- function(p1, p2,
                          bottom_label = "Balanced accuracy",
                          left_label = "Count",
                          plot_margin = c(0.5, 0.5, 1, 0.3),
                          label_fontsize = 14,
                          label_fontface = "bold",
                          widths_subplots = c(1, 1),
                          final_widths = c(2, 0.3)) {
  
  # Load required packages
  require(ggplot2)
  require(gridExtra)
  require(grid)
  
  # Adjust individual plots
  p1_mod <- p1 +
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.margin = unit(plot_margin, "cm"))
  
  p2_mod <- p2 +
    theme(legend.position = "none",
          axis.title = element_blank(),
          plot.margin = unit(plot_margin, "cm"))
  
  # Arrange the two plots side-by-side
  combined_plots <- arrangeGrob(p1_mod, p2_mod, 
                                ncol = 2, 
                                widths = widths_subplots)
  
  # Add axis labels using arrangeGrob (the bottom and left text grobs)
  combined_with_axes <- arrangeGrob(
    combined_plots,
    bottom = textGrob(bottom_label, 
                      gp = gpar(fontsize = label_fontsize, fontface = label_fontface), 
                      vjust = -1.5),
    left   = textGrob(left_label, 
                      rot = 90, 
                      gp = gpar(fontsize = label_fontsize, fontface = label_fontface))
  )
  
  # Finally, arrange the whole thing with additional spacing if needed
  final_plot <- grid.arrange(
    combined_with_axes,
    ncol = 2,
    widths = final_widths
  )
  
  return(final_plot)
}

# Example usage:
# final_plot <- combine_plots(site_distrib, isl_distrib)

