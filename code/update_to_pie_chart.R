library(magick)
library(pdftools)

# add percentages
df_counts <- df_counts |>
  dplyr::mutate(
    sigm_cat = dplyr::recode(
      sigm_cat,
      "offs↑ only" = "External data",
      "diag↑ only" = "Local data",
      "both↑"      = "Both"
    )
  )

values <- df_counts$n

labels <- paste0(df_counts$sigm_cat,
                 "\n",
                 df_counts$n,
                 " (", df_counts$pct, ")")

pdf(
  file   = "results/paper_figs/pie_chart.pdf",
  width  = 6,    # inches
  height = 6,
  family = "Helvetica"   # or another installed font
)
# Draw pie without labels
pie(values,
    labels = NA,
    col = my_cols,
    border = "white")

# # Add labels
fractions <- values / sum(values)
cum_fractions <- cumsum(fractions)
mid_angles <- 1.6 * pi * (cum_fractions - fractions / 2)

label_radius <- rep(1.3, length(values))
i <- which(df_counts$sigm_cat == "Both")
label_radius[i] <- 1.4   # move this one further out
j <- which(df_counts$sigm_cat == "Local data")
label_radius[j] <- 0.9
k <- which(df_counts$sigm_cat == "External data")
label_radius[k] <- 0.82

text_colors <- c("salmon", "plum3", "lightsteelblue")

text(label_radius * cos(mid_angles),
     label_radius * sin(mid_angles),
     labels = labels,
     adj = ifelse(cos(mid_angles) > 0, 0, 1),
     col = text_colors,
     cex = 1.2,
     xpd = TRUE)

dev.off()     # close the file


# add to fig 3 at the end of the script
# Read the pie chart PDF as an image (first page)
img <- magick::image_read_pdf("results/paper_figs/pie_chart.pdf", density = 300)

image_info(img)

img_cropped <- image_crop(img, geometry = "1800x1800+0+280") # "WIDTHxHEIGHT+LEFT+TOP"
new_width  <- 1800 - 0 # how much to crop from right
new_height <- 1800 - 550 # how much to crop from bottom

img_cropped <- image_crop(img_cropped, 
                          geometry = paste0(new_width, "x", new_height, "+0+0"))


# Convert to grob
pie_grob <- cowplot::ggdraw() + cowplot::draw_image(img_cropped)

bottom_row <- plot_grid(
  pie_grob,
  final_plot,
  rel_widths = c(0.9, 1.2),
  labels = c("(b)", "(c)"),
  label_size = 15
)


fig3 <- plot_grid(map_missing_links + theme(plot.margin = unit(c(0.8,0.2,0.2,0.2), "cm")), 
                  bottom_row, labels = c('(a)', ''), 
                  ncol = 1, rel_heights = c(1.1, 0.7), label_size = 15)


pdf(file   = "results/paper_figs/missing_interactions_degree2.pdf",
    width  = 13,    # inches
    height = 11,
    family = "Helvetica"   # or another installed font
)
fig3
dev.off()

