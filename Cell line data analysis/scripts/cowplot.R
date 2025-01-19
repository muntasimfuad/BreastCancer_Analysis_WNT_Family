# Load cowplot
library(cowplot)

# Combine plots
figure_ccl <- plot_grid(wnt2_plot, wnt7b_plot, wnt11_plot,
                        align = "v",
                        ncol = 1,label_fontfamily = "Open Sans", 
                        label_fontface = "plain",
                        labels = "AUTO",
                        label_size = 40,
                        label_x = -0.01,
                        label_y = 1.02)
# Export the plot
save_plot(filename = "figures/CCL.png",
       plot = figure_ccl,
       ncol = 1,
       base_height = 18,
       base_width = 22)
