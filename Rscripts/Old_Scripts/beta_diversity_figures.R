library(cowplot)

setwd("~/projects/DenoiseCompare/Revised_figures/")


#blueberry beta diversity figure


Blueberry_grid_plot <- plot_grid(Unifrac_box, Bray_Curtis_box, unifrac_plot, bray_curtis_plot, labels="AUTO", rel_heights = c(1,2))
Blueberry_grid_plot



#biscuit beta diversity figure


#done
biscuit_grid_plot <- plot_grid(Unifrac_box, recorded_Bray_Curtis_box, noPos_unifrac_plot, bray_curtis_plot, labels="AUTO", rel_heights = c(1,2))
biscuit_grid_plot
#done


#Exercise beta diversity figure
Exercise_grid_plot <- plot_grid(Unifrac_box, Bray_Curtis_box, unifrac_plot, bray_curtis_plot, labels="AUTO", rel_heights = c(1,2))
Exercise_grid_plot



#supplemental unwieghted figures

Unweighted_grid_plot <- plot_grid(blueberry_U_Unifrac_box, biscuit_unweighted_Unifrac_box, Exercise_U_Unifrac_box,
                                  Blueberry_unweighted_unifrac, biscuit_unweighted_plot, Exercise_U_unifrac_plot, labels="AUTO", rel_heights = c(1,2))
Unweighted_grid_plot


#supplemental deblur no postive filtering figure
NoPos_grid_plot_blueberry <- plot_grid(noPos_weighted_unifrac, noPos_unweighted, noPos_Bray_curtis, nrow=3, labels="AUTO")
NoPos_grid_plot_blueberry
