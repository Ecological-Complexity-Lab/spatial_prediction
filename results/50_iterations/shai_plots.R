# for fidelity

make_facet_scatter_plot <- function(data,                                     
                                    evaluator = "f1_score",                                      
                                    pivot_cols = c("avg_sorensen_plants", "avg_sorensen_pollinators"),                                     
                                    names_to = "jaccard_type",                                     
                                    values_to = "jaccard_value",                                     
                                    x_lab = "Jaccard similarity",                                     
                                    y_lab = "F1 score",                                     
                                    plot_title = NULL,                                     
                                    facet_scales = "free_x") {      
  # Reshape data from wide to long format for the specified pivot columns   
  df_long <- data %>%      
    pivot_longer(cols = all_of(pivot_cols),                   
                 names_to = names_to,                   
                 values_to = values_to)      
  
  # For each facet (jaccard_type), compute correlation between the evaluator and jaccard_value   
  cor_table <- df_long %>%     
    group_by(!!sym(names_to)) %>%     
    summarise(       
      cor_value = cor(.data[[evaluator]], .data[[values_to]], use = "complete.obs", method = "pearson"),       
      p_value   = cor.test(.data[[evaluator]], .data[[values_to]], method = "pearson")$p.value     
    ) %>%     
    ungroup()      
  
  # Create annotations with formatted correlation coefficients and p-values   
  cor_table_annot <- cor_table %>%     
    mutate(       
      r_fmt = formatC(cor_value, format = "f", digits = 2),       
      p_fmt = ifelse(         
        p_value < 0.001,         
        formatC(p_value, format = "e", digits = 2),  # Scientific notation         
        formatC(p_value, format = "f", digits = 3)   # Otherwise       
      ),       
      label_text = paste0("r = ", r_fmt, ", p = ", p_fmt)     
    )      
  
  # Facet labels (renaming)   
  facet_labels <- c(     
    avg_sorensen_plants = "Plants",     
    avg_sorensen_pollinators = "Pollinators"   
  )      
  
  # Construct the faceted scatter plot   
  plot <- ggplot(df_long, aes_string(x = values_to, y = evaluator, color = names_to)) +  # map color by facet type
    geom_point(alpha = 0.6, size = 3) +     
    geom_smooth(method = "lm", se = FALSE, color = "steelblue") +     
    facet_wrap(as.formula(paste("~", names_to)),                 
               scales = facet_scales,                
               labeller = as_labeller(facet_labels)) +     
    scale_x_continuous(labels = number_format(accuracy = 0.02)) +     
    # Place the correlation annotation in the upper-right corner of each facet     
    geom_text(data = cor_table_annot,               
              aes(label = label_text),               
              x = Inf,               
              y = Inf,               
              hjust = 1.1,               
              vjust = 1.2,               
              size = 5,               
              color = "black") +     
    labs(x = x_lab, y = y_lab, title = plot_title) +     
    theme_minimal() +     
    # Add spacing between panels here:    
    theme(
      panel.spacing = unit(3, "lines"),  # increase spacing between facets
      strip.text = element_text(size = 20),  # facet titles bigger
      axis.title.y = element_text(size = 20),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.ticks = element_line(color = "black") + tme
    ) +
    # Set manual colors for points by facet type:
    scale_color_manual(values = c(
      avg_sorensen_plants = "palegreen3",
      avg_sorensen_pollinators = "thistle"
    ))
  
  return(plot)
}

f1_fidelity <- make_facet_scatter_plot(data = working_df_offs, 
                                       evaluator = "f1_score",
                                       pivot_cols = c("avg_sorensen_plants", "avg_sorensen_pollinators"),
                                       x_lab = NULL,
                                       y_lab = "F1 score",
                                       facet_scales = "free_x") + theme(legend.position = "none") + tme
f1_fidelity

rmse_fidelity <- make_facet_scatter_plot(data = working_df_offs, 
                                       evaluator = "rmse",
                                       pivot_cols = c("avg_sorensen_plants", "avg_sorensen_pollinators"),
                                       x_lab = "Similarity",
                                       y_lab = "RMSE",
                                       facet_scales = "free_x") + theme(legend.position = "none", strip.text = element_blank() ) + tme
rmse_fidelity

library(patchwork)

final_combined_plot <- f1_fidelity / rmse_fidelity
final_combined_plot
