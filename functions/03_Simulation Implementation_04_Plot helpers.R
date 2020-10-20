
###'######################################################################
###'
###' Category: Define functions
###' 
###' Task: Define plot helper functions
###'
###' Date: 2020-03-27
###' 
###' Author: JoonHo Lee (joonho@berkeley.edu)
###' 
###'

###'######################################################################
###'
###' Load necessary packages
###'
###'

library(ggplot2)
library(cowplot)



###'######################################################################
###'
###' Common settings for theme and temporary labels
###'
###'

### Theme settings
theme_trend <- 
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())


theme_preset <- theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())


### Temporary labels
temp_labels <- labs(title = "Enter title here", 
                    subtitle = "Enter subtitle here", 
                    caption = "Enter caption here", 
                    y = "Enter ylabel here",  
                    x = "Enter xlabel here")


### Define manual palettes
color_palette <- c("firebrick1", "dodgerblue1", "forestgreen", "darkorchid1", "darkgoldenrod1", 
                   "blue", "green", "purple", "gold", "red")    

shape_palette <- c(16, 17, 15, 18, 1, 2, 0, 5, 6, 4, 3, 8, 10, 7, 9) 



###'######################################################################
###'
###' plot_compare_true_obs()
###' 
###' - Generate a plot comparing true and observed distribution of theta
###'
###'

### Define a function to generate plot
plot_compare_true_obs <- function(df_G){
  
  # Compare the true distribution (theta) to the observed data (Y)
  df_plot <- df_G %>%
    select(theta, Y) %>% 
    rename(True = theta, Observed = Y) %>%
    gather(key = variable, value = value)
  
  p1 <- ggplot(data = df_plot, aes(x = value, group = variable, fill = variable)) +
    geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
    geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
    labs(title = bquote("True vs. Observed " ~ theta ~" Distributions"), 
         x = expression(theta)) + theme_preset
  
  # Compare point estimates by sigma_k
  p2 <- ggplot(data = df_G %>% mutate(gap = abs(Y - theta)), aes(x = sigma2, y = gap)) + 
    geom_point(aes(size = sigma2)) + 
    labs(title = bquote("Gap between True and Observed vs. " ~ sigma^2), 
         y = bquote("Gap between True and Observed " ~ theta), 
         x = expression(sigma^2)) + theme_preset
  
  # Combine the two plots
  p_grid <- plot_grid(p1, p2, labels = "AUTO")
  return(p_grid)
}



###'######################################################################
###'
###' plot_lines_grp()
###'
###' - Generate a plot comparing groups
###'
###'

plot_lines_grp <- function(dataframe, 
                           x, 
                           y,
                           group, 
                           yline = NULL, 
                           ylim = NULL, 
                           sprintf = "%.2f", 
                           hjust = 0.5, vjust = 2.0){
  
  ###' Enquote x, y, and group variables
  ###' Renamed variables because scale::comma() didn't work with !!yvar
  xvar <- enquo(x)
  yvar <- enquo(y)
  groupvar <- enquo(group)
  dataframe <- dataframe %>%
    rename(xvar = !!xvar, yvar = !!yvar, groupvar = !!groupvar)
  
  ### Assign data and aesthetic mappings
  p <- ggplot(dataframe) + 
    aes(x = xvar, y = yvar, group = groupvar)
  
  ###' Add point, path, and value label layers
  p <- p + geom_point(aes(shape = groupvar, color = groupvar), size = 3.0) + 
    geom_path(aes(linetype = groupvar, color = groupvar), size = 1.0) + 
    geom_text(aes(label = sprintf(sprintf, yvar)), size = 3, hjust = hjust, vjust = vjust)
  
  ### Add vertical line layer
  if (!is.null(yline)){
    p <- p + geom_vline(aes(xintercept = yline), color = "red", linetype = "dashed")
  }
  
  # ### Scales
  # p <- p + 
  #   scale_x_continuous(breaks = seq(min(dataframe$xvar), max(dataframe$xvar), 
  #                                   by = xinterval)) + 
  #   scale_y_continuous(labels = comma, limits = ylim)
  
  ### Themes, temporary labels, and manual colors
  p + theme_trend + temp_labels + 
    scale_color_manual(values = rev(color_palette[seq(unique(dataframe$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(dataframe$groupvar))]))
}



###'######################################################################
###'
###' four_dim_plot()
###' 
###' - Define a function to automate faceted plotting
###'
###'

four_dim_plot <- function(df_plot, 
                          y, 
                          x, 
                          group, 
                          facet_row, 
                          facet_col, 
                          set_scale = "free_y", 
                          yline = NULL, 
                          ylim = NULL, 
                          sprintf = "%.2f", 
                          hjust = 0.5, vjust = 2.0, 
                          set_expand = c(0.1, 0.1)){
  
  ### Enquote and rename x, y, group, facet_row, and facet_col variables
  yvar <- enquo(y)
  xvar <- enquo(x)
  groupvar <- enquo(group)
  facet_row_var <- enquo(facet_row)
  facet_col_var <- enquo(facet_col)
  
  df_plot <- df_plot %>%
    rename(yvar = !!yvar, xvar = !!xvar, groupvar = !!groupvar, 
           facet_row_var = !!facet_row_var, facet_col_var = !! facet_col_var)
  
  
  ### Assign data and aesthetic mappings
  p <- ggplot(df_plot) + 
    aes(x = xvar, y = yvar, group = groupvar)
  
  
  ### Add point, path, and value label layers
  p <- p + geom_point(aes(shape = groupvar, color = groupvar), size = 3.0) + 
    geom_path(aes(linetype = groupvar, color = groupvar), size = 1.0) + 
    geom_text(aes(label = sprintf(sprintf, yvar)), size = 3, hjust = hjust, vjust = vjust)
  
  ### Add vertical line layer
  if (!is.null(yline)){
    p <- p + geom_vline(aes(xintercept = yline), color = "red", linetype = "dashed")
  }
  
  
  ### Faceting
  p <- p + facet_grid(facet_row_var ~ facet_col_var, scales = set_scale)
  
  
  
  ### Themes, temporary labels, and manual colors
  p + theme_trend + temp_labels + 
    scale_y_continuous( expand = set_expand ) + 
    scale_color_manual(values = rev(color_palette[seq(unique(df_plot$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(df_plot$groupvar))]))
}



###'######################################################################
###'
###' lolipop_3dim()
###' 
###' - Define a function to plot lolipop chart
###'
###'

lolipop_3dim <- function(df_plot, 
                         x, 
                         y, 
                         group, 
                         facet_row, 
                         facet_col,
                         set_scale = "fixed", 
                         yline = 0, 
                         ylim = NULL, 
                         sprintf = "%.0f", 
                         set_expand = c(0.2, 0.2)){
  
  
  ### Enquote and rename x, y, group, facet_row, and facet_col variables
  xvar <- enquo(x)
  yvar <- enquo(y)
  groupvar <- enquo(group)
  facet_row_var <- enquo(facet_row)
  facet_col_var <- enquo(facet_col)
  
  df_plot <- df_plot %>%
    rename(yvar = !!yvar, xvar = !!xvar, groupvar = !!groupvar, 
           facet_row_var = !!facet_row_var, 
           facet_col_var = !! facet_col_var)
  
  
  ### Assign data and aesthetic mappings
  p <- ggplot(df_plot) + 
    aes(x = xvar, y = yvar, group = groupvar)
  
  
  ### Add segment, point, and value label layers
  p <- p + geom_segment(aes(x = xvar, xend = xvar, 
                            y = 0, yend = yvar, 
                            color = groupvar)) + 
    geom_point(aes(x = xvar, y = yvar, color = groupvar), 
               size = 3) + 
    geom_text(aes(label = round(yvar, 0), vjust = ifelse(yvar >= 0, -1.0, 1.5)), 
              color = "black")
  
  ### Add horizontal line layer
  if (!is.null(yline)){
    p <- p + geom_hline(aes(yintercept = yline), color = "black", linetype = "solid")
  }
    
  ### Faceting
  p <- p + facet_grid(facet_row_var ~ facet_col_var + groupvar, 
                      scales = set_scale)
  
  
  ### Themes, temporary labels, and manual colors
  p + theme_trend + temp_labels + 
    scale_y_continuous( expand = set_expand ) + 
    scale_color_manual(values = rev(color_palette[seq(unique(df_plot$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(df_plot$groupvar))]))
}


