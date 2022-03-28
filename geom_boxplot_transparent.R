# Gist for question http://stackoverflow.com/questions/38305785/ggplot-adding-alpha-value-to-a-whole-layer
# Works for R 3.3.1, ggplot2_2.1.0

library(ggplot2)
`%||%` <- ggplot2:::`%||%`
ggname <- ggplot2:::ggname

GeomBoxplot <- ggproto("GeomBoxplot", Geom,
                       setup_data = function(data, params) {
                         data$width <- data$width %||%
                           params$width %||% (resolution(data$x, FALSE) * 0.9)
                         
                         if (!is.null(data$outliers)) {
                           suppressWarnings({
                             out_min <- vapply(data$outliers, min, numeric(1))
                             out_max <- vapply(data$outliers, max, numeric(1))
                           })
                           
                           data$ymin_final <- pmin(out_min, data$ymin)
                           data$ymax_final <- pmax(out_max, data$ymax)
                         }
                         
                         # if `varwidth` not requested or not available, don't use it
                         if (is.null(params) || is.null(params$varwidth) || !params$varwidth || is.null(data$relvarwidth)) {
                           data$xmin <- data$x - data$width / 2
                           data$xmax <- data$x + data$width / 2
                         } else {
                           # make `relvarwidth` relative to the size of the largest group
                           data$relvarwidth <- data$relvarwidth / max(data$relvarwidth)
                           data$xmin <- data$x - data$relvarwidth * data$width / 2
                           data$xmax <- data$x + data$relvarwidth * data$width / 2
                         }
                         data$width <- NULL
                         if (!is.null(data$relvarwidth)) data$relvarwidth <- NULL
                         
                         data
                       },
                       
                       draw_group = function(data, panel_scales, coord, fatten = 2,
                                             outlier.colour = NULL, outlier.shape = 19,
                                             outlier.size = 1.5, outlier.stroke = 0.5,
                                             notch = FALSE, notchwidth = 0.5, varwidth = FALSE) {
                         
                         common <- data.frame(
                           colour = data$colour,
                           size = data$size,
                           linetype = data$linetype,
                           fill = alpha(data$fill, data$alpha),
                           group = data$group,
                           alpha = data$alpha,
                           stringsAsFactors = FALSE
                         )
                         
                         whiskers <- data.frame(
                           x = data$x,
                           xend = data$x,
                           y = c(data$upper, data$lower),
                           yend = c(data$ymax, data$ymin),
                           common,
                           stringsAsFactors = FALSE
                         )
                         
                         box <- data.frame(
                           xmin = data$xmin,
                           xmax = data$xmax,
                           ymin = data$lower,
                           y = data$middle,
                           ymax = data$upper,
                           ynotchlower = ifelse(notch, data$notchlower, NA),
                           ynotchupper = ifelse(notch, data$notchupper, NA),
                           notchwidth = notchwidth,
                           common,
                           stringsAsFactors = FALSE
                         )
                         
                         if (!is.null(data$outliers) && length(data$outliers[[1]] >= 1)) {
                           outliers <- data.frame(
                             y = data$outliers[[1]],
                             x = data$x[1],
                             colour = outlier.colour %||% data$colour[1],
                             shape = outlier.shape %||% data$shape[1],
                             size = outlier.size %||% data$size[1],
                             stroke = outlier.stroke %||% data$stroke[1],
                             fill = NA,
                             alpha = data$alpha,
                             stringsAsFactors = FALSE
                           )
                           outliers_grob <- GeomPoint$draw_panel(outliers, panel_scales, coord)
                         } else {
                           outliers_grob <- NULL
                         }
                         
                         ggname("geom_boxplot", grobTree(
                           outliers_grob,
                           GeomSegment$draw_panel(whiskers, panel_scales, coord),
                           GeomCrossbar$draw_panel(box, fatten = fatten, panel_scales, coord)
                         ))
                       },
                       
                       draw_key = draw_key_boxplot,
                       
                       default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                                         alpha = NA, shape = 19, linetype = "solid"),
                       
                       required_aes = c("x", "lower", "upper", "middle", "ymin", "ymax")
)

GeomCrossbar <- ggproto("GeomCrossbar", Geom,
                        setup_data = function(data, params) {
                          GeomErrorbar$setup_data(data, params)
                        },
                        
                        default_aes = aes(colour = "black", fill = NA, size = 0.5, linetype = 1,
                                          alpha = NA),
                        
                        required_aes = c("x", "y", "ymin", "ymax"),
                        
                        draw_key = draw_key_crossbar,
                        
                        draw_panel = function(data, panel_scales, coord, fatten = 2.5, width = NULL) {
                          middle <- transform(data, x = xmin, xend = xmax, yend = y, size = size * fatten, alpha = alpha)
                          
                          has_notch <- !is.null(data$ynotchlower) && !is.null(data$ynotchupper) &&
                            !is.na(data$ynotchlower) && !is.na(data$ynotchupper)
                          
                          if (has_notch) {
                            if (data$ynotchlower < data$ymin  ||  data$ynotchupper > data$ymax)
                              message("notch went outside hinges. Try setting notch=FALSE.")
                            
                            notchindent <- (1 - data$notchwidth) * (data$xmax - data$xmin) / 2
                            
                            middle$x <- middle$x + notchindent
                            middle$xend <- middle$xend - notchindent
                            
                            box <- data.frame(
                              x = c(
                                data$xmin, data$xmin, data$xmin + notchindent, data$xmin, data$xmin,
                                data$xmax, data$xmax, data$xmax - notchindent, data$xmax, data$xmax,
                                data$xmin
                              ),
                              y = c(
                                data$ymax, data$ynotchupper, data$y, data$ynotchlower, data$ymin,
                                data$ymin, data$ynotchlower, data$y, data$ynotchupper, data$ymax,
                                data$ymax
                              ),
                              alpha = data$alpha,
                              colour = data$colour,
                              size = data$size,
                              linetype = data$linetype, fill = data$fill,
                              group = seq_len(nrow(data)),
                              stringsAsFactors = FALSE
                            )
                          } else {
                            # No notch
                            box <- data.frame(
                              x = c(data$xmin, data$xmin, data$xmax, data$xmax, data$xmin),
                              y = c(data$ymax, data$ymin, data$ymin, data$ymax, data$ymax),
                              alpha = data$alpha,
                              colour = data$colour,
                              size = data$size,
                              linetype = data$linetype,
                              fill = data$fill,
                              group = seq_len(nrow(data)), # each bar forms it's own group
                              stringsAsFactors = FALSE
                            )
                          }
                          
                          ggname("geom_crossbar", gTree(children = gList(
                            GeomPolygon$draw_panel(box, panel_scales, coord),
                            GeomSegment$draw_panel(middle, panel_scales, coord)
                          )))
                        }
)

GeomPolygon <- ggproto("GeomPolygon", Geom,
                       draw_panel = function(data, panel_scales, coord) {
                         n <- nrow(data)
                         if (n == 1) return(zeroGrob())
                         
                         munched <- coord_munch(coord, data, panel_scales)
                         # Sort by group to make sure that colors, fill, etc. come in same order
                         munched <- munched[order(munched$group), ]
                         
                         # For gpar(), there is one entry per polygon (not one entry per point).
                         # We'll pull the first value from each group, and assume all these values
                         # are the same within each group.
                         first_idx <- !duplicated(munched$group)
                         first_rows <- munched[first_idx, ]
                         
                         ggname("geom_polygon",
                                polygonGrob(munched$x, munched$y, default.units = "native",
                                            id = munched$group,
                                            gp = gpar(
                                              col = alpha(first_rows$colour, first_rows$alpha),
                                              fill = alpha(first_rows$fill, first_rows$alpha),
                                              lwd = first_rows$size * .pt,
                                              lty = first_rows$linetype
                                            )
                                )
                         )
                       },
                       
                       default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1,
                                         alpha = NA),
                       
                       handle_na = function(data, params) {
                         data
                       },
                       
                       required_aes = c("x", "y"),
                       
                       draw_key = draw_key_polygon
)

assignInNamespace("GeomBoxplot", GeomBoxplot, "ggplot2")
assignInNamespace("GeomCrossbar", GeomCrossbar, "ggplot2")
assignInNamespace("GeomPolygon", GeomPolygon, "ggplot2")

ggplot(mtcars, aes(factor(cyl), mpg)) + 
  geom_boxplot(aes(fill = factor(cyl), color = factor(cyl)), alpha = 0.4, size = 1.4) 