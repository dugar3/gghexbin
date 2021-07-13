#' Something about what hex.R does
#' 
#' @param x independent variable
#' @param y dependent variable
#' @param xmin minimum value of x, default value is min(x)
#' @param xmax maximum value of x, default value is max(x)
#' @param ymin minimum value of y, default value is min(y)
#' @param ymax maximum value of y, default value is max(y)
#' @param pct_of_plot_used maximum value of y, default value is max(y)
#' @param binning_relation a value between 0 and 1 that determines the relationship between count and hexagon area proportion. 
#' A value of 0.5 means that area and count have a 1:1 ratio, A value of 0 means that all hexagons are the same size, regardless of count
#' @param lattice logical that indicates whether the lattice information is returned
#' @param plotted logical that indicates whether the ggproto information is returned
#' 
#' @import data.table ggplot2
#' 
#' @return Returns lattice information, \code{lt}, that includes hexagons' circumcenter coordinates, radii, and counts; 
#' Also returns ggproto object, \code{hexplot}, which contains all of the geom_poly() generated hexagons for plotting.
#' 
#' @examples
#' library(ggplot2)
#' library(data.table)
#' 
#' data(iris)
#' x=iris$Sepal.Length
#' y=iris$Petal.Length
#' 
#' hex_object <- hex(x=x, y=y, plotted = TRUE, lattice = TRUE, binning_relation = 0.5)
#' lin_reg <- lm(y~x, data.table(x,y))
#' ggplot() + hex_object[[2]] + ggtitle("iris Data") + xlab("Sepal Length") + ylab("Petal Length") + geom_abline(slope = lin_reg$coefficients[2], intercept = lin_reg$coefficients[1]) + coord_fixed()  + guides(fill=guide_legend(title="Counts")) + scale_fill_gradient(low="green", high="red")
#' 
#' @rdname hex
#' 
#' @export hex

hex = function(x, y, xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y), pct_of_plot_used = 20, binning_relation = 0.5, lattice = TRUE, plotted = TRUE){
  # Inputs: n (legnth), x and y (data vectors), 
  # pct_of_plot_used (scaling parameter < IDIM), 
  # binning_relation (scaling of the hexagons with 
  # respect to the amount of points contained)
  # binning_relation = 2 means that the areas of
  # the hexagons are proportional to the number of
  # points that they contain
  #
  # Inputs: xmin, xmax, ymin, ymax are data extremes
  # Outputs: nl (length), xl and yl 
  # (lattice centers), counts (nonzero counts, nl <= n)
  
  
  # Binning relation is the norm-relationship between the relative size of the hexagons and the number of points they contain
  
  
  if(length(x) != length(y)){
    stop("x and y data must be the same length")
  }
  
  
  # Inputs
  size = pct_of_plot_used
  if (pct_of_plot_used < 0 & pct_of_plot_used > 100){
    stop("pct_of_plot_used must be between 0 and 100 (inclusive)")
  }
  
  
  n = length(x)
  counts = rep(0, n)
  xl = rep(0, n)
  yl = rep(0, n)
  
  #Lattice Dimensions
  IDIM = 100
  JDIM = floor(IDIM/sqrt(3)+1)
  
  
  lat1 = matrix(data = 0, nrow = IDIM, ncol = JDIM)
  lat2 = matrix(data = 0, nrow = IDIM, ncol = JDIM)
  
  # Scaling constants
  xr = xmax - xmin
  yr = ymax - ymin
  c1 = size/xr
  c2 = size/(xr*sqrt(3))
  c3 = 1/c1
  c4 = 1/c2
  
  # Binning
  for (k in 1:n){
    sx = c1*(x[k] - xmin)
    sy = c2*(y[k] - ymin)
    i1 = floor(sx + 0.5)
    j1 = floor(sy + 0.5)
    i2 = floor(sx)
    j2 = floor(sy)
    if (((sx-i1)^2 + 3*(sy-j1)^2) < ((sx-i2-0.5)^2 + 3*(sy - j2 - 0.5)^2)){
      lat1[i1, j1] = lat1[i1, j1] + 1
    }
    else{
      lat2[i2, j2] = lat2[i2, j2] + 1
      
    }
  }
  
  # Write Lattice Points and Counts
  
  imax = floor(size)
  jmax = floor(size/sqrt(3))
  nl = 0
  
  for(j in 1:jmax){
    for(i in 1:imax){
      if(lat1[i, j] > 0){
        nl = nl + 1
        xl[nl] = c3*i + xmin
        yl[nl] = c4*j + ymin
        counts[nl] = lat1[i, j]
      }
    }
    for(i in 1:imax){
      if(lat2[i, j] > 0){
        nl = nl+1
        xl[nl] = c3*(i+0.5) + xmin
        yl[nl] = c4*(j+0.5) + ymin
        counts[nl] = lat2[i, j]
      }
    }
  }
  max_circumradius = 0.5*xr/size
  counts2 = counts[counts != 0]
  
  if (binning_relation >= 0 & binning_relation <= 1){
    radii = (counts2/max(counts2))^(binning_relation)*max_circumradius
  }
  else{
    stop("binning_relation must be between 0 and 1 (inclusive)")
  }
  
  lt = list(xl = xl[counts != 0], yl = yl[counts != 0], counts = counts2, 
            radii =  radii
  )
  
  poly = c()
  
  fill <- group <- NULL # To get rid of R CMD check NOTEs - See option 2 here - https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  for (i in 1:length(lt$xl)){
    
    centerx = lt$xl[i]
    centery = lt$yl[i]
    radius = lt$radii[i]
    
    poly_hex <- data.table( x=c(centerx, centerx - sqrt(3)*radius/2, centerx - sqrt(3)*radius/2, centerx , centerx + sqrt(3)*radius/2, centerx + sqrt(3)*radius/2 ), y=c(centery + radius,  centery + 1/2*radius, centery - 1/2*radius, centery-radius , centery - 1/2*radius, centery + 1/2*radius ), 
                            fill = rep(lt$counts[i],6),
                            group = rep(i, 6) # This separates the hexagons
    )
    
    poly = rbind(poly, poly_hex)
    
  }
  # Setting color equal to fill will give the hexagons the same color for outline and fill
  hexplot <- list(geom_polygon(data=poly, aes(x=x, y=y, fill = fill, group = group)))
  
  # Returns
  
  if(lattice & plotted){
    return(list(lt, hexplot))
  }
  else if (lattice){
    return(list(lt, NULL))
  }
  else if (plotted){
    return(list(NULL, hexplot))
  }
  
}
