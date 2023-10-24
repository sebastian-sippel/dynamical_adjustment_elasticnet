### This script is to pruduce different color maps for plotting

## 14.05.2014
## Sebastian Sippel

frenchcolormap <- function(n = 1024) {
    return(colorRampPalette(c("red", "ghostwhite", "blue"))(n))
}


frenchcolormap.fancy <- function(n = 1024) {
  return(colorRampPalette(c("orange","red", "ghostwhite", "blue", "darkgreen"))(n))
}


# frenchcolormap_veg: from redbrown to green:
frenchcolormap.vegetation <- function(n = 1024) {
  return(colorRampPalette(c("brown", "ghostwhite", "darkgreen"))(n))
}


# frenchcolormap_veg: from redbrown to darkblue:
frenchcolormap.precipitation <- function(n = 1024) {
  return(colorRampPalette(c("brown", "ghostwhite", "darkblue"))(n))
}

# make transparent color for polygons:
make.transparent.color <- function(x, alpha = 80) {  
  temp = c(col2rgb(col = x))
  return(rgb(red=temp[1], green=temp[2], blue=temp[3], alpha=alpha, maxColorValue=256))
}


frenchcolormap.greyscale <- function(n = 1024) {
  return(colorRampPalette(c("red", "gray", "blue"))(n))
}


# frenchcolormap_veg: from redbrown to green:
colormap.vegetation <- function(n = 1024) {
  return(colorRampPalette(c("brown", "darkgreen"))(n))
}


colormap.2col.transient <- function(n = 1024, col1, col2) {
  return(colorRampPalette(c(col1, col2))(n))
}




