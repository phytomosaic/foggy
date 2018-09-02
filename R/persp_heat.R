### helper function to add heatmap colors to a perspective plot
`persp_heat` <- function(z, pal='heat', ...){
     nrz <- nrow(z) ; ncz <- ncol(z)
     zfacet <- z[-1,-1] + z[-1,-ncz] + z[-nrz,-1] + z[-nrz,-ncz]
     if(pal=='redblue'){
          cc <- colorRampPalette(c('red','transparent','blue'))(111)[
               as.numeric(cut(zfacet,breaks=111))]
     }
     if(pal=='heat'){
          cc <- heat.colors(100)[cut(zfacet, 100)]
     }
     persp(z, col=cc, ...)
}
