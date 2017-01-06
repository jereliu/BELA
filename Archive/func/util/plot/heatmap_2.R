heatmap_2 <- 
  function(obj, ...){
    colors <-
      c(seq(-1,-0.3,length=100),
        seq(-0.3,0.3,length=100),
        seq(0.3,1,length=100)) %>% unique %>% sort
    my_palette <- 
      colorRampPalette(
        c("blue", "white", "red"))(n = (length(colors)-1))
    
    heatmap.2(obj, 
              Rowv = FALSE, Colv = FALSE,
              dendrogram = "none", trace = "none",
              labRow = "", labCol = "", 
              col = my_palette, breaks = colors,
              ...)
  }