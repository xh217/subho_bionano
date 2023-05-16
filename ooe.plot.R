## Plotting the visual representation of the order of complex events
ooe.plot <- function(res) {
  library(bnlearn)
  nodes <- c()
  for (i in 1:length(res$nodes)) {
    nodes <- c(nodes, res$nodes[[i]]$mb)
    nodes <- unique(nodes)
  }
  
  inter_nodes <- c()
  for (i in 1:length(nodes)) {
    if (strsplit(nodes[i], "\\.")[[1]][1] != strsplit(nodes[i], "\\.")[[1]][3])
    {
      inter_nodes <- c(inter_nodes, nodes[i])
    }
  }
  if (length(inter_nodes) > 0) {
    highlight_list <- list(nodes = inter_nodes, fill = "orange")
    graphviz.plot(res, shape = "rectangle", highlight = highlight_list)
  } else NULL
}
