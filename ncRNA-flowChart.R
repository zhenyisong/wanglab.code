library(DiagrammeR)
library(tidyverse)

graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of ncRNA analysis") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    set_global_graph_attrs("node", "color", "grey75") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(5) %>%
    select_nodes_by_id(c(1:5)) %>% 
    set_node_attrs_ws("shape", "retangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 2->4 3->5 4->5", "black") %>%
    set_node_attrs("label",c( 'New tuxedo\nanalysis pipeline',
                              'Ab initio\ntranscriptome assembly',
                              'Transcript comparison',"ad hoc filtering\nSebnif",
                              'Qunatification\ndifferential expression')) %>%
    set_node_attrs("fontsize",10) %>%
    set_edge_attrs("arrowsize", 1)

render_graph(graph)

sessionInfo()