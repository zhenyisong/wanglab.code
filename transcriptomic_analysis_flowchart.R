
# https://github.com/rich-iannone/DiagrammeR/issues/126
library(DiagrammeR)
library(magrittr)

graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of Mature gene selection") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    set_global_graph_attrs("node", "color", "gray75") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(6) %>%
    select_nodes_by_id(c(1:6)) %>% 
    set_node_attrs_ws("shape", "rectangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 3->4 3->5 5->6", "black") %>%
    set_node_attrs("label",c( 'Genome mapping','reads',
                              'Mapping to genome','Transcript identification &\ncounting',
                              'Transcript discovery &\ncounting','Functional annoation')) %>%
    set_node_attrs("fontsize",12) %>%
    set_edge_attrs("arrowsize", 1) %>%
    select_edges_by_node_id(nodes = 1) %>%
    set_edge_attrs_ws("arrowsize", 3)

render_graph(graph)



graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of Mature gene selection") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    set_global_graph_attrs("node", "color", "gray75") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(4) %>%
    select_nodes_by_id(c(1:4)) %>% 
    set_node_attrs_ws("shape", "rectangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 3->4", "black") %>%
    set_node_attrs("label",c( 'Transcriptome mapping','reads',
                              'Mapping to transcriptome','Transcript identification &\ncounting')) %>%
    set_node_attrs("fontsize",12) %>%
    set_edge_attrs("arrowsize", 1) %>%
    select_edges_by_node_id(nodes = 1) %>%
    set_edge_attrs_ws("arrowsize", 3)

render_graph(graph)



graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of Mature gene selection") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    set_global_graph_attrs("node", "color", "gray75") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(6) %>%
    select_nodes_by_id(c(1:6)) %>% 
    set_node_attrs_ws("shape", "rectangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 3->4 4->5 5->6", "black") %>%
    set_node_attrs("label",c( 'Reference free assembly','reads',
                              'Assembly into transcripts','Map reads back',
                              'Counting','Functional annotation')) %>%
    set_node_attrs("fontsize",12) %>%
    set_edge_attrs("arrowsize", 1) %>%
    select_edges_by_node_id(nodes = 1) %>%
    set_edge_attrs_ws("arrowsize", 3)

render_graph(graph)



graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of ncRNA analysis") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    #set_global_graph_attrs("graph", "fixedsize", "true") %>%
    set_global_graph_attrs("node", "color", "black") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(5) %>%
    select_nodes_by_id(c(1:5)) %>% 
    set_node_attrs_ws("shape", "ellipse") %>%
    
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 2->4 3->5 4->5", "black") %>%
    set_node_attrs("label",c( 'New tuxedo\nanalysis pipeline',
                              'Ab initio\ntranscriptome assembly',
                              'transcript comparison',"ad hoc filtering\nSebnif",
                              'qunatification\ndifferential expression')) %>%
    set_node_attrs("fontsize",10) %>%
    set_edge_attrs("arrowsize", 1)

render_graph(graph)


graph <-
    create_graph() %>%
    set_graph_name("Flow Chart of lncRNA analysis plus novel lncRNA prediction") %>%
    set_global_graph_attrs("graph", "overlap", "true") %>%
    #set_global_graph_attrs("graph", "fixedsize", "true") %>%
    set_global_graph_attrs("node", "color", "grey75") %>%
    set_global_graph_attrs("node", "fontname", "Helvetica") %>%
    add_n_nodes(4) %>%
    select_nodes_by_id(c(1:4)) %>% 
    set_node_attrs_ws("shape", "retangle") %>%
    set_node_attrs_ws("style", "filled") %>%
    clear_selection %>%
    add_edges_w_string(
      "1->2 2->3 3->4", "black") %>%
    set_node_attrs("label",c( 'Hisat2\nmap reads with reference genome assemble',
                              'Stringtie\ntranscript assemble with reference GTF',
                              'iSeeRNA\nnovel lncRNA prediction',
                              'Limma\nqunatification\ndifferential expression')) %>%
    set_node_attrs("fontsize",10) %>%
    set_edge_attrs("arrowsize", 1)

render_graph(graph)
