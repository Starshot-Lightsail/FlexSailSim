function edge_length = edgeLength(edge_idx)

global e_na e_nb 

edge_length = distanceNodes(e_na(edge_idx),e_nb(edge_idx));

end