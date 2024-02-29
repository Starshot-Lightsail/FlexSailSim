function edge_elongation = edgeElong(edge_idx)

global e_na e_nb e_l0

edge_length = distanceNodes(e_na(edge_idx),e_nb(edge_idx));

edge_elongation = edge_length - e_l0(edge_idx);

end