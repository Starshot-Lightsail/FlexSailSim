function edge_strain = edgeStrain(edge_idx)

global e_na e_nb e_l0

edge_length = distanceNodes(e_na(edge_idx),e_nb(edge_idx));

edge_strain = edge_length / e_l0(edge_idx) - 1;

end