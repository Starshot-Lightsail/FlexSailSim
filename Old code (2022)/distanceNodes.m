function dist = distanceNodes(node1_idx, node2_idx)
global n_x n_y n_z

   dx = n_x(node2_idx) - n_x(node1_idx);
   dy = n_y(node2_idx) - n_y(node1_idx);
   dz = n_z(node2_idx) - n_z(node1_idx);
   
   dist = sqrt( dx.^2 + dy.^2 + dz.^2 );
end
