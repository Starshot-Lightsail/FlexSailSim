function area = triArea(tri_idx)
%returns area of the specified triangle

global n_x n_y n_z t_na t_nb t_nc

p1 = [n_x(t_na(tri_idx)) n_y(t_na(tri_idx)) n_z(t_na(tri_idx))];
p2 = [n_x(t_nb(tri_idx)) n_y(t_nb(tri_idx)) n_z(t_nb(tri_idx))];
p3 = [n_x(t_nc(tri_idx)) n_y(t_nc(tri_idx)) n_z(t_nc(tri_idx))];

cp = cross( p2-p1, p3-p1 );

area = 0.5 * sqrt(dot(cp,cp));

end
