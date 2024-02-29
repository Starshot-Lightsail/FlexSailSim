% This simple script is to assign global variables to temporary variables
% for the core simulation script of force calculation

% May 17th, 2023
% Ramon Gao

 % Assign global variales to temporary variables
n_x_tmp = n_x;
n_y_tmp = n_y;
n_z_tmp = n_z;

% Update triangle centroids
t_cx_tmp = ( n_x_tmp(t_na) + n_x_tmp(t_nb) + n_x_tmp(t_nc) ) / 3;
t_cy_tmp = ( n_y_tmp(t_na) + n_y_tmp(t_nb) + n_y_tmp(t_nc) ) / 3;
t_cz_tmp = ( n_z_tmp(t_na) + n_z_tmp(t_nb) + n_z_tmp(t_nc) ) / 3;

% Update time step
tt_tmp = tt;

% Update edge and triangle temperatures
e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;
t_t = ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) ./ 3;

% Assign node, edge and triangle temperatures
n_t_tmp = n_t;
e_t_tmp = e_t;
t_t_tmp = t_t;

