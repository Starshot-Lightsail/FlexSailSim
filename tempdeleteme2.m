% M_nxyz2tev1 = zeros(numnodes,numtris);
% M_nxyz2tev2 = zeros(numnodes,numtris);
% 
% for nt = 1:length(t_na)
%     M_nxyz2tev1( t_na(nt), nt) = -1;
%     M_nxyz2tev1( t_nb(nt), nt) = 1;
%     M_nxyz2tev2( t_na(nt), nt) = -1;
%     M_nxyz2tev2( t_nc(nt), nt) = 1;
% end


numtries = 1000;

A = [ n_x(t_nb) - n_x(t_na) ; n_y(t_nb) - n_y(t_na) ; n_z(t_nb) - n_z(t_na) ];
            B = [ n_x(t_nc) - n_x(t_na) ; n_y(t_nc) - n_y(t_na) ; n_z(t_nc) - n_z(t_na) ];
            
tic;

for n=1:numtries
    A = [ n_x(t_nb) - n_x(t_na) ; n_y(t_nb) - n_y(t_na) ; n_z(t_nb) - n_z(t_na) ];
    B = [ n_x(t_nc) - n_x(t_na) ; n_y(t_nc) - n_y(t_na) ; n_z(t_nc) - n_z(t_na) ];
end
toc


C = A;
D = B;

n_xyz = [n_x; n_y; n_z];

tic;
for n=1:numtries
%     C(1,:) = n_x * M_nxyz2tev1;
%     C(2,:) = n_y * M_nxyz2tev1;
%     C(3,:) = n_z * M_nxyz2tev1;
%     D(1,:) = n_x * M_nxyz2tev2;
%     D(2,:) = n_y * M_nxyz2tev2;
%     D(3,:) = n_z * M_nxyz2tev2;
C = [n_x; n_y; n_z]* M_nxyz2tev1;
    D = [n_x; n_y; n_z] * M_nxyz2tev2;
end
toc

sum(sum(abs(C-A))) + sum(sum(abs(D-B)))


