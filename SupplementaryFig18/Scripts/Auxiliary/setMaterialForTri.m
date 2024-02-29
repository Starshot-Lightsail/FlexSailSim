function setMaterialForTri( tri_idx, mat_idx)
% Assigns material properties such as thickness, modulus, etc. for the specified triangle tri_idx, from the specified
% material data structure mat_idx.  Because this is written as a function, it uses evalin('base',...) to execute all
% these commands, which makes it exceptionally slow.  Perhaps using globals would be faster, but that sounds ghetto to
% me.  Anyway, mapping material properties for large grids is really slow, but not so slow that I'm going to fix it.
if mat_idx > 0
    %global MATERIALS
                evalin('base',['t_thick(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').thickness;']);
                evalin('base',['t_dens(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').density;']);
                evalin('base',['t_ymod(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').Ymod;']);
                evalin('base',['t_emod(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').Emod;']);
                evalin('base',['t_nu(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').nu;']);
                evalin('base',['t_thermcd(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').Thermcondmm;']);
                evalin('base',['t_heatCap(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').heatCap;']);
                evalin('base',['t_CTE(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').CTE;']);
                evalin('base',['t_emis(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').Emissivity;']);
                evalin('base',['t_Iabs(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').Iabs;']);
                evalin('base',['t_Irefl(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').Irefl;']);
                evalin('base',['t_failtemp(' num2str(tri_idx) ') = MATERIALS(' num2str(mat_idx) ').failtemp;']);
                evalin('base',['t_mat(' num2str(tri_idx) ') = ' num2str(mat_idx) ';' ]);
end
end


% 
%  t_dens(tri_idx) = MATERIALS(mat_idx).density;
%                 t_ymod(tri_idx) = MATERIALS(mat_idx).Ymod;
%                 t_emod(tri_idx) = MATERIALS(mat_idx).Emod;
%                 t_nu(tri_idx) = MATERIALS(mat_idx).nu;
%                 t_thermcd(tri_idx) = MATERIALS(mat_idx).Thermcondmm;
%                 t_heatCap(tri_idx) = MATERIALS(mat_idx).heatCap;
%                 t_CTE(tri_idx) = MATERIALS(mat_idx).CTE;
%                 t_emis(tri_idx) = MATERIALS(mat_idx).Emissivity;
%                 t_Iabs(tri_idx) = MATERIALS(mat_idx).Iabs;
%                 t_Irefl(tri_idx) = MATERIALS(mat_idx).Irefl;