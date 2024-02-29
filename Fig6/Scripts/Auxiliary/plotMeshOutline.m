%% plotMeshOutline subscript
%
% Draws a line around the outer perimeter of the sail mesh, to help with visualization of of the sail shape in 3D plots.
%
% For non-broken sails (with all edges still intact) we use the "constrainedEdges" vector, which was created during 
% initial meshing to identify the edges comprising the outer perimeter.  This is reasonably fast, as we don't need to 
% calculate which edges form the perimeter, and we can use a single call to plot3() to render the perimeter.  
%
% However, once the edges start breaking, we instead look at all the edges to determine which ones ones appear only once
% in the list of remaining triangle edges, which indicates that they're on an outer edge of the remaining strucutre.  We
% then plot each of those outer edges with consecutive calls to plot3(), making this a slower rendering process, but the
% only one I've found that works with partially broken lightsail meshes.  


%input paramter:  plotMeshOutlineZ0.  This should be zero, except when the outline needs to line up with mesh faces that 
% were plotted with scewed z-axis scaling (i.e., when zreliefmag has a value other than 1) within "snapshot" figures, which
% use arbitrary z offsets to present renderings of the sail at various stages of simulation within a single compact figure.  

hold on;
max_mesh_tris=100000;
if any(t_broken)
    TRI3 = sort(TRI((~t_broken),:),2);
    EDG3 = [ max_mesh_tris*TRI3(:,1)+TRI3(:,2); max_mesh_tris*TRI3(:,2)+TRI3(:,3); max_mesh_tris*TRI3(:,1)+TRI3(:,3) ];
    EDG3 = sort(EDG3)';
    EDG4 = [EDG3 0 0];
    EDG5 = [0 EDG3 0];
    EDG6 = [0 0 EDG3];
    EDG_UNIQUE = ( EDG4 ~= EDG5) & ( EDG5 ~= EDG6 );
    EDG7 = EDG5(EDG_UNIQUE);
    EDGB = mod(EDG7,max_mesh_tris);
    EDGA = (EDG7-EDGB)/max_mesh_tris;

    for n=1:length(EDGA)
        plot3(n_x([EDGA(n) EDGB(n)]), n_y([EDGA(n) EDGB(n)]), zreliefmag*(n_z([EDGA(n) EDGB(n)])-plotMeshOutlineZ0), 'k','LineWidth',0.4);
    end
else
    plot3([n_x(constrainedEdges(:,1)); n_x(constrainedEdges(end,2))],[n_y(constrainedEdges(:,1)); n_y(constrainedEdges(end,2))],zreliefmag*([n_z(constrainedEdges(:,1)); n_z(constrainedEdges(end,2))]-plotMeshOutlineZ0),'k') 
end
plot3(n_x(1), n_y(1), zreliefmag*(n_z(1)-plotMeshOutlineZ0), 'k.');

