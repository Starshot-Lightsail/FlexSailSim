function edge_idx = findedge(node1_idx, node2_idx)

% Finds index of edge that starts and stops at the two given node indices.
% The order of the two nodes doesn't matter.  If there's an edge between
% them, this will return its index.  If not, it returns zero.

% Uses global e_na and e_nb arrays.  Iterates to length(e_na) but bails
% early if it finds an edge starting at node zero.

%This function should be improved if it's goign to be used in the main
%loop.  OK for meshing.  Will need to replace call to length(e_na) for
%C/C++ version...

global e_na e_nb 

edge_idx = 0;
for ne=1:length(e_na)
    if (e_na(ne) == node1_idx)
        if ( e_nb(ne) == node2_idx)
            edge_idx = ne;
            return;
        end
    elseif ( e_na(ne) == node2_idx )
        if ( e_nb(ne) == node1_idx )
            edge_idx = ne;
            return;
        end
    elseif ( e_na(ne) == 0 )
        return;
    end
end

end
