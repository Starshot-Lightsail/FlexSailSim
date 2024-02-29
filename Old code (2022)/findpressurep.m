function [pressure] = findpressurep(pressure_TE_collection_small_tables, pressure_TM_collection_small_tables, p,angle0)
% Function to find a simulated tilt angle that is closed to the
% numerically evolved angle by means of minimum Euclidean distance

%pressure_TE_collection_small_tables = evalin('base','pressure_TE_collection_small_tables');
%size(pressure_TE_collection_small_tables)

%pressure_TM_collection_small_tables = evalin('base','pressure_TM_collection_small_tables');

last_indices_tables = [2581          1785          3585          2390          3585          6180          6180          6180          6180          6180          6180          6180          3585          2390          3585          1785          2670 ];
border_angles = [-0.541052068118  -0.270526034059  -0.211184839491  -0.176278254451  -0.123918376892  -0.088139127226  -0.053232542186  -0.018325957146  0.016580627894  0.051487212934  0.086393797974  0.121300383014  0.172787595947  0.207694180987  0.260054058547  0.514872129338 ];
num_small_tables = 17;

euler_seq = 123;

% if euler_seq == 312
%     
%     R31 = cos(angle0(3)).*sin(angle0(2)) + cos(angle0(2)).*sin(angle0(1)).*sin(angle0(3));
%     R32 = -cos(angle0(2)).*cos(angle0(3)).*sin(angle0(1)) + sin(angle0(2)).*sin(angle0(3));
%     R33 = cos(angle0(2)).*cos(angle0(1));
%     % R21 = -cos(angle0(1)).*sin(angle0(3));
%     % R11 = cos(angle0(2)).*cos(angle0(3)) - sin(angle0(2)).*sin(angle0(1)).*sin(angle0(3));
%     
%     theta123 = asin(R31);
%     phi123 = -atan2(R32./cos(theta123),R33./cos(theta123));
%     % psi123 = -atan2(R21./cos(theta123),R11./cos(theta123));
%     
% elseif euler_seq == 123
    
    theta123 = angle0(2);
    phi123 = angle0(1);
    
% end

if strcmpi(p,'px_R1')
    if phi123 < border_angles(1)
        idx_R1 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(1),1:2,1), [phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R1,3,1);
    elseif phi123 >= border_angles(end)
        idx_R1 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables),[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R1,3,num_small_tables);
    elseif (border_angles(1) <= phi123) && (phi123 < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= phi123) & ...
            (phi123 < border_angles(2:end)) == 1);
        idx_R1 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R1,3,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'px_R2')
    if (-phi123) < border_angles(1)
        idx_R2 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(1),1:2,1), -[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R2,3,1);
    elseif (-phi123) >= border_angles(end)
        idx_R2 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables), -[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R2,3,num_small_tables);
    elseif (border_angles(1) <= (-phi123)) && ((-phi123) < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= (-phi123)) & ...
            ((-phi123) < border_angles(2:end)) == 1);
        idx_R2 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),-[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R2,3,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'px_R3')
    if phi123 < border_angles(1)
        idx_R3 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(1),1:2,1), [phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R3,3,1);
    elseif phi123 >= border_angles(end)
        idx_R3 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables),[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R3,3,num_small_tables);
    elseif (border_angles(1) <= phi123) && (phi123 < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= phi123) & ...
            (phi123 < border_angles(2:end)) == 1);
        idx_R3 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R3,3,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'px_R4')
    if (-phi123) < border_angles(1)
        idx_R4 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(1),1:2,1), -[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R4,3,1);
    elseif (-phi123) >= border_angles(end)
        idx_R4 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables), -[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R4,3,num_small_tables);
    elseif (border_angles(1) <= (-phi123)) && ((-phi123) < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= (-phi123)) & ...
            ((-phi123) < border_angles(2:end)) == 1);
        idx_R4 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),-[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R4,3,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'pz_R1')
    if phi123 < border_angles(1)
        idx_R1 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(1),1:2,1), [phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R1,4,1);
    elseif phi123 >= border_angles(end)
        idx_R1 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables),[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R1,4,num_small_tables);
    elseif (border_angles(1) <= phi123) && (phi123 < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= phi123) & ...
            (phi123 < border_angles(2:end)) == 1);
        idx_R1 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R1,4,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'pz_R2')
    if (-phi123) < border_angles(1)
        idx_R2 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(1),1:2,1), -[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R2,4,1);
    elseif (-phi123) >= border_angles(end)
        idx_R2 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables), -[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R2,4,num_small_tables);
    elseif (border_angles(1) <= (-phi123)) && ((-phi123) < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= (-phi123)) & ...
            ((-phi123) < border_angles(2:end)) == 1);
        idx_R2 = knnsearch(pressure_TE_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),-[phi123, theta123]);
        pressure = pressure_TE_collection_small_tables(idx_R2,4,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'pz_R3')
    if phi123 < border_angles(1)
        idx_R3 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(1),1:2,1), [phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R3,4,1);
    elseif phi123 >= border_angles(end)
        idx_R3 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables),[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R3,4,num_small_tables);
    elseif (border_angles(1) <= phi123) && (phi123 < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= phi123) & ...
            (phi123 < border_angles(2:end)) == 1);
        idx_R3 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R3,4,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
elseif strcmpi(p,'pz_R4')
    if (-phi123) < border_angles(1)
        idx_R4 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(1),1:2,1), -[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R4,4,1);
    elseif (-phi123) >= border_angles(end)
        idx_R4 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(end),1:2,num_small_tables), -[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R4,4,num_small_tables);
    elseif (border_angles(1) <= (-phi123)) && ((-phi123) < border_angles(end))
        idx_small_table = find((border_angles(1:end-1) <= (-phi123)) & ...
            ((-phi123) < border_angles(2:end)) == 1);
        idx_R4 = knnsearch(pressure_TM_collection_small_tables(1:last_indices_tables(idx_small_table+1),1:2,idx_small_table+1),-[phi123, theta123]);
        pressure = pressure_TM_collection_small_tables(idx_R4,4,idx_small_table+1);
    else
        warning('Error: angle out of range!');
    end
end

end
