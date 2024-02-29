%% Calculate Optical Forces vs Translation and Rotation
% For composite round lightsails, Mark 1e, Si3N4 at 1064 nm

% Compare calculated results with those in Mathematica

% Ramon Gao, February 19, 2022

%% Initialization
% Close all figures, clear command window, clear all stored variables

close all
clear
clc
clear textprogressbar

%% General parameters and settings

% Plot parameters
fontsize = 16;

% Global parameters
c0 = 299792458; % [c0] = m/s, speed of light
I0 = 1;     % [I0] = 1, normalized light intensity
I0_actual = 1e9; % in W/m^2
D = 1;      % [D] = 1, normalized diameter of lightsail
D_actual = 1; % in m
m = 1;      % [m] = 1, normalized mass

% Beam width, determines stability
w = 0.4*D;

% Orientation of gratings in their local frames
beta = [0, pi, pi/2, 3*pi/2];

% Position of boundary separating unpatterned and patterned regions
sx = 0*1/10; % [sx] = 1, unitless, actual position is [sx*D] = m
    
% Define radial boundaries of all regions/patterns (for angular boundaries,
% see below, since these are material-specific)
r_start_R1 = sx.*D;
r_start_R2 = r_start_R1;
r_start_R3 = r_start_R1;
r_start_R4 = r_start_R1;

r_stop_R1 = D./2;
r_stop_R2 = r_stop_R1;
r_stop_R3 = r_stop_R1;
r_stop_R4 = r_stop_R1;
    
% Define boundaries of all regions/patterns for a round lightsail
varphi_start_R1 = -pi/6;
varphi_stop_R1 = -varphi_start_R1;
varphi_start_R2 = pi - abs(varphi_stop_R1);
varphi_stop_R2 = pi + abs(varphi_stop_R1);
varphi_start_R3 = varphi_stop_R1;
varphi_stop_R3 = varphi_start_R2;
varphi_start_R4 = varphi_stop_R2;
varphi_stop_R4 = 2*pi - abs(varphi_stop_R1);

%% Load pressures

% Load collection of smaller look-up tables of pressures versus roll (phi)
% and phi (theta) angles for TE (region A) and TM polarization (region B)

global pressure_TE_collection_small_tables
global pressure_TM_collection_small_tables

% Load splitted smaller look-up tables of pressures vs roll and pitch
% angles for TE- and TM-polarized light. All the small look-up tables for
% TE (TM) polarization combined together yield the full look-up table for
% TE (TM) polarization. The full look-up table is being splitted into
% smaller ones to increase speed of finding pressures within the table
% later in the code

file_name_collection_small_tables = ...
    '08232021_MEPS_Si3N4_Mark1e_1064_PitchRoll_Pressures_CollectionTables_';

load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])

pressure_TE_collection_small_tables = pressure_TE_small_tables;
pressure_TM_collection_small_tables = pressure_TM_small_tables;

%% Assign important boundary phi angles

% Assign pressures in local x and z direction to vectors

dim_collection_small_tables = size(pressure_TE_small_tables);

collection_phi_angles_tmp = zeros(dim_collection_small_tables(3), ...
    dim_collection_small_tables(1));

global last_indices_tables
global num_small_tables

num_small_tables = dim_collection_small_tables(3);

last_indices_tables = zeros(1,num_small_tables);

for ii = 1:1:num_small_tables
    
    tmp = [];
    last_indices_tables(ii) = find(pressure_TE_small_tables(:,1,ii),1,'last');
    tmp = unique(pressure_TE_small_tables(1:last_indices_tables(ii),1,ii));
    collection_phi_angles_tmp(ii,1:length(tmp)) = tmp;
    
end

max_len_collection_phi_angles = 0;

for ii = 1:1:dim_collection_small_tables(3)
    if find(collection_phi_angles_tmp(ii,:),1,'last') > max_len_collection_phi_angles
        max_len_collection_phi_angles = find(collection_phi_angles_tmp(ii,:),1,'last');
    end
end

collection_phi_angles = zeros(num_small_tables, max_len_collection_phi_angles);

last_indices_angles = zeros(1,num_small_tables);

for ii = 1:1:num_small_tables
    
    last_indices_angles(ii) = find(collection_phi_angles_tmp(ii,:),1,'last');
    collection_phi_angles(ii,1:last_indices_angles(ii)) = ...
        collection_phi_angles_tmp(ii,1:last_indices_angles(ii));
  
end

global border_angles

first_angles = collection_phi_angles(:,1);
step_angles = abs(collection_phi_angles(:,1) - collection_phi_angles(:,2));

last_angles = zeros(size(first_angles));
for ii = 1:1:length(last_indices_angles)
    last_angles(ii) = collection_phi_angles(ii,last_indices_angles(ii));
end

border_angles = (last_angles(1:end-1) + first_angles(2:end))/2;

%% Definition of integrated auxiliary intensity functions

% Euler rotation sequence 1-2-3 (x-y'-z'')

integrandForce = @(r,varphi,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*((xc+r.*cos(psi).*cos(theta).*cos(varphi)+ ...
    (-1).*r.*cos(theta).*sin(psi).*sin(varphi)).^2+(r.*cos(varphi).*( ...
    sin(phi).*sin(psi)+(-1).*cos(phi).*cos(psi).*sin(theta))+r.*(cos( ...
    psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*sin(varphi)).^2+( ...
    yc+r.*cos(varphi).*(cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin( ...
    theta))+r.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin( ...
    theta)).*sin(varphi)).^2)).*r;
    
%% Definition of body-frame forces and torques
   
FxBodyFrame = @(xc,yc,psi,theta,phi) cos(theta).*cos(phi).*I0.* ...
    (c0./(I0.*D.^2)).* ...
    ( cos(beta(1)).*findpressure('px_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1) + ...
    cos(beta(2)).*findpressure('px_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2) + ...
    cos(beta(3)).*findpressure('px_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3) + ...
    cos(beta(4)).*findpressure('px_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4) );

FyBodyFrame = @(xc,yc,psi,theta,phi) cos(theta).*cos(phi).*I0.* ...
    (c0./(I0.*D.^2)).* ...
    ( sin(beta(1)).*findpressure('px_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1) + ...
    sin(beta(2)).*findpressure('px_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2) + ...
    sin(beta(3)).*findpressure('px_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3) + ...
    sin(beta(4)).*findpressure('px_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4) );

FzBodyFrame = @(xc,yc,psi,theta,phi) cos(theta).*cos(phi).*I0.* ...
    (c0./(I0.*D.^2)).* ...
    (findpressure('pz_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1) + ...
    findpressure('pz_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2) + ...
    findpressure('pz_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3) + ...
    findpressure('pz_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandForce(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4) );

%% Definition of direction cosine matrix (rotation matrix)

% Define function that describes the equations of the system in vector form

DRM123_11 = @(psi,theta,phi) cos(psi).*cos(theta);
DRM123_12 = @(psi,theta,phi) cos(psi)*sin(theta).*sin(phi) + sin(psi).*cos(phi);
DRM123_13 = @(psi,theta,phi) -cos(psi).*sin(theta).*cos(phi) + sin(psi).*sin(phi);
DRM123_21 = @(psi,theta,phi) -sin(psi).*cos(theta);
DRM123_22 = @(psi,theta,phi) -sin(psi).*sin(theta).*sin(phi) + cos(psi).*cos(phi);
DRM123_23 = @(psi,theta,phi) sin(psi)*sin(theta).*cos(phi) + cos(psi).*sin(phi);
DRM123_31 = @(psi,theta,phi) sin(theta);
DRM123_32 = @(psi,theta,phi) -cos(theta).*sin(phi);
DRM123_33 = @(psi,theta,phi) cos(theta).*cos(phi);

%% Definition of lab-frame forces

Fx = @(xc,yc,psi,theta,phi) ...
    FxBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_11(psi,theta,phi) + ...
    FyBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_21(psi,theta,phi) + ...
    FzBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_31(psi,theta,phi);

Fy = @(xc,yc,psi,theta,phi) ...
    FxBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_12(psi,theta,phi) + ...
    FyBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_22(psi,theta,phi) + ...
    FzBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_32(psi,theta,phi);

Fz = @(xc,yc,psi,theta,phi) ...
    FxBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_13(psi,theta,phi) + ...
    FyBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_23(psi,theta,phi) + ...
    FzBodyFrame(xc,yc,psi,theta,phi) .* ...
    DRM123_33(psi,theta,phi);

%% Calculate forces

start_theta = -10;
stop_theta = -start_theta;
step_theta = 0.1;

theta_list = start_theta:step_theta:stop_theta;
theta_list = deg2rad(theta_list);
num_theta = length(theta_list);

phi_list = theta_list;
num_phi = length(phi_list);

x_start = -0.8;
x_stop = -x_start;
x_step = 0.01;

x_list = x_start:x_step:x_stop;
num_x = length(x_list);

y_list = x_list;
num_y = length(y_list);

fx = zeros(num_theta*num_x,3);
fy = zeros(num_phi*num_y,3);

kk = 1;
for ii = 1:1:num_x
    for jj = 1:1:num_theta
        
        fx(kk,:) = [x_list(ii), rad2deg(theta_list(jj)), ...
            Fx(x_list(ii),0,0,theta_list(jj),0)];
        
        kk = kk + 1;
    end
    disp(ii)
end

kk = 1;
for ii = 1:1:num_y
    for jj = 1:1:num_phi
        
        fy(kk,:) = [x_list(ii), rad2deg(theta_list(jj)), ...
            Fy(0,y_list(ii),0,0,phi_list(jj))];
        
        kk = kk + 1;
    end
    disp(ii)
end

formatOut = 'mmddyy';

% Define file base name
save_name = ['data_Si3N4_',datestr(now,formatOut),'_M1e_1064nm_w',num2str(round(w,1)),'D'];

% Replace '-' with 'm' and '.' with 'p'
save_name = strrep(save_name,'-','m');
save_name = strrep(save_name,'.','p');

save(fullfile('Data',[save_name,'_generalFx.mat']),'fx');
save(fullfile('Data',[save_name,'_generalFy.mat']),'fy');


%% Auxiliary functions

function [pressure] = findpressure(p,angle0)
% Function to find a simulated tilt angle that is closed to the 
% numerically evolved angle by means of minimum Euclidean distance

global pressure_TE_collection_small_tables
global pressure_TM_collection_small_tables

global last_indices_tables
global border_angles
global num_small_tables

theta123 = angle0(2);
phi123 = angle0(1);

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