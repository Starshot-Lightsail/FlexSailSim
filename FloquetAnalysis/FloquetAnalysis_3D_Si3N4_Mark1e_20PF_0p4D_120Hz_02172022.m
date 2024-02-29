%% Stability analysis for spinning lightsails

% Do Floquet analysis to analyze stability of spinning lightsails by 
% (1) calculating the Jacobian, 
% (2) the state transition matrix by solving the differential equation based on the Jacobian, 
% (3) eigenvalues of the state transition matrix at t = T, where T is the period

% Approximate partial derivatives numerically for composite 3D MEPS made of
% high-stress stoichiometric silicon nitride (Si3N4)

% Unstable case: Mark 1e, but resonators 20% farther away, w = 0.3D, f = 120 Hz

% Ramon Gao, 02/17/2022

close all
clear
clc

%% Global parameters and settings

% Set material of lightsail to be 'SOI' or 'SiNx'
material = 'Si3N4';

% Sail shape: string, 'round' or 'square'
shape = 'round';

% Global parameters
c0 = 299792458; % [c0] = m/s, speed of light
I0 = 1;     % [I0] = 1, normalized light intensity
I0_actual = 1e9; % in W/m^2
D = 1;      % [D] = 1, normalized diameter of lightsail
D_actual = 1; % in m
m = 1;      % [m] = 1, normalized mass

% Beam width, determines stability
if strcmpi(material,'SOI')
    w = 2*D;    % [w] = 1, normalized beam width
elseif strcmpi(material,'SiNx2')
    w = 0.1*D;
elseif strcmpi(material,'Si3N4')
    w = 0.4*D;
end

% Set spinning frequency
spinning_freq = 120; % in Hz

% Orientation of gratings in their local frames
beta = [0, pi, pi/2, 3*pi/2];

% Position of boundary separating unpatterned and patterned regions
if strcmpi(shape,'round')
    sx = 0*1/10; % [sx] = 1, unitless, actual position is [sx*D] = m
elseif strcmpi(shape,'square')
    sx = 1/4; % [sx] = 1, unitless, actual position is [sx*D] = m
end

% Densities
rho_Si3N4 = 2700*1e-6; % [rho] = kg/m^2/um, density of high-stress Si3N4
rho_SiNx = 3100*1e-6; % [rho] = kg/m^2/um, density of low-stress SiNx
rho_SiO2 = 2200*1e-6; % [rho] = kg/m^2/um, density of silica
rho_Si = 2300*1e-6; % [rho] = kg/m^2/um, density of silicon

% Distribution of regions 1 and 2 on lightsail area
n = 4;  % number of wedges/pie slices
na = 2; % number of wedges/pie slices being region A
nb = n - na; % number of wedges/pie slices being region B
    
% Define radial boundaries of all regions/patterns (for angular boundaries,
% see below, since these are material-specific)
if strcmpi(shape,'round')
    
    r_start_R1 = sx.*D;
    r_start_R2 = r_start_R1;
    r_start_R3 = r_start_R1;
    r_start_R4 = r_start_R1;
    
    r_stop_R1 = D./2;
    r_stop_R2 = r_stop_R1;
    r_stop_R3 = r_stop_R1;
    r_stop_R4 = r_stop_R1;
    
elseif strcmpi(shape,'square')
    
    x_start_R1 = sx.*D;
    x_start_R2 = -D./2;
    x_start_R3 = -sx.*D;
    x_start_R4 = -sx.*D;
    
    y_start_R1 = -D./2;
    y_start_R2 = -D./2;
    y_start_R3 = 0;
    y_start_R4 = -D./2;

    x_stop_R1 = D./2;
    x_stop_R2 = -sx.*D;
    x_stop_R3 = sx.*D;
    x_stop_R4 = sx.*D;
    
    y_stop_R1 = D./2;
    y_stop_R2 = D./2;
    y_stop_R3 = D./2;
    y_stop_R4 = 0;
    
end

%% Define and calculate material-specific parameters

if strcmpi(material,'SOI')
    
    % Define densities of substrate and resonating elements
    rho_s = rho_SiO2; % [rho] = kg/m^2/um, density of silica
    rho_e = rho_Si; % [rho] = kg/m^2/um, density of silicon
    
    % Define design parameters of unit cells A and B 
    dA = 1.8; % [dA] = um, size of unit cell A
    dB = 1.775; % [dB] = um, size of unit cell B
    w1A = 0.15*dA; % [w1A] = um, width of first resonator in region A
    w2A = 0.35*dA; % [w2A] = um, width of second resonator in region A
    w1B = 0.125*dB; % [w1B] = um, width of first resonator in region B
    w2B = 0.25*dB; % [w2B] = um, width of second resonator in region B    

    tA = 0.5; % [ts] = um, thickness of substrate in region A
    hA = 0.5; % [hr] = um, height of resonators in region A
    tB = tA; % [ts] = um, thickness of substrate in region B
    hB = 0.45; % [hr] = um, height of resonators in region B

    % Define the mass per unit area of each region (unpatterned, A or B)
    chi0 = rho_s*tA;
    chiA = rho_s*tA + rho_e*hB*(w1A + w2A)/dA;
    chiB = rho_s*tA + rho_e*hB*(w1B + w2B)/dB;
    
    % Define boundaries of all regions for a round lightsail
    if strcmpi(shape,'round')
        varphi_start_R1 = -pi/4;
        varphi_stop_R1 = -varphi_start_R1;
        varphi_start_R2 = pi - abs(varphi_stop_R1);
        varphi_stop_R2 = pi + abs(varphi_stop_R1);
        varphi_start_R3 = varphi_stop_R1;
        varphi_stop_R3 = varphi_start_R2;
        varphi_start_R4 = varphi_stop_R2;
        varphi_stop_R4 = 2*pi - abs(varphi_stop_R1);
        
        % Define the total mass of the round lightsail, normalized by D^2
        mass = (pi/4)*( (4*abs(varphi_start_R1)/(2*pi))*(chiA - chiB) + chiB );  
        
    else
        % Define the total mass of the round lightsail, normalized by D^2
        mass = (1/n)*pi*(1/4 - sx^2)*(chiA*na + chiB*nb) + chi0*pi*sx^2;

    end
elseif strcmpi(material,'SiNx2')
        
    % Define density of substrate and resonating elements
    rho_s = rho_SiNx; % [rho] = kg/m^2/um
    rho_e = rho_SiNx; % [rho] = kg/m^2/um

    % Define design parameters of unit cells A and B 
    w1A = 0.10; % [w1A] = um, width of first resonator in region A
    w2A = 0.24; % [w2A] = um, width of second resonator in region A
    w1B = 0.12; % [w1B] = um, width of first resonator in region B
    w2B = 0.30; % [w2B] = um, width of second resonator in region B

    dUC_A = 0.73; % [dUC_A] = um, period of unit cell of region A
    dUC_B = 0.72; % [dUC_B] = um, period of unit cell of region B

    tA = 0.1; % [ts] = um, thickness of substrate
    hA = 0.2; % [hr] = um, height of resonators
    tB = tA; % [ts] = um, thickness of substrate
    hB = hA; % [hr] = um, height of resonators

    % Define the mass per unit area of each region (unpatterned, A or B)
    chi0 = rho_s*tA;
    chiA = rho_s*tA + rho_e*hA*(w1A + w2A)/dUC_A;
    chiB = rho_s*tB + rho_e*hB*(w1B + w2B)/dUC_B;

    if strcmpi(shape,'round')

        % Define boundaries of all regions/patterns for a round lightsail
        varphi_start_R1 = -pi/6;
        varphi_stop_R1 = -varphi_start_R1;
        varphi_start_R2 = pi - abs(varphi_stop_R1);
        varphi_stop_R2 = pi + abs(varphi_stop_R1);
        varphi_start_R3 = varphi_stop_R1;
        varphi_stop_R3 = varphi_start_R2;
        varphi_start_R4 = varphi_stop_R2;
        varphi_stop_R4 = 2*pi - abs(varphi_stop_R1);
        
        % Define the total mass of the round lightsail, normalized by D^2
        mass = (pi/4)*( (4*abs(varphi_start_R1)/(2*pi))*(chiA - chiB) + chiB );
        
    elseif strcmpi(shape,'square')
        
        % Define the total mass of the round lightsail, normalized by D^2
        mass = 2*(chiA*(1/2-sx) + chiB*sx);
        
    end
    
    
elseif strcmpi(material,'Si3N4')
    
    % Define density of substrate and resonating elements
    rho_s = rho_Si3N4; % [rho] = kg/m^2/um
    rho_e = rho_Si3N4; % [rho] = kg/m^2/um

    % Define design parameters of unit cells A and B 
    w1A = 0.2; % [w1A] = um, width of first resonator in region A
    w2A = 0.6; % [w2A] = um, width of second resonator in region A
    w1B = 0.2; % [w1B] = um, width of first resonator in region B
    w2B = 0.52; % [w2B] = um, width of second resonator in region B
    
    dUC_A = 1.60; % [dUC_A] = um, period of unit cell of region A
    dUC_B = 1.35; % [dUC_B] = um, period of unit cell of region B

    tA = 0.2; % [ts] = um, thickness of substrate
    hA = 0.4; % [hr] = um, height of resonators
    tB = tA; % [ts] = um, thickness of substrate
    hB = hA; % [hr] = um, height of resonators

    % Define the mass per unit area of each region (unpatterned, A or B)
    chi0 = rho_s*tA;
    chiA = rho_s*tA + rho_e*hA*(w1A + w2A)/dUC_A;
    chiB = rho_s*tB + rho_e*hB*(w1B + w2B)/dUC_B;
    
    if strcmpi(shape,'round')

        % Define boundaries of all regions/patterns for a round lightsail
        varphi_start_R1 = -pi/6;
        varphi_stop_R1 = -varphi_start_R1;
        varphi_start_R2 = pi - abs(varphi_stop_R1);
        varphi_stop_R2 = pi + abs(varphi_stop_R1);
        varphi_start_R3 = varphi_stop_R1;
        varphi_stop_R3 = varphi_start_R2;
        varphi_start_R4 = varphi_stop_R2;
        varphi_stop_R4 = 2*pi - abs(varphi_stop_R1);
        
        % Define the total mass of the round lightsail, normalized by D^2
        mass = (pi/4)*( (4*abs(varphi_start_R1)/(2*pi))*(chiA - chiB) + chiB );
        
    elseif strcmpi(shape,'square')
        
        % Define the total mass of the round lightsail, normalized by D^2
        mass = 2*(chiA*(1/2-sx) + chiB*sx);
        
    end

end

if strcmpi(shape,'round')
    
    % Define integrands for mass integration to get moments of inertia of a
    % round lightsail
    integrand_Ix = @(varphi,r) r.^2 .* ((sin(varphi)).^2) .* r;
    integrand_Iy = @(varphi,r) r.^2 .* ((cos(varphi)).^2) .* r;
    
    Ix_eval = chi0 * integral2(integrand_Ix,0,2*pi,0,sx.*D) + ...
        chiA .* ( ...
        integral2(integrand_Ix,varphi_start_R1,varphi_stop_R1,r_start_R1,r_stop_R1) + ...
        integral2(integrand_Ix,varphi_start_R2,varphi_stop_R2,r_start_R2,r_stop_R2)) + ...
        chiB .* ( ...
        integral2(integrand_Ix,varphi_start_R3,varphi_stop_R3,r_start_R3,r_stop_R3) + ...
        integral2(integrand_Ix,varphi_start_R4,varphi_stop_R4,r_start_R4,r_stop_R4));
    
    Iy_eval = chi0 * integral2(integrand_Iy,0,2*pi,0,sx.*D) + ...
        chiA .* ( ...
        integral2(integrand_Iy,varphi_start_R1,varphi_stop_R1,r_start_R1,r_stop_R1) + ...
        integral2(integrand_Iy,varphi_start_R2,varphi_stop_R2,r_start_R2,r_stop_R2)) + ...
        chiB .* ( ...
        integral2(integrand_Iy,varphi_start_R3,varphi_stop_R3,r_start_R3,r_stop_R3) + ...
        integral2(integrand_Iy,varphi_start_R4,varphi_stop_R4,r_start_R4,r_stop_R4));

elseif strcmpi(shape,'square')
    
    Ix_eval = (chiA.*(1/2 - sx) + chiB.*sx) .* D.^4/6;
    Iy_eval = (chiA./12 + 2.*(chiB - chiA).*sx.^3/3) .* D.^4;    
    
end

% Define prefactors for moments of inertias to be multiplied with m*D^2
prefactor_Ix = Ix_eval / mass;
prefactor_Iy = Iy_eval / mass;
prefactor_Iz = prefactor_Ix + prefactor_Iy;

% Define moments of inertia in normalzed units, i.e. in terms of m and D
Ix = prefactor_Ix .* m .* D.^2;
Iy = prefactor_Iy .* m .* D.^2;
Iz = prefactor_Iz .* m .* D.^2;

%% Load pressures

% Load splitted smaller look-up tables of pressures vs roll and pitch
% angles for TE- and TM-polarized light. All the small look-up tables for
% TE (TM) polarization combined together yield the full look-up table for
% TE (TM) polarization. The full look-up table is being splitted into
% smaller ones to increase speed of finding pressures within the table
% later in the code

if strcmpi(material,'SOI')

    file_name_collection_small_tables = ...
        '04032021_MEPS_SOI_Oggy_ReducedPitchRoll_Pressures_CollectionTables_';
    
    load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
    load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])                    

elseif strcmpi(material,'SiNx1')
    
    file_name_collection_small_tables = ...
        '03222021_MEPS_SiNx_Mark6e1_ReducedPitchRoll_Pressures_CollectionTables_';
    
    load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
    load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])
    
elseif strcmpi(material,'SiNx2')
    
    file_name_collection_small_tables = ...
        '04052021_MEPS_SiNx_Mark6e2_ReducedPitchRoll_Pressures_CollectionTables_';
    
    load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
    load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])
    
elseif strcmpi(material,'Si3N4')
    
    file_name_collection_small_tables = ...
        '11052021_MEPS_Si3N4_M1e_20PF_1064_PitchRoll_Pressures_CollectionTables_';

    load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
    load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])
    
end

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

%% Select subset of table entries from LUT

theta_start = deg2rad(-0.1);
phi_start = theta_start;

theta_stop = -theta_start;
phi_stop = -phi_start;

size_LUT = size(pressure_TE_small_tables);

for ii = 1:1:size_LUT(3)
    if min(pressure_TE_small_tables(:,1,ii)) > theta_start
        idx_page_LUT = ii - 1;
        break
    end
end

pressure_TE_table = pressure_TE_small_tables(:,:,idx_page_LUT);
pressure_TM_table = pressure_TM_small_tables(:,:,idx_page_LUT);

pressure_table_angles = pressure_TE_table(:,1:2);

indices_small_phi_angles = (pressure_table_angles(:,1) >= phi_start) & ...
    (pressure_table_angles(:,1) <= phi_stop);

indices_small_theta_angles = (pressure_table_angles(:,2) >= theta_start) & ...
    (pressure_table_angles(:,2) <= theta_stop);

indices_small_angles = find(indices_small_phi_angles + ...
    indices_small_theta_angles == 2);

pressure_TE_small_table = pressure_TE_table(indices_small_angles,:);
pressure_TM_small_table = pressure_TM_table(indices_small_angles,:);

%% Prepare data as grids for plotting & further processing

small_angle_list = unique(pressure_TE_small_table(:,1));

[phi_small_angles,theta_small_angles] = meshgrid(small_angle_list,small_angle_list);

px_TE_R1 = griddata(pressure_TE_small_table(:,1), pressure_TE_small_table(:,2), ...
    pressure_TE_small_table(:,3),phi_small_angles,theta_small_angles);

pz_TE_R1 = griddata(pressure_TE_small_table(:,1), pressure_TE_small_table(:,2), ...
    pressure_TE_small_table(:,4),phi_small_angles,theta_small_angles);

px_TM_R3 = griddata(pressure_TM_small_table(:,1), pressure_TM_small_table(:,2), ...
    pressure_TM_small_table(:,3),phi_small_angles,theta_small_angles);

pz_TM_R3 = griddata(pressure_TM_small_table(:,1), pressure_TM_small_table(:,2), ...
    pressure_TM_small_table(:,4),phi_small_angles,theta_small_angles);

%% Plot pressures at small angles
if 0
    
figure(1)
surf(phi_small_angles,theta_small_angles,px_TE_R1)
title('p_{x} for Region 1 (TE)')

figure(2)
surf(phi_small_angles,theta_small_angles,pz_TE_R1)
title('p_{z} for Region 1 (TE)')

figure(3)
surf(phi_small_angles,theta_small_angles,px_TM_R3)
title('p_{x} for Region 3 (TM)')

figure(4)
surf(phi_small_angles,theta_small_angles,pz_TM_R3)
title('p_{z} for Region 3 (TM)')

figure(5)
imagesc(([phi_start phi_stop]), ([theta_start theta_stop]), px_TE_R1);
xlabel('\phi (\circ)');
ylabel('\theta (\circ)');
colorbar;

hold on;
plot(rad2deg(pressure_TE_small_table(:,1)), ...
    rad2deg(pressure_TE_small_table(:,2)), 'k.', 'MarkerSize', 1)

figure(6)
imagesc(([phi_start phi_stop]), ([theta_start theta_stop]), pz_TE_R1);
xlabel('\phi (\circ)');
ylabel('\theta (\circ)');
colorbar;

hold on;
plot(rad2deg(pressure_TE_small_table(:,1)), ...
    rad2deg(pressure_TE_small_table(:,2)), 'k.', 'MarkerSize', 1)

figure(7)
imagesc(([phi_start phi_stop]), ([theta_start theta_stop]), px_TM_R3);
xlabel('\phi (\circ)');
ylabel('\theta (\circ)');
colorbar;

hold on;
plot(rad2deg(pressure_TE_small_table(:,1)), ...
    rad2deg(pressure_TE_small_table(:,2)), 'k.', 'MarkerSize', 1)

figure(8)
imagesc(([phi_start phi_stop]), ([theta_start theta_stop]), pz_TM_R3);
xlabel('\phi (\circ)');
ylabel('\theta (\circ)');
colorbar;

hold on;
plot(rad2deg(pressure_TE_small_table(:,1)), ...
    rad2deg(pressure_TE_small_table(:,2)), 'k.', 'MarkerSize', 1)

end

%% Define functional form of pressures to be fitted

interp_fun_px_TE_R1 = @(xi,yi) interp2(phi_small_angles,theta_small_angles,px_TE_R1,xi,yi);
interp_fun_pz_TE_R1 = @(xi,yi) interp2(phi_small_angles,theta_small_angles,pz_TE_R1,xi,yi);
interp_fun_px_TM_R3 = @(xi,yi) interp2(phi_small_angles,theta_small_angles,px_TM_R3,xi,yi);
interp_fun_pz_TM_R3 = @(xi,yi) interp2(phi_small_angles,theta_small_angles,pz_TM_R3,xi,yi);

%% Definition of integrated auxiliary intensity functions

if strcmpi(shape,'round')

integrandForce = @(r,varphi,w,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*((xc+r.*cos(psi).*cos(theta).*cos(varphi)+ ...
    (-1).*r.*cos(theta).*sin(psi).*sin(varphi)).^2+(r.*cos(varphi).*( ...
    sin(phi).*sin(psi)+(-1).*cos(phi).*cos(psi).*sin(theta))+r.*(cos( ...
    psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*sin(varphi)).^2+( ...
    yc+r.*cos(varphi).*(cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin( ...
    theta))+r.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin( ...
    theta)).*sin(varphi)).^2)).*r;
    
integrandTorqueX = @(r,varphi,w,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*((xc+r.*cos(psi).*cos(theta).*cos(varphi)+ ...
    (-1).*r.*cos(theta).*sin(psi).*sin(varphi)).^2+(r.*cos(varphi).*( ...
    sin(phi).*sin(psi)+(-1).*cos(phi).*cos(psi).*sin(theta))+r.*(cos( ...
    psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*sin(varphi)).^2+( ...
    yc+r.*cos(varphi).*(cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin( ...
    theta))+r.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin( ...
    theta)).*sin(varphi)).^2)).*r.^2.*sin(varphi);

integrandTorqueY = @(r,varphi,w,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*((xc+r.*cos(psi).*cos(theta).*cos(varphi)+ ...
    (-1).*r.*cos(theta).*sin(psi).*sin(varphi)).^2+(r.*cos(varphi).*( ...
    sin(phi).*sin(psi)+(-1).*cos(phi).*cos(psi).*sin(theta))+r.*(cos( ...
    psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*sin(varphi)).^2+( ...
    yc+r.*cos(varphi).*(cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin( ...
    theta))+r.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin( ...
    theta)).*sin(varphi)).^2)).*r.^2.*cos(varphi);

aux_force_R1 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandForce(r,varphi,w,xc,yc,psi,theta,phi),r_start_R1,r_stop_R1,...
    varphi_start_R1,varphi_stop_R1);

aux_force_R2 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandForce(r,varphi,w,xc,yc,psi,theta,phi),r_start_R2,r_stop_R2,...
    varphi_start_R2,varphi_stop_R2);

aux_force_R3 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandForce(r,varphi,w,xc,yc,psi,theta,phi),r_start_R3,r_stop_R3,...
    varphi_start_R3,varphi_stop_R3);

aux_force_R4 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandForce(r,varphi,w,xc,yc,psi,theta,phi),r_start_R4,r_stop_R4,...
    varphi_start_R4,varphi_stop_R4);

aux_torqueX_R1 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueX(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1);

aux_torqueX_R2 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueX(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2);

aux_torqueX_R3 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueX(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3);

aux_torqueX_R4 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueX(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4);

aux_torqueY_R1 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueY(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1);

aux_torqueY_R2 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueY(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2);

aux_torqueY_R3 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueY(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3);

aux_torqueY_R4 = @(w,xc,yc,psi,theta,phi) integral2(@(r,varphi) ...
    integrandTorqueY(r,varphi,w,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4);

elseif strcmpi(shape,'square')

integrandForce = @(w,xb,yb,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*(((-1).*yb.*cos(theta).*sin(phi)+xb.*sin( ...
    theta)).^2+(xc+xb.*cos(psi).*cos(theta)+yb.*(cos(phi).*sin(psi)+ ...
    cos(psi).*sin(phi).*sin(theta))).^2+(yc+(-1).*xb.*cos(theta).*sin( ...
    psi)+yb.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin(theta)) ...
    ).^2));

integrandTorqueX = @(w,xb,yb,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*(((-1).*yb.*cos(theta).*sin(phi)+xb.*sin( ...
    theta)).^2+(xc+xb.*cos(psi).*cos(theta)+yb.*(cos(phi).*sin(psi)+ ...
    cos(psi).*sin(phi).*sin(theta))).^2+(yc+(-1).*xb.*cos(theta).*sin( ...
    psi)+yb.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin(theta)) ...
    ).^2)).*yb;

integrandTorqueY = @(w,xb,yb,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*(((-1).*yb.*cos(theta).*sin(phi)+xb.*sin( ...
    theta)).^2+(xc+xb.*cos(psi).*cos(theta)+yb.*(cos(phi).*sin(psi)+ ...
    cos(psi).*sin(phi).*sin(theta))).^2+(yc+(-1).*xb.*cos(theta).*sin( ...
    psi)+yb.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin(theta)) ...
    ).^2)).*xb;

end

%% Definition of body-frame forces and torques

fx_BF = @(w,x,y,psi,theta,phi) cos(theta).*cos(phi).*I0.*(c0./(I0.*D.^2)).*( ...
    cos(beta(1)).*interp_fun_px_TE_R1(phi,theta).*aux_force_R1(w,x,y,psi,theta,phi) + ...
    cos(beta(2)).*interp_fun_px_TE_R1(-phi,-theta).*aux_force_R2(w,x,y,psi,theta,phi) + ...
    cos(beta(3)).*interp_fun_px_TM_R3(phi,theta).*aux_force_R3(w,x,y,psi,theta,phi) + ...
    cos(beta(4)).*interp_fun_px_TM_R3(-phi,-theta).*aux_force_R4(w,x,y,psi,theta,phi));
    
fy_BF = @(w,x,y,psi,theta,phi) cos(theta).*cos(phi).*I0.*(c0./(I0.*D.^2)).*(...
    sin(beta(1)).*interp_fun_px_TE_R1(phi,theta).*aux_force_R1(w,x,y,psi,theta,phi) + ...
    sin(beta(2)).*interp_fun_px_TE_R1(-phi,-theta).*aux_force_R2(w,x,y,psi,theta,phi) + ...
    sin(beta(3)).*interp_fun_px_TM_R3(phi,theta).*aux_force_R3(w,x,y,psi,theta,phi) + ...
    sin(beta(4)).*interp_fun_px_TM_R3(-phi,-theta).*aux_force_R4(w,x,y,psi,theta,phi));
    
fz_BF = @(w,x,y,psi,theta,phi) cos(theta).*cos(phi).*I0.*(c0./(I0.*D.^2)).*(...
    interp_fun_pz_TE_R1(phi,theta).*aux_force_R1(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TE_R1(-phi,-theta).*aux_force_R2(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TM_R3(phi,theta).*aux_force_R3(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TM_R3(-phi,-theta).*aux_force_R4(w,x,y,psi,theta,phi));

tx = @(w,x,y,psi,theta,phi) cos(theta).*cos(phi).*I0.*(c0./(I0.*D.^2)).*(...
    interp_fun_pz_TE_R1(phi,theta).*aux_torqueX_R1(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TE_R1(-phi,-theta).*aux_torqueX_R2(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TM_R3(phi,theta).*aux_torqueX_R3(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TM_R3(-phi,-theta).*aux_torqueX_R4(w,x,y,psi,theta,phi) );

ty = @(w,x,y,psi,theta,phi) -cos(theta).*cos(phi).*I0.*(c0./(I0.*D.^2)).*(...
    interp_fun_pz_TE_R1(phi,theta).*aux_torqueY_R1(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TE_R1(-phi,-theta).*aux_torqueY_R2(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TM_R3(phi,theta).*aux_torqueY_R3(w,x,y,psi,theta,phi) + ...
    interp_fun_pz_TM_R3(-phi,-theta).*aux_torqueY_R4(w,x,y,psi,theta,phi) );

tz = @(w,x,y,psi,theta,phi) cos(theta).*cos(phi).*I0.*(c0./(I0.*D.^2)).*( ...
    (sin(beta(1)).*interp_fun_px_TE_R1(phi,theta).*aux_torqueY_R1(w,x,y,psi,theta,phi) - ...
    cos(beta(1)).*interp_fun_px_TE_R1(phi,theta).*aux_torqueX_R1(w,x,y,psi,theta,phi)) + ...
    (sin(beta(2)).*interp_fun_px_TE_R1(-phi,-theta).*aux_torqueY_R2(w,x,y,psi,theta,phi) - ...
    cos(beta(2)).*interp_fun_px_TE_R1(-phi,-theta).*aux_torqueX_R2(w,x,y,psi,theta,phi)) + ...
    (sin(beta(3)).*interp_fun_px_TM_R3(phi,theta).*aux_torqueY_R3(w,x,y,psi,theta,phi) - ...
    cos(beta(3)).*interp_fun_px_TM_R3(phi,theta).*aux_torqueX_R3(w,x,y,psi,theta,phi)) + ...
    (sin(beta(4)).*interp_fun_px_TM_R3(-phi,-theta).*aux_torqueY_R4(w,x,y,psi,theta,phi) - ...
    cos(beta(4)).*interp_fun_px_TM_R3(-phi,-theta).*aux_torqueX_R4(w,x,y,psi,theta,phi)) );

%% Define elements of direction cosine matrix

DRM123_11 = @(psi,theta,phi) cos(psi).*cos(theta);
DRM123_12 = @(psi,theta,phi) cos(psi)*sin(theta).*sin(phi) + sin(psi).*cos(phi);
DRM123_13 = @(psi,theta,phi) -cos(psi).*sin(theta).*cos(phi) + sin(psi).*sin(phi);
DRM123_21 = @(psi,theta,phi) -sin(psi).*cos(theta);
DRM123_22 = @(psi,theta,phi) -sin(psi).*sin(theta).*sin(phi) + cos(psi).*cos(phi);
DRM123_23 = @(psi,theta,phi) sin(psi)*sin(theta).*cos(phi) + cos(psi).*sin(phi);
DRM123_31 = @(psi,theta,phi) sin(theta);
DRM123_32 = @(psi,theta,phi) -cos(theta).*sin(phi);
DRM123_33 = @(psi,theta,phi) cos(theta).*cos(phi);

%% Define lab-frame forces

fx = @(w,x,y,psi,theta,phi) ...
    fx_BF(w,x,y,psi,theta,phi) .* DRM123_11(psi,theta,phi) + ...
    fy_BF(w,x,y,psi,theta,phi) .* DRM123_21(psi,theta,phi) + ...
    fz_BF(w,x,y,psi,theta,phi) .* DRM123_31(psi,theta,phi);

fy = @(w,x,y,psi,theta,phi) ...
    fx_BF(w,x,y,psi,theta,phi) .* DRM123_12(psi,theta,phi) + ...
    fy_BF(w,x,y,psi,theta,phi) .* DRM123_22(psi,theta,phi) + ...
    fz_BF(w,x,y,psi,theta,phi) .* DRM123_32(psi,theta,phi);

%% Define "total" torques

fphi = @(w,x,y,psi,theta,phi,wx,wy,wz) ((Iy - Iz).*wy.*wz + tx(w,x,y,psi,theta,phi))/Ix;
ftheta = @(w,x,y,psi,theta,phi,wx,wy,wz) ((Iz - Ix).*wx.*wz + ty(w,x,y,psi,theta,phi))/Iy;
fpsi = @(w,x,y,psi,theta,phi,wx,wy,wz) ((Ix - Iy).*wx.*wy + tz(w,x,y,psi,theta,phi))/Iz;

%% Define equilibrium

x_eq = 0;
y_eq = 0;
psi_eq = 0;
theta_eq = 0;
phi_eq = 0;
wx_eq = 0;
wy_eq = 0;
wz_eq = 1;

% Calculate time constant
t0 = sqrt((mass*(D_actual^2))*c0/I0_actual/D_actual);

% Calculate normalized angular frequency equivalent to 
omega_z_0 = 2*pi*spinning_freq*t0;

psi = deg2rad(0);
wz = omega_z_0;

%% Calculate partial derivatives

diff_step = 1E-6;

% Simple numeric derivatives
dfdx = @(f,w,x,y,psi,theta,phi) (f(w,x + diff_step,y,psi,theta,phi) - f(w,x,y,psi,theta,phi)) ./ diff_step;
dfdy = @(f,w,x,y,psi,theta,phi) (f(w,x,y + diff_step,psi,theta,phi) - f(w,x,y,psi,theta,phi)) ./ diff_step;
dfdtheta = @(f,w,x,y,psi,theta,phi) (f(w,x,y,psi,theta + diff_step,phi) - f(w,x,y,psi,theta,phi)) ./ diff_step;
dfdphi = @(f,w,x,y,psi,theta,phi) (f(w,x,y,psi,theta,phi + diff_step) - f(w,x,y,psi,theta,phi)) ./ diff_step;
dfdpsi = @(f,w,x,y,psi,theta,phi) (f(w,x,y,psi + diff_step,theta,phi) - f(w,x,y,psi,theta,phi)) ./ diff_step;

dfadx = @(f,w,x,y,psi,theta,phi,wx,wy,wz) (f(w,x + diff_step,y,psi,theta,phi,wx,wy,wz) - f(w,x,y,psi,theta,phi,wx,wy,wz)) ./ diff_step;
dfady = @(f,w,x,y,psi,theta,phi,wx,wy,wz) (f(w,x,y + diff_step,psi,theta,phi,wx,wy,wz) - f(w,x,y,psi,theta,phi,wx,wy,wz)) ./ diff_step;
dfadtheta = @(f,w,x,y,psi,theta,phi,wx,wy,wz) (f(w,x,y,psi,theta + diff_step,phi,wx,wy,wz) - f(w,x,y,psi,theta,phi,wx,wy,wz)) ./ diff_step;
dfadphi = @(f,w,x,y,psi,theta,phi,wx,wy,wz) (f(w,x,y,psi,theta,phi + diff_step,wx,wy,wz) - f(w,x,y,psi,theta,phi,wx,wy,wz)) ./ diff_step;
dfadpsi = @(f,w,x,y,psi,theta,phi,wx,wy,wz) (f(w,x,y,psi + diff_step,theta,phi,wx,wy,wz) - f(w,x,y,psi,theta,phi,wx,wy,wz)) ./ diff_step;

% fxx = dfdx(fx,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fxy = dfdy(fx,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fxtheta = dfdtheta(fx,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fxphi = dfdphi(fx,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fxpsi = dfdpsi(fx,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq);
% 
% fyx = dfdx(fy,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fyy = dfdy(fy,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fytheta = dfdtheta(fy,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fyphi = dfdphi(fy,w,x_eq,y_eq,psi,theta_eq,phi_eq);
% fypsi = dfdpsi(fy,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq);
% 
% fthetax = dfadx(ftheta,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fthetay = dfady(ftheta,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fthetatheta = dfadtheta(ftheta,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fthetaphi = dfadphi(ftheta,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fthetapsi = dfadpsi(ftheta,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
% 
% fphix = dfadx(fphi,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fphiy = dfady(fphi,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fphitheta = dfadtheta(fphi,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fphiphi = dfadphi(fphi,w,x_eq,y_eq,psi,theta_eq,phi_eq,wx_eq,wy_eq,wz);
% fphipsi = dfadpsi(fphi,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
% 
% fpsix = dfadx(fpsi,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
% fpsiy = dfady(fpsi,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
% fpsitheta = dfadtheta(fpsi,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
% fpsiphi = dfadphi(fpsi,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
% fpsipsi = dfadpsi(fpsi,w,x_eq,y_eq,psi_eq,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);

fxx = @(t) dfdx(fx,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fxy = @(t) dfdy(fx,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fxtheta = @(t) dfdtheta(fx,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fxphi = @(t) dfdphi(fx,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fxpsi = @(t) dfdpsi(fx,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);

fyx = @(t) dfdx(fy,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fyy = @(t) dfdy(fy,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fytheta = @(t) dfdtheta(fy,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fyphi = @(t) dfdphi(fy,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);
fypsi = @(t) dfdpsi(fy,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq);

fthetax = @(t) dfadx(ftheta,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fthetay = @(t) dfady(ftheta,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fthetatheta = @(t) dfadtheta(ftheta,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fthetaphi = @(t) dfadphi(ftheta,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fthetapsi = @(t) dfadpsi(ftheta,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);

fphix = @(t) dfadx(fphi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fphiy = @(t) dfady(fphi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fphitheta = @(t) dfadtheta(fphi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fphiphi = @(t) dfadphi(fphi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz);
fphipsi = @(t) dfadpsi(fphi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);

fpsix = @(t) dfadx(fpsi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
fpsiy = @(t) dfady(fpsi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
fpsitheta = @(t) dfadtheta(fpsi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
fpsiphi = @(t) dfadphi(fpsi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);
fpsipsi = @(t) dfadpsi(fpsi,w,x_eq,y_eq,wz.*t,theta_eq,phi_eq,wx_eq,wy_eq,wz_eq);

%% Calculate state transition matrix

% Define initial conditions being the identity matrix
I8 = eye(8);
u0 = I8(:);

% Define vectorial form of differential equation for Phi' = J.Phi
f = @(t,u) [ ...
    u(17);
    u(18);
    u(19);
    u(20);
    u(21);
    u(22);
    u(23);
    u(24);
    u(25);
    u(26);
    u(27);
    u(28);
    u(29);
    u(30);
    u(31);
    u(32);
    u(1).*fxx(t) + u(9).*fxy(t) + u(33).*fxtheta(t) + u(41).*fxphi(t);
    u(2).*fxx(t) + u(10).*fxy(t) + u(34).*fxtheta(t) + u(42).*fxphi(t);
    u(3).*fxx(t) + u(11).*fxy(t) + u(35).*fxtheta(t) + u(43).*fxphi(t);
    u(4).*fxx(t) + u(12).*fxy(t) + u(36).*fxtheta(t) + u(44).*fxphi(t);
    u(5).*fxx(t) + u(13).*fxy(t) + u(37).*fxtheta(t) + u(45).*fxphi(t);
    u(6).*fxx(t) + u(14).*fxy(t) + u(38).*fxtheta(t) + u(46).*fxphi(t);
    u(7).*fxx(t) + u(15).*fxy(t) + u(39).*fxtheta(t) + u(47).*fxphi(t);
    u(8).*fxx(t) + u(16).*fxy(t) + u(40).*fxtheta(t) + u(48).*fxphi(t);
    u(1).*fyx(t) + u(9).*fyy(t) + u(33).*fytheta(t) + u(41).*fyphi(t);
    u(2).*fyx(t) + u(10).*fyy(t) + u(34).*fytheta(t) + u(42).*fyphi(t);
    u(3).*fyx(t) + u(11).*fyy(t) + u(35).*fytheta(t) + u(43).*fyphi(t);
    u(4).*fyx(t) + u(12).*fyy(t) + u(36).*fytheta(t) + u(44).*fyphi(t);
    u(5).*fyx(t) + u(13).*fyy(t) + u(37).*fytheta(t) + u(45).*fyphi(t);
    u(6).*fyx(t) + u(14).*fyy(t) + u(38).*fytheta(t) + u(46).*fyphi(t);
    u(7).*fyx(t) + u(15).*fyy(t) + u(39).*fytheta(t) + u(47).*fyphi(t);
    u(8).*fyx(t) + u(16).*fyy(t) + u(40).*fytheta(t) + u(48).*fyphi(t);
    u(57).*cos(wz.*t) + u(49).*sin(wz.*t);
    u(58).*cos(wz.*t) + u(50).*sin(wz.*t);
    u(59).*cos(wz.*t) + u(51).*sin(wz.*t);
    u(60).*cos(wz.*t) + u(52).*sin(wz.*t);
    u(61).*cos(wz.*t) + u(53).*sin(wz.*t);
    u(62).*cos(wz.*t) + u(54).*sin(wz.*t);
    u(63).*cos(wz.*t) + u(55).*sin(wz.*t);
    u(64).*cos(wz.*t) + u(56).*sin(wz.*t);
    u(49).*cos(wz.*t) - u(57).*sin(wz.*t);
    u(50).*cos(wz.*t) - u(58).*sin(wz.*t);
    u(51).*cos(wz.*t) - u(59).*sin(wz.*t);
    u(52).*cos(wz.*t) - u(60).*sin(wz.*t);
    u(53).*cos(wz.*t) - u(61).*sin(wz.*t);
    u(54).*cos(wz.*t) - u(62).*sin(wz.*t);
    u(55).*cos(wz.*t) - u(63).*sin(wz.*t);
    u(56).*cos(wz.*t) - u(64).*sin(wz.*t);
    u(1).*fphix(t) + u(9).*fphiy(t) + u(33).*fphitheta(t) + u(41).*fphiphi(t) - u(57).*wz;
    u(2).*fphix(t) + u(10).*fphiy(t) + u(34).*fphitheta(t) + u(42).*fphiphi(t) - u(58).*wz;
    u(3).*fphix(t) + u(11).*fphiy(t) + u(35).*fphitheta(t) + u(43).*fphiphi(t) - u(59).*wz;
    u(4).*fphix(t) + u(12).*fphiy(t) + u(36).*fphitheta(t) + u(44).*fphiphi(t) - u(60).*wz;
    u(5).*fphix(t) + u(13).*fphiy(t) + u(37).*fphitheta(t) + u(45).*fphiphi(t) - u(61).*wz;
    u(6).*fphix(t) + u(14).*fphiy(t) + u(38).*fphitheta(t) + u(46).*fphiphi(t) - u(62).*wz;
    u(7).*fphix(t) + u(15).*fphiy(t) + u(39).*fphitheta(t) + u(47).*fphiphi(t) - u(63).*wz;
    u(8).*fphix(t) + u(16).*fphiy(t) + u(40).*fphitheta(t) + u(48).*fphiphi(t) - u(64).*wz;
    u(1).*fthetax(t) + u(9).*fthetay(t) + u(33).*fthetatheta(t) + u(41).*fthetaphi(t) + u(49).*wz;
    u(2).*fthetax(t) + u(10).*fthetay(t) + u(34).*fthetatheta(t) + u(42).*fthetaphi(t) + u(50).*wz;
    u(3).*fthetax(t) + u(11).*fthetay(t) + u(35).*fthetatheta(t) + u(43).*fthetaphi(t) + u(51).*wz;
    u(4).*fthetax(t) + u(12).*fthetay(t) + u(36).*fthetatheta(t) + u(44).*fthetaphi(t) + u(52).*wz;
    u(5).*fthetax(t) + u(13).*fthetay(t) + u(37).*fthetatheta(t) + u(45).*fthetaphi(t) + u(53).*wz;
    u(6).*fthetax(t) + u(14).*fthetay(t) + u(38).*fthetatheta(t) + u(46).*fthetaphi(t) + u(54).*wz;
    u(7).*fthetax(t) + u(15).*fthetay(t) + u(39).*fthetatheta(t) + u(47).*fthetaphi(t) + u(55).*wz;
    u(8).*fthetax(t) + u(16).*fthetay(t) + u(40).*fthetatheta(t) + u(48).*fthetaphi(t) + u(56).*wz ...
    ];

% Define interval for numerical integration
tstart = 0;
tend = 1/(spinning_freq*t0);

% Numerically integrate differential equations
tic
[t,u] = ode45(f,[tstart tend],u0);
toc

% Date format for file name
formatOut = 'mmddyy';

% Define file base name
save_name = ['MEPS_',datestr(now,formatOut),'_3D_Si3N4_M1e_1064nm_',num2str(round(w,1)),'D_',num2str(spinning_freq),'Hz'];
    
% Save results
save(fullfile('Data',[save_name,'_STM_0_to_T.mat']),'u');

% Load saved results
load(fullfile('Data',[save_name,'_STM_0_to_T.mat']));

% State transition matrix at t = T
transposedPhi = reshape(u(end,:),8,8);
Phi = transposedPhi';

abs(eig(Phi))
