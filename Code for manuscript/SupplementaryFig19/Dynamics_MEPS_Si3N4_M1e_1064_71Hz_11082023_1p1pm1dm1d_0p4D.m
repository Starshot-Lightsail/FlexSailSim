%% Numerical evolution of coupled differential equations for rigid and
% composite round lightsails, Mark 1e, Si3N4 at 1064 nm
% Equations of motions of 2nd order are rewritten in vector form
% Simulate 5 s for paper with initial conditions:
% x0 = 0.025D, y0 = 0.025D, theta0 = -1 deg, phi0 = -1 deg, w = 0.4D


% Change (02/11/22): Add power ramp function from flexible lightsail code
% to equations of motion

% Change (02/10/22): Spin the opposite direction
% Change (02/05/22): Changed geometrical parameters to correctly represent
% unit cells of Mark 1e design

% Change (07/04/21): Corrected mistake of calculating position vector
% in argument of intensity exponetial distribution, i.e., rc + H_B^I * rb
% instead of rc + H_I^B * rb. This has the consequence that the correct
% transformation of forces from body frame to lab frame is given by the
% intuitively logical formula: F = H_B^I * F^{BF}

% Ramon Gao, October 12, 2023

%% Initialization
% Close all figures, clear command window, clear all stored variables

close all
clear
clc
clear textprogressbar

% Evolve equations of motion numerically
do_simulate = 1;

% Save simulation results (1) or not (0)
dosave = 1;

% Set material of lightsail to be 'SOI', 'SiNx2' or 'Si3N4'
material = 'Si3N4';

% Sail shape: string, 'round' or 'square'
shape = 'round';

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
if strcmpi(material,'SOI')
    w = 2*D;    % [w] = 1, normalized beam width
elseif strcmpi(material,'SiNx2')
    w = 0.1*D;
elseif strcmpi(material,'Si3N4')
    w = 0.4*D;
end

% Set spinning frequency
spinning_freq = 71; % in Hz

% Set initial conditions
x0 = 0.01*D;
y0 = 0.01*D;
theta0 = deg2rad(-1);
phi0 = deg2rad(-1);

% Set simulation duraction
tend = (2/5) .* 310.12028577; % in units of t0, or: normalized by t0

%% General parameters and settings

% Densities
rho_Si3N4 = 2700*1e-6; % [rho] = kg/m^2/um, density of high-stress Si3N4
rho_SiNx = 3100*1e-6; % [rho] = kg/m^2/um, density of low-stress SiNx
rho_SiO2 = 2200*1e-6; % [rho] = kg/m^2/um, density of silica
rho_Si = 2300*1e-6; % [rho] = kg/m^2/um, density of silicon

% Orientation of gratings in their local frames
beta = [0, pi, pi/2, 3*pi/2];

% Position of boundary separating unpatterned and patterned regions
if strcmpi(shape,'round')
    sx = 0*1/10; % [sx] = 1, unitless, actual position is [sx*D] = m
elseif strcmpi(shape,'square')
    sx = 1/4; % [sx] = 1, unitless, actual position is [sx*D] = m
end

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
    
elseif strcmpi(material,'Si3N4')
    
    % Define density of substrate and resonating elements
    rho_s = rho_Si3N4; % [rho] = kg/m^2/um
    rho_e = rho_Si3N4; % [rho] = kg/m^2/um

    % Define design parameters of unit cells A and B 
    w1A = 0.6; % [w1A] = um, width of first resonator in region A
    w2A = 0.2; % [w2A] = um, width of second resonator in region A
    w1B = 0.52; % [w1B] = um, width of first resonator in region B
    w2B = 0.2; % [w2B] = um, width of second resonator in region B
    
    dUC_A = 1.6; % [dUC_A] = um, period of unit cell of region A
    dUC_B = 1.35; % [dUC_B] = um, period of unit cell of region B

    tA = 0.2; % [ts] = um, thickness of substrate
    hA = 0.4; % [hr] = um, height of resonators
    tB = tA; % [ts] = um, thickness of substrate
    hB = hA; % [hr] = um, height of resonators

    % Define the mass per unit area of each region (unpatterned, A or B)
    chi0 = rho_s*tA;
    chiA = rho_s*tA + rho_e*hA*(w1A + w2A)/dUC_A;
    chiB = rho_s*tB + rho_e*hB*(w1B + w2B)/dUC_B;
    
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
        
end

    
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

% Define prefactors for moments of inertias to be multiplied with m*D^2
prefactor_Ix = Ix_eval / mass;
prefactor_Iy = Iy_eval / mass;
prefactor_Iz = prefactor_Ix + prefactor_Iy;

% Define moments of inertia in normalzed units, i.e. in terms of m and D
global Ix
global Iy
global Iz
Ix = prefactor_Ix .* m .* D.^2;
Iy = prefactor_Iy .* m .* D.^2;
Iz = prefactor_Iz .* m .* D.^2;

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

if strcmpi(material,'SOI')

    file_name_collection_small_tables = ...
        '04032021_MEPS_SOI_Oggy_ReducedPitchRoll_Pressures_CollectionTables_';
    
    load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
    load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])                    

elseif strcmpi(material,'Si3N4')
    
    file_name_collection_small_tables = ...
        '08232021_MEPS_Si3N4_Mark1e_1064_PitchRoll_Pressures_CollectionTables_';

    load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
    load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])
        
end

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
    
integrandTorqueX = @(r,varphi,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*((xc+r.*cos(psi).*cos(theta).*cos(varphi)+ ...
    (-1).*r.*cos(theta).*sin(psi).*sin(varphi)).^2+(r.*cos(varphi).*( ...
    sin(phi).*sin(psi)+(-1).*cos(phi).*cos(psi).*sin(theta))+r.*(cos( ...
    psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*sin(varphi)).^2+( ...
    yc+r.*cos(varphi).*(cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin( ...
    theta))+r.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin( ...
    theta)).*sin(varphi)).^2)).*r.^2.*sin(varphi);

integrandTorqueY = @(r,varphi,xc,yc,psi,theta,phi) ...
    exp(1).^((-2).*w.^(-2).*((xc+r.*cos(psi).*cos(theta).*cos(varphi)+ ...
    (-1).*r.*cos(theta).*sin(psi).*sin(varphi)).^2+(r.*cos(varphi).*( ...
    sin(phi).*sin(psi)+(-1).*cos(phi).*cos(psi).*sin(theta))+r.*(cos( ...
    psi).*sin(phi)+cos(phi).*sin(psi).*sin(theta)).*sin(varphi)).^2+( ...
    yc+r.*cos(varphi).*(cos(phi).*sin(psi)+cos(psi).*sin(phi).*sin( ...
    theta))+r.*(cos(phi).*cos(psi)+(-1).*sin(phi).*sin(psi).*sin( ...
    theta)).*sin(varphi)).^2)).*r.^2.*cos(varphi);

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

funtorquex = @(xc,yc,psi,theta,phi) cos(theta).*cos(phi).*I0.* ...
    (c0./(I0.*D.^2)).* ...
    (findpressure('pz_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1) + ...
    findpressure('pz_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2) + ...
    findpressure('pz_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3) + ...
    findpressure('pz_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4) );

funtorquey = @(xc,yc,psi,theta,phi) -cos(theta).*cos(phi).*I0.* ...
    (c0/(I0.*D.^2)).* ...
    (findpressure('pz_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1) + ...
    findpressure('pz_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2) + ...
    findpressure('pz_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3) + ...
    findpressure('pz_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4) );

funtorquez = @(xc,yc,psi,theta,phi) cos(theta).*cos(phi).*I0.* ...
    (c0/(I0.*D.^2)).* ...
    ((sin(beta(1)).*findpressure('px_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1) - ...
    cos(beta(1)).*findpressure('px_R1',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R1,r_stop_R1,varphi_start_R1,varphi_stop_R1)) + ...
    (sin(beta(2)).*findpressure('px_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2) - ...
    cos(beta(2)).*findpressure('px_R2',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R2,r_stop_R2,varphi_start_R2,varphi_stop_R2)) + ...
    (sin(beta(3)).*findpressure('px_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3) - ...
    cos(beta(3)).*findpressure('px_R3',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R3,r_stop_R3,varphi_start_R3,varphi_stop_R3)) + ...
    (sin(beta(4)).*findpressure('px_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueY(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4) - ...
    cos(beta(4)).*findpressure('px_R4',[phi,theta,psi]).* ...
    integral2(@(r,varphi) integrandTorqueX(r,varphi,xc,yc,psi,theta,phi), ...
    r_start_R4,r_stop_R4,varphi_start_R4,varphi_stop_R4)) );

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

%% Main code

% Calculate time constant
t0 = sqrt((mass*(D_actual^2))*c0/I0_actual/D_actual);

% Calculate normalized angular frequency equivalent to 
omega_z_0 = 2*pi*spinning_freq*t0;

% Power ramp parameters
global t_on_dur
global t_ramp_delay
global t_ramp_dur
t_on_dur = tend;
t_ramp_delay = 0/t0;
t_ramp_dur = 1e-6/t0;

% Define function that describes the equations of the system in vector form

f = @(t,u) [u(7);
    u(8);
    u(9);
    -cos(u(4)).*tan(u(5)).*u(10) + sin(u(4)).*tan(u(5)).*u(11) + u(12);
    sin(u(4)).*u(10) + cos(u(4)).*u(11);
    cos(u(4)).*sec(u(5)).*u(10) - sin(u(4)).*sec(u(5)).*u(11);
    Fx(u(1),u(2),u(4),u(5),u(6));
    Fy(u(1),u(2),u(4),u(5),u(6));
    Fz(u(1),u(2),u(4),u(5),u(6));
    ((Iy - Iz).*u(11).*u(12) + funtorquex(u(1),u(2),u(4),u(5),u(6)))./Ix;
    ((Iz - Ix).*u(10).*u(12) + funtorquey(u(1),u(2),u(4),u(5),u(6)))./Iy;
    ((Ix - Iy).*u(10).*u(11) + funtorquez(u(1),u(2),u(4),u(5),u(6)))./Iz;
    ];    

% Define initial conditions according to
% '1' : moderate initial condition
% '2' : extreme initial condition
% '3' : sweep over different initial conditions

initialcond = '1';

if initialcond == '1'
    u0 = [x0; y0; 0; 0; theta0; phi0; 0; 0; 0; 0; 0; -omega_z_0];
elseif initialcond == '2' 
    
    numInitialConditions = 3;
    
    u0mat = zeros(numInitialConditions,12);
    
    u0mat(1,:) = [0.02*D; 0.02*D; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2*pi*1*t0];
    u0mat(2,:) = [0.05*D; 0.05*D; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2*pi*1*t0];
    u0mat(3,:) = [0.05*D; 0.05*D; 0; 0; 0; 0; 0; 0; 0; 0; 0; 2*pi*50*t0];

end

tstart = 0;

if initialcond == '1'
    tend = tend;
elseif initialcond == '2'
    tend = 1;
end 

%% Calculate force-torque pairs as in figure 4c in [OIlic2019]

if 0
    
theta_list = -15:1:15;
theta_list = deg2rad(theta_list);
phi_list = theta_list;

fx = zeros(size(theta_list));
fy = zeros(size(theta_list));
tx = zeros(size(theta_list));
ty = zeros(size(theta_list));

for ii = 1:1:length(theta_list)
    fx(ii) = cos(theta_list(ii)) .* FxBodyFrame(0,0,0,theta_list(ii),0) + ...
        sin(theta_list(ii)) .* FzBodyFrame(0,0,0,theta_list(ii),0);
    fy(ii) = cos(phi_list(ii)) .* FyBodyFrame(0,0,0,0,phi_list(ii)) - ...
        sin(phi_list(ii)) .* FzBodyFrame(0,0,0,0,phi_list(ii));
    tx(ii) = funtorquex(0,0,0,0,phi_list(ii));
    ty(ii) = funtorquey(0,0,0,theta_list(ii),0);
end

figure(1)
plot(100*tan(theta_list),fx)
hold on
plot(100*tan(theta_list),10*ty)
xlim([-21, 21])
% ylim([-0.62, 0.62])

figure(2)
plot(100*tan(theta_list),fy)
hold on
plot(100*tan(theta_list),10*tx)
xlim([-21, 21])
% ylim([-0.62, 0.62])

end

%% Single Simulation
if do_simulate
    
rel_tol = 1e-8;
abs_tol = rel_tol;

% options = odeset('MaxStep',0.01,'OutputFcn',@odeprog,'Events',@odeabort);
options = odeset('RelTol',rel_tol,'AbsTol',abs_tol,'OutputFcn',@odetpbar);
% options = odeset('OutputFcn',@odetpbar);
% options = odeset('OutputFcn',@odetpbar);

tic
% [t,u] = ode45(f,[tstart tend],u0,options);
[t,u] = ode45(@(t,u) odefunMEPS(t,u,Fx,Fy,Fz,funtorquex,funtorquey,funtorquez), ...
    [tstart tend],u0,options);
toc

linewidth = 2;

figure(1)
% plot(t(:,1),100.*tan(u(:,5)),'Linewidth',linewidth)
% title('tan(\theta)')
plot(t(:,1),rad2deg(u(:,5)),'Linewidth',linewidth)
xlabel('\it{t}/\it{t}_{0}')
ylabel('\theta (°)')
title('Pitch')

figure(2)
plot(t(:,1),rad2deg(u(:,6)),'Linewidth',linewidth)
xlabel('\it{t}/\it{t}_{0}')
ylabel('\phi (°)')
title('Roll')
 
% figure(4)
% plot(t(:,1),rad2deg(mod(u(:,4),2*pi)),'Linewidth',linewidth)
% xlabel('\it{t}/\it{t}_{0}')
% ylabel('\psi (°)')
% title('Yaw')

tstepmax = max(t(2:end) - t(1:end-1));
fprintf('Maximum time step taken by ode: %f\n',tstepmax);

statevec = [ t u];

% Date format for file name
formatOut = 'mmddyy';

% Define file base name
save_name = ['Si3N4_',datestr(now,formatOut),'_M1e_1064nm_',num2str(100*x0),'p', ...
    num2str(100*y0),'p',num2str(rad2deg(theta0)),'d',num2str(rad2deg(phi0)), ...
    'd_',num2str(spinning_freq),'Hz_t',num2str(round(tend,2)),'_w',num2str(round(w,1)),'D'];

% Replace '-' with 'm' and '.' with 'p'
save_name = strrep(save_name,'-','m');
save_name = strrep(save_name,'.','p');

if initialcond == '1' && dosave == 1
    save(fullfile('Data',[save_name,'_states.mat']),'statevec');
elseif initialcond == '2' && dosave == 1
    save('MEPS_02232021_SiNx_Mark3_3D_5p5p_50Hz_states.mat','statevec');
end

figure(3)
plot(u(:,1),u(:,2),'Linewidth',linewidth)
xlabel('\it{x/D}')
ylabel('\it{y/D}')
% axis([-1.5*max(u0(1),u0(2)) 1.5*max(u0(1),u0(2)) ...
%     -1.5*max(u0(1),u0(2)) 1.5*max(u0(1),u0(2))])
title(['\it{f} = ',num2str(spinning_freq),' Hz, \it{x}_0 = ',num2str(x0), ...
    '\it{D}, \it{y}_0 = ', num2str(y0),'\it{D}, \theta = ', ...
    num2str(rad2deg(theta0)),'°, \phi = ',num2str(rad2deg(phi0)),'°'])

saveas(gcf,fullfile('Figures',[save_name,'_traj.png']));


end

%% Auxiliary functions

function dudt = odefunMEPS(t,u,Fx,Fy,Fz,funtorquex,funtorquey,funtorquez)
% Time-dependent equations of motion for rigid lightsail, taking into
% account a power ramp of incident laser power

global t_ramp_delay
global t_ramp_dur
global t_on_dur

global Ix
global Iy
global Iz

pwrRampFcn = @(tt) 0.5*double(tt > t_ramp_delay).*double(tt <= t_ramp_delay + t_ramp_dur) .* ...
    (1-cos(pi*(tt - t_ramp_delay)./(t_ramp_dur))) + double(tt > t_ramp_delay + t_ramp_dur) .* ...
    double(tt <= t_ramp_delay + t_ramp_dur + t_on_dur) + ...
    0.5*double(tt > t_ramp_delay + t_ramp_dur + t_on_dur) .* ...
    double(tt <= t_ramp_delay + t_ramp_dur + t_on_dur + t_ramp_dur) .* ...
    (1 + cos(pi*(tt - t_ramp_delay - t_ramp_dur - t_on_dur)./t_ramp_dur));

dudt = zeros(12,1);

dudt(1) = u(7);
dudt(2) = u(8);
dudt(3) = u(9);
dudt(4) = -cos(u(4)).*tan(u(5)).*u(10) + sin(u(4)).*tan(u(5)).*u(11) + u(12);
dudt(5) = sin(u(4)).*u(10) + cos(u(4)).*u(11);
dudt(6) = cos(u(4)).*sec(u(5)).*u(10) - sin(u(4)).*sec(u(5)).*u(11);
dudt(7) = pwrRampFcn(t).*Fx(u(1),u(2),u(4),u(5),u(6));
dudt(8) = pwrRampFcn(t).*Fy(u(1),u(2),u(4),u(5),u(6));
dudt(9) = pwrRampFcn(t).*Fz(u(1),u(2),u(4),u(5),u(6));
dudt(10) = ((Iy - Iz).*u(11).*u(12) + pwrRampFcn(t).*funtorquex(u(1),u(2),u(4),u(5),u(6)))./Ix;
dudt(11) = ((Iz - Ix).*u(10).*u(12) + pwrRampFcn(t).*funtorquey(u(1),u(2),u(4),u(5),u(6)))./Iy;
dudt(12) = ((Ix - Iy).*u(10).*u(11) + pwrRampFcn(t).*funtorquez(u(1),u(2),u(4),u(5),u(6)))./Iz;

end


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