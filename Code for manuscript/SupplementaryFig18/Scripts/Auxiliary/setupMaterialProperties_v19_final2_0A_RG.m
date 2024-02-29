%% setupMaterialProperties_v17
% (c) 2018-2021 Michael Kelzenberg, Ramon Gao -- California Institute of Technology
% Sets up material properties variables... oh, and also physical constants used in simulations
% 
% These properties are set by this subscript in hopes of ensuring better consistency between various versions of the 
% simulator code.

global MATERIALS

%% Physical constants
% Residents of dimension C137 shouldn't need to change these values:
c0 = 299792458; % m/sec, speed of light.
c0mm = c0*1000; % mm/sec, speed of light
SBC = 5.67e-8; % W/m^2/K^4; boltzmann constant, for radiative cooling.
SBCmm = SBC / 1e6; % W/mm^2/K^4

%% Material selector
% String variable helps us select between different sets of properties and behaviors throughout the simulations.
%material = 'silicon' ;

disp(['Setting up physical constants and material properties.' ]);

MATERIALS = [];

MATERIALS(1).name = 'Si43nm';
MATERIALS(1).Iabs = 1e-6; % absorptivity
MATERIALS(1).Irefl = 0.5; % reflectivity
MATERIALS(1).Emissivity = 0.8; % emissivity
    %lambda = 1550e-6; %mm, wavelength
MATERIALS(1).ior_n = 3.52; % refractive index
MATERIALS(1).ior_k = 4e-5;
MATERIALS(1).density = 2.3e-3; % g/mm3
MATERIALS(1).Ymod = 169e9;  % Pa, N/m^2
MATERIALS(1).Emod = 1.20e9;  % Pa, tensile strength
MATERIALS(1).nu = 0.25; % Poisson's ratio
MATERIALS(1).Thermcondmm = 0.157; % W/mm/degC
MATERIALS(1).heatCap = 0.7; % J/g/degC
MATERIALS(1).CTE = 2.33; % PPM/K, coefficient of thermal expansion
MATERIALS(1).failtemp = 1300; % Kelvin, failure temperature
MATERIALS(1).thickness = 43e-6; % mm

MATERIALS(2).name = 'SiNx_0521';
%I0 = 1000;%1000; % W/mm^2, peak incident power density (1000 = 1GW/m2)  %MK: it goes faster at 10000
MATERIALS(2).Iabs = 1e-5; % absorptivity
MATERIALS(2).Irefl = 0.35; % reflectivity
MATERIALS(2).Emissivity = 0.1; % emissivity
%lambda = 514e-6; % mm, wavelength
MATERIALS(2).ior_n = 2.2;  % refractive index
MATERIALS(2).ior_k = 0.001;
MATERIALS(2).density = 3.1e-3; % g/mm3
MATERIALS(2).Ymod = 270e9;  % Pa, N/m^2
%warning('Using fake value of modulus!');
%Ymod = Ymod/50;
MATERIALS(2).Emod = 5.5e9;  % Pa, tensile strength
%warning('Using fake value of tensile strength!');
%Emod = Emod/100;
MATERIALS(2).nu = 0.25; % Poisson's ratio

MATERIALS(2).Thermcondmm = 0.02; % W/mm/degC
MATERIALS(2).heatCap = 0.75; % J/g/degC
MATERIALS(2).CTE = 2.4; % PPM/K, coefficient of thermal expansion
MATERIALS(2).thickness = 176e-6; % mm
MATERIALS(2).failtemp = 1100; % Kelvin, failure temperature


MATERIALS(3).name = 'Si20nm';
MATERIALS(3).Iabs = 1e-6*20000/43; % absorptivity
MATERIALS(3).Irefl = 0.5; % reflectivity
MATERIALS(3).Emissivity = 1; % emissivity
%lambda = 1550e-6; %mm, wavelength
MATERIALS(3).ior_n = 3.52; % refractive index
MATERIALS(3).ior_k = 4e-5;
MATERIALS(3).density = 2.3e-3; % g/mm3
MATERIALS(3).Ymod = 169e9;  % Pa, N/m^2
MATERIALS(3).Emod = 1.20e9;  % Pa, tensile strength
MATERIALS(3).nu = 0.25; % Poisson's ratio
MATERIALS(3).Thermcondmm = 0.157; % W/mm/degC
MATERIALS(3).heatCap = 0.7; % J/g/degC
MATERIALS(3).CTE = 2.33; % PPM/K, coefficient of thermal expansion
MATERIALS(3).failtemp = 1300; % Kelvin, failure temperature
MATERIALS(3).thickness = 20000e-6; % mm


% This section is dedicated to calculated the thickness of a plain Si3N4
% with the same mass (density) as that of the TE or TM metagrating design.
% I figured that this is the best way to keep values up to date, as the
% unit cell designs need to be entered here!

% IMPORTANT Si3N4 MATERIAL PROPERTIES DEFINED HERE
Ymod_Si3N4 = 270e9;  % in Pa = N/m^2, [DJWilson2009]
rho_Si3N4 = 2700*1e-6; % [rho] = kg/m^2/um = g/mm^3, [DJWilson2009]
nu_Si3N4 = 0.27; % [DJWilson2009]
Emod_Si3N4 = 6.4e9; % in Pa, [AKaushik2005]
Thermcondmm_Si3N4 = 3*1e-3; % in W/mm/degC, [HFtouni2015]
heatCap_Si3N4 = 0.8; % in J/g/degC, [HFtouni2015]
CTE_Si3N4 = 2.3; % in PPM/K, [DJWilson2012_Thesis]
% IMPORTANT: Factor 2 (previously omitted) needed to account for frontside
% and backside contribution to radiated power
Emissivity_Si3N4 = 2*0.1; % not directly from literature, calculations based on [CZhang2020] would give 0.1161/0.1171 for TE/TM at 300 K
n_Si3N4 = 2; % [MKaruza2013]
k_Si3N4 = 2e-6; % [MKaruza2013]

% Densities
rho_s = rho_Si3N4; % [rho] = kg/m^2/um
rho_e = rho_Si3N4; % [rho] = kg/m^2/um
    
% Define design parameters of unit cells A and B 
w1_TE = 0.6; % [w1A] = um, width of first resonator in region A
w2_TE = 0.2; % [w2A] = um, width of second resonator in region A
w1_TM = 0.52; % [w1B] = um, width of first resonator in region B
w2_TM = 0.2; % [w2B] = um, width of second resonator in region B
    
dUC_TE = 1.6; % [dUC_A] = um, period of unit cell of region A
dUC_TM = 1.35; % [dUC_B] = um, period of unit cell of region B

t_TE = 0.2; % [ts] = um, thickness of substrate
h_TE = 0.4; % [hr] = um, height of resonators
t_TM = t_TE; % [ts] = um, thickness of substrate
h_TM = h_TE; % [hr] = um, height of resonators

% Define the mass per unit area of each region (unpatterned, A or B)
chi_TE = rho_s*t_TE + rho_e*h_TE*(w1_TE + w2_TE)/dUC_TE;
chi_TM = rho_s*t_TM + rho_e*h_TM*(w1_TM + w2_TM)/dUC_TM;
    
% Define boundaries of all regions/patterns for a round lightsail
varphi_start_R1 = -pi/6;

% Define the total mass of the round lightsail, normalized by D^2
mass = (pi/4)*( (4*abs(varphi_start_R1)/(2*pi))*(chi_TE - chi_TM) + chi_TM );

% Calculate mass of round lightsail fully patterned with TE metagratings,
% normalized by D^2
mass_allTE = (pi/4) * chi_TE;
mass_allTM = (pi/4) * chi_TM;

% Calculate thickness [m] of equivalent plain film with mass of lightsail
% having all-TE/TM gratings
radiusmm = 500; % DOUBLE CHECK WITH VALUE in generateMesh!!
t_TE = mass_allTE .* (2*radiusmm*1e3).^2 ./ (pi*(rho_Si3N4*1e6)*(radiusmm*1e3)^2 );
t_TM = mass_allTM .* (2*radiusmm*1e3).^2 ./ (pi*(rho_Si3N4*1e6)*(radiusmm*1e3)^2 );

% From thicknesses, calculate absorptivity
Iabs_Si3N4_TE = 1 - exp(-4*pi*k_Si3N4*t_TE / 1064e-9 );
Iabs_Si3N4_TM = 1 - exp(-4*pi*k_Si3N4*t_TM / 1064e-9 );

% Assign defined and calculated properties

MATERIALS(4).name = 'Si3N4_TE';
MATERIALS(4).Iabs = 0*Iabs_Si3N4_TE; % absorptivity
MATERIALS(4).Irefl = 0.15; % reflectivity
MATERIALS(4).Emissivity = 0*Emissivity_Si3N4; % emissivity
MATERIALS(4).ior_n = n_Si3N4;  % refractive index
MATERIALS(4).ior_k = k_Si3N4;
MATERIALS(4).density = rho_Si3N4; % g/mm3
MATERIALS(4).Ymod = Ymod_Si3N4;  % Pa, N/m^2
MATERIALS(4).Emod = Emod_Si3N4;  % Pa, tensile strength
MATERIALS(4).nu = nu_Si3N4; % Poisson's ratio

MATERIALS(4).Thermcondmm = Thermcondmm_Si3N4; % W/mm/degC, [HFtouni2015]
MATERIALS(4).heatCap = heatCap_Si3N4; % J/g/degC, [HFtouni2015]
MATERIALS(4).CTE = CTE_Si3N4; % PPM/K, coefficient of thermal expansion
MATERIALS(4).thickness = t_TE*1e3; % mm
MATERIALS(4).failtemp = 1400 + 273.15; % Kelvin, failure temperature

MATERIALS(5).name = 'Si3N4_TM';
MATERIALS(5).Iabs = 0*Iabs_Si3N4_TM; % absorptivity
MATERIALS(5).Irefl = 0.15; % reflectivity
MATERIALS(5).Emissivity = 0*Emissivity_Si3N4; % emissivity
MATERIALS(5).ior_n = n_Si3N4;  % refractive index
MATERIALS(5).ior_k = k_Si3N4;
MATERIALS(5).density = rho_Si3N4; % g/mm3
MATERIALS(5).Ymod = Ymod_Si3N4;  % Pa, N/m^2
MATERIALS(5).Emod = Emod_Si3N4;  % Pa, tensile strength
MATERIALS(5).nu = nu_Si3N4; % Poisson's ratio

MATERIALS(5).Thermcondmm = Thermcondmm_Si3N4; % W/mm/degC, [HFtouni2015]
MATERIALS(5).heatCap = heatCap_Si3N4; % J/g/degC, [HFtouni2015]
MATERIALS(5).CTE = CTE_Si3N4; % PPM/K, coefficient of thermal expansion
MATERIALS(5).thickness = t_TM*1e3; % mm
MATERIALS(5).failtemp = 1400 + 273.15; % Kelvin, failure temperature


% defaults:
material = 'silicon';
   %I0 = 4000; % W/mm^2, peak incident power density (1000 = 1GW/m2)
Iabs = 1e-6; % absorptivity
Irefl = 0.5; % reflectivity
Emissivity = 0.8; % emissivity
lambda = 1550e-6; %mm, wavelength
ior_n = 3.52; % refractive index
ior_k = 4e-5;
density = 2.3e-3; % g/mm3
Ymod = 169e9;  % Pa, N/m^2
Emod = 1.20e9;  % Pa, tensile strength
nu = 0.25; % Poisson's ratio

Thermcondmm = 0.157; % W/mm/degC
heatCap = 0.7; % J/g/degC
CTE = 2.33; % PPM/K, coefficient of thermal expansion
failtemp = 1300; %kelvin
thickness = 43e-6; % mm 
    
tensile_strain_limit = Emod/Ymod;
thicknessnm = thickness.*1e6;

Ymodmm = Ymod_Si3N4./1e6;

