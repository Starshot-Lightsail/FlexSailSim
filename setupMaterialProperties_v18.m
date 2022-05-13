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

MATERIALS(1).name = 'SiNx37nm';
MATERIALS(1).Iabs = 1e-6; % absorptivity
MATERIALS(1).Irefl = 0.091; % reflectivity
MATERIALS(1).Emissivity = 0.1; % emissivity
    %lambda = 1550e-6; %mm, wavelength
MATERIALS(1).ior_n = 2; % refractive index
MATERIALS(1).ior_k = 4e-6;
MATERIALS(1).density = 2.7e-3; % g/mm3
MATERIALS(1).Ymod = 270e9;  % Pa, N/m^2
MATERIALS(1).Emod = 6.4e9;  % Pa, tensile strength
MATERIALS(1).nu = 0.27; % Poisson's ratio
MATERIALS(1).Thermcondmm = 0.010; % W/mm/degC
MATERIALS(1).heatCap = 0.9; % J/g/degC
MATERIALS(1).CTE = 2.0; % PPM/K, coefficient of thermal expansion
MATERIALS(1).failtemp = 1300; % Kelvin, failure temperature
MATERIALS(1).thickness = 37e-6; % mm

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
    
    
    MATERIALS(3).name = 'Si111_43nm';
MATERIALS(3).Iabs = 1e-6*20000/43; % absorptivity
MATERIALS(3).Irefl = 0.51; % reflectivity
MATERIALS(3).Emissivity = 0.3; % emissivity
    %lambda = 1550e-6; %mm, wavelength
MATERIALS(3).ior_n = 3.48; % refractive index
MATERIALS(3).ior_k = 4e-6;
MATERIALS(3).density = 2.33e-3; % g/mm3
MATERIALS(3).Ymod = 169e9;  % Pa, N/m^2
MATERIALS(3).Emod = 2.1e9;  % Pa, tensile strength
MATERIALS(3).nu = 0.262; % Poisson's ratio
MATERIALS(3).Thermcondmm = 0.160; % W/mm/degC
MATERIALS(3).heatCap = 0.67; % J/g/degC
MATERIALS(3).CTE = 2.33; % PPM/K, coefficient of thermal expansion
MATERIALS(3).failtemp = 1600; % Kelvin, failure temperature
MATERIALS(3).thickness = 43e-6; % mm

    
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
Ymodmm = Ymod./1e6;

