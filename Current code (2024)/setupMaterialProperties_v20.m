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

%% universal values
lambda = 1550e-6; %mm, wavelength (of laser propulsion)

%% Default material:
% these values are used throughout the sail unless modified by material assignments in the mesh generator
% Note: You can't use advanced optical models (LUTs, TMM, ...) with the default material!

material = 'silicon'; % just a display name, to prevent confusion

%Optical:
Iabs = 3e-10; % absorptivity
Irefl = 0.5; % reflectivity
    %note:  default material doesn't have refractive index because it is incompatible with the TMM code
    %To use specular TMM, assign a MATERIAL to the sail.

%Thermal:
Emissivity = 0.8; % emissivity
Emissivity_fwd = 0.8;    Thermcondmm = 0.157; % W/mm/degC
heatCap = 0.7; % J/g/degC
CTE = 2.33; % PPM/K, coefficient of thermal expansion
failtemp = 600; %kelvin

%Mechanical:
density = 2.3e-3; % g/mm3
Ymod = 169e9;  % Pa, N/m^2
Emod = 1.20e9;  % Pa, tensile strength
nu = 0.25; % Poisson's ratio
thickness = 43e-6; % mm 
    
disp(['Setting up physical constants and material properties.' ]);

MATERIALS = [];

% new (2023) settings per material
%   .Emissivity_fwd -- effective hemispherical emissivity of the forward-facing (i.e., away from the laser) surface of the lightsail, in case that's
%       different than the emissivity on the laser-facing side
%   .failtemp -- temperature (K) above which the material is assumed to fail

MATERIALS(1).name = 'SiNx37nm';
MATERIALS(1).Iabs = 1e-6; % absorptivity
MATERIALS(1).Irefl = 0.091; % reflectivity
MATERIALS(1).Emissivity = 0.1; % emissivity
MATERIALS(1).Emissivity_fwd = 0.1; % forward emissivity (facing away from laser, in case it's different?)
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
    MATERIALS(1).Emissivity_fwd = 0.1; % forward emissivity (facing away from laser, in case it's different?)
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
    
    
    MATERIALS(3).name = 'Si111_43nm_old';
    MATERIALS(3).Iabs = 7e-7; % absorptivity
    MATERIALS(3).Irefl = 0.51; % reflectivity
    MATERIALS(3).Emissivity = 0.8; % emissivity
    MATERIALS(3).Emissivity_fwd = 0.8; % forward emissivity (in case it's different?)
        %lambda = 1550e-6; %mm, wavelength
    MATERIALS(3).ior_n = 3.48; % refractive index
    MATERIALS(3).ior_k = 1.2e-6;
    MATERIALS(3).density = 2.33e-3; % g/mm3
    MATERIALS(3).Ymod = 169e9;  % Pa, N/m^2
    MATERIALS(3).Emod = 2.1e9;  % Pa, tensile strength
    MATERIALS(3).nu = 0.262; % Poisson's ratio
    MATERIALS(3).Thermcondmm = 0.160; % W/mm/degC
    MATERIALS(3).heatCap = 0.67; % J/g/degC
    MATERIALS(3).CTE = 2.33; % PPM/K, coefficient of thermal expansion
    MATERIALS(3).failtemp = 1000; % Kelvin, failure temperature
    MATERIALS(3).thickness = 43e-6; % mm

    MATERIALS(4).name = 'Si_43nm_450k';
    MATERIALS(4).Iabs = 20*7.70242e-09; % absorptivity
    MATERIALS(4).Irefl = 0.454941; % reflectivity
    MATERIALS(4).Emissivity = 0.1; % emissivity  %1.01271e-05 is correct for 450k
    MATERIALS(4).Emissivity_fwd = 0.1; % forward emissivity (in case it's different?)
        %lambda = 1550e-6; %mm, wavelength
    MATERIALS(4).ior_n = 3.48636; % refractive index
    MATERIALS(4).ior_k = 1.29905e-07;
    MATERIALS(4).density = 2.33e-3; % g/mm3
    MATERIALS(4).Ymod = 169e9;  % Pa, N/m^2
    MATERIALS(4).Emod = 2.1e9;  % Pa, tensile strength
    MATERIALS(4).nu = 0.262; % Poisson's ratio
    MATERIALS(4).Thermcondmm = 0.160; % W/mm/degC
    MATERIALS(4).heatCap = 0.67; % J/g/degC
    MATERIALS(4).CTE = 2.33; % PPM/K, coefficient of thermal expansion
    MATERIALS(4).failtemp = 1200; % Kelvin, failure temperature
    MATERIALS(4).thickness = 43e-6; % mm
    %MATERIALS(4).optDBtemps = [300 450 600];
    %MATERIALS(4).optDBn = { 'emis/si_n.csv' 'emis/si_n.csv' 'emis/si_n.csv' };
    %MATERIALS(4).optDBn = { 'emis/si_absco_300k.csv' 'emis/si_absco_450k.csv' 'emis/si_absco_600k.csv' };
    %MATERIALS(4).optDB_emisoverride = { 0.5 0.5 0.5 };

    


