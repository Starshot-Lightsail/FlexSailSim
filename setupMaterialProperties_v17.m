%% setupMaterialProperties_v17
% (c) 2018-2021 Michael Kelzenberg, Ramon Gao -- California Institute of Technology
% Sets up material properties variables... oh, and also physical constants used in simulations
% 
% These properties are set by this subscript in hopes of ensuring better consistency between various versions of the 
% simulator code.

%% Physical constants
% Residents of dimension C137 shouldn't need to change these values:
c0 = 299792458; % m/sec, speed of light.
c0mm = c0*1000; % mm/sec, speed of light
SBC = 5.67e-8; % W/m^2/K^4; boltzmann constant, for radiative cooling.
SBCmm = SBC / 1e6; % W/mm^2/K^4

%% Material selector
% String variable helps us select between different sets of properties and behaviors throughout the simulations.
material = 'silicon' ;

disp(['Setting up physical constants and material properties.  Material = ''' material '''' ]);

if strcmpi(material,'silicon')
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
    thickness = 43e-6; % mm

elseif strcmpi(material,'nitride')
    %I0 = 1000;%1000; % W/mm^2, peak incident power density (1000 = 1GW/m2)  %MK: it goes faster at 10000
    Iabs = 1e-5; % absorptivity
    Irefl = 0.35; % reflectivity
    Emissivity = 0.1; % emissivity
    lambda = 514e-6; % mm, wavelength
    ior_n = 2.2;  % refractive index
    ior_k = 0.001;
    density = 3.1e-3; % g/mm3
    Ymod = 270e9;  % Pa, N/m^2
    %warning('Using fake value of modulus!');
    %Ymod = Ymod/50;
    Emod = 5.5e9;  % Pa, tensile strength
    %warning('Using fake value of tensile strength!');
    %Emod = Emod/100;
    nu = 0.25; % Poisson's ratio
    
    Thermcondmm = 0.02; % W/mm/degC
    heatCap = 0.75; % J/g/degC
    CTE = 2.4; % PPM/K, coefficient of thermal expansion
    thickness = 176e-6; % mm
end

tensile_strain_limit = Emod/Ymod;
thicknessnm = thickness.*1e6;
Ymodmm = Ymod./1e6;

