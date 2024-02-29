%% STARSHOT LIGHTSAIL SIMULATOR
% (c) 2018-2023 Michael Kelzenberg, Ramon Gao -- California Institute of Technology
% Largely well commented, with external documentation now in development.  Ad Astra!

% This is a modified version 19, which looks up simulated optical pressures
% using Ramon's slow approach. The reason for this is due to discrepancies
% between the interpolated pressures from loadLUTs_v17.m and the actual
% tabulated pressures, especially at 0 deg. Rigidified simulations with
% Ramon's slow LUT approach agree much better with numerical solutions from
% solving the equations of motion (ODEs).

% 05/19/2023: Implemented advanced numerical integration/time stepping such
% as Midpoint or RK4 by outsourcing sections to separate auxiliary script

disp('Flexible lightsail simulator ');
disp('(c) 2018-2023 Michael D. Kelzenberg and Ramon Gao, California Institute of Technology');
disp('Version: 2023-05-04 (v20)');

%% Initialization block:
%clear  % clear the workspace to make sure the code will run in a fresh workspace, and prevent carry-over of prior
        % simulation data.  We often comment this out during development though, to preserve analysis/support variables.

close all hidden  % close all windows from prior simulations or whatever.  The 'hidden' part closes message boxes too.

% Limit maximum number of cores to ~ of cores on ahura
maxNumCompThreads(30)

% Set this to 0 if using Ramon's LUT approach, or 1 for Mike's LUT approach
read_LUT = 0;

simdesc = 'flexM1e80HzFig6RK4'; % Use this if you want to add a short description of why you're running the simulation.  
              % It will be appended the output file names.  You can leave it blank ('').  Or, write something 
              % meaningful if you are performing a specific study.
              % I recommend limiting this to 10 characters, safe for file names (e.g., no '/').  

try %close prior video recorder.  I will be left open if the prior simulation errored out.
    if exist('V1')
        close(V1);
    end
catch ME
    warning('Closing the prior simulation''s video recorder produced some sort of error.');
    %whatever... we tried
end

startscripttic = tic;  %so we can record the elapsed time for the entire script

% Each simulation has a unique identifier number, which should automatically change between simulations.
% This value is used in the creation of all the output files resulting from each simulation.  That way, each time a
% simulation is run, a unique set of output files are produced, and it is easy to figure out which output files
% correspond to which simulation.
% The identifier is called 'movieNum' and each time the simulation script is run, we create a new movieNum based on the
% date and time, down to the second so hopefully we don't ever devise a duplicate movieNum.  If the simulation 
% variable iteration loops are used ("outer loops" 1 and 2 below) the movienum increments by 1 between each iteration.
% That way we'll know which simulations are grouped together by their common filename bases.  
ungt = datetime;
ungts = sprintf('%02d%02d%02d%02d%02d%02d00', ungt.Year-100*round(ungt.Year/100), ungt.Month, ungt.Day, ungt.Hour, ungt.Minute, round(ungt.Second));
movieNum = str2num(ungts);  %for example, 21042202550700.  Note:  if you're going to iterate through more than 100 conditions, add an extra 0 to the initial movieNum.
if ~isempty(simdesc) 
    filebasename = sprintf('sim%d_%s', movieNum, simdesc);
else
    filebasename = sprintf('sim%d', movieNum );
end

% Make general folder for figures, if it does not already exist
if not(isfolder('Figures'))
    mkdir('Figures')
end

% Make general folder for workspaces, if it does not already exist
if not(isfolder('Workspaces'))
    mkdir('Workspaces')
end

% Make general folder for results, if it does not already exist
if not(isfolder('Results'))
    mkdir('Results')
end

% Save path to parent directory of simulation
parent_dir = pwd;

% Create subfolders for figures for better organization & remember the path
cd('Figures')
if not(isfolder(filebasename))
    mkdir(filebasename)
    cd(filebasename)
    path_figures = pwd;
end

% Create subfolders for scripts for better organization & remember the path
cd(parent_dir)
cd('Scripts')
if not(isfolder(filebasename))
    mkdir(filebasename)
    cd(filebasename)
    path_scripts = pwd;
end

% Create subfolders for results for better organization & remember the path
cd(parent_dir)
cd('Workspaces')
if not(isfolder(filebasename))
    mkdir(filebasename)
    cd(filebasename)
    path_workspaces = pwd;
end

% Create subfolders for results for better organization & remember the path
cd(parent_dir)
cd('Results')
if not(isfolder(filebasename))
    mkdir(filebasename)
    cd(filebasename)
    path_results = pwd;
end

cd(parent_dir)

% Here we save a copy of the exact simulation code used for each simulation session, so we can go back and deduce
% the specific settings and algorithms that were used for any recorded simulation
copyfile([ mfilename '.m'], fullfile(path_scripts,[filebasename '.script']));

% We also save the terminal output so we can see what was displayed during the simulation.
diary(fullfile(path_results,[filebasename '.diary'])); % Log command window text to file


%% Choose numerical integration (time stepping) method
% Currently, available methods are semi-implicit Euler ('SIEuler'), velocity Verlet ('VVerlet') and Runge-Kutta 4 ('RK4')[
listIntegrationMethods = {'SIEuler', 'MP','VVerlet', 'RK4'};
namesIntegrationMethods = {'Semi-implicit Euler', 'Midpoint','velocity Verlet', 'Runge-Kutta 4'};

% Set index for method: '1' -> 'SIEuler', '2' -> 'VVerlet, '3' -> 'RK4'
integrationMethodIdx = 4;
integrationMethod = listIntegrationMethods{integrationMethodIdx};
nameIntegrationMethod = namesIntegrationMethods{integrationMethodIdx};

disp(['Choosing the ', nameIntegrationMethod, ' method for numerical integration.']);


%% Set up material properties  (external script)
% Defines basic material properties and physical constants.  As of this writing, supported materials are 'silicon' and 
% 'nitride.'  For various reasons, we define certain lightsail specific parameters here as well, such as thickness.  If
% you want to change these values iteratively during a batch of simulations, make sure that such changes are applied
% prior to mesh generation.
setupMaterialProperties_v19_final2_RG;

%% Import optical response of MEPS textures as lookup tables (LUTs)  (external script)
% MK (April 6-8 2021):  This is a new approach to loading and accessing the MEPS data from Ramon's simulations.  This
% approach can also be used to load any other optical properties of the lightsail regions, such as from theory or from
% measurements.  See the load function called below for uncharasterically thorough comments!
%
% Since we don't generally change the source of the optical data between iterative simulations, we load it here at the
% beginning of the script to save time and prevent redundant outputs.
%
% IMPORTANT:  As of this writing, All LUTs must use the same values for LUTPsiStepDeg and LUTThetaStepDeg!
%  The values retained in LUT(n).pstep and .tstep are not used elsewhere in the simulation.
%  Using the same LUT sizes reduces the complexity of indexing the tables during the time-domain simulations.  

LUTPsiStepDeg = 0.20;    %  These specify the angular resolution to use for the LUTs.  The values might be overwritten 
LUTThetaStepDeg = 0.20;  %  in setupLUTs depending on the code configuration.  But you need to set a positive nonzero 
                         %  value for these two settings, even if not using LUTs, otherwise the simulation will error
                         %  out.

% If not using LUTs, set LUTs = [] here instead of calling setupLUTs().                         

if read_LUT
    setupLUTs_v17;  % now moved to subscript.  
end

% Note:  If you want to change LUTs between simulations, you can modify the setupLUTs script to change its behavior
% based on 'outerloop' variables, then re-invoke it within each loop, e.g., near generateMesh.  Or I often find it
% easier to load all the LUTs for all designs here at the beginning, then assign whichever LUTs are actually used for 
% each structure during generateMesh.                                     

%% Import optical response of specular (non-MEPS) surfaces as lookup lists (LULs)
% Unlike complex nanophotonic grating structures (MEPS), the optical resonse of simple thin-films membranes, and even 
% 1D stacks such as bragg reflectors, can be calculated as a function of incidence angle via closed-form expressions.
% Despite this, in practice, it's almost always faster to use a look-up approach in this case as well, in which we 
% obtain approximate values for reflectance and absorption by retrieving pre-calculated values from memory rather than 
% calculating them on the spot.  For specular surfaces, this can be accomplished by a 1D look-up-table.  For the sake of
% clarity, we call these "look-up lists" (LULs) to differentiate them from the 2D LUTs used for MEPS.
%
% Currently, we don't support having multiple "regions" of different specular optical response, but I will probably add
% this support in the future, in the same way that we support multiple LUT (texture) region mapping for MEPS. For now, 
% we assume that all specular lightsail surfaces are made from the same material and thickness(es), thus a global 
% specular response value/LUL is used.  This specular behavior is used for all triangles with "t_tex" values of zero. 
%
%There are two modes supported for the specular response calculation:
%
%   useSpecularTMM = 0:  Uses constant scalar values for the the reflectance and absorption of the lightsail surface,
%                        regardless of the incidence angle.  The values are Irefl and Iabs, and an example value for 
%                        each is 0.98 and 0.02, which would reasoanbly describe a metallic silver lightsail surface 
%                        (e.g., a solar sail).  This is a reasonable assumption for such a sail, because metalic silver
%                        has similary high reflectance across all incidence angles.  This simple approximation might 
%                        also be useful for conceptual studies or thought experiments.  However, as we know, metallic 
%                        reflectors are unsuitable for "starshot" illumination levels, because they would be destroyed.
%                        Dielectric membranes can achieve much lower absorption, but the reflectance and absorption 
%                        of such membranes can vary considerably with incidence angle (and also polarization, which we
%                        don't yet support sorry.)  Thus for dielectric membranes, set UseSpecularTMM to 1...
%
%                        Specular behavior determined by:  Irefl (reflection coefficient) and Iabs (absorption)
%
%   useSpecularTMM = 1:  Uses a pre-calculated list of values to determine the specular reflectance (and absorption) as 
%                        a function of incidence angle.  When I wrote this functionality, I was assuming that such 
%                        values would be calculated via the transfer matrix method (TMM), hence the name.  But you can 
%                        generate the list of of R and A values however you want.  The supplied calculator works for the
%                        simple case of a single dielectric membrane with complex refractive index ior = ior_n - i*ior_k
%                        (Yes, the complex ior has negative complex value, because that's the convention used in the
%                        reference I followed when writing the code... deal with it.  Do not specify negative values for
%                        ior_k, just trust the code and be sure to drink your ovaltine.  
%
%                        Specular behavior determined by:  Values stored in LUL structure (A,R) spaced by LULAngleStep
useSpecularTMM = 0;
LULAngleStepDeg = 0.2; % angular resolution of LUL.  Don't comment this out, even if not using SpecularTMM!
LULAngleStep = deg2rad(LULAngleStepDeg);
% setupSpecularLULs_v18; % Ramon: Didn't get this file from you, Mike

%% Import and set up simulated optical pressures from COMSOL simulation

% Load all small look-up tables of pressures vs roll and pitch angles for
% TE- and TM-polarized light. All the small look-up tables for TE (TM)
% polarization combined together yield the full look-up table for TE (TM)
% polarization. The full look-up table is being splitted into smaller ones
% to increase speed of finding pressures within the table later in the code

if ~read_LUT
    
file_name_collection_small_tables = ...
    '08232021_MEPS_Si3N4_Mark1e_1064_PitchRoll_Pressures_CollectionTables_';

load([fullfile('Data',file_name_collection_small_tables) 'TE.mat'])
load([fullfile('Data',file_name_collection_small_tables) 'TM.mat'])

% global last_indices_tables
% global num_small_tables

dim_collection_small_tables = size(pressure_TE_small_tables);
num_small_tables = dim_collection_small_tables(3);

collection_phi_angles_tmp = zeros(num_small_tables, dim_collection_small_tables(1));

last_indices_tables = zeros(1,num_small_tables);

for ii = 1:1:num_small_tables
%     tmp = [];
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

first_angles = collection_phi_angles(:,1);
step_angles = abs(collection_phi_angles(:,1) - collection_phi_angles(:,2));

last_angles = zeros(size(first_angles));
for ii = 1:1:length(last_indices_angles)
    last_angles(ii) = collection_phi_angles(ii,last_indices_angles(ii));
end

% global border_angles

border_angles = (last_angles(1:end-1) + first_angles(2:end))/2;

end

%% Definition of direction cosine matrix (rotation matrix)

% Define function that describes the equations of the system in vector form

DRM123_11 = @(psi,theta,phi) cos(psi).*cos(theta);
DRM123_12 = @(psi,theta,phi) cos(psi).*sin(theta).*sin(phi) + sin(psi).*cos(phi);
DRM123_13 = @(psi,theta,phi) -cos(psi).*sin(theta).*cos(phi) + sin(psi).*sin(phi);
DRM123_21 = @(psi,theta,phi) -sin(psi).*cos(theta);
DRM123_22 = @(psi,theta,phi) -sin(psi).*sin(theta).*sin(phi) + cos(psi).*cos(phi);
DRM123_23 = @(psi,theta,phi) sin(psi).*sin(theta).*cos(phi) + cos(psi).*sin(phi);
DRM123_31 = @(psi,theta,phi) sin(theta);
DRM123_32 = @(psi,theta,phi) -cos(theta).*sin(phi);
DRM123_33 = @(psi,theta,phi) cos(theta).*cos(phi);

%% WHAT'S THE DIFFERENCE BETWEEN "SPECULAR/LUL" AND "MEPS/LUT" BEHAVIOR?  
%  In short, specular behavior is calculated by 1D lookups (or 0D for fixed-reflectance surfaces), whereas MEPS behavior
%  is calculated by 2D lookups, as a function of incidence angle.  ONLY specular/LUL lightsail surfaces can be used with
%  RayTracingMode=1 multi-reflection optics, whereas both specular/LUL and MEPS/LUT surfaces can only be used with 
%  RaytracingMode=0 single-surface optics.  
%
%  - MEPS/LUT behavior is used to simulate complex nanophotonic structures with non-specular behavor, such as gratings, 
%      which can produce lateral optical forces and torques in almost any direction depending on the design and incidence
%      angle. To model the behavior of these regions, we keep track of the orientation of the MEPS "texture" vector 
%      within each triangle, and use a 2D lookup table to determine the optical response based on how the incident light 
%      interacts with nanophotonic surface at each specific incidence angle and orientation of the surface texture.  
%      The lookup tables produce net optical forces and absorption for any incidence angle, but give us no insight 
%      into the direction(s) that any reflected, transmitted, or scattered light will go.  Thus, we can't determine if 
%      reflected light would impact any other region of the lightsail.  Thus, MEPS/LUT behavior is primarly used to
%      model flat or near-flat lightsails, in which stability is obtained via the nanophotonic texturing and layout of 
%      the sail.  
%
%  - SPECULAR/LUL behavior is used to simulate specular optical surfaces such as dielectric films or 1D stacks thereof.
%      Specular surfaces transmit and reflect light in a specific, known direction, and do not scatter or diffract light
%      in other directions.  Furthermore their reflectance and transmittance coefficients do not depend on the specific
%      in-plane orienation of the material surface.  Thus, we only need to know the scalar values of reflectance and
%      absorption to calcualte the optical forces for such surfaces, and these values depend only on a single value 
%      (incidence angle).  (Note:  the situation becomes more complex for polarized light, but our simulator doesn't yet
%      handle polarization for ANY optical response model.)  Because specular surfaces reflect light in a single known
%      direction, we can determine if the reflected light will impact another part of the lightsail, and thus reasonably
%      aporoximate the behavior of curved lightsail shapes in which the incident light is reflected off the surface more
%      than one time (see below).

%% Single vs multi-surface optics mode  (RayTracingMode)
%  A perflectly flat lightsail, or a slightly distorted flat lightsail, can be adequately modeled by a single set of  
%  surface interactions with the incident beam.  The incident light arrives, and is either absorbed, or is reflected, 
%  transmitted, or scattered away from the lightsail, producing a certain amount of photon pressure and absorption 
%  heating on that portion of the sail.  But we don't have to worry about any of the reflected/transmitted/scattered 
%  light impacting any other portion of the lightsail.  
%
%  But for substantially convex lightsails, e.g., curved shapes such as parabolods or cones, we might expect some of the
%  reflected, transmitted, or scattered light to impact other portions of the sail, resulting in additional optical
%  forces that we must consider before we can expect our simulations to stand up to the scrutiny of peer-review.  
%
%  We have no hope or intent of modeling the behavior of every possible combination of surface behavior and lightsail 
%  shape, but one particular class of lightsails is of partular relevance to our current work:  We wish to study open, 
%  convex, upright lightsails, such as paraboloids or cones, made from dielectric membranes having specular surface 
%  behavior.  We consider this such shapes of specular lightsails to be one of the two primary competing design 
%  approaches to address the challenges of the Breakthrough Starshot Lightsail program.  (The other being flat 
%  lightsails with nanophotonic texturing.)  
%
%  So, to support modeling these convex shaped, specular lightsails, we have implemented a crude yet effective first-
%  order multi-reflection optics calculation method, which considers a finite number of specular reflections within the 
%  interrior 'bottom' surface of the lightsail only.  
%
%  Limitations and requirements for RayTracingMode:
%  - The lightsail can only contain specular surfaces.  MEPS/LUT REGIONS CANNOT BE USED.
%  - All reflections occur only between discrete triangles, without interpolation or dithering; thus the accuracy of 
%      the approximation degrades considerably for coarse meshes with few triangles.  
%  - Ray-tracing begins with the calculation of the initial optical interaction between the incident beam and each
%      triangle, evaluated at the centroid of each triangle, just as in the non-raytracing scenario.  Following this,
%      the reflected beam is calculated for each illuminated triangle, and iterative raytracing proceeds. 
%  - Each reflected beam originates at the centroid of the source triangle, and impacts either zero or one other triangles.  
%      - if the reflected beam intersects zero triangles, then this branch of raytracing terminates.
%      - if the reflected beam intersects one other triangle, then the optical forces and absorption are calculated for 
%          the impacted triangle as if the entire beam impacted the centroid of that triangle (distributed equally to 
%          each of the three impacted nodes, without any second-order corrections based on the proximity of impact to 
%          any node or adjacent triangle.  Then, the power and direction of the reflected beam is calculated, and 
%          raytracing continues, with the new reflected beam originating from the centroid of the impacted triangle, 
%          again, regardless of actual point of impact.
%      - if the reflected beam intersects more than one other triangle, one and only one triangle will be selected for
%          calculation of optical interactions as described above.  I do not know which triangle will be selected, and 
%          the process may well not be deterministic.  In general, multiple triangle intersections should occur only
%          in the case of concave lightsail shapes, which fall outside the scope of interest for our current work.
%   - Arithmetic errors and numerical tolerances may result in a small percentage (hopefully <<1%) of raytracing 
%          calculations terminating prematurely or incorrectly.   I haven't studied this in detail, but believe the
%          effect to be entirely negligable in typical applications.  
%   - Raytracing does not occur for transmitted light, as no further interactions are expected for convex upright 
%          lightsails.  
RayTracingMode = 0;   % whether or not to raytracing stuff
MAX_REFL = 2;     % maximum number of ray reflections

if RayTracingMode
   % convert LULs to gpuArrays for RayTracing mode. 
end

init_time = toc(startscripttic);

%%  Outer sweeps setup
% The simulator can perform a 2D parameter sweep via two "outer loops".  Between each iteration of the outer
% loops, various simulation settings can be changed.  An example application would be to study the stability of
% acceleration as a function of the initial beam offset (t0t_x0) and beam diameter (Irad).
%
% Specify values for the inner and outer sweep variables here.  For a single simulation, enter a single (scalar) value
% for each outerseepvals array.  If you want to simulate a bunch of various configurations affecting multiple variables, 
% changes, you can implement that as a list with a state variable as the index, so long as you modify the initialization
% code appropriately.  
%
% Implementing the outer sweeps will require some modification of the simulation code below.  For example, if you are
%  sweeping the values of t0t_x0 and Irad, find and change the definition of each to :
%   t0t_x0 = outersweepvals1(nouterloop1)
%   Irad = outersweepvals2(nouterloop2)
%
% Note:  outer loop #1 is the outermost of the two outer loops.
%
% Note from Ramon:  Mike's definition of the beam width is two times that of Oggy's. More importantly, this value is
% multiplied by the radius, so another factor of 2 is needed to account for beam width in terms of the diameter

beam_width = 0.4; % in units of sail's diameter [D];
spin_frequency = 80; % in Hz

outersweepvals2 = [ [beam_width, spin_frequency] ]; % sweeping beam width & spin_frequency simultaneously
outersweepvals1 = [ 0 ]; %currently used for:  raytrace_mode
burstTestVals = [];

outerloopstable = ''; %this initializes a string array to store the output table for all subsequent simulations
multiloop = (length(outersweepvals1) + length(outersweepvals2)) > 2;
nloop = 0;

for nouterloop1 = 1:length(outersweepvals1)
    for nouterloop2 = 1:size(outersweepvals2,1)
        startlooptic = tic;
        nloop = nloop + 1;
        outersweepval1 = outersweepvals1(nouterloop1);
        outersweepval2 = outersweepvals2(nouterloop2,:);
        
        fprintf('\n');
        if (length(outersweepvals2) > 1) || (length(outersweepvals1) > 1)
            disp(['OuterLoop [' num2str(nouterloop2) ',' num2str(nouterloop1) '] '  num2str(outersweepvals2(nouterloop2)) ' ' num2str(outersweepvals1(nouterloop1)) ] );
        end
        
        RayTracingMode = outersweepvals1;

        %% Generate mesh
        cd(fullfile(parent_dir,'Scripts','Auxiliary'))
        generateMesh_v20_medium_RG; % This helper script sets up the mesh.  Modify it to change mesh settings such as shape, patterning, and material properties.
        cd(parent_dir)
        
        if ~read_LUT
           % Get indices of which triangles correspond to which texture region
            % indices_Region0 = find(t_tex == 0);
            indices_Region1 = find(t_tex == 1);
            indices_Region2 = find(t_tex == 2);
            indices_Region3 = find(t_tex == 3);
            indices_Region4 = find(t_tex == 4);

            % Pre-allocate memory for pressures & forces on mesh triangles
            t_press_x = zeros(size(t_a));
            t_press_x_tmp = zeros(size(t_press_x));
            t_press_y = zeros(size(t_a));
            t_press_y_tmp = zeros(size(t_press_y));
            t_press_z = zeros(size(t_a)); 
            t_press_z_tmp = zeros(size(t_press_z));
        end
        
        saveas(gcf,fullfile(path_figures,[filebasename '_mesh.fig']));  % save a rendering of the initial mesh, which is created as a
        saveas(gcf,fullfile(path_figures,[filebasename '_mesh.png']));  % figure by the generateMesh script
                
        %% Illumination mode and power
        %  As of this writing, only 'gaussian' illumination profiles are supported.  We plan to add arbitrary profiles via LUTs
        %  in the future.  Beam is always centered at (0,0).  Starting sail offsets have now been moved to a different
        %  section.
        illum_mode = 'gaussian';  %Illumination mode selector.  Not yet implemented in the code (FIXME).
        I0 = 1000; % W/mm^2, peak incident power density (1000 = 1GW/m2).
        
        %% Illumination settings for gaussian beam profile.
        Irad = (outersweepval2(1)) * sqrt(2) * radiusmm; % mm, beam waist radius, for gaussian illumination
        
        %% Illumination settings for other beam profile modes.  To be added in the future.
        if ~isequal(illum_mode,'gaussian')
            error(['Attempted to use non-defined illumination mode ''' illum_mode '''']);
        end
                
        %% Spin-stabilization settings
        rps_target = outersweepval2(2);%outersweepval2; % Rotations per second for spin stabilization.  Can be negative to spin the other way.
        spinup_gs = max(rps_target/100,10);  % How quickly to spin-up the lightsail.  As a time-domain simulator, the lightsail must begin at rest.
        %  We spin it up to a desired rotation speed during a brief time prior to the propulsion phase.  I'm not sure what
        %  the units are for this parameter, but it works out that we get 1808.8 rotations per second per simulated second of
        %  spin-up for a value of spinup_gs = 1.  Depending on material properties, too much angular acceleration can 
        %  result in unintended pressure waves ("shock") developing in the sail during spin-up.  Too rapid a spin-up 
        %  will also reduce the number of time steps over which spin-up occurs, leading to slight slight errors in the
        %  achieved rotational speed vs. rps_target, due to the discrete stepping.  On the other hand, spinning up too 
        %  slowly wastes time.  simulation time.  
        % NOTE:  Don't make this negative.
        
        % this is the resulting duration for the spinup phase.  It is calculated, do not change:
        spinup_time = abs(rps_target)/1801.8/abs(spinup_gs); % How long to apply the spinup forces, to achieve the desired RPS        
        
        %% Simulation start-up and duration
        % The lightsail launch sequence begins t=0.  If we want a spin-stabilized sail, it needs to be already rotating 
        % by this point in time.  Due to the specific implemenation of our time domain simulator, we can't instantaniously 
        % accelerate the sail from rest to the desired speed at t=0.  Instead, we must allow this spin-up process to 
        % occur over a certain duration of time, starting at t = -spinup_time, and ending at t=0, at which point the 
        % lightsail should be spinning about its axis at the desired rate.
        %
        % Note:  at t=0, the sail can be moved or rotated to simulate a position or tilt error prior to acceleration.
        % See below.
        %
        % Then, optionally, the start of illumination can be delayed by a further amount t_ramp_delay, to allow
        % us to observe the behavior of the spinning sail prior to acceleration.
        %
        % At t=t_ramp_delay, we begin turning on the acceleration beam.  In our time-domain simulation approach, it is useful to
        % turn on the beam gradually, so we don't induce shockwaves in the lightsail due to abrupt application of high forces.
        % This is accomplished over a period of t_ramp_dur, during which the illumination profile is varied from 0 to 100%
        % intensity following the profile of a half-phase sinusiod function.  Specifically the intensity weighting function is:
        % intensity weight factor  = 1 - 0.5*cos( pi * t' / t_ramp_dur), where t' is the elapsed time since starting the beam
        % turn-on.  The beam thus reaches full power at t = t_ramp_delay + t_ramp_dur seconds.
        %
        % After reaching full power, the beam remains on continuiously for t_on_dur seconds.  During this phase we can observe
        % the sail's acceleration and stabililty, or just say we did while browsing social media instead.
        %
        % Then, the beam is switched off via the reverse of the turn-on ramp (over the course of t_ramp_dur seconds), and 
        % the simulation proceeds for an additional t_coast_dur seconds, so that we might observe the final trajectory of the sail.
        %
        % Thus the simulation terminates at t_ramp_delay + 2*t_ramp_dur + t_on_dur + t_coast_dur.
        %
        % Note that simulations can be automatically terminated eariler if the sail veers off course or breaks apart (see below.)
        t_ramp_delay = 0*1e-3; % (seconds), time delay before starting beam turn-on, formerly t_start_ramp
        t_ramp_dur = 1e-6; %20;%3e-3; % (seconds), time duration of beam turn-on. DO NOT SET TO ZERO!  (Use a small number e.g., 1e-8 for instant turn-on.)        
        t_on_dur = 1.01; %00e-3; % (seconds) duration of acceleration phase, after which beam is turned off
        t_coast_dur = 00e-3; % (seconds) duration of time after beam turn-off for observing the final sail trajectory.
        
        dt = rec_dt; % seconds, the time step for the simulation.  rec_dt is the "recommended" calculated value from the mesh generator.  I do not recommend changing this.
        if rps_target  % reduce time step further if needed for rapidly spinning structures.  
            if dt > (1/rps_target/360)
                dt = (1/rps_target/360);
                warning(['Reducing time step from mesh-calculated value due to high spin speed.  dt = ' num2str(dt)]);
            end
        end

        % This is the calculated value for the total simulation duration after t=0; do not change:
        sim_dur = t_ramp_delay + 2*t_ramp_dur + t_on_dur + t_coast_dur; % seconds; note this doesn't include the t<0 spinup time.
        
        % Here is the actual power ramp function.  Do not change.
        pwrRampFcn = @(tt) 0.5*double(tt > t_ramp_delay).*double(tt <= t_ramp_delay+t_ramp_dur).*(1-cos(pi*(tt-t_ramp_delay)./(t_ramp_dur))) + double(tt > t_ramp_delay+t_ramp_dur).*double(tt <= t_ramp_delay+t_ramp_dur+t_on_dur) + ...
            0.5*double(tt > t_ramp_delay+t_ramp_dur+t_on_dur).*double(tt <= t_ramp_delay+t_ramp_dur+t_on_dur+t_ramp_dur).*(1+cos(pi*(tt-t_ramp_delay-t_ramp_dur-t_on_dur)./t_ramp_dur));
        
        %% Initial sail orientation and offsets  (t=0 transform)
        % The initial mesh and the beam are both always centered at (0,0).  So, to simulate a beam offset, tilt offset,
        % or a specific starting angular orientation for the sail, we need to translate or rotate the mesh after
        % spin-up, prior to beginning the acceleration.  (Spin-up forces are centered around (0,0) as well.)
        %
        % This is done as an instantanious rigid transformation at t=0, in violation of physical continuity, but 
        % ultimately this appraoch made the most sense.  
        t0t_changeSpinAngle = 1; % whether or not to change the rotational position of the sail, useful for spinning sails
        t0t_s0 = psi0;  % the desired rotational position, if enabled.  
                                      %Typically we use psi0 from the mesh generator -- the initial angle of the mesh.
                                      %For a different rotational position, use desiredValue-psi0
        
        t0t_x0 = 0.1*radiusmm;
        t0t_y0 = 0.1*radiusmm;
        t0t_t0_x = deg2rad(-2); % (rad), initial tilt of sail about the x axis
        t0t_t0_y = deg2rad(-2); % (rad), initial tilt of sail about the y axis
        t0t_t0_z = 0*deg2rad(45); %  (rad), initial tilt of sail about the z axis

        %% Simulation video recording mode and rendering options:
        % In general it is nice to plot and record as much data as possible, but that can slow down the simulations.  Use
        % plotFast=1 to render a simplified plot of the lightsail for the video recorder, or plotFast=2 to disable plotting and
        % video recording alltogether.
        videoMode = 2;  %0 = full plots; 1 = simplified plot (slightly faster movie rendering); 2 = no plots (no video; fastest option)
        zreliefmag = 1; %increase z-values of mesh points in 3D plots by this amount to accentuate shape distortions.  Should be at least 1, do not set to zero.
                            
        %% Output data logging frequency, video frame rates, etc.
        % The time step used in our calcuations is automatically calculated, based on the material properites and scaling of
        % the simulation mesh, to ensure fidelity of the results within the capabilities of floating point arithmatic.
        % With this approach, typical simulations produce useful results over a simulated time period of several seconds of
        % acceleration.
        %
        % Ideally, we would record the entire state of the simulation mesh at every time step, so we could go back and study
        % what was happening anywhere at any point in time.  But in practice, this generates too much data for our
        % computers to
        % handle, and we rarely need such fine time resolution (exception:  FFT analysis, see below).  Also, we rarely need to
        % know the entire mesh state, instead we can monitor a smaller subset of data ('key state variables').  For example, 
        % instead of recording the location of every mesh node, we can record the location of the center of mass. 
        %
        % For saving key state variables, we downsample the output monitors by the following factor.  For example, if
        % sim_output_downsample is 4, we will record key state variables only every 4th time step.  This saves time and memory,
        % but if complete time resolution is desired, set sim_output_downsample to 1.  NOTE: during simulations, the
        % timestep iteration index is 'nt', whereas the downsampled output monitor index is 'nto'.  After the time
        % domain simulation ends, we set nt = nto, so that nt values correspond to the downsampled data sets.
        sim_output_downsample = 8;  % typical recommended value: 8
        sim_angles_downsample = 32;  % typical recommended value: 8
        
        % This parameter lets us skip re-calculation of optical forces during most time steps.  The optical forces on
        % each triangle tend to change fairly slowly compared to the mechanical forces.  We can make simulations run
        % faster, with little loss of accuracy, by reducing the frequency of optical force calculation.  However, this
        % does fundamentally reduce the accuracy of simulations.  Value of 2 to 16 seem to provide reasonalbe speedup
        % for our simulations, without substantially affecting the results.  
        opt_skip_downsample = 1;  % number of time steps between optical force calculation.  Do not set to zero!
                                % Set to 1 for maximal accuracy of optical forces.  
        
        % monitor_mode:  if monitor_mode is set (1), we will record a handful of state variables at each iteration of the
        % loop, with no downsampling, so we can perform FFT analysis later.  This doesn't affect the key state
        % variables, which are still stored at downsampled intervals.  
        fft_mode = 0;  % set to 1 to record node position and forces at all times, then perform FFT analysis
        
        % Our code also produces a video recording of the lightsail's behavior, rendered via a MATLAB figure window.  Rendering
        % each video frame is computationally intensive, typically requiring 100-10,000x longer than the calculation time for each
        % simulation time step.  Furthermore the visual changes to the mesh between consecutive time steps are impercievably
        % small.  Thus we can render video frames at a substantially reduced frequency, specified here by the value
        % frame_interval_s (in seconds).  Typical values of interest range from ~1e-5 seconds for recording extreme acceleration
        % events or catastrophic shape failures, to 1e-3 seconds to illustrate the trajectory and steady-state behavior of 
        % stabilized lightsail designs.
        frame_interval_s = 2*1e-4;     %seconds
        
        % frame rate reduction for spinning sails.  If you want the rotation of the sail to be clearly percievable in
        % the video recording, uncomment this line to make sure the frame rate is slow enough to capture it.  
        % if (rps_target)  frame_interval_s = min(frame_interval_s, 1/rps_target/36); end % make sure there are at least 36 frames per revolution
                
        % Skip rendering the spin-up sequence?  This prevents video rendering during spin-up, which can save time.
        skip_spinup_video = 0;  
        
        % When simulating unstable lightsail designs resulting in structural failure, it is sometimes of interest to visulize
        % their failure with greater time resolution, to help understand the cause.  This variable can be used to reduce the
        % video rendering time step upon the onset of structural failure:
        breakup_slowmo_interval_s = max( frame_interval_s/10, 5e-6);   %seconds; increases frame rate 10x for failures
        
        %% Snapshot display settings
        % Our simulator code also produces a composite 'snapshot' rendering of the evolution of the lightsail shape and
        % (transverse) position during the simulation, overlaying snapshots of the sail at fixed time intervals upon a single
        % axis, with each consecutive snapshot offset vertically above the prior one, accompanied by a timestamp annotation.
        % Because only a handful of shape renderings can be effectively displayed in a single figure, it is difficult to set
        % this value correctly in advance.  At a minimum, the sail's initial state and final state are rendered in this figure,
        % and with judicious choice of the snapshot interval, the evoloution from the former to the latter can be effectively
        % conveyed within a single static figure.  Good luck.
        snapshot_interval_s = 0.1;%0.020;%.250;
%         snapshot_interval_s = 0.25;%0.020;%.250;
        
        save_first_snapshot = 1; % Auxiliary variable to save snapshot at t = 0 too!
        
        % Our simulator detects the structural failure of a lightsail as the moment when the stress/strain in any single mesh
        % edge exceeds the tensile strength of the material, after which the failed edge is omitted (mechanically) in future
        % shape calculations.  It is useful to visualize the shape of the sail at and shorly after the onset of such failures.
        % The following settings determine when extra "failure" snapshots are rendered in the snapshot window.
        nextExplosionSnapshot = 5;  % force a snapshot when number of edge failures reaches this threshold
        explosionSnapshotRatio = 10;   %and thereafter increase the forced-snapshot threshold by this ratio
        
        %% Simulation termination conditions
        % We don't want to waste time simulating lightsails that have veered off course or failed structurally.  Here are the
        % early termination criteria that will stop a simulation immediately:
        flyaway_ratio = 50*2.5;  % by how many times the initial sail radius the sail can be displaced (from the beam center) before the simulation stops.  Set to inf to disable.
        explosion_ratio = 6;  % by how many times the initial sail radius the tattered remains of the sail can extend (radially) before the simulation stops.  Set to inf to disable.
        failure_ratio = .00001; %the fraction of triangles broken before termination.  Set >1 to disable
        
        % We can also opt to terminate the simulation after a certain time delay following certain situations.  The delay lets
        % us continue to visualize the mode of failure for a brief time before stopping the simulation.  
        terminate_on_upsidedown = 0;  %initiates a delayed termination if any of the triangles are upside down, which usually means the launch has failed.
        terminate_on_breakup = 0; %initiates a delayed termination upon the first tensile failure of the mesh.
        termination_delay_s = 0.050; %seconds of additional simulation time before any auto-shutoff conditions above terminate the simulation
        
        %% END OF COMMONLY CHANGED SETTINGS. 
        % We don't recommend changing values beyond this point unless you know what you're doing.
        
        %% Set up misc. constants, control variables, and initial temps...  you shouldn't need to modify these values.
        frame_interval = max(floor(frame_interval_s/dt),1);  %number of loops between movie frames
        snapshot_intvl = max(floor(snapshot_interval_s/dt),1);  % number of loops between snapshots
        sim_dur_nt = ceil( (sim_dur + spinup_time) / dt );  %total number of iterations (time steps)
        sim_dur_nto = ceil(sim_dur_nt / sim_output_downsample);
        
        sim_dur_nto2 = ceil(sim_dur_nt / sim_angles_downsample);

        dnto = sim_output_downsample * dt;
        dnto2 = sim_angles_downsample * dt;

        tt = -spinup_time;  % initial physical time counter
        if tt > -dt
            tt=-dt;  % give us one frame to trigger t0t before beginning sim
        end
        uhohs = 0;  % keeps track of number of edge failures
        tensilefailures = 0;
        thermalfailures = 0;
        pf = 1;  % this variable is used to decide when to plot frames for video export
        nf = 1;  % number of video frames rendered
        ps = 1;  % this variable is used to decide when to plot the mesh shape in the 'snapshot' window.
        ns = 1;  % number of snapshots included in the snapshot display
        
        tscale_min = t0-1;  %keeps track of minimum temperature for temperature colormaps
        tscale_max = t0+1;  %keeps track of maximum temperature for temperature colormaps
        max_strain = 0;  % keeps track of max strain for strain colormaps
        min_strain = 0;  % keeps track of min strain for strain colormaps
        max_abs_pwr_dens = 0; % gives consistent color range for the absorbed power density colormap
                
        %% Display settings (for log/diary)
        dispvar2('I0','Irad');
        dispvar3('spinup_gs','spinup_time','videoMode');
        dispvar3('t_ramp_delay','t_ramp_dur','t_on_dur');
        dispvar3('dt','frame_interval','snapshot_intvl');
        dispvar2('sim_dur','sim_dur_nto');
        if (sim_dur_nto > 1e7)
            error('Number of output steps exceeds 10M.  Probably mistake?');
        end
        
        %% Setup video window
        if (videoMode < 2)
%             V1 = VideoWriter([ filebasename '_video' ],'MPEG-4');  %alternate format for Ramon:  'Archival'
%             V1.Quality = 100;

            V1 = VideoWriter(fullfile(path_results,[ filebasename '_video' ]),'Archival');  %alternate format for Ramon:  'Archival'

            open(V1);
        end
        if (videoMode==1)
            fSimVid = figure('pos',[20 50 1200 800]);
        end
        if videoMode==0
            fSimVid = figure('pos',[100 100 2000 800]);
        end
        if videoMode<2
            set(gcf,'color','w');     
        end
        
        %% Set up simulation abort button
         FS = stoploop({'Press this button' 'to terminate simulation.'}) ;
        
        %% Setup snapshot window
        fSnapshot = figure('pos',[200 100 400 800]);  set(gcf,'color','w');
        xbox = 1.10*(radiusmm + sqrt(t0t_x0^2+t0t_y0^2));
        zbox = 0.2*xbox + xbox*zAspectRatio/2; %(0.8*radiusmm*max(1/(pD/radialRings),0.8)) * (1+1/radialRings);
        
        if (pD == 0) 
            zbox = xbox*0.8; 
        end
        
        plotSnapshot;
        hSnapshotLight2 = lightangle(00,20);
        hSnapshotLight = lightangle(180,40);
        axis equal;
%         colormap jet;

        % Scientific color map from https://www.fabiocrameri.ch/colourmaps/
        load('lajolla.mat');
        colormap(flipud(lajolla));

        colorbar ;
        
        %% Initialize monitors:
        comxs = zeros(1,sim_dur_nto);          % COM of whole (remaining) structure
        comys = zeros(1,sim_dur_nto);
        comzs = zeros(1,sim_dur_nto);
        comrs = zeros(1,sim_dur_nto);          % Distance between beam centerline and the COM (mm)
        velxs = zeros(1,sim_dur_nto);          % Velocity of COM
        velys = zeros(1,sim_dur_nto);
        velzs = zeros(1,sim_dur_nto);
        acczs = zeros(1,sim_dur_nto);          % Z-accel of COM
        accts = zeros(1,sim_dur_nto);          % Overall acceleration of COM (magnitude)
        maxvs = zeros(1,sim_dur_nto);          % Velocity of the node having greatest relative velocity to the COM
        maxts = zeros(1,sim_dur_nto);          % Temperature of hottest node
        mints = zeros(1,sim_dur_nto);
        avgts = zeros(1,sim_dur_nto);          % Average temperature of remaining structure (mass weighted)
        vertxs = zeros(1,sim_dur_nto);         % X position of vertex node
        vertys = zeros(1,sim_dur_nto);         % Y position of vertex node
        vertzs = zeros(1,sim_dur_nto);
        maxstrains = zeros(1,sim_dur_nto);     % Strain of most-tensioned edge
        minstrains = zeros(1,sim_dur_nto);     % Strain of most-compressed edge
        maxstrainsnorm = zeros(1,sim_dur_nto); % Strain of most-tensioned edge, normalized to tensile failure limit
        potenergies = zeros(1,sim_dur_nto);    % Total potential energy of the mesh (sum of 1/2*k*dx^2 over all remaining edges)
        kinenergies = zeros(1,sim_dur_nto);    % Total kinetic energy of all vertices (including broken ones?)
        kinzenergies = zeros(1,sim_dur_nto);   % z-direction energy of COM (1/2 * mass * com_velocity^2
        kinoenergies = zeros(1,sim_dur_nto);   % Difference between above two KEs = rotational and vibrational KE.
        ttlinpwr = zeros(1,sim_dur_nto);       % Total input power, sum of all power incident on all non-broken triangles (including upside down ones)
        ttlabspwr = zeros(1,sim_dur_nto);      % Total absorbed power, sum across non-broken triangles
        ttlradpwr = zeros(1,sim_dur_nto);      % Total radiated power, sum across non-broken triangles
        graphuos = zeros(1,sim_dur_nto);       % Number or broken edges
        graphtensfails = zeros(1,sim_dur_nto);  % num of tensile failures
        graphthermfails = zeros(1,sim_dur_nto);  % num of thermal failures 
        graphts = zeros(1,sim_dur_nto);        % simulation time at each recording
        graphts2 = zeros(1,sim_dur_nto2);        % simulation time at each recording for angles
        upsidedowntris = zeros(1,sim_dur_nto);  % number of upside-down triangles (surface normal has negative Z)
        mintrinormz = zeros(1,sim_dur_nto);    % Z component of most downward-facing triangle surface norm
        areas = zeros(1,sim_dur_nto);          % Sum of all non-broken triangle areas
        areasxy = zeros(1,sim_dur_nto);        % Sum of all non-broken triangles projected area area in illumination plane (typically in xy plane unless light is tilted).  Upside down triangles count as negative area, I think.
        areasrightsideupxy = zeros(1,sim_dur_nto);  % Same as above but excluding triangles having surface norms with negative z ("upside down")
        monitorqty = zeros(1,sim_dur_nto);     % misc. monitor vector, for various debug / tracking
        monitorcs = zeros(length(meshcs),sim_dur_nto);      % misc monitor for monitoring a mesh variable along a cross section of the mesh.  meshcs is defined in generateMesh.
        vertextemps = zeros(1,sim_dur_nto);    % Temperature of vertex node
        xradpttemps = zeros(1,sim_dur_nto);    % Temperature of last node (which is always somewhere on the edge)
        xradptxyzs = zeros(3, sim_dur_nto);     % XYZ coordiantes of the last node
        xradptvxyzs = zeros(3, sim_dur_nto);    % velocity of the last node (vx vy vz)
        yradpttemps = zeros(1,sim_dur_nto);    % Temperature of last node (which is always somewhere on the edge)
        yradptxyzs = zeros(3, sim_dur_nto);     % XYZ coordiantes of the last node
        yradptvxyz = zeros(3, sim_dur_nto);    % velocity of the last node (vx vy vz)
        nxradpttemps = zeros(1,sim_dur_nto);    % Temperature of last node (which is always somewhere on the edge)
        nxradptxyzs = zeros(3, sim_dur_nto);     % XYZ coordiantes of the last node
        nxradptvxyz = zeros(3, sim_dur_nto);    % velocity of the last node (vx vy vz)
        nyradpttemps = zeros(1,sim_dur_nto);    % Temperature of last node (which is always somewhere on the edge)
        nyradptxyzs = zeros(3, sim_dur_nto);     % XYZ coordiantes of the last node
        nyradptvxyz = zeros(3, sim_dur_nto);    % velocity of the last node (vx vy vz)
        hoop12xyz = zeros(3, sim_dur_nto);     % XZY coordinates of the hoop node furthest from the hoop COM
        hoop6xyz = zeros(3, sim_dur_nto);      % XYZ coord of the node opposite the hoop12 node (just by adding 1/2 to the ring index)
        hoop9xyz = zeros(3, sim_dur_nto);      % XYZ coord of node opposite hoop3 node
        hoop3xyz = zeros(3, sim_dur_nto);      % XYZ coord of the hoop node closest to hoop COM.  Note that hoop3 and hoop12 aren't always spaced 90 degrees apart on the ring!
        hooprads = zeros(4, sim_dur_nto);      % This is the "radius" (distance to hoop COM) of all 4
        hoopcoms = zeros(3, sim_dur_nto);      % This is the hoop COM
        hoop12vxyz = zeros(3, sim_dur_nto);    % Velocity of hoop12 node
        hoop3vxyz = zeros(3, sim_dur_nto);     % Vel of hoop3 node
        angVels = zeros(3, sim_dur_nto);       % angular velocity of mesh structure
        angVelsRigid = zeros(3, sim_dur_nto);
        torques = zeros(3, sim_dur_nto);       % torque applied to mesh structure (rigid sims only)
        spinAngs = zeros(1,sim_dur_nto);       % keeps track of relative rotation of spin-stabilized sail.  Polar angle of a particular point on edge of sail. 
        sailNorms = zeros(3,sim_dur_nto);      % center-of-sail normal vectors, average of norms of all tris bordering center node.
        sailPitches = zeros(1,sim_dur_nto);    % 'pitch' of center-of-sail, defined as zero being upright, and 90 degrees being edge-on.  
        ofxs = zeros(1,sim_dur_nto);           % total optical force in X dir
        ofys = zeros(1,sim_dur_nto);           % total optical force in y dir
        ofzs = zeros(1,sim_dur_nto);           % total optical force in z dir
        nofxs = zeros(1,sim_dur_nto);           % normalized optical force in x dir (for comparison to Ramon's RB code)
        nofys = zeros(1,sim_dur_nto);           % normalized optical force in x dir 
        nofzs = zeros(1,sim_dur_nto);           % normalized optical force in x dir 
        tqxs = zeros(1,sim_dur_nto);           % rigid body torque about x
        tqys = zeros(1,sim_dur_nto);           % rigid body torque about y
        tqzs = zeros(1,sim_dur_nto);           % rigid body torque about z
        ntqxs = zeros(1,sim_dur_nto);           % normalized rigid body torque about x (for comparison to Ramon's RB code)
        ntqys = zeros(1,sim_dur_nto);           % rigid body torque about y
        ntqzs = zeros(1,sim_dur_nto);           % rigid body torque about z

        tab_pitch = zeros(max(size(t_a)),sim_dur_nto2);  % histogrammed pitch angles
        tab_roll = zeros(max(size(t_a)),sim_dur_nto2);  % histogrammed pitch angles
        tab_yaw = zeros(max(size(t_a)),sim_dur_nto2);  % histogrammed pitch angles

        tab_n_x = zeros(length(n_x),sim_dur_nto2);
        tab_n_y = zeros(length(n_y),sim_dur_nto2);
        tab_n_z = zeros(length(n_z),sim_dur_nto2);

        tab_t_t = zeros(length(t_t),sim_dur_nto2);

        db_bm_timer_gpu = 0;
        db_bm_timer_vec = 0;
        db_bm_timer_itr = 0;
        
        monitor_depth = floor(5e8/length(n_x)/2)*2;  % Number of full-grid monitor records to try to fit in memory.
        if fft_mode
            monitorfzs = zeros(length(n_z),monitor_depth);     %Mech. force on each node, in z direction
            monitorrfs = zeros(length(n_z),monitor_depth);     %mech force on each node, projected to
            monitortfs = zeros(length(n_z),monitor_depth);
            monitorxs = zeros(length(n_x),monitor_depth);      %Complete position monitor,
            monitorys = zeros(length(n_y),monitor_depth);      % All nodes
            monitorzs = zeros(length(n_z),monitor_depth);      % This is where all the memory has gone!
        end
        
        %% Initialize buffers:
        t_optFcn_r2 = (zeros(size(t_na)));  % intermediate variable for calculating gaussian illumination
        t_optFcn_r2_tmp = (zeros(size(t_optFcn_r2)));

        t_ev1 = zeros(3,numtris); % edge vector for first triangle edge
        t_ev1_tmp = zeros(size(t_ev1));

        t_ev2 = zeros(3,numtris);  % edge vector for second edge
        t_ev2_tmp = zeros(size(t_ev2));

        t_ev1n = t_ev1; % normalized (unit) vector for #1 triangle edge
        t_ev1n_tmp = t_ev1n;

        t_ev2n = t_ev2; % normalized (unit) vector for #2 triangle edge
        t_ev2n_tmp = t_ev2n;

        t_norms = (zeros(3,length(t_na)));   % triangle normal vectors, normalized
        t_norms_tmp = zeros(size(t_norms));

        t_IDir = ([ zeros(size(t_na)); zeros(size(t_na)); ones(size(t_na)) ]);   % light incidence direction at each triangle centroid
         % note:  This initial value is not ever changed in the current version of the simulation, because we always assume light
         % oriented on the positive Z axis.  
       
        t_IMag = (zeros(size(t_na)));  % magnitude of light intensity at each triangle centroid (watts per area)
        t_IMag_tmp = zeros(size(t_IMag));

        t_Ipwr = (zeros(size(t_na)));  % incident light power on each triangle (watts)
        t_Ipwr_tmp = zeros(size(t_Ipwr));

        t_aproj = t_a0xy;  % triangle area in light incidence plane, i.e., projected area in xy plane
        t_aproj_tmp = zeros(size(t_aproj));

        t_specA = zeros(size(t_na));  % specular absorbed power
        t_specA_tmp = zeros(size(t_specA));

        t_specR = zeros(size(t_na));  % specular reflection
        t_specR_tmp = zeros(size(t_specR));

        t_abs = (zeros(size(t_na)));  %absorbed optical power in each triangle (watts)
        t_abs_tmp = zeros(size(t_abs));

        t_abs_mr = (zeros(size(t_na)));  %absorbed optical power in each triangle (watts)
        t_abs_mr_tmp = zeros(size(t_abs_mr));

        t_reflp = (zeros(size(t_na)));  %reflected optical power in each triangle (watts)
        t_reflp_tmp = zeros(size(t_reflp));

        t_refln = zeros(3,length(t_na)); % reflected light direction (unit vector)
        t_refln_tmp = zeros(size(t_refln));

        t_oth = zeros(3,length(t_na));  % optical thrust vectors (force, newtons)
        t_oth_tmp = zeros(size(t_oth));

        t_oth_mr = zeros(3,length(t_na)); % optical thust vectors, multi-reflection (force, newtons)
        t_oth_mr_tmp = zeros(size(t_oth_mr));

        t_texn = (zeros(3,length(t_na)));  % unit vector within each triangle plane, indicating orientation axis of texture ("parallel" direction)
        t_texn_tmp = zeros(size(t_texn));

        t_textn = (zeros(3,length(t_na)));  % unit vector within each triangle plane, perpendicular to the texture orientation ("transverse" direction)
        t_textn_tmp = zeros(size(t_textn));

        t_press_n = (zeros(size(t_na)));  %LUT photon pressure coefficient, normal to triangle
        t_press_n_tmp = zeros(size(t_press_n));

        t_press_p = (zeros(size(t_na)));  %LUT photon pressure coefficient, parallel to triangle
        t_press_p_tmp = zeros(size(t_press_p));

        t_press_t = (zeros(size(t_na)));  %LUT photon pressure coefficient, transverse to triangle
        t_press_t_tmp = zeros(size(t_press_t));

        t_lut_a = (zeros(size(t_na)));    %LUT absorption coefficient for each triangle
        t_lut_a_tmp = zeros(size(t_lut_a));

        e_dl = (zeros(1,length(e_na)));  % mm, 
        e_dl_tmp = zeros(size(e_dl));

        e_nl = (zeros(3,length(e_na)));  % normalized (unit) direction vector of each edge
        e_nl_tmp = zeros(size(e_nl));

        e_l = (zeros(size(e_na))); %mm, current edge length
        e_l_tmp = zeros(size(e_l));
        
        e_l0t = zeros(size(e_na)); %mm, current unstrained edge length, including thermal expansion
        e_l0t_tmp = zeros(size(e_l0t));

        % already initialized in generateMesh:  n_a = (zeros(size(n_x)));  % mm2, effective radiation area for each node, for radiative cooling calcs         
        n_abs = (n_a);  % absorbed power (watts), from propulsion beam, heating input to each node
        n_abs_tmp = zeros(size(n_abs));

        n_ems = (n_a);  % emitted power (watts), from radiative cooling, heating output from each node
        n_ems_tmp = zeros(size(n_ems));

        n_hf = (n_a);  % heat flow (watts), from thermal conduction, heat into or out of each node
        n_hf_tmp = zeros(size(n_hf));
        
        %% Initialize various simulation variables, state variables, and flags prior to beginning time-domain simulation
        % There are no user-configurable settings in this section.  Just initialization stuff.
        nto = 0;
        nto2 = 0; % for recording angles
        rot0 = 0;
        accz = 0; acct = 0;
        anybroken = 0;
        nto_firstBroken = sim_dur_nto;
        nto_startAccel = 0; % index of first output (downsampled) data frame for acceleration phase.  This is a monitor, not a setting.
        nt_startAccel = 0;  % iteration step when light first turns on.  This is a monitor, not a setting. 
        time_frstBrkn = inf;
        pwrRampVal = 0;

        % Potential & kinetic energy
        PE = 0; KE = 0; KE_COM = 0; KEo = 0;
        end_reason = 'launched, simulation timeout';
      
        % Needed to decide when to switch from explicit Euler integration to
        % velocity Verlet algorithm, and immediately after that, skipping the first
        % time calculation of the velocities
        %MK:  This is currently unused, but hopefully we will re-enable Verlet stepping in the future.
        first_switch_integrator = 1;
        
        usercancel = 0;  %keeps track of the user's cancelling of the simulation via the dialog box.
        
        timers = zeros(1,15);  % this is where we keep track of the time spent on various simulation calcuation steps, for
        % performance/optimizaiton purposes.  Here are the bins:
        timerBins = { 'spin/incd', ...              1
                      'tri/spec', ...          2
                      'LUT', ...             3
                      'T2N', ...         4
                      'edge_calcs', ...       5
                      'E2N', ...        6
                      'break_test', ...           7
                      'iterate', ...             8
                      'monitor', ...            9
                      'gatherMsh', ...           10 
                      'vidRender', ....        11
                      'Snapshot', ...        12
                      'Termination', ...    13
                      'RayTrace'               ,... 14
                      'waitGPU' };            % 15
        
        terminate_t = 0;  % (seconds) value for delayed termination
        lastdispnt = 0;
        lastoutofboundtswarn = -1; % set to -1 to disable console warnings about upside down triangles
        
        debug_nt_for_pause = -1;
        
        %% Render / record first video frame
        if videoMode < 2
            figure(fSimVid);
        end

        plotMesh_v18_RG;
        hold on
        
        if (videoMode < 2)
            writeVideo(V1,getframe(gcf));
        end
       
        %% Start timers for simulation loop
        tic;
        startsimtic = tic;
        lastdisptic = startsimtic;
        setup_time = toc(startlooptic);
        lastbuttoncheck = tic;
        
        %debug
        specialtimer1 = 0;
        specialtimer2 = 0;
        
        %error check
        if max(t_tex) > 0
            if RayTracingMode
                error('At least one triangle has a MEPS/LUT texture assignment, which is incompatible with RayTracingMode.  Either disable raytracing, or assign specular behavior throughout the mesh.');
            end
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Main time-domain simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nt = 1;

        while nt <= sim_dur_nt 
            %% Notes
            %  In addition to state variables, the following derived variables should be correctly established at the
            %  beginning of each iteration:
            %       inertiaT -- intertia tensor
            %       angVel -- mesh-calculated angular velocity (broken)
            %       t_cx...z -- triangle centroid positions
            %       comx...z -- center of mass
            %       n_dx...z -- node positions relative to COM
            %  disp(['nt=' num2str(nt) ]);

            dooptics = ~mod(nt-1, opt_skip_downsample); % if opt_skip_downsample = 1, dooptics is always 1!
            
            %% Numerical integration
            
            % I am commenting out the rigid body calculations for now, as I
            % am not sure about the correct way for more advanced numerical
            % integration

%             for nrb = 1:num_rb
% 
%                 bodymask = n_rb == nrb;
%                 
%                 rb_force(1,nrb) = sum(n_fx(bodymask));
%                 rb_force(2,nrb) = sum(n_fy(bodymask));
%                 rb_force(3,nrb) = sum(n_fz(bodymask));
%                 
%                 rb_torque(1,nrb) = sum( n_rb_dy(bodymask).*n_fz(bodymask) - n_rb_dz(bodymask).*n_fy(bodymask) );
%                 rb_torque(2,nrb) = sum( n_rb_dz(bodymask).*n_fx(bodymask) - n_rb_dx(bodymask).*n_fz(bodymask) );
%                 rb_torque(3,nrb) = sum( n_rb_dx(bodymask).*n_fy(bodymask) - n_rb_dy(bodymask).*n_fx(bodymask) );
%                 
%                 rb_angMom(:,nrb) = rb_angMom(:,nrb) + 1e6.* rb_torque(:,nrb) .* dt;
%                 rb_angVel(:,nrb) = inv(rb_inertiaT(:,:,nrb)) * rb_angMom(:,nrb);
%                 
%                 % iterate the velocity and COM position
%                 rb_vel(:,nrb) = rb_vel(:,nrb) + 1e6 .* rb_force(:,nrb) ./ rb_mass(nrb) .* dt;
%                 rb_COM(:,nrb) = rb_COM(:,nrb) + dt.*rb_vel(:,nrb);
%                 
%                 dangle = rb_angVel(:,nrb) .* dt; % the angle of the dangle?
%                 
%                 ihaterotationmatrices = ... % gpuArray( ...   % does this need to be promoted to gpuArray?
%                     [ cos(dangle(3)) -sin(dangle(3))  0; ...
%                     sin(dangle(3))  cos(dangle(3))  0; ...
%                     0               0               1 ] * ...
%                     ...
%                     [ cos(dangle(2))  0               sin(dangle(2)) ; ...
%                     0               1               0 ; ...
%                     -sin(dangle(2)) 0               cos(dangle(2)) ] * ...
%                     ...
%                     [ 1               0               0 ; ...
%                     0               cos(dangle(1)) -sin(dangle(1)) ; ...
%                     0               sin(dangle(1))  cos(dangle(1)) ] ;  %) ;
%             
%                 % Let's try Oggy's convention for rotation matrices ...
%             %     ihaterotationmatrices = ...
%             %         [   cos(dangle(3))  sin(dangle(3))  0; ...
%             %             -sin(dangle(3)) cos(dangle(3))  0; ...
%             %             0               0               1 ] * ...
%             %         ...
%             %         [   cos(dangle(2))  0               -sin(dangle(2)) ; ...
%             %             0               1               0 ; ...
%             %             sin(dangle(2)) 0               cos(dangle(2)) ] * ...
%             %         ...
%             %         [   1               0               0 ; ...
%             %             0               cos(dangle(1))  sin(dangle(1)) ; ...
%             %             0               -sin(dangle(1)) cos(dangle(1)) ] ;
%             
%                 n_rb_dxyz_new = ihaterotationmatrices * [ n_rb_dx(bodymask); n_rb_dy(bodymask); n_rb_dz(bodymask) ];
%                 
%                 n_rb_dx(bodymask) = n_rb_dxyz_new(1,:);  %save these so we don't have to recalculate rb_dxyz later
%                 n_rb_dy(bodymask) = n_rb_dxyz_new(2,:);
%                 n_rb_dz(bodymask) = n_rb_dxyz_new(3,:);
%                 
%                 %Now we calculate the new node relative velocities due to rotational motion.  These will be added to the
%                 %body COM velocity to set the new node velocity shortly
%                 n_vx_rot = ( rb_angVel(2,nrb).*n_rb_dxyz_new(3,:) - rb_angVel(3,nrb).*n_rb_dxyz_new(2,:) );
%                 n_vy_rot = -( rb_angVel(1,nrb).*n_rb_dxyz_new(3,:) - rb_angVel(3,nrb).*n_rb_dxyz_new(1,:) );
%                 n_vz_rot = ( rb_angVel(1,nrb).*n_rb_dxyz_new(2,:) - rb_angVel(2,nrb).*n_rb_dxyz_new(1,:) );
%                 
%                 % calculate the new node absolute positions for this body
%             
%                 n_x_new = rb_COM(1,nrb) + n_rb_dxyz_new(1,:);
%                 n_y_new = rb_COM(2,nrb) + n_rb_dxyz_new(2,:);
%                 n_z_new = rb_COM(3,nrb) + n_rb_dxyz_new(3,:);
%                 
%                 n_vx_new = rb_vel(1,nrb) + n_vx_rot;
%                 n_vy_new = rb_vel(2,nrb) + n_vy_rot;
%                 n_vz_new = rb_vel(3,nrb) + n_vz_rot;
%             
%                 n_x(bodymask) = n_x_new;
%                 n_y(bodymask) = n_y_new;
%                 n_z(bodymask) = n_z_new;
%                     
%                 n_vx(bodymask) = n_vx_new;
%                 n_vy(bodymask) = n_vy_new;
%                 n_vz(bodymask) = n_vz_new;
%                 
%             end

            % Flexible dynamics
            if any_flex
                
                % Save original positions
                n_x_orig = n_x;
                n_y_orig = n_y;
                n_z_orig = n_z;

                % Save original time
                tt_orig = tt;

                % Save original node temperatures
                n_t_orig = n_t;

                if ~strcmp(integrationMethod,'VVerlet')
    
                    % Calculate forces at tt
                    core_simulate_first_v1;
                    
                end

                % NOTE about units of velocities and positions: positions
                % are in units of mm. Given that masses are in units of mg,
                % time in units of s and forces in units of N, in order for
                % velocity to have a unit of mm/s, acceleration will need
                % to be multiplied by a factor of 1e3 * 1e3 = 1e6

                % Semi-implicit (symplectic) Euler method
                if strcmp(integrationMethod,'SIEuler')

                    % Iterate the node velocities
                    n_vx(n_flexible) = n_vx(n_flexible) + 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                    n_vy(n_flexible) = n_vy(n_flexible) + 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                    n_vz(n_flexible) = n_vz(n_flexible) + 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;
                
                    % Iterate the node positions
                    n_x(n_flexible) = n_x(n_flexible) + n_vx(n_flexible) .* dt;
                    n_y(n_flexible) = n_y(n_flexible) + n_vy(n_flexible) .* dt;
                    n_z(n_flexible) = n_z(n_flexible) + n_vz(n_flexible) .* dt;

                    % Iterate node, edge and triangle temperatures
                    n_t = n_t + (n_hf + n_abs - n_ems).*dt ./ n_heatmass;   % todo:  need to support variable heat cap across the array

                    % Calculate edge temperature as averaged temperature of
                    % node temperatures forming the edge
                    e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;

                    % Calculate triangle temperature as average of
                    % temperatures of nodes making up the triangle
                    t_t  =  ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) ./ 3;

                    % Update centroids
                    t_cx = ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
                    t_cy = ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
                    t_cz = ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;


               elseif strcmp(integrationMethod,'MP')

                    % NOT DONE YET!

                    n_fx = n_fx_tmp;
                    n_fy = n_fy_tmp;
                    n_fz = n_fz_tmp;

                    % Calculate velocity at midpoint tt + dt/2

                    k2_x(n_flexible) = n_vx(n_flexible) + 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* (dt/2);
                    k2_y(n_flexible) = n_vy(n_flexible) + 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* (dt/2);
                    k2_z(n_flexible) = n_vz(n_flexible) + 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* (dt/2);

                    % Calculate position at tt + dt/2
                    n_x(n_flexible) = n_x(n_flexible) + n_vx(n_flexible) .* (dt/2);
                    n_y(n_flexible) = n_y(n_flexible) + n_vy(n_flexible) .* (dt/2);
                    n_z(n_flexible) = n_z(n_flexible) + n_vz(n_flexible) .* (dt/2);

                    % Calculate forces at midpoint tt + dt/2 with position
                    % (and velocity) at tt + dt/2
                    tt = tt_orig + dt/2;

                    assign_variables_for_core_simulate;

                    core_simulate_v1;

                    get_variables_from_core_simulate;

                    n_fx = n_fx_tmp;
                    n_fy = n_fy_tmp;
                    n_fz = n_fz_tmp;

                    % Compute position at tt + dt with velocity at tt + dt/2
                    n_x(n_flexible) = n_x(n_flexible) + k2_x(n_flexible) .* dt;
                    n_y(n_flexible) = n_y(n_flexible) + k2_y(n_flexible) .* dt;
                    n_z(n_flexible) = n_z(n_flexible) + k2_z(n_flexible) .* dt;

                    % Compute velocity at tt + dt with acceleration at tt + dt/2
                    n_vx(n_flexible) = n_vx(n_flexible) + 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                    n_vy(n_flexible) = n_vy(n_flexible) + 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                    n_vz(n_flexible) = n_vz(n_flexible) + 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;

                elseif strcmp(integrationMethod,'VVerlet')

                    % NOT DONE YET!

                    % Calculate forces at tt
                    if nt == 1

                        core_simulate_v1;

                        n_fx = n_fx_tmp;
                        n_fy = n_fy_tmp;
                        n_fz = n_fz_tmp;

                    elseif nt > 1

                        n_fx = n_fx_next;
                        n_fy = n_fy_next;
                        n_fz = n_fz_next;
                        
                    else
                        error('Something went wrong with the time stepping!')
                    end

                    n_x(n_flexible) = n_x(n_flexible) + n_vx(n_flexible) .* dt + ...
                        0.5 .* (dt).^2 * n_fx(n_flexible)./n_m(n_flexible);
                    n_y(n_flexible) = n_y(n_flexible) + n_vy(n_flexible) .* dt + ...
                        0.5 .* (dt).^2 * n_fy(n_flexible)./n_m(n_flexible);
                    n_z(n_flexible) = n_z(n_flexible) + n_vz(n_flexible) .* dt + ...
                        0.5 .* (dt).^2 * n_fz(n_flexible)./n_m(n_flexible);

                    % Calculate forces at tt + dt
                    tt = tt_orig + dt;

                    assign_variables_for_core_simulate;

                    core_simulate_v1;

                    n_fx_next = n_fx_tmp;
                    n_fy_next = n_fy_tmp;
                    n_fz_next = n_fz_tmp;

                    n_vx(n_flexible) = n_vx(n_flexible) + 0.5 .* dt .* ...
                        (n_fx(flexible)./n_m(n_flexible) + n_fx_next(flexible)./n_m(n_flexible));
                    n_vy(n_flexible) = n_vy(n_flexible) + 0.5 .* dt .* ...
                        (n_fy(flexible)./n_m(n_flexible) + n_fy_next(flexible)./n_m(n_flexible));
                    n_vz(n_flexible) = n_vz(n_flexible) + 0.5 .* dt .* ...
                        (n_fz(flexible)./n_m(n_flexible) + n_fz_next(flexible)./n_m(n_flexible));

                elseif strcmp(integrationMethod,'RK4')

                    % Finishing first RK4 evaluation

                    k1_vx = 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                    k1_vy = 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                    k1_vz = 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;

                    k1_x = n_vx(n_flexible) .* dt;
                    k1_y = n_vy(n_flexible) .* dt;
                    k1_z = n_vz(n_flexible) .* dt;

                    k1_n_t = (n_hf + n_abs - n_ems).*dt ./ n_heatmass;

                    % 2ND RK4 EVALUATION

                    % Step time by half of original time step
                    tt = tt_orig + dt/2;
                    dt_tmp = dt/2;

                    % Update position at this new time step
                    % Note: no need to step velocity just now, as force
                    % calculations are independent of velocity
                    n_x(n_flexible) = n_x_orig(n_flexible) + k1_x/2;
                    n_y(n_flexible) = n_y_orig(n_flexible) + k1_y/2;
                    n_z(n_flexible) = n_z_orig(n_flexible) + k1_z/2;

                    % Update node, edge and triangle temperatures
                    n_t = n_t_orig + k1_n_t/2;                    

                    % Calculate forces at new time step with updated
                    % positions
                    assign_variables_for_core_simulate;
                    core_simulate_v1;
                    get_variables_from_core_simulate;

                    n_fx = n_fx_tmp;
                    n_fy = n_fy_tmp;
                    n_fz = n_fz_tmp;

                    k2_vx = 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                    k2_vy = 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                    k2_vz = 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;

                    k2_x = (n_vx(n_flexible) + 0.5 .* k1_vx(n_flexible)) .*dt;
                    k2_y = (n_vy(n_flexible) + 0.5 .* k1_vy(n_flexible)) .*dt;
                    k2_z = (n_vz(n_flexible) + 0.5 .* k1_vz(n_flexible)) .*dt;
                    
                    k2_n_t = (n_hf_tmp + n_abs_tmp - n_ems_tmp).*dt_tmp ./ n_heatmass;

                    % 3RD RK4 EVALUATION

                    % Step time by half of original time step (superfluous,
                    % but just to be safe ...)
                    tt = tt_orig + dt/2;
                    dt_tmp = dt/2;

                    % Update position at this new time step
                    % Note: no need to step velocity just now, as force
                    % calculations are independent of velocity
                    n_x(n_flexible) = n_x_orig(n_flexible) + k2_x/2;
                    n_y(n_flexible) = n_y_orig(n_flexible) + k2_y/2;
                    n_z(n_flexible) = n_z_orig(n_flexible) + k2_z/2;

                    % Update node, edge and triangle temperatures
                    n_t = n_t_orig + k2_n_t/2;  

                    % Calculate forces at new time step with updated
                    % positions
                    assign_variables_for_core_simulate;
                    core_simulate_v1;
                    get_variables_from_core_simulate;

                    n_fx = n_fx_tmp;
                    n_fy = n_fy_tmp;
                    n_fz = n_fz_tmp;

                    k3_vx = 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                    k3_vy = 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                    k3_vz = 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;

                    k3_x = (n_vx(n_flexible) + 0.5 .* k2_vx(n_flexible)) .*dt;
                    k3_y = (n_vy(n_flexible) + 0.5 .* k2_vy(n_flexible)) .*dt;
                    k3_z = (n_vz(n_flexible) + 0.5 .* k2_vz(n_flexible)) .*dt;

                    k3_n_t = (n_hf_tmp + n_abs_tmp - n_ems_tmp).*dt_tmp ./ n_heatmass;

                    % 4TH RK4 EVALUATION

                    % Step time by original time step
                    tt = tt_orig + dt;
                    dt_tmp = dt;

                    % Update position at this new time step
                    % Note: no need to step velocity just now, as force
                    % calculations are independent of velocity
                    n_x(n_flexible) = n_x_orig(n_flexible) + k3_x;
                    n_y(n_flexible) = n_y_orig(n_flexible) + k3_y;
                    n_z(n_flexible) = n_z_orig(n_flexible) + k3_z;

                    % Update node, edge and triangle temperatures
                    n_t = n_t_orig + k3_n_t;

                    % Calculate forces at new time step with updated
                    % positions
                    assign_variables_for_core_simulate;
                    core_simulate_v1;
                    get_variables_from_core_simulate;

                    n_fx = n_fx_tmp;
                    n_fy = n_fy_tmp;
                    n_fz = n_fz_tmp;

                    k4_vx = 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                    k4_vy = 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                    k4_vz = 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;

                    k4_x = (n_vx(n_flexible) + k3_vx(n_flexible)) .*dt;
                    k4_y = (n_vy(n_flexible) + k3_vy(n_flexible)) .*dt;
                    k4_z = (n_vz(n_flexible) + k3_vz(n_flexible)) .*dt;

                    k4_n_t = (n_hf_tmp + n_abs_tmp - n_ems_tmp).*dt_tmp ./ n_heatmass;

                    % Update velocities & positions at new time step
                    n_vx(n_flexible) = n_vx(n_flexible) + (1/6) .* ...
                        (k1_vx + 2.*k2_vx + 2.*k3_vx + k4_vx);
                    n_vy(n_flexible) = n_vy(n_flexible) + (1/6) .* ...
                        (k1_vy + 2.*k2_vy + 2.*k3_vy + k4_vy);
                    n_vz(n_flexible) = n_vz(n_flexible) + (1/6) .* ...
                        (k1_vz + 2.*k2_vz + 2.*k3_vz + k4_vz);

                    n_x(n_flexible) = n_x_orig(n_flexible) + (1/6) .* ...
                        (k1_x + 2.*k2_x + 2.*k3_x + k4_x);
                    n_y(n_flexible) = n_y_orig(n_flexible) + (1/6) .* ...
                        (k1_y + 2.*k2_y + 2.*k3_y + k4_y);
                    n_z(n_flexible) = n_z_orig(n_flexible) + (1/6) .* ...
                        (k1_z + 2.*k2_z + 2.*k3_z + k4_z);

                    % Update node, edge and triangle temperatures at new
                    % time step
                    n_t = n_t_orig + (1/6) .* ...
                        (k1_n_t + 2.*k2_n_t + 2.*k3_n_t + k4_n_t);
                    e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;
                    t_t = ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) ./ 3;

                    % Update centroids
                    t_cx = ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
                    t_cy = ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
                    t_cz = ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;

                end
                
            end

            % Go back to original time step
            tt = tt_orig;
            
            % update overall COM:
            coms_x = n_x .* n_m;
            coms_y = n_y .* n_m;
            coms_z = n_z .* n_m;
            
            comx = sum(coms_x)./totalmass_nodes ;
            comy = sum(coms_y)./totalmass_nodes ;
            comz = sum(coms_z)./totalmass_nodes ;


            %% Advance time step
            tt = tt + dt;
            

            %% Apply initial sail offest and tilts at end of spin-up
            if (tt >= 0) && ( (tt-dt < 0) ) % note: tt has already been incremented at this point...
                disp('Doing initial offset stuff');
                
                rb_COM = zeros(3,num_rb);
                
                if t0t_changeSpinAngle
                    approxSpinAng = atan2(n_y(xradidx)-comy, n_x(xradidx)-comx);
                    t0ta = t0t_s0 - approxSpinAng;
                    
                    n_x_new = cos(t0ta).*n_x - sin(t0ta).*n_y;
                    n_y = sin(t0ta).*n_x + cos(t0ta).*n_y;
                    n_x = n_x_new;
                    
                    n_vx_new = cos(t0ta).*n_vx - sin(t0ta).*n_vy;
                    n_vy = sin(t0ta).*n_vx + cos(t0ta).*n_vy;
                    n_vx = n_vx_new;
                    
                    if num_rb ~= 0
                        rb_COM(1,nrb) = sum(coms_x(bodymask))./rb_mass(nrb) ;
                        rb_COM(2,nrb) = sum(coms_y(bodymask))./rb_mass(nrb) ;
                        rb_COM(3,nrb) = sum(coms_z(bodymask))./rb_mass(nrb) ;

                        temp_dx = n_x(bodymask) - rb_COM(1,nrb);
                        temp_dy = n_y(bodymask) - rb_COM(2,nrb);
                        temp_dz = n_z(bodymask) - rb_COM(3,nrb);
                        
                        n_rb_dx(bodymask) = temp_dx;  % save these because we'll need them for torque calculation
                        n_rb_dy(bodymask) = temp_dy;
                        n_rb_dz(bodymask) = temp_dz;
                    end

                end
                
                if t0t_t0_x ~= 0 || t0t_t0_y ~= 0 || t0t_t0_z ~= 0

                    % Evaluate them first to save time
                    DRM123_11_eval = DRM123_11(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_12_eval = DRM123_12(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_13_eval = DRM123_13(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_21_eval = DRM123_21(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_22_eval = DRM123_22(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_23_eval = DRM123_23(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_31_eval = DRM123_31(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_32_eval = DRM123_32(t0t_t0_z,t0t_t0_y,t0t_t0_x);
                    DRM123_33_eval = DRM123_33(t0t_t0_z,t0t_t0_y,t0t_t0_x);

                    % Transform node coordinates 
                    n_x_new = n_x .*  DRM123_11_eval + ...
                        n_y .* DRM123_21_eval + n_z .* DRM123_31_eval;

                    n_y_new = n_x .* DRM123_12_eval + ...
                        n_y .* DRM123_22_eval + n_z .* DRM123_32_eval;

                    n_z_new = n_x .* DRM123_13_eval + ...
                        n_y .* DRM123_23_eval + n_z .* DRM123_33_eval;

                    n_x = n_x_new;
                    n_y = n_y_new;
                    n_z = n_z_new;

                    % Transform center-of-mass coordinates
                    comx_new = comx .* DRM123_11_eval + ...
                        comy .* DRM123_21_eval + comz .* DRM123_31_eval;

                    comy_new = comx .* DRM123_12_eval + ...
                        comy .* DRM123_22_eval + comz .* DRM123_32_eval;

                    comz_new = comx .* DRM123_13_eval + ...
                        comy .* DRM123_23_eval + comz .* DRM123_33_eval;

                    comx = comx_new;
                    comy = comy_new;
                    comz = comz_new;

                    % Transform velocities
                    n_vx_new = n_vx .*  DRM123_11_eval + ...
                        n_vy .* DRM123_21_eval + n_vz .* DRM123_31_eval;

                    n_vy_new = n_vx .* DRM123_12_eval + ...
                        n_vy .* DRM123_22_eval + n_vz .* DRM123_32_eval;

                    n_vz_new = n_vx .* DRM123_13_eval + ...
                        n_vy .* DRM123_23_eval + n_vz .* DRM123_33_eval;

                    n_vx = n_vx_new;
                    n_vy = n_vy_new;
                    n_vz = n_vz_new;
                    
                    if num_rb ~= 0
                        % Transform center-of-mass coordinates of rigid bodies
                        rb_COM_new(1,:) = rb_COM(1,:) .*  DRM123_11_eval + ...
                            rb_COM(2,:) .* DRM123_21_eval + rb_COM(3,:) .* DRM123_31_eval;

                        rb_COM_new(2,:) = rb_COM(1,:) .*  DRM123_12_eval + ...
                            rb_COM(2,:) .* DRM123_22_eval + rb_COM(3,:) .* DRM123_32_eval;

                        rb_COM_new(3,:) = rb_COM(1,:) .*  DRM123_13_eval + ...
                            rb_COM(2,:) .* DRM123_23_eval + rb_COM(3,:) .* DRM123_33_eval;

                        rb_COM = rb_COM_new;

                        % Transform velocities of rigid bodies
                        rb_vel_new(1,:) = rb_vel(1,:) .*  DRM123_11_eval + ...
                            rb_vel(2,:) .* DRM123_21_eval + rb_vel(3,:) .* DRM123_31_eval;

                        rb_vel_new(2,:) = rb_vel(1,:) .*  DRM123_12_eval + ...
                            rb_vel(2,:) .* DRM123_22_eval + rb_vel(3,:) .* DRM123_32_eval;

                        rb_vel_new(3,:) = rb_vel(1,:) .*  DRM123_13_eval + ...
                            rb_vel(2,:) .* DRM123_23_eval + rb_vel(3,:) .* DRM123_33_eval;

                        rb_vel = rb_vel_new;
                    
                    end

                end
                
                if t0t_x0 ~= 0
                    n_x = n_x + t0t_x0;
                    comx = comx + t0t_x0;
%                     for nrb = 1:num_rb % Ramon: I don't think the
%                     for-loop is needed
                        rb_COM(1,:) = rb_COM(1,:) + t0t_x0;
%                     end
                end

                if t0t_y0 ~= 0
                    n_y = n_y + t0t_y0;
                    comy = comy + t0t_y0;
%                     for nrb = 1:num_rb % Ramon: I don't think the
%                     for-loop is needed
                        rb_COM(2,:) = rb_COM(2,:) + t0t_x0;
%                     end
                end
                
%                 if t0t_t0 ~= 0
%                     error('I haven''t yet written the code for tilting the sail at t0');
%                     %Note to Ramon: When pasting in your tilt code here, you'll need to deal with the rigid bodies as
%                     %well, if you want rigid bodies to work.  Not quite sure how to do that!
%                 end
                
            end
              
            %% Update total body velocity, KE, spin speed...
            % Calculate kinetic energy as KE = (1/2) * m * v^2
            KE = 0.5e-9 * sum( n_m .* (n_vx.^2 + n_vy.^2 + n_vz.^2) ); % Joules
            
            meanvx = sum(n_vx.*n_m)./totalmass_nodes;
            meanvy = sum(n_vy.*n_m)./totalmass_nodes;
            meanvz = sum(n_vz.*n_m)./totalmass_nodes;
            maxvrel = max( sqrt( (n_vx - meanvx).^2 + (n_vy - meanvy).^2 + ...
                (n_vz - meanvz).^2 ) );
            
            KE_COM = 0.5e-9 * totalmass_nodes * norm([meanvz meanvy meanvx]).^2;
            KEo = KE - KE_COM;
            
            COMVel = [meanvx; meanvy; meanvz];
      
            % calculate spin speed during spin-up
            if tt<0
                if any(n_af(1,:))
                    n_rps = sqrt(n_vx.^2+n_vy.^2)./n_arad./pi/2;
                    %n_rps(1)=0; % Note:  angular velocity of first (center/vertex) node and last (tethered payload
                    %node) might be NaN, but we don't need to zero them out here, just omit them from averaging on the
                    %following line:
                    rot0 = mean(n_rps(2:end-1));  % get spin speed as average node rotational velocities for all nodes except first and last
                end
            end
            
            %% Update inertia tensors for rigid bodies
            % Calculate the local positions of nodes (relative to COM) and the
            % inertia tensor matrix
            for nrb=1:num_rb
                bodymask = n_rb == nrb;

                temp_dx2 = n_rb_dx(bodymask).^2;
                temp_dy2 = n_rb_dy(bodymask).^2;
                temp_dz2 = n_rb_dz(bodymask).^2;
                temp_dr = sqrt( temp_dx2 + temp_dy2 + temp_dz2);
    
                rb_inertiaT(:,:,nrb) = [  sum( n_m(bodymask) .* ( temp_dy2 + temp_dz2 ) )    ...
                    sum( -n_m(bodymask) .*  temp_dx .* temp_dy  )         ...
                    sum( -n_m(bodymask) .* temp_dx .* temp_dz );  ...
                    sum( -n_m(bodymask) .* temp_dy .* temp_dx )            ...
                    sum( n_m(bodymask) .* (temp_dx2 + temp_dz2) )  ...
                    sum( -n_m(bodymask) .* temp_dy .* temp_dz );   ...
                    sum( -n_m(bodymask) .* temp_dz .* temp_dx )            ...
                    sum( -n_m(bodymask) .* temp_dz .* temp_dy )        ...
                    sum( n_m(bodymask) .* (temp_dx2 + temp_dy2) ) ] ;

%                     n_rb_dx(bodymask) = temp_dx;  % save these because we'll need them for torque calculation
%                     n_rb_dy(bodymask) = temp_dy;
%                     n_rb_dz(bodymask) = temp_dz;
    
            end
            
            %% Update triangle centroids
            
            % Calculate triangle centroids, where the x, y or z coordinate of a
            % centroid is calculated as the arithmetic mean of the respective x, y
            % or z vertex/nodes coordinates
            t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
            t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
            t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;
            
            
            %% FFT analysis -- record monitor data.
            if fft_mode
                n_xrelcom = n_dx;
                n_yrelcom = n_dy;
                n_rrelcom = sqrt(n_xrelcom.^2+n_yrelcom.^2);  % this is NOT a duplicate of n_dr-- this one is projected in xy plane
                n_xrelcom = n_xrelcom./n_rrelcom;
                n_yrelcom = n_yrelcom./n_rrelcom;
                % n_tan_norm = [ -n_yrelcom ; n_xrelcom ];
                % projet a onto b:  a dot bnorm
                %n_rmf = n_mf(1,:).*n_xrelcom + n_mf(2,:).*n_yrelcom;
                %n_tmf = n_mf(1,:).*(-n_yrelcom) + n_mf(2,:).*n_xrelcom;
                %n_rmf = n_pec(1,:).*n_xrelcom + n_pec(2,:).*n_yrelcom;
                %n_tmf = n_pec(1,:).*(-n_yrelcom) + n_pec(2,:).*n_xrelcom;
                n_rmf = n_mf(1,:).*n_xrelcom + n_mf(2,:).*n_yrelcom;
                n_tmf = n_mf(1,:).*(-n_yrelcom) + n_mf(2,:).*n_xrelcom;
            end
            
            timers(8) = timers(8) + toc;
            tic;
            
            
          %  gather(g_n_x(1));  %this should force us to wait until all GPU calculations have terminated...  not sure how far back things might be queued
            
            timers(15) = timers(15) + toc;
            tic;
            
            %% Record standard output monitors (downsampled)  
            % Note that video/snapshot outputs occur only during monitor output frames, to avoid rendudant gather() steps.  
             if ( ~mod(nt-1, sim_angles_downsample) )
                 
                nto2 = nto2 + 1;
                
                tab_pitch(:,nto2) = pitch_angles;
                tab_roll(:,nto2) = roll_angles;
                tab_yaw(:,nto2) = yaw_angles;

                tab_n_x(:,nto2) = n_x;
                tab_n_y(:,nto2) = n_y;
                tab_n_z(:,nto2) = n_z;

                tab_t_t(:,nto2) = t_t;
                
                graphts2(nto2) = tt;
                
            end
            
            if ( ~mod(nt-1, sim_output_downsample) )  %TODO:  let's use nt-1 so we catch the first loop output?
                
                nto = nto + 1;
                
                % Calculate some values here that don't need to be calculated during the loop:
                % Determine maximum and minimum strain on any of the edges in the mesh (not including broken edges)
                max_strain = max(e_s.*(e_notbroken));
                min_strain = min(e_s.*(e_notbroken));
                max_strain_norm = max(e_s(e_notbroken)./e_tensile_strain_limit(e_notbroken));
            
                
        
                comxs(nto) = comx; comys(nto) = comy; comzs(nto) = comz;
                comrs(nto) = sqrt( (comxs(nto)).^2 + comys(nto).^2);
                velxs(nto) = meanvx; velys(nto) = meanvy; velzs(nto) =meanvz; 
                
                vertxs(nto) = n_x(1);
                vertys(nto) = n_y(1);
                vertzs(nto) = n_z(1);
              
                if (num_rb>0)
                    angVelsRigid(:,nto) = (rb_angVel(:,1));
                else
                    angVelsRigid(:,nto) = angVel;
                end

                angVels(:,nto) = (angVel);
                
                torques(:,nto) = (externalTorque);
                spinAngs(nto) = (atan2(n_dy(xradidx), n_dx(xradidx)));  %disp([num2str(rad2deg(approxSpinAng)) ' sam_wait=' num2str(sam_wait) ] );
                sailNorms(:,nto) = (mean(t_norms(:,n_centerTris),2));
                sailPitches(:,nto) = atan2( sqrt(sailNorms(1,nto).^2 + sailNorms(2,nto).^2), sailNorms(3,nto) );
                
                ofxs(nto) = externalForceX;
                ofys(nto) = externalForceY;
                ofzs(nto) = externalForceZ;
                
                nofxs(nto) = ofxs(nto) / (I0 * ((2*radiusmm).^2) / c0);
                nofys(nto) = ofys(nto) / (I0 * ((2*radiusmm).^2) / c0);
                nofzs(nto) = ofzs(nto) / (I0 * ((2*radiusmm).^2) / c0);
                
                tqxs(nto) = externalTorqueX;
                tqys(nto) = externalTorqueY;
                tqzs(nto) = externalTorqueZ;

                ntqxs(nto) = tqxs(nto) / (I0 * ((2*radiusmm).^3) / c0);
                ntqys(nto) = tqys(nto) / (I0 * ((2*radiusmm).^3) / c0);
                ntqzs(nto) = tqzs(nto) / (I0 * ((2*radiusmm).^3) / c0);
                graphuos(nto) = ( uhohs );
                graphtensfails(nto) = tensilefailures;
                graphthermfails(nto) = thermalfailures;
                graphts(nto) = tt;
                avgts(nto) = (sum(t_t.*t_m))/totalmass_tris;
                maxts(nto) = max(t_t);
                mints(nto) = min(t_t);
                maxstrains(nto) = max_strain;
                maxstrainsnorm(nto) = max_strain_norm; 
                minstrains(nto) = min_strain; 
                potenergies(nto) = PE;
                kinenergies(nto) = KE;
                kinzenergies(nto) = KE_COM;
                kinoenergies(nto) = KEo;
                ttlinpwr(nto)  = (sum(t_Ipwr.*(t_notbroken)));
                ttlabspwr(nto) = (sum(t_abs(t_notbroken)) + sum(t_abs_mr(t_notbroken)));
                ttlradpwr(nto) = (sum(n_ems));
                areas(nto) = (sum(t_a.*(t_notbroken)));
                areasxy(nto) = (sum(t_aproj.*(t_notbroken)));
%                 upsidedowntris(nto) = (upsidedowntris);
                areasrightsideupxy(nto) = (sum( t_aproj.*(t_aproj>0).*(t_notbroken)));
                mintrinormz(nto) = (min( t_norms(3,:) ));
                % monitorqty(nto) = mean(n_rps(2:end));
                xradpttemps(nto) = n_t(end-1);
                vertextemps(nto) = n_t(1);
                
                if(nto==1)
                    accz = velzs(nto) ./ dt;
                    acct = norm( [velxs(nto) velys(nto) velzs(nto) ] ) ./dt;
                else
                    accz = (velzs(nto) - velzs(nto-1))./(dt*sim_output_downsample);
                    acct = norm( [velxs(nto) velys(nto) velzs(nto)] - [velxs(nto-1) velys(nto-1) velzs(nto-1) ] )./(dt*sim_output_downsample);
                end
                acczs(nto)=accz;
                accts(nto)=acct;
                
                
                % monitor evolution of any mesh variable, as evaluated along a specific (typically radial line)
                % across the mesh here.
                monitorcs(:,nto) = (n_z(meshcs)-comz);
                
                % monitor evolution any other variable here:
                monitorqty(nto) = (rot0);
                
                % elongation ratio monitors  (these are all gpu arrays now)
                ringcomx = sum(coms_x(outerringstartidx:outerringendidx))./sum(n_m(outerringstartidx:outerringendidx));
                ringcomy = sum(coms_y(outerringstartidx:outerringendidx))./sum(n_m(outerringstartidx:outerringendidx));
                ringcomz = sum(coms_z(outerringstartidx:outerringendidx))./sum(n_m(outerringstartidx:outerringendidx));
                ringcom = [ringcomx; ringcomy; ringcomz];
                
                ringradii = sqrt( ...
                    (n_x(outerringstartidx:outerringendidx) - ringcomx).^2 + ...
                    (n_y(outerringstartidx:outerringendidx) - ringcomy).^2 + ...
                    (n_z(outerringstartidx:outerringendidx) - ringcomz).^2   );
                
                [furthestrad, furthestidx] = max(ringradii);
                [closestrad, closestidx] = min(ringradii);
                
                
                idx6 = mod(furthestidx + round(numringpts/2) -1, numringpts ) + 1;
                idx3 = mod(furthestidx + round(numringpts*0.25) -1 , numringpts) + 1;
                idx9 = mod(closestidx + round(numringpts/2) -1 , numringpts) + 1;
                
                hoopcoms(:,nto) = (ringcom);
                hooprads(:,nto) = ([ringradii(furthestidx); ringradii(closestidx); ringradii(idx6); ringradii(idx9) ]);
                                
                hoop12xyz(:,nto) = ([ n_x(outerringstartidx+furthestidx-1); n_y(outerringstartidx+furthestidx-1); n_z(outerringstartidx+furthestidx-1) ]);
                hoop6xyz(:,nto)  = ([ n_x(outerringstartidx+idx6-1);        n_y(outerringstartidx+idx6-1);        n_z(outerringstartidx+idx6-1)        ]);
                hoop3xyz(:,nto)  = ([ n_x(outerringstartidx+closestidx-1);  n_y(outerringstartidx+closestidx-1);  n_z(outerringstartidx+closestidx-1)  ]);
                hoop9xyz(:,nto)  = ([ n_x(outerringstartidx+idx9-1);        n_y(outerringstartidx+idx9-1);        n_z(outerringstartidx+idx9-1)        ]);
                
                hoop12vxyz(:,nto) = ([ n_vx(outerringstartidx+furthestidx-1); n_vy(outerringstartidx+furthestidx-1); n_vz(outerringstartidx+furthestidx-1) ] );
                hoop3vxyz(:,nto)  = ([ n_vx(outerringstartidx+closestidx-1);  n_vy(outerringstartidx+closestidx-1);  n_vz(outerringstartidx+closestidx-1) ] );
                
                xradptxyzs(:,nto)  = ([n_x(end);  n_y(end);  n_z(end)] );
                xradptvxyzs(:,nto) = ([n_vx(end); n_vy(end); n_vz(end)] );
            
            
                % keep track of time spent recording the monitor data
                timers(9) = timers(9) + toc;
                tic;

                %% 17.  Bring back all state vectors from gpu for plots, if needed
                % here we are only bringing back the data -- plot functions follow.
                if (nt - pf >= frame_interval - 1)
                    
                    % note that this is probably much more than we will actually use, we can probably comment out the unused
                    % data if the retrieve operations take too long.
                    
                    % gather things needed from gpu
                   % ringradii  = gather(g_ringradii);
                    
                    % node vectors and arrays
%                     n_x = gather(g_n_x);
%                     n_y = gather(g_n_y);
%                     n_z = gather(g_n_z);
%                     n_t = gather(g_n_t);
%                     n_a = gather(g_n_a);
%                     n_af = gather(g_n_af);
%                     n_dr = gather(g_n_dr);
%                     n_dx = gather(g_n_dx);
%                     n_dy = gather(g_n_dy);
%                     n_dz = gather(g_n_dz);
%                     n_dx2 = gather(g_n_dx2);
%                     n_dy2 = gather(g_n_dy2);
%                     n_dz2 = gather(g_n_dz2);
%                     n_dxyz = gather(g_n_dxyz);
%                     n_fx = gather(g_n_fx);
%                     n_fy = gather(g_n_fy);
%                     n_fz = gather(g_n_fz);
%                     n_mf = gather(g_n_mf);
%                     n_of = gather(g_n_of);
%                     n_vx = gather(g_n_vx);
%                     n_vy = gather(g_n_vy);
%                     n_vz = gather(g_n_vz);
%                     comx = gather(comx);
%                     comy = gather(comy);
%                     comz = gather(comz);
%                     n_abs = gather(g_n_abs);  % absorbed power (watts), from propulsion beam, heating input to each node
%                     n_ems = gather(g_n_ems);  % emitted power (watts), from radiative cooling, heating output from each node
%                     n_hf = gather(g_n_hf);  % heat flow (watts), from thermal conduction, heat into or out of each node
%                     
%                     % edge vectors and arrays
%                     e_t = gather(g_e_t);
%                     e_broken = gather(g_e_broken);
%                     e_notbroken = gather(g_e_notbroken);
%                     e_dl = gather(g_e_dl);  % mm,
%                     e_nl = gather(g_e_nl);  % normalized (unit) direction vector of each edge
%                     e_l = gather(g_e_l); %mm, current edge length
%                     e_l0t = gather(g_e_l0t ); %mm, current unstrained edge length, including thermal expansion
%                     
%                     % triangle vectors and arrays
%                     t_cx = gather(g_t_cx);
%                     t_cy = gather(g_t_cy);
%                     t_cz = gather(g_t_cz);
%                     t_t = gather(g_t_t);  %temperature
%                     t_broken = gather(g_t_broken);
%                     t_notbroken = gather(g_t_notbroken);
%                     t_optFcn_r2 = gather(g_t_optFcn_r2);  % intermediate variable for calculating gaussian illumination
%                     t_ev1 = gather(g_t_ev1); % edge vector for first triangle edge
%                     t_ev2 = gather(g_t_ev2);  % edge vector for second edge
%                     t_ev1n = gather(g_t_ev1n); % normalized (unit) vector for #1 triangle edge
%                     t_ev2n = gather(g_t_ev2n); % normalized (unit) vector for #2 triangle edge
%                     t_norms = gather(g_t_norms);   % triangle normal vectors, normalized
%                     t_IDir = gather(g_t_IDir);   %light incidence direction at each triangle centroid
%                     t_IMag = gather(g_t_IMag);  % magnitude of light intensity at each triangle centroid (watts per area)
%                     t_Ipwr = gather(g_t_Ipwr);  % incident light power on each triangle (watts)
%                     t_aproj = gather(g_t_aproj);  % triangle area in light incidence plane, i.e., projected area in xy plane
%                     t_specA = gather(g_t_specA);  % specular absorption
%                     t_specR = gather(g_t_specR);  % specular reflection
%                     t_abs = gather(g_t_abs);  %absorbed optical power in each triangle (watts)
%                     t_abs_mr = gather(g_t_abs_mr);  %absorbed optical power in each triangle (watts)
%                     t_reflp = gather(g_t_reflp);  %reflected optical power in each triangle (watts)
%                     t_refln = gather(g_t_refln); % reflected light direction (unit vector)
%                     t_oth = gather(g_t_oth);  % optical thrust vectors (force, newtons)
%                     t_oth_mr = gather(g_t_oth_mr); % optical thust vectors, multi-reflection (force, newtons)
%                     t_texn = gather(g_t_texn);  % unit vector within each triangle plane, indicating orientation axis of texture ("parallel" direction)
%                     t_textn = gather(g_t_textn);  % unit vector within each triangle plane, perpendicular to the texture orientation ("transverse" direction)
%                     t_press_n = gather(g_t_press_n);  %LUT photon pressure coefficient, normal to triangle
%                     t_press_p = gather(g_t_press_p);  %LUT photon pressure coefficient, parallel to triangle
%                     t_press_t = gather(g_t_press_t);  %LUT photon pressure coefficient, transverse to triangle
%                     t_lut_a = gather(g_t_lut_a);    %LUT absorption coefficient for each triangle
%                     
                    %anybroken = gather(anybroken);
                    
                    timers(10) = timers(10) + toc;
                    tic;
                    
                    %% PLOT VIDEO OUTPUT FRAME
                    
                    if (nt - pf >= frame_interval - 1)
                        pf = nt + 1;  % time step counter of prior frame
                        nf = nf + 1;  % number of movie frames
                        
                        %       [~,idx_max_theta] = max(abs(pitch_angles312));
                        %        [~,idx_max_phi] = max(abs(roll_angles312));
                        %       [~,idx_max_psi] = max(yaw_angles312);
                        
                        
                        stringtodisplay = ['n = ' num2str(nt) ',  nto = ' num2str(nto) ',  nf = ' num2str(nf) ',  tt = ' num2str(tt)  ...
                            ... %                  ',  max_theta = ' num2str(rad2deg(pitch_angles312(idx_max_theta)))  ',  max_phi = ' num2str(rad2deg(roll_angles312(idx_max_phi)))  ',  max_psi = ' num2str(rad2deg(yaw_angles312(idx_max_psi)))   ...
                            ',  x = ' num2str(comxs(nto)) ',  y = ' num2str(comys(nto)) ',  z = ' num2str(comzs(nto)) ...
                            ',  max_str = ' num2str(maxstrainsnorm(nto)) ',  maxtemp = ' num2str(maxts(nto)) ',  maxv = ' num2str((maxvrel)) ...
                            ',  uh-ohs: ' num2str(graphuos(nto)) ',  vz = ' num2str(velzs(nto)) ',  vx = ' num2str(velxs(nto)) ',  vy = ' num2str(velys(nto)) ...
                            ',  PE = ' num2str(potenergies(nto)) ',  KE = ' num2str(kinenergies(nto)) ...
                            ... %   '      tavg=' num2str(tempavg) ' tmax=' num2str(temppk) ' qin=' num2str(qin) ' qout=' num2str(qout) ...
                            ...  %   '      mincomp=' num2str(mincomp) ' maxcomp=' num2str(maxcomp) ' ET=' num2str(toc) 's' ]);
                            ',  CPS(sim) = ' num2str((nt-lastdispnt)/toc(lastdisptic)) ',  ET = ' num2str(toc(startsimtic)) 's' ...
                            ',  timeRatio = ' num2str((tt+spinup_time)/toc(startsimtic)) ];
                        
                        
                        % fprintf('%-4d  Node 5:  pos=%+.3e,%+.3e,%+.3e     vel=%+.3e,%+.3e,%+.3e      mf=%+.3e,%+.3e,%+.3e       of=%+.3e,%+.3e,%+.3e \n', nt, n_x(5), n_y(5), n_z(5), n_vx(5), n_vy(5), n_vz(5), n_mf(1,5), n_mf(2,5), n_mf(3,5), n_of(1,5), n_of(2,5), n_of(3,5) );
                        
                        if (videoMode<2)
                            if tt > 0 || ( tt < 0 && ~skip_spinup_video )
                                try
                                    if (gcf ~= fSimVid)
                                        figure(fSimVid);
                                    end
                                    plotMesh_v18_RG;
                                    open(V1);
                                    currentFrame = getframe(gcf);
%                                     size(currentFrame.cdata)
%                                     cdata = print(gcf, '-RGBImage');
                                    writeVideo(V1,currentFrame);
                                catch the_last_few_minutes_of_an_interesting_tv_show
                                    warning('Something went wrong with the video recorder!  No longer recording video output.');
                                    videoMode = 2;
                                    try
                                        close(V1);
                                        close(fSimVid);
                                    catch light_of_the_fact_that_youve_failed_once_again
                                        disp('Fuck it dude, let''s go bowling...');
                                    end
                                end
                            end
                        end
                        
                        disp( [ stringtodisplay ', CPS(ttl) = ' num2str((nt-lastdispnt)/toc(lastdisptic)) ] );
                        
                        lastdispnt = nt;
                        lastdisptic = tic;
                    end
            
                    if (videoMode < 2) && (graphuos(nto) > 0)
                        frame_interval = max(floor(breakup_slowmo_interval_s/dt),1); %slow down for breakup!
                    end
            
                    % keep track of time spent plotting video frames
                    timers(11) = timers(11) + toc;
                    tic;
            
                    %% 18. plot to snapshot window
                    if (tt < 0)
                        ps = nt;  % Do not start recording snapshots during spin-up 
                    end

                    if ((nt-ps >= snapshot_intvl-1 ) && (ns < 50) ) || (graphuos(nto) >= nextExplosionSnapshot) || ( (tt >= 0) && save_first_snapshot == 1)
                                                
                        ps=nt+1;
                        ns=ns+1;

                        tscale_max = max(tscale_max, maxts(nto));
                        tscale_min = min(tscale_min, maxts(nto));
                        if (graphuos(nto) >= nextExplosionSnapshot)
                            nextExplosionSnapshot = nextExplosionSnapshot * explosionSnapshotRatio;
                        end
                        figure(fSnapshot);
                        plotSnapshot;
                        
                        % Save workspace at each snapshot 
                        disp('Saving workspace');
                        allvars_snapshot = whos;
                        tosave_snapshot = cellfun(@isempty, regexp({allvars_snapshot.class}, 'matlab.(graphics|ui)'));
%                         tosave_snapshot = cellfun(@isempty, regexp({allvars_snapshot.class}, '^matlab\.(ui|graphics)\.'));
        
                        % Pass these variable names to save
                        save(fullfile(path_workspaces,[filebasename '_workspace',num2str(ns-1),'.mat']),  allvars_snapshot(tosave_snapshot).name, '-v7.3');
                        
                        % Setting to save snapshot at t = 0
                        if save_first_snapshot == 1
                            save_first_snapshot = 0;
                        end

                    end

                    % keep track of time spent plotting snapshots
                    timers(12) = timers(12) + toc;
                    tic;
                    
                end
                
            end
            
            %% 19. Check for user interaction with dialog boxes
            % Check at least once per 200 steps, two frame updates, or 5 seconds (whichever happens first)
            if (~mod(nt,min(200,2*frame_interval))) || (toc(lastbuttoncheck)>5)
                lastbuttoncheck = tic;
                usercancel = FS.Stop();
                % This is no longer used, because we no longer have a global rigid vs. flexible mode setting.  Rigid
                % bodes and flexible bodies are now allowed to coexist with the updated code.
%                 if MODEBOX.Stop()
%                     rigid_mode = ~rigid_mode;
%                     if rigid_mode
%                         MODEBOX = stoploop({'Currently using rigid dynamics.' 'Press OK to change to flexible.'}) ;
%                         % calcualte inertia tensor, which is otherwise skipped in flex mode
% %                         inertiaT = [  sum( n_m .* ( n_dy2 + n_dz2 ) )        ...
% %                             sum( -n_m .* n_dx .* n_dy)         ...
% %                             sum( -n_m .* n_dx .* n_dz );  ...
% %                             sum( -n_m .* n_dy .* n_dx )            ...
% %                             sum( n_m .* (n_dx2 + n_dz2) )      ...
% %                             sum( -n_m .* n_dy .* n_dz );  ...
% %                             sum( -n_m .* n_dz .* n_dx )            ...
% %                             sum( -n_m .* n_dz .* n_dy )        ...
% %                             sum( n_m .* (n_dx2 + n_dy2) ) ];
%                             g_inertiaT = [  ...
%                                 sum( g_n_m .* (g_n_dy2 + g_n_dz2) )        ...
%                                 sum( -g_n_m .* g_n_dx .* g_n_dy)         ...
%                                 sum( -g_n_m .* g_n_dx .* g_n_dz );  ...
%                                 sum( -g_n_m .* g_n_dy .* g_n_dx )            ...
%                                 sum( g_n_m .* (g_n_dx2 + g_n_dz2) )      ...
%                                 sum( -g_n_m .* g_n_dy .* g_n_dz );  ...
%                                 sum( -g_n_m .* g_n_dz .* g_n_dx )            ...
%                                 sum( -g_n_m .* g_n_dz .* g_n_dy )        ...
%                                 sum( g_n_m .* (g_n_dx2 + g_n_dy2) ) ]';
%                     else
%                         MODEBOX = stoploop({'Currently using flexible dynamics.' 'Press OK to change to rigid.'}) ;
%                     end
%                 end
           
            %% 20.  Check for termination conditions
            % record reason, save final snapshot, then break simulation loop

            xextent = (max(n_dr));
            xextentratio = xextent / (radiusmm);
            failureratio = sum(t_broken) / length(t_na);
            COMBeamDistNorm = (sqrt( (comx)^2 + (comy)^2))/(radiusmm);
            %nto>749 ||
            if( usercancel || ((max(n_vx)) > 1e12) || (xextentratio > explosion_ratio ) || ...
                    (COMBeamDistNorm > flyaway_ratio) || (failureratio > failure_ratio) || ...
                    ( terminate_t && (tt > terminate_t) )    )
                
%                 % gather things needed from gpu
%                     ringradii  = gather(g_ringradii);
%                     
%                     % node vectors and arrays
%                     n_x = gather(g_n_x);
%                     n_y = gather(g_n_y);
%                     n_z = gather(g_n_z);
%                     n_t = gather(g_n_t);
%                     n_a = gather(g_n_a);
%                     n_af = gather(g_n_af);
%                     n_dr = gather(g_n_dr);
%                     n_dx = gather(g_n_dx);
%                     n_dy = gather(g_n_dy);
%                     n_dz = gather(g_n_dz);
%                     n_dx2 = gather(g_n_dx2);
%                     n_dy2 = gather(g_n_dy2);
%                     n_dz2 = gather(g_n_dz2);
%                     n_dxyz = gather(g_n_dxyz);
%                     n_fx = gather(g_n_fx);
%                     n_fy = gather(g_n_fy);
%                     n_fz = gather(g_n_fz);
%                     n_mf = gather(g_n_mf);
%                     n_of = gather(g_n_of);
%                     n_vx = gather(g_n_vx);
%                     n_vy = gather(g_n_vy);
%                     n_vz = gather(g_n_vz);
%                     comx = gather(comx);
%                     comy = gather(comy);
%                     comz = gather(comz);
%                     n_abs = gather(g_n_abs);  % absorbed power (watts), from propulsion beam, heating input to each node
%                     n_ems = gather(g_n_ems);  % emitted power (watts), from radiative cooling, heating output from each node
%                     n_hf = gather(g_n_hf);  % heat flow (watts), from thermal conduction, heat into or out of each node
%                     
%                     % edge vectors and arrays
%                     e_t = gather(g_e_t);
%                     e_broken = gather(g_e_broken);
%                     e_notbroken = gather(g_e_notbroken);
%                     e_dl = gather(g_e_dl);  % mm,
%                     e_nl = gather(g_e_nl);  % normalized (unit) direction vector of each edge
%                     e_l = gather(g_e_l); %mm, current edge length
%                     e_l0t = gather(g_e_l0t ); %mm, current unstrained edge length, including thermal expansion
%                     
%                     % triangle vectors and arrays
%                     t_cx = gather(g_t_cx);
%                     t_cy = gather(g_t_cy);
%                     t_cz = gather(g_t_cz);
%                     t_t = gather(g_t_t);  %temperature
%                     t_broken = gather(g_t_broken);
%                     t_notbroken = gather(g_t_notbroken);
%                     t_optFcn_r2 = gather(g_t_optFcn_r2);  % intermediate variable for calculating gaussian illumination
%                     t_ev1 = gather(g_t_ev1); % edge vector for first triangle edge
%                     t_ev2 = gather(g_t_ev2);  % edge vector for second edge
%                     t_ev1n = gather(g_t_ev1n); % normalized (unit) vector for #1 triangle edge
%                     t_ev2n = gather(g_t_ev2n); % normalized (unit) vector for #2 triangle edge
%                     t_norms = gather(g_t_norms);   % triangle normal vectors, normalized
%                     t_IDir = gather(g_t_IDir);   %light incidence direction at each triangle centroid
%                     t_IMag = gather(g_t_IMag);  % magnitude of light intensity at each triangle centroid (watts per area)
%                     t_Ipwr = gather(g_t_Ipwr);  % incident light power on each triangle (watts)
%                     t_aproj = gather(g_t_aproj);  % triangle area in light incidence plane, i.e., projected area in xy plane
%                     t_specA = gather(g_t_specA);  % specular absorption
%                     t_specR = gather(g_t_specR);  % specular reflection
%                     t_abs = gather(g_t_abs);  %absorbed optical power in each triangle (watts)
%                     t_abs_mr = gather(g_t_abs_mr);  %absorbed optical power in each triangle (watts)
%                     t_reflp = gather(g_t_reflp);  %reflected optical power in each triangle (watts)
%                     t_refln = gather(g_t_refln); % reflected light direction (unit vector)
%                     t_oth = gather(g_t_oth);  % optical thrust vectors (force, newtons)
%                     t_oth_mr = gather(g_t_oth_mr); % optical thust vectors, multi-reflection (force, newtons)
%                     t_texn = gather(g_t_texn);  % unit vector within each triangle plane, indicating orientation axis of texture ("parallel" direction)
%                     t_textn = gather(g_t_textn);  % unit vector within each triangle plane, perpendicular to the texture orientation ("transverse" direction)
%                     t_press_n = gather(g_t_press_n);  %LUT photon pressure coefficient, normal to triangle
%                     t_press_p = gather(g_t_press_p);  %LUT photon pressure coefficient, parallel to triangle
%                     t_press_t = gather(g_t_press_t);  %LUT photon pressure coefficient, transverse to triangle
%                     t_lut_a = gather(g_t_lut_a);    %LUT absorption coefficient for each triangle                     

                if (max(n_vx) > 1e12)
                    end_reason = 'sim-diverged';
                elseif (xextentratio > explosion_ratio )
                    end_reason = 'exploded';
                elseif (failureratio > failure_ratio ) 
                    end_reason = 'substantial-failure';
                elseif (COMBeamDistNorm > flyaway_ratio)
                    if ~any(e_broken)
                        if (tt > (t_ramp_delay+2*t_ramp_dur+t_on_dur) )
                            end_reason = 'launched-off_course';
                        else
                            end_reason = 'flyaway';
                        end
                    else
                        end_reason = 'brokeaway';
                    end
                elseif ( terminate_t && (tt>terminate_t) )
                    end_reason = 'delay';
                    if any(upsidedowntris)
                        end_reason = [end_reason '_flip'];
                    end
                    if graphuos(nto)
                        end_reason = [end_reason '_fail'];
                    end
                else
                    end_reason = 'cancelled';
                end
                if (isfinite(snapshot_intvl))
                    ns=ns+1;
                    ps=ps+1;
                    figure(fSnapshot);
                    plotSnapshot;
                end
                timers(12) = timers(12) + toc;
                tic;
                break;
            end
            timers(12) = timers(12) + toc;
            tic;
            
            end
  
            nt=nt+1; %completion of time-domain loop (originally was a for loop, now while with end increment.
            
        end
        
        %% Display simulation timing statistics
        fprintf('SIMULATION ENDED:  ');
        disp(end_reason);
        allvars = whos;
        disp(['Workspace memory usage: ' num2str(sum([allvars(:).bytes])./1024./1024) ' MB']);
        disp(['Elapsed time: ' num2str(toc(startlooptic)) ]);
        disp(['Sim time: ' num2str(toc(startsimtic)) ' (time ratio ' num2str((tt+spinup_time)/toc(startsimtic)) ')']);
        ttltimers = sum(timers);
        et=toc(startsimtic);
        rate_oa = nt/et;
        sim_time = sum(timers([1:10 13 14 15]));
        render_time = sum(timers([ 11 12] ));
        raytrace_time = sum(timers(14));
        rate_loop = nt/sim_time;
        rate_rend = (nf+ns)/render_time;
        rate_rt = max((nt - nt_startAccel),0) / raytrace_time;
        disp(['Timer results:    Init time   ' sprintf('%10.3f',init_time)        ' sec.     ']);
        init_time=0;  % clear for successive looped sims
        disp(['                  Setup time      ' sprintf('%10.3f',setup_time)       ' sec.     ']);
        disp(['                  Sim time        ' sprintf('%10.3f',et)               ' sec.     Rate = ' num2str(rate_oa) ' loops/sec']);
        disp(['                  MainLoop time   ' sprintf('%10.3f',sim_time)         ' sec.     Rate = ' num2str(rate_loop) ' loops/sec']);
        disp(['                  Render time     ' sprintf('%10.3f',render_time)      ' sec.     Rate = ' num2str(rate_rend) ' frames/sec']);
        disp(['                  Raytrace time   ' sprintf('%10.3f',raytrace_time)    ' sec.     Rate = ' num2str(rate_rt) ' frames/sec']);
        fprintf('Timer:   ');
        %         for ii=1:length(timerBins)    fprintf('%-12s\t', timerBins{ii}); end
        %         fprintf(newline);
        %         fprintf('Percent: ');
        %         for ii=1:length(timerBins)    fprintf('%% %-10g\t', timers(ii)./ttltimers*100 ); end
        %         fprintf(newline);
        %         fprintf('Seconds: ');
        %         for ii=1:length(timerBins)    fprintf('s %-10g\t', timers(ii)); end
        %         fprintf(newline);
        for ii=1:length(timerBins)    fprintf('%-11s\t', timerBins{ii}); end
        fprintf(newline);
        fprintf('Percent: ');
        for ii=1:length(timerBins)    fprintf('%-11.2f\t', timers(ii)./ttltimers*100 ); end
        fprintf(newline);
       
        %% Complete video recording, close video window
        if (videoMode < 2)
            try
                close(V1);
            catch a_spinning_frisbee
                % Nope, fumbled it...
            end
            try
                close(fSimVid);
            catch covid_nineteen
                % shit, better drink some bleach.
            end
        end
        FS.Clear();
        %MODEBOX.Clear();
        
        %% Follow-up calculations and data handling
        nto_firstBroken = min(nto_firstBroken,nto+1);
        
        longaxis_length = vecnorm(hoop12xyz-hoop6xyz);
        shortaxis_length = vecnorm(hoop3xyz-hoop9xyz);
        
        nt=nto;   % change meaning of nt, now it means the number of monitor frames.
        plot_time = graphts(1:nt);
        
        
        if fft_mode
            if (monitor_depth > nto)
                monitorrfs = monitorrfs(:,1:nto);
                monitortfs = monitortfs(:,1:nto);
                %monitorfxs = monitorfxs(:,1:nto);
                %monitorfys = monitorfys(:,1:nto);
                monitorfzs = monitorfzs(:,1:nto);  %trim these down now...
                monitorxs = monitorxs(:,1:nto);
                monitorys = monitorys(:,1:nto);
                monitorzs = monitorzs(:,1:nto);
            else
                warning('Monitor buffer overflow/wrap!!! fix me!!!');
            end
        end
        
        %% Prepare / save main simulation output window
        disp('Generating output plot.');
%         hOutputPlot = figure('pos',[100 100 2100 1200]);
        hOutputPlot = figure('pos',[50 100 1800 800]);
        set(gcf,'color','w');
        
        hafp1 = subplot(4,5,1);
        plot(plot_time, velxs(1:nt));
        title('X & Y velocity (mm/s)');
        xlabel('Time');
        hold on
        plot(plot_time, velys(1:nt));
        legend({'X' 'Y'});
        
        
        hafp2 = subplot(4,5,2);
        plot(plot_time, velzs(1:nt));
        title('Z velocity (mm/s)');
        
        hafp3 = subplot(4,5,3);
        plot(plot_time, areas(1:nt));
        title('Area (mm^2)');
        hold on
        plot(plot_time, areasxy(1:nt));
        plot(plot_time, areasrightsideupxy(1:nt), '--');
        Ylims = get(gca,'YLim');
        set(gca,'Ylim',[ Ylims(1) Ylims(2)+0.3*(Ylims(2)-Ylims(1)) ] );
        legend({'Surface area', 'Projected area', 'Upright projected' } );
        
        
        hafp4 = subplot(4,5,4);
        plot(plot_time, minstrains(1:nt));
        title('Strain min/max');
        hold on;
        plot(plot_time, maxstrains(1:nt));
        plot(plot_time([1 end]), [0 0], 'k');
        
        hafp5 = subplot(4,5,5);
        plot(plot_time, graphuos(1:nt));
        title('Number of tensile failures');
        set(gca,'YScale','log')
        
        hafp6 = subplot(4,5,6 );
        plot(plot_time, comxs(1:nt)-((0)));
        title({'Position (mm)' '(COM rel to beam ctr)'});
        hold on;
        plot(plot_time, comys(1:nt)-0);
        plot(plot_time, comrs(1:nt), '--');
        legend({'X' 'Y' 'R'});
        
        
        hafp7 = subplot(4,5,7);
        plot(plot_time, comzs(1:nt));
        title('Z position (mm)');
        
        hafp8 = subplot(4,5,8);
        plot(plot_time, kinoenergies(1:nt));
        title('Non-ballistic KE (J)');
        
        hafp9 = subplot(4,5,9);
        plot(plot_time, potenergies(1:nt));
        title('Potential energy (J)');
        
        
        %subplot(3,4,9);
        %plot(plot_time, maxvs(1:nt));
        %title('Max velocity rel. to COM');
        
        hafp10 = subplot(4,5,10);
        plot(plot_time, rad2deg( asin(mintrinormz(1:nt))));
        title('Steepest triangle angle (\circ)');
        hold on
        plot(plot_time([1 end]), [0 0], 'k');
        
        
        
        hafp11 = subplot(4,5,11);
        plot(plot_time, ttlabspwr(1:nt));
        title('Power absorbed / emitted');
        hold on;
        plot(plot_time, ttlradpwr(1:nt));
        legend({ 'Abs', 'Ems' });
        
        hafp12 = subplot(4,5,12);
        plot(plot_time, acczs(1:nt));
        hold on;
        plot(plot_time, accts(1:nt));
        title('Acceleration (mm/s/s)');
        legend({'Z' 'Total'});
        
        hafp13 = subplot(4,5,13);
        plot(plot_time, avgts(1:nt));
        hold on
        plot(plot_time, vertextemps(1:nt));
        plot(plot_time, xradpttemps(1:nt));
        plot(plot_time, maxts(1:nt));
        title('Temperature (\circ K)');
        legend({ 'Avg' 'Vertex' 'Edge point' 'Peak' } );
        
        hafp14 = subplot(4,5,14);
        plot(plot_time(2:end), diff(avgts(1:nt))./dnto);
        hold on
        %plot(plot_time(2:end), diff(maxts(1:nt))./dnto);  % this one isn't really valid as the point index shifts all the time...
        plot(plot_time(2:end), diff(vertextemps(1:nt))./dnto);
        plot(plot_time(2:end), diff(xradpttemps(1:nt))./dnto);
        title('dT/dt (\circK/s)');
        plot(plot_time([1 end]), [0 0], 'k');
        legend({ 'Avg' 'Vertex' 'Edge point' } );
        
        hafp15 = subplot(4,5,15);
        plot(graphts(1:nto_firstBroken-1), longaxis_length(1:nto_firstBroken-1));
        hold on;
        title({'Elongation' '[shape base diameter (mm)]'});
        plot(graphts(1:nto_firstBroken-1), shortaxis_length(1:nto_firstBroken-1));
        plot(plot_time([1 end]), 2*[radActualx0mm radActualx0mm], 'k--');
        legend({ 'Major dia' 'Minor dia' 'Rest dia'} );
        
        hafp16 = subplot(4,5,16);
        plot(plot_time, rad2deg(sailPitches(1:nt)) );
        hold on;
        plot(plot_time, rad2deg(atan2( sailNorms(2,1:nt), sailNorms(3,nt) )) )
        plot(plot_time, -rad2deg(atan2( sailNorms(1,1:nt), sailNorms(3,nt) )) )
        title('Sail tilt (\circ)' );
        legend({ 'total off Z'  'X' 'Y'});
        
        hafp17 = subplot(4,5,17);
        plot(plot_time, abs(angVels(3,1:nt)./(2*pi)) );
        hold on;
        plot(plot_time, vecnorm(angVels(:,1:nt))./(2*pi) );
        title('Angular velocity (Hz)' );
        legend({'Along Z (abs)' 'Total'});
        
        hafp18 = subplot(4,5,18);
        plot(plot_time, angVels(1,1:nt)./(2*pi) );
        hold on;
        plot(plot_time, angVels(2,1:nt)./(2*pi) );
        title('X and Y angular velocity (Hz)' );
        legend({ 'X' 'Y' } );
        
        hafp19 = subplot(4,5,19);
        plot(plot_time, ofxs(1:nt) );
        hold on;
        plot(plot_time, ofys(1:nt) );
        title('Lateral optical force' );
        legend({ 'X' 'Y' } );
        
        hafp20 = subplot(4,5,20);
        plot(plot_time, tqxs(1:nt) );
        hold on;
        plot(plot_time, tqys(1:nt) );
        plot(plot_time, tqzs(1:nt) );
        title('Effective torque on sail' );
        legend({ 'X' 'Y' 'Z'} );
        
        
        linkaxes([hafp1, hafp2, hafp3, hafp4, hafp5, hafp6, hafp7, hafp8, hafp9, hafp10, hafp11, hafp12, hafp13, hafp14, hafp15, hafp16, hafp17, hafp18, hafp19, hafp20 ], 'x');
        
        annotation('textbox', [.02 .95 .95 .04], 'String', ...
            sprintf(['Trajectory for %s       \n' ...
            'pD=%g\t      radius=%g mm\t     rot0=%d Hz\t   ' ...
            'I0=%d   \t     x0=%g mm\t     y0=%g mm\t     beamRadius=%g mm\n' ], ...
            filebasename,  pD, radiusmm, rot0, I0, t0t_x0, t0t_y0, Irad), ...
            'EdgeColor', 'none', 'Interpreter', 'none' );
        
        saveas(gcf,fullfile(path_figures,[filebasename '.fig']));
        saveas(gcf,fullfile(path_figures,[filebasename '.png']));
        
        %% Save snapshot figure
        disp('Generating snapshot plot.');
        figure(fSnapshot);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        set(gca,'ZTick',[]);
        set(gca,'Box','on');
        title(filebasename, 'Interpreter', 'none');
        plot3([((0)) ((0))], [((0)) ((0))], get(gca,'ZLim'),'m');
        
        saveas(gcf,fullfile(path_figures,[filebasename '_snapshot.fig']));
        saveas(gcf,fullfile(path_figures,[filebasename '_snapshot.png']));
        
        %% Generate / save trajectory figure
        disp('Generating trajectory plot.');
        if multiloop
            if nloop == 1
                hTrajectoryFig = figure;
            else
                try
                    figure(hTrajectoryFig);
                catch dontcare
                    hTrajectoryFig = figure;
                end
            end
        else
            hTrajectoryFig = figure;
        end
        plotTrajectory;
        hold on;
        saveas(gcf,fullfile(path_figures,[filebasename '_trajxy.fig']));
        saveas(gcf,fullfile(path_figures,[filebasename '_trajxy.png']));
       
        %% Generate / save restoring force figure
        if nto_startAccel && (nto_firstBroken > nto_startAccel)
            disp('Generating restoring force analysis plot.');
            plotRestoringForceAnalysis;
            saveas(gcf,fullfile(path_figures,[filebasename '_forcetorque.fig']));
            saveas(gcf,fullfile(path_figures,[filebasename '_forcetorque.png']));
        end
        
        %% Generate/display data table entry
        %        disp('Data table entry:')
        
        peakstrain = max(maxstrains);
        peakt = max(maxts);
        peakaccz = max(acczs);
        peakpwr = max(ttlinpwr);
        
        spacer = 0;  %just used to fill in missing columns in the table
%         dispvars = { 'movieNum',...
%             'radiusmm','zAspectRatio','totalarea','totalmass','sum(t_a0xy)','rot0', ...
%             'numnodes','numedges','numtris',...
%             'thickness','density','Ymod',...
%             'Emod','Thermcondmm','heatCap',... %   'cte0','cte1','cts','ctw',  ...
%             't0t_x0','t0t_y0','t0t_s0','t0t_t0_x','t0t_t0_y','t0t_t0_z', ...
%             'I0','Irad','spacer', ...
% %             'Iabs','Irefl','Emissivity', ...
%             'spinup_gs','spinup_time','videoMode', ...
%             't_ramp_delay','t_ramp_dur','t_on_dur', ...
%             'dt','tt','nt','frame_interval','snapshot_intvl', ...
%             'uhohs','comx','comy','comz', ...
%             'meanvx','meanvy','meanvz', ...
%             'peakt','peakaccz','peakpwr','peakstrain', ...
%             'rate_oa', 'rate_loop', 'rate_rend', ...
%             'spacer','time_frstBrkn','et' };
        dispvars = { 'movieNum',...
            'radiusmm','zAspectRatio','totalarea','totalmass','sum(t_a0xy)','rot0', ...
            'numnodes','numedges','numtris',...
            'thickness','density','Ymod',...
            'Emod','Thermcondmm','heatCap',... %   'cte0','cte1','cts','ctw',  ...
            't0t_x0','t0t_y0','t0t_s0','t0t_t0_x','t0t_t0_y','t0t_t0_z', ...
            'I0','Irad','spacer', ...
            'spinup_gs','spinup_time','videoMode', ...
            't_ramp_delay','t_ramp_dur','t_on_dur', ...
            'dt','tt','nt','frame_interval','snapshot_intvl', ...
            'uhohs','comx','comy','comz', ...
            'meanvx','meanvy','meanvz', ...
            'peakt','peakaccz','peakpwr','peakstrain', ...
            'rate_oa', 'rate_loop', 'rate_rend', ...
            'spacer','time_frstBrkn','et' };
        dispstrs = {'illum_mode', 'zMode', 'xyMode', 'curvatureMode', 'end_reason'}; 
        
        %      disptablerow( dispvars, dispstrs );  %this is fucking obnoxious
        
        if nouterloop1==1 && nouterloop2==1
            outerloopstable = getdisptableheader( dispvars, dispstrs);
        end
        
        outerloopstable = [ outerloopstable newline ...
            getdisptablerow( dispvars, dispstrs ) ];
        
        burstTestVals(end+1) = pwrRampVal * I0
        
        %% Save the whole workspace (except figures)
        disp('Saving workspace');
        allvars = whos;
        % Identify the variables that ARE NOT graphics handles. This uses a regular
        % expression on the class of each variable to check if it's a graphics object
%         tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
        tosave = cellfun(@isempty, regexp({allvars.class}, 'matlab.(graphics|ui)'));
        % Pass these variable names to save
        save(fullfile(path_workspaces,[filebasename '_workspace.mat']),  allvars(tosave).name, '-v7.3');
        %  Note:  I don't need to worry about freeing unused buffer space before saving the workspace, because somehow
        %  MATLAB compresses the workspace when saving to disk.
        
        %% Wrap up, prepare for next simulation
        movieNum = movieNum+1;
        if ~isempty(simdesc)
            filebasename = sprintf('sim%d_%s', movieNum, simdesc);
        else
            filebasename = sprintf('sim%d', movieNum );
        end
        
        % End of outerloop2
        if (usercancel)
            break;
        end
        
        % close figures if perfomring multiple sims
        if multiloop
            try
                close(hOutputPlot);
                close(fSnapshot);
                close(hRestoringForce);
                close(hMeshPreview);
            catch nothing
            end
        end
        
    end
    
    % End of outerloop1
    if (usercancel)
        break;
    end
    
    if length(outersweepvals1)>1
        try
            close(hTrajectoryFig)
        catch nothing
        end
    end
end

if multiloop && ~usercancel
    fprintf('\n\n\n');
    disp(['Completed simulation script, elapsed time = ' num2str(toc(startscripttic)) ' sec.']);
    fprintf('\n\n');
    disp(outerloopstable);
end

diary off;
