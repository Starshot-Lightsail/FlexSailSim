%% STARSHOT LIGHTSAIL SIMULATOR
% (c) 2018-2023 Michael Kelzenberg, Ramon Gao -- California Institute of Technology
% Largely well commented, with external documentation now in development.  Ad Astra!

disp('Flexible lightsail simulator ');
disp('(c) 2018-2023 Michael D. Kelzenberg and Ramon Gao, California Institute of Technology');
disp('Version: 2023-10-02 (v20)');

%% Initialization block:
%clear  % clear the workspace to make sure the code will run in a fresh workspace, and prevent carry-over of prior
        % simulation data.  We often comment this out during development though, to preserve analysis/support variables.
maxNumCompThreads(8);
close all hidden  % close all windows from prior simulations or whatever.  The 'hidden' part closes message boxes too.

simdesc = ''; % Use this if you want to add a short description of why you're running the simulation.  
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

% Here we save a copy of the exact simulation code used for each simulation session, so we can go back and deduce
% the specific settings and algorithms that were used for any recorded simulation
copyfile([ mfilename '.m'], [filebasename '.script']);

% We also save the terminal output so we can see what was displayed during the simulation.
diary([filebasename '.diary']); % Log command window text to file

%% Set up material properties  (external script)
% Defines basic material properties and physical constants.  As of this writing, supported materials are 'silicon' and 
% 'nitride.'  For various reasons, we define certain lightsail specific parameters here as well, such as thickness.  If
% you want to change these values iteratively during a batch of simulations, make sure that such changes are applied
% prior to mesh generation.
setupMaterialProperties_v20;  


%% Select whether to use MEPS LUTs, and pre-load the data if needed (external script)
% MK (April 6-8 2021):  This is a new approach to loading and accessing the MEPS data from Ramon's simulations.  This
% approach can also be used to load any other optical properties of the lightsail regions, such as from theory or from
% measurements.  See the load function called below for comments!
%
% Since we don't generally change the source of the optical data between iterative simulations, we load it here at the
% beginning of the script to save time and prevent redundant outputs.  In the future, we could increase the dimensionality
% of the lookup tables, for example, to enable properties to vary as a function of temperature as well.  
%
% IMPORTANT:  As of this writing, All LUTs must use the same values for LUTPsiStepDeg and LUTThetaStepDeg!
%  The values retained in LUT(n).pstep and .tstep are not used elsewhere in the simulation.
%  Using the same LUT sizes reduces the complexity of indexing the tables during the time-domain simulations.  
%
% ALSO:  Note that LUTs only affect the optical response of triangles for which you've assigned a non-zero texture value!
%
% FINALLY:  useMepsLUTs can't be used in conjunction with useRayTraceMode

useMepsLUTs = 0;  % if enabled, optical response for non-zero texture regions is calcualted from 2D LUTs, replacing whatever
                  % response was calculated based on specular (constant or TMM) model.
                  % IF DISABLED:  Doesn't use MEPS LUTS even if you specify textures for mesh regions

if useMepsLUTs
    LUTPsiStepDeg = 0.20;    %  These specify the angular resolution to use for the LUTs.  The values might be overwritten 
    LUTThetaStepDeg = 0.20;  %  in setupLUTs depending on the code configuration.  But you need to set a positive nonzero 
                             %  value for these two settings, even if not using LUTs, otherwise the simulation will error
                             %  out.
    setupLUTs_v17;  % now moved to subscript.  

    % Note:  If you want to change LUTs between consecutive simulations, you can modify the setupLUTs script to change its 
    % behavior based on 'outerloop' variables, then re-invoke it within each loop, e.g., near generateMesh.  Or I often find 
    % it easier to load all the LUTs for all designs here at the beginning, then assign whichever LUTs are actually used for 
    % each structure during mapping in generateMesh.  
else
    LUTs = [];
end
                         



%% Select whether to use TMM specular reflectance model, and pre-calculate optical response if used
% Unlike complex nanophotonic grating structures (MEPS), the optical resonse of simple thin-films membranes, and even 
% 1D stacks such as bragg reflectors, can be calculated as a function of incidence angle via closed-form expressions.
% Despite this, in practice, it's almost always faster to use a look-up approach in this case as well, in which we 
% obtain approximate values for reflectance and absorption by retrieving pre-calculated values from memory rather than 
% calculating them on the spot.  For specular surfaces, this can be accomplished by a 1D look-up-table.  For the sake of
% clarity, we call these "look-up lists" (LULs) to differentiate them from the 2D LUTs used for MEPS.
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
%                        Specular behavior determined by:  Lookup list structure calculated from material ior_n and ior_k.
%
% IMPORTANT:  the specular TMM LUL approach is only used for regions with a non-zero material assignment!
%   The default (0) material allows for fixed reflectance and fixed absorption only.
useSpecularTMM = 0;
if useSpecularTMM
    LULAngleStepDeg = 0.5; % angular resolution of LUL
    LULAngleStep = deg2rad(LULAngleStepDeg);
    setupSpecularLULs_v20;
else
    LULs = [];
    LUL = [];
end



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
RayTracingMode = 0;   % whether or not to raytracing stuff   %% moved to below (inside sweep loop)
MAX_REFL = 4;     % maximum number of ray reflections

%% Reverse raytrace mode, to see if a paraboloid sail could be used as a reflector for comms?
ReverseRayTrace = 0;  % whether or not reverse raytracing is used.  Only produces valid results for a tethered payload!


%% Radiative heat transfer
do_rht = 0; %flag for enabling RHT during simulations

%% LATERAL STABILITY VISUALIZATION
% I find this helpful when pondering the stability of spinning sails.  These precess instead of tipping in response
% to applied torque. Non-spinning sails tend to tip themselves back towards the beam in the process of stabilizing...
% but the spinning sails stay misaligned.  So it is important that they have a lateral equilibrium point over a 
% wide range of tip angles.  Or at least that's how I think of it...
%
% Anyway, this is a special visualization tool that runs right before the simulator starts.  It takes the whole sail
% and moves it around in X and Y, over a range of tip angles, to look at the lateral forces... and eventually prepare
% a plot.  Then, it puts the sail back where it found it and tiptoes out, so the simulation can begin normally.
doStabilityVisualization = 0;
 % 0:  Do not do the stability visualization thing
 % 1:  Do the visualization, then proceed with the simulation
 % 2:  Do the visualization only, then cancel the simulation before it starts


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
outersweepvals2 = 500;%1760 ; %fliplr(50 * 10.^[1.6:0.025:2]);  % currently used for:  nothing
outersweepvals1 = ([0 ]); 
burstTestVals = [];

outerloopstable = ''; %this initializes a string array to store the output table for all subsequent simulations
multiloop = (length(outersweepvals1) + length(outersweepvals2)) > 2;
nloop = 0;

for nouterloop1 = 1:length(outersweepvals1)
    for nouterloop2 = 1:length(outersweepvals2)
        startlooptic = tic;
        nloop = nloop + 1;
        outersweepval1 = outersweepvals1(nouterloop1);
        outersweepval2 = outersweepvals2(nouterloop2);
        fprintf('\n');
        if (length(outersweepvals2) > 1) || (length(outersweepvals1) > 1)
            disp(['OuterLoop [' num2str(nouterloop2) ',' num2str(nouterloop1) '] '  num2str(outersweepvals2(nouterloop2)) ' ' num2str(outersweepvals1(nouterloop1)) ] );
        end
        
       % RayTracingMode = 0;

        %% Generate mesh
        generateMesh_v20; % This helper script sets up the mesh.  Modify it to change mesh settings such as shape, patterning, and material properties.
        saveas(gcf,[filebasename '_mesh.fig']);  % save a rendering of the initial mesh, which is created as a
        saveas(gcf,[filebasename '_mesh.png']);  % figure by the generateMesh script
        
        %% Rigid/flexible mechanics:
        %  This is no longer used, we can define any number of rigid bodies in the mesher.  Do not explicitly set a
        %  value for rigid_mode!
        %rigid_mode = 0;  % set to 1 to use rigid body mechanics, or 0 for flexible mechanics.  This can be changed dymaically
        % throughout the simulation, but must be specified here for initial mesh generation.
        
        %% Illumination mode and power
        %  As of this writing, only 'gaussian' illumination profiles are supported.  We plan to add arbitrary profiles via LUTs
        %  in the future.  Beam is always centered at (0,0).  Starting sail offsets have now been moved to a different
        %  section.
        illum_mode = 'gaussian';  %Illumination mode selector.  Not yet implemented in the code (FIXME).
        %DEPRECIATED:  Ictrx = -radiusmm * 0.1;%outersweepvals1(nouterloop1);  % mm, beam offset in x direction, relative to origin of initial mesh.
        %DEPRECIATED:  Ictry = -radiusmm * 0.1;%outersweepvals1(nouterloop1);  % mm, beam offset in y direction, relative to origin of initial mesh.
        I0 = 5000; % W/mm^2, peak incident power density (1000 = 1GW/m2).
        
        %% Illumination settings for gaussian beam profile.
        Irad = radiusmm*sqrt(2); % * 0.5 * sqrt(2); %radiusmm * outersweepvals2(nouterloop2);  % mm, beam waist radius, for gaussian illumination
        
        %% Illumination settings for other beam profile modes.  To be added in the future.
        if ~isequal(illum_mode,'gaussian')
            error(['Attempted to use non-defined illumination mode ''' illum_mode '''']);
        end
        
                
        %% Spin-stabilization settings
        rps_target = 40;%outersweepval2; % Rotations per second for spin stabilization.  Can be negative to spin the other way.
        spinup_gs = max(abs(rps_target/100),2);  % How quickly to spin-up the lightsail.  As a time-domain simulator, the lightsail must begin at rest.
        %  We spin it up to a desired rotation speed during a brief time prior to the propulsion phase.  I'm not sure what
        %  the units are for this parameter, but it works out that we get 1808.8 rotations per second per simulated second of
        %  spin-up for a value of spinup_gs = 1.  Depending on material properties, too much angular acceleration can 
        %  result in unintended pressure waves ("shock") developing in the sail during spin-up.  Too rapid a spin-up 
        %  will also reduce the number of time steps over which spin-up occurs, leading to slight slight errors in the
        %  achieved rotational speed vs. rps_target, due to the discrete stepping.  On the other hand, spinning up too 
        %  slowly wastes time.    
        %NOTE:  Don't make spinup_gs negative!  Make rps_target negative to spin backwards. 
        
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
        t_ramp_delay = 0.00000; % (seconds), time delay before starting beam turn-on, formerly t_start_ramp
        t_ramp_dur = 0.001;%20;%3e-3; % (seconds), time duration of beam turn-on. DO NOT SET TO ZERO!  (Use a small number e.g., 1e-8 for instant turn-on.)
        t_on_dur = 1;%00e-3; % (seconds) duration of acceleration phase, after which beam is turned off
        t_coast_dur = 0.0; % (seconds) duration of time after beam turn-off for observing the final sail trajectory.
        
        dt = rec_dt; %seconds, the time step for the simulation.  rec_dt is the "recommended" calculated value from the mesh generator.  I do not recommend changing this.
        if rps_target  %reduce time step further if needed for rapidly spinning structures.  
            if dt > (1/rps_target/360)
                dt = (1/rps_target/360);
                warning(['Reducing time step from mesh-calculated value due to high spin speed.  dt = ' num2str(dt)]);
            end
        end

        %This is the calculated value for the total simulation duration after t=0; do not change:
        sim_dur = t_ramp_delay + 2*t_ramp_dur + t_on_dur+t_coast_dur; % seconds; note this doesn't include the t<0 spinup time.
        
        %Here is the actual power ramp function.  Do not change.
        pwrRampFcn = @(tt) 0.5*double(tt > t_ramp_delay).*double(tt <= t_ramp_delay+t_ramp_dur).*(1-cos(pi*(tt-t_ramp_delay)./(t_ramp_dur))) + double(tt > t_ramp_delay+t_ramp_dur).*double(tt <= t_ramp_delay+t_ramp_dur+t_on_dur) + ...
            0.5*double(tt > t_ramp_delay+t_ramp_dur+t_on_dur).*double(tt <= t_ramp_delay+t_ramp_dur+t_on_dur+t_ramp_dur).*(1+cos(pi*(tt-t_ramp_delay-t_ramp_dur-t_on_dur)./t_ramp_dur));
        
        %% Initial sail orientation and offsets  (t=0 transform)
        % The initial mesh and the beam are both always centered at (0,0).  So, to simulate a beam offset, tilt offset,
        % or a specific starting angular orientation for the sail, we need to translate or rotate the mesh after
        % spin-up, prior to beginning the acceleration.  (Spin-up forces are centered around (0,0) as well.)
        %
        % This is done as an instantanious rigid transformation at t=0, in violation of physical continuity, but 
        % ultimately this appraoch made the most sense.  
        t0t_changeSpinAngle = 0; % whether or not to change the rotational position of the sail, useful for spinning sails, so it always has the same rotational orientation when the light first turns on.
        t0t_s0 = psi0;  %the desired rotational position, if enabled.  
                                      %Typically we use psi0 from the mesh generator -- the initial angle of the mesh.
                                      %For a different rotational position, use desiredValue-psi0
        t0t_x0 = 0.1*radiusmm;
        t0t_y0 = 0.1*radiusmm;
        t0t_t0_x = 0; % (rad) initial tilt of sail (misalignment) about x-axis, performed at t=0
        t0t_t0_y = 0; % (rad) initial tilt of sail (misalignment) about y-axis, performed at t=0
        t0t_t0_z = 0; % (rad) additional t=0 rotation about z-axis, I'm not sure why we would ever do this, but it looks like the code supports it.
        
        
        %% Simulation video recording mode and rendering options:
        % In general it is nice to plot and record as much data as possible, but that can slow down the simulations.  Use
        % plotFast=1 to render a simplified plot of the lightsail for the video recorder, or plotFast=2 to disable plotting and
        % video recording alltogether.
        videoMode = 1;  %0 = full plots; 1 = simplified plot (slightly faster movie rendering); 2 = no plots (no video; fastest option)
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
        
        % This parameter lets us skip re-calculation of optical forces during most time steps.  The optical forces on
        % each triangle tend to change fairly slowly compared to the mechanical forces.  We can make simulations run
        % faster, with little loss of accuracy, by reducing the frequency of optical force calculation.  However, this
        % does fundamentally reduce the accuracy of simulations.  Value of 2 to 16 seem to provide reasonalbe speedup
        % for our simulations, without substantially affecting the results.  This affects all calculations of optical 
        % forces, as well as radiative heat transfer, raytracing, absorption heating, etc.  It does not affect mechanics,
        % thermal conduction, spinup, etc.
        %
        % WARNING:  Don't use huge opt_skip values for spin-stabilized sails.  The stale optical forces will apply a torque
        %           that slows down the rotation!
        opt_skip_downsample=4;  % number of time steps between optical force calculation.  Do not set to zero!
                                % Set to 1 for maximal accuracy of optiacal forces.  
        
        % monitor_mode:  if monitor_mode is set (1), we will record a handful of state variables at each iteration of the
        % loop, with no downsampling, so we can perform FFT analysis later.  This doesn't affect the key state
        % variables, which are still stored at downsampled intervals.  
        fft_mode = 0;  % set to 1 to record node position and forces at all times, then perform FFT analysis
        fft_start_time = 0.5; %seconds
        fft_dur = 0.5; %seconds

        % post-video mode:  this stores full node spatial information 
        save_mesh_mode = 1; % stores certain full-mesh information so we can post-render videos from saved data...
        save_mesh_dt = 0.0001; % seconds, time step between mesh frames
        
        % Our code also produces a video recording of the lightsail's behavior, rendered via a MATLAB figure window.  Rendering
        % each video frame is computationally intensive, typically requiring 100-10,000x longer than the calculation time for each
        % simulation time step.  Furthermore the visual changes to the mesh between consecutive time steps are impercievably
        % small.  Thus we can render video frames at a substantially reduced frequency, specified here by the value
        % frame_interval_s (in seconds).  Typical values of interest range from ~1e-5 seconds for recording extreme acceleration
        % events or catastrophic shape failures, to 1e-3 seconds to illustrate the trajectory and steady-state behavior of 
        % stabilized lightsail designs.
        frame_interval_s = 0.5e-3;     %seconds
        
        % frame rate reduction for spinning sails.  If you want the rotation of the sail to be clearly percievable in
        % the video recording, uncomment this line to make sure the frame rate is slow enough to capture it.  
        % if (rps_target)  frame_interval_s = min(frame_interval_s, 1/rps_target/36); end % make sure there are at least 36 frames per revolution
                
        % Skip rendering the spin-up sequence?  This prevents video rendering during spin-up, which can save time.
        skip_spinup_video = 0;  
        
        % When simulating unstable lightsail designs resulting in structural failure, it is sometimes of interest to visulize
        % their failure with greater time resolution, to help understand the cause.  This variable can be used to reduce the
        % video rendering time step upon the onset of structural failure:
        breakup_slowmo_interval_s = max( frame_interval_s/4, 5e-6);   %seconds; increases frame rate 10x for failures
        %breakup_slowmo_interval_s = frame_interval_s;  % use this value to disable slowmo for lightsail structure failure
        
        %% Snapshot display settings
        % Our simulator code also produces a composite 'snapshot' rendering of the evolution of the lightsail shape and
        % (transverse) position during the simulation, overlaying snapshots of the sail at fixed time intervals upon a single
        % axis, with each consecutive snapshot offset vertically above the prior one, accompanied by a timestamp annotation.
        % Because only a handful of shape renderings can be effectively displayed in a single figure, it is difficult to set
        % this value correctly in advance.  At a minimum, the sail's initial state and final state are rendered in this figure,
        % and with judicious choice of the snapshot interval, the evoloution from the former to the latter can be effectively
        % conveyed within a single static figure.  Good luck.
        snapshot_interval_s = 0.008;%0.020;%.250;
        
        % Our simulator detects the structural failure of a lightsail as the moment when the stress/strain in any single mesh
        % edge exceeds the tensile strength of the material, after which the failed edge is omitted (mechanically) in future
        % shape calculations.  It is useful to visualize the shape of the sail at and shorly after the onset of such failures.
        % The following settings determine when extra "failure" snapshots are rendered in the snapshot window.
        nextExplosionSnapshot = 5000;  % force a snapshot when number of edge failures reaches this threshold
        explosionSnapshotRatio = 10;   %and thereafter increase the forced-snapshot threshold by this ratio
        
        %% Simulation termination conditions
        % We don't want to waste time simulating lightsails that have veered off course or failed structurally.  Here are the
        % early termination criteria that will stop a simulation immediately:
        flyaway_ratio = 2.5;  % by how many times the initial sail radius the sail can be displaced (from the beam center) before the simulation stops.  Set to inf to disable.
        explosion_ratio = 6;  % by how many times the initial sail radius the tattered remains of the sail can extend (radially) before the simulation stops.  Set to inf to disable.
        failure_ratio = .4; %the fraction of triangles broken before termination.  Set >1 to disable
        
        % We can also opt to terminate the simulation after a certain time delay following certain situations.  The delay lets
        % us continue to visualize the mode of failure for a brief time before stopping the simulation.  
        terminate_on_upsidedown = 0;  %initiates a delayed termination if any of the triangles are upside down, which usually means the launch has failed.
        terminate_on_breakup = 1; %initiates a delayed termination upon the first tensile failure of the mesh.
        termination_delay_s = 0.10; %seconds of additional simulation time before any auto-shutoff conditions above terminate the simulation
        
        
        %% END OF COMMONLY CHANGED SETTINGS. 
        % We don't recommend changing values beyond this point unless you know what you're doing.
        
        
        %% Set up misc. constants, control variables, and initial temps...  you shouldn't need to modify these values.
        frame_interval = max(floor(frame_interval_s/dt),1);  %number of loops between movie frames
        snapshot_intvl = max(floor(snapshot_interval_s/dt),1);  % number of loops between snapshots
        %kick_interval = max(floor(kick_interval_s/dt),1);
        sim_dur_nt = ceil( (sim_dur + spinup_time) / dt );  %total number of iterations (time steps)
        sim_dur_nto = ceil(sim_dur_nt / sim_output_downsample);
        dnto = sim_output_downsample * dt;
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
        dispvar3('rigid_mode','I0','Irad');
        dispvar3('Iabs','Irefl','Emissivity');
        dispvar3('spinup_gs','spinup_time','videoMode');
        dispvar3('t_ramp_delay','t_ramp_dur','t_on_dur');
        dispvar3('dt','frame_interval','snapshot_intvl');
        dispvar2('sim_dur','sim_dur_nto');
        if (sim_dur_nto > 1e7)
            error('Number of output steps exceeds 10M.  Probably mistake?');
        end
        
        

        %logic check for reverse ray tracing
        if ReverseRayTrace 
            if ~usetether
                warning('Reverse ray tracing only works with a tethered payload enabled.  Disabling RRT!');
                ReverseRayTrace = 0;
            end
        end
        
        %% Initialize monitors:
        
        if ReverseRayTrace
            RRT_ang_bins = 201;
            RRT_bins = zeros(RRT_ang_bins, sim_dur_nto);
            RRT_ax = zeros(nt, sim_dur_nto);
            RRT_ay = zeros(nt, sim_dur_nto);
            RRT_a = RRT_ay;
        end
        
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
        maxstrainsnorm = zeros(1,sim_dur_nto);     % Strain of most-tensioned edge, normalized to tensile failure limit
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
        I0s = zeros(1,sim_dur_nto);            %Beam gate value (0 to 1)
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
        
        db_bm_timer_gpu = 0;
        db_bm_timer_vec = 0;
        db_bm_timer_itr = 0;
        
        % %hooprl = zeros(1, sim_dur_nto);
        % %hooprs = zeros(1,sim_dur_nto);
        % %hoopa = zeros(1,sim_dur_nto);
        % %monitorcs = zeros(length(mycs),sim_dur_nto);    %this can track any property along a cross section (across the membrane along Y)
        % %monitorfxs = zeros(length(n_x),sim_dur_nto);
        % %monitorfys = zeros(length(n_y),sim_dur_nto);
        
        fft_depth = min( fft_dur./dt, floor(5e8/length(n_x)/2)*2 );  % Number of full-grid monitor records to try to fit in memory.
        fft_num=0;
        if fft_mode
            monitorfzs = zeros(length(n_z),fft_depth);     %Mech. force on each node, in z direction
            monitorrfs = zeros(length(n_z),fft_depth);     %mech force on each node, projected to
            monitortfs = zeros(length(n_z),fft_depth);
            monitorxs = zeros(length(n_x),fft_depth);      %Complete position monitor,
            monitorys = zeros(length(n_y),fft_depth);      % All nodes
            monitorzs = zeros(length(n_z),fft_depth);      % This is where all the memory has gone!
            disp(['Note:  FFT mode analysis is enabled; buffer = ' num2str(fft_depth) ' frames (' num2str(6 * 8 * length(n_z) * fft_depth/1024/1024 ) ' MB)']);
            disp(['       FFT start time:  ' num2str(fft_start_time) '   duration: ' num2str(fft_dur) ]);
        end
        savemeshdepth = min( ceil(sim_dur / save_mesh_dt ), floor(1e8/length(n_z)/2)*2 );
        savemesh_num = 0;
        savemesh_last = -inf;
        if save_mesh_mode  % currrently I have this set up to study radiative heat transfer
            savemeshn_x = zeros(length(n_z), savemeshdepth);
            savemeshn_y = zeros(length(n_z), savemeshdepth);
            savemeshn_z = zeros(length(n_z), savemeshdepth);
            savemeshn_t = zeros(length(n_z), savemeshdepth);
            savemeshn_hf = zeros(length(n_z), savemeshdepth);
            %savemesht_t = zeros(length(t_na), savemeshdepth);
            savemesht_abs = zeros(length(t_na), savemeshdepth); 
            savemesht_abs_mr = zeros(length(t_na), savemeshdepth); 
            savemesht_ems = zeros(length(t_na), savemeshdepth); 
            %savemesht_cnd = zeros(length(t_na), savemeshdepth); 
            savemesht_rht = zeros(length(t_na), savemeshdepth); 
            savemesht_a = zeros(length(t_na), savemeshdepth); 
            savemesht_maxstr = zeros(length(t_na), savemeshdepth); 
            savemesh_ntos = zeros(1,savemeshdepth);
            disp(['Note:  Mesh saving mode for post-video rendering is enabled; frame depth = ' num2str(savemeshdepth) '  (' num2str(8 * savemeshdepth * ( 5*length(n_z) + 6*length(t_na) + 1 )/1024/1024 ) ' MB)']);
            disp(['       Mesh memory buffer duration: ' num2str( save_mesh_dt * savemeshdepth) ' s']);
        end
        
        %% Initialize buffers:
        t_optFcn_r2 = (zeros(size(t_na)));  % intermediate variable for calculating gaussian illumination
        t_ev1 = zeros(3,numtris); % edge vector for first triangle edge
        t_ev2 = zeros(3,numtris);  % edge vector for second edge
        t_ev1n = t_ev1; % normalized (unit) vector for #1 triangle edge
        t_ev2n = t_ev2; % normalized (unit) vector for #2 triangle edge
        t_norms = (zeros(3,length(t_na)));   % triangle normal vectors, normalized
        t_IDir = ([ zeros(size(t_na)); zeros(size(t_na)); ones(size(t_na)) ]);   %light incidence direction at each triangle centroid
         % note:  This initial value is not ever changed in the current version of the simulation, because we always assume light
         % oriented on the positive Z axis.  
        t_IMag = (zeros(size(t_na)));   % magnitude of light intensity at each triangle centroid (watts per area)
        t_Ipwr = (zeros(size(t_na)));   % incident light power on each triangle (watts)
        t_aproj = t_a0xy;  % triangle area in light incidence plane, i.e., projected area in xy plane
        t_specA = zeros(size(t_na));   % specular absorbed power
        t_specR = zeros(size(t_na));   % specular reflection
        t_abs = (zeros(size(t_na)));   % heat: absorbed optical power in each triangle (watts)
        t_abs_mr = (zeros(size(t_na)));% heat: absorbed optical power in each triangle (watts) from multiple reflections (not including first)
        t_ems = (zeros(size(t_na)));   % heat: emitted power from triangle
        t_cnd = (zeros(size(t_na)));   % heat: conducted power from (to?) each tringle 
        t_rht = (zeros(size(t_na)));   % heat: absorbed radiative heat from other triangles 
        t_reflp = (zeros(size(t_na)));  %reflected optical power in each triangle (watts)
        t_refln = zeros(3,length(t_na)); % reflected light direction (unit vector)
        t_oth = zeros(3,length(t_na));  % optical thrust vectors (force, newtons)
        t_oth_mr = zeros(3,length(t_na)); % optical thust vectors, multi-reflection (force, newtons)
        %t_vel = (zeros(3,length(t_na)));  % velocity of each triangle   ---->  not used?  Will delete soon.
        t_texn = (zeros(3,length(t_na)));  % unit vector within each triangle plane, indicating orientation axis of texture ("parallel" direction)
        t_textn = (zeros(3,length(t_na)));  % unit vector within each triangle plane, perpendicular to the texture orientation ("transverse" direction)
        t_press_n = (zeros(size(t_na)));  %LUT photon pressure coefficient, normal to triangle
        t_press_p = (zeros(size(t_na)));  %LUT photon pressure coefficient, parallel to triangle
        t_press_t = (zeros(size(t_na)));  %LUT photon pressure coefficient, transverse to triangle
        t_lut_a = (zeros(size(t_na)));    %LUT absorption coefficient for each triangle
        e_dl = (zeros(1,length(e_na)));  % mm, 
        e_nl = (zeros(3,length(e_na)));  % normalized (unit) direction vector of each edge
        e_l = (zeros(size(e_na))); %mm, current edge length
        e_l0t = zeros(size(e_na)); %mm, current unstrained edge length, including thermal expansion
        % already initilaized in generateMesh:  n_a = (zeros(size(n_x)));  % mm2, effective radiation area for each node, for radiative cooling calcs         
        n_abs = (n_a);  % absorbed power (watts), from propulsion beam, heating input to each node
        n_ems = (n_a);  % emitted power (watts), from radiative cooling, heating output from each node
        n_hf = (n_a);  % heat flow (watts), from thermal conduction, heat into or out of each node
        
        
        
        %% Initialize various simulation varilables, state variables, and flags prior to beginning time-domain simulation
        % There are no user-configurable settings in this section.  Just initialization stuff.
        %angVel = [0; 0; 0];
        %COMVel = [0; 0; 0];
        nt = 0;
        nto = 0;
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
                      'waitGPU' };            %15
        
        terminate_t = 0;  % (seconds) value for delayed termination
        lastdispnt = 0;
        lastoutofboundtswarn = -1; % set to -1 to disable console warnings about upside down triangles
      
        
        debug_nt_for_pause = -1;
        
        
       
        %% start timers for simulation loop
        tic;
        startsimtic = tic;
        lastdisptic = startsimtic;
        setup_time = toc(startlooptic);
        lastbuttoncheck = tic;
        
        %debug
        specialtimer1 = 0;
        specialtimer2 = 0;
        
        %error check
        if useMepsLUTs && max(t_tex) > 0
            if RayTracingMode
                error('RayTracingMode and useMepsLUTs options cannot be used at the same time.');
            end
        end

        % need these for t=0 tilt rotations:
        DRM123_11 = @(psi,theta,phi) cos(psi).*cos(theta);
        DRM123_12 = @(psi,theta,phi) cos(psi).*sin(theta).*sin(phi) + sin(psi).*cos(phi);
        DRM123_13 = @(psi,theta,phi) -cos(psi).*sin(theta).*cos(phi) + sin(psi).*sin(phi);
        DRM123_21 = @(psi,theta,phi) -sin(psi).*cos(theta);
        DRM123_22 = @(psi,theta,phi) -sin(psi).*sin(theta).*sin(phi) + cos(psi).*cos(phi);
        DRM123_23 = @(psi,theta,phi) sin(psi).*sin(theta).*cos(phi) + cos(psi).*sin(phi);
        DRM123_31 = @(psi,theta,phi) sin(theta);
        DRM123_32 = @(psi,theta,phi) -cos(theta).*sin(phi);
        DRM123_33 = @(psi,theta,phi) cos(theta).*cos(phi);
        
        %% (V18 new) move everything to the GPU for the simulation
%         disp('Moving stuff to GPU');
% %         % mesh stuff & constant arrays
% %         
% %         % node mesh constants
% %             %not needed on gpu: n_ang0, n_r0
% %             g_n_ea = gpuArray(n_ea);
% %             g_n_ne = gpuArray(n_ne);
% %             g_n_ta = gpuArray(n_ta);
% %             g_n_nt = gpuArray(n_nt);
% %             g_n_m = gpuArray(n_m);
% %         % edge mesh constants
% %             %Not needed on gpu:  e_f0   e_l1 e_l2  e_reg   e_t1 e_t2  e_thk1,2  e_w1,2+b    
% %             g_e_al = gpuArray(e_al);
% %             g_e_k = gpuArray(e_k);
% %             g_e_l0 = gpuArray(e_l0);
% %             g_e_m = gpuArray(e_m);
% %             g_e_na = gpuArray(e_na);
% %             g_e_nb = gpuArray(e_nb);
% %         % triangle mesh constants
% %             %Not needed on gpu:  g_t_a0zy  t_centroid_ang,r  t_x0,y0,z0
% %             %t_f0      t_reg(obsolete) t_thick
% %             g_t_m = gpuArray(t_m);
% %             g_thissemsdumb = gpuArray(thisseemsdumb);
% %             g_t_a0 = gpuArray(t_a0);
% %             g_t_costexa = gpuArray(t_costexa);
% %             g_t_sintexa = gpuArray(t_sintexa);
% %             g_t_e1 = gpuArray(t_e1);
% %             g_t_e2 = gpuArray(t_e2);
% %             g_t_e3 = gpuArray(t_e3);
% %             g_t_na = gpuArray(t_na);
% %             g_t_nb = gpuArray(t_nb);
% %             g_t_nc = gpuArray(t_nc);
% %             g_t_tex = gpuArray(t_tex);
% %             g_t_texa = gpuArray(t_texa);
% %             
% %         
% %         
% %         % node vectors and arrays
% %         g_n_x = gpuArray(n_x);
% %         g_n_y = gpuArray(n_y);
% %         g_n_z = gpuArray(n_z);
% %         g_n_t = gpuArray(n_t);
% %         g_n_a = gpuArray(n_a);
% %         g_n_af = gpuArray(n_af);
% %         g_n_dr = gpuArray(n_dr);
% %         g_n_dx = gpuArray(n_dx);
% %         g_n_dy = gpuArray(n_dy);
% %         g_n_dz = gpuArray(n_dz);
% %         g_n_dx2 = gpuArray(n_dx2);
% %         g_n_dy2 = gpuArray(n_dy2);
% %         g_n_dz2 = gpuArray(n_dz2);
% %         g_n_dxyz = gpuArray(n_dxyz);
% %         g_n_fx = gpuArray(n_fx);
% %         g_n_fy = gpuArray(n_fy);
% %         g_n_fz = gpuArray(n_fz);
% %         g_n_mf = gpuArray(n_mf);
% %         g_n_of = gpuArray(n_of);
% %         g_n_vx = gpuArray(n_vx);
% %         g_n_vy = gpuArray(n_vy);
% %         g_n_vz = gpuArray(n_vz);
% %         g_coms_x = gpuArray(coms_x);
% %         g_coms_y = gpuArray(coms_y);
% %         g_coms_z = gpuArray(coms_z);
% %         g_n_abs = gpuArray(n_abs);  % absorbed power (watts), from propulsion beam, heating input to each node
% %         g_n_ems = gpuArray(n_ems);  % emitted power (watts), from radiative cooling, heating output from each node
% %         g_n_hf = gpuArray(n_hf);  % heat flow (watts), from thermal conduction, heat into or out of each node
% %         
% %         % edge vectors and arrays
% %         g_e_t = gpuArray(e_t);
% %         g_e_broken = gpuArray(e_broken);
% %         g_e_notbroken = gpuArray(e_notbroken);
% %         g_e_dl = gpuArray(e_dl);  % mm, 
% %         g_e_nl = gpuArray(e_nl);  % normalized (unit) direction vector of each edge
% %         g_e_l = gpuArray(e_l); %mm, current edge length
% %         g_e_l0t = gpuArray(e_l0t ); %mm, current unstrained edge length, including thermal expansion
% %         
% %         % triangle vectors and arrays
% %         g_t_cx = gpuArray(t_cx);
% %         g_t_cy = gpuArray(t_cy);
% %         g_t_cz = gpuArray(t_cz);
% %         g_t_t = gpuArray(t_t);  %temperature
% %         g_t_broken = gpuArray(t_broken);
% %         g_t_notbroken = gpuArray(t_notbroken);
% %         g_t_optFcn_r2 = gpuArray(t_optFcn_r2);  % intermediate variable for calculating gaussian illumination
% %         g_t_ev1 = gpuArray(t_ev1); % edge vector for first triangle edge
% %         g_t_ev2 = gpuArray(t_ev2);  % edge vector for second edge
% %         g_t_ev1n = gpuArray(t_ev1n); % normalized (unit) vector for #1 triangle edge
% %         g_t_ev2n = gpuArray(t_ev2n); % normalized (unit) vector for #2 triangle edge
% %         g_t_norms = gpuArray(t_norms);   % triangle normal vectors, normalized
% %         g_t_IDir = gpuArray(t_IDir);   %light incidence direction at each triangle centroid
% %         g_t_IMag = gpuArray(t_IMag);  % magnitude of light intensity at each triangle centroid (watts per area)
% %         g_t_Ipwr = gpuArray(t_Ipwr);  % incident light power on each triangle (watts)
% %         g_t_aproj = gpuArray(t_aproj);  % triangle area in light incidence plane, i.e., projected area in xy plane
% %         g_t_specA = gpuArray(t_specA);  % specular absorption
% %         g_t_specR = gpuArray(t_specR);  % specular reflection
% %         g_t_abs = gpuArray(t_abs);  %absorbed optical power in each triangle (watts)
% %         g_t_abs_mr = gpuArray(t_abs_mr);  %absorbed optical power in each triangle (watts)
% %         g_t_reflp = gpuArray(t_reflp);  %reflected optical power in each triangle (watts)
% %         g_t_refln = gpuArray(t_refln); % reflected light direction (unit vector)
% %         g_t_oth = gpuArray(t_oth);  % optical thrust vectors (force, newtons)
% %         g_t_oth_mr = gpuArray(t_oth_mr); % optical thust vectors, multi-reflection (force, newtons)
% %         g_t_texn = gpuArray(t_texn);  % unit vector within each triangle plane, indicating orientation axis of texture ("parallel" direction)
% %         g_t_textn = gpuArray(t_textn);  % unit vector within each triangle plane, perpendicular to the texture orientation ("transverse" direction)
% %         g_t_press_n = gpuArray(t_press_n);  %LUT photon pressure coefficient, normal to triangle
% %         g_t_press_p = gpuArray(t_press_p);  %LUT photon pressure coefficient, parallel to triangle
% %         g_t_press_t = gpuArray(t_press_t);  %LUT photon pressure coefficient, transverse to triangle
% %         g_t_lut_a = gpuArray(t_lut_a);    %LUT absorption coefficient for each triangle
% % 
% %         
% %         
% %         % helper sparse matrices:
% %         g_M_e2n = gpuArray(M_e2n);
% %         g_M_t2n = gpuArray(M_t2n);
% %         g_M_nxyz2tev1 = gpuArray(M_nxyz2tev1);
% %         g_M_nxyz2tev2 = gpuArray(M_nxyz2tev2);
% %         
% %         
% %         % other stuff I want to leave on the GPU
% %         g_inertiaT = gpuArray(inertiaT);
% %         
% %         %LULs
% %         if useSpecularTMM
% %             g_LUL = [];
% %             g_LUL.R = gpuArray(LUL.R);
% %             g_LUL.A = gpuArray(LUL.A);
% %         end
% %         
% %         %LUTs
% %         g_LUTs = [];
% %         for n=1:length(LUTs)
% %             g_LUTs(n).n = LUTs(n).n;
% %             g_LUTs(n).p = LUTs(n).p;
% %             g_LUTs(n).t = LUTs(n).t;
% %             g_LUTs(n).a = LUTs(n).a;
% %         end
% %         %asdf = gather(g_t_lut_a);  % need to let GPU catch up with array creation, otherwise we crash
% %         
%         disp('Done moving stuff to GPU, starting loop');
            
        RT_bmechmark_count = 0;
        lutmask = t_tex>0;
        anyluts = any(lutmask);

        %these are shared wiht plotting subscripts:
        xbox = 1.10*(radiusmm + sqrt(t0t_x0^2+t0t_y0^2));% 
        zbox = 0.2*xbox + xbox*zAspectRatio/2; %(0.8*radiusmm*max(1/(pD/radialRings),0.8)) * (1+1/radialRings);
        if (pD == 0) zbox = xbox*0.8; end


        %% LATERAL STABILITY VISUALIZER
        % perform optional lateral stability vs. tilt analysis thing
        if doStabilityVisualization
            tiltStabilityAnalysis_v20;
            if doStabilityVisualization>1
                return;
            else
                close(gcf);
            end
        end


        %% Setup video window
        if (videoMode < 2)
            V1 = VideoWriter([ filebasename '_video' ],'Archival');  %alternate format for Ramon:  'Archival'
            %V1.Quality = 100;
            open(V1);
        end
        if (videoMode==1)
            %fSimVid = figure('pos',[100 100 1600 1200]);
            fSimVid = figure('pos',[50        50        1280         1024]);  % for my laptop!
        end
        if videoMode==0
            fSimVid = figure('pos',[100 100 2000 800]);
        end
        if videoMode<2
            set(gcf,'color','w');     
        end
        
        %% set up rigid/flexible mode selector, and simulation abort button
        % No longer use global rigid_mode flag, we now support multiple rigid bodies 
%         if rigid_mode
%             MODEBOX = stoploop({'Currently using rigid dynamics.' 'Press OK to change to flexible.'}) ;
%         else
%             MODEBOX = stoploop({'Currently using flexible dynamics.' 'Press OK to change to rigid.'}) ;
%         end
         FS = stoploop({'Press this button' 'to terminate simulation.'}) ;
        
        %% Setup snapshot window
        fSnapshot = figure('pos',[200 100 400 800]);  set(gcf,'color','w');

        plotSnapshot;
        hSnapshotLight2 = lightangle(00,20);
        hSnapshotLight = lightangle(180,40);
        axis equal;
        colormap jet;
        colorbar ;

        %% Render / record first video frame
        if videoMode < 2
            figure(fSimVid);
        end
        plotMesh_v17;
        if (videoMode < 2)
            writeVideo(V1,getframe(gcf));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Main time-domain simulation loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nt=1;
        while nt<=sim_dur_nt 
            %% Notes
            %  In addition to state variables, the following derived variables should be correctly established at the
            %  beginning of each iteration:
            %       inertiaT -- intertia tensor
            %       angVel -- mesh-calculated angular velocity (broken)
            %       t_cx...z -- triangle centroid positions
            %       comx...z -- center of mass
            %       n_dx...z -- node positions relative to COM
%  disp(['nt=' num2str(nt) ]);

    dooptics = ~mod(nt-1, opt_skip_downsample);
            
            %% 0. Calculate spin-up forces
            if  (tt < 0) && rps_target  % spin-up forces are applied for tt<0.  
                n_arad = sqrt( n_x.^2 + n_y.^2 );
                n_af(1,:) = sign(rps_target)*spinup_gs/88.2.*n_m.*(n_y) ;%%% .* (n_arad < 0.6 * radiusmm);
                n_af(2,:) = sign(rps_target)*-spinup_gs/88.2*n_m.*(n_x) ;%%%.* (n_arad < 0.6 * radiusmm);
                %n_af(3,:) = spinup_gs/88.2.*n_m.*(n_y) ;%%% .* (n_arad < 0.6 * radiusmm);
                %n_af(2,:) = -spinup_gs/88.2*n_m.*(n_z) ;%%%.* (n_arad < 0.6 * radiusmm);
            end
            if ( tt >= 0) || anybroken  % stop spin-up at t=0 (or if any edges broke during spinup)
                n_af(:)=0;
            end
  %disp('done with spinup forces');
            
            %% 1.  Calculate beam intensity profile and direction
            % here we set t_IMag and t_IDir, which are the magnitude and direction of the incident light at each
            % triangle centroid.  Normally t_IDir is (0,0,1) so we never change it here, but in theory you could set
            % whatever value you want.  Note:  t_IDir must be unit vectors.  
            pwrRampVal = pwrRampFcn(tt);
%             if pwrRampVal < 0.9999
%                 dooptics = 1;
%             end
            if (pwrRampVal == 0)
                t_IMag = zeros(size(t_cx));
            else
                t_optFcn_r2 = (t_cx - 0).^2 + (t_cy - 0).^2;
                t_IMag = pwrRampVal * I0 * exp( - t_optFcn_r2 ./ (Irad^2) );  % This is in W/mm2.  Divide by c0 to get newtons / mm2 radiation pressure (N/mm2)
                %t_iz(:) = 1; can skip this since we never change incidence vector
                %for now...
                if nto_startAccel==0
                    nto_startAccel = nto;
                    nt_startAccel = nt;
                end

            end
            
            % keep track of time spent on spin-up force and light intensity calculations
            timers(1) = timers(1) + toc;
            tic;
    %asdf = gather(t_optFcn_r2);
  % disp('Done with intensity');
            

            %% 2. Calculate triangle area, normal vectors, and texture orientation vectors
            if isbig  %choose between two calculation mechanisms for optimal performance on my home computer (MK)
                t_ev1(1,:) = n_x(t_nb) - n_x(t_na);
                t_ev1(2,:) = n_y(t_nb) - n_y(t_na);
                t_ev1(3,:) = n_z(t_nb) - n_z(t_na);
                t_ev2(1,:) = n_x(t_nc) - n_x(t_na);
                t_ev2(2,:) = n_y(t_nc) - n_y(t_na);
                t_ev2(3,:) = n_z(t_nc) - n_z(t_na);
            else
                t_ev1 = [n_x; n_y; n_z]* M_nxyz2tev1;
                t_ev2 = [n_x; n_y; n_z]* M_nxyz2tev2;
            end
            t_norms = cross( t_ev1, t_ev2 );  %triangle normal vectors, not yet normalized; we will normalize shortly
   % asdf = gather(t_norms);
%  disp('got norms');
            
            upsidedowntris = sum( t_norms(3,:) < 0);  % we determine if a triangle is upside down based on the z 
                % component of its normal vector.  Here there's an implicit assumption that light is incident along the 
                % positive Z axis, so if you ever set t_IDir to something other than (0,0,1), you might need to
                % determine a new way to detect upside down triangles.  
            
            if upsidedowntris
                if ~terminate_t && (terminate_on_upsidedown)  %initiate delayed termination if triangles are upside down
                    terminate_t = tt + termination_delay_s;
                    disp('Initiating delayed termination due to upside down triangles!');
                end
            end
 %asdf = gather(upsidedowntris);
 % disp('tested norms');
            
            % Calculate length of the non-normalized normal vectors, which happen to also be 2x the triangle area
            t_a = max(vecnorm(t_norms), 1e-8);  %initially this is 2x the triangle area
            
            % Normalize normal vectors to each triangle
            t_norms = t_norms ./ t_a ; %ok, now t_norms contains proper unit normal vectors
  % asdf = gather(t_norms);
%disp('normalized norms');            
            
            t_a = 0.5 .* t_a; % ok, now t_a contains the actual area of each triangle
            
            if anybroken
            % If any triangles are 'broken', we replace their now-meanlingless calculated area with their original area
                t_a(t_broken) = t_a0(t_broken);  % Don't worry, other calcs will continue ignoring broken triangles.  
            end
            
%disp('set area');            
            % Normalize the #1 edge vectors
            %t_ev1n = t_ev1./(max(vecnorm(t_ev1),1e-8));
            %test =;
%disp('did vecnorm');
            t_ev1n = t_ev1./ vecnorm(t_ev1);%test;
%disp('got first edge');            
            
            % During mesh generation, we recorded a value 't_texa' for each triangle, which is the in-plane angle
            % between the texture orientation vector and the #1 edge vector.  We use that now to calculate the 
            % new texture orientation vectors for the current mesh state... well, actually, we use the cached values 
            % of the sin and cosine of this angle to hopefully save time (t_costexa, t_sintexa).
            t_texn = t_ev1n .* (t_costexa) + ...
                cross(t_norms, t_ev1n ) .* (t_sintexa)  + ...
                t_norms .* dot(t_norms, t_ev1n ) .* (1 - cos(t_texa)) ;  %TODO:  Can I safely replace the cos() with t_costexa???

% disp('got texn')
    
            % Calculate cos(theta), where theta corresponds to the light incidence angle at each triangle.              
            t_costheta = dot(t_norms, t_IDir);
            
            % Here we calculate the projected area of non-broken triangles in the X-Y plane.  This is essentially the
            % aperture area of our sail in the plane of beam incidence.  Note again the implicit assumption that
            % illumination is along (0,0,1).  
            t_aproj = t_a .* t_costheta .* (t_notbroken);
            
 %  disp('Done with triangle stuff');
            
            %% 3 Calculate reflection and absorption for specular regions
            %  Currently we use simple scalar values for reflection and absorption.  If you were so inclined, you could
            %  implement different physics to calculate how much light should be reflected and absorbed by specular sail
            %  surfaces, such as using fresnel equations.  However, to implement more complex optical behavior, I 
            %  recommend using look-up tables rather than trying to calculate values here.
            %
            %  Available inputs for calculating specular response:
            %  t_costheta: cosine of angle of incidence
            %  t_t: local temperature
            %  t_x0, t_y0, t_z0: Initial coordinates of this triangle (prior to simulation beginning)
            %  t_thick: thickness of film here 
            %  t_m: mass of this triangle
            %
            %  Required outputs from this section: 
            %  t_specA, t_specR, absorption and reflection coefficients (coefficients--not power!  Values should be 0 to 1)
                
            if dooptics
                t_specA(:) = t_Iabs;   % default specular reflection
                t_specR(:) = t_Irefl;
                if useSpecularTMM  % will be overwritten by specularTMM if enabled and has material
                    t_theta = acos(t_costheta);
                    if min(t_theta) < 0 
                        warning('Negative theta in useSpecularTMM code');
                        t_theta(t_theta < 0) = -t_theta(t_theta < 0);  % I don't think this should ever happen
                    end
                    if max(t_theta) > pi/2
                    %    warning('Reverse illuminated triangles in specularTMM code');
                        t_theta(t_theta > pi/2) = pi/2 - t_theta(t_theta > pi/2); % flip the angle so we get the correct LUL value
                    end
                    t_LULi = max(min(floor(t_theta ./ LULAngleStep)+1,LUL.num),1);
                    t_specA = LUL.A(  sub2ind(LUL.sz, t_mat, t_LULi) );
                    t_specR = LUL.R( sub2ind(LUL.sz, t_mat, t_LULi) ) ;

                    % will presumably have to iterate over materials
                end

     %  disp('Done with specular');


                %% 4. Calculate per-triangle incident power, specular optical forces, & absorption heat-loading
                t_Ipwr = abs(t_aproj.*t_IMag);  %incident power on each triangle, watts (includes costheta and broken factors)
                if anybroken
                    t_Ipwr(t_broken)=0;
                end

                % Calculate specular absorbed power in each triangle.  This gives the absorption heat input.  
                % Note:  Radiatiave cooling occurs later at the node level.  That's why triangle area is distributed back 
                %  to the nodes during each simulation step.  
                t_abs = t_specA .* t_Ipwr; % Watt
                % Note:  This value will be overwitten in textured regions via the LUT code

                % Calculate specular reflected power.
                t_reflp = t_Ipwr .* t_specR;

                % calculate optical force ("optical thrust") for each triangle
                t_oth = 2 * t_reflp ./ c0 .* t_norms .* t_costheta  + ... %  this is the reflection force
                    t_abs       ./ c0 .* t_IDir ;    % and this is the absorption force
                % Note: this will be overwritten in textured regions with LUT values.
            end

            
            
%             if pwrRampVal > 0
%                  t_refln = t_IDir - 2 .* t_costheta .* t_norms;
%                  t_vert0 = [n_x(t_na); n_y(t_na); n_z(t_na)];
%                % figure
%                % plotMesh_v17;
%                % hold on;
%                % quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
%                 return;
%             end
              
            % keep track of time spent doing triangle-based calculations (and specular optics)
            timers(2) = timers(2) + toc;
            tic;
            timerStartRT = tic;
            
            
            %% 5b. Do RayTracingMode optics  
            if dooptics && RayTracingMode && (pwrRampVal > 0)
                
                t_refln = t_IDir - 2 .* t_costheta .* t_norms;
                t_vert0 = [n_x(t_na); n_y(t_na); n_z(t_na)];

                
                processMR_GPU_v17;  % call raytracing subscript.  Generates:
                      %  t_oth_mr :  additional optical thrust from multiple reflections
                      %  t_abs_mr :  additional absorption from multiple reflections.
                      
                      
%              db_bm_timer_gpu = db_bm_timer_gpu + db_bm_timer_gpu_val;
%              db_bm_timer_vec = db_bm_timer_vec + db_bm_timer_vec_val;
%              db_bm_timer_itr = db_bm_timer_itr + db_bm_timer_itr_val;
             
                timers(14) = timers(14) + toc(timerStartRT); % the RT code sometimes calls tic() on its own, so need to explicitly reference the toc() to the correct time.
                tic;
                
%                 RT_bmechmark_count = RT_bmechmark_count + 1;
%                 if RT_bmechmark_count == 301    % 123456789012                          Benchmark timers:  
%                     fprintf('Benchmark timers:      GPU    \t        VEC    \t   ITR   (numtris=%d)  \n                   ', numtris);
%                     fprintf('%11f \t', [ db_bm_timer_gpu db_bm_timer_vec db_bm_timer_itr]./300    );
%                     return;
%                 end
            end

            if dooptics && useMepsLUTs && (pwrRampVal>0) && anyluts
            %% 5b. Calculate optical response of non-specular ('textured') regions via look-up tables (LUTs).  
            % Look-up tables are used to specify arbitrary optical properties of textured regions of the mesh.  For
            % those regions, we look up pre-calculated values for absorption and optical thrust, based on the 
            % incidence angle relative to the triangle surface and texture orientation.  
            % For these calculations we are relying on the following vectors as the starting point. 
            %     t_norms :  unit normal vectors for each triangle
            %     t_texn :   unit vectors aligned with the grating ('parallel' axis)
            %     t_Idir:  unit vectors describing light incidence direction
            
                % First we generate the transverse unit vectors for each triangle
                t_textn = cross(t_norms, t_texn);  % this is the transverse unit vector

                % sanity checks, we can probably delete this, but it's nice to know if we've messed up:
                if mod(nt,1000)==0  %to save time, only do this check once every 1000 loops
                    %      disp([ 'This should be zero: ' sprintf('%e', max(abs( dot( t_norms, t_texn) ) ) ) ] );
                    %check for orthagonolity of t_norms and t_texn
                    if (max(abs( dot( t_norms, t_texn) ) ) >1e-15)
                        warning('Error with vector math in optical force LUT code...  the triangle normals aren''t orthogonal to the texture vectors!');
                    end
                    %      disp([ 'This should be zero: ' sprintf('%e', max(abs( dot( t_textn, t_textn ) - 1) ) ) ]);
                    %check that t_textn comprises unit vectors
                    if (max(abs( dot( t_textn, t_textn ) - 1) ) > 1e-14)
                        warning('Error with vector math in optical force LUT code...  the transverse texture vector isn''t normalized.  Something is wrong.');
                    end
                end

                % now we're going to project the incidence vector into the parallel and transverse grating planes.

                t_incsignt = dot( t_IDir, t_textn );  %projection mag along transverse vector
                t_incprojp = t_IDir - t_textn .* t_incsignt ; % inc vector projected into parallel plane
                t_incprojpn = t_incprojp ./ vecnorm(t_incprojp) ; % now normalized

                t_incsignp = dot( t_IDir, t_texn );  % projection mag along parallel vector
                t_incprojt = t_IDir - t_texn .* t_incsignp ; % inc vector projected into transverse plane
                t_incprojtn = t_incprojt ./ vecnorm(t_incprojt) ; % now normalized

                t_incsignt = sign(t_incsignt);  % get sign for angle unwrap from acos function
                t_incsignp = sign(t_incsignp);  % ...

                t_psi = t_incsignp .* acos( max(min(dot( t_incprojpn, t_norms ), 1),-1) );  % the min/max confine us to real limits for acos().  
                t_theta = t_incsignt .* acos( max(min(dot( t_incprojtn, t_norms ), 1),-1) );% they shouldn't be needed, but I don't trust floating
                                                                                            % point arithmatic quite enough to omit them

                t_ipsi = round((t_psi+pi/2)./deg2rad(LUTPsiStepDeg)); % Psi index for LUTs
                t_itheta = round((t_theta+pi/2)./deg2rad(LUTThetaStepDeg)); % Theta index for LUTs

                % check for (and correct) out-of-bounds psi and theta values, which occur for reverse-illuminated triangles
                outofboundspsilow = t_ipsi < 1;
                outofboundspsihigh = t_ipsi > (180/LUTPsiStepDeg+1);
                outofboundsthetalow = t_itheta < 1;
                outofboundsthetahigh = t_itheta > (180/LUTThetaStepDeg+1);
                noobpl = sum(outofboundspsilow);
                noobph = sum(outofboundspsihigh);
                noobtl = sum(outofboundsthetalow);
                noobth = sum(outofboundsthetahigh);
                outofboundts = noobpl + noobph + noobtl + noobth;
                if (outofboundts)
    %                 if (outofboundts ~= lastoutofboundtswarn) && (lastoutofboundtswarn > 0)
    %                     lastoutofboundtswarn = outofboundts;  
    %                     warning(['There are ' num2str(outofboundts) ' reverse-illuminated triangles at nt=' num2str(nt)  '!  ' ...
    %                         '(psi_low = ' num2str(noobpl) ',  psi_high = ' num2str(noobph) ...
    %                         ',  theta_low = ' num2str(noobtl) ',  theta_high = ' num2str(noobth) ')' ]);
    %                 end
                    if noobpl t_ipsi(outofboundspsilow) = 1; end
                    if noobph t_ipsi(outofboundspsihigh) = 180/LUTPsiStepDeg+1; end
                    if noobtl t_itheta(outofboundsthetalow) = 1; end
                    if noobth t_itheta(outofboundsthetahigh) = 180/LUTThetaStepDeg+1; end
                end

                % actually perform the look-ups:
                for nlut = 1:length(LUTs)
                    lutinds = sub2ind(size(LUTs(nlut).n), t_ipsi(t_tex == nlut),t_itheta(t_tex == nlut));
                    if any(lutinds)
                        t_press_n(t_tex == nlut) = LUTs(nlut).n(lutinds);
                        t_press_t(t_tex == nlut) = LUTs(nlut).t(lutinds);
                        t_press_p(t_tex == nlut) = LUTs(nlut).p(lutinds);
                        t_lut_a( t_tex == nlut ) = LUTs(nlut).a(lutinds);
                    end
                end

                % calculate optical thrust and absorbed power based on LUT values
                
                if any(lutmask)
                    t_oth(:,lutmask) = t_Ipwr(lutmask) .*  ( ...
                        t_norms(:,lutmask) .* t_press_n(lutmask)  + ... this is the normal force
                        t_textn(:,lutmask) .* t_press_t(lutmask)  + ... this is the transverse force
                        t_texn(:,lutmask) .* t_press_p(lutmask)); % this is the parallel/axial force
						%technically, I should also add absorption (force) here but it's always negligable
                    t_abs(lutmask) = t_lut_a(lutmask) .* t_Ipwr(lutmask);
                end
            

                % keep track of time spent on LUT calculations
                timers(3) = timers(3) + toc;
                tic;
            end

            %% Calculate emission and radiative heat transfer
            %% Calculate thermal emission (radiatively emitted power)
            % Calculate radiated power using Stefan-Boltzmann's law given by P_rad
            % = A * emissivity * sigma * ( T^4 - T_ref^4). Here, the radiated power
            % at each node is calculated.
            %n_ems = SBCmm .* n_a .* Emissivity .* (n_t.^4); % Watt
            % UPDATE:  We don't radiate from nodes, we radiate from triangles, and have already distributed that to n_ems above
            if tt>t_ramp_delay   % also, I disable radiative heat treansfer / cooling until we start turning on the beam
                t_t4 = t_t.^4;
                t_ems = SBCmm .* t_a .* (t_emis + t_emis_fr) .* (t_t4); %watt, thermal radiated power from each triangle in all directions

                if dooptics && do_rht  % if doing radiative heat transfer
                    % Note: this would be an optional time to re-calculate the RHT conductance array if the geometry or optical properties have changed
                    % significantly during simulation... for now, I just use the t=0 RHT matrix, under the assumption the shape doesn't change much.
    
                    % We need to calculate heat ARRIVING to each tringle, from all other triangles
                    
                    % here we calculate the heat radiated by each triangle towards all other triangles
                    M_rht = (t_t4) .* t_rht_aemvfs;
                    %as sanity check, we calculate heat EMITTED from each triangle which impinges any other triangle (this should be the same)
    
                    t_rht = SBCmm .* sum(M_rht,2)';
                    %t_rht_out = sum(M_rht,1);
                    if anybroken
                        t_rht(t_broken) = 0;
                    end
                end
                if anybroken
                        t_ems(t_broken) = 0;
                end
            end
                
            
            %% 6 Distribute optical forces (& emission & radiative heat transfer) to the nodes
            if dooptics
                if RayTracingMode
                    if isbig  %this selects between two calculation methods, based on performance thresholds on my computer (MK).
                        n_of(1,:) = (t_oth(1,:) + t_oth_mr(1,:)) * M_t2n ./ 3;
                        n_of(2,:) = (t_oth(2,:) + t_oth_mr(2,:)) * M_t2n ./ 3;
                        n_of(3,:) = (t_oth(3,:) + t_oth_mr(3,:)) * M_t2n ./ 3;
                    else
                        n_of = (t_oth + t_oth_mr) * M_t2n ./ 3;
                    end
                else
                    if isbig  %this selects between two calculation methods, based on performance thresholds on my computer (MK).
                        n_of(1,:) = t_oth(1,:) * M_t2n ./ 3;
                        n_of(2,:) = t_oth(2,:) * M_t2n ./ 3;
                        n_of(3,:) = t_oth(3,:) * M_t2n ./ 3;
                    else
                        n_of = t_oth * M_t2n ./ 3;
                    end
                end


                % distribute absorbed power and effective radiator area to each node
                % update:  we no longer distribute radiator area to the nodes, we distribute radiated (and absorbed) POWER to each node

                t_p2n = t_abs - t_ems;
                if do_rht
                    t_p2n = t_p2n + t_rht;
                end
                if (RayTracingMode)
                    t_p2n = t_p2n + t_abs_mr;
                end
                
                n_ht_flux = t_p2n * M_t2n ./ 3;
                

                % quick sanity check, if we've messed up a calulation, terminate the simulation immediately so we can debug
                % before the error corrupts the entire mesh.
                if (isnan(sum(n_abs)))
                    error('Found NaN!!! (in optical forces)');
                end
            end
            % keep track of time spent on triangle-to-node calculations
            timers(4) = timers(4) + toc;
            tic;
            
                        

            
            
            %% 8. Calculate mechanical forces & heat conduction
            
            % Matrix containing all edge direction vectors with x, y, z components as columns
            e_nl = [ n_x(e_nb) - n_x(e_na) ; n_y(e_nb) - n_y(e_na) ; ...
                n_z(e_nb) - n_z(e_na) ];
            
            % Calculate current lengths in mm of all edges
            e_l = vecnorm( e_nl);
            
            % Normalize all edge direction vectors. Here 'e_nl' stands for edge, normalized length
            e_nl = e_nl ./ e_l;   %unit vectors
            
            % Calculate linear expansion due to temperature difference e_t minus initial temperature t0 for each edge 
            % multiplied by the coefficient of thermal expansion.
            % Amount of thermal expansion can be described by material strain
            % epsilon_thermal = (L_final - L_initial)/L_initial, where the strain
            % is due change in temperature proportional to the coefficient of
            % thermal expansion, epsilon_thermal = -alpha_L * (T_final - T_initial)
            % This can be rewritten in terms of length change dL / L_initial =
            % 1 - (1 + alpha_L * (T_final - T_initial) )
            
            e_tex = 1 + (e_t - t0).*(CTE*1e-6);  % edge length ratio due to thermal expansion
            e_l0t = e_l0 .* e_tex;  % resting (unstrained) edge length, including thermal expansion
            e_dl = e_l - e_l0t ;    % mm, current edge length difference vs. unstrained edge length, positive for tension
            
            % Calculate strain based on current vs. resting/unstrained length
            e_s = e_l ./ (e_l0t) - 1;   % strain %/100.   No stress -> 0 strain.  Elongate to 2x length -> 1.0 strain
            

            % Calculate mechanical force on edge as the the edge length difference
            % times the edge stiffness/spring constant being equal to the 2D
            % Young's modulus x area / length^2 (see generateMesh.m for more info).
            % The sign of the mechanical force on a specific edge is determined by
            % the corresponding edge direction vector 'e_nl'
            e_mf = (e_notbroken & e_flexible) .* (e_k) .* e_dl .* e_nl;  % Newtons
            % Note that we zero the mechanical force if the edge is fully within a rigid body.  
            
            % calculate potential energy due to edge strain
            e_mf_mag2 = e_k .* e_dl.^2 .* (e_notbroken);  
            PE = 5e-4 * sum(e_mf_mag2);  % should be in Joules, I think...
            
            % Calculate temperature difference across all edges due to different temperatures on every node
            e_dt = n_t(e_nb) - n_t(e_na);   % kelvin
            
            % Thermal conduction, i.e. heat flow power on edges ('e_hf') according to formula kappa * (A/l) * (T_hot - T_cold)
            e_hf = e_heatCondc .* e_dt;  % Watts
            
            % Terminate script if we've picked up a math error somehow, so we can debug before it corrupts the mesh.
            if isnan(sum(e_hf))
                error('Found NaN!  (in heatflow calc)')
            end
            
            % Keep track of time spent doing edge-based calculations
            timers(5) = timers(5) + toc;
            tic;
        
        
            %% 9. Assign mechanical forces to nodes (from edges) and also heat flow
            %if ~rigid_mode   %we can skip the rigid body mask if the entire mesh is flexible...
                n_mf(1,:) = e_mf(1,:) * M_e2n;
                n_mf(2,:) = e_mf(2,:) * M_e2n;
                n_mf(3,:) = e_mf(3,:) * M_e2n;
            %else  % Actually, e_mf will be zero for rigid edges now, so no need for if-statement here.. 
                

            % Calculate node heat flow from edge heat flow
            n_hf = e_hf * M_e2n;

            %DELETE ME:  This is for debugging: 
            t_cnd = n_hf(t_na) ./ n_nt(t_na) + n_hf(t_nb) ./ n_nt(t_nb) + n_hf(t_nc) ./ n_nt(t_nc);
            
            % keep track of time spent doing edge-to-node calculations
            timers(6) = timers(6) + toc;
            tic;

        
            %% 10. Evaluate for tensile or thermal failure
            % DISCLAIMER:  This simulator does not attempt to accurately model the evolution of tensile failure, other
            % than to detect the onset of such failure whenever the strain on a mesh edge exceeds the tensile limit of
            % the material.  It would be prudent to immediately terminate the simulation at the first such failure...
            %
            % but...
            %
            % Letting the simulator continue to run after tensile failure turns out to produce fun videos of the
            % lightsail tearing itself apart and flailing around in pieces!  So, to keep ourselves entertained, we can
            % crudely (and with no intent or claim of accuracy) continue to simulate the partially failed lightsail, by
            % essentially omitting the failed edges and triangles from optical and mechanical calculations, as if they
            % weren't there in the first place.  We denote such failed elements as 'broken'.  Here's where we detect
            % and keep track of broken edges and triangles:
            if any_flex  % we can skip tensile failure analysis if the entire mesh is rigid...
    
                e_tensilefailure =  ( e_flexible &  e_notbroken) & ( e_s > e_tensile_strain_limit )  ;  %all edges that are now broken from tensile failure
                e_thermalfailure = ( e_notbroken .* (e_t > e_failtemp ) );
                newtensilefailures = sum(e_tensilefailure);
                newthermalfailures = sum(e_thermalfailure);
                e_broken = e_broken | (e_tensilefailure > 0) | e_thermalfailure;
                %e_broken = e_nowbroken;
                
                % Just for fun, let's take the potential energy stored in recently severed edges, and dump it into
                % kinetic (recoil) energy of the adjoined nodes!  This is new as of v16, and clearly I'm off task again,
                % but some men just want to watch the world burn.
                % Potential energy due to edge strain = 5e-4 * e_mf_mag2.   One half of PE will go to each adjoining node.
                % From the edge's frame of reference, the adjoining nodes have approximately zero kinetic energy, so we
                %  can just use 1/2PE = deltaKE=1/2 m deltaV^2  -->  deltaV = sqrt( PE / m) to determine the velocity 
                %  delta to add to the adjoining nodes.  We'll use 1/2 e_m to get the approximate mass of the adjoining 
                %  nodes, which we can't access directly, so this will assume both nodes have equal mass.  
                %  We will subtract the velocity deltas in the distribution expressions, since the M_e2n matrix is set 
                %  up for tensile (inward) forces, and I'm pretty sure the recoil impulse should occur in the outward direction.
%                 if any(e_newbroken)
%                     
%                     n_vx = n_vx - (e_newbroken.*sqrt(2.5e+2*e_mf_mag2./e_m).*e_nl(1,:)) * M_e2n;
%                     n_vy = n_vy - (e_newbroken.*sqrt(2.5e+2*e_mf_mag2./e_m).*e_nl(2,:)) * M_e2n;
%                     n_vz = n_vz - (e_newbroken.*sqrt(2.5e+2*e_mf_mag2./e_m).*e_nl(3,:)) * M_e2n;
%                     
%                     %debug
%                     disp(['Edge failure (' num2str(sum(e_newbroken)) '), imparted velocities:  ' num2str(sqrt(2.5e+2*e_mf_mag2(e_newbroken)./e_m(e_newbroken))) ]);
%                 end
                % Edit:  Nope, this did not have the desired entertainment effect.  It is stupid and pointless.  Will
                % delete in next version.
                
                % Count number of broken bonds/edges and call it 'uhohs'
                tensilefailures = tensilefailures + newtensilefailures;
                thermalfailures = thermalfailures + newthermalfailures;
                
                uhohs = sum(e_broken);
                
                % Determine which triangles are broken... we used to delete triangles if any single edge failed, ...
                t_broken =  ( e_broken(t_e1) | e_broken(t_e2) | e_broken(t_e3) );  %and that's still the best i could come up with
                
                % Convenience variables so we don't have to perform negation every time we want the unbroken edges/tris
                e_notbroken = ~e_broken;
                t_notbroken = ~t_broken;
                
                % 'any' command returns '1' if any of the vector elements is nonzero.
                % 'anybroken' is a number that is '1' if there is at least one broken
                % bond in the mesh
                anybroken = anybroken || (any(t_broken));  %whether or not there has ever been a failure
                
                % Record time and time-index if the first breakage has just occurred
                if (anybroken) && (nto_firstBroken == sim_dur_nto)
                    nto_firstBroken = nto;
                    time_frstBrkn = tt;
                    
                    if terminate_on_breakup
                        if ~terminate_t
                            terminate_t = tt + termination_delay_s;
                            if ~terminate_t
                                terminate_t = dt;  %just in case the first breakage would try to schedule termination at exactly tt=0, we'll change it to a non-zero value, because zero is a sentinel value that prevents termination...
                            end
                        end
                    end
                end
                
                % keep track of time spent on edge breakage detection and bookkeeping
                timers(7) = timers(7) + toc;
                tic;
            end
            
            %% 11. Mechanical calculations and time stepping
            
            % We calculate external forces (not including elastic mesh forces) as an intermediate step
            n_rfx = n_of(1,:) + n_af(1,:);
            n_rfy = n_of(2,:) + n_af(2,:);
            n_rfz = n_of(3,:) + n_af(3,:);
            externalForceX = sum( n_rfx );
            externalForceY = sum( n_rfy );
            externalForceZ = sum( n_rfz );
            externalForce = [ externalForceX; externalForceY; externalForceZ];
            
            %We can calculate an effective torque on the whole mesh (as if it were rigid).  This is not actually used in
            %subsequent calculations...
            externalTorqueX = sum( n_dy.*n_rfz - n_dz.*n_rfy );
            externalTorqueY = sum( n_dz.*n_rfx - n_dx.*n_rfz );
            externalTorqueZ = sum( n_dx.*n_rfy - n_dy.*n_rfx );
            externalTorque = [ externalTorqueX; externalTorqueY; externalTorqueZ ];
            
            % sum mechanical and optical forces on each node
            n_fx = n_mf(1,:) + n_rfx;
            n_fy = n_mf(2,:) + n_rfy;
            n_fz = n_mf(3,:) + n_rfz;
            % Note that n_mf doesn't include edge spring contributions from within rigid bodies
            
            if rigid_mode  % include rigid body dynamics
                
                for nrb=1:num_rb
                    bodymask = n_rb == nrb;
                    
                    rb_force(1,nrb) = sum(n_fx(bodymask));
                    rb_force(2,nrb) = sum(n_fy(bodymask));
                    rb_force(3,nrb) = sum(n_fz(bodymask));
                    
                    rb_torque(1,nrb) = sum( n_rb_dy(bodymask).*n_fz(bodymask) - n_rb_dz(bodymask).*n_fy(bodymask) );
                    rb_torque(2,nrb) = sum( n_rb_dz(bodymask).*n_fx(bodymask) - n_rb_dx(bodymask).*n_fz(bodymask) );
                    rb_torque(3,nrb) = sum( n_rb_dx(bodymask).*n_fy(bodymask) - n_rb_dy(bodymask).*n_fx(bodymask) );
                    
                    rb_angMom(:,nrb) = rb_angMom(:,nrb) + 1e6.* rb_torque(:,nrb) .* dt;
                    rb_angVel(:,nrb) = inv(rb_inertiaT(:,:,nrb)) * rb_angMom(:,nrb);
%                 %if switching from flexible to rigid dynamics, angMomRigid will be 0.  If that's the case, we'll use the
%                 % approximate angular momentum we calculated from the mesh as our starting rigid angular momentum.
%                 % Otherwise, we'll use the previous value of angMomRigid.  Note:  my calculation of angular velocity is
%                 % probably wrong, this is pretty sketch.
%                 if (all(angMomRigid == [0; 0; 0]) )
%                     angMomRigid = inertiaT * angVel;
%                 end
                
                % iterate the angular velocity by calculating angular momentum using the inertia tensor, applying the
                % impulse of torque over the time step interval, then calculating the new angular velocity via the
                % inverse inertia tensor:
                %angMom = inertiaT * angVelRigid;
                %fprintf('%-5d ', nt);
                %fprintf('%g  ', angMomRigid);
               % angMomRigid = angMomRigid + 1e6.* externalTorque .* dt;
               % newAngVel = inv(inertiaT) * angMomRigid;
%                 fprintf('  |  ');
%                 fprintf('%g  ', angMomRigid);
%                 fprintf('  |  ');
%                 fprintf('%g  ', rigidTorque);
%                 fprintf('  |  ');
%                 fprintf('%g  ', newAngVel);
%                 fprintf('\n');
                %delAngVel = newAngVel - angVelRigid;
                %angVelRigid = newAngVel;
                
                % iterate the velocity and COM position
                rb_vel(:,nrb) = rb_vel(:,nrb) + 1e6 .* rb_force(:,nrb) ./ rb_mass(nrb) .* dt;
                rb_COM(:,nrb) = rb_COM(:,nrb) + dt.*rb_vel(:,nrb);
                
                dangle = rb_angVel(:,nrb) .* dt; % the angle of the dangle?
                
                ihaterotationmatrices = ...% gpuArray( ...   % does this need to be promoted to gpuArray?
                    [ cos(dangle(3)) -sin(dangle(3))  0; ...
                    sin(dangle(3))  cos(dangle(3))  0; ...
                    0               0               1 ] * ...
                    ...
                    [ cos(dangle(2))  0               sin(dangle(2)) ; ...
                    0               1               0 ; ...
                    -sin(dangle(2)) 0               cos(dangle(2)) ] * ...
                    ...
                    [ 1               0               0 ; ...
                    0               cos(dangle(1)) -sin(dangle(1)) ; ...
                    0               sin(dangle(1))  cos(dangle(1)) ] ;  %) ;
                
                n_rb_dxyz_new = ihaterotationmatrices * [ n_rb_dx(bodymask); n_rb_dy(bodymask); n_rb_dz(bodymask) ];
                
                n_rb_dx(bodymask) = n_rb_dxyz_new(1,:);  %save these so we don't have to recalculate rb_dxyz later
                n_rb_dy(bodymask) = n_rb_dxyz_new(2,:);
                n_rb_dz(bodymask) = n_rb_dxyz_new(3,:);
                
                %Now we calculate the new node relative velocities due to rotational motion.  These will be added to the
                %body COM velocity to set the new node velocity shortly
                n_vx_rot = ( rb_angVel(2,nrb).*n_rb_dxyz_new(3,:) - rb_angVel(3,nrb).*n_rb_dxyz_new(2,:) );
                n_vy_rot = -( rb_angVel(1,nrb).*n_rb_dxyz_new(3,:) - rb_angVel(3,nrb).*n_rb_dxyz_new(1,:) );
                n_vz_rot = ( rb_angVel(1,nrb).*n_rb_dxyz_new(2,:) - rb_angVel(2,nrb).*n_rb_dxyz_new(1,:) );
                
                % calculate the new node absolute positions for this body
                
                
                n_x_new = rb_COM(1,nrb) + n_rb_dxyz_new(1,:);
                n_y_new = rb_COM(2,nrb) + n_rb_dxyz_new(2,:);
                n_z_new = rb_COM(3,nrb) + n_rb_dxyz_new(3,:);
                
                n_vx_new = rb_vel(1,nrb) + n_vx_rot;
                n_vy_new = rb_vel(2,nrb) + n_vy_rot;
                n_vz_new = rb_vel(3,nrb) + n_vz_rot;
                

                n_x(bodymask) = n_x_new;
                n_y(bodymask) = n_y_new;
                n_z(bodymask) = n_z_new;
                    
                n_vx(bodymask) = n_vx_new;
                n_vy(bodymask) = n_vy_new;
                n_vz(bodymask) = n_vz_new;
                end
                    
            end
            if any_flex % flexible dynamics
                
                
                
                % iterate the node velocities
                n_vx(n_flexible) = n_vx(n_flexible) + 1e6 .* n_fx(n_flexible)./n_m(n_flexible) .* dt;
                n_vy(n_flexible) = n_vy(n_flexible) + 1e6 .* n_fy(n_flexible)./n_m(n_flexible) .* dt;
                n_vz(n_flexible) = n_vz(n_flexible) + 1e6 .* n_fz(n_flexible)./n_m(n_flexible) .* dt;
                
                % iterate the node positions
                n_x(n_flexible) = n_x(n_flexible) + n_vx(n_flexible) .* dt;
                n_y(n_flexible) = n_y(n_flexible) + n_vy(n_flexible) .* dt;
                n_z(n_flexible) = n_z(n_flexible) + n_vz(n_flexible) .* dt;
                
            end
            
            % update overall COM:
            coms_x = n_x .* n_m;
            coms_y = n_y .* n_m;
            coms_z = n_z .* n_m;
            
            comx = sum(coms_x)./totalmass_nodes ;
            comy = sum(coms_y)./totalmass_nodes ;
            comz = sum(coms_z)./totalmass_nodes ;
            
                                    
            %% 12. Advance time step
            tt = tt + dt;
            
            %% Apply initial sail offest and tilts at end of spin-up
            if (tt >= 0) && ( (tt-dt < 0) ) % note: tt has already been incremented at this point...
                disp('Doing initial offset stuff');
                
                if t0t_changeSpinAngle  % this just "resets" the spinning sail to a particular point in its rotation, for consistency between sims
                    approxSpinAng = atan2(n_y(xradidx)-comy, n_x(xradidx)-comx);
                    t0ta = t0t_s0 - approxSpinAng;
                    
                    n_x_new = cos(t0ta).*n_x - sin(t0ta).*n_y;
                    n_y = sin(t0ta).*n_x + cos(t0ta).*n_y;
                    n_x = n_x_new;
                    
                    n_vx_new = cos(t0ta).*n_vx - sin(t0ta).*n_vy;
                    n_vy = sin(t0ta).*n_vx + cos(t0ta).*n_vy;
                    n_vx = n_vx_new;
                    
                    for nrb = 1:num_rb
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
                 %   rb_COM(1,:) = rb_COM(1,:) + t0t_x0;
                end
                if t0t_y0 ~= 0
                    n_y = n_y + t0t_y0;
                    comy = comy + t0t_y0;
                 %   rb_COM(2,:) = rb_COM(2,:) + t0t_x0;
                end
                
                
                
            end

            
            %% 13. Update temperatures
            % Heat capacity links the temperature rise of a body to the required
            % amount of energy. The unit of heat capacity is [J/g/degC], which
            % describes the amount of energy required to raise the temperature of 1
            % gram of an object by 1�C (or 1 Kelvin). Energy is equal to power
            % times time, such that the temperature rise is given by energy / mass
            % / heat capacity, where energy = (input power - output power) * dt.
            % Input power is given by the opticaly absorbed power and power from
            % heat flow/thermal conduction, whereas output power is due to
            % radiative cooling.
            
            n_t = n_t + (n_hf + n_ht_flux).*dt ./ n_heatmass;  
            
            % Calculate the edge temperature as the averaged temperature of the sum
            % of the temperatures of the two nodes that form the edge
            e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;
            
            % Calculate temperature of triangle as the average of the temperatures
            % of the nodes that make up the triangle
            t_t  =  ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) ./ 3;
            
            
            %% 14.  Update total body velocity, KE, spin speed...
            % Calculate kinetic energy as KE = (1/2) * m * v^2, which lags behind
            % during propulsion by one time step, i.e. dt, due to velocity Verlet
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
            

            %% 15:  Update inertia tensors for rigid bodies
            %calculate the local positions of nodes (relative to COM) and the
            %inertia tensor matrix
            if rigid_mode
                for nrb=1:num_rb
                    bodymask = n_rb == nrb;
                    %I can skip calculating COM from the node positions, because the rigid body COM was already used to
                    %set the node positions above...
%                     rb_COM(1,nrb) = sum(coms_x(bodymask))./rb_mass(nrb) ;
%                     rb_COM(2,nrb) = sum(coms_y(bodymask))./rb_mass(nrb) ;
%                     rb_COM(3,nrb) = sum(coms_z(bodymask))./rb_mass(nrb) ;
        
                    % I can skip calculating relative positions too, since I retained them during the node position
                    % update above...
%                     temp_dx = n_x(bodymask) - rb_COM(1,nrb);
%                     temp_dy = n_y(bodymask) - rb_COM(2,nrb);
%                     temp_dz = n_z(bodymask) - rb_COM(3,nrb);
    
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

%         n_rb_dx(bodymask) = temp_dx;  % save these because we'll need them for torque calculation
%         n_rb_dy(bodymask) = temp_dy;
%         n_rb_dz(bodymask) = temp_dz;
        
            
                end
            end
            
            
            %% This is a crappy attempt to calcualate overall angular velocity of the mesh, it doesn't really work and should probably 
            % be deleted...
            
            n_dx = n_x - comx;
            n_dy = n_y - comy;
            n_dz = n_z - comz;
            n_dxyz = [n_dx; n_dy; n_dz];
            
            n_dx2 = n_dx.^2;
            n_dy2 = n_dy.^2;
            n_dz2 = n_dz.^2;
            n_dr = sqrt( n_dx2 + n_dy2 + n_dz2);

            %rotational position and angular velocity
            %approxSpinAng = ...  moved this to monitor output, don't need it each loop
            
            %make sure we don'r risk dividing by zero when calculating the body's angular velocity:
            fpmsucksthreshold = 0.00001*radiusmm;
            fpmsucksmask = ( n_dx2 > fpmsucksthreshold ) & ( n_dy2 > fpmsucksthreshold ) & ( n_dr > fpmsucksthreshold );
            
            n_dxyzav = [n_dx(fpmsucksmask); n_dy(fpmsucksmask); n_dz(fpmsucksmask)];
            n_dxav2 =  n_dx2(fpmsucksmask);
            n_dyav2 =  n_dy2(fpmsucksmask);
            n_dzav2 =  n_dz2(fpmsucksmask);
            n_drav = n_dr(fpmsucksmask);
            
            n_vav = [ n_vx(fpmsucksmask)-meanvx; n_vy(fpmsucksmask)-meanvy; n_vz(fpmsucksmask)-meanvz ];
            angVelTemp = cross( n_dxyzav, n_vav );
            angVelTemp(1,:) = angVelTemp(1,:) ./ (n_dzav2 + n_dyav2 );
            angVelTemp(2,:) = angVelTemp(2,:) ./ (n_dxav2 + n_dzav2 );
            angVelTemp(3,:) = angVelTemp(3,:) ./ (n_dxav2 + n_dyav2 ); %./ n_drav.^2;%
            
            angVel = mean(angVelTemp,2);
            
            % inertia tensor
%             if rigid_mode
%                 inertiaT = [  ...
%                     sum( n_m .* (n_dy2 + n_dz2) )        ...
%                     sum( -n_m .* n_dx .* n_dy)         ...
%                     sum( -n_m .* n_dx .* n_dz );  ...
%                     sum( -n_m .* n_dy .* n_dx )            ...
%                     sum( n_m .* (n_dx2 + n_dz2) )      ...
%                     sum( -n_m .* n_dy .* n_dz );  ...
%                     sum( -n_m .* n_dz .* n_dx )            ...
%                     sum( -n_m .* n_dz .* n_dy )        ...
%                     sum( n_m .* (n_dx2 + n_dy2) ) ]';
%             else
%                 angVelRigid = [0; 0; 0 ];
%             end
            
            % Calculate triangle centroids, where the x, y or z coordinate of a
            % centroid is calculated as the arithmetic mean of the respective x, y
            % or z vertex/nodes coordinates
            t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
            t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
            t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;
            
            
%             %% 15.  FFT analysis -- record monitor data.
%             if fft_mode
%                 n_xrelcom = n_dx;
%                 n_yrelcom = n_dy;
%                 n_rrelcom = sqrt(n_xrelcom.^2+n_yrelcom.^2);  % this is NOT a duplicate of n_dr-- this one is projected in xy plane
%                 n_xrelcom = n_xrelcom./n_rrelcom;
%                 n_yrelcom = n_yrelcom./n_rrelcom;
%                 % n_tan_norm = [ -n_yrelcom ; n_xrelcom ];
%                 % projet a onto b:  a dot bnorm
%                 %n_rmf = n_mf(1,:).*n_xrelcom + n_mf(2,:).*n_yrelcom;
%                 %n_tmf = n_mf(1,:).*(-n_yrelcom) + n_mf(2,:).*n_xrelcom;
%                 %n_rmf = n_pec(1,:).*n_xrelcom + n_pec(2,:).*n_yrelcom;
%                 %n_tmf = n_pec(1,:).*(-n_yrelcom) + n_pec(2,:).*n_xrelcom;
%                 n_rmf = n_mf(1,:).*n_xrelcom + n_mf(2,:).*n_yrelcom;
%                 n_tmf = n_mf(1,:).*(-n_yrelcom) + n_mf(2,:).*n_xrelcom;
%             end
            
            timers(8) = timers(8) + toc;
            tic;
            
            
          %  gather(g_n_x(1));  %this should force us to wait until all GPU calculations have terminated...  not sure how far back things might be queued
            
            timers(15) = timers(15) + toc;
            tic;
            
            %% 16.  Record standard output monitors (downsampled)  
            % Note that video/snapshot outputs occur only during monitor output frames, to avoid rendudant gather() steps.  
            if ( ~mod(nt-1, sim_output_downsample) )  %TODO:  let's use nt-1 so we catch the first loop output?
                
                nto = nto + 1;
                
                if ReverseRayTrace && usetether
                   %RRT_tri_ctrs = [n_x(t_na) + n_x(t_nb) + n_x(t_nc); ...
                   %                n_y(t_na) + n_y(t_nb) + n_z(t_nc); ...
                   %               n_z(t_na) + n_z(t_nb) + n_z(t_nc)] ./ 3;
                   
                   RRT_idirs = [t_cx; t_cy; t_cz] - [ n_x(end); n_y(end); n_z(end) ];
                   RRT_idirs = RRT_idirs ./ vecnorm(RRT_idirs);
                   RRT_costhetas = dot(t_norms, RRT_idirs);
                   RRT_odirs = RRT_idirs - 2 .* RRT_costhetas .* t_norms; 
                   
                   if (max(abs( vecnorm(RRT_odirs) - 1 ) > 0.00001 ))
                       error('You screwed up.');
                   end
                   
                   RRT_ax(:, nto) = rad2deg(asin(RRT_odirs(1,:)));
                   RRT_ay(:, nto) = rad2deg(asin(RRT_odirs(2,:)));
                   RRT_a(:,nto) = rad2deg(asin(sqrt( RRT_odirs(1,:).^2 + RRT_odirs(2,:).^2 )));
                    
                end
                
                
                
                %calculate some values here that don't need to be calculated during the loop:
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
%                 if (rigid_mode)
%                     angVels(:,nto) = angVelsRigid(:,nto);
%                 else
                    angVels(:,nto) = (angVel);
%                 end
                
                torques(:,nto) = (externalTorque);
                spinAngs(nto) = (atan2(n_dy(xradidx), n_dx(xradidx)));  %disp([num2str(rad2deg(approxSpinAng)) ' sam_wait=' num2str(sam_wait) ] );
                sailNorms(:,nto) = (mean(t_norms(:,t_centerTris),2));  % this is an approximation of the sail tilt direction -- the average norms of the center triangles comprising the sail
                sailPitches(:,nto) = atan2( sqrt(sailNorms(1,nto).^2 + sailNorms(2,nto).^2), sailNorms(3,nto) ); % this is the tilt angle in degrees(from vertical) of the approximate sail norm
                
                ofxs(nto) = externalForceX;
                ofys(nto) = externalForceY;
                ofzs(nto) = externalForceZ;
                
                nofxs(nto) = ofxs(nto) / (I0 * ((2*radiusmm).^2) / c0);
                nofys(nto) = ofys(nto) / (I0 * ((2*radiusmm).^2) / c0);
                nofzs(nto) = ofzs(nto) / (I0 * ((2*radiusmm).^2) / c0);
                
                tqxs(nto) = externalTorqueX;
                tqys(nto) = externalTorqueY;
                tqzs(nto) = externalTorqueZ;
%                 tqxs(nto) = rigidTorqueX;
%                 tqys(nto) = rigidTorqueY;
%                 tqzs(nto) = rigidTorqueZ;
                ntqxs(nto) = tqxs(nto) / (I0 * ((2*radiusmm).^3) / c0);
                ntqys(nto) = tqys(nto) / (I0 * ((2*radiusmm).^3) / c0);
                ntqzs(nto) = tqzs(nto) / (I0 * ((2*radiusmm).^3) / c0);
                graphuos(nto) = ( uhohs );
                graphtensfails(nto) = tensilefailures;
                graphthermfails(nto) = thermalfailures;
                graphts(nto) = tt;
                I0s(nto) = pwrRampVal;
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
                ttlradpwr(nto) = (sum(t_ems));
                areas(nto) = (sum(t_a.*(t_notbroken)));
                areasxy(nto) = (sum(t_aproj.*(t_notbroken)));
                upsidedowntris(nto) = (upsidedowntris);
                areasrightsideupxy(nto) = (sum( t_aproj.*(t_aproj>0).*(t_notbroken)));
                mintrinormz(nto) = (min( t_norms(3,:) ));
                % monitorqty(nto) = mean(n_rps(2:end));
                xradpttemps(nto) = n_t(end-1);
                vertextemps(nto) = n_t(1);
                
%skipping these due to disuse:
%                 xradptxyzs(:,nto) = [n_x(xradidx); n_y(xradidx); n_z(xradidx)];
%                 yradptxyzs(:,nto) = [n_x(yradidx); n_y(yradidx); n_z(yradidx)];
%                 nxradptxyzs(:,nto) = [n_x(nxradidx); n_y(nxradidx); n_z(nxradidx)];
%                 nyradptxyzs(:,nto) = [n_x(nyradidx); n_y(nyradidx); n_z(nyradidx)];
%                 %monitorfxs(:,nto) = n_mf(1,:);
%                 %monitorfys(:,nto) = n_mf(2,:);
                if fft_mode && tt >= fft_start_time  && tt <= (fft_start_time + fft_dur) && fft_num < fft_depth
                    fft_num = fft_num+1;
                    monitorfzs(:,fft_num) = n_mf(3,:);
                    monitorrfs(:,fft_num) = n_rmf;
                    monitortfs(:,fft_num) = n_tmf;
                    monitorxs(:,fft_num) = n_x;
                    monitorys(:,fft_num) = n_y;
                    monitorzs(:,fft_num) = n_z;
                end
                if save_mesh_mode && savemesh_num < savemeshdepth
                    if (tt - savemesh_last) > save_mesh_dt  
                        savemesh_last = tt;
                        savemesh_num = savemesh_num+1;
                                    
                        savemeshn_x(:,savemesh_num) = n_x;
                        savemeshn_y(:,savemesh_num) = n_y; 
                        savemeshn_z(:,savemesh_num) = n_z; 
                        %savemesht_t(:,savemesh_num) = t_t; 
                        savemeshn_t(:, savemesh_num) = n_t;
                        savemeshn_hf(:,savemesh_num) = n_hf;
                        savemesht_abs(:, savemesh_num) = t_abs;
                        savemesht_abs_mr(:, savemesh_num) = t_abs_mr;
                        savemesht_ems(:,savemesh_num) = t_ems;
                        savemesht_rht(:, savemesh_num) = t_rht;
                        %savemesht_cnd(:, savemesh_num) = n_hf(t_na) ./ n_nt(t_na) + ...
                        %                                 n_hf(t_nb) ./ n_nt(t_nb) + ...
                        %                                 n_hf(t_nc) ./ n_nt(t_nc);
                        savemesht_a(:, savemesh_num) = t_a; % note we can use this to display strain, usually we do sqrt(ta./t_a0)-1
                        savemesht_maxstr(:,savemesh_num) = max( max( e_s(t_e1), e_s(t_e2) ), e_s(t_e3) );  % alternate view of strain
                        savemesh_ntos(savemesh_num) = nto;  % for convenience of synchronizing index to standard output data
                    end
                end
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

                                
                
                    
                %% PLOT VIDEO OUTPUT FRAME
                
                if (nt - pf >= frame_interval - 1)
                    pf = nt + 1;  % time step counter of prior frame
                    nf = nf + 1;  % number of movie frames
                    
                    %       [~,idx_max_theta] = max(abs(pitch_angles312));
                    %        [~,idx_max_phi] = max(abs(roll_angles312));
                    %       [~,idx_max_psi] = max(yaw_angles312);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
                    if (RayTracingMode)
                        pause(0.5);  % this is stupid and terrible, my computer is overheating and I'm too lazy to fix it.
                    end
                    timers(10) = timers(10) + toc;
                    tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
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
                            %try
                                if (gcf ~= fSimVid)
                                    figure(fSimVid);
                                end
                                plotMesh_v17;
                                writeVideo(V1,getframe(gcf));
                            %catch the_last_few_minutes_of_an_interesting_tv_show
%                                     warning('Something went wrong with the video recorder!  No longer recording video output.');
%                                     videoMode = 2;
%                                     try
%                                         close(V1);
%                                         close(fSimVid);
%                                     catch light_of_the_fact_that_youve_failed_once_again
%                                         disp('Fuck it dude, let''s go bowling...');
%                                     end
                            %end
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

                if ((nt-ps >= snapshot_intvl-1 ) && (ns < 10) ) || (graphuos(nto) >= nextExplosionSnapshot)
                    ps=nt+1;
                    ns=ns+1;

                    tscale_max = max(tscale_max, maxts(nto));
                    tscale_min = min(tscale_min, maxts(nto));
                    if (graphuos(nto) >= nextExplosionSnapshot)
                        nextExplosionSnapshot = nextExplosionSnapshot * explosionSnapshotRatio;
                    end
                    figure(fSnapshot);
                    plotSnapshot;
                end

                % keep track of time spent plotting snapshots
                timers(12) = timers(12) + toc;
                tic;
                
            end
            
            
            %% 19. check for user interaction with dialog boxes
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
%                     
%                 
                
                
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
        
        %% display simulation timing statistics
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
        
        
        %% complete video recording, close video window
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
        
        
        %% follow-up calculations and data handling
        nto_firstBroken = min(nto_firstBroken,nto+1);
        hoop12xyz = hoop12xyz(:,1:nto);
        hoop6xyz = hoop6xyz(:,1:nto);
        hoop3xyz = hoop3xyz(:,1:nto);
        hoop9xyz = hoop9xyz(:,1:nto);

        longaxis_length = vecnorm(hoop12xyz-hoop6xyz);
        shortaxis_length = vecnorm(hoop3xyz-hoop9xyz);
        
        nt=nto;   % change meaning of nt, now it means the number of monitor frames.
        plot_time = graphts(1:nt);
        
        
        if fft_mode
            if (fft_depth > fft_num)
                monitorrfs = monitorrfs(:,1:fft_num);
                monitortfs = monitortfs(:,1:fft_num);
                %monitorfxs = monitorfxs(:,1:nto);
                %monitorfys = monitorfys(:,1:nto);
                monitorfzs = monitorfzs(:,1:fft_num);  %trim these down now...
                monitorxs = monitorxs(:,1:fft_num);
                monitorys = monitorys(:,1:fft_num);
                monitorzs = monitorzs(:,1:fft_num);
            else
                warning('Monitor buffer overflow/wrap!!! fix me!!!');
            end
        end
        if save_mesh_mode 
            if savemesh_num < savemeshdepth
        
                savemeshn_x(:,savemesh_num+1:end) = [];
                savemeshn_y(:,savemesh_num+1:end) = [];
                savemeshn_z(:,savemesh_num+1:end) = [];
%                savemesht_t(:,savemesh_num+1:end) = [];
                savemesht_abs(:,savemesh_num+1:end) = [];
                savemesht_abs_mr(:,savemesh_num+1:end) = [];
                savemesht_ems(:,savemesh_num+1:end) = [];
                savemesht_rht(:,savemesh_num+1:end) = [];

            %    savemesht_cnd(:,savemesh_num+1:end) = [];
                savemesht_a(:,savemesh_num+1:end) = [];
                savemesht_maxstr(:,savemesh_num+1:end) = [];
            end
        end
        
        
        %% prepare / save main simulation output window
        disp('Generating output plot.');
        hOutputPlot = figure('pos',[100 100 2100 1200]);
        set(gcf,'color','w');
        
        hafp1 = subplot(4,5,1);
        velxs = velxs(1:nt);
        plot(plot_time, velxs);
        title('X & Y velocity (mm/s)');
        xlabel('Time');
        hold on
        velys = velys(1:nt);
        plot(plot_time, velys);
        legend({'X' 'Y'});
        
        
        hafp2 = subplot(4,5,2);
        velzs = velzs(1:nt);
        plot(plot_time, velzs);
        title('Z velocity (mm/s)');
        
        hafp3 = subplot(4,5,3);
        areas = areas(1:nt);
        plot(plot_time, areas);
        title('Area (mm^2)');
        hold on
        areasxy = areasxy(1:nt);
        plot(plot_time, areasxy);
        areasrightsideupxy = areasrightsideupxy(1:nt);
        plot(plot_time, areasrightsideupxy, '--');
        Ylims = get(gca,'YLim');
        set(gca,'Ylim',[ Ylims(1) Ylims(2)+0.3*(Ylims(2)-Ylims(1)) ] );
        legend({'Surface area', 'Projected area', 'Upright projected' } );
        
        
        hafp4 = subplot(4,5,4);
        minstrains = minstrains(1:nt);
        plot(plot_time, minstrains);
        title('Strain min/max');
        hold on;
        maxstrains = maxstrains(1:nt);
        plot(plot_time, maxstrains(1:nt));
        plot(plot_time([1 end]), [0 0], 'k');
        
        hafp5 = subplot(4,5,5);
        graphuos = graphuos(1:nt);
        plot(plot_time, graphuos(1:nt));
        title('Number of tensile failures');
        set(gca,'YScale','log')
        
        hafp6 = subplot(4,5,6 );
        comxs = comxs(1:nt);
        plot(plot_time, comxs(1:nt)-((0)));
        title({'Position (mm)' '(COM rel to beam ctr)'});
        hold on;
        comys = comys(1:nt);
        plot(plot_time, comys(1:nt)-0);
        comrs = comrs(1:nt);
        plot(plot_time, comrs(1:nt), '--');
        legend({'X' 'Y' 'R'});
        
        
        hafp7 = subplot(4,5,7);
        comzs = comzs(1:nt);
        plot(plot_time, comzs(1:nt));
        title('Z position (mm)');
        
        hafp8 = subplot(4,5,8);
        kinoenergies = kinoenergies(1:nt);
        plot(plot_time, kinoenergies(1:nt));
        title('Non-ballistic KE (J)');
        
        hafp9 = subplot(4,5,9);
        potenergies = potenergies(1:nt);
        plot(plot_time, potenergies(1:nt));
        title('Potential energy (J)');
        
        
        %subplot(3,4,9);
        %plot(plot_time, maxvs(1:nt));
        %title('Max velocity rel. to COM');
        
        hafp10 = subplot(4,5,10);
        mintrinormz = mintrinormz(1:nt);
        plot(plot_time, rad2deg( asin(mintrinormz(1:nt))));
        title('Steepest triangle angle (\circ)');
        hold on
        plot(plot_time([1 end]), [0 0], 'k');
        
        
        
        hafp11 = subplot(4,5,11);
        plot(plot_time, ttlabspwr(1:nt));
        ttlabspwr = ttlabspwr(1:nt);
        title('Power absorbed / emitted');
        hold on;
        ttlradpwr = ttlradpwr(1:nt);
        plot(plot_time, ttlradpwr(1:nt));
        legend({ 'Abs', 'Ems' });
        
        hafp12 = subplot(4,5,12);
        acczs = acczs(1:nt);
        plot(plot_time, acczs(1:nt));
        hold on;
        accts = accts(1:nt);
        plot(plot_time, accts(1:nt));
        title('Acceleration (mm/s/s)');
        legend({'Z' 'Total'});
        
        hafp13 = subplot(4,5,13);
        avgts = avgts(1:nt);
        plot(plot_time, avgts(1:nt));
        hold on
        vertextemps = vertextemps(1:nt);
        plot(plot_time, vertextemps(1:nt));
        xradpttemps = xradpttemps(1:nt);
        plot(plot_time, xradpttemps(1:nt));
        maxts = maxts(1:nt);
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
        sailPitches = sailPitches(1:nt);
        plot(plot_time, rad2deg(sailPitches(1:nt)) );
        hold on;
        sailNorms = sailNorms(:,1:nt);
        plot(plot_time, rad2deg(atan2( sailNorms(2,1:nt), sailNorms(3,nt) )) )
        plot(plot_time, -rad2deg(atan2( sailNorms(1,1:nt), sailNorms(3,nt) )) )
        title('Sail tilt (\circ)' );
        legend({ 'total off Z'  'X' 'Y'});
        
        hafp17 = subplot(4,5,17);
        angVels = angVels(:,1:nt);
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
        ofxs = ofxs(1:nt);
        plot(plot_time, ofxs(1:nt) );
        hold on;
        ofys = ofys(1:nt);
        plot(plot_time, ofys(1:nt) );
        title('Lateral optical force' );
        legend({ 'X' 'Y' } );
        
        hafp20 = subplot(4,5,20);
        tqxs = tqxs(1:nt);
        plot(plot_time, tqxs(1:nt) );
        hold on;
        tqys = tqys(1:nt);
        tqzs = tqzs(1:nt);
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
        
        saveas(gcf,[filebasename '.fig']);
        saveas(gcf,[filebasename '.png']);
        
        %% save snapshot figure
        disp('Generating snapshot plot.');
        figure(fSnapshot);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        set(gca,'ZTick',[]);
        set(gca,'Box','on');
        title(filebasename, 'Interpreter', 'none');
        plot3([((0)) ((0))], [((0)) ((0))], get(gca,'ZLim'),'m');
        
        saveas(gcf,[filebasename '_snapshot.fig']);
        saveas(gcf,[filebasename '_snapshot.png']);
        
        %% generate / save trajectory figure
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
        saveas(gcf,[filebasename '_trajxy.fig']);
        saveas(gcf,[filebasename '_trajxy.png']);
        
        
        %% generate / save restoring force figure
%         if nto_startAccel && (nto_firstBroken > nto_startAccel)
%             disp('Generating restoring force analysis plot.');
%             plotRestoringForceAnalysis;
%             saveas(gcf,[filebasename '_forcetorque.fig']);
%             saveas(gcf,[filebasename '_forcetorque.png']);
%         end
        
        %% RRT stuff
        if (ReverseRayTrace && ~usercancel)
            RRT_thing;
            RRTfom = sqrt(x_fom1.^2+y_fom1.^2);
        else
            RRTfom = 0;
        end
        %% generate/display data table entry
        %        disp('Data table entry:')
        
        
        
        peakstrain = max(maxstrains);
        peakt = max(maxts);
        peakaccz = max(acczs);
        peakpwr = max(ttlinpwr);
        
        spacer = 0;  %just used to fill in missing columns in the table
        dispvars = { 'movieNum', 'RRTfom'...
            'radiusmm','zAspectRatio','totalarea','totalmass','sum(t_a0xy)','rot0', ...
            'numnodes','numedges','numtris',...
            'thickness','density','Ymod',...
            'Emod','Thermcondmm','heatCap',... %   'cte0','cte1','cts','ctw',  ...
            't0t_x0','t0t_y0','t0t_t0_x','t0t_t0_y' , ...
            'I0','Irad','spacer', ...
            'Iabs','Irefl','Emissivity', ...
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
        
        burstTestVals(end+1) = pwrRampVal * I0;
        
        
        
        %% save the whole workspace (except figures)
        disp('Saving workspace');
        allvars = whos;
        % Identify the variables that ARE NOT graphics handles. This uses a regular
        % expression on the class of each variable to check if it's a graphics object
        tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));
        % Pass these variable names to save
        save([filebasename '_workspace.mat'],  allvars(tosave).name);
        %  Note:  I don't need to worry about freeing unused buffer space before saving the workspace, because somehow
        %  MATLAB compresses the workspace when saving to disk.
        disp(['Workspace memory usage: ' num2str(sum([allvars(:).bytes])./1024./1024) ' MB']);
        
        %% wrap up, prepare for next simulation
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
