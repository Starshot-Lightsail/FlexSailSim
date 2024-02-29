clear
clc
close all


% Include MEPS?
includeMEPS = 0;

% Stop simulation if lasersail is displaced by more than this threshold
stopSimulation_Displacement = 1; % in radius*x0

% Outermost loop over different beam waist radii, if value is 1, then the
% beam waist radii is equal to the radius of the lightsail, if round
outersweepvals2 = [2.00]; 

% Outer loop over displacement of beam center from origin
outersweepvals1 = 1*0.2; %0.55:0.05:0.8;
outerloopstable = ''; % What does this variable represent? And what does it do?
for nouterloop2 = 1:length(outersweepvals2)
for nouterloop1 = 1:length(outersweepvals1)

    disp(['OuterLoop [' num2str(nouterloop2) ',' num2str(nouterloop1) '] '  num2str(outersweepvals2(nouterloop2)) ' ' num2str(outersweepvals1(nouterloop1)) ] );

% simulate
% close(V1);
movieNum = 1; % Iterator for simulation run?  
% movieNum = movieNum + 1; % Need to un-comment this if first data set is
%already created
diary(['sim' num2str(movieNum) '.diary']); % Log command window text to file

%% generate mesh & set up constants 
generateMesh_02092021_RG; % This is another script/function being called here

% Needed to decide when to switch from explicit Euler integration to
% velocity Verlet algorithm, and immediately after that, skipping the first
% time calculation of the velocities
first_switch_integrator = 1;

I0 = 1000; % W/mm^2, peak incident power density (1000 = 1GW/m2)
Irad = radius*x0*outersweepvals2(nouterloop2);  % mm, beam waist radius
Ictrx = x0*radius*outersweepvals1(nouterloop1);
spinup_gs = 1;  % Arbitrary units.  I'm not sure what the scale factor was here.  But it works out to 3181.7 RPS/spinup_seconds.  
rps_target = 80; % Rotations per second
spinup_time = rps_target/1801.8/spinup_gs; % Why 1801.8?
t_start_ramp = 1e-3;%4 / rps_target;%5e-3;%5e-3;  
t_ramp_dur = 5e-3;%1 / rps_target;%5e-3;%5e-3;
t_on_dur = 4; %10;
% sim_dur = 4 + t_start_ramp + t_ramp_dur;
sim_dur = t_start_ramp + t_ramp_dur + t_on_dur; % Should it be this, Mike?

sim_output_downsample = 8;
frame_interval_s = 5*1e-4;    
snapshot_interval_s = 0.01;%.250;
breakup_slowmo_interval_s = max( frame_interval_s/8, 5e-6);
kick_interval_s = 2e-6;
monitor_depth = floor(5e8/length(n_x)/2)*2;  % Number of full-grid monitor records to try to fit in memory.  

pwrRampFcn = @(tt) 0.5*double(tt > t_start_ramp).*double(tt <= t_start_ramp+t_ramp_dur).*(1-cos(pi*(tt-t_start_ramp)./(t_ramp_dur))) + double(tt > t_start_ramp+t_ramp_dur).*double(tt <= t_start_ramp+t_ramp_dur+t_on_dur) + ...
    0.5*double(tt > t_start_ramp+t_ramp_dur+t_on_dur).*double(tt <= t_start_ramp+t_ramp_dur+t_on_dur+t_ramp_dur).*(1+cos(pi*(tt-t_start_ramp-t_ramp_dur-t_on_dur)./t_ramp_dur));

% here are some important constants
c0 = 3e8; % m/sec, speed of light.  I am confused as to whether or not this should be converted to mm for my calculations.  
c0mm = c0*1000; % mm/sec, speed of light
% To get rad presure, use I0 / c0.  This gives Pa / mm2
SBC = 5.67e-8; % W/m^2/K^4  for radiative cooling.  
SBCmm = SBC / 1e6; % W/mm^2/K^4

%% Material constants, e.g. for silicon nitride or silicon

% Iabs = 1e-5; % absorptivity
Iabs = 1e-5; % absorptivity
Irefl = 0.15; % reflectivity
Emissivity = 0.1; % emissivity

%% More memory allocation
n_accx = zeros(1,numnodes); % x-coordinate of each node
n_accy = zeros(1,numnodes); % y-coordinate of each node
n_accz = zeros(1,numnodes); % z-coordinate of each node

%% more variables, settings, and setup
plotFast = 0;
monitor_mode = 0;
%dt = 4*5e-7;
dt = rec_dt;
frame_interval = max(floor(frame_interval_s/dt),1);  %number of loops between movie frames
snapshot_interval = max(floor(snapshot_interval_s/dt),1);  % number of loops between snapshots
kick_interval = max(floor(kick_interval_s/dt),1);
sim_dur_nt = ceil( (sim_dur + spinup_time) / dt );
sim_dur_nto = ceil(sim_dur_nt / sim_output_downsample);
dnto = sim_output_downsample * dt;
tt = -spinup_time;
uhohs = 0;
pf = 1;
nf = 1;
ps = 1;
ns = 1;
nextExplosionSnapshot = 10;  % force a snapshot when #uhohs reaches this value
explosionSnapshotRatio = 10;   %and thereafter increase the threshold by this ratio
tscale_min = t0-1;
tscale_max = t0+1;
sscale_max = 0.0002;
max_strain = 0;
min_strain = 0;

%% Import and set up simulated optical pressures from COMSOL simulation

% resizePressureLookUpTables_v3; % call Ramon's script

% Load full look-up table of pressures vs roll (phi) and pitch (theta)
% angles for TE- and TM-polarized light.
load(fullfile('Data','02012021_MEPS_SiNx_Mark4_ReducedPitchRoll_Pressures_TE_Table.mat'))
load(fullfile('Data','02012021_MEPS_SiNx_Mark4_ReducedPitchRoll_Pressures_TM_Table.mat'))

% Load splitted smaller look-up tables of pressures vs roll and pitch
% angles for TE- and TM-polarized light. All the small look-up tables for
% TE (TM) polarization combined together yield the full look-up table for
% TE (TM) polarization. The full look-up table is being splitted into
% smaller ones to increase speed of finding pressures within the table
% later in the code
for ii = 1:10                    
    load([fullfile('Data','02012021_MEPS_SiNx_Mark4_ReducedPitchRoll_Pressures_TE_Table') ...
        num2str(ii) '.mat'])
    load([fullfile('Data','02012021_MEPS_SiNx_Mark4_ReducedPitchRoll_Pressures_TM_Table') ...
        num2str(ii) '.mat'])
end

% Assign simulated roll and pitch angles to n x 2 matrix
phi_angles = pressure_TE_table(:,1);
phi_angles = deg2rad(phi_angles);
theta_angles = pressure_TE_table(:,2);
theta_angles = deg2rad(theta_angles);
rollpitch_angles = [phi_angles, theta_angles];

% Speed of light in m/s
c0 = 299792458;

% Assign pressures in local x and z direction to vectors
phi_angles_1 = unique(pressure_TE_table1(:,1));
phi_angles_2 = unique(pressure_TE_table2(:,1));
phi_angles_3 = unique(pressure_TE_table3(:,1));
phi_angles_4 = unique(pressure_TE_table4(:,1));
phi_angles_5 = unique(pressure_TE_table5(:,1));
phi_angles_6 = unique(pressure_TE_table6(:,1));
phi_angles_7 = unique(pressure_TE_table7(:,1));
phi_angles_8 = unique(pressure_TE_table8(:,1));
phi_angles_9 = unique(pressure_TE_table9(:,1));
phi_angles_10 = unique(pressure_TE_table10(:,1));

% Assign important boundary phi angles
coarse_first_angle = deg2rad(phi_angles_1(1));
coarse_last_angle = deg2rad(phi_angles_1(end));
coarse_angle_step = deg2rad(abs(phi_angles_1(2)-phi_angles_1(1)));

finer_first_angle = deg2rad(phi_angles_2(1));
finer_last_angle = deg2rad(phi_angles_2(end));
finer_angle_step = deg2rad(abs(phi_angles_2(2)-phi_angles_2(1)));

even_finer_first_angle = deg2rad(phi_angles_3(1));
even_finer_last_angle = deg2rad(phi_angles_3(end));
even_finer_angle_step = deg2rad(abs(phi_angles_3(2)-phi_angles_3(1)));

finest_first_angle = deg2rad(phi_angles_4(1));
finest_last_angle = -finest_first_angle;
finest_angle_step = deg2rad(abs(phi_angles_4(2)-phi_angles_4(1)));

finest_intermediate_angle1 = deg2rad(phi_angles_4(end));
finest_intermediate_angle2 = deg2rad(phi_angles_5(end));

% Get indices of which triangles correspond to which texture region
indices_Region0 = find(t_tex == 0);
indices_Region1 = find(t_tex == 1);
indices_Region2 = find(t_tex == 2);
indices_Region3 = find(t_tex == 3);
indices_Region4 = find(t_tex == 4);
indices_Region5 = find(t_tex == 5);
indices_Region6 = find(t_tex == 6);

% Pre-allocate memory for pressures & forces on mesh triangles
t_press_x = zeros(size(t_a));
t_press_y = zeros(size(t_a));
t_press_z = zeros(size(t_a));

t_force_x_BF = zeros(size(t_a));
t_force_y_BF = zeros(size(t_a));
t_force_z_BF = zeros(size(t_a));
t_force_x = zeros(size(t_a));
t_force_y = zeros(size(t_a));
t_force_z = zeros(size(t_a));

%% display settings (for log)
dispvar3('I0','Irad','Ictrx');
dispvar3('Iabs','Irefl','Emissivity');
dispvar3('spinup_gs','spinup_time','plotFast');
dispvar3('t_start_ramp','t_ramp_dur','t_on_dur');
dispvar3('dt','frame_interval','snapshot_interval');
dispvar2('sim_dur','sim_dur_nto');
if (sim_dur_nto > 1e7)
    error('Number of output steps exceeds 10M.  Probably mistake?');
end

%% setup video window

if (plotFast < 2)
    
    V1 = VideoWriter([ 'sim' num2str(movieNum) '_video' ],'Archival');

%     V1 = VideoWriter([ 'sim' num2str(movieNum) '_video' ],'MPEG-4');
%     V1.Quality = 100;
    open(V1);
end

FS = stoploop({'Press this button' 'to terminate simulation.'}) ;
xbox = 1.0*radius*x0;% + Ictrx;
zbox = (0.5*radius*x0*max(1/(pD/radius),0.8)) * (1+1/radius);

close all;

if (plotFast)
    fSimVid = figure('pos',[100 100 1400 1000]);
else
    fSimVid = figure('pos',[100 100 2000 800]);
end

set(gcf,'color','w');

%% setup snapshot window
fSnapshot = figure('pos',[200 100 600 1000]);  set(gcf,'color','w');
plotSnapshot;
hSnapshotLight2 = lightangle(00,20);
hSnapshotLight = lightangle(180,40);
axis equal;
colormap jet;
colorbar ;

%% set up cross section for hex
mycs1 = ringttl(1:radius-1) + round(1.5*(1:radius-1)+1);
mycs2 = ringttl(1:radius-1) + round(4.5*(1:radius-1)+3)+1;
mycs = [ fliplr(mycs2) 1 mycs1 ];

%% initialize monitors:
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
avgts = zeros(1,sim_dur_nto);          % Average temperature of remaining structure (mass weighted)
vertxs = zeros(1,sim_dur_nto);         % X position of vertex node
vertys = zeros(1,sim_dur_nto);         % Y position of vertex node
maxstrains = zeros(1,sim_dur_nto);     % Strain of most-tensioned edge
minstrains = zeros(1,sim_dur_nto);     % Strain of most-compressed edge
potenergies = zeros(1,sim_dur_nto);    % Total potential energy of the mesh (sum of 1/2*k*dx^2 over all remaining edges)
kinenergies = zeros(1,sim_dur_nto);    % Total kinetic energy of all vertices (including broken ones?)
kinzenergies = zeros(1,sim_dur_nto);   % z-direction energy of COM (1/2 * mass * com_velocity^2
kinoenergies = zeros(1,sim_dur_nto);   % Difference between above two KEs = rotational and vibrational KE.
ttlinpwr = zeros(1,sim_dur_nto);       % Total input power, sum of all power incident on all non-broken triangles (including upside down ones)
ttlabspwr = zeros(1,sim_dur_nto);      % Total absorbed power, sum across non-broken triangles
ttlradpwr = zeros(1,sim_dur_nto);      % Total radiated power, sum across non-broken triangles
graphuos = zeros(1,sim_dur_nto);       % Number or broken edges
graphts = zeros(1,sim_dur_nto);        % simulation time at each recording
graphdts = zeros(1,sim_dur_nto);       % time steps per output (obsolete - simulation now uses fixed dt)
upsidedowntris = zeros(1,sim_dur_nt);  % number of upside-down triangles (surface normal has negative Z)
mintrinormz = zeros(1,sim_dur_nto);    % Z component of most downward-facing triangle surface norm
areas = zeros(1,sim_dur_nto);          % Sum of all non-broken triangle areas
areasxy = zeros(1,sim_dur_nto);        % Sum of all non-broken triangles projected area onto xy plane
areasrightsideupxy = zeros(1,sim_dur_nto);  % Same as above but excluding triangles having surface norms with negative z ("upside down")
monitorqty = zeros(1,sim_dur_nto);  % misc. monitor vector, for various debug / tracking
vertextemps = zeros(1,sim_dur_nto);    % Temperature of vertex node
edgepttemps = zeros(1,sim_dur_nto);    % Temperature of last node (which is always somewhere on the edge)
lastptxyz = zeros(3, sim_dur_nto);     % XYZ coordiantes of the last node
lastptvxyz = zeros(3, sim_dur_nto);    % velocity of the last node (vx vy vz)
hoop12xyz = zeros(3, sim_dur_nto);     % XZY coordinates of the hoop node furthest from the hoop COM
hoop6xyz = zeros(3, sim_dur_nto);      % XYZ coord of the node opposite the hoop12 node (just by adding 1/2 to the ring index)
hoop9xyz = zeros(3, sim_dur_nto);      % XYZ coord of node opposite hoop3 node
hoop3xyz = zeros(3, sim_dur_nto);      % XYZ coord of the hoop node closest to hoop COM.  Note that hoop3 and hoop12 aren't always spaced 90 degrees apart on the ring!
hooprads = zeros(4, sim_dur_nto);      % This is the "radius" (distance to hoop COM) of all 4 
hoopcoms = zeros(3, sim_dur_nto);      % This is the hoop COM 
hoop12vxyz = zeros(3, sim_dur_nto);    % Velocity of hoop12 node
hoop3vxyz = zeros(3, sim_dur_nto);     % Vel of hoop3 node
% %hooprl = zeros(1, sim_dur_nto);
% %hooprs = zeros(1,sim_dur_nto);
% %hoopa = zeros(1,sim_dur_nto);
% %monitorcs = zeros(length(mycs),sim_dur_nto);    %this can track any property along a cross section (across the membrane along Y)
% %monitorfxs = zeros(length(n_x),sim_dur_nto);
% %monitorfys = zeros(length(n_y),sim_dur_nto);
if monitor_mode
    monitorfzs = zeros(length(n_z),monitor_depth);     %Mech. force on each node, in z direction
    monitorrfs = zeros(length(n_z),monitor_depth);     %mech force on each node, projected to 
    monitortfs = zeros(length(n_z),monitor_depth);
    monitorxs = zeros(length(n_x),monitor_depth);      %Complete position monitor, 
    monitorys = zeros(length(n_y),monitor_depth);      % All nodes
    monitorzs = zeros(length(n_z),monitor_depth);      % This is where all the memory has gone!
end

%% initialize buffers:
t_optFcn_r2 = zeros(size(t_na));
t_p1_x = zeros(size(t_na));
t_p1_y = zeros(size(t_na));
t_p1_z = zeros(size(t_na));
t_p2_x = zeros(size(t_na));
t_p2_y = zeros(size(t_na));
t_p2_z = zeros(size(t_na));
t_p3_x = zeros(size(t_na));
t_p3_y = zeros(size(t_na));
t_p3_z = zeros(size(t_na));
t_norms = zeros(3,length(t_na));
t_ix = zeros(size(t_na)); 
t_iy = zeros(size(t_na));
t_iz = zeros(size(t_na));
t_Ipwr = zeros(size(t_na));
t_abs = zeros(size(t_na));
t_vel = zeros(3,length(t_na));
%t_mf = zeros(3,length(t_na));

e_dl = zeros(3,length(e_na));   %bring in edge length and tension separately...
e_nl = zeros(3,length(e_na));
e_l = zeros(size(e_na));
e_t = zeros(size(e_na));
                               % " asdfasdfasdf ";
n_a = zeros(size(n_x));        % "Potential Energy" using matrix operator.  Fast, but some risk of failure.  

n_pec = zeros(3,length(n_x));  %"Potential energy" counter... not sure if this will work.
n_abs = n_a;
n_ems = n_a;
n_hf = n_a;

figure(fSimVid);
nt = 0;
nto = 0;
accz = 0; acct = 0;
anybroken = 0;
nto_firstBroken = sim_dur_nto;
time_firstBroken = inf;
pwrRampVal = 0;
% Potential & kinetic energy
PE = 0; KE = 0; KE_COM = 0; KEo = 0; 
end_reason = 'timeout';
plotMesh;

if (plotFast < 2)
    writeVideo(V1,getframe(gcf));
end

usercancel = 0;

% Vectors of simulation times for specific functions
timers = zeros(1,10);

startsimtic = tic;
lastdisptic = startsimtic;
lastdispnt = 0;

% Loop over 
for nt = 1:sim_dur_nt
       
    %% 1. Reset forcing counters on each node
    %n_mf(:) = 0;
    %n_of(:) = 0;
    %n_af(:) = 0;
    %n_a(:) = 0;
    %n_abs(:) = 0;
    %n_hf(:) = 0;
    
    
    %%  New 1.5.   Get centroid of each tri
%     if (nt == 1)
%         t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
%         t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
%         t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;
%         t_t  =  ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) / 3;
%     end
    %hold on;
    %plot3(t_cx, t_cy, t_cz, '+');
    %return
    
    %% 2.  Get light intensity at each node and store it there
    %  Note that I hard coded a gaussian here.  Trying to keep things fast!
    %     pwrRampVal = pwrRampFcn(tt);
    %     n_ix(:) = 0;
    %     n_iy(:) = 0;
    %     if (pwrRampVal == 0)
    %         n_iz(:) = 0;
    %     else
    %         n_optFcn_r2 = n_x.^2 + n_y.^2;
    %         n_iz = pwrRampVal * I0 * exp( -n_optFcn_r2 ./ (Irad^2) );
    %     end
    
    %% New 2. Get light intensity at each centroid and store it there!
    %  Note that I hard coded a gaussian here.  Trying to keep things fast!
    %     tt=1.5e-3;
    pwrRampVal = pwrRampFcn(tt);
    t_ix(:) = 0; 
    t_iy(:) = 0;
    if (pwrRampVal == 0)
        t_iz(:) = 0;
    else
        t_optFcn_r2 = (t_cx - Ictrx).^2 + t_cy.^2;
        t_iz = pwrRampVal * I0 * exp( - t_optFcn_r2 ./ (Irad^2) );  % This is in W/mm2.  Divide by c0 to get newtons / mm2 radiation pressure (N/mm2)
    end
    
    %warning:  hard-coded normalized vector used here:
    t_IDir = [ zeros(1,length(t_na)); zeros(1,length(t_na)); ones(1,length(t_na)); ];
    t_IMag = t_iz;
    %     hold off
    %     trisurf(TRI,n_x,n_y,n_z, t_iz);
    %     axis equal
    %     colorbar
    %     return
    
    %% 2.5 Calculate spin-up force
    if (tt <= (0 + dt))
        %n_ang = angle(n_x+1i*n_y);
        n_arad = sqrt( n_x.^2 + n_y.^2 );
        %n_af(1,:) = -9.8*spinup_gs/x0/radius.*n_m.*n_arad.*sin(n_ang);
        n_af(1,:) = -spinup_gs/88.2.*n_m.*(n_y) ;%%% .* (n_arad < 0.6 * radius * x0);
        %n_af(2,:) = 9.8*spinup_gs/x0/radius.*n_m.*n_arad.*cos(n_ang);
        n_af(2,:) = spinup_gs/88.2*n_m.*(n_x) ;%%%.* (n_arad < 0.6 * radius * x0);
        %%%n_af(3,:) = n_m .* (n_arad < 0.6 * radius * x0);
        n_rps = sqrt(n_vx.^2+n_vy.^2)./n_arad./pi/2;
        n_rps(1)=0;
        rot0 = mean(n_rps(2:end));
    end
    if ( tt >= 0) || anybroken
        n_af(:)=0;
    end
        % if myx < 0
                   %     mya = mya+pi;
                   % elseif myy < 0
                   %     mya = mya + 2*pi;
                   % end
    
    
    %% 3. Calculate triangle area and normal.... 
    
    % Assign x-coordinates of every triangle node/vertex
    t_p1_x = n_x(t_na); %this should be little penalty due to copy-on-write
    t_p2_x = n_x(t_nb);
    t_p3_x = n_x(t_nc);

    % Assign y-coordinates of every triangle node/vertex
    t_p1_y = n_y(t_na);
    t_p2_y = n_y(t_nb);
    t_p3_y = n_y(t_nc);

    % Assign z-coordinates of every triangle node/vertex
    t_p1_z = n_z(t_na);
    t_p2_z = n_z(t_nb);
    t_p3_z = n_z(t_nc);
    
    % Calculate normal to each triangle using cross product
    t_cps = cross( [ t_p2_x - t_p1_x ; t_p2_y - t_p1_y ; t_p2_z - t_p1_z ], ...
        [ t_p3_x - t_p1_x ; t_p3_y - t_p1_y ; t_p3_z - t_p1_z ] );
    %upsidedowntris(nt) = sum( t_cps(3,:) < 0);
    
    % Calculate length of normal vectors to each triangle, but why this way
    % instead of just norm()? And is the max command there to prevent zero
    % length and thus division by 0 later on?
    % Note that the norm of the cross-product a x b is also the area of the
    % parallegram being spanned by vectors a and b. Therefore, the area of
    % the triangle spanned by a and b is |a x b|/2. The missing factor 1/2
    % will be accounted for later
    t_a = max(sqrt(dot(t_cps,t_cps)), 1e-8);
%     t_a = max(vecnorm(t_cps), 1e-8);
    
    % Normalize normal vectors to each triangle
    t_norms = t_cps ./ t_a ;
    
    % Here is where the area of each triangle is correctly calculated by
    % including the factor 1/2 missing above
    t_a = 0.5 .* t_a;
    
    % Calculate size of each triangle regardless of whether it is broken or
    % not. If broken, use its original area calculated in rest position,
    % i.e. during mesh generation
    t_a = (t_notbroken).*t_a + t_broken.*t_a0;
                
    % Create matrix that contains one of two vectors spanning up each
    % triangle in the current mesh configuration
    t_edg1 = [  n_x(t_nb) - n_x(t_na)   ;
                n_y(t_nb) - n_y(t_na)   ;
                n_z(t_nb) - n_z(t_na)   ];
    
    % Normalize matrix containing one of two vectors spanning up each
    % triangle in current mesh configuration
    t_edg1n = t_edg1./(max(vecnorm(t_edg1),1e-8));
    
    % Noting that t_texa corresponds to angle between one spanning vector
    % of triangle and projected vector based on characteristic directional
    % vector for given texture, cos(t_texa) corresponds to the dot product
    % between the normalized projected vector within the texture and one of
    % the normalized edge vectors of the triangle
    t_texn = t_edg1n .* cos(t_texa) + ...
             cross(t_norms, t_edg1n ) .* sin(t_texa)  + ...
             t_norms .* dot(t_norms, t_edg1n ) .* (1 - cos(t_texa)) ;
         
    %t_texn = t_texn .* (t_tex>0);   % get rid of texture vectors on
    %unpatterned regions...?
    
    % Calculate cos(theta), where theta corresponds to the angle between
    % the normalized normal to each triangle and t_IDir, the directional
    % vector of the incident beam
    
    t_costheta = dot(t_norms, t_IDir);
    
    % Calculate reduced incident vector by first projecting normalized
    % normal vectors to each triangle, t_norms, to directional incident
    % vector, t_IDir, via cos(theta) * t_norms
    
    t_incProj = t_IDir - t_costheta .* t_norms;
%     t_incProj = t_IDir - dot(  t_norms, t_IDir ) .* t_norms;

    % Normalize projected incident vector
    t_incProjn = t_incProj ./ max(vecnorm( t_incProj ),1e-8);
    
    % Calculate cos(angle) between texture vector and projected incident
    % vector
    t_yaw_cos = dot( t_texn, t_incProjn );
    t_yaw = acos(t_yaw_cos)   .* ...
        -sign( dot( cross( t_incProjn, t_texn ), t_norms ) );
    
    % Calculate area of each triangle in xy plane from area in xyz via
    % calculated angle between incident vector and normal vector to
    % triangle. If this angle is zero, then cos(theta) = 1 and the area in
    % xy and xyz are the same. Consider only triangles with unbroken edges
    t_axy = t_a .* t_costheta .* (t_notbroken);
    
    %% 3.5 Here is where we can calculate reflection and absorption
    %  Available inputs:
    %  t_costheta: cosine of angle of incidence  
    %  t_t: local temperature
    %  t_x0, t_y0, t_z0: Initial coordinates of this triangle (prior to
    %  simulation beginning)
    %  t_thick: thickness of film here
    %  t_m: mass of this triangle
    %    
    % REQUIRED OUTPUTS: t_myA, t_myR, absorption and reflection
    % COEFFICIENTS, not amounts of power
    
    % Assign absorptivity and reflectivity to new variables
    t_myA = Iabs;
    t_myR = Irefl;
    
    % Calculate transmissivity from absorptivity and reflectivity
    t_myT = 1 - t_myA - t_myR;
    
    %% 4. Calculate optical forces & absorption heat-loading
    %  This gives the absorption heat input.  Radiatiave cooling occurs
    %  later at the node level.  That's why triangle area is assigned back
    %  to each node during the simulation.  
    
    % t_Ipwr = t_a .* dot(t_norms, [ t_ix; t_iy; t_iz ] );  % W, total incident laser power on this triangle, including costheta and area factors
    
    % Calculate total incident laser power on triange for given power
    % density t_IMag and projected triangle area into xy plane t_axy
    if includeMEPS
        t_Ipwr(indices_Region0) = abs(t_axy(indices_Region0) .* ...
            t_IMag(indices_Region0));  % Watt
    else
        t_Ipwr = abs(t_axy.*t_IMag);
    end
    % note some of these are negative...  I correct for that later.
    
    % Calculate absorbed power in each triangle
    t_abs = abs(t_myA .* (abs(t_axy.*t_IMag))); % Watt
    
    % note:  I am using the m/s value of c0 here, since conversion to mm
    % and 
    t_oth = 2 * t_myR * t_Ipwr ./ c0 .* t_norms .* t_costheta  + ... %  this is the reflection force
                   t_abs       ./ c0 .* t_IDir ;    % and this is the absorption force
    
            
    t_Ipwr = abs(t_Ipwr); %this corrects the incident power to be always positive.  It needs to be negative in the line above
    
    % Count amount of time needed for this calculation
    timers(1) = timers(1) + toc;
    tic;
    
    %% 4b. Calculate roll, pitch and yaw angles of each triangle
    
    if includeMEPS
    
    % Calculate individual elements from the rotation matrices that rotate
    % each mesh triangle into its current state, with formula specified in
    % Ramon's notes
    
    rotmat_31 = temp_tex_proj_norms(1,:).*t_texn(3,:) + ...
        (-temp_tex_proj_norms(3,:) .* temp_tri_normvec(2,:) + ...
        temp_tex_proj_norms(2,:) .* temp_tri_normvec(3,:)) .* ...
        (-t_texn(2,:).*t_norms(1,:) + t_texn(1,:).*t_norms(2,:)) + ...
        temp_tri_normvec(1,:).* t_norms(3,:);
    
    rotmat_32 = temp_tex_proj_norms(2,:).*t_texn(3,:) + ...
        (temp_tex_proj_norms(3,:) .* temp_tri_normvec(1,:) - ...
        temp_tex_proj_norms(1,:) .* temp_tri_normvec(3,:)) .* ...
        (-t_texn(2,:) .* t_norms(1,:) + t_texn(1,:) .* t_norms(2,:)) + ...
        temp_tri_normvec(2,:) .* t_norms(3,:);
    
    rotmat_33 = temp_tex_proj_norms(3,:).*t_texn(3,:) + ...
        (-temp_tex_proj_norms(2,:).*temp_tri_normvec(1,:) + ...
        temp_tex_proj_norms(1,:).*temp_tri_normvec(2,:)) .* ...
        (-t_texn(2,:) .* t_norms(1,:) + t_texn(1,:) .* t_norms(2,:)) + ...
        temp_tri_normvec(3,:) .* t_norms(3,:);
        
    rotmat_11 = temp_tex_proj_norms(1,:).*t_texn(1,:) + ...
        temp_tri_normvec(1,:) .* t_norms(1,:) + ...
        (-temp_tex_proj_norms(3,:).*temp_tri_normvec(2,:) + ...
        temp_tex_proj_norms(2,:).*temp_tri_normvec(3,:)) .* ...
        (-t_texn(3,:) .* t_norms(2,:) + t_texn(2,:) .* t_norms(3,:));
        
    rotmat_21 = temp_tex_proj_norms(1,:).*t_texn(2,:) + ...
        temp_tri_normvec(1,:) .* t_norms(2,:) + ...
        (-temp_tex_proj_norms(3,:).*temp_tri_normvec(2,:) + ...
        temp_tex_proj_norms(2,:).*temp_tri_normvec(3,:)) .* ...
        (t_texn(3,:) .* t_norms(1,:) - t_texn(1,:) .* t_norms(3,:));
        
    % Calculate pitch, roll and yaw angles based on formulas specified in
    % Ramon's notes. Note that technically, the yaw angles are not
    % important, as we assume the polarization to rotate synchronously with
    % the sail and thus its MEPS. However, calculating the yaw_angles
    % allows us to track the effect of spinning the sail

    pitch_angles = -asin(rotmat_31);
    roll_angles = atan2(rotmat_32./cos(pitch_angles), rotmat_33./cos(pitch_angles));
    yaw_angles = atan2(rotmat_21./cos(pitch_angles), rotmat_11./cos(pitch_angles));
    
    % Correct for negative yaw angles by adding 360� or 2*pi
    yaw_angles(yaw_angles < 0) = yaw_angles(yaw_angles < 0) + 2*pi; 
    
    % Allocate memory for current roll and pitch angles of each triangle
    specific_angles = zeros(max(size(t_tex)),2);

    % Assemble correct definition of angles for different regions. Here,
    % regions 3 & 4 are mirror-symmetric with respect to region 1 & 6, and
    % region 5 is mirror-symmetric with respect to region 2. Therefore,
    % roll and pitch angles for regions 3, 4 and 5 need to have opposite
    % signs compared to the roll and pitch angles for regions 1, 2 and 6,
    % since only regions 1 and 2 have been simulated (region 6 is equal to
    % region 1, region 3 is equal to region 4)
    
    specific_angles(indices_Region2,:) = ...
        [roll_angles(indices_Region2)', pitch_angles(indices_Region2)'];
    specific_angles(indices_Region3,:) = ...
        [roll_angles(indices_Region3)', pitch_angles(indices_Region3)'];
    specific_angles(indices_Region4,:) = ...
        [roll_angles(indices_Region4)', pitch_angles(indices_Region4)'];
    
    specific_angles(indices_Region1,:) = ...
        [-roll_angles(indices_Region1)', -pitch_angles(indices_Region1)'];
    specific_angles(indices_Region5,:) = ...
        [-roll_angles(indices_Region5)', -pitch_angles(indices_Region5)'];
    specific_angles(indices_Region6,:) = ...
        [-roll_angles(indices_Region6)', -pitch_angles(indices_Region6)'];
    
    
    % Get indices indicating in which range of roll angles every angle set
    % (roll, pitch) for each triangle can be found. This is needed later to
    % retrieve the pressures from the correct sub-look-up table
    
    specific_angles_indices1 = specific_angles(:,1) < ...
        (coarse_last_angle + finer_first_angle)/2;
    specific_angles_indices2 = specific_angles(:,1) >= ...
        (coarse_last_angle + finer_first_angle)/2 & specific_angles(:,1) < ...
        (finer_last_angle + even_finer_first_angle)/2;
    specific_angles_indices3 = specific_angles(:,1) >= ...
        (finer_last_angle + even_finer_first_angle)/2 & specific_angles(:,1) < ...
        (even_finer_last_angle + finest_first_angle)/2;
    specific_angles_indices4 = specific_angles(:,1) >= ...
        (even_finer_last_angle + finest_first_angle)/2 & specific_angles(:,1) < ...
        (finest_intermediate_angle1 + finest_angle_step/2);
    specific_angles_indices5 = specific_angles(:,1) >= ...
        (finest_intermediate_angle1 + finest_angle_step/2) & specific_angles(:,1) ...
        < (finest_intermediate_angle2 + finest_angle_step/2);
    specific_angles_indices6 = specific_angles(:,1) >= ...
        (finest_intermediate_angle2 + finest_angle_step/2) & specific_angles(:,1) ...
        < (-finest_intermediate_angle1 + finest_angle_step/2);
    specific_angles_indices7 = specific_angles(:,1) >= ...
        (-finest_intermediate_angle1 + finest_angle_step/2) & specific_angles(:,1) ...
        < (-finest_first_angle - even_finer_last_angle)/2;
    specific_angles_indices8 = specific_angles(:,1) >= ...
        (-finest_first_angle - even_finer_last_angle)/2 & specific_angles(:,1) < ...
        (-even_finer_first_angle - finer_last_angle)/2;
    specific_angles_indices9 = specific_angles(:,1) >= ...
        (-even_finer_first_angle - finer_last_angle)/2 & specific_angles(:,1) < ...
        (-finer_first_angle - coarse_last_angle)/2;
    specific_angles_indices10 = specific_angles(:,1) >= ...
        (-finer_first_angle - coarse_last_angle)/2;
    
    % Find pair of roll-pitch angles from look-up table of pressures that
    % are closest to the calculated ones using Matlab's built-in function.
    % Apparently, this is already fast, but still, this line IS the reason
    % for the significant slow-down of the simulation. I have alleviated
    % this a bit by removing some more extreme angles and the pressures at
    % these angles from the look-up table
    
    % Then, assign to each triangle its corresponding x, y and z pressure
    % given the definition of the respective grating orientation
    
    % Note that pressures calculated in COMSOL are normalized by the speed
    % of light, therefore, one has to multiply these by c0 to get the
    % actual pressures. However, to get the forces, the pressures need to
    % be weighted by 1/c again, "canceling" the effect of multiplying the
    % pressures with c here, and giving the correct order of magnitude of
    % the restoring forces compared to the forces from specular reflection
    % on the central region of the lightsail

    selected_pressures_x_TE = zeros(size(t_press_x));
    selected_pressures_z_TE = zeros(size(t_press_x));
    selected_pressures_x_TM = zeros(size(t_press_x));
    selected_pressures_z_TM = zeros(size(t_press_x));

    if ~isempty(specific_angles(specific_angles_indices1))
        [indices_angles1, ~] = knnsearch(pressure_TE_table1(:,1:2), ...
            specific_angles(specific_angles_indices1,:));
        selected_pressures_x_TE(specific_angles_indices1) = ...
            c0*pressure_TE_table1(indices_angles1,3);
        selected_pressures_z_TE(specific_angles_indices1) = ...
            c0*pressure_TE_table1(indices_angles1,4);
        selected_pressures_x_TM(specific_angles_indices1) = ...
            c0*pressure_TM_table1(indices_angles1,3);
        selected_pressures_z_TM(specific_angles_indices1) = ...
            c0*pressure_TM_table1(indices_angles1,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices2)) 
        [indices_angles2, ~] = knnsearch(pressure_TE_table2(:,1:2), ...
            specific_angles(specific_angles_indices2,:));
        selected_pressures_x_TE(specific_angles_indices2) = ...
            c0*pressure_TE_table2(indices_angles2,3);
        selected_pressures_z_TE(specific_angles_indices2) = ...
            c0*pressure_TE_table2(indices_angles2,4);
        selected_pressures_x_TM(specific_angles_indices2) = ...
            c0*pressure_TM_table2(indices_angles2,3);
        selected_pressures_z_TM(specific_angles_indices2) = ...
            c0*pressure_TM_table2(indices_angles2,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices3))
        [indices_angles3, ~] = knnsearch(pressure_TE_table3(:,1:2), ...
            specific_angles(specific_angles_indices3,:));
        selected_pressures_x_TE(specific_angles_indices3) = ...
            c0*pressure_TE_table3(indices_angles3,3);
        selected_pressures_z_TE(specific_angles_indices3) = ...
            c0*pressure_TE_table3(indices_angles3,4);
        selected_pressures_x_TM(specific_angles_indices3) = ...
            c0*pressure_TM_table3(indices_angles3,3);
        selected_pressures_z_TM(specific_angles_indices3) = ...
            c0*pressure_TM_table3(indices_angles3,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices4))
        [indices_angles4, ~] = knnsearch(pressure_TE_table4(:,1:2), ...
            specific_angles(specific_angles_indices4,:));
        selected_pressures_x_TE(specific_angles_indices4) = ...
            c0*pressure_TE_table4(indices_angles4,3);
        selected_pressures_z_TE(specific_angles_indices4) = ...
            c0*pressure_TE_table4(indices_angles4,4);
        selected_pressures_x_TM(specific_angles_indices4) = ...
            c0*pressure_TM_table4(indices_angles4,3);
        selected_pressures_z_TM(specific_angles_indices4) = ...
            c0*pressure_TM_table4(indices_angles4,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices5))
        [indices_angles5, ~] = knnsearch(pressure_TE_table5(:,1:2), ...
            specific_angles(specific_angles_indices5,:));
        selected_pressures_x_TE(specific_angles_indices5) = ...
            c0*pressure_TE_table5(indices_angles5,3);
        selected_pressures_z_TE(specific_angles_indices5) = ...
            c0*pressure_TE_table5(indices_angles5,4);
        selected_pressures_x_TM(specific_angles_indices5) = ...
            c0*pressure_TM_table5(indices_angles5,3);
        selected_pressures_z_TM(specific_angles_indices5) = ...
            c0*pressure_TM_table5(indices_angles5,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices6))
        [indices_angles6, ~] = knnsearch(pressure_TE_table6(:,1:2), ...
            specific_angles(specific_angles_indices6,:));
        selected_pressures_x_TE(specific_angles_indices6) = ...
            c0*pressure_TE_table6(indices_angles6,3);
        selected_pressures_z_TE(specific_angles_indices6) = ...
            c0*pressure_TE_table6(indices_angles6,4);
        selected_pressures_x_TM(specific_angles_indices6) = ...
            c0*pressure_TM_table6(indices_angles6,3);
        selected_pressures_z_TM(specific_angles_indices6) = ...
            c0*pressure_TM_table6(indices_angles6,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices7))
        [indices_angles7, ~] = knnsearch(pressure_TE_table7(:,1:2), ...
            specific_angles(specific_angles_indices7,:));
        selected_pressures_x_TE(specific_angles_indices7) = ...
            c0*pressure_TE_table7(indices_angles7,3);
        selected_pressures_z_TE(specific_angles_indices7) = ...
            c0*pressure_TE_table7(indices_angles7,4);
        selected_pressures_x_TM(specific_angles_indices7) = ...
            c0*pressure_TM_table7(indices_angles7,3);
        selected_pressures_z_TM(specific_angles_indices7) = ...
            c0*pressure_TM_table7(indices_angles7,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices8))
        [indices_angles8, ~] = knnsearch(pressure_TE_table8(:,1:2), ...
            specific_angles(specific_angles_indices8,:));
        selected_pressures_x_TE(specific_angles_indices8) = ...
            c0*pressure_TE_table8(indices_angles8,3);
        selected_pressures_z_TE(specific_angles_indices8) = ...
            c0*pressure_TE_table8(indices_angles8,4);
        selected_pressures_x_TM(specific_angles_indices8) = ...
            c0*pressure_TM_table8(indices_angles8,3);
        selected_pressures_z_TM(specific_angles_indices8) = ...
            c0*pressure_TM_table8(indices_angles8,4);
    end
    
    if ~isempty(specific_angles(specific_angles_indices9))
        [indices_angles9, ~] = knnsearch(pressure_TE_table9(:,1:2), ...
            specific_angles(specific_angles_indices9,:));
        selected_pressures_x_TE(specific_angles_indices9) = ...
            c0*pressure_TE_table9(indices_angles9,3);
        selected_pressures_z_TE(specific_angles_indices9) = ...
            c0*pressure_TE_table9(indices_angles9,4);
        selected_pressures_x_TM(specific_angles_indices9) = ...
            c0*pressure_TM_table9(indices_angles9,3);
        selected_pressures_z_TM(specific_angles_indices9) = ...
            c0*pressure_TM_table9(indices_angles9,4);
        end
    
    if ~isempty(specific_angles(specific_angles_indices10))
        [indices_angles10, ~] = knnsearch(pressure_TE_table10(:,1:2), ...
            specific_angles(specific_angles_indices10,:));
        selected_pressures_x_TE(specific_angles_indices10) = ...
            c0*pressure_TE_table10(indices_angles10,3);
        selected_pressures_z_TE(specific_angles_indices10) = ...
            c0*pressure_TE_table10(indices_angles10,4);
        selected_pressures_x_TM(specific_angles_indices10) = ...
            c0*pressure_TM_table10(indices_angles10,3);
        selected_pressures_z_TM(specific_angles_indices10) = ...
            c0*pressure_TM_table10(indices_angles10,4);
    end
    
    % Transform local frames of individual regions to global body frame.
    % Note that the lightsail is flipped upside down for this simulation
    % compared to those in COMSOL, thus requiring appropriate sign changes
    % for respective pressures
    
    t_press_x(indices_Region1) = selected_pressures_x_TE(indices_Region1);
    t_press_y(indices_Region1) = 0;
    t_press_z(indices_Region1) = -selected_pressures_z_TE(indices_Region1);

    t_press_x(indices_Region2) = 0;
    t_press_y(indices_Region2) = selected_pressures_x_TM(indices_Region2);
    t_press_z(indices_Region2) = -selected_pressures_z_TM(indices_Region2);

    t_press_x(indices_Region3) = -selected_pressures_x_TE(indices_Region3);
    t_press_y(indices_Region3) = 0;
    t_press_z(indices_Region3) = -selected_pressures_z_TE(indices_Region3);

    t_press_x(indices_Region4) = -selected_pressures_x_TE(indices_Region4);
    t_press_y(indices_Region4) = 0;
    t_press_z(indices_Region4) = -selected_pressures_z_TE(indices_Region4);

    t_press_x(indices_Region5) = 0;
    t_press_y(indices_Region5) = -selected_pressures_x_TM(indices_Region5);
    t_press_z(indices_Region5) = -selected_pressures_z_TM(indices_Region5);
 
    t_press_x(indices_Region6) = selected_pressures_x_TE(indices_Region6);
    t_press_y(indices_Region6) = 0;
    t_press_z(indices_Region6) = -selected_pressures_z_TE(indices_Region6);

    % Calculate body-frame forces by discretizing the double integral into a
    % multiplication by the area of a mesh element (see Ramon's ppt
    % documentation for details), while accounting for the projected area
    % Note that the unit of the intensity [t_IMag] is W/(mm^2), while the unit
    % of the triangle areas [t_a] is mm^2, and the speed of light [c0] being
    % m/s, we get the correct unit for the forces, namely W*s/m = J/m = N
    % Note that cos(pitch) * cos(roll) accounts for the projected area

    t_force_x_BF = cos(pitch_angles) .* cos(roll_angles) .* t_a .* ...
        (t_IMag./c0) .* t_press_x;
    t_force_y_BF = cos(pitch_angles) .* cos(roll_angles) .* t_a .* ...
        (t_IMag./c0) .* t_press_y;
    t_force_z_BF = cos(pitch_angles) .* cos(roll_angles) .* t_a .* ...
        (t_IMag./c0) .* t_press_z;

    % Transform forces on triangle in local body frame to global laser frame by
    % multiplying body-frame force vector with direction cosine matrix H_I^B
    % according to OIlic2019
    
    t_force_x = cos(pitch_angles) .* t_force_x_BF + ...
        sin(pitch_angles) .* t_force_z_BF;
    t_force_y = -sin(pitch_angles) .* sin(roll_angles) .* t_force_x_BF + ...
        cos(roll_angles) .* t_force_y_BF + ...
        cos(pitch_angles) .* sin(roll_angles) .* t_force_z_BF;
    t_force_z = -sin(pitch_angles) .* cos(roll_angles) .* t_force_x_BF - ...
        sin(roll_angles) .* t_force_y_BF + ...
        cos(pitch_angles) .* cos(roll_angles) .* t_force_z_BF;
   
    end
        
    %% 5. Distribute optical forces (& absorbed heat) to the nodes
    
    if includeMEPS
        % Sum up forces from MEPS/gratings (first terms) and forces from
        % specular region at the center of the lightsail (second terms)
        
        t_oth(1,:) =  t_oth(1,:) + t_force_x;
        t_oth(2,:) =  t_oth(2,:) + t_force_y;
        t_oth(3,:) =  t_oth(3,:) + t_force_z;
        
        n_of(1,:) = t_oth(1,:) * M_t2n ./ 3;
        n_of(2,:) = t_oth(2,:) * M_t2n ./ 3;
        n_of(3,:) = t_oth(3,:) * M_t2n ./ 3; 
    else
        n_of(1,:) = t_oth(1,:) * M_t2n ./ 3;
        n_of(2,:) = t_oth(2,:) * M_t2n ./ 3;
        n_of(3,:) = t_oth(3,:) * M_t2n ./ 3; 
    end
    n_a = t_a * M_t2n ./ 3;
    n_abs = t_abs * M_t2n ./ 3;
    
%     for ntr=1:length(t_na)
%         if (t_broken(ntr))
%            % n_a( [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) = n_a( [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) + [t_a0(ntr) t_a0(ntr) t_a0(ntr)] ./ 3;
%             n_a(   t_na(ntr) )  = n_a( t_na(ntr) ) + t_a0(ntr);
%             n_a(   t_nb(ntr) )  = n_a( t_nb(ntr) ) + t_a0(ntr);
%             n_a(   t_nc(ntr) )  = n_a( t_nc(ntr) ) + t_a0(ntr);
%         else
%         %n_a(  [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] )  = n_a( [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) + [ t_a(ntr)  t_a(ntr) t_a(ntr)]./3;
%         n_a(   t_na(ntr) )  = n_a( t_na(ntr) ) + t_a(ntr);
%         n_a(   t_nb(ntr) )  = n_a( t_nb(ntr) ) + t_a(ntr);
%         n_a(   t_nc(ntr) )  = n_a( t_nc(ntr) ) + t_a(ntr);
%         
%         n_of(:,[ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) = n_of(:, [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) + t_oth(:,[ ntr ntr ntr ]);
%         %n_of(:,t_na(ntr)) = n_of(:,t_na(ntr)) + t_oth(:,ntr);
%         %n_of(:,t_nb(ntr)) = n_of(:,t_nb(ntr)) + t_oth(:,ntr);
%         %n_of(:,t_nc(ntr)) = n_of(:,t_nc(ntr)) + t_oth(:,ntr);
%         
%         %n_abs( [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) = n_abs( [ t_na(ntr) t_nb(ntr) t_nc(ntr) ] ) + t_abs(ntr)./3;
%         n_abs( t_na(ntr) ) = n_abs(  t_na(ntr)) + t_abs(ntr);
%         n_abs( t_nb(ntr) ) = n_abs(  t_nb(ntr)) + t_abs(ntr);
%         n_abs( t_nc(ntr) ) = n_abs(  t_nc(ntr)) +
%         t_abs(ntr);trtrrtrrrrrrrrtrrttb 
%         
%         end
%         %t_vel(:,ntr) = n_vx(t_na(ntr)) + n_vy(t_nb(ntr)) + n_vz(t_nc(ntr));
%     end
%     n_a   = n_a   ./ 3;
%     n_of  = n_of  ./ 3;
%     n_abs = n_abs ./ 3;
    
    if (isnan(sum(n_abs)))
        error('Found NaN!!!');
    end
    
    % Count amount of time needed for this calculation
    timers(2) = timers(2) + toc;
    tic;
    

    %% 5.5 Calculate thermal emission (radiatively emitted power)
    % Calculate radiated power using Stefan-Boltzmann's law given by P_rad
    % = A * emissivity * sigma * ( T^4 - T_ref^4). Here, the radiated power
    % at each node is calculated.
    % Mike, for sake of correctness, should the temperature of space as the
    % reference temperature be included here?
    n_ems = SBCmm .* n_a .* Emissivity .* (n_t.^4); % Watt
    
    % I suggest this to be correct?
%     t_space = 2.7; % K, black body temperature of background radiation in space
%     n_ems = SBCmm .* n_a * Emissivity .* (n_t.^4 - t_space^4);  
    
    
    %% 6. Calculate mechanical forces & heat conduction
    %  mincomp = 1e6;
    %  maxcomp = 0;
    
    % Matrix containing all edge direction vectors with x, y, z components
    % as columns
    e_nl = [ n_x(e_nb) - n_x(e_na) ; n_y(e_nb) - n_y(e_na) ; ...
        n_z(e_nb) - n_z(e_na) ];
    
    % Calculate current lengths in mm of all edges 
    e_l = sqrt( dot( e_nl, e_nl ));
    
    % Normalize all edge direction vectors. Here 'e_nl' stands for edge,
    % normalized length
    e_nl = e_nl ./ e_l;   %normlizing the edge direction vector
    
    % Here, calculate linear expansion due to temperature difference e_t
    % minus initial temperature t0 for each edge multiplied by the
    % coefficient of thermal expansion
    % Amount of thermal expansion can be described by material strain
    % epsilon_thermal = (L_final - L_initial)/L_initial, where the strain
    % is due change in temperature proportional to the coefficient of 
    % thermal expansion, epsilon_thermal = -alpha_L * (T_final - T_initial)
    % This can be rewritten in terms of length change dL / L_initial = 
    % 1 - (1 + alpha_L * (T_final - T_initial) )
    
    e_tex = 1 + (e_t - t0).*(CTE*1e-6);
    e_dl = e_l - (e_l0 .* e_tex) ;    % mm, edge length difference
    
    % Calculate strain on every edge due to thermal expansion
    e_s = e_l ./ e_l0 - 1;   % strain %/100.   No stress -> 0 strain.  Elongate to 2x length -> 1.0 strain
    
    % Determine maximum and minimum strain on any of the edges in the mesh
    max_strain = max(e_s.*(e_notbroken));
    min_strain = min(e_s.*(e_notbroken));
    
    % Calculate mechanical force on edge as the the edge length difference
    % times the edge stiffness/spring constant being equal to the 2D
    % Young's modulus x area / length^2 (see generateMesh.m for more info).
    % The sign of the mechanical force on a specific edge is determined by
    % the corresponding edge direction vector 'e_nl'
    e_mf = (e_notbroken) .* (Ymodmm .* e_al) .* e_dl .* e_nl;  % Newtons
    
    e_mf_mag2 = Ymodmm * e_al .* e_dl.^2 .* (e_notbroken);
    
    %PE = 5e-4 * Ymodmm * sum( e_al .* e_dl.^2 .* (e_notbroken) );  % Joules maybe? 
    PE = 5e-4 * sum(e_mf_mag2);
    
    % Calculate temperature difference across all edges due to different
    % temperatures on every node
    e_dt = n_t(e_nb) - n_t(e_na);   %convention: add to a
    
    % Thermal conduction, i.e. heat flow power on edges ('e_hf') according
    % to formula kappa * (A/l) * (T_hot - T_cold)
    e_hf = Thermcondmm .* e_al .* e_dt;  % Watt
    
    % What is this for?
    if isnan(sum(e_hf))
        error('Found it!')
    end

    % Count amount of time needed for this calculation
    timers(3) = timers(3) + toc;
    tic;
    
    %% 7. Assign mechanical forces to nodes (iterate) and also heat flow
    n_mf(1,:) = e_mf(1,:) * M_e2n;
    n_mf(2,:) = e_mf(2,:) * M_e2n;
    n_mf(3,:) = e_mf(3,:) * M_e2n;   
    n_pec(1,:) = e_mf_mag2 .* abs(e_nl(1,:)) * abs(M_e2n);
    n_pec(2,:) = e_mf_mag2 .* abs(e_nl(2,:)) * abs(M_e2n);
    n_pec(3,:) = e_mf_mag2 .* abs(e_nl(3,:)) * abs(M_e2n);
    
    % Convert thermal conduction (power) on edges to thermal conduction
    % (power) on nodes
    n_hf = e_hf * M_e2n;
    
%     for nedg=1:length(e_na)
%         if e_broken(nedg)
%         
%         else
%         %n_mf(:,[e_na(nedg) e_nb(nedg)]) = n_mf(:,[e_na(nedg) e_nb(nedg)]) + e_mf(:,[nedg nedg]).*[1 -1] ;
%         n_mf(:,e_na(nedg)) = n_mf(:,e_na(nedg)) + e_mf(:,nedg) ;
%         n_mf(:,e_nb(nedg)) = n_mf(:,e_nb(nedg)) - e_mf(:,nedg) ;
%         %n_hf([e_na(nedg) e_nb(nedg)]) = n_hf([e_na(nedg) e_nb(nedg)]) - e_hf([nedg nedg]).*[-1 1] ;
%         n_hf(e_na(nedg)) = n_hf(e_na(nedg)) - e_hf(nedg);
%         n_hf(e_nb(nedg)) = n_hf(e_nb(nedg)) + e_hf(nedg);
%         end
%     end

    % Count amount of time needed for this calculation
    timers(4) = timers(4) + toc;
    tic;
   
    
    %% 8. Sum forces
    % Total forces = optical + mechanical + acceleration/spin-up forces
    
    n_fx = n_mf(1,:) + n_of(1,:) + n_af(1,:);
    n_fy = n_mf(2,:) + n_of(2,:) + n_af(2,:);
    n_fz = n_mf(3,:) + n_of(3,:) + n_af(3,:);
    
    % wave excitation kick:
%    if (nt<kick_interval)
%        n_fy(end-2*radius:end-radius) = n_fy(end-2*radius:end-radius) -10/kick_interval*sin(pi*nt/kick_interval); %P-wave
%      %  n_fx(n_y < (origmeshymin + x0)) = n_fx(n_y < (origmeshymin + x0)) -2/kick_interval*sin(pi*nt/kick_interval); % S-wave
%    end
    
    
    %% 9. Evaluate for brokenness
    
    % Use element-by-element logical OR operation to get broken bonds. If
    % an edge is broken, it will be indexed as '1' in e_broken, otherwise
    % '0'. 'e_s' is a vector describing the strain for each edge. If the
    % strain is greater than the tensile strain limit, then this edge is
    % evaluated to be '1' due to the logical OR, thus becoming broken. Any
    % edge that is already broken will remain broken forever due to the
    % logical OR operation.
    e_broken =  ( e_broken | ( e_s > tensile_strain_limit ));
    
    % Count number of broken bonds/edges and call it 'uhohs'
    uhohs = sum(e_broken);
    
    % Vector indexing for each triangle whether it is broken or not by
    % containing at least one broken bond/edge
    t_broken =  ( e_broken(t_e1) | e_broken(t_e2) | e_broken(t_e3) );
    
    % Vector indexing unbroken edges/bonds as '1' and broken ones as '0'
    e_notbroken = ~e_broken;
    
    % Vector indexing unbroken triangles as '1' and broken ones as '0'
    t_notbroken = ~t_broken;
    
    % 'any' command returns '1' if any of the vector elements is nonzero.
    % 'anybroken' is a number that is '1' if there is at least one broken
    % bond in the mesh
    anybroken = anybroken || any(t_broken);
    
    % Determine step and time at which first breaking occurs
    if (anybroken) && (nto_firstBroken == sim_dur_nto)
        nto_firstBroken = nto;
        time_firstBroken = tt;
    end
        
    
    %% 10. Allow time step, update positions
    
    % Implementation of velocity Verlet algorithm as integrator for the
    % propulsion phase, while the explicit Euler method's is employed for
    % the spin-up phase prior to release & propulsion. This distinction is
    % necessary because the velocities calculated with the Verlet algorithm
    % will lag behind the positions (and accelerations) by one time step.
    % While this is OK during propulsion, because only the positions are
    % needed to calculate the forces, with velocities only needed to
    % calculate the kinetic energy and for logging purposes, during the
    % spin-up phase (see section 2), the velocities will actually be needed
    % in 'real-time'. Note that it is assumed that the accelerations only
    % depend on positions via the forces, i.e. a = F/m
    
    % Compared to the explicit Euler's method, the velocity Verlet
    % integration is of second order, thus being an order better than the
    % former in terms of dt. The lagging velocities 
    
    % It seems that Mike assumes the node masses to remain fixed during the
    % whole simulation after the mesh generation, but doesn't deformation
    % of the sail result in changes in triangle areas and thus changes in
    % node masses?
    
    if (tt <= (0 + dt)) % explicit Euler's method
        
        % Remember that [n_m] = g, so this accounts for 1e3
        
        n_vx = n_vx + 1e6 .* (n_fx./n_m) .* dt;
        n_vy = n_vy + 1e6 .* (n_fy./n_m) .* dt;
        n_vz = n_vz + 1e6 .* (n_fz./n_m) .* dt;

        n_x = n_x + n_vx .* dt;
        n_y = n_y + n_vy .* dt;
        n_z = n_z + n_vz .* dt;
    
    else % velocity Verlet integration
        if first_switch_integrator
            
            % Do this only once, namely at the beginning, i.e., after the
            % spin-up phase upon release/propulsion. This is to 'make' the
            % velocities lag behind the positions, or in other words, to
            % avoid calculating the velocities during this step
            
            n_accx = 1e6 .* (n_fx ./ n_m);
            n_accy = 1e6 .* (n_fy ./ n_m);
            n_accz = 1e6 .* (n_fz ./ n_m);
        
            n_x = n_x + n_vx .* dt + 0.5 .* n_accx .* (dt).^2;
            n_y = n_y + n_vy .* dt + 0.5 .* n_accy .* (dt).^2;
            n_z = n_z + n_vz .* dt + 0.5 .* n_accz .* (dt).^2;
        
            first_switch_integrator = 0; % Such that never land here again
            
            dt = 2*dt; % Accounting for second-order integration method
            
        else
            
            n_vx = n_vx + 0.5 .* (n_accx + 1e6 .* (n_fx ./ n_m) ) .* dt;
            n_vy = n_vy + 0.5 .* (n_accy + 1e6 .* (n_fy ./ n_m) ) .* dt;
            n_vz = n_vz + 0.5 .* (n_accz + 1e6 .* (n_fz ./ n_m) ) .* dt;
       
            n_accx = 1e6 .* (n_fx ./ n_m);
            n_accy = 1e6 .* (n_fy ./ n_m);
            n_accz = 1e6 .* (n_fz ./ n_m);
       
            n_x = n_x + n_vx .* dt + 0.5 .* n_accx .* (dt).^2;
            n_y = n_y + n_vy .* dt + 0.5 .* n_accy .* (dt).^2;
            n_z = n_z + n_vz .* dt + 0.5 .* n_accz .* (dt).^2;
            
        end
    end
    
    % Start calculating center-of-mass of sail as mass-weighted sum of
    % every x/y/z coordinate divided by total mass of sail, by calculating
    % the numerator here first
    
    coms_x = n_x .* n_m;
    coms_y = n_y .* n_m;
    coms_z = n_z .* n_m;
    
    % Calculate temperature
    % Heat capacity links the temperature rise of a body to the required
    % amount of energy. The unit of heat capacity is [J/g/degC], which
    % describes the amount of energy required to raise the temperature of 1
    % gram of an object by 1�C (or 1 Kelvin). Energy is equal to power
    % times time, such that the temperature rise is given by energy / mass
    % / heat capacity, where energy = (input power - output power) * dt.
    % Input power is given by the opticaly absorbed power and power from
    % heat flow/thermal conduction, whereas output power is due to
    % radiative cooling.
    
    n_t = n_t + (n_hf + n_abs - n_ems).*dt ./ heatCap ./ n_m;
    
    % Calculate the edge temperature as the averaged temperature of the sum
    % of the temperatures of the two nodes that form the edge
    e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;
    
    % Advance time by step dt
    tt = tt + dt;
    
    % Calculate kinetic energy as KE = (1/2) * m * v^2, which lags behind
    % during propulsion by one time step, i.e. dt, due to velocity Verlet
    
    KE = 0.5e-9 * sum( n_m .* (n_vx.^2 + n_vy.^2 + n_vz.^2) ); % Joules
    
    
    %% 11. Update mean velocity, center of mass...
    
    meanvx = sum(n_vx.*n_m)./totalmass_nodes;
    meanvy = sum(n_vy.*n_m)./totalmass_nodes;
    meanvz = sum(n_vz.*n_m)./totalmass_nodes;
    maxvrel = max( sqrt( (n_vx - meanvx).^2 + (n_vy - meanvy).^2 + ...
        (n_vz - meanvz).^2 ) );

    KE_COM = 0.5e-9 * totalmass_nodes * norm([meanvz meanvy meanvx]).^2;
    KEo = KE - KE_COM;
    
    % Complete calculating center-of-mass of sail as mass-weighted sum of
    % every x/y/z coordinate divided by total mass of sail, by taking the
    % already calculated numerators and dividing by total mass of nodes
    
    comx = sum(coms_x)./totalmass_nodes ;
    comy = sum(coms_y)./totalmass_nodes ;
    comz = sum(coms_z)./totalmass_nodes ;
    
    % Calculate triangle centroids, where the x, y or z coordinate of a
    % centroid is calculated as the arithmetic mean of the respective x, y
    % or z vertex/nodes coordinates

    t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
    t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
    t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;
    
    % Calculate temperature of triangle as the average of the temperatures
    % of the nodes that make up the triangle
    t_t  =  ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) / 3;
    
    %% Decompose mf to find radial and tangential forces
    if monitor_mode
        n_xrelcom = n_x - comx;
        n_yrelcom = n_y - comy;
        n_rrelcom = sqrt(n_xrelcom.^2+n_yrelcom.^2);
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
    
    
    %% Save all state variables and outputs for future plotting vs. time
    if ( ~mod(nt, sim_output_downsample) )
        nto = nto + 1;
    
        comxs(nto) = comx;
        comys(nto) = comy;
        comzs(nto) = comz;
        comrs(nto) = sqrt( (comx-Ictrx).^2 + comy.^2);
        velxs(nto) = meanvx;
        velys(nto) = meanvy;
        velzs(nto) = meanvz;
        maxvs(nto) = maxvrel;
        vertxs(nto) = n_x(1);
        vertys(nto) = n_y(1);
        graphuos(nto) = uhohs;
        graphts(nto) = tt;
        avgts(nto) = sum(t_t.*t_m)./sum(t_m);
        maxts(nto) = max(t_t);
        maxstrains(nto) = max_strain;
        minstrains(nto) = min_strain;
        potenergies(nto) = PE;
        kinenergies(nto) = KE;
        kinzenergies(nto) = KE_COM;
        kinoenergies(nto) = KEo;
        ttlinpwr(nto)  = sum(t_Ipwr.*(t_notbroken));
        ttlabspwr(nto) = sum(t_abs.*(t_notbroken));
        ttlradpwr(nto) = sum(n_ems);
        areas(nto) = sum(t_a.*(t_notbroken));
        areasxy(nto) = sum(t_axy.*(t_notbroken));
        upsidedowntris(nto) = sum( t_cps(3,:) < 0);
        areasrightsideupxy(nto) = sum( t_axy.*(t_axy>0).*(t_notbroken));
        mintrinormz(nto) = min( t_norms(3,:) );
       % monitorqty(nto) = mean(n_rps(2:end));
        edgepttemps(nto) = n_t(end);
        vertextemps(nto) = n_t(1);
        %monitorfxs(:,nto) = n_mf(1,:);
        %monitorfys(:,nto) = n_mf(2,:);
        if monitor_mode
            monitorfzs(:,mod(nto-1,monitor_depth)+1) = n_mf(3,:);
            monitorrfs(:,mod(nto-1,monitor_depth)+1) = n_rmf;
            monitortfs(:,mod(nto-1,monitor_depth)+1) = n_tmf;
            monitorxs(:,mod(nto-1,monitor_depth)+1) = n_x;
            monitorys(:,mod(nto-1,monitor_depth)+1) = n_y;
            monitorzs(:,mod(nto-1,monitor_depth)+1) = n_z;
        end
        if(nto==1)
            accz = meanvz ./ dt;
            acct = norm( [meanvx meanvy meanvz ] ) ./dt;
        else
            accz = (meanvz - velzs(nto-1))./(dt*sim_output_downsample);
            acct = norm( [meanvx meanvy meanvz] - [velxs(nto-1) velys(nto-1) velzs(nto-1) ] )./(dt*sim_output_downsample);
        end
        acczs(nto)=accz;
        accts(nto)=acct;
%        lastptxyz(:,nto) = [n_x(end); n_y(end); n_z(end)];
%        lastpt

        %monitorcs(:,nto) = n_fx(mycs);
        monitorqty(nto) = rot0;
        
        %  elongation ratio monitors
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
        %idx3 = mod(furthestidx + round(numringpts*0.25) -1 , numringpts) + 1;
        idx9 = mod(closestidx + round(numringpts/2) -1 , numringpts) + 1;
        
        hoopcoms(:,nto) = ringcom;
        hooprads(:,nto) = [ringradii(furthestidx); ringradii(closestidx); ringradii(idx6); ringradii(idx9) ];
                
        hoop12xyz(:,nto) = [ n_x(outerringstartidx+furthestidx-1); n_y(outerringstartidx+furthestidx-1); n_z(outerringstartidx+furthestidx-1) ];
        hoop6xyz(:,nto) = [ n_x(outerringstartidx+idx6-1); n_y(outerringstartidx+idx6-1); n_z(outerringstartidx+idx6-1) ];
        hoop3xyz(:,nto) = [ n_x(outerringstartidx+closestidx-1); n_y(outerringstartidx+closestidx-1); n_z(outerringstartidx+closestidx-1) ];
        hoop9xyz(:,nto) = [ n_x(outerringstartidx+idx9-1); n_y(outerringstartidx+idx9-1); n_z(outerringstartidx+idx9-1) ];
        
        hoop12vxyz(:,nto) = [ n_vx(outerringstartidx+furthestidx-1); n_vy(outerringstartidx+furthestidx-1); n_vz(outerringstartidx+furthestidx-1) ];
        hoop3vxyz(:,nto) = [ n_vx(outerringstartidx+closestidx-1); n_vy(outerringstartidx+closestidx-1); n_vz(outerringstartidx+closestidx-1) ];
        
        lastptxyz(:,nto) = [n_x(end); n_y(end); n_z(end)];
        lastptvxyz(:,nto) = [n_vx(end); n_vy(end); n_vz(end)];
    end
    
    timers(5) = timers(5) + toc;
    tic;
    
    %display plot for movie:
    if (nt - pf >= frame_interval - 1)
        pf = nt + 1;
        nf = nf + 1;
        
        
                disp(['n=' num2str(nt) '  nto = ' num2str(nto) '  nf=' num2str(nf) '  tt=' num2str(tt) '  max_str=' num2str(max_strain) '    maxv=' num2str(maxvrel) ...
            '      uh-ohs: ' num2str(uhohs) '   vz=' num2str(meanvz) '  vx=' num2str(meanvx) '  vy=' num2str(meanvy) ...
             '      z=' num2str(comz) '  x=' num2str(comx) '  y=' num2str(comy) ...
             '  PE=' num2str(PE) '  KE=' num2str(KE) ...
         ... %   '      tavg=' num2str(tempavg) ' tmax=' num2str(temppk) ' qin=' num2str(qin) ' qout=' num2str(qout) ...
        ...  %   '      mincomp=' num2str(mincomp) ' maxcomp=' num2str(maxcomp) ' ET=' num2str(toc) 's' ]);
        ' CPS=' num2str((nt-lastdispnt)/toc(lastdisptic)) ' ET=' num2str(toc(startsimtic)) 's' ]);
        lastdispnt = nt;
        
       % fprintf('%-4d  Node 5:  pos=%+.3e,%+.3e,%+.3e     vel=%+.3e,%+.3e,%+.3e      mf=%+.3e,%+.3e,%+.3e       of=%+.3e,%+.3e,%+.3e \n', nt, n_x(5), n_y(5), n_z(5), n_vx(5), n_vy(5), n_vz(5), n_mf(1,5), n_mf(2,5), n_mf(3,5), n_of(1,5), n_of(2,5), n_of(3,5) );
        
        if (gcf ~= fSimVid) 
            figure(fSimVid);
        end
        
        % Call self-written external function
        plotMesh;
        if (plotFast < 2)
            writeVideo(V1,getframe(gcf));
        end
        
        timers(6) = timers(6) + toc;   
        lastdisptic = tic;
    else
        timers(6) = timers(6) + toc;
        tic;
    end
    
    
    if (tt < 0)
        ps = nt;
    end
    
    if (plotFast < 2) && (uhohs > 0)
        frame_interval = max(floor(breakup_slowmo_interval_s/dt),1); %slow down for breakup!
    end
    
    if ((nt-ps >= snapshot_interval-1 ) && (ns < 10) ) || (uhohs >= nextExplosionSnapshot)
        ps=nt+1;
        ns=ns+1;
        if (uhohs >= nextExplosionSnapshot)
            nextExplosionSnapshot = nextExplosionSnapshot * explosionSnapshotRatio;
        end
        figure(fSnapshot);
        plotSnapshot;
    end
    
    timers(7) = timers(7) + toc;
    tic;
    
    xextent = max(abs(n_x-comx));
    xextentratio = xextent / (radius*x0);
    COMBeamDistNorm = sqrt( (comx - Ictrx)^2 + comy^2)./(radius*x0);
        
    if ~mod(nt,min(500,5*frame_interval))
        usercancel = FS.Stop();
    end

    if(usercancel || (max(n_vx) > 1e12) || (xextentratio > 5 ) || ...
            (COMBeamDistNorm > stopSimulation_Displacement))
        if (max(n_vx) > 1e12)
            end_reason = 'diverged';
        elseif (xextentratio > 5 )
            end_reason = 'exploded';
        elseif (COMBeamDistNorm > stopSimulation_Displacement)
            if ~any(e_broken)
                if (tt > (t_start_ramp+2*t_ramp_dur+t_on_dur) )
                    end_reason = 'launched';
                else
                    end_reason = 'flyaway';
                end
            else
                end_reason = 'brokeaway';
            end
        else
            end_reason = 'cancelled';
        end
        if (isfinite(snapshot_interval))
            ns=ns+1;
            ps=ps+1;
            figure(fSnapshot);
            plotSnapshot;
        end
        break;
    end
    timers(8) = timers(8) + toc;
    tic;
    
end

fprintf('SIMULATION ENDED:  ');
disp(end_reason);

ttltimers = sum(timers);
et=toc(startsimtic);
rate_oa = nt/et;
rate_loop = nt/(sum(timers(1:5))+timers(8));
rate_rend = (nf+ns)/(timers(6)+timers(7));
disp(['Timer results:    Total time  ' num2str(et)               ' sec.     Rate = ' num2str(rate_oa) ' loops/sec']);
disp(['                  Calc time   ' num2str(sum(timers([1:5 8]))) ' sec.      Rate = ' num2str(rate_loop) ' loops/sec']);
disp(['                  Render time ' num2str(sum(timers(6:7)))        ' sec.      Rate = ' num2str(rate_rend) ' frames/sec']);
disp(  ['Breakdown:  Optical      Tri_iter     Mechanical   Edge_iter    Evolve       PlotMovie    Snapshots   EndCalcs']);
fprintf('      '); 
fprintf('   %10.5f', timers(1:8)./et );  
fprintf('\n');

if (plotFast < 2)
close(V1);
FS.Clear();
end

nto_firstBroken = min(nto_firstBroken,nto+1);

longaxis_length = vecnorm(hoop12xyz-hoop6xyz);
shortaxis_length = vecnorm(hoop3xyz-hoop9xyz);

nt=nto;   % go back to original definition of nt!
plot_time = graphts(1:nt);

if monitor_mode
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
        warning('The monitors wrapped!!! fix me!!!');
    end
end

figure('pos',[100 100 2100 900]);
set(gcf,'color','w');


hafp1 = subplot(3,5,1);
plot(plot_time, velxs(1:nt));
title('X & Y velocity (mm/s)');
xlabel('Time');
hold on
plot(plot_time, velys(1:nt));
legend({'X' 'Y'});


hafp2 = subplot(3,5,2);
plot(plot_time, velzs(1:nt));
title('Z velocity (mm/s)');

hafp3 = subplot(3,5,3);
plot(plot_time, areas(1:nt));
title('Area (mm^2)');
hold on
plot(plot_time, areasxy(1:nt));
plot(plot_time, areasrightsideupxy(1:nt), '--');
Ylims = get(gca,'YLim');
set(gca,'Ylim',[ Ylims(1) Ylims(2)+0.3*(Ylims(2)-Ylims(1)) ] );
legend({'Surface area', 'Projected area', 'Upright projected' } );


hafp4 = subplot(3,5,4);
plot(plot_time, minstrains(1:nt));
title('Strain min/max');
hold on;
plot(plot_time, maxstrains(1:nt));
plot(plot_time([1 end]), [0 0], 'k');

hafp5 = subplot(3,5,5);
plot(plot_time, graphuos(1:nt));
title('Number of tensile failures');
set(gca,'YScale','log')

hafp6 = subplot(3,5,6 );
plot(plot_time, comxs(1:nt)-Ictrx);
title({'Position (mm)' '(COM rel to beam ctr)'});
hold on;
plot(plot_time, comys(1:nt));
plot(plot_time, comrs(1:nt), '--');
legend({'X' 'Y' 'R'});


hafp7 = subplot(3,5,7);
plot(plot_time, comzs(1:nt));
title('Z position (mm)');

hafp8 = subplot(3,5,8);
plot(plot_time, kinoenergies(1:nt));
title('Non-ballistic KE (J)');

hafp9 = subplot(3,5,9);
plot(plot_time, potenergies(1:nt));
title('Potential energy (J)');


%subplot(3,4,9);
%plot(plot_time, maxvs(1:nt));
%title('Max velocity rel. to COM');

hafp10 = subplot(3,5,10);
plot(plot_time, rad2deg( asin(mintrinormz(1:nt))));
title('Steepest triangle angle (\circ)');
hold on
plot(plot_time([1 end]), [0 0], 'k');



hafp11 = subplot(3,5,11);
plot(plot_time, ttlabspwr(1:nt));
title('Power absorbed / emitted');
hold on;
plot(plot_time, ttlradpwr(1:nt));

hafp12 = subplot(3,5,12);
plot(plot_time, acczs(1:nt));
hold on;
plot(plot_time, accts(1:nt));
title('Acceleration (mm/s/s)');
legend({'Z' 'Total'});

hafp13 = subplot(3,5,13);
plot(plot_time, avgts(1:nt));
hold on
plot(plot_time, vertextemps(1:nt));
plot(plot_time, edgepttemps(1:nt));
plot(plot_time, maxts(1:nt));
title('Temperature (\circ K)');
legend({ 'Avg' 'Vertex' 'Edge point' 'Peak' } );

hafp14 = subplot(3,5,14);
plot(plot_time(2:end), diff(avgts(1:nt))./dnto);
hold on
%plot(plot_time(2:end), diff(maxts(1:nt))./dnto);  % this one isn't really valid as the point index shifts all the time...
plot(plot_time(2:end), diff(vertextemps(1:nt))./dnto);
plot(plot_time(2:end), diff(edgepttemps(1:nt))./dnto);
title('dT/dt (\circK/s)');
plot(plot_time([1 end]), [0 0], 'k');
legend({ 'Avg' 'Vertex' 'Edge point' } );

hafp15 = subplot(3,5,15);
plot(graphts(1:nto_firstBroken-1), longaxis_length(1:nto_firstBroken-1));
hold on;
title({'Elongation' '[shape base diameter (mm)]'});
plot(graphts(1:nto_firstBroken-1), shortaxis_length(1:nto_firstBroken-1));
plot(plot_time([1 end]), 2*[radActual0 radActual0], 'k--');
legend({ 'Major dia' 'Minor dia' 'Rest dia'} );

linkaxes([hafp1, hafp2, hafp3, hafp4, hafp5, hafp6, hafp7, hafp8, hafp9, hafp10, hafp11, hafp12, hafp13, hafp14, hafp15 ], 'x');

saveas(gcf,['sim' num2str(movieNum) '.fig']);
saveas(gcf,['sim' num2str(movieNum) '.png']);

figure(fSnapshot);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'ZTick',[]);
set(gca,'Box','on');
plot3([Ictrx Ictrx], [0 0], get(gca,'ZLim'),'m');

saveas(gcf,['sim' num2str(movieNum) '_snapshot.fig']);
saveas(gcf,['sim' num2str(movieNum) '_snapshot.png']);

%figure('pos',[100 100 2100 900]);
figure;
set(gcf,'color','w');
hline = plot3(comxs(1:nto_firstBroken-1),comys(1:nto_firstBroken-1),plot_time(1:nto_firstBroken-1),'DisplayName',num2str(movieNum));
if (nto_firstBroken < nto)
    hold on;
    plot3(comxs(nto_firstBroken:nto),comys(nto_firstBroken:nto),plot_time(nto_firstBroken:nto),'Color',get(hline,'Color'),'LineStyle',':');
    hline2=plot3(comxs(nto_firstBroken),comys(nto_firstBroken),plot_time(nto_firstBroken),'*');
    set(hline2,'Color',get(hline,'Color'));
end
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
box on;
view(0,90);

saveas(gcf,['sim' num2str(movieNum) '_trajxy.fig']);
saveas(gcf,['sim' num2str(movieNum) '_trajxy.png']);

copyfile('simulate_02092021.m', ['sim' num2str(movieNum) '_sim.script']);
copyfile('generateMesh_v2.m', ['sim' num2str(movieNum) '_gen.script']);

disp('Data table entry:')
spacer = 0;
numnodes = length(n_x);
numedges = length(e_na);
numtris = length(t_na);
peakstrain = max(maxstrains);
peakt = max(maxts);
peakaccz = max(acczs);
peakpwr = max(ttlinpwr);
disptablerow( { 'movieNum',...
	'radius','x0','pD','totalarea','totalmass','sum(t_a0xy)','rot0', ...
    'numnodes','numedges','numtris',...
	'thickness','density','Ymod',...
	'Emod','Thermcondmm','heatCap',... %   'cte0','cte1','cts','ctw',  ...
    'spacer','spacer','spacer','spacer', ...
	'I0','Irad','Ictrx', ...
	'Iabs','Irefl','Emissivity', ...
	'spinup_gs','spinup_time','plotFast', ...
	't_start_ramp','t_ramp_dur','t_on_dur', ...
	'dt','tt','nt','frame_interval','snapshot_interval', ...
    'uhohs','comx','comy','comz', ...
    'meanvx','meanvy','meanvz', ...
    'peakt','peakaccz','peakpwr','peakstrain', ...
    'rate_oa', 'rate_loop', 'rate_rend', ...
    'spacer','time_firstBroken','et' }, ...
    {'meshtype', 'end_reason'} );
outerloopstable = [ outerloopstable ...
    getdisptablerow( { 'movieNum',...
	'radius','x0','pD','totalarea','totalmass','sum(t_a0xy)','rot0', ...
    'numnodes','numedges','numtris',...
	'thickness','density','Ymod',...
	'Emod','Thermcondmm','heatCap',... %   'cte0','cte1','cts','ctw',  ...
    'spacer','spacer','spacer','spacer', ...
	'I0','Irad','Ictrx', ...
	'Iabs','Irefl','Emissivity', ...
	'spinup_gs','spinup_time','plotFast', ...
	't_start_ramp','t_ramp_dur','t_on_dur', ...
	'dt','tt','nt','frame_interval','snapshot_interval', ...
    'uhohs','comx','comy','comz', ...
    'meanvx','meanvy','meanvz', ...
    'peakt','peakaccz','peakpwr','peakstrain', ...
    'rate_oa', 'rate_loop', 'rate_rend', ...
    'spacer','time_firstBroken','et' }, ...
    {'meshtype', 'end_reason'} ) ];

diary off;


%% end of outerloop1
if (usercancel)
    break;
end
end
if (usercancel)
    break;
end
end
disp(outerloopstable);
