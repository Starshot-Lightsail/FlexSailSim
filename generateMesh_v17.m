% Lightsail mesh generator, v17 (2021 05 03) 
% (c) 2019-2021 Michael D. Kelzenberg, Ramon Gao, California Institute of Technology
% 
% Generates mesh structures for use with our simulator code, for the purpose of simulating lightsails having various
% shapes, mechanical properties, optical properties, nanophotonic texturing, etc.  It can produce a handful of curved
% or flat open membrane shapes with general symmetry about the Z axis.  Topologically, all meshes must represent open 
% 2D surfaces, and furthermore the initial mesh must project into the XY plane without self-intersection or singularites
% in the projected area of each mesh triangle.  Thus for example, closed surfaces such as a complete sphere can't be 
% meshed.  At most a full hemisphere could be meshed, but even this seems inadvisable wihtout substantial improvements 
% to the optics simulation code to deal with multiple reflections and reverse-illuminated surfaces.  Ultimately the 
% simulation is currently suitable only for flat or modestly curved membrane shapes.  An implicit assumption is that 
% lightsail membranes are infinitesimally thin and thus have negligable (zero) out-of-plane bending stiffness.  If you 
% wanted to simulate a structure with out-of plane stiffness (other than a rigid body), you could implement appropriate 
% mechanics in the membrane simulation code, or alternatly, develop a 3D meshing topology to simulate multi-layer meshes
% while also altering the optical and thermal physics calculations to deal with the new topology.  
%
% Tip:  This generator script can be run on its own and will produce a plot of the resulting mesh.  This is
% the best way to figure out which settings to use.  Make sure to run the setupMaterialProperties script first.

%% Save a copy of the mesh generation script each time it's run...
% So we can figure out what settings we used when trying to make sense of the results later...
%warning('fixme:  Commented out important part of generate_mesh')
copyfile([ mfilename  '.m'], [filebasename '_gen.script']) % save a copy of this mesh generation script, under the 
 % current simulation ID (movieNum) via filebasename, so we know what mesh settings were used for this simulation!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Geometry Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radiusmm = 2*450;%710/2;  % (mm) the nominal radius (or inner half-width) of the lightsail.  NOTE:  due to the way I've 
  % written the node generators, this is NOT the actual radial extent of the mesh for non-flat structures!  Imagine 
  % taking a flat lightasil disc with the desired radius (radiusmm) and draping it over a mandrel with the desired 
  % 3D profile.  That is kind of how my node generators work.  Higher aspect ratio shapes will have reduced actual 
  % radius in the xy plane.  This is in attempt to produce uniform edge lengths and triangle size, regardless of 
  % aspect ratio.  

radialRings = 20   ; % Number of meshed ring layers spanning radially from the origin to the outer boundary of the lightsail.  
 % Increasing this parameter will increase the number of simulation mesh elements, reduce the size of the mesh
 % triangles, and increase overall fidelity of the simulations.  
 % MUST BE POSITIVE NONZERO INTEGER.  Reasonable values are 2-5 for rigid bodies, 20-50 for flexible bodies.  Values
 % above ~20 start to get really slow in simulations!
    radialRings = round(abs(radialRings));  %make sure it's an integer
    radius=radialRings;  %duplicate as 'radius' for legacy reasons.  Do not change.

zMode = 'paraboloid';  % Can be 'flat', 'paraboloid', 'cone', or 'sphere'.  Determines the z-profile 
                        % of the sail mesh.  I use terms for three-dimensional strucutres, but this setting more
                        % accurately describes the 2D cross section of the mesh in the principal vertical plane.(i.e.,
                        % flat line, parabola, sloped line, or circle).  Setting "curvatureMode" to 'smooth' will
                        % produce the described shape profiles.  
                       
zAspectRatio = 1.5;%.2;  % a relative measure of the height of the sail strucutre produced by the non-flat node generators. 
                   % This is ignored for 'flat' zMode.  A value of zero gives 'flat' zProfileMode structures even for 
                   % non-flat zProfileModes. Due to the convoluted and evolving ways I wrote the meshing code, this value 
                   % doesn't actually specify the numerical aspect ratio (height to width ratio) for any given strucutre,
                   % but for any given structure, a higher zAspectRatio will produce a higher aspect ratio for the
                   % generated mesh, hopefully in a way that scales sensibly with changes to the radiusmm and 
                   % radialRings values.  If you want a specific shape, you will have to study the code below to 
                   % figure out what the hell is going on.   Certain ranges of zAspectRatio are invalid for certain 
                   % zMode node generators, but most invalid combinations should at least produce an error message.  
                   % Good luck!
                        
curvatureMode = 'faceted'; % can be 'smooth' or 'faceted'.  Determines how the out-of-plane curvature is applied for 
                          % non-flat non-round sail geometries.  'smooth' curvature will produce axisymmetric meshes, i.e.
                          % surfaces of revolution, but will leave 'jagged' bottom edges for meshes with non-circular
                          % footprints.  'faceted' curvature will produce structures with discrete rotational 
                          % symmetry about the Z axis, but with constant Z-height for each mesh ring.  The profile of 
                          % each facet conforms to the desired zProfileMode contour only at the orthagonal radial plane
                          % which bisects the 'facet'.  
                          % Note:  This setting is effectively ignored for "round" xyMode and for "flat"
                          % zProfileMode, which do not allow for faceted geometry.

xyMode = 'round'; % Can be 'hex', 'round', or 'square'
 % this is essentially the in-plane node generation method, and determines the projected shape of the sail.   
 % 'hex': produces lightsails with hexagon footprints.  Hex meshing produces equilateral triangles
 %    throughout the mesh (except for distortions due to surface curvature) which yield the most faithfully isotropic
 %    membrane behavior in simulation.  
 % 'round':  produces membranes with circular footprints.  Mesh triangulation ranges from nearly equilateral along the 
 %    30/90/150 degree radials, to nearly 45-degree isosceles along 0/60/120 radials, and thus can produce 6-fold 
 %    symmetric artifacts in simulation due to underlying non-isotropic mesh cells. 
 % 'square':  Produces membranes with square footprints.  All mesh triangles are essentially 45-degree isosceles,
 %    so this is the least isotropic of the node placement methods.  
 %
 % NOTE:  the square mesher interprets the radius settings to specify the minor radius, i.e., half-width of the square.  
 %        the hex mesher interprets the radius settings to specify the major radius, i.e., vertex-to-center distance. 
 
z0 = 0;  % mm initial height of vertex or center of the sail
vz0 = 0; % mm/sec initial z velocity
t0 = 300; % K initial temperature

psi0 = deg2rad(30);  % Rotates the whole mesh by this amount (around Z) after node generation, prior to texture/region 
                    % assignment.  This is useful so we can rotate round meshes by 30 degrees so the triangles align 
                    % nicely with Ramon's nanophotonic design as of early April 2021, which has 60 degree wide pizza slices
                    % in a bowtie configuration along the X axis, so setting psi0 to 30 degrees prevents jagged
                    % texture boundaries.
                    
rescaleRadius = 1;  % whether or not to rescale the mesh to achieve the desired value of radiusmm.  For legacy reasons, 
                    % non-flat meshes will have their radius reduced to account for the curvature of the surface,
                    % seeking to produce somewhat consistent edge lengths throughout the mesh.  The actual radius (or
                    % more accurately, the outer extent of the mesh along the x axis) will be stored as radActualx0mm.
                    % Setting rescaleRadius will re-scale the n_x, n_y, and n_z points prior to meshing so that 
                    % radActualx0mm is equal to radiusmm.  
                    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
%% Texture (region) mapper settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Texture' is used to identify the specific optical properties of each
% triangle comprising the lightsail.  By assigning a texture to each 
% triangle, we can simulate complex lightsail designs to probe concepts 
% such as self-stabilization via metasurface/nanophotonic patterning.
%
% Assigning texture is analogous to picking patches of fabric for a quilt.
% Thanks, grandma!
%
% The array t_tex stores the texture value for each triangle.  The texture
% value is a non-negative integer.  Each value describes a different
% type of optical behavior. As of this writing, there are two classes of
% texture supported:
%
%  - a texture value of zero means the simulator should apply simple
%    specular optical physics to calculate that triangle's response to
%    incident light (via scalar constants for reflection and absorption in
%    the main simulator)
%
%  - a texture value 1,2,... specifies that a lookup table (LUT) is to be 
%    used to determine optical forces.  The number specifies which lookup 
%    table to use.  The simulator must be separately configured with the 
%    correct lookup tables for the desired lightsail materials.  
%
% Also defined in this section is the "orientation" of each texture, via
% the value t_texa, which stores the angle between the triangle's first 
% edge and the texture's primary axis.  In the code below, the vector
% this_tex_dir can be used to specify the texture's primary axis, from
% which the correct value of t_texa will automatically be calculated.
% this_tex_dir is projected into the plane of the triangle and normalized,
% so it needn't necessarily be a unit vector or lie within the triangle
% plane to start with, but obviously the code won't work if this_tex_dir is
% parallel to the normal vector, and projecting out-of-plane vectors onto
% curved surfaces (e.g., projecting 'verti' orientation on paraboloids) 
% gives slightly distorted texture directions when viewed from above.  
%
% For convinence, a handful of commonly used texture orientations are 
% defined here, such as 'horiz' (along x axis) and 'verti' (along y axis).
%
% Note that the same texture (lookup table) can be used in multiple
% orientations througout the lightsail.  For example, if you're simulating
% a mirror symmetric membrane with two regions, you can use the same
% texture lookup table for both halves -- just use an opposite texture
% direction vector for the left vs. right half.


%% texture mapping mode:  
% I have implemented three general-purpose texture mappers, or you can specify your own.
% 'pizza': divides the mesh into radial slices centered around the origin.
% 'stuffedcrust':  works just like the pizza mapper, but each slice can
%    also be divided into two texture regions, split at a certain radius. Different
%    crust radii, crust texture, and crust texture orientation, can be specified for 
%    each slice.
% 'brownies':  divides the mesh into regions split by horizontal and
%    vertical dividing lines.  
% 'user':  Map the regions yourself, near the bottom of this script.
texture_map_modes = { 'pizza' };  %must be a cell array of strings
% V17 update:  texture_map_mode is now a cell array of strings.  This lets us apply consecutive texture mappers to our
% mesh.  Mappers are applied consecutively, so the last one listed will have the final word.  To leave the texture
% assignment from a prior mapper, specify a texture value of -1 in the mapper settings for that region.  Currently, 
% each texture mapper has only a single set of settings, so to enable two passes of the same mapper, define a second 
% version of that mapper mode (e.g., 'stuffedcrust2') and define settings for that new mode (e.g., pizza_crust_texs2, 
% pizza_crust_texs_dirs2, pizza_crust_texs_radii2), then duplicate duplicate the mapping iterator code to support that 
% texture mode as well.  Then, set the texture map mode to {'stuffedcrust', ... 'stuffedcrust2'}


%% texture orientation convience vectors (do not change)
% The following vectors are defined for convienice of specifying cardinal directions and 30/45-degree incremented 
% directions when defining texture orientation, plus a few specialty kludgy vector flags.  Resist the urge to change the
% pre-defined vectors, but feel free to add your own here.
horiz = [1 0 0]';  % the x axis
verti = [0 1 0]';  % the y axis
vec0=horiz;
vec30 = [sqrt(3)/2 0.5 0]';
vec45 = [sqrt(2)/2 sqrt(2)/2 0]';
vec60 = [0.5 sqrt(3)/2 0]'; % 60� between diagonal and x > 0 axis
vec90=verti;
vec120 = [-0.5 sqrt(3)/2 0]'; % 120� between diagonal and x > 0 axis
vec135 = [-sqrt(2)/2 sqrt(2)/2 0]';
vec150 = [-sqrt(3)/2 0.5 0]';
%The above vectors can be combined via arithmetic operations in specifying 
%the texture orientations, e.g., sqrt(2)*(horiz+verti).  (Or just use vec45).
%However, the below values are special flags that can't be manipulated except
%for negation:
vec_radial = [0 0 123]'; % this is a special (kludgy) flag vector.
%  It will be replaced with the local radial direction vector
 % during the mapping.  Its actual value (123z) is not used -- it's
 % just a flag value for convience.  You can use -vec_radial to specify an
 % inward-pointing direction, but no other manipulation will work, since
 % it's just a flag, not a direction vector.  
vec_tang = [0 0 456]'; % this is a special (kludgey) flag vector. 
 % It will be replaced with the local tangential/circumferential 
 % direction vector during mapping (facing clockwise).
 % Its actual value (456z) is not used -- it's just a flag value for 
 % convience.  You can use -vec_tang to specify the clockwise direction, 
 % but no other manipulation will work, since it's just a flag, not an
 % actual direction vector.  

 
%% PIZZA slice texture mapper
% One way of assigning texture regions is using a pizza slicer approach.  
% Any number of regions can be defined by their radial position (angle)
% relative to the origin.  
pizza_angles = deg2rad([0 45 135 225 315 360.1]); %boundary angles of each slice. 
%pizza_angles = deg2rad([0 30 150 210 330 360.1]); %boundary angles of each slice. 
%pizza_angles = deg2rad([0 90 270 360.1]); %boundary angles of each slice. 
 % this can be arbitrarily long, but it should be monotonic, start with 0,
 % and end at 2pi_delta (I use 360.1 degrees, the extra 0.1 is because sometime
 % floating point errors happen).  Specifying n pizza angles will define
 % n-1 regions.  
pizza_texs = [ 1 2 1 2 1 ]+2;  % tex (LUT index) values for each pizza slice.
 % this should be one element shorter than pizza_angles.
pizza_dirs = [ verti -horiz -verti horiz  verti]; % the texture
 % directions for each pizza slice.  (these are vectors.) 

% % alternate set of pizza settings for quick change to specular behavior 
% if outersweepval1 == 2
%     pizza_texs = [1 1 1 1 1 ]+2;
% end
% if outersweepval1 == 3
%     pizza_texs = [1 2 1 2 1 ];
% end
% if outersweepval1 == 4
%     pizza_texs = [1 1 1 1 1];
% end
pizza_angles = deg2rad([0 360.1]);
pizza_texs = [ 0 ];
pizza_dirs = [vec_tang];
 
 
%% STUFFED CRUST texture mapper
% uses the pizza mapper values above, but also lets each slice be split
% into an inner and an outer region.  The outer region settings are:
pizza_crust_texs = [ 3 4 3 4 3 ];  %texture assignments for outer crust of each slice.  Must match size of pizza_texs;
pizza_crust_dirs = [ vec_radial vec_tang vec_tang vec_radial vec_tang vec_tang vec_radial]; % must match dims of pizza_dirs
pizza_crust_radii = radiusmm*[ .4 .8 .4 .8 .4 .8 .4 ];   % mm?  

% this alternate set of definitions converts a pizza slize mapping to have a non-textured inner region:
% pizza_crust_texs = pizza_texs; % copy inner region texture assignments to outer region
% pizza_texs = zeros(size(pizza_texs));  % replace inner region texture assignments with zeros (specular)
% pizza_crust_dirs = pizza_dirs; % copy inner region texture orientations to outer region
% pizza_crust_radii = 0.5 * radiusmm * ones(size(pizza_texs));   % set value for region boundary radius


%% BROWNIES texture mapper
% assigns textures by splitting the mesh along a grid of horizontal and
% vertical lines.
brownies_slices_x = [ -0.5*radiusmm 0 0.5*radiusmm ];
brownies_slices_y = [ -0.6*radiusmm -0.3*radiusmm 0 0.3*radiusmm 0.6*radiusmm ];
%if brownies_slices_x contains m values, and brownies_slices_y n, then the
%brownies_texs matrix should be (m+1) cols   x   (n+1) rows  in size.  
brownies_texs = [ 1 2 1 2 ; 2 1 2 1 ; ...
                  0 3 4 0 ; 0 4 3 0 ; ...
                  1 2 1 2 ; 2 1 2 1 ]';
%similarly, brownies_dirs should be a (n+1)-element cell array.  Each
%element is a (m+1) long array of vectors.
brownies_dirs = { [ verti horiz -verti -horiz ], [ vec60 vec120 -vec120 -vec60], ...
                  [ vec_radial vec_radial -vec_radial -vec_radial ], [vec_tang vec_tang -vec_tang -vec_tang], ...
                  [ vec45 vec135 -vec45 -vec135 ] , [ vec30 vec150 -vec150 -vec30] };
 
 
% This was my attempt to create oggy's tiling on the fly during a call with Ramon:              
%brownies_slices_x = x0*radialRings*[ -.7 0 .7 ] ;
%brownies_slices_y = x0*radialRings*[0];

%brownies_texs = [ [1 2 2 1] ; [ 1 2 2 1 ] ]';
%brownies_dirs = { [-verti -vec150 vec30 verti ] , [-verti -vec30 -vec150 verti] };


% this set of brownie texturing can be used to define a '+' specular region atop a prior pizza mapping, if specified as
% the final course of the meal.
brownies_slices_x = radiusmm .* [-0.1 0.1];
brownies_slices_y = radiusmm .* [-0.1 0.1];
brownies_texs = [ -1 0 -1; ...
                  0  0  0; ...
                  -1 0 -1; ];
brownies_dirs = { [verti verti verti], [verti verti verti], [verti verti verti] };
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rigid bodies (future feature -- not yet implemented!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lets us define certain nodes as belonging to rigid bodies, approximating support structures or payloads for the 
% lightsail.  I was working on this, but had to switch gears for a while.  I'm leaving the framework up here for now.

nrb = 0;  % number of rigid bodies

rb_types = { };  % rigid body types for each rigid body (strings).  
% I plan to support the following types of rigid bodies:
%  'hoop':  the outermost nodes of the structure.  Params = [width (mm)].  Width from outer radius.
%  'disc':  The innermost nodes of the structure, defined by a cutoff radius.  Params = [radius (mm)]
%  'cross':  Nodes lying within a certain distance of two orthogonal axes.  Params = [ halfwidth (mm), phi (rad)]
%  '5pt':  Similar to 'cross', but only includes nodes at the center and ends of each arm.  Same params. 
%  '4pt':  Similar to '5pt', but lacking the center region.  Same params.

rb_params = { }; % extra params for rigid body definition.

rb_mass = [ ]; % mass of each rigid body.  I will figure out how to distribute the mass and calculate the intertia 
   % tensors based on this value.  
   
rb_texoverride = [ ]; % overrides texture assignment for triangles included in the rigid body footprint.  Set to -1 to 
                      % leave original texture assignment.  If there's no lightsail in the rigid body's footprint,
                      % configure a LUT with all-zero values and assign its value here. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of commonly changed user settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
%% calculated values, do not change
x0 = radiusmm / radialRings; % the nominal edge length between nodes (mm).  This is calculated; do not change.  
                             % My node generators generally attempt to produce consistent edge length throughout the array
                             % by applying first-order corrections for the curvature of the surfaces.  However,
                             % only flat hex meshes are perfectly isotropic (i.e. all triangles are equilateral).
                             
% Calculate parameters for zMode
if isequal(zMode, 'flat')
    pD = 0;
    zAspectRatio = 0;
else
    if zAspectRatio == 0
        zMode = 'flat';
        warning('Changed mesh generator shape zMode to ''flat'' because zAspectRatio is zero!');
        pD = 0;
    else
        pD = abs(radialRings/zAspectRatio); % parabolic diameter @ latus rectum term, larger values = less curvature, ...
    end                               
end

disp('Sharshot lightsail simulator mesh generator version 17 (20210503) (c) 2019-2021 Michael Kelzenberg, Ramon Gao, California Institute of Technology');
fprintf('%s\n',datestr(now)); % Date & time
dispvar = @(varname1                     ) evalin('base',['fprintf('' %+12s: %-12g\n'', ''' varname1 ''', ' varname1 ');'] );
dispvar2 = @(varname1, varname2          ) evalin('base',['fprintf('' %+12s: %-12g  %+12s: %-12g\n'', ''' varname1 ''', ' varname1 ', ''' varname2 ''', ' varname2 ');'] );
dispvar3 = @(varname1, varname2, varname3) evalin('base',['fprintf('' %+12s: %-12g  %+12s: %-12g  %+12s: %-12g\n'', ''' varname1 ''', ' varname1 ', ''' varname2 ''', ' varname2 ', ''' varname3 ''', ' varname3 ');'] );
dispvar3s = @(varname1, varname2, varname3) evalin('base',['fprintf('' %+12s: %-12s  %+12s: %-12s  %+12s: %-12s\n'', ''' varname1 ''', ' varname1 ', ''' varname2 ''', ' varname2 ', ''' varname3 ''', ' varname3 ');'] );
dispvar3s('zMode','curvatureMode','xyMode');
dispvar3('radiusmm','radialRings','zAspectRatio');
dispvar3('z0', 'vz0', 'psi0' );
dispvar3('thicknessnm','density','Ymod');
dispvar3('Emod','Thermcondmm','heatCap');

%% error checking
if ~( strcmp(zMode,'flat') || strcmp(zMode,'paraboloid') || strcmp(zMode,'cone') || strcmp(zMode,'sphere') )
    error(['Invalid zMode value "' zMode '"   Valid values are "flat" "paraboloid" "cone" or "sphere."']);
end
if ~( strcmp(xyMode,'hex') || strcmp(xyMode,'round') || strcmp(xyMode,'square') )
    error(['Invalid xyMode value "' xyMode '"   Valid values are "hex" "round" or "square"']);
end
if ~( strcmp(curvatureMode,'smooth') || strcmp(curvatureMode,'faceted') )
    error(['Invalid curvatureMode value "' curvatureMode '"   Valid values are "smooth" or "faceted"']);
end
if radialRings <= 0
    error('radialRings must be a positive integer.  Try again.');
end
if radialRings > 50
    error('radialRings is absurdely large.');
end
flipzslater = 0;
if zAspectRatio<0
    flipzslater=1;
    zAspectRatio=abs(zAspectRatio);
end
if strcmp(zMode,'sphere') && (zAspectRatio > 0.80)
    error(['The Z aspect ratio for sphere zMode can''t exceed 0.80 (was ' num2str(zAspectRatio) ')']);
end
if strcmp(zMode,'sphere') && (zAspectRatio > 0.60)
    warning(['Z aspect ratios above ~0.6 for spherical sails produce weird triangulation near the base. (val=' num2str(zAspectRatio) ')']);
end
if strcmp(zMode,'cone') && (zAspectRatio > 1)
    %warning(['Z aspect ratios above 1.0 for spherical sails are not accurately modeled by this simulator -- they require multiple-reflection optical calculations. (val=' num2str(zAspectRatio) ')']);
end

%% enums / preprocessor stuff
region_main = 0; % 0 indicates that node is in the interior of the mesh
region_hoop = 1; % 1 indicates that node is on the edge of the mesh
region_payload = 2;
region_ctr = 3;
region_test = 4;

MAX_EDGES_PER_NODE = 10; % what is this useful for? Seems like number of edges per node will be at most 6
    %MK:  Yes, this is for memory allocation purposes only.   For all meshing strategies currently employed, the max 
    % number of edges should be 6.  But the future is hard to predict, so let's arbitrarily add a bit of margin to the
    % buffer size.  If we're lucky, it will confuse someone in the future. 

global n_x n_y n_z e_na e_nb e_l0 t_na t_nb t_nc
% since we're working in the base workspace, making these global isn't necessary, but it makes them display in a
% different font color in my editor, which helps my brain sometimes.  

%% Node generation
% Here are the node-based mesh variables:
%  ea:  (const) edges array (reference to edges) (zero padded array)
%  ta:  (const) triangles array (reference to triangles) (zero padded array)
%  nt:  (const) number of triangles bordering the node
%  ne:  (const) number of edges connected to the node
% Here are the key static node creation-frame variables
%  r0 a0 x0...z0:  original radius, polar position (rad), and xyz position in the resting structure, prior to any transformations
%  m:  mass
% Here are some of the key node simulation variables (dynamics state variables):
%  x y z:  actual position, mm (scalar)
%  vx vy vz:  actual velocity, mm/s (scalar)
%  t    temperature (scalar), K, currently set to t0 everywhere

%
% We start by generating the x,y,z coordinates of each mesh vertex (a.k.a. nodes)
% n_x, n_y, and n_z comprise one set of state vectors for our simulation mesh. 
% Currently, our mesher only works for open convex 2D surface membranes from a topological standpoint.  In other words, 
% simulatable membranes are things you could imagine crafting from a single piece of of foil or seran wrap without 
% punching holes or taping edges together.  We further requre that all simulation meshes project into the XY plane 
% without self-intersection.  Of in other words, that we can see the entire top side of our lightsail surface when 
% viewing from directly above, along the Z axis.  Once again, all membrane meshes must represent 2D toplogocial 
% surfaces; out-of-plane connectivity is vollstaendigly verboten.  
%
% Presently all node generators follow an "onion peel" or stepped outward spiral node placement algorithm.  In other
% words, the first node is placed at the origin, then additional sets of nodes are placed in a circle/hexagonal/square
% rings surrounding the origin, and so forth until the desired number of concentric node rings (radialRings) is reached. 
%
% Note:  The node generation code is in sad need of being refactored due to piecemeal 

disp('Generating nodes...')
tic;  % for timing metrics

%sphereical surface calcs:
    if strcmp(zMode,'sphere')
        sphrad = radiusmm / (pi/2 - asin(1 - zAspectRatio));
        sphdalpha = (pi/2 - asin(1-zAspectRatio))/radialRings;
    end
%conical surface calcs:
    if strcmp(zMode,'cone')
        conex0 = x0 / sqrt( 1 + zAspectRatio.^2 );
    end

if strcmp(xyMode,'hex') || strcmp(xyMode,'round')   %  6x nodes per ring
    
    ringttl = @(x) 1 + 3.*x + 3.*x.^2; % function to calculate number of nodes for this meshing strategy, with x being number of concetric rings. 
      %Ramon's notes:  Note that the code currently only works with this function, i.e. changing anything about it will yield errors - why?
      %MK answer:  Because this function is used throughout the algorithm to determine the node indices of the first and last
      % node in each concentric ring.  Changing the function will yield incorrect indices while building the mesh, and
      % is likely to result in index-exceeds-bounds type errors, and will certainly result in invalid/unintended mesh
      % structures.  
    numnodes = ringttl(radialRings); % total number of nodes for this mesh
    
%     constrainedEdges = [ ( (ringttl(radius-1)+1):(ringttl(radius)-1) )' ... MK's version
%         ( (ringttl(radius-1)+2):(ringttl(radius)) )' ; ...
%          (ringttl(radius)) (ringttl(radius-1)+1);];
    
    constrainedEdges = [ ( (ringttl(radius-1)+1):(ringttl(radius)-1) )' ...
        ( (ringttl(radius-1)+2):(ringttl(radius)) )' ; ...
         (ringttl(radius)) (ringttl(radius-1)+1);];
    %^^ this is a list of edges (defined as pairs of nodes) comprising the outer perimeter of the mesh.  It is required
    %to know which edges form the perimeter during the meshing (edge creation) step, so we can delete superfluous
    %perimeter triangles.
    
    
    n_x = zeros(1,numnodes); % initialize vector for x-coordinate of each node
    n_y = zeros(1,numnodes); % initialize vector for y-coordinate of each node
    n_z = zeros(1,numnodes); % initialize vector for z-coordinate of each node 
    n_ang0 = zeros(1,numnodes); %initialize vector for initial angular position of each node (i.e., polar coordinates)
    n_r0 = zeros(1,numnodes); %initialize vector for initial radial position of each node (i.e., polar coordinates)
    n_reg = zeros(1,numnodes); %initialize vector for region assignment for each node (this is NOT the same as texture assignment and is largely no longer in use!)

    hs = x0; % nominal horizontal spacing between nodes for hex grid
    hu = hs*sqrt(3)/2; % nominal vertical spacing between nodes for hex grid
    ups = [ hu 0 -hu -hu 0 hu ]; % vertical (y) displacements between consecutive nodes within each hextant for our algorithm.
    rights = [ -hs/2 -hs -hs/2 hs/2 hs hs/2 ]; % horizontal (x) displacements between consecutive nodes within each hextant for our algorithm. 
    
    np = 1;  %node (point) index for iterative construction of the node position arrays
    n_z(np) = z0; % assign initial height of first node (the vertex node).  (Its X and Y values are always zero.)
    
    myxx = 0; %radial x-axis position (mm) for node construction
    
    
    
    for nr = 1:radialRings % sweep over concentric rings
        % increment radial x position for next concentric ring of nodes, reducing x step to account for surface curvature:
        if ( ( pD == 0)) % flat surface, no correction
            myxx = myxx + x0;
        else % non-flat surface, reduce x step according to local steepness of the curve, so as to preserve approximately
            if strcmp(zMode,'paraboloid')                                            % uniform edge length througout the mesh.
                myxx = myxx + x0^2 * pD / ( (pD*x0)^2 + 2 * (myxx)^2).^0.5;        
            elseif strcmp(zMode,'cone')
                myxx = myxx +  conex0;
            elseif strcmp(zMode,'sphere')
                myxx = sphrad * sin( nr * sphdalpha ); %oof
            else
                error('Go fuck yourself, San Diego!');
            end
        end
                
        myx = myxx;
        mynormx = myx / x0;
        
        % set z-height for this ring of nodes (faceted curve mode)
        if ( (pD == 0))  % flat surface, no change to z height
            myz = z0;
        else
            if strcmp(zMode,'paraboloid') % parabolic surface, change z height accordingly
                myz = z0 - (x0/pD)*(mynormx^2); %I can't believe this works.  Should probably be redfined as a function of nr for sanity...
            elseif strcmp(zMode,'cone')
                myz = z0 - nr * conex0 * zAspectRatio;
            elseif strcmp(zMode,'sphere')
                myz = z0 - sphrad*( 1 - cos(nr*sphdalpha) );
            end       
        end
        
        if strcmp(xyMode,'hex')                              % for hex meshes, we'll adjust the vertical and horizontal
            myups = ups * mynormx / nr;                      % node position offsets in step with the already-adjusted 
            myrights = rights * mynormx / nr;                % x0 offset to account for surface curvature, then populate
                                                             % the node ring by following the vertical ("ups") and
                                                             % horizontal ("rights") node offset vectors in each hextant
            myy = 0; % Starting at the vertex along the x axis...  first point will be in the first quadrant... last point will be here
            
            for ns = 1:6  % there are 6 hextants 
                for nh = 1:nr  % the number of "hops" per hextant grows with our radial index...
                    %increment the local node position appropriately for this hextant:
                    myx = myx + myrights(ns);
                    myy = myy + myups(ns);
                    myr = sqrt(myx^2 + myy^2);
                    %calculate the angular position of this node (i.e., for polar coordinates)
                    mya = atan(myy/myx);
                    if myx < 0  %correct atan ouput to give actual angle in [0 2pi] range
                        mya = mya+pi;
                    elseif myy < 0
                        mya = mya + 2*pi;
                    end
                    mya2 = atan2(myy,myx);  % todo:  check if this gives the same results!

                    np = np + 1;     % increment to next point
                    n_x(np) = myx;   % assign next point x
                    n_y(np) = myy;   % assign next point y
                    if strcmp(curvatureMode,'smooth')  % assign next point z for smooth vs. faceted paraboloid
                        if (~pD) 
                            n_z(np) = z0;
                        else
                            if strcmp(zMode,'paraboloid')
                                n_z(np) = z0-((myx/x0)^2+(myy/x0)^2)/pD*x0    ;%   +0.1*sin(2*pi*nh/nr)*nr;
                            elseif strcmp(zMode,'cone')
                                n_z(np) = z0 - zAspectRatio*myr;
                            elseif strcmp(zMode,'sphere')
                                n_z(np) = z0 - sphrad*(1-sqrt(1-(myr/sphrad)^2));
                            end
                        end
                    else  %faceted
                        n_z(np) = myz                     ;%   +0.1*sin(2*pi*nh/nr)*nr;
                    end
                    n_r0(np) = myr;    % retain initial radial position for this node
                    n_ang0(np) = mya;    % retain initial angular position for this node
                    if (nr == radialRings)  %assign n_reg value, I forgot why I do this.
                        n_reg(np) = region_hoop;
                    else
                        n_reg(np) = region_main;
                    end
                end
            end
        else %'round' mesh; similar approach, but use sin() and cos() to assign positions instead of stepping along six hexagonal directions.  
            for na = 0:(6*nr-1)
                mya = 2*pi*(na)/(6*nr);
                myr = myx;
                np = np + 1;
                n_x(np) = myx * cos(mya);
                n_y(np) = myx * sin(mya);
                %n_x(np) = (myx * cos(mya)*cos(psi) - myx * sin(mya)*sin(psi));
                %n_y(np) = (myx * cos(mya)*sin(psi) + myx * sin(mya)*cos(psi));
                n_z(np) = myz;
                n_r0(np) = myx;
                n_ang0(np) = mya;
                
                if (nr == radialRings)
                    n_reg(np) = region_hoop;
                else
                    n_reg(np) = region_main;
                end
            end
        end
    end
    %% Set up cross section for round meshes
    %  This is so I can easily plot any mesh property along an initially defined cross section of the mesh...
    if strcmp(xyMode,'round')
        mycs1 = ringttl(0:radialRings-1) + 1;%round(1.5*(1:radius-1)+1);
        mycs2 = ringttl(0:radialRings-1) + 3*(0:radialRings-1) + 3 + 1; %round(4.5*(1:radius-1)+3)+1;
        meshcs = [ fliplr(mycs2) 1 mycs1 ];  % a list of node indices falling on the initial horizontal axis cross section of the mesh
        
        % these are the node indices of the outermost nodes along the x and y axes.  I might use these later to determine
        % the approximate pitch and orientation of the membrane.  
        xradidx = ringttl(radialRings-1)+1;    
        yradidx = ringttl(radialRings-1)+floor(1.5*radialRings)+1;
        nxradidx = ringttl(radialRings-1)+3*radialRings+1;
        nyradidx = ringttl(radialRings-1)+floor(4.5*radialRings)+1;
    end
    %% Set up cross section for hex meshes
    if strcmp(xyMode,'hex')
        mycs1 = ringttl(1:radialRings) + 0;%round(1.5*(1:radius-1)+1);
        mycs2 = ringttl(0:radialRings-1) + 3*(0:radialRings-1) + 3 + 0; %round(4.5*(1:radius-1)+3)+1;
        meshcs = [ fliplr(mycs2) 1  mycs1 ];  % a list of node indices falling on the initial horizontal axis cross section of the mesh
            
        % these are the node indices of the outermost nodes along the x and y axes.  I might use these later to determine
        % the approximate pitch and orientation of the membrane.  
        xradidx = ringttl(radialRings);    
        yradidx = ringttl(radialRings-1)+floor(1.5*radialRings);
        nxradidx = ringttl(radialRings-1)+3*radialRings;
        nyradidx = ringttl(radialRings-1)+floor(4.5*radialRings);
    end
    
        

elseif strcmp(xyMode,'square')  % square lightsail
        
    ringttl = @(x) (2*x + 1).^2; % function to calculate number of nodes with x being number of concetric rings. Note that this is different than the hex/round version.
    numnodes = ringttl(radialRings); % number of nodes in this mesh
    constrainedEdges = [ ((ringttl(radialRings-1)+1):(ringttl(radialRings)-1)) ...
        ringttl(radialRings) ; ...
        ((ringttl(radialRings-1)+2):(ringttl(radialRings))) ringttl(radialRings-1)+1]';
    %^^ this is a list of edges (defined as pairs of nodes) comprising the outer perimeter of the mesh.  It is required
    %to know which edges form the perimeter during the meshing (edge creation) step, so we can delete superfluous
    %perimeter triangles.
    

    
    n_x = zeros(1,numnodes); % initialize vector for x-coordinate of each node
    n_y = zeros(1,numnodes); % initialize vector for y-coordinate of each node
    n_z = zeros(1,numnodes); % initialize vector for z-coordinate of each node 
    n_ang0 = zeros(1,numnodes); %initialize vector for initial angular position of each node (i.e., polar coordinates)
    n_r0 = zeros(1,numnodes); %initialize vector for initial radial position of each node (i.e., polar coordinates)
    n_reg = zeros(1,numnodes); %initialize vector specifying whether node is in interior region (0) or on the edge (1)
    %n_reg is not the same as the texture/region assignment for MEPS.  It
    %is mostly left over from a previous version where I was varying the
    %mass density in various regions of the mesh.
   
    ups = x0.*[ 1 0 -1 0 ];  % same approach as hex node generator, but with only 4 sections to traverse.
    rights = x0.*[ 0 -1 0 1 ];
    
    np = 1;
    n_z(np) = z0; % assign initial height of vertex to first element of this vector
    
    myxx = 0;
    
    for nr = 1:radialRings % sweep over concentric rings
        % get ring seed position (myxx), adjusting lateral grid density to account for surface curvature :
        if ( ( pD == 0))
            myxx = myxx + x0;
        else %non-flat surface, adjust ring seed X position for surface curvature
            if strcmp(zMode,'paraboloid')
                myxx = myxx + x0^2 * pD / ( (pD*x0)^2 + 2 * (myxx)^2).^0.5;
            elseif strcmp(zMode,'cone')
                myxx = myxx +  conex0;
            elseif strcmp(zMode,'sphere')
                myxx = sphrad * sin( nr * sphdalpha ); 
            else
                error('Go fuck yourself, San Diego!');
            end
        end
        
        myx = myxx;
        mynormx = myx / x0;
        
        % get ring seed z height:
        if ( (pD == 0))  % flat surface, no change to z height
            myz = z0;
        else
            if strcmp(zMode,'paraboloid') 
                myz = z0 - (x0/pD)*(mynormx^2); 
            elseif strcmp(zMode,'cone')
                myz = z0 - nr * conex0 * zAspectRatio;
            elseif strcmp(zMode,'sphere')
                myz = z0 - sphrad*( 1 - cos(nr*sphdalpha) );
            end
        end
        
            myups = ups * mynormx / nr;
            myrights = rights * mynormx / nr;
            myy = -myx;% * last point will be here, in the lower right hand corner.  First point will be just above it.
            
            for ns = 1:4
                for nh = 1:2*nr
                
                    myx = myx + myrights(ns);
                    myy = myy + myups(ns);
                    mya = atan(myy/myx);
                    if myx < 0
                        mya = mya+pi;
                    elseif myy < 0
                        mya = mya + 2*pi;
                    end
                    myr = sqrt(myx^2 + myy^2);

                    np = np + 1;
                    n_x(np) = myx;
                    n_y(np) = myy;
                    %testing only -- disable these two lines
                    %  n_x(np) = myx                          + cos(mya)*nr/radius/5;
                    %  n_y(np) = myy                          + sin(mya)*nr/radius/5;
                    if strcmp(curvatureMode,'smooth')  % smooth paraboloid
                        if (~pD) 
                            n_z(np) = z0;
                        else
                            if strcmp(zMode,'paraboloid')
                                n_z(np) = z0-((myx/x0)^2+(myy/x0)^2)/pD*x0    ;%   +0.1*sin(2*pi*nh/nr)*nr;
                            elseif strcmp(zMode,'cone')
                                n_z(np) = z0 - zAspectRatio*myr;
                            elseif strcmp(zMode,'sphere')
                                n_z(np) = z0 - sphrad*(1-sqrt(1-(myr/1.42/sphrad)^2)); %z0 - sphrad*(1-cos(nr*sphdalpha));           
                            end
                        end
                    else  % faceted 
                        n_z(np) = myz                     ;%   +0.1*sin(2*pi*nh/nr)*nr;
                    end
                    n_r0(np) = sqrt(myx^2 + myy^2);
                    n_ang0(np) = mya;
                    if (nr == radialRings)
                        n_reg(np) = region_hoop;
                    else
                        n_reg(np) = region_main;
                    end
                end
            end
    end
    %% Set up cross section for square
    mycs1 = ringttl(1:radialRings) ; %round(1.5*(1:radius-1)+1);
    mycs2 = ringttl(1:radialRings) - 4*(1:radialRings) + 0;%round(4.5*(1:radius-1)+3)+1;
    meshcs = [ fliplr(mycs2) 1 mycs1 ];  % a list of node indices falling on the initial horizontal axis cross section of the mesh
    
    % these are the node indices of the four outermost nodes lying along the x and y axes (min and max).  I might use 
    % these later to determine the approximate pitch or orientation of the membrane.  
    xradidx = ringttl(radialRings);
    yradidx = ringttl(radialRings-1)+2*radialRings;
    nxradidx = ringttl(radialRings-1)+4*radialRings;
    nyradidx = ringttl(radialRings-1)+6*radialRings;
else
    error('Unrecognized mesh type');
end

if flipzslater
    n_z = z0 - (n_z-z0);
end

radActualx0mm = max(n_x);    % actual xy-plane radius of mesh, along x axis, major or minor depending on shape type...

if rescaleRadius
    if radActualx0mm ~= radiusmm
        n_x = n_x * radiusmm/radActualx0mm;
        n_y = n_y * radiusmm/radActualx0mm;
        n_z = (n_z-z0)*radiusmm/radActualx0mm + z0;
        radActualx0mm = max(n_x);
    end
end

%%Apply initial rotation (this is for static positioning, not spin stabilizing)
if psi0
    n_x_new = cos(psi0).*n_x - sin(psi0).*n_y;
    n_y = sin(psi0).*n_x + cos(psi0).*n_y;
    n_x = n_x_new;
    clear('n_x_new');
end



outerringstartidx = ringttl(radialRings-1)+1; % index of first node in outermost ring
outerringendidx = ringttl(radialRings); % index of last node in outermost ring
    % Mistake? Missing a " +1 "   %MK: probably?  your version looks right.
numringpts = outerringendidx - outerringstartidx + 1;


%  NOTE:  n_a0 is used here as the initial angular position of the node,
%  for calculating initial velocity of each node in the rotating body.
%  Later, it will be re-purposed to mean the initial "area" of the node, 
%  which is used in calculating radiative cooling.
    
% initial velocity , temperature, and other properties
n_vz = (zeros(1,numnodes)) + vz0; % accounting for initial velocity 
n_vx = (zeros(1,numnodes));%-pi*rot0*n_r0.*sin(n_a0); % initial velocity in x direction
n_vy = (zeros(1,numnodes));%pi*rot0*n_r0.*cos(n_a0); % initial velocity in y direction
  %Note:  it doesn't work to set vx and vy here for the initial rotational
  %velocity of the sail, because the mesh is created at rest.  Subjecting 
  %it to full rotational velocity instantaneously is unrealistic, and the 
  %resulting shock loading will destroy it.  Instead we need to spin it up
  %from rest in the time domain.  
n_t = (zeros(1,numnodes)) + t0; % temperature at every node
% n_isr = ones(1,numnodes); % ?
% n_isf = ones(1,numnodes); % ?
% n_a = (zeros(1,numnodes)); % node "area" is defined as 1/3 the sum of the areas of all surrounding triangles
% n_r = ones(1,numnodes); % ?
% n_e = ones(1,numnodes); % ?
% n_sh = zeros(1,numnodes); % ?

n_mf = (zeros(3,numnodes)); % Matrix for mechanical forces on each node
n_of = (zeros(3,numnodes)); % Matrix for optical forces on each node
n_af = (zeros(3,numnodes)); % Matrix for spin-up forces on each node

n_fx = (zeros(1,numnodes)); % Total force in x-direction on each node
n_fy = (zeros(1,numnodes)); % Total force in y-direction on each node
n_fz = (zeros(1,numnodes)); % Total force in z-direction on each node

n_m = (zeros(1,numnodes)); % Vector with masses of each node
%n_af = zeros(1,numnodes);
%n_ar = zeros(1,numnodes);

%n_ix = zeros(1,numnodes);
%n_iy = zeros(1,numnodes);
%n_iz = zeros(1,numnodes);

n_ta = (zeros(MAX_EDGES_PER_NODE,numnodes));   % node triangle array.  For each node (column), which triangles border that node?
                                               % Note: this is zero padded, since the number of tris per node can vary.
                                               % Use n_nt to determine how many triangles are in each column.
%Note:  Due to the way we populate this array, the triangle indices appear in increasing order (by row)                                               

% Vector specifying for each node (index) how many triangles it belongs to (row)
n_nt = (zeros(1,numnodes));

% Matrix specifying for each node (column) which edges are connected to it
% (rows)
n_ea = (zeros(MAX_EDGES_PER_NODE,numnodes));

% Vector specifying for each node (index) how many edges are connected to
% it (entry)
n_ne = (zeros(1,numnodes)); 

% Built-in MATLAB function, with first argument specifying the 2D points,
% and the second argument the edge constraints, where each row of this
% matrix defines the start and end vertex IDs of a constrained edge. See
% also https://www.mathworks.com/help/matlab/ref/delaunaytriangulation.html#d122e296596
% Each row in DT.ConnectivityList from DT represents a triangle in the
% triangulation, where in each row, the three points/nodes that make up the
% triangle are listed
DT = delaunayTriangulation([n_x' n_y'],constrainedEdges); 

edgesOK = isInterior(DT); % Test if triangles are in the interior of a 2-D constrained Delaunay triangulation, with constrained edges in triangulation defining the boundaries

TRI = DT.ConnectivityList(edgesOK, :); % TRI contains the connectivity list, i.e. each row containing three nodes that connect to a triangle
TRI2 = DT.ConnectivityList(~edgesOK, :); % These are the bullshit zero-area perimeter triangles.  
   %MK note to Ramon:  TRI2 will not be empty, at least in the case of hex meshing, along the non-orthonormal edges of the
   %mesh.  This is due to floating point arithmetic errors for irrational numbers... the nodes don't lie perfectly on a line.
   %So unless we exclude these bullshit perimeter triangles via the constrainedEdges / isInterrior calculations, we'll wind
   %up with bullshit triangles in our simulation mesh.  It is important not to include the entire DT.ConnectivityList
   %with our triangulation stored in TRI -- we must filter it with edgesOK!

t_na = TRI(:,1)'; % vector containing first row of TRI, i.e. first set of nodes belong to triangles. Its size is the total number of triangles
t_nb = TRI(:,2)'; % vector containing second row of TRI, i.e. second set of nodes belong to triangles. Its size is the total number of triangles
t_nc = TRI(:,3)'; % vector containing third row of TRI, i.e. third set of nodes belong to triangles. Its size is the total number of triangles
t_a0 = zeros(size(t_na)); % vector with areas of triangles in 3D
t_a0xy = zeros(size(t_na)); % vector with areas of triangles projected into 2D plane
t_m = zeros(size(t_na)); % vector with masses of each triangle
t_thick = zeros(size(t_na)); % vector with thicknesses of triangles
t_e1 = zeros(size(t_na)); % vector specifying for each triangle first edge
t_e2 = zeros(size(t_na)); % vector specifying for each triangle second edge
t_e3 = zeros(size(t_na)); % vector specifying for each triangle third edge
t_reg = zeros(size(t_na)); % vector defining 'special' region
t_broken = false(size(t_na));
t_notbroken = true(size(t_na));
t_cx = zeros(size(t_na));
t_cy = zeros(size(t_na));
t_cz = zeros(size(t_na));
t_t = zeros(size(t_na)) + t0;
numtris = length(t_na);
        isbig = numtris > 500;
        

% triplot(TRI,n_x,n_y); % visualize 2D Delaunay triangulation
% axis equal

% voronoi(n_x,n_y,TRI); % plot Voronoi diagram. Voronoi polygons are
% defined by boundaries such that all points in them are closer to the node
% than any other point in the set

%% End of node generation
et = toc;
disp(['Generated ' num2str(length(n_x)) ' vertices in ' num2str(et) ' sec.']);

tic;
lastET = 3; % Seconds after which progress of edge generation is displayed
disp('Generating edges...');

nodes = [];
edges = [];
triangles = [];

% meanvx = 0;
% meanvy = 0;
% meanvz = vz0;

%% Mesher  -  generate edges

% Cannot use Delaunay edges because they include the crap triangles 
% But what are "crap" triangles? For the specific case of 3 concentric
% rings, the Delaunay edges are the same as the edge array built by Mike?

% edges = DT.edges;
% e_na = edges(:,1);
% e_nb = edges(:,2);
% e_broken = zeros(1,length(e_na));

% So I need to build my own edge array, assign edge IDs, and
% cross-reference everything on my own

ne = 0; % index counting number of new edges to be created

% Note that these vectors are intentionally larger than to be expected at
% the end, with 3*length(t_na) being the total number of iterations of the
% two following for loops -> buffers

e_na = zeros(1,3*length(t_na)); % vector for first node of all edges
e_nb = zeros(1,3*length(t_na)); % vector for second node of all edges

% Buffers for each edge in every triangle, with e_l0 for original edge to
% be conidered, i.e., edge 'c' as in [AVanGelder2004], whereas e_l1 and
% e_l2 represent edges 'a' and 'b'
e_l0 = zeros(1,3*length(t_na)); % vector containing length of edges
e_l1 = zeros(1,3*length(t_na)); % vector containing length of edges
e_l2 = zeros(1,3*length(t_na)); % vector containing length of edges

e_t1 = zeros(1,3*length(t_na)); % vector specifying triangle (entry) for each edge (row)
e_t2 = zeros(1,3*length(t_na)); % vector specifying triangle (entry) for each edge (row)
e_w1 = zeros(1,3*length(t_na)); 
e_w2 = zeros(1,3*length(t_na));

% Additional buffers to calculate edge stiffnesses using Van [phi,theta]'s
% expansive form accounting for non-zero Poisson's ratio, containing the
% second term
e_w1b = zeros(1,3*length(t_na));
e_w2b = zeros(1,3*length(t_na));

e_thk1 = zeros(1,3*length(t_na)); % vector with thickness of each edge
e_thk2 = zeros(1,3*length(t_na)); % vector with thickness of each edge



permutationMat = [1 2 3; 2 3 1; 3 1 2];

% Loop over all triangles in Delaunay triangulation
for nt = 1:length(t_na) 
    
    myarea = triArea(nt); % Calculate area of specific triangle 
    myprojarea = triAreaXYProj(nt); % Calculate area of specific triangle in 2D, or projected into 2D
    t_a0(nt) = myarea; % Assign triangle area to vector component
    t_a0xy(nt) = myprojarea; % Assign projected triangle area in 2D
    
    % Define thickness of triangle... we used to have code here to vary the thickness, but it fell out of use.  
    %  if you want to vary the thickness with position, this would be a good palce to do it.
    mythick = thickness;
        
    % Calculate the mass of specific triangle
    mymass = density * myarea * mythick; % should be in grams
    
    t_m(nt) = mymass; % Assign mass of specific triangle to vector
    t_thick(nt) = mythick; % Assign thickness of specific triangle to vector
    
    % Give each node a reference to the triangle
    % Calculate number of triangles each node belongs to for n_nt
    % Assign indices of triangles to each node (column) for n_ta
    n_nt(t_na(nt)) = n_nt(t_na(nt)) + 1;
    n_ta(n_nt(t_na(nt)),t_na(nt)) = nt;
    n_nt(t_nb(nt)) = n_nt(t_nb(nt)) + 1;
    n_ta(n_nt(t_nb(nt)),t_nb(nt)) = nt;
    n_nt(t_nc(nt)) = n_nt(t_nc(nt)) + 1;
    n_ta(n_nt(t_nc(nt)),t_nc(nt)) = nt;
    
    % And distribute mass to each node by assuming that each node in a
    % triangle has 1/3 of the triangle's mass, and adding up the masses
    % from all triangles that share the same node
    n_m(t_na(nt)) = n_m(t_na(nt)) + mymass/3;
    n_m(t_nb(nt)) = n_m(t_nb(nt)) + mymass/3;
    n_m(t_nc(nt)) = n_m(t_nc(nt)) + mymass/3;
    
    % Figure out if I have a special region - NOT DONE, ask Mike about it
    tempreg = n_reg(t_na(nt));
    if (tempreg > 0) % if node of that specific triangle is on the edge 
       if (tempreg == n_reg(t_nb(nt))) % if second node of specific triangle is on the edge
           if (tempreg == n_reg(t_nc(nt))) % if third node of specific triangle is on the edge
               t_reg(nt) = tempreg; % Fill corresponding entry for that triangle with 1
           end
       end
    end
    
    % For each edge, find it or add it anew to the edge list, then give the
    % edge a reference to this triangle and add this triangle's "width" to
    % that edge for later calculation of edge stiffness.  If creating a new
    % edge, calculate its initial length, and give a reference to the edge
    % to its two nodes.
    
    startnodes = [t_na(nt) t_nb(nt) t_nc(nt)]; % Nodes of specific triangle
    endnodes =   [t_nb(nt) t_nc(nt) t_na(nt)]; % Permuted nodes of specific triangle
    
    % Loop three times over permutations of pairs from 3 nodes of triangle
    for doit3times = 1:3
        
        startnode = startnodes(doit3times); % Take one of 3 nodes from specific triangle
        endnode = endnodes(doit3times); % Take second of 3 nodes from specific triangle that is different from 'startnode' by definition of vector endnodes
        
        startnode2 = startnodes(permutationMat(doit3times,2));
        endnode2 = endnodes(permutationMat(doit3times,2));
        
        startnode3 = startnodes(permutationMat(doit3times,3));
        endnode3 = endnodes(permutationMat(doit3times,3));
        
        % Self-written function to find edge connecting two nodes. Note
        % that since e_na and e_nb are vectors containing 0 only at the
        % beginning, for the first iterations, findedge will only return
        % existingEdgeIdx = 0
        existingEdgeIdx = findedge(startnode,endnode);
        
        if (existingEdgeIdx == 0) % new edge, create and add to list, back-reference from nodes.
            
            ne = ne + 1;
            
            e_na(ne) = startnode;
            e_nb(ne) = endnode;
            
            e_l0(ne) = distanceNodes(startnode,endnode);
            e_l1(ne) = distanceNodes(startnode2,endnode2);
            e_l2(ne) = distanceNodes(startnode3,endnode3);
            
            % e_t1, e_w1, e_thk1 with '1' (one, NOT small L) are variables
            % reserved for newly created edges
            e_t1(ne) = nt; % reference this as first triangle
            
%             e_w1(ne) = (2/3) * myarea / e_l0(ne); % Mike, why 2/3?  
            e_w1(ne) = myarea / e_l0(ne); % according to [AVanGelder2004]

            % [AVanGelder2004]: expansive form for stiffness accounting for
            % non-zero Poisson's ratio
            e_w1b(ne) = abs( (e_l0(ne)).^2 + abs(e_l1(ne)).^2 - ...
                abs(e_l2(ne)).^2 ) ./ (8 .* myarea); 
            
            e_thk1(ne) = t_thick(nt); % assign thickness to edge
            
            % Fill vector n_ne at column 'startnode' for that specific node
            % with numbers of edges that are connected to it iteratively by
            % continuously adding 1
            n_ne(startnode) = n_ne(startnode) + 1;  % back-reference new edge from node a
            
            % Fill matrix n_ea at column 'startnode' and row
            % 'n_ne(startnode)' with edge 'ne' being connected to this node
            n_ea( n_ne(startnode), startnode ) = ne; 
            
            % Fill vector n_ne at column 'endnode' for that specific node
            % with numbers of edges that are connected to it iteratively by
            % continuously adding 1
            n_ne(endnode) = n_ne(endnode) + 1;  % back-reference new edge from node b
            
            % Fill matrix n_ea at column 'endnode' and row
            % 'n_ne(endnode)' with edge 'ne' being connected to this node
            n_ea( n_ne(endnode), endnode ) = ne;
            
            % Assign to variable the current edge index 'ne'
            existingEdgeIdx = ne;
            
        else
            
            if (e_t2(existingEdgeIdx) > 0)
                warning(['Found more than 2 triangles at edge ' num2str(existingEdgeIdx) ]);
            else
                % e_t2, e_w2, e_thk2 with '2' are variables reserved for
                % already created edges
                e_t2(existingEdgeIdx) = nt;
%                 e_w2(existingEdgeIdx) = (2/3) * myarea / e_l0(existingEdgeIdx);
                e_w2(existingEdgeIdx) = myarea / e_l0(existingEdgeIdx);
                e_w2b(existingEdgeIdx) = abs( (e_l0(existingEdgeIdx)).^2 + ...
                    abs(e_l1(existingEdgeIdx)).^2 - ...
                    abs(e_l2(existingEdgeIdx)).^2 ) ./ (8 .* myarea);  

                e_thk2(existingEdgeIdx) = t_thick(nt);
            end
            
        end
        
        % Make case where to assign edge to, t_e1, t_e2 or t_e3 depending
        % on which permutation pair is currently considered
        if (doit3times == 1)
            t_e1(nt) = existingEdgeIdx;
        elseif (doit3times == 2)
            t_e2(nt) = existingEdgeIdx;
        elseif (doit3times == 3)
            t_e3(nt) = existingEdgeIdx;
        end
        
    end
    
    et = toc;
    if (et - lastET) > 1
        % Display progress of edge generation to user in %
        lastET = et;
        disp([ '    ' num2str(100*nt/length(t_na)) '%...']);
    end
end

% Assign vector with 3D triangle areas to new variable - why? This to keep
% the original triangle areas stored in t_a0, which will be used to
% evaluate the total area of a semi-broken sail
t_a = t_a0;

% Initialize vector to be of same size as vector containing nodes
n_a = zeros(size(n_x));

% Loop over number of triangles to assign to n_a three times the area of
% the triangle to which the node belongs. Note that n_a has as many
% elements as there are nodes, with each element stating 3x the triangle
% area to which this node belongs
for ntr = 1:length(t_na)
        n_a(   t_na(ntr) )  = n_a( t_na(ntr) ) + t_a(ntr);
        n_a(   t_nb(ntr) )  = n_a( t_nb(ntr) ) + t_a(ntr);
        n_a(   t_nc(ntr) )  = n_a( t_nc(ntr) ) + t_a(ntr);
end

% Divide calculated areas by 3 to account for fact that triangles have been
% counted thrice, since there are 3 nodes per triangle
n_a = n_a./3;

%  Re-purposing n_a0 from initial angular position to mean the initial
%  "area" of the node, which is used in calculating radiative cooling.
n_ang0 = n_a;

%% Now trim the edge arrays to minimum size
% This is because in the initial initialization of these vectors above,
% their sizes were assumed to be 3*length(t_na), i.e. since each triangle
% has three edges, the maximum number of edges will be 3 times the number
% of triangles. However, most edges are shared between two triangles, and
% edges on the boundary of the mesh even only by one. Therefore, the actual
% sizes of these vectors is smaller than 3*length(t_na)

% Alternative way to calculate actual number of edges is multiply number of
% triangles times 3, then divide by 2 to account for double counting of
% triangles, but then add number of edges/triangles on the boundary divided
% by 2 

e_na = e_na(1:ne);
e_nb = e_nb(1:ne);
e_l0 = e_l0(1:ne);
e_t1 = e_t1(1:ne);
e_t2 = e_t2(1:ne);
e_w1 = e_w1(1:ne);
e_w2 = e_w2(1:ne);

e_w1b = e_w1b(1:ne);
e_w2b = e_w2b(1:ne);

e_thk1 = e_thk1(1:ne);
e_thk2 = e_thk2(1:ne);
e_al = zeros(size(e_w2)); % 'aspect length', representing effective thickness/cross-section of each edge
e_broken = false(size(e_w2));
e_notbroken = true(size(e_w2));
e_reg = e_al;
e_m = e_al;  %average mass of the two nodes at either end of the edge... not used for practical purposes, but having this
  % lets us crudely approximate the recoil effect of edge failure.  Does not conserve mass.  Ignore this, it is only for
  % fun.
edgeregionhist = [0 0 0 0 0]; % vector with 1st element = # boundary edges

%% Generate edge-to-node matrix
% Need this for assigning mechanical forces and heat flow to nodes
% Matrix is constructed such that each column denotes a specific node, and
% its entries are nonzero if the corresponding row indexing a specific node
% contains this node.  This lets us transfer edge calculations (mech. forces, heat
% flow, ...) to nodes via matrix multiplication, rather than iteration, which 
% makes the code a billion times faster.

M_e2n = (zeros(ne,numnodes));

% Needed for program to work efficiently

for ne = 1:length(e_na)
    % Assign direction A -> B with nodes A & B and edge -> to value 1, and
    % the opposite direction B -> A with the same nodes and edge to -1
    M_e2n(ne,e_na(ne)) = 1; % test by switching this line below, should "exponentially blow up"
    M_e2n(ne,e_nb(ne)) = -1;
%     disp([ne,e_na(ne)])
%     disp([ne,e_nb(ne)])
    e_m(ne) = 0.5 * n_m(e_na(ne)) + 0.5 * n_m(e_nb(ne)) ;
end
M_e2n = sparse(M_e2n); % squeeze out 0 elements to save memory.  This is critical for speed.
    
%% Generate triangle-to-node matrix
% Need this for distributing optical forces and absorbed heat to nodes
% Matrix is constructed such that each column denotes a specific node, and
% its entries are 1 if the corresponding row indexing a specific triangle
% contains this node, and otherwise 0.  This lets us assign triangle calcs
% (e.g., optical force) to nodes via matrix multiplication, rather than 
% a for loop, which makes the code a billion times faster.

M_t2n = (zeros(length(t_na),numnodes));
for nt = 1:length(t_na)
    M_t2n(nt,t_na(nt)) = 1;
    M_t2n(nt,t_nb(nt)) = 1;
    M_t2n(nt,t_nc(nt)) = 1;
end
M_t2n = sparse(M_t2n); % squeeze out 0 elements to save memory

M_nxyz2tev1 = zeros(numnodes,numtris);
M_nxyz2tev2 = zeros(numnodes,numtris);

for nt = 1:length(t_na)
    M_nxyz2tev1( t_na(nt), nt) = -1;
    M_nxyz2tev1( t_nb(nt), nt) = 1;
    M_nxyz2tev2( t_na(nt), nt) = -1;
    M_nxyz2tev2( t_nc(nt), nt) = 1;
end
M_nxyz2tev1 = sparse(M_nxyz2tev1);
M_nxyz2tev2 = sparse(M_nxyz2tev2);

%% Calculate the spring constants and maybe also thermal conductivity 

% Auxiliary counting variable, this will count number of edges on the
% boundary of the mesh
numouter = 0;

% Loop over all edges
for ne = 1:length(e_na)
    if (e_t2(ne) == 0)
        numouter = numouter + 1;
%         e_w2(ne)=e_w1(ne);
%         e_thk2(ne)=e_thk1(ne);
    end
    
    % Assign to variable 'tempreg' whether node belong to edge under
    % consideration is on the mesh boundary (1) or in the interior (0)
    tempreg = n_reg(e_na(ne));
    
    % If we deal with a boundary node, and the other node that together
    % with it forms an edge is on the boundary of the mesh too, then we
    % have an boundary edge, and we note this in 'e_reg', which indexes
    % each edge in the mesh, and notes whether it is on the mesh boundary
    % (1) or not (0) -> the total number of ones in 'e_reg' should be equal
    % to the number of edges on the mesh boundary. Not sure about
    % 'edgeregionhist', whose first element seems to display the total
    % number of edges on the mesh boundary (but why is it a 1x5 vector?)
    if (tempreg > 0) && (tempreg == n_reg(e_nb(ne)) )
        e_reg(ne) = tempreg;
        edgeregionhist(tempreg) = edgeregionhist(tempreg) + 1;
    end
end

disp(['   found... ' num2str(numouter) ' outer edges.']);

% According to [AVanGelder2004], the spring stiffness coefficient of an
% edge for triangulated spring meshes varies as triangle-area over
% edge-length squared, i.e. k = E2D * sum_{e} area(T_e) / |c|^2, where the
% sum is over (nominally two) triangles T_e incident upon edge c, |c| is
% the length of the edge and E2D the two-dimensional Young's modulus of the
% membrane to be simulated. If the membrane has a constant thickness of t,
% then E2D = E * t, with E being the Young's modulus. Note that equation 27
% depicts a more general expression for the spring constant, but needs to
% be treated with care for positive Poisson ratios, which could result in
% negative spring constants. Therefore, [AVanGelder2004] assumes nu = 0
% (zero Poisson ratio), where the stiffness varies as area over length
% squared, with proportionality factor E2D. It remains to be seen whether
% this is a valid assumption for our purposes!


% According to [AVanGelder2004] expansive form for stiffness, accounting
% for non-zero Poisson coefficient, see equation (29)

e_al = (e_thk1.*e_w1 + e_thk2.*e_w2) ./ e_l0;
% e_al = (e_thk1.*e_w1 + e_thk2.*e_w2) ./ e_l0 ./ (1 + nu) + ...
%     (nu ./ ( 1 - nu.^2)) .* (e_thk1.*e_w1b + e_thk2.*e_w2b);

if sum(e_al < 0) > 0
    warning('Negative edge stiffness found, consider setting Poisson ratio to 0');
end

e_k = e_al .* (Ymodmm);

% To be written ...
% al ... area & length, describes relative aspect ratio, length to
% thickness and cross-sectional area, approximate each edge as rectangular
% material slab to first order
% e_tc = e_al .* Thermcondmm;

% Frequency of each edge is calculated as the square root of the edge's
% spring constant divided by the mass as the sum of the two node masses
% that constitute the edge 
e_f0 = sqrt( e_k ./ ((n_m(e_na) + n_m(e_nb)) .*1e-6) ) /(2*pi);

% Frequency of each triangle is calculated as the averaged sum of the
% individual edge frequencies that make up the triangle
t_f0 = (e_f0(t_na) + e_f0(t_nb) + e_f0(t_nc) )./3;

% Calculate the (initial) edge temperature as the averaged temperature of 
% the sum of the temperatures of the two nodes that form the edge
e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;

et = toc;
disp(['Generated ' num2str(length(e_na)) ' edges in ' num2str(et) ' sec.']);
numedges = length(e_na);

% Start new timer
tic;

%% Calculate the triangle centroids
% The triangle centroid is the point of intersection of its medians, i.e.
% the lines joining each vertex with the midpoint of the opposite site
% [German: geometrischer Schwerpunkt]
% For a triangle, the x, y or z coordinate of the centroid is calculated as
% the arithmetic mean of the respective x, y or z vertex/nodes coordinates

t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;

%t_cx,y,z will get updated with each simulation iteration.  

%Here, we'll store the initial triangle centroids in a separate array, in
%case we are ever curious about the triangle's original position at a later
%time.

t_cx0 = t_cx;
t_cy0 = t_cy;
t_cz0 = t_cz;

%% centerTris is a quick list of the triangles surrounding the center node.
% I average the norms of the center triangles to approximate the pitch of the sail during the simulation
n_centerTris = n_ta(:,1);
n_centerTris = n_centerTris(n_centerTris>0);



              
%% Region mapping
% assign values for texture and texture orientation for each triangle
t_texa = zeros(1,length(t_na));
t_tex = t_texa;
t_costexa = ones(1,length(t_na));
t_sintexa = t_texa;

% Create two matrices containing the vectors that span up each triangle in
% the mesh (note that two vectors are sufficient for every triangle)
temp_tri_edg1 = [   n_x(t_nb) - n_x(t_na) ;
                    n_y(t_nb) - n_y(t_na) ;
                    n_z(t_nb) - n_z(t_na) ];
temp_tri_edg2 = [   n_x(t_nc) - n_x(t_na) ;
                    n_y(t_nc) - n_y(t_na) ;
                    n_z(t_nc) - n_z(t_na) ];

% Normalize each diagonal spanning up each triangle in the mesh
temp_tri_edg1_norm = temp_tri_edg1./vecnorm(temp_tri_edg1);
temp_tri_edg2_norm = temp_tri_edg2./vecnorm(temp_tri_edg2);

% calculate and normalize the triangle normal vectors
temp_tri_normvec = cross(temp_tri_edg1_norm,temp_tri_edg2_norm);
temp_tri_normvec = temp_tri_normvec ./ vecnorm( temp_tri_normvec);

% Initialize matrix for normalized vectors projected into corresponding
% triangles using characteristic directional vector
temp_tex_proj_norms = zeros(3,length(t_na));

% Auxiliary counting variable, counts number of textured triangles (i.e.
% excluding those in the interior region)
numTexTris = 0;

% Calculate angle of centroids of all triangles
t_centroid_ang = atan2(t_cy0,t_cx0);
t_centroid_ang(t_centroid_ang < 0) = t_centroid_ang(t_centroid_ang < 0) + 2*pi; 

t_centroid_r = sqrt( t_cx0.^2 + t_cy0.^2 );

for ntm = 1:length(texture_map_modes)
    texture_map_mode = texture_map_modes{ntm};
    % Loop over all triangles in mesh
    for nt = 1:length(t_na)
        %Here, we assign a value of 'texture' (i.e., a region) to each triangle.
        
        this_existingtex = t_tex(nt);    
        dotex = 0;
        
        if isequal(texture_map_mode, 'pizza')
            radi = sum( t_centroid_ang(nt) >= pizza_angles );
            if pizza_texs( radi ) >= 0
                t_tex(nt) = pizza_texs( radi );
                this_tex_dir = pizza_dirs(:,radi);
                dotex = 1;
            end
            
        elseif isequal(texture_map_mode, 'stuffedcrust')
            radi = sum( t_centroid_ang(nt) >= pizza_angles );
            if (t_centroid_r(nt) >  pizza_crust_radii(radi) )
                if pizza_crust_texs(radi) >= 0
                    t_tex(nt) = pizza_crust_texs(radi);
                    this_tex_dir = pizza_crust_dirs(:,radi);
                    dotex = 1;
                end
            else
                if pizza_texs( radi) >= 0
                    t_tex(nt) = pizza_texs( radi);
                    this_tex_dir = pizza_dirs(:,radi);
                    dotex = 1;
                end
            end
            
        elseif isequal(texture_map_mode, 'brownies')
            bxi = sum( t_cx(nt) > brownies_slices_x )+1;
            byi = sum( t_cy(nt) > brownies_slices_y )+1;
            mydirs = brownies_dirs{byi};
            if brownies_texs(bxi,byi) >= 0
                t_tex(nt) = brownies_texs(bxi,byi);
                this_tex_dir = mydirs(:,bxi);
                dotex = 1;
            end
            
        elseif isequal(texture_map_mode,'user')
            error('User texture mapping is not yet specified.');
        else
            error('Invalid texture mapping mode in mesh generator');
        end
        
        if (~this_existingtex && t_tex(nt)) numTexTris = numTexTris+1; end
        if (this_existingtex && ~t_tex(nt)) numTexTris = numTexTris-1; end
        
        if dotex
            if isequal(this_tex_dir, vec_radial)
                this_tex_dir = [t_cx(nt); t_cy(nt); 0]./t_centroid_r(nt);
            elseif isequal(this_tex_dir, -vec_radial)
                this_tex_dir = -[t_cx(nt); t_cy(nt); 0]./t_centroid_r(nt);
            elseif isequal(this_tex_dir, vec_tang)
                this_tex_dir = [-t_cy(nt); t_cx(nt); 0]./t_centroid_r(nt);
            elseif isequal(this_tex_dir, -vec_tang)
                this_tex_dir = [t_cy(nt); -t_cx(nt); 0]./t_centroid_r(nt);
            end

            temp_tex_proj = this_tex_dir - dot(temp_tri_normvec(:,nt),this_tex_dir) * temp_tri_normvec(:,nt) ;
            temp_tex_proj_norm = temp_tex_proj ./ norm(temp_tex_proj);
            temp_tex_proj_norms(:,nt) = temp_tex_proj_norm;
            %... get
            temp_tex_adot = dot(  temp_tri_edg1_norm(:,nt), temp_tex_proj_norm );
            t_debug(nt) = temp_tex_adot;
            temp_tex_adot = max(min( temp_tex_adot, 1), -1);
            %the min(...,1) here is necessary because sometimes floating point
            %division error causes the dot product to be slightly more than 1.
            % (like ... 1 + 2e-16) which makes the angle array complex :(
            t_texa(nt) = acos( temp_tex_adot ) ;

            if dot( cross(temp_tex_proj_norm, temp_tri_edg1_norm(:,nt)), temp_tri_normvec(:,nt) ) > 0
                t_texa(nt) = t_texa(nt) *-1;
            end
            t_costexa(nt) = cos(t_texa(nt));
            t_sintexa(nt) = sin(t_texa(nt));
        end
    end
end

% Assign temporary matrix with normalized projected vectors specific to
% each texcture region to another matrix variable
t_texn = temp_tex_proj_norms;

% Calculate total mass of nodes by summing up masses from every node
totalmass_nodes = sum(n_m);

% Calculate total mass of triangles by summing up masses from every
% triangle in mesh. As a sanity check, this should be identical to the
% total mass of nodes, since to calculate the mass of each node, the mass
% of a triangle was equally distributed to its nodes
totalmass_tris = sum(t_m);

% Define total mass of sail as total mass of nodes
totalmass = totalmass_nodes;

% Note that the total mass should roughly be equal to the theoretical mass
% calculated as density * thickness * area, where for a round sail, area is
% equal to pi * (x0 * radius)^2. The relative error becomes smaller the
% finer the mesh is
totalmass_theory_round = density * thickness * pi * (x0 * radialRings)^2;

% Calculate total area of sail. This should be roughly equal to the
% theoretical sail area as pi * (x0 * radius)^2. The relative error becomes
% smaller the finer the mesh is
totalarea = sum(t_a0);   
totalarea_theory_round = pi * (x0 * radialRings)^2;

% Calculate total mean x, y and z velocity components from averaging over
% all velocities of every node in x, y and z. If no spin-stabilization is
% considered, and there is no initial velocity, then these mean velocity
% components will be 0
meanvx = mean(n_vx);
meanvy = mean(n_vy);
meanvz = mean(n_vz);

% Calculate maximum velocity determined by maximum velocity among all nodes
maxv = max( sqrt( n_vx.^2 + n_vy.^2 + n_vz.^2 ) );

% Calculate center-of-mass of sail as mass-weighted sum of every x/y/z
% coordinate divided by total mass of sail
coms_x = n_x .* n_m;
coms_y = n_y .* n_m;
coms_z = n_z .* n_m;
comx = sum(coms_x)./totalmass_nodes ;
comy = sum(coms_y)./totalmass_nodes ;
comz = sum(coms_z)./totalmass_nodes ;

angVel = [0; 0; 0];
angVelRigid = [0; 0; 0];
angMomRigid = [0; 0; 0];

%if rigid_mode
    
    n_dx = n_x - comx;
    n_dy = n_y - comy;
    n_dz = n_z - comz;
    n_dxyz = [ n_dx; n_dy; n_dz ];
    
    n_dx2 = n_dx.^2;
    n_dy2 = n_dy.^2;
    n_dz2 = n_dz.^2;
    n_dr = sqrt( n_dx2 + n_dy2 + n_dz2);
    
    inertiaT = [  sum( n_m .* ( n_dy2 + n_dz2 ) )    ...
        sum( -n_m .*  n_dx .* n_dy  )         ...
        sum( -n_m .* n_dx .* n_dz );  ...
        sum( -n_m .* n_dy .* n_dx )            ...
        sum( n_m .* (n_dx2 + n_dz2) )  ...
        sum( -n_m .* n_dy .* n_dz );   ...
        sum( -n_m .* n_dz .* n_dx )            ...
        sum( -n_m .* n_dz .* n_dy )        ...
        sum( n_m .* (n_dx2 + n_dy2) ) ] ;
    
    approxSpinAng = atan2(n_dy(xradidx), n_dx(xradidx));
%end

% Display information about generated mesh and associated properties to
% user in command window/console
disp(['Num nodes: ' num2str(length(n_x)) ]);
disp(['Num edges: ' num2str(length(e_na)) ]);
disp(['Num triangles: ' num2str(length(t_na)) ]);
disp(['Mass: ' num2str(totalmass_nodes) ' / ' num2str(totalmass_nodes) ' grams']); 
disp(['Area: ' num2str(totalarea) ' mm^2   (' num2str(sum(t_a0xy)) ' mm^2 aperture)' ]);
disp(['AspectRatio: ' num2str((max(n_z)-min(n_z))/(((max(n_x)-min(n_x))+(max(n_y)-min(n_y)))/2)) ] );
disp(   ['Spatial extent:  X (mm)      Y (mm)      Z(mm) ']);
fprintf('                 %-8g    %-8g    %-8g\n',max(n_x)-min(n_x),max(n_y)-min(n_y),max(n_z)-min(n_z))
disp(['Avg speed of propagation: ' num2str(mean(e_f0.*e_l0/1000)) ' m/s']);

% Time step to be used for simulation is 10x smaller than that recommended
% maximum time step given the highest resonant frequency
% In [XProvot1995], it is stated that numerically solving linear equations
% as these may become ill-conditioned if the time step is greater than the
% natural period of the system. This is why the highest resonant frequency
% of the membrane depicts an upper bound for the time step.
% rec_dt = (1/max(e_f0)/10); 
rec_dt = (1/max(e_f0)/(10*2*pi)); 

disp(['Highest resonant frequency: ' num2str(max(e_f0))  ' Hz    (recommended max time step ' num2str(rec_dt) ' sec)']);

% Used in simulate.m, but what for? It corresponds to the most negative y
% value of any node in the mesh
origmeshymin = min(n_y);

% trimesh(TRI,n_x,n_y,n_z, n_a0);
% hold on;

n_x0 = n_x;
n_y0 = n_y;
n_z0 = n_z;

%% this seems dumb, but...
thisseemsdumb = 1:numtris;      %this precalculated thing helps with raytracing
tri_broken2orig_map = thisseemsdumb;

%% Translate and rotate lightsail
% no longer used -- use the t0t settings in the main simulator instead!
% if 0
% psi_ex = deg2rad(130);
% theta_ex = deg2rad(2.1);
% phi_ex = deg2rad(3.5);
% 
% psi_ex = deg2rad(0);
% theta_ex = deg2rad(0);
% phi_ex = deg2rad(0);
% 
% DRM123_11 = @(psi,theta,phi) cos(psi).*cos(theta);
% DRM123_12 = @(psi,theta,phi) cos(psi)*sin(theta).*sin(phi) + sin(psi).*cos(phi);
% DRM123_13 = @(psi,theta,phi) -cos(psi).*sin(theta).*cos(phi) + sin(psi).*sin(phi);
% DRM123_21 = @(psi,theta,phi) -sin(psi).*cos(theta);
% DRM123_22 = @(psi,theta,phi) -sin(psi).*sin(theta).*sin(phi) + cos(psi).*cos(phi);
% DRM123_23 = @(psi,theta,phi) sin(psi)*sin(theta).*cos(phi) + cos(psi).*sin(phi);
% DRM123_31 = @(psi,theta,phi) sin(theta);
% DRM123_32 = @(psi,theta,phi) -cos(theta).*sin(phi);
% DRM123_33 = @(psi,theta,phi) cos(theta).*cos(phi);
% 
% DRM_mat = [cos(psi_ex).*cos(theta_ex), cos(psi_ex)*sin(theta_ex).*sin(phi_ex) + sin(psi_ex).*cos(phi_ex), -cos(psi_ex).*sin(theta_ex).*cos(phi_ex) + sin(psi_ex).*sin(phi_ex); ...
%     -sin(psi_ex).*cos(theta_ex), -sin(psi_ex).*sin(theta_ex).*sin(phi_ex) + cos(psi_ex).*cos(phi_ex), sin(psi_ex)*sin(theta_ex).*cos(phi_ex) + cos(psi_ex).*sin(phi_ex); ...
%     sin(theta_ex), -cos(theta_ex).*sin(phi_ex), cos(theta_ex).*cos(phi_ex)];
% 
% n_x_tmp = n_x .* DRM123_11(psi_ex,theta_ex,phi_ex) + ...
%     n_y .* DRM123_12(psi_ex,theta_ex,phi_ex) + ...
%     n_z .* DRM123_13(psi_ex,theta_ex,phi_ex);
% 
% n_y_tmp = n_x .* DRM123_21(psi_ex,theta_ex,phi_ex) + ...
%     n_y .* DRM123_22(psi_ex,theta_ex,phi_ex) + ...
%     n_z .* DRM123_23(psi_ex,theta_ex,phi_ex);
% 
% n_z_tmp = n_x .* DRM123_31(psi_ex,theta_ex,phi_ex) + ...
%     n_y .* DRM123_32(psi_ex,theta_ex,phi_ex) + ...
%     n_z .* DRM123_33(psi_ex,theta_ex,phi_ex);
% 
% n_x_BF = n_x;
% n_y_BF = n_y;
% n_z_BF = n_z;
% 
% n_x = n_x_tmp;
% n_y = n_y_tmp;
% n_z = n_z_tmp;
% 
%                 comx = sum(coms_x)./totalmass_nodes ;
%                 comy = sum(coms_y)./totalmass_nodes ;
%                 comz = sum(coms_z)./totalmass_nodes ;
%                 
%                 n_dx = n_x - comx;
%                 n_dy = n_y - comy;
%                 n_dz = n_z - comz;
%                 n_dxyz = [n_dx; n_dy; n_dz];
%                 
%                 n_dx2 = n_dx.^2;
%                 n_dy2 = n_dy.^2;
%                 n_dz2 = n_dz.^2;
%                 n_dr = sqrt( n_dx2 + n_dy2 + n_dz2);
%                 
%                 approxSpinAng = atan2(n_dy(xradidx), n_dx(xradidx));
%                 
%                 inertiaT = [  sum( n_m .* ( n_dy2 + n_dz2 ) )        ...
%                         sum( -n_m .* n_dx .* n_dy)         ...
%                         sum( -n_m .* n_dx .* n_dz );  ...
%                         sum( -n_m .* n_dy .* n_dx )            ...
%                         sum( n_m .* (n_dx2 + n_dz2) )      ...
%                         sum( -n_m .* n_dy .* n_dz );  ...
%                         sum( -n_m .* n_dz .* n_dx )            ...
%                         sum( -n_m .* n_dz .* n_dy )        ...
%                         sum( n_m .* (n_dx2 + n_dy2) ) ]';
%                     
%                 t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
%                 t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
%                 t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;
%                 
%                 t_ev1 = [n_x; n_y; n_z]* M_nxyz2tev1;
%                 t_ev2 = [n_x; n_y; n_z]* M_nxyz2tev2;
%                 t_norms = cross( t_ev1, t_ev2 );
%                 
%                 t_a = max(vecnorm(t_norms), 1e-8);
%             
%             % Normalize normal vectors to each triangle
%             t_norms = t_norms ./ t_a ; %ok now t_norms contains proper normal vectors (normalized)
%             t_a = 0.5 .* t_a;
%             t_ev1n = t_ev1./(max(vecnorm(t_ev1),1e-8));
%             
%             t_texn = t_ev1n .* (t_costexa) + ...
%             cross(t_norms, t_ev1n ) .* (t_sintexa)  + ...
%             t_norms .* dot(t_norms, t_ev1n ) .* (1 - cos(t_texa)) ;
% end
            
% Plot triangulated mesh using Delaunay triangulation with different
% texture regions
hMeshPreview = figure('Position', [50   600    1200   600]);
trisurf(TRI,n_x,n_y,n_z, t_tex);
axis equal
hc = colorbar ;
hc.Label.String = 'Texture region';
hold on;
set(gca,'fontsize',16)
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
ctrlptidxs = [xradidx yradidx nxradidx nyradidx];


plot3(n_x(meshcs), n_y(meshcs), n_z(meshcs), '-mo', 'MarkerSize',10 )
plot3(n_x(ctrlptidxs), n_y(ctrlptidxs), n_z(ctrlptidxs), 'co', 'MarkerSize', 10)
plot3(n_x(xradidx), n_y(xradidx), n_z(xradidx), 'b*', 'MarkerSize', 10)
plot3(n_x(end), n_y(end), n_z(end), 'kx', 'MarkerSize', 10)
if exist('filebasename')
    title(['Simulation mesh: ' filebasename], 'Interpreter', 'none');
end
%colormap(gca,'jet');


% Matrix 'temp_tex_ax_norm' contains first spanning vectors of each
% triangle in every texture region rotated by angle calculated from
% projected vector (see above) around respective normal on triangles 
% for nt = 1:length(t_na)
%     temp_tex_ax_norm(:,nt) = rodrigues_rot(  temp_tri_edg1_norm(:,nt), temp_tri_normvec(:,nt),   t_texa(nt)  );
% end

% This is just another way to calculate rotated texture axes
temp_tex_ax_norm = temp_tri_edg1_norm .* cos(t_texa) + ...
                   cross(temp_tri_normvec, temp_tri_edg1_norm ) .* sin(t_texa)  + ...
                   temp_tri_normvec .* dot(temp_tri_normvec, temp_tri_edg1_norm ) .* (1 - cos(t_texa)) ;


quiver3( t_cx, t_cy, t_cz, t_texn(1,:), t_texn(2,:), t_texn(3,:),'r' ); 
%quiver3( t_cx, t_cy, t_cz, temp_tex_ax_norm(1,:), temp_tex_ax_norm(2,:), temp_tex_ax_norm(3,:),'r' ); 
% quiver3( t_cx, t_cy, t_cz, temp_tex_ax_norm(1,:), temp_tex_ax_norm(2,:), temp_tex_ax_norm(3,:),'r','LineWidth',2); 
%quiver3( t_cx, t_cy, t_cz, temp_tex_proj_norms(1,:), temp_tex_proj_norms(2,:), temp_tex_proj_norms(3,:), 'r'); 
%quiver3( t_cx, t_cy, t_cz, temp_tex_ax_norm(1,:), temp_tex_ax_norm(2,:), temp_tex_ax_norm(3,:), 'r'); 
hold off
hold off;