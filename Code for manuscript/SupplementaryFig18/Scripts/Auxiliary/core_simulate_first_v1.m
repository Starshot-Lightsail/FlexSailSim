%% Core Computation Script for STARSHOT LIGHTSAIL SIMULATOR
% (c) 2019-2023 Michael D. Kelzenberg, Ramon Gao, California Institute of Technology

% Do time-dependent optical, mechanical and thermal calculations here
% Outsourced from main code for all versions 19 and earlier to allow for
% use of more advanced numerical integration methods

% Version 1, May 14th, 2023

%% Step 0: Calculate spin-up forces
if  (tt < 0) && rps_target  % spin-up forces are applied for tt<0.  
    n_arad = sqrt( n_x.^2 + n_y.^2 );
    n_af(1,:) = sign(rps_target)*spinup_gs/88.2.*n_m.*(n_y) ;%%% .* (n_arad < 0.6 * radiusmm);
                n_af(2,:) = sign(rps_target)*-spinup_gs/88.2*n_m.*(n_x) ;%%%.* (n_arad < 0.6 * radiusmm);
end

if ( tt >= 0) || anybroken  % stop spin-up at t=0 (or if any edges broke during spinup)
    n_af(:) = 0;
end


%% Step 1: Calculate beam intensity profile and direction
% Here, we set t_IMag and t_IDir, which are the magnitude and direction of the incident light at each
% triangle centroid. Normally, t_IDir is (0, 0, 1), so we never change it here, but in theory, you could set
% whatever value you want. Note: t_IDir must be unit vectors.

pwrRampVal = pwrRampFcn(tt);

if (pwrRampVal == 0)
    t_IMag = zeros(size(t_cx));
else
    t_optFcn_r2 = (t_cx - 0).^2 + (t_cy - 0).^2;
    t_IMag = pwrRampVal * I0 * exp( - t_optFcn_r2 ./ (Irad^2) );  % This is in W/mm2.  Divide by c0 to get newtons / mm2 radiation pressure (N/mm2)
  
    if nto_startAccel == 0
        nto_startAccel = nto;
        nt_startAccel = nt;
    end
end

% Keep track of time spent on spin-up force and light intensity calculations
timers(1) = timers(1) + toc;
tic;


%% Step 2: Calculate triangle area, normal vectors, and texture orientation vectors
if isbig  % choose between two calculation mechanisms for optimal performance on my home computer (MK)
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

if read_LUT

    t_norms = cross( t_ev1, t_ev2 );  % triangle normal vectors, not yet normalized; we will normalize shortly

    % We determine if a triangle is upside down based on the z component of
    % its normal vector.  Here, there's an implicit assumption that light is incident along the 
    % positive Z axis, so if you ever set t_IDir to something other than (0, 0, 1), you might need to
    % determine a new way to detect upside down triangles.
    upsidedowntris = sum( t_norms(3,:) < 0);   

elseif ~read_LUT

    t_norms = cross( t_ev2, t_ev1 );  % for propulsion in -z - sorry, I (Ramon) still struggle to figure out correct calculation of forces in +z
    upsidedowntris(nt) = sum( t_norms(3,:) > 0); 

end

if upsidedowntris
    if ~terminate_t && (terminate_on_upsidedown)  % initiate delayed termination if triangles are upside down
        terminate_t = tt + termination_delay_s;
        disp('Initiating delayed termination due to upside down triangles!');
    end
end

% Calculate length of the non-normalized normal vectors, which happen to also be 2x the triangle area
t_a = max(vecnorm(t_norms), 1e-8);  % initially, this is 2x the triangle area

% Normalize normal vectors to each triangle
t_norms = t_norms ./ t_a ; % ok, now t_norms contains proper unit normal vectors

t_a = 0.5 .* t_a; % ok, now t_a contains the actual area of each triangle

% If any triangles are 'broken', we replace their now-meanlingless calculated area with their original area
if anybroken
    t_a(t_broken) = t_a0(t_broken);  % Don't worry, other calculations will continue ignoring broken triangles.  
end

% Normalize the #1 edge vectors
t_ev1n = t_ev1./ vecnorm(t_ev1);%test;

% During mesh generation, we recorded a value 't_texa' for each triangle, which is the in-plane angle
% between the texture orientation vector and the #1 edge vector.  We use that now to calculate the 
% new texture orientation vectors for the current mesh state... well, actually, we use the cached values 
% of the sin and cosine of this angle to hopefully save time (t_costexa, t_sintexa).
t_texn = t_ev1n .* (t_costexa) + ...
    cross(t_norms, t_ev1n ) .* (t_sintexa)  + ...
    t_norms .* dot(t_norms, t_ev1n ) .* (1 - cos(t_texa)) ;  %TODO: Can I safely replace the cos() with t_costexa???

% Calculate cos(theta), where theta corresponds to the light incidence angle at each triangle.              
t_costheta = dot(t_norms, t_IDir);

% Here we calculate the projected area of non-broken triangles in the X-Y plane.  This is essentially the
% aperture area of our sail in the plane of beam incidence.  Note again the implicit assumption that
% illumination is along (0,0,1).  
t_aproj = t_a .* t_costheta .* (t_notbroken);


%% Step 3: Calculate reflection and absorption for specular regions
%  Currently we use simple scalar values for reflection and absorption.  If you were so inclined, you could
%  implement different physics to calculate how much light should be reflected and absorbed by specular sail
%  surfaces, such as using Fresnel equations.  However, to implement more complex optical behavior, I 
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
    if useSpecularTMM
        t_theta = acos(t_costheta);
        if min(t_theta) < 0 
            warning('Negative theta in useSpecularTMM code');
            t_theta(t_theta < 0) = 0;  % I don't think this should ever happen
        end
        if max(t_theta) > pi/2
        %    warning('Reverse illuminated triangles in specularTMM code');
            t_theta(t_theta > pi/2) = pi/2 - t_theta(t_theta > pi/2); % flip the angle so we get the correct LUL value
        end
        t_LULi = max(min(floor(t_theta ./ LULAngleStep)+1,LUL.num),1);

        t_specA = LUL.A(t_LULi);
        t_specR = LUL.R(t_LULi);
    else
        t_specA(:) = t_Iabs;
        t_specR(:) = t_Irefl;
    end
         
    %% Calculate per-triangle incident power, specular optical forces, & absorption heat-loading
    t_Ipwr = abs(t_aproj.*t_IMag);  % incident power on each triangle, watts (includes costheta and broken factors)

    % Calculate specular absorbed power in each triangle.  This gives the absorption heat input.  
    % Note: Radiatiave cooling occurs later at the node level. That's why triangle area is distributed back 
    % to the nodes during each simulation step.  
    t_abs = t_specA .* t_Ipwr; % Watt
    % Note:  This value will be overwitten in textured regions via the LUT code

    % Calculate specular reflected power
    t_reflp = t_Ipwr .* t_specR;

    % calculate optical force ("optical thrust") for each triangle
    t_oth = 2 * t_reflp ./ c0 .* t_norms .* t_costheta  + ...  %  this is the reflection force
        t_abs ./ c0 .* t_IDir ;    % and this is the absorption force
    % Note: this will be overwritten in textured regions with LUT values.
end
  
% keep track of time spent doing triangle-based calculations (and specular optics)
timers(2) = timers(2) + toc;
tic;
timerStartRT = tic;


%% Step 4: Do ray tracing and/or calculate non-specular optics  
if RayTracingMode && (pwrRampVal > 0)
    if dooptics
    
    t_refln = t_IDir - 2 .* t_costheta .* t_norms;
    t_vert0 = [n_x(t_na); n_y(t_na); n_z(t_na)];
    
    processMR_GPU_v17;  % call raytracing subscript.  Generates:
          %  t_oth_mr :  additional optical thrust from multiple reflections
          %  t_abs_mr :  additional absorption from multiple reflections
    
%     db_bm_timer_gpu = db_bm_timer_gpu + db_bm_timer_gpu_val;
%     db_bm_timer_vec = db_bm_timer_vec + db_bm_timer_vec_val;
%     db_bm_timer_itr = db_bm_timer_itr + db_bm_timer_itr_val;
 
    timers(14) = timers(14) + toc(timerStartRT); % the RT code sometimes calls tic() on its own, so need to explicitly reference the toc() to the correct time.
    tic;
    
%                 RT_bmechmark_count = RT_bmechmark_count + 1;
%                 if RT_bmechmark_count == 301    % 123456789012                          Benchmark timers:  
%                     fprintf('Benchmark timers:      GPU    \t        VEC    \t   ITR   (numtris=%d)  \n                   ', numtris);
%                     fprintf('%11f \t', [ db_bm_timer_gpu db_bm_timer_vec db_bm_timer_itr]./300    );
%                     return;
%                 end
    end
          
else
    if dooptics
        if read_LUT

            %% 4b. Calculate optical response of non-specular ('textured') regions via look-up tables (LUTs).  
            % Look-up tables are used to specify arbitrary optical properties of textured regions of the mesh.  For
            % those regions, we look up pre-calculated values for absorption and optical thrust, based on the 
            % incidence angle relative to the triangle surface and texture orientation.  
            % For these calculations we are relying on the following vectors as the starting point. 
            %     t_norms :  unit normal vectors for each triangle
            %     t_texn :   unit vectors aligned with the grating ('parallel' axis)
            %     t_Idir:  unit vectors describing light incidence direction
            
            % First, we generate the transverse unit vectors for each triangle
            t_textn = cross(t_norms, t_texn);  % this is the transverse unit vector
            
            % Sanity checks, we can probably delete this, but it's nice to know if we've messed up:
            if mod(nt,1000)==0  % to save time, only do this check once every 1000 loops 
%                 disp([ 'This should be zero: ' sprintf('%e', max(abs( dot( t_norms, t_texn) ) ) ) ] );
                % Check for orthagonolity of t_norms and t_texn
                if (max(abs( dot( t_norms, t_texn) ) ) > 1e-15)
                    warning('Error with vector math in optical force LUT code ... the triangle normals are not orthogonal to the texture vectors!');
                end
    
%                 disp([ 'This should be zero: ' sprintf('%e', max(abs( dot( t_textn, t_textn ) - 1) ) ) ]);
                % Check that t_textn comprises unit vectors
                if (max(abs( dot( t_textn, t_textn ) - 1) ) > 1e-14)
                    warning('Error with vector math in optical force LUT code...  the transverse texture vector is not normalized. Something is wrong.');
                end
            end
            
            % Now we're going to project the incidence vector into the parallel and transverse grating planes.
            
            t_incsignt = dot( t_IDir, t_textn );  % projection mag along transverse vector
            t_incprojp = t_IDir - t_textn .* t_incsignt ; % inc vector projected into parallel plane
            t_incprojpn = t_incprojp ./ vecnorm(t_incprojp) ; % now normalized
        
            t_incsignp = dot( t_IDir, t_texn );  % projection mag along parallel vector
            t_incprojt = t_IDir - t_texn .* t_incsignp ; % inc vector projected into transverse plane
            t_incprojtn = t_incprojt ./ vecnorm(t_incprojt) ; % now normalized
        
            t_incsignt = sign(t_incsignt);  % get sign for angle unwrap from acos function
            t_incsignp = sign(t_incsignp);  % ...
        
            t_psi = t_incsignp .* acos( max(min(dot( t_incprojpn, t_norms ), 1),-1) );   % the min/max confine us to real limits for acos().  
            t_theta = t_incsignt .* acos( max(min(dot( t_incprojtn, t_norms ), 1),-1) ); % they shouldn't be needed, but I don't trust floating
                                                                                         % point arithmatic quite enough to omit them
        
            t_ipsi = round((t_psi+pi/2)./deg2rad(LUTPsiStepDeg)); % Psi index for LUTs
            t_itheta = round((t_theta+pi/2)./deg2rad(LUTThetaStepDeg)); % Theta index for LUTs
            
            % Check for (and correct) out-of-bounds psi and theta values, which occur for reverse-illuminated triangles
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
                if noobpl 
                    t_ipsi(outofboundspsilow) = 1; 
                end
                if noobph 
                    t_ipsi(outofboundspsihigh) = 180/LUTPsiStepDeg+1; 
                end
                if noobtl 
                    t_itheta(outofboundsthetalow) = 1; 
                end
                if noobth 
                    t_itheta(outofboundsthetahigh) = 180/LUTThetaStepDeg+1; 
                end
            end
            
            % Actually perform the look-ups:
            for nlut = 1:length(LUTs)
                lutinds = sub2ind(size(LUTs(nlut).n), t_ipsi(t_tex == nlut),t_itheta(t_tex == nlut));
                if any(lutinds)
                    t_press_n(t_tex == nlut) = LUTs(nlut).n(lutinds);
                    t_press_t(t_tex == nlut) = LUTs(nlut).t(lutinds);
                    t_press_p(t_tex == nlut) = LUTs(nlut).p(lutinds);
                    t_lut_a( t_tex == nlut ) = LUTs(nlut).a(lutinds);
                end
            end
            
            % Calculate optical thrust and absorbed power based on LUT values
            lutmask = t_tex>0;
            if any(lutmask)
                t_oth(:,lutmask) = t_Ipwr(lutmask) .*  ( ...
                    t_norms(:,lutmask) .* t_press_n(lutmask)  + ... this is the normal force
                    t_textn(:,lutmask) .* t_press_t(lutmask)  + ... this is the transverse force
                    t_texn(:,lutmask) .* t_press_p(lutmask)); % this is the parallel/axial force
                t_abs(lutmask) = t_lut_a(lutmask) .* t_Ipwr(lutmask);
            end

        elseif ~read_LUT

            %% 4c. Calculate optical response of non-specular ('textured') regions with Ramon's explicit, but slower approach  
    
            % Pre-calculate terms used multiple times
            term_rot_cross_1 = ...
                t_texn(2,:) .* t_norms(3,:) - t_texn(3,:) .* t_norms(2,:);
            
            term_rot_cross_2 = ...
                t_texn(3,:) .* t_norms(1,:) - t_texn(1,:) .* t_norms(3,:);
            
            term_rot_cross_3 = ...
                t_texn(1,:) .* t_norms(2,:) - t_texn(2,:) .* t_norms(1,:);
            
            term_cross_1 = ...
                temp_tex_proj_norms(2,:) .* temp_tri_normvec(3,:) - ...
                temp_tex_proj_norms(3,:) .* temp_tri_normvec(2,:);
            
            term_cross_2 = ...
                temp_tex_proj_norms(3,:) .* temp_tri_normvec(1,:) - ...
                temp_tex_proj_norms(1,:) .* temp_tri_normvec(3,:);                
            
            term_cross_3 = ...
                temp_tex_proj_norms(1,:) .* temp_tri_normvec(2,:) - ...
                temp_tex_proj_norms(2,:) .* temp_tri_normvec(1,:);  
                  
            % Calculate selected entries of rotation matrix
            rotmat_11 = ...
                t_texn(1,:) .* temp_tex_proj_norms(1,:) + ...
                t_norms(1,:) .* temp_tri_normvec(1,:)  + ...
                term_rot_cross_1 .* term_cross_1;
            
            rotmat_12 = ...
                t_texn(1,:) .* temp_tex_proj_norms(2,:) + ...
                t_norms(1,:) .* temp_tri_normvec(2,:)  + ...
                term_rot_cross_1 .* term_cross_2;                
            
            rotmat_13 = ...
                t_texn(1,:) .* temp_tex_proj_norms(3,:) + ...
                t_norms(1,:) .* temp_tri_normvec(3,:)  + ...
                term_rot_cross_1 .* term_cross_3;                   
            
            rotmat_23 = ...
                t_texn(2,:) .* temp_tex_proj_norms(3,:) + ...
                t_norms(2,:) .* temp_tri_normvec(3,:)  + ...
                term_rot_cross_2 .* term_cross_3;                 

            rotmat_33 = ...
                t_texn(3,:) .* temp_tex_proj_norms(3,:) + ...
                t_norms(3,:) .* temp_tri_normvec(3,:)  + ...
                term_rot_cross_3 .* term_cross_3;
                        
            % Calculate pitch, roll and yaw angles based on formulas specified in
            % Ramon's notes. Note that technically, the yaw angles are not
            % important, as we assume the polarization to rotate synchronously with
            % the sail and thus its MEPS. However, calculating the yaw_angles
            % allows us to track the effect of spinning the sail
        
            % Note that the additional minus signs for the pitch and roll angles
            % are there to account for Oggy's convention of sign of incidence angles
        
            pitch_angles = asin(rotmat_13);
            roll_angles = -atan2(rotmat_23, rotmat_33);
            yaw_angles = -atan2(rotmat_12, rotmat_11);
    
            % Allocate memory for current roll and pitch angles of each triangle
            specific_angles = zeros(max(size(t_tex)),2);
            specific_angles_indices = zeros(max(size(t_tex)),num_small_tables);
        
            % Assemble correct definition of angles for different regions. Here,
            % regions 3 & 4 are mirror-symmetric with respect to region 1 & 6, and
            % region 5 is mirror-symmetric with respect to region 2. Therefore,
            % roll and pitch angles for regions 3, 4 and 5 need to have opposite
            % signs compared to the roll and pitch angles for regions 1, 2 and 6,
            % since only regions 1 and 2 have been simulated (region 6 is equal to
            % region 1, region 3 is equal to region 4)
        
            specific_angles(indices_Region1,:) = ...
                [roll_angles(indices_Region1)', pitch_angles(indices_Region1)'];
            specific_angles(indices_Region2,:) = ...
                [roll_angles(indices_Region2)', pitch_angles(indices_Region2)'];
        
            specific_angles(indices_Region3,:) = ...
                [-roll_angles(indices_Region3)', -pitch_angles(indices_Region3)'];
            specific_angles(indices_Region4,:) = ...
                [-roll_angles(indices_Region4)', -pitch_angles(indices_Region4)'];
        
            % Get indices indicating in which range of roll angles every angle set
            % (roll, pitch) for each triangle can be found. This is needed later to
            % retrieve the pressures from []the correct sub-look-up table
        
            for ii = 1:1:num_small_tables
                if ii == 1
                    specific_angles_indices(:,ii) =  ...
                        specific_angles(:,1) < border_angles(1);
                elseif ii == num_small_tables
                    specific_angles_indices(:,ii) = ...
                        specific_angles(:,1) >= border_angles(end); 
                elseif (ii > 1) && (ii < num_small_tables)
                    specific_angles_indices(:,ii) = ...
                        (specific_angles(:,1) >= border_angles(ii-1)) & ...
                        (specific_angles(:,1) < border_angles(ii));
                else
                    warning('Error: angle out of range!');
                end
            end
        
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
        
            for ii = 1:1:num_small_tables
                
                logical_specific_angles_indices = ...
                    logical(specific_angles_indices(:,ii));
                    
                if  ~isempty(specific_angles(logical_specific_angles_indices))
        
                    [indices_angles, ~] = knnsearch(pressure_TE_small_tables(1:last_indices_tables(ii),1:2,ii), ...
                        specific_angles(logical_specific_angles_indices,:));
        
                    selected_pressures_x_TE(logical(specific_angles_indices(:,ii))) = ...
                        c0*pressure_TE_small_tables(indices_angles,3,ii);
                    selected_pressures_z_TE(logical(specific_angles_indices(:,ii))) = ...
                        c0*pressure_TE_small_tables(indices_angles,4,ii);
                    selected_pressures_x_TM(logical(specific_angles_indices(:,ii))) = ...
                        c0*pressure_TM_small_tables(indices_angles,3,ii);
                    selected_pressures_z_TM(logical(specific_angles_indices(:,ii))) = ...
                        c0*pressure_TM_small_tables(indices_angles,4,ii);
        
                end
            end
        
            % Transform local frames of individual regions to global body frame.
            % Note that the lightsail is flipped upside down for this simulation
            % compared to those in COMSOL, thus requiring appropriate sign changes
            % for respective pressures
        
            t_press_x(indices_Region1) = selected_pressures_x_TE(indices_Region1);
            t_press_y(indices_Region1) = 0;
            t_press_z(indices_Region1) = selected_pressures_z_TE(indices_Region1);
        
            t_press_x(indices_Region2) = 0;
            t_press_y(indices_Region2) = selected_pressures_x_TM(indices_Region2);
            t_press_z(indices_Region2) = selected_pressures_z_TM(indices_Region2);
        
            t_press_x(indices_Region3) = -selected_pressures_x_TE(indices_Region3);
            t_press_y(indices_Region3) = 0;
            t_press_z(indices_Region3) = selected_pressures_z_TE(indices_Region3);
        
            t_press_x(indices_Region4) = 0;
            t_press_y(indices_Region4) = -selected_pressures_x_TM(indices_Region4);
            t_press_z(indices_Region4) = selected_pressures_z_TM(indices_Region4);
        
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
        
            % Transform forces on triangle in local body frame to global laser
            % frame by multiplying body-frame force vector with the direction
            % cosine matrix H_I^B(apparently NOT the transpose, although I don't
            % understand?) according to OIlic2019
        
            t_force_x = t_force_x_BF .* DRM123_11(yaw_angles,pitch_angles,roll_angles) + ...
                t_force_y_BF .* DRM123_21(yaw_angles,pitch_angles,roll_angles) + ...
                t_force_z_BF .* DRM123_31(yaw_angles,pitch_angles,roll_angles);
        
            t_force_y = t_force_x_BF .* DRM123_12(yaw_angles,pitch_angles,roll_angles) + ...
                t_force_y_BF .* DRM123_22(yaw_angles,pitch_angles,roll_angles) + ...
                t_force_z_BF .* DRM123_32(yaw_angles,pitch_angles,roll_angles);
        
            t_force_z = t_force_x_BF .* DRM123_13(yaw_angles,pitch_angles,roll_angles) + ...
                t_force_y_BF .* DRM123_23(yaw_angles,pitch_angles,roll_angles) + ...
                t_force_z_BF .* DRM123_33(yaw_angles,pitch_angles,roll_angles);
        
            t_oth(1,:) =  t_force_x;
            t_oth(2,:) =  t_force_y;
            t_oth(3,:) =  t_force_z;
        
        end

        % Keep track of time spent on LUT calculations
        timers(3) = timers(3) + toc;
        tic;       
    end
end


%% Step 5: Distribute optical forces (and absorbed heat) to the nodes
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
    n_a = t_a * M_t2n ./ 3;
    if (RayTracingMode)
        n_abs = (t_abs+t_abs_mr) * M_t2n ./ 3;
    else
        n_abs = t_abs * M_t2n ./ 3;
    end

    % quick sanity check, if we've messed up a calulation, terminate the simulation immediately so we can debug
    % before the error corrupts the entire mesh.
    if (isnan(sum(n_abs)))
        error('Found NaN!!!');
    end
end

% Keep track of time spent on triangle-to-node calculations
timers(4) = timers(4) + toc;
tic;


%% Step 6: Calculate thermal emission (radiatively emitted power)
% Calculate radiated power using Stefan-Boltzmann's law given by P_rad
% = A * emissivity * sigma * ( T^4 - T_ref^4). Here, the radiated power
% at each node is calculated.

% n_ems = SBCmm .* n_a .* Emissivity .* (n_t.^4); % Watt
n_ems = SBCmm .* n_a .* Emissivity_Si3N4 .* (n_t.^4); % Watt


%% Step 7: Calculate mechanical forces & heat conduction

% Matrix containing all edge direction vectors with x, y, z components as columns
e_nl = [ n_x(e_nb) - n_x(e_na) ; n_y(e_nb) - n_y(e_na) ; n_z(e_nb) - n_z(e_na) ];

% Calculate current lengths in mm of all edges
e_l = vecnorm( e_nl);

% Normalize all edge direction vectors. Here, 'e_nl' stands for edge, normalized length
e_nl = e_nl ./ e_l;   % unit vectors

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
e_s = e_l ./ (e_l0t) - 1;   % strain %/100.   No stress -> 0 strain. Elongate to 2x length -> 1.0 strain

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
e_dt = n_t(e_nb) - n_t(e_na);   % Kelvin

% Thermal conduction, i.e. heat flow power on edges ('e_hf') according to formula kappa * (A/l) * (T_hot - T_cold)
e_hf = e_heatCondc .* e_dt;  % Watts

% Terminate script if we've picked up a math error somehow, so we can debug before it corrupts the mesh.
if isnan(sum(e_hf))
    error('Found NaN!')
end

% Keep track of time spent doing edge-based calculations
timers(5) = timers(5) + toc;
tic;


%% Step 8: Assign mechanical forces to nodes (from edges) and also heat flow
n_mf(1,:) = e_mf(1,:) * M_e2n;
n_mf(2,:) = e_mf(2,:) * M_e2n;
n_mf(3,:) = e_mf(3,:) * M_e2n;

% Calculate node heat flow from edge heat flow
n_hf = e_hf * M_e2n;

% Keep track of time spent doing edge-to-node calculations
timers(6) = timers(6) + toc;
tic;


%% Step 9: Evaluate for tensile or thermal failure
% DISCLAIMER: This simulator does not attempt to accurately model the evolution of tensile failure, other
% than to detect the onset of such failure whenever the strain on a mesh edge exceeds the tensile limit of
% the material. It would be prudent to immediately terminate the simulation at the first such failure...
%
% but...
%
% Letting the simulator continue to run after tensile failure turns out to produce fun videos of the
% lightsail tearing itself apart and flailing around in pieces! So, to keep ourselves entertained, we can
% crudely (and with no intent or claim of accuracy) continue to simulate the partially failed lightsail, by
% essentially omitting the failed edges and triangles from optical and mechanical calculations, as if they
% weren't there in the first place.  We denote such failed elements as 'broken'. Here's where we detect
% and keep track of broken edges and triangles:

if any_flex  % we can skip tensile failure analysis if the entire mesh is rigid...

    e_tensilefailure =  ( e_flexible &  e_notbroken) & ( e_s > e_tensile_strain_limit )  ;  %all edges that are now broken from tensile failure
    e_thermalfailure = ( e_notbroken .* (e_t > e_failtemp ) );
    newtensilefailures = sum(e_tensilefailure);
    newthermalfailures = sum(e_thermalfailure);
    e_broken = e_broken | (e_tensilefailure > 0) | e_thermalfailure;
                    
    % Count number of broken bonds/edges and call it 'uhohs'
    tensilefailures = tensilefailures + newtensilefailures;
    thermalfailures = thermalfailures + newthermalfailures;
    
    uhohs = sum(e_broken);
    
    % Determine which triangles are broken ... we used to delete triangles if any single edge failed, ...
    t_broken =  ( e_broken(t_e1) | e_broken(t_e2) | e_broken(t_e3) );  % and that's still the best I could come up with
    
    % Convenience variables so we don't have to perform negation every time we want the unbroken edges/tris
    e_notbroken = ~e_broken;
    t_notbroken = ~t_broken;
    
    % 'any' command returns '1' if any of the vector elements is nonzero.
    % 'anybroken' is a number that is '1' if there is at least one broken
    % bond in the mesh
    anybroken = anybroken || (any(t_broken));  % whether or not there has ever been a failure
    
    % Record time and time-index if the first breakage has just occurred
    if (anybroken) && (nto_firstBroken == sim_dur_nto)
        nto_firstBroken = nto;
        time_frstBrkn = tt;
        
        if terminate_on_breakup
            if ~terminate_t
                terminate_t = tt + termination_delay_s;
                if ~terminate_t
                    terminate_t = dt;  % just in case the first breakage would try to schedule termination at exactly tt=0, we'll change it to a non-zero value, because zero is a sentinel value that prevents termination...
                end
            end
        end
    end
    
    % Keep track of time spent on edge breakage detection and bookkeeping
    timers(7) = timers(7) + toc;
    tic;

end


%% Step 10: Mechanical calculations and time stepping

% We calculate external forces (not including elastic mesh forces) as an intermediate step
n_rfx = n_of(1,:) + n_af(1,:);
n_rfy = n_of(2,:) + n_af(2,:);
n_rfz = n_of(3,:) + n_af(3,:);

externalForceX = sum( n_rfx );
externalForceY = sum( n_rfy );
externalForceZ = sum( n_rfz );

externalForce = [ externalForceX; externalForceY; externalForceZ];

% We can calculate an effective torque on the whole mesh (as if it were rigid).
% This is not actually used in subsequent calculations...
externalTorqueX = sum( n_dy.*n_rfz - n_dz.*n_rfy );
externalTorqueY = sum( n_dz.*n_rfx - n_dx.*n_rfz );
externalTorqueZ = sum( n_dx.*n_rfy - n_dy.*n_rfx );

externalTorque = [ externalTorqueX; externalTorqueY; externalTorqueZ ];

% Sum mechanical and optical forces on each node
n_fx = n_mf(1,:) + n_rfx;
n_fy = n_mf(2,:) + n_rfy;
n_fz = n_mf(3,:) + n_rfz;
% Note that n_mf doesn't include edge spring contributions from within rigid bodies


%% Step 11: Update triangle centroids
t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;


%% Step 12: Update temperatures
% Heat capacity links the temperature rise of a body to the required
% amount of energy. The unit of heat capacity is [J/g/degC], which
% describes the amount of energy required to raise the temperature of 1
% gram of an object by 1 deg C (or 1 Kelvin). Energy is equal to power
% times time, such that the temperature rise is given by energy / mass
% / heat capacity, where energy = (input power - output power) * dt.
% Input power is given by the opticaly absorbed power and power from
% heat flow/thermal conduction, whereas output power is due to
% radiative cooling.

% MOVED TO MAIN SCRIPT

% n_t = n_t + (n_hf + n_abs - n_ems).*dt ./ n_heatmass;   % todo:  need to support variable heat cap across the array
% 
% % Calculate the edge temperature as the averaged temperature of the sum
% % of the temperatures of the two nodes that form the edge
% e_t = (n_t(e_na) + n_t(e_nb)) ./ 2;
% 
% % Calculate temperature of triangle as the average of the temperatures
% % of the nodes that make up the triangle
% t_t  =  ( n_t(t_na) + n_t(t_nb) + n_t(t_nc) ) ./ 3;

    
