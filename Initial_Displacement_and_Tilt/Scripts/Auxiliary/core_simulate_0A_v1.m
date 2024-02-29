%% Core Computation Script for STARSHOT LIGHTSAIL SIMULATOR
% (c) 2019-2023 Michael D. Kelzenberg, Ramon Gao, California Institute of Technology

% Do time-dependent optical, mechanical and thermal calculations here
% Outsourced from main code for all versions 19 and earlier to allow for
% use of more advanced numerical integration methods

% Version 1, May 14th, 2023

%% Step 0: Calculate spin-up forces
if  (tt_tmp < 0) && rps_target  % spin-up forces are applied for tt<0.  
    n_arad = sqrt( n_x_tmp.^2 + n_y_tmp.^2 );
    n_af_tmp(1,:) = sign(rps_target)*spinup_gs/88.2.*n_m.*(n_y_tmp) ;%%% .* (n_arad < 0.6 * radiusmm);
    n_af_tmp(2,:) = sign(rps_target)*-spinup_gs/88.2*n_m.*(n_x_tmp) ;%%%.* (n_arad < 0.6 * radiusmm);
end

if ( tt_tmp >= 0) || anybroken  % stop spin-up at t=0 (or if any edges broke during spinup)
    n_af_tmp(:) = 0;
end


%% Step 1: Calculate beam intensity profile and direction
% Here, we set t_IMag and t_IDir, which are the magnitude and direction of the incident light at each
% triangle centroid. Normally, t_IDir is (0, 0, 1), so we never change it here, but in theory, you could set
% whatever value you want. Note: t_IDir must be unit vectors.

pwrRampVal_tmp = pwrRampFcn(tt_tmp);

if (pwrRampVal_tmp == 0)
    t_IMag_tmp = zeros(size(t_cx_tmp));
else
    t_optFcn_r2_tmp = (t_cx_tmp - 0).^2 + (t_cy_tmp - 0).^2;
    t_IMag_tmp = pwrRampVal_tmp * I0 * exp( - t_optFcn_r2_tmp ./ (Irad^2) );  % This is in W/mm2.  Divide by c0 to get newtons / mm2 radiation pressure (N/mm2)
  
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
    t_ev1_tmp(1,:) = n_x_tmp(t_nb) - n_x_tmp(t_na);
    t_ev1_tmp(2,:) = n_y_tmp(t_nb) - n_y_tmp(t_na);
    t_ev1_tmp(3,:) = n_z_tmp(t_nb) - n_z_tmp(t_na);
    t_ev2_tmp(1,:) = n_x_tmp(t_nc) - n_x_tmp(t_na);
    t_ev2_tmp(2,:) = n_y_tmp(t_nc) - n_y_tmp(t_na);
    t_ev2_tmp(3,:) = n_z_tmp(t_nc) - n_z_tmp(t_na);
else
    t_ev1_tmp = [n_x_tmp; n_y_tmp; n_z_tmp]* M_nxyz2tev1;
    t_ev2_tmp = [n_x_tmp; n_y_tmp; n_z_tmp]* M_nxyz2tev2;
end

if read_LUT

    t_norms_tmp = cross( t_ev1_tmp, t_ev2_tmp );  % triangle normal vectors, not yet normalized; we will normalize shortly

    % We determine if a triangle is upside down based on the z component of
    % its normal vector.  Here, there's an implicit assumption that light is incident along the 
    % positive Z axis, so if you ever set t_IDir to something other than (0, 0, 1), you might need to
    % determine a new way to detect upside down triangles.
    upsidedowntris_tmp = sum( t_norms_tmp(3,:) < 0);   

elseif ~read_LUT

    t_norms_tmp = cross( t_ev2_tmp, t_ev1_tmp );  % for propulsion in -z - sorry, I (Ramon) still struggle to figure out correct calculation of forces in +z
    upsidedowntris_tmp(nt) = sum( t_norms_tmp(3,:) > 0); 

end

if upsidedowntris_tmp
    if ~terminate_t_tmp && (terminate_on_upsidedown)  % initiate delayed termination if triangles are upside down
        terminate_t_tmp = tt_tmp + termination_delay_s;
        disp('Initiating delayed termination due to upside down triangles!');
    end
end

% Calculate length of the non-normalized normal vectors, which happen to also be 2x the triangle area
t_a_tmp = max(vecnorm(t_norms_tmp), 1e-8);  % initially, this is 2x the triangle area

% Normalize normal vectors to each triangle
t_norms_tmp = t_norms_tmp ./ t_a_tmp ; % ok, now t_norms contains proper unit normal vectors

t_a_tmp = 0.5 .* t_a_tmp; % ok, now t_a contains the actual area of each triangle

% If any triangles are 'broken', we replace their now-meanlingless calculated area with their original area
if anybroken
    t_a_tmp(t_broken) = t_a0(t_broken);  % Don't worry, other calculations will continue ignoring broken triangles.  
end

% Normalize the #1 edge vectors
t_ev1n_tmp = t_ev1_tmp./ vecnorm(t_ev1_tmp);%test;

% During mesh generation, we recorded a value 't_texa' for each triangle, which is the in-plane angle
% between the texture orientation vector and the #1 edge vector.  We use that now to calculate the 
% new texture orientation vectors for the current mesh state... well, actually, we use the cached values 
% of the sin and cosine of this angle to hopefully save time (t_costexa, t_sintexa).
t_texn_tmp = t_ev1n_tmp .* (t_costexa) + ...
    cross(t_norms_tmp, t_ev1n_tmp ) .* (t_sintexa)  + ...
    t_norms_tmp .* dot(t_norms_tmp, t_ev1n_tmp ) .* (1 - t_costexa) ;  %TODO: Can I safely replace the cos() with t_costexa???

% Calculate cos(theta), where theta corresponds to the light incidence angle at each triangle.              
t_costheta_tmp = dot(t_norms_tmp, t_IDir);

% Here we calculate the projected area of non-broken triangles in the X-Y plane.  This is essentially the
% aperture area of our sail in the plane of beam incidence.  Note again the implicit assumption that
% illumination is along (0,0,1).  
t_aproj_tmp = t_a_tmp .* t_costheta_tmp .* (t_notbroken);


%% Step 3: Calculate reflection and absorption for specular regions
%  Currently we use simple scalar values for reflection and absorption.  If you were so inclined, you could
%  implement different physics to calculate how much light should be reflected and absorbed by specular sail
%  surfaces, such as using Fresnel equations.  However, to implement more complex optical behavior, I 
%  recommend using look-up tables rather than trying to calculate values here.
%
%  Available inputs for calculating specular response:
%  t_costheta: cosine of angle of incidence
%  t_t_tmp: local temperature
%  t_x0, t_y0, t_z0: Initial coordinates of this triangle (prior to simulation beginning)
%  t_thick: thickness of film here 
%  t_m: mass of this triangle
%
%  Required outputs from this section: 
%  t_specA, t_specR, absorption and reflection coefficients (coefficients--not power!  Values should be 0 to 1)
    
if dooptics
    if useSpecularTMM
        t_theta_tmp = acos(t_costheta_tmp);
        if min(t_theta_tmp) < 0 
            warning('Negative theta in useSpecularTMM code');
            t_theta_tmp(t_theta_tmp < 0) = 0;  % I don't think this should ever happen
        end
        if max(t_theta_tmp) > pi/2
        %    warning('Reverse illuminated triangles in specularTMM code');
            t_theta_tmp(t_theta_tmp > pi/2) = pi/2 - t_theta_tmp(t_theta_tmp > pi/2); % flip the angle so we get the correct LUL value
        end
        t_LULi_tmp = max(min(floor(t_theta_tmp ./ LULAngleStep)+1,LUL.num),1);

        t_specA_tmp = LUL.A(t_LULi_tmp);
        t_specR_tmp = LUL.R(t_LULi_tmp);
    else
        t_specA_tmp(:) = t_Iabs;
        t_specR_tmp(:) = t_Irefl;
    end
         
    %% Calculate per-triangle incident power, specular optical forces, & absorption heat-loading
    t_Ipwr_tmp = abs(t_aproj_tmp.*t_IMag_tmp);  % incident power on each triangle, watts (includes costheta and broken factors)

    % Calculate specular absorbed power in each triangle.  This gives the absorption heat input.  
    % Note: Radiatiave cooling occurs later at the node level. That's why triangle area is distributed back 
    % to the nodes during each simulation step.  
    t_abs_tmp = t_specA_tmp .* t_Ipwr_tmp; % Watt
    % Note:  This value will be overwitten in textured regions via the LUT code

    % Calculate specular reflected power
    t_reflp_tmp = t_Ipwr_tmp .* t_specR_tmp;

    % calculate optical force ("optical thrust") for each triangle
    t_oth_tmp = 2 * t_reflp_tmp ./ c0 .* t_norms_tmp .* t_costheta_tmp  + ...  %  this is the reflection force
        t_abs_tmp ./ c0 .* t_IDir ;    % and this is the absorption force
    % Note: this will be overwritten in textured regions with LUT values.
end
  
% keep track of time spent doing triangle-based calculations (and specular optics)
timers(2) = timers(2) + toc;
tic;
timerStartRT = tic;


%% Step 4: Do ray tracing and/or calculate non-specular optics  
if RayTracingMode && (pwrRampVal_tmp > 0)
    if dooptics
    
    t_refln_tmp = t_IDir - 2 .* t_costheta_tmp .* t_norms_tmp;
    t_vert0_tmp = [n_x_tmp(t_na); n_y_tmp(t_na); n_z_tmp(t_na)];
    
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
            t_textn_tmp = cross(t_norms_tmp, t_texn_tmp);  % this is the transverse unit vector
            
            % Sanity checks, we can probably delete this, but it's nice to know if we've messed up:
            if mod(nt,1000)==0  % to save time, only do this check once every 1000 loops 
%                 disp([ 'This should be zero: ' sprintf('%e', max(abs( dot( t_norms, t_texn) ) ) ) ] );
                % Check for orthagonolity of t_norms and t_texn
                if (max(abs( dot( t_norms_tmp, t_texn_tmp) ) ) > 1e-15)
                    warning('Error with vector math in optical force LUT code ... the triangle normals are not orthogonal to the texture vectors!');
                end
    
%                 disp([ 'This should be zero: ' sprintf('%e', max(abs( dot( t_textn, t_textn ) - 1) ) ) ]);
                % Check that t_textn comprises unit vectors
                if (max(abs( dot( t_textn_tmp, t_textn_tmp ) - 1) ) > 1e-14)
                    warning('Error with vector math in optical force LUT code...  the transverse texture vector is not normalized. Something is wrong.');
                end
            end
            
            % Now we're going to project the incidence vector into the parallel and transverse grating planes.
            
            t_incsignt_tmp = dot( t_IDir, t_textn_tmp );  % projection mag along transverse vector
            t_incprojp_tmp = t_IDir - t_textn_tmp .* t_incsignt_tmp ; % inc vector projected into parallel plane
            t_incprojpn_tmp = t_incprojp_tmp ./ vecnorm(t_incprojp_tmp) ; % now normalized
        
            t_incsignp_tmp = dot( t_IDir, t_texn_tmp );  % projection mag along parallel vector
            t_incprojt_tmp = t_IDir - t_texn_tmp .* t_incsignp_tmp ; % inc vector projected into transverse plane
            t_incprojtn_tmp = t_incprojt_tmp ./ vecnorm(t_incprojt_tmp) ; % now normalized
        
            t_incsignt_tmp = sign(t_incsignt_tmp);  % get sign for angle unwrap from acos function
            t_incsignp_tmp = sign(t_incsignp_tmp);  % ...
        
            t_psi_tmp = t_incsignp_tmp .* acos( max(min(dot( t_incprojpn_tmp, t_norms_tmp ), 1),-1) );   % the min/max confine us to real limits for acos().  
            t_theta_tmp = t_incsignt_tmp .* acos( max(min(dot( t_incprojtn_tmp, t_norms_tmp ), 1),-1) ); % they shouldn't be needed, but I don't trust floating
                                                                                         % point arithmatic quite enough to omit them
        
            t_ipsi_tmp = round((t_psi_tmp+pi/2)./deg2rad(LUTPsiStepDeg)); % Psi index for LUTs
            t_itheta_tmp = round((t_theta_tmp+pi/2)./deg2rad(LUTThetaStepDeg)); % Theta index for LUTs
            
            % Check for (and correct) out-of-bounds psi and theta values, which occur for reverse-illuminated triangles
            outofboundspsilow = t_ipsi_tmp < 1;
            outofboundspsihigh = t_ipsi_tmp > (180/LUTPsiStepDeg+1);
            outofboundsthetalow = t_itheta_tmp < 1;
            outofboundsthetahigh = t_itheta_tmp > (180/LUTThetaStepDeg+1);
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
                    t_ipsi_tmp(outofboundspsilow) = 1; 
                end
                if noobph 
                    t_ipsi_tmp(outofboundspsihigh) = 180/LUTPsiStepDeg+1; 
                end
                if noobtl 
                    t_itheta_tmp(outofboundsthetalow) = 1; 
                end
                if noobth 
                    t_itheta_tmp(outofboundsthetahigh) = 180/LUTThetaStepDeg+1; 
                end
            end
            
            % Actually perform the look-ups:
            for nlut = 1:length(LUTs)
                lutinds_tmp = sub2ind(size(LUTs(nlut).n), t_ipsi_tmp(t_tex == nlut),t_itheta_tmp(t_tex == nlut));
                if any(lutinds_tmp)
                    t_press_n_tmp(t_tex == nlut) = LUTs(nlut).n(lutinds_tmp);
                    t_press_t_tmp(t_tex == nlut) = LUTs(nlut).t(lutinds_tmp);
                    t_press_p_tmp(t_tex == nlut) = LUTs(nlut).p(lutinds_tmp);
                    t_lut_a_tmp( t_tex == nlut ) = LUTs(nlut).a(lutinds_tmp);
                end
            end
            
            % Calculate optical thrust and absorbed power based on LUT values
            lutmask = t_tex>0;
            if any(lutmask)
                t_oth_tmp(:,lutmask) = t_Ipwr_tmp(lutmask) .*  ( ...
                    t_norms_tmp(:,lutmask) .* t_press_n_tmp(lutmask)  + ... this is the normal force
                    t_textn_tmp(:,lutmask) .* t_press_t_tmp(lutmask)  + ... this is the transverse force
                    t_texn_tmp(:,lutmask) .* t_press_p_tmp(lutmask)); % this is the parallel/axial force
                t_abs_tmp(lutmask) = t_lut_a_tmp(lutmask) .* t_Ipwr_tmp(lutmask);
            end

        elseif ~read_LUT

            %% 4c. Calculate optical response of non-specular ('textured') regions with Ramon's explicit, but slower approach  
    
            % Pre-calculate terms used multiple times
            term_rot_cross_1_tmp = ...
                t_texn_tmp(2,:) .* t_norms_tmp(3,:) - t_texn_tmp(3,:) .* t_norms_tmp(2,:);
            
            term_rot_cross_2_tmp = ...
                t_texn_tmp(3,:) .* t_norms_tmp(1,:) - t_texn_tmp(1,:) .* t_norms_tmp(3,:);
            
            term_rot_cross_3_tmp = ...
                t_texn_tmp(1,:) .* t_norms_tmp(2,:) - t_texn_tmp(2,:) .* t_norms_tmp(1,:);
            
            term_cross_1_tmp = ...
                temp_tex_proj_norms(2,:) .* temp_tri_normvec(3,:) - ...
                temp_tex_proj_norms(3,:) .* temp_tri_normvec(2,:);
            
            term_cross_2_tmp = ...
                temp_tex_proj_norms(3,:) .* temp_tri_normvec(1,:) - ...
                temp_tex_proj_norms(1,:) .* temp_tri_normvec(3,:);                
            
            term_cross_3_tmp = ...
                temp_tex_proj_norms(1,:) .* temp_tri_normvec(2,:) - ...
                temp_tex_proj_norms(2,:) .* temp_tri_normvec(1,:);  
                  
            % Calculate selected entries of rotation matrix
            rotmat_11_tmp = ...
                t_texn_tmp(1,:) .* temp_tex_proj_norms(1,:) + ...
                t_norms_tmp(1,:) .* temp_tri_normvec(1,:)  + ...
                term_rot_cross_1_tmp .* term_cross_1_tmp;
            
            rotmat_12_tmp = ...
                t_texn_tmp(1,:) .* temp_tex_proj_norms(2,:) + ...
                t_norms_tmp(1,:) .* temp_tri_normvec(2,:)  + ...
                term_rot_cross_1_tmp .* term_cross_2_tmp;                
            
            rotmat_13_tmp = ...
                t_texn_tmp(1,:) .* temp_tex_proj_norms(3,:) + ...
                t_norms_tmp(1,:) .* temp_tri_normvec(3,:)  + ...
                term_rot_cross_1_tmp .* term_cross_3_tmp;                   
            
            rotmat_23_tmp = ...
                t_texn_tmp(2,:) .* temp_tex_proj_norms(3,:) + ...
                t_norms_tmp(2,:) .* temp_tri_normvec(3,:)  + ...
                term_rot_cross_2_tmp .* term_cross_3_tmp;                 

            rotmat_33_tmp = ...
                t_texn_tmp(3,:) .* temp_tex_proj_norms(3,:) + ...
                t_norms_tmp(3,:) .* temp_tri_normvec(3,:)  + ...
                term_rot_cross_3_tmp .* term_cross_3_tmp;

            % Calculate pitch, roll and yaw angles based on formulas specified in
            % Ramon's notes. Note that technically, the yaw angles are not
            % important, as we assume the polarization to rotate synchronously with
            % the sail and thus its MEPS. However, calculating the yaw_angles
            % allows us to track the effect of spinning the sail
        
            % Note that the additional minus signs for the pitch and roll angles
            % are there to account for Oggy's convention of sign of incidence angles
        
            pitch_angles_tmp = asin(rotmat_13_tmp);
            roll_angles_tmp = -atan2(rotmat_23_tmp, rotmat_33_tmp);
            yaw_angles_tmp = -atan2(rotmat_12_tmp, rotmat_11_tmp);
    
            % Allocate memory for current roll and pitch angles of each triangle
            specific_angles_tmp = zeros(max(size(t_tex)),2);
            specific_angles_indices_tmp = zeros(max(size(t_tex)),num_small_tables);
        
            % Assemble correct definition of angles for different regions. Here,
            % regions 3 & 4 are mirror-symmetric with respect to region 1 & 6, and
            % region 5 is mirror-symmetric with respect to region 2. Therefore,
            % roll and pitch angles for regions 3, 4 and 5 need to have opposite
            % signs compared to the roll and pitch angles for regions 1, 2 and 6,
            % since only regions 1 and 2 have been simulated (region 6 is equal to
            % region 1, region 3 is equal to region 4)
        
            specific_angles_tmp(indices_Region1,:) = ...
                [roll_angles_tmp(indices_Region1)', pitch_angles_tmp(indices_Region1)'];
            specific_angles_tmp(indices_Region2,:) = ...
                [roll_angles_tmp(indices_Region2)', pitch_angles_tmp(indices_Region2)'];
        
            specific_angles_tmp(indices_Region3,:) = ...
                [-roll_angles_tmp(indices_Region3)', -pitch_angles_tmp(indices_Region3)'];
            specific_angles_tmp(indices_Region4,:) = ...
                [-roll_angles_tmp(indices_Region4)', -pitch_angles_tmp(indices_Region4)'];
        
            % Get indices indicating in which range of roll angles every angle set
            % (roll, pitch) for each triangle can be found. This is needed later to
            % retrieve the pressures from []the correct sub-look-up table
        
            for ii = 1:1:num_small_tables
                if ii == 1
                    specific_angles_indices_tmp(:,ii) =  ...
                        specific_angles_tmp(:,1) < border_angles(1);
                elseif ii == num_small_tables
                    specific_angles_indices_tmp(:,ii) = ...
                        specific_angles_tmp(:,1) >= border_angles(end); 
                elseif (ii > 1) && (ii < num_small_tables)
                    specific_angles_indices_tmp(:,ii) = ...
                        (specific_angles_tmp(:,1) >= border_angles(ii-1)) & ...
                        (specific_angles_tmp(:,1) < border_angles(ii));
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
        
            selected_pressures_x_TE_tmp = zeros(size(t_press_x_tmp));
            selected_pressures_z_TE_tmp = zeros(size(t_press_x_tmp));
            selected_pressures_x_TM_tmp = zeros(size(t_press_x_tmp));
            selected_pressures_z_TM_tmp = zeros(size(t_press_x_tmp));
        
            for ii = 1:1:num_small_tables
                
                logical_specific_angles_indices_tmp = ...
                    logical(specific_angles_indices_tmp(:,ii));
                    
                if  ~isempty(specific_angles_tmp(logical_specific_angles_indices_tmp))
        
                    [indices_angles_tmp, ~] = knnsearch(pressure_TE_small_tables(1:last_indices_tables(ii),1:2,ii), ...
                        specific_angles_tmp(logical_specific_angles_indices_tmp,:));
        
                    selected_pressures_x_TE_tmp(logical(specific_angles_indices_tmp(:,ii))) = ...
                        c0*pressure_TE_small_tables(indices_angles_tmp,3,ii);
                    selected_pressures_z_TE_tmp(logical(specific_angles_indices_tmp(:,ii))) = ...
                        c0*pressure_TE_small_tables(indices_angles_tmp,4,ii);
                    selected_pressures_x_TM_tmp(logical(specific_angles_indices_tmp(:,ii))) = ...
                        c0*pressure_TM_small_tables(indices_angles_tmp,3,ii);
                    selected_pressures_z_TM_tmp(logical(specific_angles_indices_tmp(:,ii))) = ...
                        c0*pressure_TM_small_tables(indices_angles_tmp,4,ii);
        
                end
            end
        
            % Transform local frames of individual regions to global body frame.
            % Note that the lightsail is flipped upside down for this simulation
            % compared to those in COMSOL, thus requiring appropriate sign changes
            % for respective pressures
        
            t_press_x_tmp(indices_Region1) = selected_pressures_x_TE_tmp(indices_Region1);
            t_press_y_tmp(indices_Region1) = 0;
            t_press_z_tmp(indices_Region1) = selected_pressures_z_TE_tmp(indices_Region1);
        
            t_press_x_tmp(indices_Region2) = 0;
            t_press_y_tmp(indices_Region2) = selected_pressures_x_TM_tmp(indices_Region2);
            t_press_z_tmp(indices_Region2) = selected_pressures_z_TM_tmp(indices_Region2);
        
            t_press_x_tmp(indices_Region3) = -selected_pressures_x_TE_tmp(indices_Region3);
            t_press_y_tmp(indices_Region3) = 0;
            t_press_z_tmp(indices_Region3) = selected_pressures_z_TE_tmp(indices_Region3);
        
            t_press_x_tmp(indices_Region4) = 0;
            t_press_y_tmp(indices_Region4) = -selected_pressures_x_TM_tmp(indices_Region4);
            t_press_z_tmp(indices_Region4) = selected_pressures_z_TM_tmp(indices_Region4);
        
            % Calculate body-frame forces by discretizing the double integral into a
            % multiplication by the area of a mesh element (see Ramon's ppt
            % documentation for details), while accounting for the projected area
            % Note that the unit of the intensity [t_IMag] is W/(mm^2), while the unit
            % of the triangle areas [t_a] is mm^2, and the speed of light [c0] being
            % m/s, we get the correct unit for the forces, namely W*s/m = J/m = N
            % Note that cos(pitch) * cos(roll) accounts for the projected area
        
            t_force_x_BF_tmp = cos(pitch_angles_tmp) .* cos(roll_angles_tmp) .* t_a_tmp .* ...
                (t_IMag_tmp./c0) .* t_press_x_tmp;
            t_force_y_BF_tmp = cos(pitch_angles_tmp) .* cos(roll_angles_tmp) .* t_a_tmp .* ...
                (t_IMag_tmp./c0) .* t_press_y_tmp;
            t_force_z_BF_tmp = cos(pitch_angles_tmp) .* cos(roll_angles_tmp) .* t_a_tmp .* ...
                (t_IMag_tmp./c0) .* t_press_z_tmp;
        
            % Transform forces on triangle in local body frame to global laser
            % frame by multiplying body-frame force vector with the direction
            % cosine matrix H_I^B(apparently NOT the transpose, although I don't
            % understand?) according to OIlic2019
        
            t_force_x_tmp = t_force_x_BF_tmp .* DRM123_11(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp) + ...
                t_force_y_BF_tmp .* DRM123_21(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp) + ...
                t_force_z_BF_tmp .* DRM123_31(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp);
        
            t_force_y_tmp = t_force_x_BF_tmp .* DRM123_12(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp) + ...
                t_force_y_BF_tmp .* DRM123_22(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp) + ...
                t_force_z_BF_tmp .* DRM123_32(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp);
        
            t_force_z_tmp = t_force_x_BF_tmp .* DRM123_13(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp) + ...
                t_force_y_BF_tmp .* DRM123_23(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp) + ...
                t_force_z_BF_tmp .* DRM123_33(yaw_angles_tmp,pitch_angles_tmp,roll_angles_tmp);
        
            t_oth_tmp(1,:) =  t_force_x_tmp;
            t_oth_tmp(2,:) =  t_force_y_tmp;
            t_oth_tmp(3,:) =  t_force_z_tmp;
        
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
            n_of_tmp(1,:) = (t_oth_tmp(1,:) + t_oth_mr(1,:)) * M_t2n ./ 3;
            n_of_tmp(2,:) = (t_oth_tmp(2,:) + t_oth_mr(2,:)) * M_t2n ./ 3;
            n_of_tmp(3,:) = (t_oth_tmp(3,:) + t_oth_mr(3,:)) * M_t2n ./ 3;
        else
            n_of_tmp = (t_oth_tmp + t_oth_mr) * M_t2n ./ 3;
        end
    else
        if isbig  %this selects between two calculation methods, based on performance thresholds on my computer (MK).
            n_of_tmp(1,:) = t_oth_tmp(1,:) * M_t2n ./ 3;
            n_of_tmp(2,:) = t_oth_tmp(2,:) * M_t2n ./ 3;
            n_of_tmp(3,:) = t_oth_tmp(3,:) * M_t2n ./ 3;
        else
            n_of_tmp = t_oth_tmp * M_t2n ./ 3;
        end
    end


    % distribute absorbed power and effective radiator area to each node
    n_a_tmp = t_a_tmp * M_t2n ./ 3;
    if (RayTracingMode)
        n_abs_tmp = (t_abs_tmp+t_abs_mr) * M_t2n ./ 3;
    else
        n_abs_tmp = t_abs_tmp * M_t2n ./ 3;
    end

    % quick sanity check, if we've messed up a calulation, terminate the simulation immediately so we can debug
    % before the error corrupts the entire mesh.
    if (isnan(sum(n_abs_tmp)))
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

% n_ems_tmp = SBCmm .* n_a .* Emissivity .* (n_t_tmp.^4); % Watt
% n_ems_tmp = SBCmm .* n_a_tmp .* Emissivity_Si3N4 .* (n_t_tmp.^4); % Watt
n_ems_tmp = SBCmm .* n_a_tmp .* 0*Emissivity_Si3N4 .* (n_t_tmp.^4); % Watt


%% Step 7: Calculate mechanical forces & heat conduction

% Matrix containing all edge direction vectors with x, y, z components as columns
e_nl_tmp = [ n_x_tmp(e_nb) - n_x_tmp(e_na) ; n_y_tmp(e_nb) - n_y_tmp(e_na) ; n_z_tmp(e_nb) - n_z_tmp(e_na) ];

% Calculate current lengths in mm of all edges
e_l_tmp = vecnorm( e_nl_tmp);

% Normalize all edge direction vectors. Here, 'e_nl' stands for edge, normalized length
e_nl_tmp = e_nl_tmp ./ e_l_tmp;   % unit vectors

% Calculate linear expansion due to temperature difference e_t minus initial temperature t0 for each edge 
% multiplied by the coefficient of thermal expansion.
% Amount of thermal expansion can be described by material strain
% epsilon_thermal = (L_final - L_initial)/L_initial, where the strain
% is due change in temperature proportional to the coefficient of
% thermal expansion, epsilon_thermal = -alpha_L * (T_final - T_initial)
% This can be rewritten in terms of length change dL / L_initial =
% 1 - (1 + alpha_L * (T_final - T_initial) )

e_tex_tmp = 1 + (e_t_tmp - t0).*(CTE*1e-6);  % edge length ratio due to thermal expansion
e_l0t_tmp = e_l0 .* e_tex_tmp;  % resting (unstrained) edge length, including thermal expansion
e_dl_tmp = e_l_tmp - e_l0t_tmp ;    % mm, current edge length difference vs. unstrained edge length, positive for tension

% Calculate strain based on current vs. resting/unstrained length
e_s_tmp = e_l_tmp ./ (e_l0t_tmp) - 1;   % strain %/100.   No stress -> 0 strain. Elongate to 2x length -> 1.0 strain

% Calculate mechanical force on edge as the the edge length difference
% times the edge stiffness/spring constant being equal to the 2D
% Young's modulus x area / length^2 (see generateMesh.m for more info).
% The sign of the mechanical force on a specific edge is determined by
% the corresponding edge direction vector 'e_nl'
e_mf_tmp = (e_notbroken & e_flexible) .* (e_k) .* e_dl_tmp .* e_nl_tmp;  % Newtons
% Note that we zero the mechanical force if the edge is fully within a rigid body.  

% calculate potential energy due to edge strain
e_mf_mag2_tmp = e_k .* e_dl_tmp.^2 .* (e_notbroken);  
PE_tmp = 5e-4 * sum(e_mf_mag2_tmp);  % should be in Joules, I think...

% Calculate temperature difference across all edges due to different temperatures on every node
e_dt_tmp = n_t_tmp(e_nb) - n_t_tmp(e_na);   % Kelvin

% Thermal conduction, i.e. heat flow power on edges ('e_hf') according to formula kappa * (A/l) * (T_hot - T_cold)
e_hf_tmp = e_heatCondc .* e_dt_tmp;  % Watts

% Terminate script if we've picked up a math error somehow, so we can debug before it corrupts the mesh.
if isnan(sum(e_hf_tmp))
    error('Found NaN!')
end

% Keep track of time spent doing edge-based calculations
timers(5) = timers(5) + toc;
tic;


%% Step 8: Assign mechanical forces to nodes (from edges) and also heat flow
n_mf_tmp(1,:) = e_mf_tmp(1,:) * M_e2n;
n_mf_tmp(2,:) = e_mf_tmp(2,:) * M_e2n;
n_mf_tmp(3,:) = e_mf_tmp(3,:) * M_e2n;

% Calculate node heat flow from edge heat flow
n_hf_tmp = e_hf_tmp * M_e2n;

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

    e_tensilefailure =  ( e_flexible &  e_notbroken) & ( e_s_tmp > e_tensile_strain_limit )  ;  %all edges that are now broken from tensile failure
    e_thermalfailure = ( e_notbroken .* (e_t_tmp > e_failtemp ) );
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
        time_frstBrkn = tt_tmp;
        
        if terminate_on_breakup
            if ~terminate_t
                terminate_t = tt_tmp + termination_delay_s;
                if ~terminate_t
                    terminate_t = dt_tmp;  % just in case the first breakage would try to schedule termination at exactly tt=0, we'll change it to a non-zero value, because zero is a sentinel value that prevents termination...
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
n_rfx_tmp = n_of_tmp(1,:) + n_af_tmp(1,:);
n_rfy_tmp = n_of_tmp(2,:) + n_af_tmp(2,:);
n_rfz_tmp = n_of_tmp(3,:) + n_af_tmp(3,:);

externalForceX_tmp = sum( n_rfx_tmp );
externalForceY_tmp = sum( n_rfy_tmp );
externalForceZ_tmp = sum( n_rfz_tmp );

externalForce_tmp = [ externalForceX_tmp; externalForceY_tmp; externalForceZ_tmp];

% We can calculate an effective torque on the whole mesh (as if it were rigid).
% This is not actually used in subsequent calculations...
externalTorqueX_tmp = sum( n_dy.*n_rfz_tmp - n_dz.*n_rfy_tmp );
externalTorqueY_tmp = sum( n_dz.*n_rfx_tmp - n_dx.*n_rfz_tmp );
externalTorqueZ_tmp = sum( n_dx.*n_rfy_tmp - n_dy.*n_rfx_tmp );

externalTorque_tmp = [ externalTorqueX_tmp; externalTorqueY_tmp; externalTorqueZ_tmp ];

% Sum mechanical and optical forces on each node
n_fx_tmp = n_mf_tmp(1,:) + n_rfx_tmp;
n_fy_tmp = n_mf_tmp(2,:) + n_rfy_tmp;
n_fz_tmp = n_mf_tmp(3,:) + n_rfz_tmp;
% Note that n_mf_tmp doesn't include edge spring contributions from within rigid bodies


%% Step 11: Update temperatures
% Heat capacity links the temperature rise of a body to the required
% amount of energy. The unit of heat capacity is [J/g/degC], which
% describes the amount of energy required to raise the temperature of 1
% gram of an object by 1 deg C (or 1 Kelvin). Energy is equal to power
% times time, such that the temperature rise is given by energy / mass
% / heat capacity, where energy = (input power - output power) * dt.
% Input power is given by the opticaly absorbed power and power from
% heat flow/thermal conduction, whereas output power is due to
% radiative cooling.

% MOVED TO MAIN SCRIPT FOR RK4

% n_t_tmp = n_t_tmp + (n_hf_tmp + n_abs_tmp - n_ems_tmp).*dt_tmp ./ n_heatmass;   % todo:  need to support variable heat cap across the array
% 
% % Calculate the edge temperature as the averaged temperature of the sum
% % of the temperatures of the two nodes that form the edge
% e_t_tmp = (n_t_tmp(e_na) + n_t_tmp(e_nb)) ./ 2;
% 
% % Calculate temperature of triangle as the average of the temperatures
% % of the nodes that make up the triangle
% t_t_tmp  =  ( n_t_tmp(t_na) + n_t_tmp(t_nb) + n_t_tmp(t_nc) ) ./ 3;


%% Step 12: Update triangle centroids

% MOVED TO MAIN SCRIPT FOR RK4

% t_cx_tmp =  ( n_x_tmp(t_na) + n_x_tmp(t_nb) + n_x_tmp(t_nc) ) / 3;
% t_cy_tmp =  ( n_y_tmp(t_na) + n_y_tmp(t_nb) + n_y_tmp(t_nc) ) / 3;
% t_cz_tmp =  ( n_z_tmp(t_na) + n_z_tmp(t_nb) + n_z_tmp(t_nc) ) / 3;
    
