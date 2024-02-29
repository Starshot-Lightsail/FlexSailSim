%tiltStabilityAnalysisSub

%% part of V2
 t_cx =  ( n_x(t_na) + n_x(t_nb) + n_x(t_nc) ) / 3;
                t_cy =  ( n_y(t_na) + n_y(t_nb) + n_y(t_nc) ) / 3;
                t_cz =  ( n_z(t_na) + n_z(t_nb) + n_z(t_nc) ) / 3;
t_ev1 = [n_x; n_y; n_z]* M_nxyz2tev1;
t_ev2 = [n_x; n_y; n_z]* M_nxyz2tev2;
t_norms = cross( t_ev1, t_ev2 );
t_a = max(vecnorm(t_norms), 1e-8);
t_norms = t_norms ./ t_a ; %ok now t_norms contains proper normal vectors (normalized)
t_a = 0.5 .* t_a;
t_a = (t_notbroken).*t_a + t_broken.*t_a0;
t_ev1n = t_ev1./(max(vecnorm(t_ev1),1e-8));
t_texn = t_ev1n .* (t_costexa) + ...
    cross(t_norms, t_ev1n ) .* (t_sintexa)  + ...
    t_norms .* dot(t_norms, t_ev1n ) .* (1 - cos(t_texa)) ;
t_textn = cross(t_norms, t_texn);  % this is the transverse unit vector

       t_optFcn_r2 = ( (t_cx) ).^2 + (t_cy ).^2;
       t_IMag = I0 * exp( - t_optFcn_r2 ./ (Irad^2) );

%% part of V1
        t_costheta = dot(t_norms, t_IDir);
        t_aproj = t_a .* t_costheta .* (t_notbroken);
        if useSpecularTMM
            t_theta = acos(t_costheta);
            if min(t_theta) < 0
                warning('Negative theta in useSpecularTMM code');
                t_theta(t_theta < 0) = 0;
            end
            if max(t_theta) > pi/2
                %    warning('Reverse illuminated triangles in specularTMM code');
                t_theta(t_theta > pi/2) = pi/2 - t_theta(t_theta > pi/2);
            end
            t_LULi = max(min(floor(t_theta ./ LULAngleStep)+1,LUL.num),1);
            
            t_specA = LUL.A(t_LULi);
            t_specR = LUL.R(t_LULi);
        else
            t_specA(:) = Iabs;
            t_specR(:) = Irefl;
        end
        t_Ipwr = abs(t_aproj.*t_IMag);
        t_abs = t_specA .* t_Ipwr;
        t_reflp = t_Ipwr .* t_specR;
        t_oth = 2 * t_reflp ./ c0 .* t_norms .* t_costheta  + ... %  this is the reflection force
            t_abs       ./ c0 .* t_IDir ;    % and this is the absorption force
        
        if RayTracingMode
            t_refln = t_IDir - 2 .* t_costheta .* t_norms;
                
            processMR_GPU_v17;  % call raytracing subscript.  Generates:
                  %  t_oth_mr :  additional optical thrust from multiple reflections
                  %  t_abs_mr :  additional absorption from multiple reflections.
                  
                if isbig  %this selects between to calculation methods, based on performance thresholds on my computer (MK).
                    n_of(1,:) = (t_oth(1,:) + t_oth_mr(1,:)) * M_t2n ./ 3;
                    n_of(2,:) = (t_oth(2,:) + t_oth_mr(2,:)) * M_t2n ./ 3;
                    n_of(3,:) = (t_oth(3,:) + t_oth_mr(3,:)) * M_t2n ./ 3;
                else
                    n_of = (t_oth + t_oth_mr) * M_t2n ./ 3;
                end
        
        else
            %% LUT stuff for MEPS
            t_incsignt = dot( t_IDir, t_textn );  %projection mag along transverse vector
            t_incprojp = t_IDir - t_textn .* t_incsignt ; % inc vector projected into parallel plane
            t_incprojpn = t_incprojp ./ vecnorm(t_incprojp) ; % now normalized

            t_incsignp = dot( t_IDir, t_texn );  % projection mag along parallel vector
            t_incprojt = t_IDir - t_texn .* t_incsignp ; % inc vector projected into transverse plane
            t_incprojtn = t_incprojt ./ vecnorm(t_incprojt) ; % now normalized

            t_incsignt = sign(t_incsignt);  % get sign for angle unwrap from acos
            t_incsignp = sign(t_incsignp);  % ...

            t_psi = t_incsignp .* acos( max(min(dot( t_incprojpn, t_norms ), 1),-1) );  % the min/max confine us to real limits for acos().
            t_theta = t_incsignt .* acos( max(min(dot( t_incprojtn, t_norms ), 1),-1) );  % they shouldn't be needed, but I'm not sure.

            t_ipsi = round((t_psi+pi/2)./deg2rad(LUTPsiStepDeg)); % Psi index for LUTs
            t_itheta = round((t_theta+pi/2)./deg2rad(LUTThetaStepDeg)); % Theta index for LUTs

            % debug
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
                if (outofboundts ~= lastoutofboundtswarn)
                    lastoutofboundtswarn = outofboundts;

                    error(['There are ' num2str(outofboundts) ' reverse-illuminated triangles at nt=' num2str(nt)  '!  ' ...
                        '(psi_low = ' num2str(noobpl) ',  psi_high = ' num2str(noobph) ...
                        ',  theta_low = ' num2str(noobtl) ',  theta_high = ' num2str(noobth) ')' ]);
                end
                if noobpl t_ipsi(outofboundspsilow) = 1; end
                if noobph t_ipsi(outofboundspsihigh) = 180/LUTPsiStepDeg+1; end
                if noobtl t_itheta(outofboundsthetalow) = 1; end
                if noobth t_itheta(outofboundsthetahigh) = 180/LUTThetaStepDeg+1; end
            end


            for nlut = 1:length(LUTs)
                lutinds = sub2ind(size(LUTs(nlut).n), t_ipsi(t_tex == nlut),t_itheta(t_tex == nlut));
                t_press_n(t_tex == nlut) = LUTs(nlut).n(lutinds);
                t_press_t(t_tex == nlut) = LUTs(nlut).t(lutinds);
                t_press_p(t_tex == nlut) = LUTs(nlut).p(lutinds);
                t_lut_a( t_tex == nlut ) = LUTs(nlut).a(lutinds);
            end

            lutmask = t_tex>0;
            t_oth(:,lutmask) = t_Ipwr(lutmask) .*  ( ...
                t_norms(:,lutmask) .* t_press_n(lutmask)  + ... this is the normal force
                t_textn(:,lutmask) .* t_press_t(lutmask)  + ... this is the transverse force
                t_texn(:,lutmask) .* t_press_p(lutmask)); % this is the parallel/axial force

            t_abs(lutmask) = t_lut_a(lutmask) .* t_Ipwr(lutmask);
            
            if isbig  %this selects between to calculation methods, based on performance thresholds on my computer (MK).
                n_of(1,:) = t_oth(1,:) * M_t2n ./ 3;
                n_of(2,:) = t_oth(2,:) * M_t2n ./ 3;
                n_of(3,:) = t_oth(3,:) * M_t2n ./ 3;
            else
                n_of = t_oth * M_t2n ./ 3;
            end
        end