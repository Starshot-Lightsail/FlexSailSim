% tilt-displacement stability visualizer v17 - 05/04/2021 (c) Michael Kelzenberg, California Institute of Technology
disp('Tilt Stability Visualizer Script, V20 (Jul 2023), Michael Kelzenberg, California Institute of Technology');

isRotationallySymmetric=1;

isam_idir_xa = 20*deg2rad( -1:.04:1  ); %1.  Sweep x incidence direction
isam_idir_ya = 20*deg2rad( -1:.04:1  ); %2.  Sweep y incidence direction
isam_disp_xo = 1.5*radActualx0mm .* (-1:.04:1); % 3.  Sweep x beam center
isam_disp_yo = 1.5*radActualx0mm .* (-1:.04:1); % 3.  Sweep y beam center

isamx_forces = ones( length(isam_disp_xo), length(isam_idir_ya));
isamy_forces = ones( length(isam_disp_yo), length(isam_idir_xa));
isamy_torques = zeros( length(isam_disp_xo), length(isam_idir_ya));
isamx_torques = zeros( length(isam_disp_yo), length(isam_idir_xa));

DRM123_11 = @(psi,theta,phi) cos(psi).*cos(theta);
DRM123_12 = @(psi,theta,phi) cos(psi)*sin(theta).*sin(phi) + sin(psi).*cos(phi);
DRM123_13 = @(psi,theta,phi) -cos(psi).*sin(theta).*cos(phi) + sin(psi).*sin(phi);
DRM123_21 = @(psi,theta,phi) -sin(psi).*cos(theta);
DRM123_22 = @(psi,theta,phi) -sin(psi).*sin(theta).*sin(phi) + cos(psi).*cos(phi);
DRM123_23 = @(psi,theta,phi) sin(psi)*sin(theta).*cos(phi) + cos(psi).*sin(phi);
DRM123_31 = @(psi,theta,phi) sin(theta);
DRM123_32 = @(psi,theta,phi) -cos(theta).*sin(phi);
DRM123_33 = @(psi,theta,phi) cos(theta).*cos(phi);

disp('Starting 2D stability analysis...');
disp(['Stability sweep 1:  X offset from ' num2str((isam_disp_xo(1))) ' to ' num2str((isam_disp_xo(end))) ' mm x ' num2str(length(isam_disp_xo)) '.']);
disp(['                    Y tilt from '  num2str(rad2deg(isam_idir_ya(1))) ' to ' num2str(rad2deg(isam_idir_ya(end))) ' deg x ' num2str(length(isam_idir_ya)) '.']);



saveStabVideo = 0;

if saveStabVideo
    figure('Position',[625         500        2048        1024]);
    VS = VideoWriter([ filebasename '_tiltDisp2' ],'MPEG-4');  %alternate format for Ramon:  'Archival'
    VS.Quality = 100;
    open(VS);
end

climmax = 0;

startisamtic = tic;
updatetic = tic;
tic;
nsamnt = length(isam_disp_xo) * length(isam_idir_ya);
nsamn = 0;

for nsamx = 1:length(isam_disp_xo)
    
    disp_xo = isam_disp_xo(nsamx);
    
    for nsamy = 1:length(isam_idir_ya)
        nsamn = nsamn+1;
        
        tilt_ya = isam_idir_ya(nsamy);
        
        phi_ex = 0;
        theta_ex = tilt_ya;
        psi_ex = 0;
        
        DRM_mat = [cos(psi_ex).*cos(theta_ex), cos(psi_ex)*sin(theta_ex).*sin(phi_ex) + sin(psi_ex).*cos(phi_ex), -cos(psi_ex).*sin(theta_ex).*cos(phi_ex) + sin(psi_ex).*sin(phi_ex); ...
            -sin(psi_ex).*cos(theta_ex), -sin(psi_ex).*sin(theta_ex).*sin(phi_ex) + cos(psi_ex).*cos(phi_ex), sin(psi_ex)*sin(theta_ex).*cos(phi_ex) + cos(psi_ex).*sin(phi_ex); ...
            sin(theta_ex), -cos(theta_ex).*sin(phi_ex), cos(theta_ex).*cos(phi_ex)];
        
        n_x = n_x0 .* DRM123_11(psi_ex,theta_ex,phi_ex) + ...
            n_y0 .* DRM123_12(psi_ex,theta_ex,phi_ex) + ...
            n_z0 .* DRM123_13(psi_ex,theta_ex,phi_ex);
        
        n_y = n_x0 .* DRM123_21(psi_ex,theta_ex,phi_ex) + ...
            n_y0 .* DRM123_22(psi_ex,theta_ex,phi_ex) + ...
            n_z0 .* DRM123_23(psi_ex,theta_ex,phi_ex);
        
        n_z = n_x0 .* DRM123_31(psi_ex,theta_ex,phi_ex) + ...
            n_y0 .* DRM123_32(psi_ex,theta_ex,phi_ex) + ...
            n_z0 .* DRM123_33(psi_ex,theta_ex,phi_ex);
        
        n_x = n_x + disp_xo;
        
        
        
        tiltStabilityAnalysisSub;
        
        if ~mod(nsamy,25) && ~mod(nsamx,5)
           if saveStabVideo
            plotMeshSAM;
            title({['X=' num2str(disp_xo) '  Ytilt = ' num2str(rad2deg(tilt_ya)) '\circ  Fx=' num2str(sum(n_of(1,:))) '  Fy=' num2str(sum(n_of(2,:))) '  Fz=' num2str(sum(n_of(3,:)))]}); 
            writeVideo(VS,getframe(gcf));
            %pause(0.01);
           end
           
        end
        
        if ~mod(nsamy,10)
            if toc(updatetic)>5
               disp([ num2str(100*nsamn/nsamnt) '%...']); updatetic = tic;
            end
        end
        
        %n_rfx = n_of(1,:) ;
        %n_rfy = n_of(2,:) ;
        %n_rfz = n_of(3,:) ;
        %rigidForceX = sum( n_of(1,:) );
        isamx_forces(nsamx, nsamy) = sum( n_of(1,:) );
        isamy_torques(nsamx, nsamy) = sum( n_dz.*n_of(1,:) - n_dx.*n_of(3,:) );
        %upsidedowntris = sum( t_norms(3,:) < 0);
        if upsidedowntris
            isamx_forces(nsamx, nsamy) = 0;
            isamy_torques(nsamx, nsamy) = 0;
        end
        %rigidForceY = sum( n_of(2,:) );
        %rigidForceZ = sum( n_of(3,:) );
        %rigidForce = [ rigidForceX; rigidForceY; rigidForceZ];
        
        %rigidTorqueX = sum( n_dy.*n_rfz - n_dz.*n_rfy );
        %rigidTorqueY = sum( n_dz.*n_rfx - n_dx.*n_rfz );
        %rigidTorqueZ = sum( n_dx.*n_rfy - n_dy.*n_rfx );
        %rigidTorque = [ rigidTorqueX; rigidTorqueY; rigidTorqueZ ];
        
        % % debug, make sure I got the axis labels right in the plot, mark a known area on the data set
        % if disp_xo/radActualx0mm >= .5  && disp_xo/radActualx0mm <= .9
        %     if rad2deg(tilt_ya) >= -3  && rad2deg(tilt_ya) <= -2
        %         isamx_forces(nsamx, nsamy) = NaN;
        %         isamy_torques(nsamx, nsamy) = NaN;
        %     end
        % end
    end
end

if ~isRotationallySymmetric
    disp(['Stability sweep 2:  Y offset from ' num2str((isam_disp_yo(1))) ' to ' num2str((isam_disp_yo(end))) ' mm x ' num2str(length(isam_disp_yo)) '.']);
    disp(['                    X tilt from '  num2str(rad2deg(isam_idir_xa(1))) ' to ' num2str(rad2deg(isam_idir_xa(end))) ' deg x ' num2str(length(isam_idir_xa)) '.']);
    tic;
    
    nsamnt = length(isam_disp_yo) * length(isam_idir_xa);
    nsamn = 0;
    
    for nsamy = 1:length(isam_disp_yo)
        
        disp_yo = isam_disp_yo(nsamy);
        
        for nsamx = 1:length(isam_idir_xa)
            nsamn = nsamn+1;
            tilt_xa = isam_idir_xa(nsamx);
            
             phi_ex = tilt_xa;
            theta_ex = 0;
            psi_ex = 0;
            
            DRM_mat = [cos(psi_ex).*cos(theta_ex), cos(psi_ex)*sin(theta_ex).*sin(phi_ex) + sin(psi_ex).*cos(phi_ex), -cos(psi_ex).*sin(theta_ex).*cos(phi_ex) + sin(psi_ex).*sin(phi_ex); ...
                -sin(psi_ex).*cos(theta_ex), -sin(psi_ex).*sin(theta_ex).*sin(phi_ex) + cos(psi_ex).*cos(phi_ex), sin(psi_ex)*sin(theta_ex).*cos(phi_ex) + cos(psi_ex).*sin(phi_ex); ...
                sin(theta_ex), -cos(theta_ex).*sin(phi_ex), cos(theta_ex).*cos(phi_ex)];
            
            n_x = n_x0 .* DRM123_11(psi_ex,theta_ex,phi_ex) + ...
                n_y0 .* DRM123_12(psi_ex,theta_ex,phi_ex) + ...
                n_z0 .* DRM123_13(psi_ex,theta_ex,phi_ex);
            
            n_y = n_x0 .* DRM123_21(psi_ex,theta_ex,phi_ex) + ...
                n_y0 .* DRM123_22(psi_ex,theta_ex,phi_ex) + ...
                n_z0 .* DRM123_23(psi_ex,theta_ex,phi_ex);
            
            n_z = n_x0 .* DRM123_31(psi_ex,theta_ex,phi_ex) + ...
                n_y0 .* DRM123_32(psi_ex,theta_ex,phi_ex) + ...
                n_z0 .* DRM123_33(psi_ex,theta_ex,phi_ex);
            
            n_y = n_y + disp_yo;
            
            tiltStabilityAnalysisSub;
            
            if ~mod(nsamx,25) && ~mod(nsamy,5)
                if saveStabVideo
                    plotMeshSAM;
                    title({['Y=' num2str(disp_yo) '  Xtilt = ' num2str(rad2deg(tilt_xa)) '\circ  Fy=' num2str(sum(n_of(2,:))) '  Fx=' num2str(sum(n_of(1,:))) '  Fz=' num2str(sum(n_of(3,:)))]}); 
                    writeVideo(VS,getframe(gcf));
                    %pause(0.05);
                end
    
            end
            
            if ~mod(nsamx,10)
                if toc(updatetic)>5
                   disp([ num2str(100*nsamn/nsamnt) '%...']); updatetic = tic;
                end
            end
    
            
            
            %n_rfx = n_of(1,:) ;
            %n_rfy = n_of(2,:) ;
            %n_rfz = n_of(3,:) ;
            %rigidForceX = sum( n_of(1,:) );
            isamy_forces(nsamy, nsamx) = sum( n_of(2,:) );
            isamx_torques(nsamy,nsamx) = sum( n_dy.*n_of(3,:) - n_dz.*n_of(2,:) );
            if upsidedowntris
                isamy_forces(nsamy, nsamx) = 0;
                isamx_torques(nsamy,nsamx) = 0;
            end
            %rigidForceY = sum( n_of(2,:) );
            %rigidForceZ = sum( n_of(3,:) );
            %rigidForce = [ rigidForceX; rigidForceY; rigidForceZ];
            
            %rigidTorqueX = sum( n_dy.*n_rfz - n_dz.*n_rfy );
            %rigidTorqueY = sum( n_dz.*n_rfx - n_dx.*n_rfz );
            %rigidTorqueZ = sum( n_dx.*n_rfy - n_dy.*n_rfx );
            %rigidTorque = [ rigidTorqueX; rigidTorqueY; rigidTorqueZ ];
            
            % % debug, make sure I got the axis labels right in the plot, mark a known area on the data set
            % if disp_yo/radActualx0mm >= -1.2  && disp_yo/radActualx0mm <= -.9
            %     if rad2deg(tilt_xa) >= 0  && rad2deg(tilt_xa) <= 4.5
            %         isamy_forces(nsamy, nsamx) = NaN;
            %         isamx_torques(nsamy, nsamx) = NaN;
            %     end
            % end
        end
    end
end

%Ok, return the mesh to where it was before we started this nonsense!
n_x = n_x0;
n_y = n_y0;
n_z = n_z0;
coms_x = n_x .* n_m;
coms_y = n_y .* n_m;
coms_z = n_z .* n_m;
comx = sum(coms_x)./totalmass_nodes ;
comy = sum(coms_y)./totalmass_nodes ;
comz = sum(coms_z)./totalmass_nodes ;
n_dx = n_x - comx;
n_dy = n_y - comy;
n_dz = n_z - comz;

if saveStabVideo
    close(VS);
    close(gcf);
end

if (~exist('forceColormap6','var'))
        load ForceColormap.mat;
end

if isRotationallySymmetric
    figure('Position',[46        600        800         1200], 'Color', 'w');
    subplot(2,1,1);
else
    figure('Position',[46        600        1400         1200], 'Color', 'w');
    subplot(2,2,1)
end

imagesc(([isam_disp_xo(1) isam_disp_xo(end)])./radActualx0mm, rad2deg([ isam_idir_ya(1) isam_idir_ya(end)]), isamx_forces', 'Interpolation', 'bilinear' );
    title(['Restoring F_X (R_{beam}/R_{sail} = ' num2str(Irad/sqrt(2) /radActualx0mm) ')'] );
    xlabel('X_{offset} / R_{sail}');
    ylabel('Y tilt (\circ)');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap6);
    hb = colorbar;
    hb.Label.String = ['Force (N)  [I_0=' num2str(I0/1e3) ' GW/m^2]'];
    set(hb,'FontSize',14);
    set(gca,'FontSize',14);

if isRotationallySymmetric
    subplot(2,1,2);
else
    subplot(2,2,3);
end

imagesc(rad2deg([ isam_idir_ya(1) isam_idir_ya(end)]), ([isam_disp_xo(1) isam_disp_xo(end)])./radActualx0mm, isamy_torques./1000, 'Interpolation', 'bilinear' );
    title(['Restoring \tau_Y  (R_{beam}/R_{sail} = ' num2str(Irad/sqrt(2) /radActualx0mm) ')'] );
    ylabel('X_{offset} / R_{sail}');
    xlabel('Y tilt (\circ)');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap6);
    hb = colorbar;
    hb.Label.String = ['Torque (N\cdotm)  [I_0 = ' num2str(I0/1e3) ' GW/m^2]'];
    set(hb,'FontSize',14);
    set(gca,'FontSize',14);


    
if ~isRotationallySymmetric    
    subplot(2,2,2)
    imagesc(([isam_disp_yo(1) isam_disp_yo(end)])./radActualx0mm, rad2deg([ isam_idir_xa(1) isam_idir_xa(end)]), isamy_forces', 'Interpolation', 'bilinear' );
    title(['Restoring F_Y (R_{beam}/R_{sail} = ' num2str(Irad/sqrt(2) /radActualx0mm) ')'] );
    %title(['Restoring $$Y$$ force       ($$\frac{R_{beam}}{R_{sail}} = ' num2str(Irad/sqrt(2) /radActualx0mm) ', I_0 = ' num2str(I0/1e3) '$$ GW/m$^2$)'])
    xlabel('Y_{offset} / R_{sail}');
    ylabel('X tilt (\circ)');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap6);
    hb = colorbar;
    hb.Label.String = ['Force (N)  [I_0=' num2str(I0/1e3) ' GW/m^2]'];
    set(hb,'FontSize',14);
    set(gca,'FontSize',14);

    subplot(2,2,4)
    imagesc( rad2deg([ isam_idir_xa(1) isam_idir_xa(end)]), ([isam_disp_yo(1) isam_disp_yo(end)])./radActualx0mm, isamx_torques./1000, 'Interpolation', 'bilinear' );
    title(['Restoring \tau_X (R_{beam}/R_{sail} = ' num2str(Irad/sqrt(2) /radActualx0mm) ')'] );
    %title(['Restoring $$X$$ torque       ($$\frac{R_{beam}}{R_{sail}} = ' num2str(Irad/sqrt(2) /radActualx0mm) ', I_0 = ' num2str(I0/1e3) '$$ GW/m$^2$)'])
    xlabel('X tilt (\circ)');
    ylabel('Y_{offset} / R_{sail}');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap6);
    hb = colorbar;
    hb.Label.String = ['Torque (N\cdotm)  [I_0=' num2str(I0/1e3) ' GW/m^2]'];
    set(hb,'FontSize',14);
    set(gca,'FontSize',14);
end

saveas(gcf,[filebasename '_latStab.fig']);
saveas(gcf,[filebasename '_latStab.png']);
                
   
                