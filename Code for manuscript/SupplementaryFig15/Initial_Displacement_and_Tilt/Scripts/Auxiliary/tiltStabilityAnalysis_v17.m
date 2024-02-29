% tilt-displacement stability visualizer v17 - 05/04/2021 (c) Michael Kelzenberg, California Institute of Technology
disp('Tilt Stability Visualizer Script, V17 (May 4 2021), Michael Kelzenberg, California Institute of Technology');

%set up simulation environment:
filebasename = 'tiltdisplacementsim';
setupMaterialProperties_v17;
%setupLUTs_v17;
setupSpecularLULs_v17;

generateMesh_v17;

MAX_REFL = 6;
useSpecularTMM = 1;
RayTracingMode = 1;
illum_mode = 'gaussian';
I0 = 1000;
Irad=0.5*radiusmm*sqrt(2);

isam_idir_xa = 20*deg2rad( -1:.05:1  ); %1.  Sweep x incidence direction
isam_idir_ya = 20*deg2rad( -1:.05:1  ); %2.  Sweep y incidence direction
isam_disp_xo = 3*radActualx0mm .* [-1:.05:1]; % 3.  Sweep x beam center
isam_disp_yo = 3*radActualx0mm .* [-1:.05:1]; % 3.  Sweep y beam center

isamx_forces = ones( length(isam_disp_xo), length(isam_idir_ya));
isamy_forces = ones( length(isam_disp_yo), length(isam_idir_xa));

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


if ~isequal(illum_mode,'gaussian')
    error(['Attempted to use non-defined illumination mode ''' illum_mode '''']);
end

savevideo = 0;

if savevideo
    figure('Position',[625         500        2048        1024]);
    V1 = VideoWriter([ 'tiltDisp_video2' ],'MPEG-4');  %alternate format for Ramon:  'Archival'
    V1.Quality = 100;
    open(V1);
end

t_specA = zeros(size(t_na));  % specular absorption
t_specR = zeros(size(t_na));  % specular reflection
t_IDir = zeros(3, length(t_na));
t_IDir(3,:) = 1;
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
        
        
       % t_optFcn_r2 = ( (t_cx - disp_xo) / cos(tilt_ya) ).^2 + (t_cy - 0).^2;
        % t_IMag = I0 * exp( - t_optFcn_r2 ./ (Irad^2) );
        
        %t_IDir(1,:) = -sin(tilt_ya);
        %t_IDir(3,:) = cos(tilt_ya);
        
        tiltStabilityAnalysisSub;
        
        if ~mod(nsamy,5) && ~mod(nsamx,25)
           if savevideo
            plotMeshSAM;
            title({['X=' num2str(disp_xo) '  Ytilt = ' num2str(rad2deg(tilt_ya)) '\circ  Fx=' num2str(sum(n_of(1,:))) '  Fy=' num2str(sum(n_of(2,:))) '  Fz=' num2str(sum(n_of(3,:)))]}); 
            writeVideo(V1,getframe(gcf));
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
        %rigidForceY = sum( n_of(2,:) );
        %rigidForceZ = sum( n_of(3,:) );
        %rigidForce = [ rigidForceX; rigidForceY; rigidForceZ];
        
        %rigidTorqueX = sum( n_dy.*n_rfz - n_dz.*n_rfy );
        %rigidTorqueY = sum( n_dz.*n_rfx - n_dx.*n_rfz );
        %rigidTorqueZ = sum( n_dx.*n_rfy - n_dy.*n_rfx );
        %rigidTorque = [ rigidTorqueX; rigidTorqueY; rigidTorqueZ ];
        
        % debug, make sure I got the axis labels right in the plot, mark a known area on the data set
%         if disp_xo/radActualx0mm >= .3  && disp_xo/radActualx0mm <= .5
%             if rad2deg(tilt_ya) >= -2.5  && rad2deg(tilt_ya) <= -2
%                 isamx_forces(nsamx, nsamy) = 0;
%             end
%         end
    end
end

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
        
         phi_ex = -tilt_xa;
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
        
        if ~mod(nsamx,5) && ~mod(nsamy,25)
            if savevideo
                plotMeshSAM;
                title({['Y=' num2str(disp_yo) '  Xtilt = ' num2str(rad2deg(tilt_xa)) '\circ  Fy=' num2str(sum(n_of(2,:))) '  Fx=' num2str(sum(n_of(1,:))) '  Fz=' num2str(sum(n_of(3,:)))]}); 
                writeVideo(V1,getframe(gcf));
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
        %rigidForceY = sum( n_of(2,:) );
        %rigidForceZ = sum( n_of(3,:) );
        %rigidForce = [ rigidForceX; rigidForceY; rigidForceZ];
        
        %rigidTorqueX = sum( n_dy.*n_rfz - n_dz.*n_rfy );
        %rigidTorqueY = sum( n_dz.*n_rfx - n_dx.*n_rfz );
        %rigidTorqueZ = sum( n_dx.*n_rfy - n_dy.*n_rfx );
        %rigidTorque = [ rigidTorqueX; rigidTorqueY; rigidTorqueZ ];
        
        % debug, make sure I got the axis labels right in the plot, mark a known area on the data set
%         if disp_yo/radActualx0mm >= -.5  && disp_yo/radActualx0mm <= -.3
%             if rad2deg(tilt_xa) >= -4.2  && rad2deg(tilt_xa) <= -4
%                 isamy_forces(nsamy, nsamx) = 0;
%             end
%         end
    end
end

n_x = n_x0;
n_y = n_y0;
n_z = n_z0;

if savevideo
    close(V1);
    close(gcf);
end

if (~exist('forceColormap5'))
        load ForceColormap_v17.mat;
end

figure('Position',[46        1105        1470         656]);
subplot(1,2,1)
imagesc(([isam_disp_xo(1) isam_disp_xo(end)])./radActualx0mm, rad2deg([ isam_idir_ya(1) isam_idir_ya(end)]), isamx_forces' );
    title(['Restoring force along X (Irad=' num2str(Irad /radActualx0mm) ')'] );
    xlabel('X displacement, normalized to radius');
    ylabel('Y tilt (\circ)');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap5)
    colorbar;
    
    subplot(1,2,2)
imagesc(([isam_disp_yo(1) isam_disp_yo(end)])./radActualx0mm, rad2deg([ isam_idir_xa(1) isam_idir_xa(end)]), isamy_forces' );
    title('Restoring force along Y');
    xlabel('Y displacement, normalized to radius');
    ylabel('X tilt (\circ)');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap5)
    colorbar;
                
   
                