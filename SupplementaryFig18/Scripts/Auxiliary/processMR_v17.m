%% processMR_v17  -- process multi-reflection forces for non-flat lightsails
%
% Applies iterative approximate raytracing appraoch, from each triangle centroid along reflection vector, to find next
% triangle hit by reflected ray, up to MAX_REFL
%
% commented lines assist with debugging.
%
%
%close all hidden;

%simulate_v17;  % for debug, add a "return" in the simulate code just after the light has been turned on, after trinorms
                % and specular optics have been calculated.

%n_xyz = [n_x' n_y' n_z'];  % node position vector

%t_cxyz = [t_cx' t_cy' t_cz'];  %ray's position vector

% ray's direction:    % t_refln'

t_vert0 = [n_x(t_na); n_y(t_na); n_z(t_na)];
t_ev1;
t_ev2;

%t_nbxyz = [n_x(t_nb); n_y(t_nb); n_z(t_nb)]';
%t_ncxyz = [n_x(t_nc); n_y(t_nc); n_z(t_nc)]';
%t_edg1 = t_ev1';
%t_edg2 = t_ev2';

%t_rofx = zeros(size(t_na));
%t_rofy = zeros(size(t_na));
%t_rofz = zeros(size(t_na));

%hold on;

mr_nhits = zeros(1,MAX_REFL);
t_mr_na = zeros(size(t_na));
t_mr_nh = zeros(size(t_na));
t_oth_mr = zeros(size(t_oth));
t_abs_mr = zeros(size(t_na));

% prettycolors = { 'r' 'g' 'c' 'm' 'y' };
% figure
% hold on;
% %    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
% %  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
%   %  vectors (can be confusing).
% 
% quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
% hold on;
% xlabel('x (mm)');
% ylabel('y (mm)');
% zlabel('z (mm)');
% axis equal;
% hitts = [];

tic;
for ti = 1:length(t_na)
    if t_broken(ti)
            continue
    end
%intersect = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
    ar_xyz = [t_cx(ti); t_cy(ti); t_cz(ti)];
    [intersect, t] = TriangleRayIntersectionMK( ar_xyz , t_refln(:,ti), t_vert0, t_ev1, t_ev2);
%  intersect = TriangleRayIntersection( t_cxyz, t_refln', t_naxyz, t_nbxyz, t_ncxyz );
    myhitts = find(intersect);
    if myhitts 
        myhittvals = t(intersect);
        rr_srct = ti;
        rr_dest = myhitts(1);
        rr_dir = t_refln(:,ti);
        rr_pwr = t_reflp(ti);
        ri = 1;
        while 1
            if t_broken(rr_dest), break, end
            ar_costheta = dot(t_norms(:,rr_dest), rr_dir);
            if ar_costheta < 0,  break, end  %warning('Negative angle!'),
            mr_nhits(ri) = mr_nhits(ri) + 1;
            t_mr_nh(ti) = t_mr_nh(ti) +1;
            t_mr_na(rr_dest) = t_mr_na(rr_dest) +1;
%            hitts(end+1) = myhittvals(1); % for debug only
%            mycolor = prettycolors{ mod(mr_nhits(1)-1, 5)+1 };
            
%             plot3( ar_xyz(1) + myhittvals(1).* rr_dir(1) .* [0 1], ...
%                    ar_xyz(2) + myhittvals(1).* rr_dir(2) .* [0 1], ...
%                    ar_xyz(3) + myhittvals(1).* rr_dir(3) .* [0 1], mycolor);
%  
            if useSpecularTMM
                ar_theta = acos(ar_costheta);
                %if min(ar_theta) < 0 
                %    warning('Negative theta in useSpecularTMM code');
                 %   t_theta(t_theta < 0) = 0;
                %end
                %if max(t_theta) > pi/2
                %    warning('Reverse illuminated triangles in specularTMM code');
                %    t_theta(t_theta > pi/2) = pi/2 - t_theta(t_theta > pi/2);
                %end
                ar_LULi = max(min(floor(ar_theta ./ LULAngleStep)+1,LUL.num),1);
                
                ar_specA = LUL.A(ar_LULi);
                ar_specR = LUL.R(ar_LULi);
                
                ar_abs = ar_specA .* rr_pwr;
                ar_rpwr = rr_pwr * ar_specR;
            else
                ar_abs = Iabs .* rr_pwr; %watt
                ar_rpwr = rr_pwr * Irefl; % watt
            end

            t_abs_mr(rr_dest) = t_abs_mr(rr_dest) + ar_abs; % Watt
 
            % calculate optical force ("optical thrust") for each triangle
            t_oth_mr(:,rr_dest) = t_oth_mr(:,rr_dest) + ...
                2 * ar_rpwr ./ c0 .* t_norms(:,rr_dest) .* ar_costheta  + ... %  this is the reflection force
                ar_abs       ./ c0 .* rr_dir ;    % and this is the absorption force
            
            if ri >= MAX_REFL  % end of loop termination
                break
            end
            
            % calculate intersect point:
            ar_xyz = [t_cx(rr_srct) + myhittvals(1).* rr_dir(1) ;   ...
                      t_cy(rr_srct) + myhittvals(1).* rr_dir(2) ;   ...
                      t_cz(rr_srct) + myhittvals(1).* rr_dir(3) ];
            
            % calculate new reflected light vector
            rr_dir = rr_dir - 2 .* ar_costheta .* t_norms(:,rr_dest);
            
            % calculate new reflected ray intersections  %[t_cx(myhitts(1)); t_cy(myhitts(1)); t_cz(myhitts(1))]
            
            [intersect, t] = TriangleRayIntersectionMK(  ar_xyz ,...
                                    rr_dir, t_vert0, t_ev1, t_ev2 );   % this version traces from centroids
                                    
            %[intersect, t] = TriangleRayIntersectionMK( ar_xyz, rr_dir', t_vert0, t_edg1, t_edg2, 'border','exclusive' );   % this version traces from point of reflection.
            
            
            myhitts = find(intersect);
            if myhitts
                myhittvals = t(intersect);
                rr_srct = rr_dest;
                rr_dest = myhitts(1);
                %rr_dir = t(refln(:,ti));
                rr_pwr = ar_rpwr;
                ri = ri+1;
            else
                break;
            end
        end
    end
      
    %fprintf('%g  ', myhitts); fprintf('     ');  fprintf('%g  ',myhittvals);  fprintf('\n');
end

%toc
% fprintf('Hits:  '); fprintf('%d  ', mr_nhits); fprintf('\n');
% 
% cdata=t_mr_na;
% trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
% hcolorbar = colorbar;
% hcolorbar.Label.String = '#ray hits per triangle';
% view([0 0]);
% axis tight;
% set(gca,'YLim',[0 2*radiusmm]);
% maxoth = max(vecnorm(t_oth));
% maxoth_mr = max(vecnorm(t_oth_mr));
% mr_ratio = maxoth_mr / maxoth;
% quiver3(t_cx, t_cy, t_cz, t_oth(1,:), t_oth(2,:), t_oth(3,:), 'm');
% quiver3(t_cx, t_cy, t_cz, t_oth_mr(1,:), t_oth_mr(2,:), t_oth_mr(3,:), mr_ratio, 'b');
