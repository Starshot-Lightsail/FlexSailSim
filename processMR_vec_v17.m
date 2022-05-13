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

showplots = 0;

if showplots
prettycolors = { 'r' 'g' 'c' 'm' 'y' };
figure
hold on;
%    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
%  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
  %  vectors (can be confusing).

quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
hold on;
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal;
end       
                
%n_xyz = [n_x' n_y' n_z'];  % node position vector

%t_cxyz = [t_cx' t_cy' t_cz'];  %ray's position vector

% ray's direction:    % t_refln'
tic
t_indxs = 1:numtris;
numtrisgood = sum(t_notbroken);
numrays = numtrisgood;
src_tris  = t_indxs(t_notbroken);
mr_nhits = zeros(1,MAX_REFL);
%t_mr_na = zeros(1, length(t_na));
t_mr_nh = zeros(size(t_na));  %how many rays reach each tri?
t_oth_mr = zeros(size(t_oth));
t_abs_mr = zeros(size(t_na));

ri = 1;

t_vert0f = [n_x(t_na(t_notbroken)); n_y(t_na(t_notbroken)); n_z(t_na(t_notbroken))];
t_ev1f = t_ev1(:,t_notbroken);
t_ev2f = t_ev2(:,t_notbroken);

t_cxyz = [t_cx(t_notbroken); t_cy(t_notbroken); t_cz(t_notbroken)];

ray_origins = t_cxyz;
ray_dirs = t_refln(:, t_notbroken);
ray_pwrs = t_reflp(:, t_notbroken);

t_vert0_exp = repmat(t_vert0f, 1, numrays);
t_ev1_exp = repmat(t_ev1f, 1, numrays );
t_ev2_exp = repmat(t_ev2f, 1, numrays );

% start loop somehow
while 1
    % ok, replicate the arrays
    numcombs = numrays*numtrisgood;
    ray_origins_exp = repelem(ray_origins, 1, numtrisgood);
    ray_dirs_exp = repelem(ray_dirs, 1, numtrisgood);
    
    %  disp(['ORI ' num2str(size(ray_origins_exp))]);
    %  disp(['DIR ' num2str(size(ray_dirs_exp))]);
    %  disp(['vert0 ' num2str(size(t_vert0_exp))]);
    %  disp(['edg1 ' num2str(size(t_ev1_exp))]);
    %  disp(['edg2 ' num2str(size(t_ev2_exp))]);
    
    [intersect, t] = TriangleRayIntersectionMultiMK( ray_origins_exp, ray_dirs_exp, t_vert0_exp(:,1:(numcombs)), t_ev1_exp(:,1:(numcombs)), t_ev2_exp(:,1:(numcombs)));
    
    intersectArr = reshape(intersect,[ numtrisgood numrays]);
    %tArr = reshape(t,[ numtrisgood numrays]);
    
    
    
    [hit_mask_ar, hitidxs_ar] = max(intersectArr);
    hit_ts = t( sub2ind([ numtrisgood  numrays ],  hitidxs_ar, 1:length(hitidxs_ar)   ) );
    %hit_ts = tArr( sub2ind([ numtrisgood numrays ],   1:length(hitidxs_ar), hitidxs_ar  ) )
    
  
    
    hit_mask_ar = (hit_mask_ar > 0);  % do I need this?  I think so, max() could reteturn values of 2 or more ...
    ray_dest_idxs = hitidxs_ar(hit_mask_ar);
    
    
    % OK do actual optics physics stuff here..
    myhittvals = t(intersect);
    %rr_srct = ti;
    %rr_dest = myhitts(1);
    %rr_dir = t_refln(:,ti);
    %rr_pwr = t_reflp(ti);
    %ri = 1;
    %while 1
    %if t_broken(rr_dest), break, end
    ar_costheta = dot(t_norms(:,ray_dest_idxs), ray_dirs(:,hit_mask_ar));
    hit_mask_ar( ar_costheta < 0 ) = 0;  % also filter out back-illumianted tris from the source hit mask
    ar_costheta = ar_costheta( ar_costheta >= 0);  % remove back-illuminated from the costheta array
    ray_dest_idxs = hitidxs_ar(hit_mask_ar);  %remove back illuminated from the destination index array
    hit_ts = hit_ts(hit_mask_ar);
    numrays = length(ray_dest_idxs);
    mr_nhits(ri) = numrays; %keep track of number or rays for this iteration
    
    %temp_mat1 = zeros(numrays,numrays);
    %temp_mat2 = zeros(numrays,numrays);
    %temp_mat3 = zeros(numrays,numrays);
    %temp_mat4 = zeros(numrays,numrays);
    
    % Because some triangles receive multiple ray hits, there may be multiple repeated index values in ray_dest_idxs
    % So we have to either construct a big array to expand our calculations, or try figuring out unique() and
    % accumarray() to get the desired behavior.  Not sure which is faster.
    %
    [ray_dest_idxs_unique,idontknow,idxUnique] = unique(ray_dest_idxs);
    t_mr_nh(ray_dest_idxs_unique) = t_mr_nh(ray_dest_idxs_unique) + accumarray(idxUnique, 1)';  % keep track of number of rays reaching each triangle
    
    ray_pwrs = ray_pwrs(hit_mask_ar);  % select only successfull rays
    ray_dirs = ray_dirs(:,hit_mask_ar);  % select only successful rays
    
    if useSpecularTMM  % now calculate optical response and assign forces/power to triangles
        ar_theta = acos(ar_costheta);
        ar_LULi = max(min(floor(ar_theta ./ LULAngleStep)+1,LUL.num),1);
        
        ar_specA = LUL.A(ar_LULi);
        ar_specR = LUL.R(ar_LULi);
        
        ar_abs = ar_specA .* ray_pwrs .* t_notbroken(ray_dest_idxs);
        ar_rpwr = ray_pwrs .* ar_specR .* t_notbroken(ray_dest_idxs);
    else
        ar_abs = Iabs .* ray_pwrs .* t_notbroken(ray_dest_idxs); %watt
        ar_rpwr = ray_pwrs .* Irefl .* t_notbroken(ray_dest_idxs); % watt
    end
    
    t_abs_mr(ray_dest_idxs_unique) = t_abs_mr(ray_dest_idxs_unique) + accumarray(idxUnique, ar_abs)'; % Watt
    
    % calculate optical force ("optical thrust") for each triangle
    t_oth_mr(1,ray_dest_idxs_unique) = t_oth_mr(1,ray_dest_idxs_unique) + ...
        accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
        ar_abs  ./ c0 .* ray_dirs(1,:)    )';    % and this is the absorption force
    t_oth_mr(2,ray_dest_idxs_unique) = t_oth_mr(2,ray_dest_idxs_unique) + ...
        accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
        ar_abs  ./ c0 .* ray_dirs(2,:)    )';    % and this is the absorption force
    t_oth_mr(3,ray_dest_idxs_unique) = t_oth_mr(3,ray_dest_idxs_unique) + ...
        accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
        ar_abs  ./ c0 .* ray_dirs(3,:)    )';    % and this is the absorption force
    
    ri = ri+1;
    
    
    if ri > MAX_REFL || numrays == 0
        break;
    end
   
    if showplots
    np2 = 0;
    for npl = 1:numrays
        if  hit_mask_ar(npl)
            np2 = np2+1;
            plot3( ray_origins(1, npl) + hit_ts(np2).* ray_dirs(1,np2) .* [0 1], ...
                   ray_origins(2, npl) + hit_ts(np2).* ray_dirs(2,np2) .* [0 1], ...
                   ray_origins(3, npl) + hit_ts(np2).* ray_dirs(3,np2) .* [0 1]);
        end
    end
    end
  
    
    
    ray_origins = ray_origins(:,hit_mask_ar) + hit_ts.* ray_dirs ;
              % [t_cx(ray_dest_idxs) + hit_ts.* ray_dirs(1) ;   ...
              %  t_cy(ray_dest_idxs) + hit_ts.* ray_dirs(2) ;   ...
              %  t_cz(ray_dest_idxs) + hit_ts.* ray_dirs(3) ];
    
    % calculate new reflected direction
    %ray_origins = t_cxyz(:,ray_dest_idxs);
    ray_pwrs = ar_rpwr;
    ray_dirs = ray_dirs - 2 .* ar_costheta .* t_norms(:,ray_dest_idxs);
    
end

toc

disp(num2str(mr_nhits))
disp(num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]));
disp(num2str(sum(t_abs_mr)));

if showplots
    %cdata=t_mr_nh
    cdata=t_abs_mr;
    trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    hcolorbar = colorbar;
    %hcolorbar.Label.String = '#ray hits per triangle';
    hcolorbar.Label.String = 'Multibounce absorbed power (W)';
    view([0 0]);
    axis tight;
    set(gca,'YLim',[0 2*radiusmm]);
    maxoth = max(vecnorm(t_oth));
    maxoth_mr = max(vecnorm(t_oth_mr));
    mr_ratio = maxoth_mr / maxoth;
    quiver3(t_cx, t_cy, t_cz, t_oth(1,:), t_oth(2,:), t_oth(3,:), 'm');
    quiver3(t_cx, t_cy, t_cz, t_oth_mr(1,:), t_oth_mr(2,:), t_oth_mr(3,:), mr_ratio, 'b');
end









%hold on;
if 1
    
if showplots
    prettycolors = { 'r' 'g' 'c' 'm' 'y' };
    figure
    hold on;
    %    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
    %  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
    %  vectors (can be confusing).
    
    quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
    hold on;
    xlabel('x (mm)');
    ylabel('y (mm)');
    zlabel('z (mm)');
    axis equal;
end

tic;

mr_nhits = zeros(1,MAX_REFL);
t_mr_na = zeros(size(t_na));
t_mr_nh = zeros(size(t_na));
t_oth_mr = zeros(size(t_oth));
t_abs_mr = zeros(size(t_na));


 hitts = [];


 for ti = 1:length(t_na)
     if t_broken(ti)
         continue
     end
     %intersect = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
     ar_xyz = [t_cx(ti); t_cy(ti); t_cz(ti)];
     [intersect, t] = TriangleRayIntersectionMK( ar_xyz, t_refln(:,ti), t_vert0, t_ev1, t_ev2);
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
             hitts(end+1) = myhittvals(1); % for debug only
             mycolor = prettycolors{ mod(mr_nhits(1)-1, 5)+1 };
             if showplots
                 plot3( ar_xyz(1) + myhittvals(1).* rr_dir(1) .* [0 1], ...
                     ar_xyz(2) + myhittvals(1).* rr_dir(2) .* [0 1], ...
                     ar_xyz(3) + myhittvals(1).* rr_dir(3) .* [0 1], mycolor);
             end
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
             ar_xyz =   ar_xyz + myhittvals(1).* rr_dir;
             %     [t_cx(rr_srct) + myhittvals(1).* rr_dir(1) ;   ...
             %      t_cy(rr_srct) + myhittvals(1).* rr_dir(2) ;   ...
             %      t_cz(rr_srct) + myhittvals(1).* rr_dir(3) ];
             
             % calculate new reflected light vector
             rr_dir = rr_dir - 2 .* ar_costheta .* t_norms(:,rr_dest);
             
             % calculate new reflected ray intersections
             %ar_xyz = [t_cx(myhitts(1)); t_cy(myhitts(1)); t_cz(myhitts(1))];
             
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

toc
%fprintf('Hits:  '); fprintf('%d  ', mr_nhits); fprintf('\n');
 
disp(num2str(mr_nhits))
disp(num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]));
disp(num2str(sum(t_abs_mr)));

if showplots
cdata=t_abs_mr;
%cdata=t_mr_na;
trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
hcolorbar = colorbar;
%hcolorbar.Label.String = '#ray hits per triangle';
hcolorbar.Label.String = 'Multibounce absorbed power (W)';
view([0 0]);
axis tight;
set(gca,'YLim',[0 2*radiusmm]);
maxoth = max(vecnorm(t_oth));
maxoth_mr = max(vecnorm(t_oth_mr));
mr_ratio = maxoth_mr / maxoth;
quiver3(t_cx, t_cy, t_cz, t_oth(1,:), t_oth(2,:), t_oth(3,:), 'm');
quiver3(t_cx, t_cy, t_cz, t_oth_mr(1,:), t_oth_mr(2,:), t_oth_mr(3,:), mr_ratio, 'b');
end

end
