
%close all hidden;

simulate_v17;  % for debug, this poops out after light has just turned on.  Not sure if 

%n_xyz = [n_x' n_y' n_z'];  % node position vector

%t_cxyz = [t_cx' t_cy' t_cz'];  %ray's position vector

% ray's direction:    % t_refln'
tic;
t_vert0 = [n_x(t_na); n_y(t_na); n_z(t_na)]';
%t_nbxyz = [n_x(t_nb); n_y(t_nb); n_z(t_nb)]';
%t_ncxyz = [n_x(t_nc); n_y(t_nc); n_z(t_nc)]';
t_edg1 = t_ev1';
t_edg2 = t_ev2';

t_rofx = zeros(size(t_na));
t_rofy = zeros(size(t_na));
t_rofz = zeros(size(t_na));

t_rof = zeros(3,length(t_na));
hold on;
MAX_REFL = 10;
mr_nhits = zeros(1,MAX_REFL);
t_mr_na = zeros(size(t_na));
t_mr_nh = zeros(size(t_na));
t_othr = zeros(size(t_oth));

prettycolors = { 'r' 'g' 'b' 'k' 'm' 'y' };

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

hitts = [];

for ti = 1:length(t_na)
    if t_broken(ti)
            continue
    end
%intersect = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
    [intersect, t, u, v, xcoor] = TriangleRayIntersectionMK( [t_cx(ti) t_cy(ti) t_cz(ti)] , t_refln(:,ti)', t_vert0, t_edg1, t_edg2, 'border','exclusive' );
%  intersect = TriangleRayIntersection( t_cxyz, t_refln', t_naxyz, t_nbxyz, t_ncxyz );
    myhitts = find(intersect);
    myhittvals = t(intersect);
    if myhitts 
        if t_broken(myhitts(1))
            continue
        end
        
        ri = 1;
        while 1
            
            mr_nhits(ri) = mr_nhits(ri) + 1;
            t_mr_nh(ti) = t_mr_nh(ti) +1;
            t_mr_na(myhitts(1)) = t_mr_na(myhitts(1)) +1;
            hitts(end+1) = myhittvals(1); % for debug only
            mycolor = prettycolors{ mod(mr_nhits(1)-1, 6)+1 };
            plot3( t_cx(ti) + myhittvals(1).* t_refln(1,ti) .* [0 1], ...
                   t_cy(ti) + myhittvals(1).* t_refln(2,ti) .* [0 1], ...
                   t_cz(ti) + myhittvals(1).* t_refln(3,ti) .* [0 1], mycolor);
               %todo:  distribute force
            ar_costheta = dot(t_norms(:,myhitvals(1)), t_refln(:,ti));
            if ar_costheta < 0, warning('Negative angle!'), end
            %ar_aproj = t_a(myhitts(1)) .* ar_costheta ;
            
            
        
        
           
           
        [intersect, t, u, v, xcoor] = TriangleRayIntersectionMK( [t_cx(myhitts(1)) t_cy(myhitts(1)) t_cz(myhitts(1))] , t_refln(:,ti)', t_vert0, t_edg1, t_edg2, 'border','exclusive' );   
        myhitts2 = find(intersect);
        myhittvals2 = t(intersect);   
        if myhitts2 
          nhits2 = nhits2 + 1;
            t_mr_nh(myhitts(1)) = t_mr_nh(myhitts(1)) +1;
            t_mr_na(myhitts2(1)) = t_mr_na(myhitts2(1)) +1;
            hitts(end+1) = myhittvals(1);
            plot3( t_cx(myhitts(1)) + myhittvals2(1).* t_refln(1,myhitts(1)) .* [0 1], ...
                   t_cy(myhitts(1)) + myhittvals2(1).* t_refln(2,myhitts(1)) .* [0 1], ...
                   t_cz(myhitts(1)) + myhittvals2(1).* t_refln(3,myhitts(1)) .* [0 1], mycolor);
        end
    end
    %fprintf('%g  ', myhitts); fprintf('     ');  fprintf('%g  ',myhittvals);  fprintf('\n');
end

toc
disp(['First hits: ' num2str(mr_nhits)]);
disp(['Second hits: ' num2str(nhits2)]);

cdata=t_mr_na;
trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
colorbar;