function [intersect, t] = TriangleRayIntersectionMK(  ori, dir, vert0, edge1, edge2 )

% determines ray-triangle intersection based on method of Möller and Trumbore ('97)
%
% ori: ray origin (one vector)
% dir: ray's direction (one vector)
% vert0: primary triangle vertices (many vectors)
% edge1, edge2:  Edge vectors spanning the triangle (many vectors)

numtris = length(vert0);
% stupidly duplicate the origin vector so we can vectorize the code
ori  = repmat(ori , 1, numtris);
% same for direction vector
dir   = repmat(dir, 1, numtris);


% disp(['ORI ' num2str(size(ori))]);
% disp(['DIR ' num2str(size(dir))]);
% disp(['vert0 ' num2str(size(vert0))]);
% disp(['edg1 ' num2str(size(edge1))]);
% disp(['edg2 ' num2str(size(edge2))]);

% floating point math sucks
eps        = 1e-5;

% initialize default output
intersect = false(1, numtris); % by default there are no intersections
t = inf+zeros(1, numtris); u=t; v=t;
%xcoor = nan+zeros(size(ori)); delete me

% find faces parallel to the ray
tvec  = ori - vert0;          % vector from vert0 to ray origin
pvec  = cross(dir, edge2);  % begin calculating determinant - also used to calculate U parameter
det   = sum(edge1.*pvec);   % determinant of the matrix M = dot(edge1,pvec)

angleOK = (abs(det)>eps); % if determinant is near zero then ray lies in the plane of the triangle
    
if all(~angleOK), disp('tacos'); disp(' '); return; end % if all parallel than no intersections

%% Different behavior depending on one or two sided triangles
det(~angleOK) = nan;              % change to avoid division by zero
u    = sum(tvec.*pvec)./det;    % 1st barycentric coordinate
% disp(['DET ' num2str(size(det))]);
% disp(['U ' num2str(size(u))]);
% disp(['TVEC ' num2str(size(tvec))]);
% disp(['AOK ' num2str(size(angleOK))]);

  v = nan+zeros(size(u)); t=v;
%  disp(['V ' num2str(size(v))]);
  ok = angleOK & u>=eps & u<=1.0-eps; % mask
%  disp(['OK ' num2str(size(ok))]);
  % if all line/plane intersections are outside the triangle than no intersections
  if ~any(ok), intersect = ok;  disp('burritos'); disp(' '); return; end

  qvec = cross(tvec(:,ok), edge1(:,ok)); % prepare to test V parameter
  v(:,ok) = sum(dir(:,ok).*qvec) ./ det(:,ok); % 2nd barycentric coordinate
  t(:,ok) = sum(edge2(:,ok).*qvec)./det(:,ok);

     
  % test if line/plane intersection is within the triangle
  ok = (ok & v>=eps & u+v<=1.0-eps);

  intersect = (ok & t>=eps); % return only thos with intersections on the correct side of origin.
