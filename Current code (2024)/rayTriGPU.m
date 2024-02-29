function [dist, flag] = rayTriGPU(P0x, P0y, P0z, ...
                        P1x, P1y, P1z, ...
                        P2x, P2y, P2z, ...
                        orx, ory, orz, ...
                        Dx, Dy, Dz)

% Ray-triangle intersection algorithm of Muller and Trumbore (1997)
% formatted for arrayfun to allow hardware acceleration

% INPUT: 
% P0x, P0y, P0z / P1x, P1y, P1z / P2x, P2y, P2z: xyz components of the
% triangle objects
% orx, ory, orz: xyz componants of the ray origin
% Dx, Dy, Dz: xyz components of the ray directional (unit) vectors

% OUTPUT:
% dist: distance from from the ray-tri intersection to the origin or nan 
% if no intersection is located

% Usage example: 
% Step 1: convert mx3 direction vectors, D = [Dx Dy Dz] to gpuarray object
% >> gD = gpuArray(D);
% Step 2: call rayTriGPU using arrayfun with scalar input formatting
% where P0, P1, P2 are the nx3 vertex lists of the triangle corner points
% and where or is the xyz coordinates of the origin
% >> [dist, flag] = arrayfun(@rayTriGPU, P0(:,1)', P0(:,2)', P0(:,3)', ...
%                             P1(:,1)', P1(:,2)', P1(:,3)', ...
%                             P2(:,1)', P2(:,2)', P2(:,3)', ...
%                             or(:,1), or(:,2), or(:,3), ...
%                             gD(:,1),gD(:,2),gD(:,3));
% Step 3: recover data
% distances = gather(dist);
% Output is an mxn array containing a the distance from the ray-tri
% intersection to the origin or nan if no intersection is located
% Implentation based upon that of Paul Peeling (originally from Jesus P.
% Mena-Chalco of FEX), MathWorks (which returns a flag but not the 
% intersection distance).

% Per ray f lags can be obtained from the output dist using the following 
% method:
% >> flagT = true(size(D,1),1); 
% >> flagT(sum(isnan(dist),2) == size(P0,1)) = false;
% This may save transfer time off the GPU

% Dependencies: requires Parallel Computing Toolbox

epsilon = 0.00001; %rounding error

% E1 = P1-P0;
E1x = P1x-P0x;
E1y = P1y-P0y;
E1z = P1z-P0z;

% E2 = P2-P0;
E2x = P2x-P0x;
E2y = P2y-P0y;
E2z = P2z-P0z;

[Qx,Qy,Qz] = scalarCross(Dx,Dy,Dz,E2x,E2y,E2z);
A = E1x*Qx + E1y*Qy + E1z*Qz;
if (A > -epsilon && A < epsilon)
    % the vector is parallel to the plane (the intersection is at infinity)
    flag = false; dist = nan;
    return
end

F = 1/A;

% S = or-P0;
Sx = orx - P0x;
Sy = ory - P0y;
Sz = orz - P0z;

% U = F.*sum(S.*Q,2);
U = F*(Sx*Qx + Sy*Qy + Sz*Qz);

% flagV(U<0.0) = 0;
if (U < 0.0)
    % the intersection is outside of the triangle
    flag = false; dist = nan;
    return
end

% R = cross(S,E1,2);
[Rx,Ry,Rz] = scalarCross(Sx,Sy,Sz,E1x,E1y,E1z);

% V = F.*sum(D.*R,2);
V = F*(Dx*Rx + Dy*Ry + Dz*Rz);

    
if (V < 0.0 || U+V > 1.0)
    % the intersection is outside of the triangle
    flag = false; dist = nan;
    return
end

% triangle intersection
flag = true;
dist = F*(E2x*Rx + E2y*Ry + E2z*Rz); 
return

% scalar cross
function [w1,w2,w3] = scalarCross(u1,u2,u3,v1,v2,v3)
w1 = u2*v3 - u3*v2;
w2 = u3*v1 - u1*v3;
w3 = u1*v2 - u2*v1;
end

end
