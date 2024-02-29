function [flag, t] = RTI_MK_GPU (ox,oy,oz,...
    dx,dy,dz,...
    p0x,p0y,p0z,...
    e1x,e1y,e1z,...
    e2x,e2y,e2z)

% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997).

    eps1 = 0.000001;
    eps2 = 1e-6;
    t=0;
    flag = false;

    % switched to using pre-calculated edge vectors as input rather than calculating them here.
    % e1 = p1-p0;
    %e1x = p1x-p0x;
    %e1y = p1y-p0y;
    %e1z = p1z-p0z;
    % e2 = p2-p0;
    %e2x = p2x-p0x;
    %e2y = p2y-p0y;
    %e2z = p2z-p0z;
        
    % q  = cross(d,e2);
    [qx,qy,qz] = gpu_cross_product(dx,dy,dz,e2x,e2y,e2z);  % this is pvec in prior version
    
    % a  = dot(e1,q); % determinant of the matrix M
    a = e1x*qx + e1y*qy + e1z*qz;  % this is 'det' in the prior version
    
    if (a>-eps1 && a<eps1) 
        % the vector is parallel to the plane (the intersection is at infinity)
        %flag = false;
        return
    end
    
    f = 1/a;
    % s = o-p0;
    sx = ox - p0x;    % this is tvec in prior code
    sy = oy - p0y;
    sz = oz - p0z;
    
    % u = f*dot(s,q);
    u = f*(sx*qx + sy*qy + sz*qz);  % this is first barycentric coord value
    
    if (u<eps2)% || (u > 1-eps2)
        % the intersection is outside of the triangle, in the -u direction.  (MK:  Can we quit for u > 1 too?  doesn't seem to matter.)
        %flag = false;
        return
    end
          
    % r = cross(s,e1);
    [rx,ry,rz] = gpu_cross_product(sx,sy,sz,e1x,e1y,e1z);  % this is qvec in prior code
    
    % v = f*dot(d,r);
    v = f*(dx*rx + dy*ry + dz*rz);  % this is second barycentric coord
    
    if (v<eps2 || u+v>1.0-eps2)
        % the intersection is outside of the triangle
        %flag = false;
        return;
    end

    t = f*(e2x*rx + e2y*ry + e2z*rz);
    if (t<eps2)  % forward rays only!  
        %flag = false;
        return;
    end
    
    flag = true;
    return;
end

function [w1,w2,w3] = gpu_cross_product(u1,u2,u3,v1,v2,v3)
w1 = u2*v3 - u3*v2;
w2 = u3*v1 - u1*v3;
w3 = u1*v2 - u2*v1;
end