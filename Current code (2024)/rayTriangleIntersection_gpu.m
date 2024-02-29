function flag = rayTriangleIntersection_gpu (ox,oy,oz,...
    dx,dy,dz,...
    p0x,p0y,p0z,...
    p1x,p1y,p1z,...
    p2x,p2y,p2z)

% Ray/triangle intersection using the algorithm proposed by Möller and Trumbore (1997).

    epsilon = 0.00001;

    % e1 = p1-p0;
    e1x = p1x-p0x;
    e1y = p1y-p0y;
    e1z = p1z-p0z;
    % e2 = p2-p0;
    e2x = p2x-p0x;
    e2y = p2y-p0y;
    e2z = p2z-p0z;
        
    % q  = cross(d,e2);
    [qx,qy,qz] = cross_product(dx,dy,dz,e2x,e2y,e2z);
    
    % a  = dot(e1,q); % determinant of the matrix M
    a = e1x*qx + e1y*qy + e1z*qz;
    
    if (a>-epsilon && a<epsilon) 
        % the vector is parallel to the plane (the intersection is at infinity)
        flag = false;
        return
    end;
    
    f = 1/a;
    % s = o-p0;
    sx = ox - p0x;
    sy = oy - p0y;
    sz = oz - p0z;
    
    % u = f*dot(s,q);
    u = f*(sx*qx + sy*qy + sz*qz);
    
    if (u<0.0)
        % the intersection is outside of the triangle
        flag = false;
        return
    end;
          
    % r = cross(s,e1);
    [rx,ry,rz] = cross_product(sx,sy,sz,e1x,e1y,e1z);
    
    % v = f*dot(d,r);
    v = f*(dx*rx + dy*ry + dz*rz);
    
    if (v<0.0 || u+v>1.0)
        % the intersection is outside of the triangle
        flag = false;
        return;
    end;
    
    flag = true;
    return;
end

function [w1,w2,w3] = cross_product(u1,u2,u3,v1,v2,v3)
w1 = u2*v3 - u3*v2;
w2 = u3*v1 - u1*v3;
w3 = u1*v2 - u2*v1;
end