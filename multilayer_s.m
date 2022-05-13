function [rho,tao,R,T] = multilayer_s(d,y,y_inc,y_sub,lambda,theta)
% This code is for s-polarization.
for m = 1:length(theta),
for i = 1:length(y),
y_s(m,i) = sqrt(real(y(i))^2-imag(-y(i))^2 - y_inc^2*sind(theta(m))^2-2*1i*real(y(i))*imag(-y(i)));
if imag(y_s(m,i))>0,
y_s(m,i) = -(y_s(m,i));
end
delta_s(m,i) = 2*pi*d(i)/lambda*y_s(m,i);
y_inc_s(m) = y_inc*cosd(theta(m));
y_sub_s(m) = sqrt(real(y_sub)^2-imag(-y_sub)^2 - y_inc^2*sind(theta(m))^2-2*1i*real(y_sub)*imag(-y_sub));
if imag(y_sub_s(m))>0,
y_sub_s(m) = -(y_sub_s(m));
end




tmp = [cos(delta_s(m,i)), 1i*sin(delta_s(m,i))/y_s(m,i);
1i*y_s(m,i)*sin(delta_s(m,i)), cos(delta_s(m,i))];
if i == 1,
M = tmp;
else
M = M*tmp;
end
end
BC = M*[1; y_sub_s(m)];
B(m) = BC(1);
C(m) = BC(2);
% rho and tao are the amplitude reflection and transmission coefficients
rho(m) = (y_inc_s(m)*B(m)-C(m))/(y_inc_s(m)*B(m)+C(m));
phi_r(m) = angle(rho(m))/pi*180; % phase
tao(m) = 2*y_inc_s(m)/(y_inc_s(m)*B(m)+C(m));
phi_t(m) = angle(tao(m))/pi*180; % phase
% R and T are the intensity reflection and transmission coefficients
% Note: for transmission, the angle changes due to refraction, thus area is different
R(m) = abs(rho(m))^2;
T(m) = 4*y_inc_s(m)*real(y_sub_s(m))/abs((y_inc_s(m)*B(m)+C(m)))^2;
end