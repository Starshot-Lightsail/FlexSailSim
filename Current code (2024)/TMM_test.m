% initial program
%clear all;
%close all;
um = 1e-6;
nm = 1e-9;
% define admittance for medium
Air = 1;
yH = 2.35;
yL = 1.35;
Glass = 1.52;
% define the reference wavelength and thickness
lambda_f = 480*nm;
dH = lambda_f/4/yH; % quarter stack
dL = lambda_f/4/yL;
% define visible region of light
lambda = linspace(350,850,501)*nm;
% define thin film structure
y_inc = Air; % incident medium admittance
y_sub = Glass; % substrate medium admittance
theta = 0; % in degree
d = [dH,dL,dH,dL,dH,dL,dH,dL,dH,dL,dH,1.2*dL,1.4*dH,1.4*dL,1.4*dH,1.4*dL,...
1.4*dH,1.4*dL,1.4*dH,1.4*dL,1.4*dH,1.4*dL,1.4*dH];
y = [yH,yL,yH,yL,yH,yL,yH,yL,yH,yL,yH,yL,yH,yL,yH,yL,...
yH,yL,yH,yL,yH,yL,yH];
for ii=1:length(lambda),
[rho(ii),tao(ii),R(ii),T(ii)] = multilayer_s(d,y,y_inc,y_sub,lambda(ii),theta);
end
% plot the result
% figure, plot(lambda/nm,R*100,'-k','LineWidth',3);
% title('A broadband reflector for the visible region');
% xlabel('Wavelength(nm)');
% ylabel('Reflectance(%)');
% axis([350 850 0 105]);
% grid on;

y_inc = 1;
lambda = 1064*1e-6; %mm
myn = 3.52 - 1i * 1e-6;
myt = 43 .* 1e-6; %mm
myt = 40 .* 1e-6;
y_sub = 1;
thetas = deg2rad(0:89);

%[ rho, tao, R, T ] = multilayer_s( myt, myn, 1, 1, lambda, thetas );

%for m = 1:length(theta),
%y_s = zeros(size(thetas));
%for i = 1:length(y),
y_s = sqrt( real(myn)^2 - imag(-myn)^2 - 1.*sin(thetas).^2 - 2*1i*real(myn)*imag(-myn) );
%if imag(y_s(m,i))>0,
y_s(imag(y_s)>0) = -y_s(imag(y_s)>0);
%end
delta_s = 2*pi*myt/lambda*y_s;
y_0_s = 1*cos(thetas);
%y_0_s = y_0_s; %sqrt(1 - sind(thetas).^2 );
%if imag(y_sub_s(m))>0,
%y_sub_s(m) = -(y_sub_s(m));
%end
y_p = myn^2./y_s;
delta_p = delta_s;
y_0_p = 1./y_0_s;
%y_sub_p(m) = y_sub^2/y_sub_s(m);
% tmp = [cos(delta_p), 1i*sin(delta_p)/y_p;
% 1i*y_p(m,i)*sin(delta_p(m,i)), cos(delta_p(m,i))];
% if i == 1,
% M = tmp;
% else
% M = M*tmp;
% end
% end
% BC = M*[1; y_sub_p(m)];

% %tmp = [cos(delta_s(m,i)), 1i*sin(delta_s(m,i))/y_s(m,i);
% %1i*y_s(m,i)*sin(delta_s(m,i)), cos(delta_s(m,i))];
%if i == 1,
%M = tmp;
%else
%M = M*tmp;
%end

Bs = cos(delta_s) + 1i*sin(delta_s)./y_s.*y_0_s;
Cs = 1i*y_s.*sin(delta_s) + cos(delta_s).*y_0_s;

Bp = cos(delta_s) + 1i*sin(delta_s)./y_p.*y_0_p;
Cp = 1i*y_p.*sin(delta_s) + cos(delta_s).*y_0_p;


%end
%BC = M*[1; y_sub_s(m)];
%B(m) = BC(1);
%C(m) = BC(2);
% rho and tao are the amplitude reflection and transmission coefficients
rho_s = (y_0_s.*Bs-Cs)./(y_0_s.*Bs+Cs);
rho_p = (y_0_p.*Bp-Cp)./(y_0_p.*Bp+Cp);
%phi_r = angle(rho(m))/pi*180; % phase
tao_s = 2*y_0_s./(y_0_s.*Bs+Cs);
tao_p = 2*y_0_s./(y_0_s.*Bp+Cp);
%tao_p = 2*y_inc_s(m)/(y_inc_s(m)*B(m)+C(m));
%phi_t(m) = angle(tao(m))/pi*180; % phase
% R and T are the intensity reflection and transmission coefficients
% Note: for transmission, the angle changes due to refraction, thus area is different
Rs = abs(rho_s).^2;
Rp = abs(rho_p).^2;
Ts = 4*y_0_s.^2 ./ abs(y_0_s.*Bs+Cs).^2;
Tp = 4*y_0_p.^2 ./ abs(y_0_p.*Bp+Cp).^2;
As = 1 - Rs - Ts;
Ap = 1 - Rp - Tp;

figure 
plot(rad2deg(thetas), Rs );
hold on
plot(rad2deg(thetas), Rp );
plot(rad2deg(thetas), 0.5*Rp + 0.5*Rs );
xlabel('\theta (\circ)');
ylabel('R');

figure 
plot(rad2deg(thetas), Ts );
hold on
plot(rad2deg(thetas), Tp );
plot(rad2deg(thetas), 0.5*Tp + 0.5*Ts );
xlabel('\theta (\circ)');
ylabel('T');

figure 
plot(rad2deg(thetas), As );
hold on
plot(rad2deg(thetas), Ap );
plot(rad2deg(thetas), 0.5*Ap + 0.5*As );
xlabel('\theta (\circ)');
ylabel('A');