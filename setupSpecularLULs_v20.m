% setupSpecularLULs
% (c) MK 2023 CIT

%here, we pre-calculate the specular reflectance of each available material 

% Inputs from workspace:
% LULAngleStep - step size for LUL (in radians)

%if not specified, use a default for LUTANgleStep
if ~exist('LULAngleStepDeg')
    LULAngleStepDeg = 0.2;
end
LULAngleStep = deg2rad(LULAngleStepDeg);

if exist('LULs')
    clear LULs;
end
if exist('LUL')
    clear LUL;
end
for luln=1:length(MATERIALS)
    LULs(luln) = getSpectralLUL(MATERIALS(luln).name,  MATERIALS(luln).thickness, MATERIALS(luln).ior_n, MATERIALS(luln).ior_k, lambda, LULAngleStepDeg);
end
if ~isempty(luln)
    LUL.R = LULs(1).R;
    LUL.A = LULs(1).A;
    LUL.angleStep = LULs(1).angleStep;
    LUL.angles = LULs(1).angles;
    LUL.num = LULs(1).num;
else
    error(['ERROR: Attempting to set up LULs for TMM specular reflection model, but there are no materials defined!  (Need to specify at least one material in setupMaterialProperties)']);
end
for luln=2:length(MATERIALS)
    LUL.R(luln,:) = LULs(luln).R;
    LUL.A(luln,:) = LULs(luln).A;
end
LUL.sz = size(LUL.R);

% moved to helper function to keep intermediate calculation variables out of the workspace
function LUL = getSpectralLUL(material, thickness, ior_n, ior_k, lambda, LULAngleStepDeg)

showLULPlots = 0;

disp(' ');
disp('setupSpecularLULs v18,  Sept 8 2021, (c) 2021 Michael Kelzenberg, California Institute of Technology');
disp(['Material:  ' material]);
disp(['Thickness:  ' num2str(thickness*1e6) ' nm']);
disp(['IOR:  n=' sprintf('%-6g', ior_n) '  k=' sprintf('%-6g',ior_k) ]);
disp(['lambda:  ' num2str(lambda*1e6) ' nm']);
disp(['Angle step: ' num2str(LULAngleStepDeg) ' deg']);


LUL = [];


LULAngles = deg2rad(  unique( [ 0:LULAngleStepDeg:90  90 ] ) );

myn = ior_n - 1i * ior_k;
myt = thickness; %mm


y_s = sqrt( real(myn)^2 - imag(-myn)^2 - 1.*sin(LULAngles).^2 - 2*1i*real(myn)*imag(-myn) );

y_s(imag(y_s)>0) = -y_s(imag(y_s)>0);

delta_s = 2*pi*myt/lambda*y_s;
y_0_s = 1*cos(LULAngles);

y_p = myn^2./y_s;
delta_p = delta_s;
y_0_p = 1./y_0_s;



Bs = cos(delta_s) + 1i*sin(delta_s)./y_s.*y_0_s;
Cs = 1i*y_s.*sin(delta_s) + cos(delta_s).*y_0_s;

Bp = cos(delta_s) + 1i*sin(delta_s)./y_p.*y_0_p;
Cp = 1i*y_p.*sin(delta_s) + cos(delta_s).*y_0_p;

% rho and tao are the amplitude reflection and transmission coefficients
rho_s = (y_0_s.*Bs-Cs)./(y_0_s.*Bs+Cs);
rho_p = (y_0_p.*Bp-Cp)./(y_0_p.*Bp+Cp);
%phi_r = angle(rho(m))/pi*180; % phase
tao_s = 2*y_0_s./(y_0_s.*Bs+Cs);
tao_p = 2*y_0_s./(y_0_s.*Bp+Cp);
%tao_p = 2*y_inc_s(m)/(y_inc_s(m)*B(m)+C(m));

% R and T are the intensity reflection and transmission coefficients

Rs = abs(rho_s).^2;
Rp = abs(rho_p).^2;
Ts = 4*y_0_s.^2 ./ abs(y_0_s.*Bs+Cs).^2;
Tp = 4*y_0_p.^2 ./ abs(y_0_p.*Bp+Cp).^2;
As = 1 - Rs - Ts;
Ap = 1 - Rp - Tp;

LUL.angles = LULAngles;
LUL.angleStepDeg = deg2rad(LULAngleStepDeg);
LUL.angleStep = deg2rad(LULAngleStepDeg);
LUL.R = 0.5*Rp + 0.5*Rs;
LUL.A = 0.5*Ap + 0.5*As;
LUL.T = 0.5*Tp + 0.5*Ts;
LUL.num = length(LULAngles);

S = whos('LUL');
disp(['Total LUL memory footprint:  ' num2str(S.bytes/1024) ' kB']);
fprintf('\n');

if showLULPlots

figure('Position',[204   236   624   970]);

subplot(3,1,1)
plot(rad2deg(LULAngles), Rs );
hold on
plot(rad2deg(LULAngles), Rp );
plot(rad2deg(LULAngles), 0.5*Rp + 0.5*Rs );
xlabel('\theta (\circ)');
ylabel('R');
legend({ 'P', 'S', 'LUL' }, 'Location', 'northwest');
hold off;
title(['SpecularTMM:  ' material], 'Interpreter','none');

subplot(3,1,2)
plot(rad2deg(LULAngles), Ts );
hold on
plot(rad2deg(LULAngles), Tp );
plot(rad2deg(LULAngles), 0.5*Tp + 0.5*Ts );
xlabel('\theta (\circ)');
ylabel('T');
hold off;

subplot(3,1,3)
plot(rad2deg(LULAngles), As );
hold on
plot(rad2deg(LULAngles), Ap );
plot(rad2deg(LULAngles), 0.5*Ap + 0.5*As );
xlabel('\theta (\circ)');
ylabel('A');
hold off;

end

end