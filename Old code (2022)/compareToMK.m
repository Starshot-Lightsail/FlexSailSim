% facilitates comparisons to MK's rigid body mesh simulator


% to be run in the workspace following a simulation.

%t_cutoff = 0.095371;

%nt1 = floor(interp1(plot_time,1:length(plot_time),0));
nt1 = 1;

%nt2 = ceil(interp1(plot_time, 1:length(plot_time),t_cutoff));
nt2 = length(t);

%close all

pnt = nt1:nt2;
pt = t.*t0;
if 1
figure(1) % position, velocity, force, accel=dv/dt
sps = [];

sps(end+1) = subplot(4,3,1);  % x vs time
plot(pt, (u(pnt,1)).*D_actual.*1000);
hold on;
title('X');
ylabel('Position (mm)');

sps(end+1) = subplot(4,3,2);  % y vs time
plot(pt, (u(pnt, 2)).*D_actual.*1000);
hold on;
title('Y');

sps(end+1) = subplot(4,3,3);  % z vs time
plot(pt, -u(pnt, 3).*D_actual.*1000);
hold on;
title('Z');



Vx_actual = u(pnt, 7).*1000.*sqrt(I0_actual*D_actual^3/(M_actual * c0));
sps(end+1) = subplot(4,3,4);  % vx vs time
plot(pt, u(pnt, 7).*1000.*sqrt(I0_actual*D_actual^3/(M_actual * c0)));
hold on;
ylabel('Vel (mm/s)');

Vy_actual = u(pnt, 8).*1000.*sqrt(I0_actual*D_actual^3/(M_actual * c0));
sps(end+1) = subplot(4,3,5);  % vy vs time
plot(pt, u(pnt, 8).*1000.*sqrt(I0_actual*D_actual^3/(M_actual * c0)));
hold on;

Vz_actual = u(pnt, 9).*1000.*sqrt(I0_actual*D_actual^3/(M_actual * c0));
sps(end+1) = subplot(4,3,6);  % vz vs time
plot(pt, -u(pnt, 9).*1000.*sqrt(I0_actual*D_actual^3/(M_actual * c0)));
hold on;


sps(end+1) = subplot(4,3,7);  % fx vs time
plot(pt, FxLabFrames(pnt).* (I0_actual.*D_actual^2./c0));
hold on;
ylabel('Force (N)');

sps(end+1) = subplot(4,3,8);  % fy vs time
plot(pt, FyLabFrames(pnt).* (I0_actual.*D_actual^2./c0));
hold on;

sps(end+1) = subplot(4,3,9);  % fz vs time
plot(pt, -FzLabFrames(pnt).* (I0_actual.*D_actual^2./c0));  
hold on;


dpt = pt(1:end-1);
dts = diff(pt);

sps(end+1) = subplot(4,3,10);  % ax vs time
plot(dpt, diff(Vx_actual)./dts);
hold on;
ylabel('dv/dt');

sps(end+1) = subplot(4,3,11);  % ay vs time
plot(dpt, diff(Vy_actual)./dts);
hold on;

sps(end+1) = subplot(4,3,12);  % az vs time
plot(dpt, -diff(Vz_actual)./dts);  
hold on;


end



if 1
figure(2) % ang position, ang velocity, torque, diff-ang-vel
sps2 = [];

sps2(end+1) = subplot(4,3,1);  % angle-x vs time
plot(pt, rad2deg(u(pnt,6)));     % phi
hold on;
title('X');
ylabel('Angle (\circ)');

sps2(end+1) = subplot(4,3,2);  % angle-y vs time
plot(pt, rad2deg(u(pnt,5)));   % theta
hold on;
title('Y');

sps2(end+1) = subplot(4,3,3);  % angle-z vs time
plot(pt, mod(rad2deg(-u(pnt,4)),360)-180);   % psi
hold on;
title('Z');



sps2(end+1) = subplot(4,3,4);  % avx vs time
plot(pt, (u(pnt,10))./t0);
hold on;
ylabel('AngVel (rad/s)');

sps2(end+1) = subplot(4,3,5);  % avy vs time
plot(pt, (u(pnt,11))./t0);
hold on;

sps2(end+1) = subplot(4,3,6);  % avz vs time
plot(pt, -(u(pnt,12))./t0);
hold on;


sps2(end+1) = subplot(4,3,7);  % tqx vs time
plot(pt, Tq1(pnt));
hold on;
ylabel('Torque (AU)');

sps2(end+1) = subplot(4,3,8);  % fy vs time
plot(pt, Tq2(pnt));
hold on;

sps2(end+1) = subplot(4,3,9);  % fz vs time
plot(pt, Tq3(pnt));  
hold on;



sps2(end+1) = subplot(4,3,10);  % davx vs time
plot(dpt, diff(u(pnt,10))./dts);
hold on;
ylabel('d\omega/dt');
xlabel('Time (s)')


sps2(end+1) = subplot(4,3,11);  % davy vs time
plot(dpt, diff(u(pnt,11))./dts);
hold on;
xlabel('Time (s)')


sps2(end+1) = subplot(4,3,12);  % davz vs time
plot(dpt, diff(u(pnt,12))./dts);  %FIXME after recording this
hold on;
xlabel('Time (s)')

end
