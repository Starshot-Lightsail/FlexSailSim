% facilitates comparisons to Ramon's rigid body simulator


% to be run in the workspace following a simulation.

t_cutoff = 1;

nt1 = floor(interp1(plot_time,1:length(plot_time),0))+1;

nt2 = ceil(interp1(plot_time, 1:length(plot_time),t_cutoff));
if isnan(nt2)
    nt2 = length(plot_time);
end

%close all
%close(figure(1));
%close(figure(2));

pnt = nt1:nt2;
pt = plot_time(nt1:nt2);

figure(1) % position, velocity, force, accel=dv/dt
sps = [];

sps(end+1) = subplot(4,3,1);  % x vs time
plot(pt, comxs(pnt));
hold on;
title('X');
ylabel('Position (mm)');

sps(end+1) = subplot(4,3,2);  % y vs time
plot(pt, comys(pnt));
hold on;
title('Y');

sps(end+1) = subplot(4,3,3);  % z vs time
plot(pt, comzs(pnt));
hold on;
title('Z');



sps(end+1) = subplot(4,3,4);  % vx vs time
plot(pt, velxs(pnt));
hold on;
ylabel('Vel (mm/s)');

sps(end+1) = subplot(4,3,5);  % vy vs time
plot(pt, velys(pnt));
hold on;

sps(end+1) = subplot(4,3,6);  % vz vs time
plot(pt, velzs(pnt));
hold on;


sps(end+1) = subplot(4,3,7);  % fx vs time
plot(pt, ofxs(pnt));
hold on;
ylabel('Force (N)');

sps(end+1) = subplot(4,3,8);  % fy vs time
plot(pt, ofys(pnt));
hold on;

sps(end+1) = subplot(4,3,9);  % fz vs time
plot(pt, ofzs(pnt)); 
hold on;


dpt = pt(1:end-1);
dts = diff(pt);

sps(end+1) = subplot(4,3,10);  % ax vs time
plot(dpt, diff(velxs(pnt))./dts);
hold on;
ylabel('dv/dt');
xlabel('Time (s)')

sps(end+1) = subplot(4,3,11);  % ay vs time
plot(dpt, diff(velys(pnt))./dts);
hold on;
xlabel('Time (s)')

sps(end+1) = subplot(4,3,12);  % az vs time
plot(dpt, diff(velzs(pnt))./dts);  
hold on;
xlabel('Time (s)')











figure(2) % ang position, ang velocity, torque, diff-ang-vel
sps2 = [];

sps2(end+1) = subplot(4,3,1);  % angle-x vs time
plot(pt, rad2deg(atan2(sailNorms(2,pnt), sailNorms(3,pnt))));
hold on;
title('X');
ylabel('Angle (\circ)');

sps2(end+1) = subplot(4,3,2);  % angle-y vs time
plot(pt, rad2deg(atan2(sailNorms(1,pnt), sailNorms(3,pnt))));
hold on;
title('Y');

sps2(end+1) = subplot(4,3,3);  % angle-z vs time
plot(pt, rad2deg(spinAngs(pnt)));
hold on;
legend({'approx'});
title('Z');



sps2(end+1) = subplot(4,3,4);  % avx vs time
plot(pt, angVels(1,pnt));
hold on;
ylabel('AngVel (rad/s)');

sps2(end+1) = subplot(4,3,5);  % avy vs time
plot(pt, angVels(2,pnt));
hold on;

sps2(end+1) = subplot(4,3,6);  % avz vs time
plot(pt, angVels(3,pnt));
hold on;


sps2(end+1) = subplot(4,3,7);  % tqx vs time
plot(pt, torques(1,pnt));
hold on;
ylabel('Torque (AU)');

sps2(end+1) = subplot(4,3,8);  % fy vs time
plot(pt, torques(2,pnt));
hold on;

sps2(end+1) = subplot(4,3,9);  % fz vs time
plot(pt, torques(3,pnt));  %FIXME after recording this
hold on;



sps2(end+1) = subplot(4,3,10);  % davx vs time
plot(dpt, diff(angVels(1,pnt))./dts);
hold on;
ylabel('d\omega/dt');
xlabel('Time (s)')

sps2(end+1) = subplot(4,3,11);  % davy vs time
plot(dpt, diff(angVels(2,pnt))./dts);
hold on;
xlabel('Time (s)')

sps2(end+1) = subplot(4,3,12);  % davz vs time
plot(dpt, diff(angVels(3,pnt))./dts);  %FIXME after recording this
hold on;
xlabel('Time (s)')



