
plot_positionxs = comxs(nto_startAccel:nto_firstBroken-1)-Ictrx;
plot_positionys = comys(nto_startAccel:nto_firstBroken-1)-Ictry;
plot_positionrs = sqrt( plot_positionxs.^2 + plot_positionys.^2);

plot_tiltxs = rad2deg(atan2(sailNorms(2,nto_startAccel:nto_firstBroken-1),sailNorms(3,nto_startAccel:nto_firstBroken-1)));
plot_tiltys = -rad2deg(atan2(sailNorms(1,nto_startAccel:nto_firstBroken-1),sailNorms(3,nto_startAccel:nto_firstBroken-1)));
plot_tiltrs = rad2deg(atan2( sqrt(sailNorms(1,nto_startAccel:nto_firstBroken-1).^2 + sailNorms(2,nto_startAccel:nto_firstBroken-1).^2), sailNorms(3,nto_startAccel:nto_firstBroken-1)));

%plot_dxs = 0;
%plot_dxs(2:length(plot_positionxs)-1) = diff(plot_positionxs,2);
plot_dxs = diff(plot_positionxs,2);
%plot_dxs(end+1) = 0;

%plot_dys = 0;
%plot_dys(2:length(plot_positionys)-1) = diff(plot_positionys,2);
plot_dys = diff(plot_positionys,2);
%plot_dys(end+1) = 0;

%plot_drs = 0;
%plot_drs(2:length(plot_positionrs)-1) = diff(plot_positionrs,2);
plot_drs = diff(plot_positionrs,2);
%plot_drs(end+1) = 0;

%plot_dtiltxs = 0;
%plot_dtiltxs(2:length(plot_tiltxs)-1) = diff(plot_tiltxs,2);
plot_dtiltxs = diff(plot_tiltxs,2);
%plot_dtiltxs(end+1)=0;

%plot_dtiltys = 0;
%plot_dtiltys(2:length(plot_tiltys)-1) = diff(plot_tiltys,2);
plot_dtiltys = diff(plot_tiltys,2);
%plot_dtiltys(end+1)=0;

%plot_dtiltrs = 0;
%plot_dtiltrs(2:length(plot_tiltrs)-1) = diff(plot_tiltrs,2);
plot_dtiltrs = diff(plot_tiltrs,2);
%plot_dtiltrs(end+1)=0;


            
 figure('Position',[500        200        1200         1600]);
            set(gcf,'color','w');
            subplot(3,2,1)
            hline = plot(plot_positionxs, ofxs(nto_startAccel:nto_firstBroken-1),'DisplayName','X forces');
            hold on
            hline = plot(plot_positionys, ofys(nto_startAccel:nto_firstBroken-1),'DisplayName','Y forces');
%             if (nto_firstBroken < nto)
%                 hold on;
%                 plot3(comxs(nto_firstBroken:nto),comys(nto_firstBroken:nto),plot_time(nto_firstBroken:nto),'Color',get(hline,'Color'),'LineStyle',':');
%                 hline2=plot3(comxs(nto_firstBroken),comys(nto_firstBroken),plot_time(nto_firstBroken),'*');
%                 set(hline2,'Color',get(hline,'Color'));
%             end
            legend            ;
            xlabel('Displacement to beam center (mm)');
            ylabel('Force');
            
            subplot(3,2,2)
            hline = plot(plot_tiltxs, tqxs(nto_startAccel:nto_firstBroken-1),'DisplayName','X Torques');
            hold on
            hline = plot(plot_tiltys, tqys(nto_startAccel:nto_firstBroken-1),'DisplayName','Y torques');
            legend;
            xlabel('Tilt about axis (\circ)')
            ylabel('Torque');
            
            
            
            
            subplot(3,2,3)
            plot(plot_positionxs(2:end-1), plot_dxs, 'DisplayName','X');
            hold on;
            plot(plot_positionys(2:end-1), plot_dys, 'DisplayName','Y');
            xlabel('Displacement to beam center (mm)');
            ylabel('Accel (A.U.)');
            legend;
            
            
                   
            
            subplot(3,2,4)
            plot(plot_tiltxs(2:end-1), plot_dtiltxs, 'DisplayName','X');
            hold on;
            plot(plot_tiltys(2:end-1), plot_dtiltys, 'DisplayName','Y');
            xlabel('Tilt about axis (\circ)')
            %ylabel('$\frac{\partial^2 tilt}{/difft^2}$','Interpreter' ,'latex');
            ylabel('\Delta^2 Tilt');
            
            myymax = max( prctile(plot_dtiltxs,99), prctile(plot_dtiltys,99) );
            myymin = min( prctile(plot_dtiltxs,1), prctile(plot_dtiltys,1) );
            
            set(gca,'YLim', 1.5*[myymin myymax] )
            legend;
            
            
            
            subplot(3,2,5)
            plot(plot_positionrs(2:end-1), plot_drs, 'DisplayName','R');
            xlabel('Radial displacement to beam center (mm)');
            ylabel('\Delta^2 Displacement (A.U.)');
                        
                              
            
            subplot(3,2,6)
            plot(plot_tiltrs(2:end-1), plot_dtiltrs, 'DisplayName','R');
            xlabel('Radial tilt about axis (\circ)')
            %ylabel('$\frac{\partial^2 tilt}{/difft^2}$','Interpreter' ,'latex');
            ylabel('\Delta^2 Tilt (A.U.)');
            
            myymax = prctile(plot_dtiltrs,99);
            myymin = prctile(plot_dtiltxs,1);
            
            set(gca,'YLim', 1.5*[myymin myymax] )