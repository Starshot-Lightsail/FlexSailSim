% plotMesh

goodTs = ~t_broken;
zmaxforplot = max( n_z(t_na(goodTs)))+2;
zminforplot = min( n_z(t_na(goodTs)))+2;
xmaxforplot = max( n_x(t_na(goodTs)));
xminforplot = min( n_x(t_na(goodTs)));
ymaxforplot = max( n_y(t_na(goodTs)));
yminforplot = min( n_y(t_na(goodTs)));

hold off;


    cdata=t_Ipwr;
showMeshCutoff = 10;
%trimesh(TRI,n_x,n_y,n_z, n_rps);
if (radius > showMeshCutoff)
    trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
   % trimesh(TRI,n_x,n_y,n_z,  n_m,'FaceColor','interp','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    hold on
    lightangle(-120,60);
    view([35.8339   22.7847]);
    plotMeshOutline;
    plot3(n_x(1), n_y(1), n_z(1), 'k.');
    
    
else
    trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    hold on;
%  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
  %  vectors (can be confusing).
end
%debug
    quiver3( t_cx, t_cy, t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
%debug

    quiver3( t_cx, t_cy, t_cz, t_IMag.*t_IDir(1,:), t_IMag.*t_IDir(2,:), t_IMag.*t_IDir(3,:), 'r', 'LineWidth',0.5, 'ShowArrowHead', 'off');

plot3(n_x(xradidx), n_y(xradidx), n_z(xradidx), 'k.'); %one dot at the edge to help us track rotation...
%if  ( (sam<5) && (sam>2) )
%        plot3([1 1].*Ictrx, [1 1].*Ictry, [comz-zbox comz+zbox], 'k');
%end
%view(0,90)

%if any(t_broken)
%    plot3(t_cx(t_broken),t_cy(t_broken),t_cz(t_broken),'k.','MarkerSize',4)
%end
plot3(n_x(end), n_y(end), n_z(end), 'k.');
axis equal

        colormap jet
        colorbar('southoutside');
        climmax = max(climmax, max(cdata(goodTs)));
        set(gca,'CLim',[0 climmax] );

set(gca,'FontSize',14);



%set(gca,'CLim',[min(n_rps) max(n_rps)]);
hold on;
quiver3(t_cx,t_cy,t_cz,t_oth(1,:),t_oth(2,:),t_oth(3,:),1,'m'); 
%if (tt < 0)
%    quiver3(n_x,n_y,n_z,n_af(1,:),n_af(2,:),n_af(3,:),'r');
%end
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal


%set(gca,'XLim',[comx-xbox comx+xbox]);
set(gca,'XLim',[min(xminforplot,comx-xbox) max(xmaxforplot,comx+xbox)]);
%set(gca,'YLim',[comy comy+xbox]);
set(gca,'YLim',[min(yminforplot,comy-xbox) max(ymaxforplot,comy+xbox)]);
set(gca,'ZLim',[min(zminforplot,comz-zbox/1.5) max(zmaxforplot,comz+zbox/1.5)]);
% plot3(comx, comy, (ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
if (tt > t_ramp_delay) && (tt < t_ramp_delay + 2*t_ramp_dur + t_on_dur)
   plot3([1 1].*0, [1 1].*0,  [comz-zbox comz+zbox], 'm');
end


%lighting gouraud;



