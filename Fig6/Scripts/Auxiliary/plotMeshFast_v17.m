% plotMesh

goodTs = ~t_broken;
zmaxforplot = max( n_z(t_na(goodTs)))+2;
zminforplot = min( n_z(t_na(goodTs)))-2;
xmaxforplot = max( n_x(t_na(goodTs)));
xminforplot = min( n_x(t_na(goodTs)));
ymaxforplot = max( n_y(t_na(goodTs)));
yminforplot = min( n_y(t_na(goodTs)));
tscale_min = min(tscale_min,min(t_t));
tscale_max = max(tscale_max,max(t_t));

hold off;


    if (~exist('straincolormap'))
        load StrainColormap.mat
    end

%trimesh(TRI,n_x,n_y,n_z, n_m);%sqrt(t_mf(1,:).^2 + t_mf(2,:).^2 + t_mf(3,:).^2));
%cdata = t_a./t_a0-1;
if colormode==0
        cdata = (sqrt(t_a./t_a0)-1)./t_tensile_strain_limit;
elseif colormode==1
        cdata = t_tex;
else
    cdata=t_t;
end
%trimesh(TRI,n_x,n_y,n_z, n_rps);
if (radialRings > showMeshCutoff)
    trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
   % trimesh(TRI,n_x,n_y,n_z,  n_m,'FaceColor','interp','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    hold on
    lightangle(-120,60);
    view([35.8339   22.7847]);
    plotMeshOutline;
    plot3(n_x(1), n_y(1), zreliefmag*n_z(1), 'k.', 'MarkerSize',10);
    
    
else
    trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    hold on;
%%    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
%  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
  %  vectors (can be confusing).
end

if usetether
    for ti = 1:nte
        plot3([n_x(end) n_x(tethered_node_idxs(ti))], [n_y(end) n_y(tethered_node_idxs(ti))], [n_z(end) n_z(tethered_node_idxs(ti))], '--bo', 'MarkerSize', 10 ) 
    end
end


plot3(n_x(xradidx), n_y(xradidx), zreliefmag*n_z(xradidx), 'k.', 'MarkerSize',12); %one dot at the edge to help us track rotation...

%view(0,90)

if any(t_broken)
    plot3(t_cx(t_broken),t_cy(t_broken),zreliefmag*t_cz(t_broken),'k.','MarkerSize',4)
end
%plot3(n_x(end), n_y(end), zreliefmag*n_z(end), 'k.');
if zreliefmag == 1
    axis equal
end
if colormode==0
    colormap( straincolormap );
    colorbar('southoutside')
    set(gca,'CLim',[-1 1] ); %[-tensile_strain_limit tensile_strain_limit]);
elseif colormode==1
    colormap(gca,'parula');
    colorbar('southoutside');
elseif colormode==2
    colormap(flipud(lajolla));
    colorbar('southoutside')
else
    colormap jet
    colorbar('southoutside');
    set(gca,'CLim',[tscale_min tscale_max] );
end
set(gca,'FontSize',14);



%set(gca,'CLim',[min(n_rps) max(n_rps)]);
hold on;


if (tt < 0)
%%    quiver3(n_x,n_y,zreliefmag*n_z,n_af(1,:),n_af(2,:),n_af(3,:),'r');
else
%%    quiver3(t_cx,t_cy,zreliefmag*t_cz,t_oth(1,:),t_oth(2,:),t_oth(3,:),1,'m'); 
    if RayTracingMode 
       max_oth_scale = max(vecnorm(t_oth));
       max_oth_mr_scale = max(vecnorm(t_oth_mr));
       mr_quiver_ratio = max_oth_mr_scale / max_oth_scale;
%%       quiver3(t_cx,t_cy,zreliefmag*t_cz,t_oth_mr(1,:),t_oth_mr(2,:),t_oth_mr(3,:),mr_quiver_ratio,'r');
    end
end
xlabel('x (mm)');
ylabel('y (mm)');
if zreliefmag~= 1
        zlabel(['z \times ' num2str(zreliefmag) ' (mm)']); 
else
        zlabel('z (mm)');
end

if (colormode==0) title({'Strain (normalized to tensile limit)'}); end
if (colormode==2) title({'Temperature'}); end%, ['t=' num2str(tt) ' s' ] });
%set(gca,'XLim',[comx-xbox comx+xbox]);
set(gca,'XLim',[min(xminforplot,comx-xbox) max(xmaxforplot,comx+xbox)]);
%set(gca,'YLim',[comy comy+xbox]);
set(gca,'YLim',[min(yminforplot,comy-xbox) max(ymaxforplot,comy+xbox)]);
if zreliefmag == 1
    %axis equal
    set(gca,'ZLim',[min(zminforplot,comz-zbox) max(zmaxforplot,comz+zbox)]);
    plot3(comx, comy, zreliefmag*(ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
else
   % pbaspect([1 1 zreliefmag])
    set(gca,'ZLim',zreliefmag*n_z(1) + (zreliefmag .* max( abs(zmaxforplot-n_z(1)), abs(n_z(1)-zminforplot)) .* [-1 1]) );
end

% plot3(comx, comy, (ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
if (tt > t_ramp_delay) && (tt < t_ramp_delay + 2*t_ramp_dur + t_on_dur)
   beamdispheight = get(gca,'ZLim');
    plot3([1 1].*((0)), [1 1].*((0)),  beamdispheight, 'm');
end


%lighting gouraud;



