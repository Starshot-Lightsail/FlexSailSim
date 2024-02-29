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

if colormode==0
        cdata = (sqrt(t_a./t_a0)-1)./t_tensile_strain_limit;
elseif colormode==1
        cdata = t_tex;
else
    cdata=t_t;
end

if (radialRings > showMeshCutoff)
    
%     max_distance_centroids_from_flat_shape = max(abs(distance_centroids_from_flat_shape));
    total_min_caxis_disp_nodes = min(distance_nodes_from_flat_shape);

    total_max_caxis_disp_nodes = max(distance_nodes_from_flat_shape);

%     trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
%     trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  ones(size(cdata(goodTs))),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    patch('Faces',TRI,'Vertices',[n_x n_y n_z],'FaceVertexCData', distance_nodes_from_flat_shape,'FaceColor','interp','LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
%     hold off;
%     lightangle(-120,60);
    view([35.8339   22.7847]);
    
%     max_limit = max(abs(total_min_caxis_disp_nodes),abs(total_max_caxis_disp_nodes));
    max_limit = max(max_distance_nodes_from_flat_shape);


    plotMeshOutline;
        hold on

%     plot3(n_x(1), n_y(1), zreliefmag*n_z(1), 'k.', 'MarkerSize',10);
    
%     delete(p)
    
else

%     trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  distance_nodes_from_flat_shape','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);

    hold on;
end


plot3(n_x(xradidx), n_y(xradidx), zreliefmag*n_z(xradidx), 'k.', 'MarkerSize',12); %one dot at the edge to help us track rotation...

if any(t_broken)
    plot3(t_cx(t_broken),t_cy(t_broken),zreliefmag*t_cz(t_broken),'k.','MarkerSize',4)
end

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
elseif colormode == 3
    colormap(vik)
    cb = colorbar('southoutside');
    colorbarPos = cb.Position;
%     cb.Position = [colorbarPos(1) 0.01 colorbarPos(3) colorbarPos(4)];
    caxis([-max_limit max_limit])

else
    colormap jet
    colorbar('southoutside');
    set(gca,'CLim',[tscale_min tscale_max] );
end
set(gca,'FontSize',14);



%set(gca,'CLim',[min(n_rps) max(n_rps)]);
% hold on;


xlabel('x (mm)');
ylabel('y (mm)');
if zreliefmag~= 1
        zlabel(['z \times ' num2str(zreliefmag) ' (mm)']); 
else
        zlabel('z (mm)');
end

if (colormode==0) title({'Strain (normalized to tensile limit)'}); end
if (colormode==2) title({'Temperature'}); end%, ['t=' num2str(tt) ' s' ] });
if (colormode==3) title({'Displacement from flat shape (mm)'}); end%, ['t=' num2str(tt) ' s' ] });

%set(gca,'XLim',[comx-xbox comx+xbox]);
set(gca,'XLim',[min(xminforplot,comx-xbox) max(xmaxforplot,comx+xbox)]);
%set(gca,'YLim',[comy comy+xbox]);
set(gca,'YLim',[min(yminforplot,comy-xbox) max(ymaxforplot,comy+xbox)]);
if zreliefmag == 1
    %axis equal
    set(gca,'ZLim',[min(zminforplot,comz-zbox) max(zmaxforplot,comz+zbox)]);
%     plot3(comx, comy, zreliefmag*(ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
else
   % pbaspect([1 1 zreliefmag])
    set(gca,'ZLim',zreliefmag*n_z(1) + (zreliefmag .* max( abs(zmaxforplot-n_z(1)), abs(n_z(1)-zminforplot)) .* [-1 1]) );
end

% plot3(comx, comy, (ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
if (tt > t_ramp_delay) && (tt < t_ramp_delay + 2*t_ramp_dur + t_on_dur)
   beamdispheight = get(gca,'ZLim');
   plot3([1 1].*((0)), [1 1].*((0)),  beamdispheight, 'm');
end


% hold off


