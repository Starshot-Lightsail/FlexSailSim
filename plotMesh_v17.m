% plotMesh

%% general plot settings
showMeshCutoff = 7;  % value of "radialRings" , above which we don't render the mesh edges. 
  % for meshes with fewer than this number of radial Rings, we draw each mesh edge with black lines
plotMeshOutlineZ0=0;
colormode = 2;  %0 use strain to color the mesh
                %1 use tex value to color the mesh
                %2 use temperature to color the mesh
                
                
if (videoMode > 1)
    return;
elseif videoMode > 0
    plotMeshFast_v17
else
    
    %% Main 3D view plot:
    
    
    goodTs = ~t_broken;
    if (~exist('straincolormap'))
        load StrainColormap.mat
    end
    
    zmaxforplot = max( n_z(t_na(goodTs)));
    zminforplot = min( n_z(t_na(goodTs)));
    
    ha1 = subplot(2,4,[1,6]);
    hold off;
    
    %trimesh(TRI,n_x,n_y,n_z, n_m);%sqrt(t_mf(1,:).^2 + t_mf(2,:).^2 + t_mf(3,:).^2));
    if colormode==0
        cdata = (sqrt(t_a./t_a0)-1)./t_tensile_strain_limit;
    elseif colormode==1
        cdata = t_tex;
    else
        cdata = t_t;
    end
%    cdata = n_fx; % for shear wave plotting
    if (radius > showMeshCutoff)
        %trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'LineStyle','none');
        trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
%        trimesh(TRI(goodTs',:),n_x,n_y,n_z,  cdata); for shear wave plotting
        lightangle(-120,60);
        hold on
        plotMeshOutline;
       
        if any(t_broken) 
            plot3(t_cx(t_broken),t_cy(t_broken),zreliefmag*t_cz(t_broken),'k.', 'MarkerSize',4)
        else
            
        end
    else
        trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
        hold on;
    % if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end
    end
     %debug
            %quiver3( t_cx, t_cy, t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead','off');
%%            quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:), t_texn(2,:), t_texn(3,:), 'b', 'LineWidth',0.5, 'ShowArrowHead','off');
    if colormode==0
        colormap(ha1, straincolormap );
        colorbar('southoutside')
    % s_max = max(s_max, max( cdata - 1 ));
    % tscale_min = min(tscale_min,min(t_t));
    % tscale_max = max(tscale_max,max(t_t));
    %sscale_max = max( max(sscale_max, max_strain), -min_strain);
    %set(gca,'CLim',[1+min(-sscale_max,min_strain) 1+max(sscale_max, max_strain)]);
        set(gca,'CLim',[1 1]); % .* tensile_strain_limit
    elseif colormode==1
        colormap(ha1,'parula');
        colorbar('southoutside')
    else
        colormap(ha1,'jet');
        colorbar('southoutside');
        set(gca,'CLim',[tscale_min tscale_max] );
    end
    hold on;
    
    if pwrRampVal > 0  % light is on, plot optical force vectors
%%        quiver3(t_cx,t_cy,zreliefmag*t_cz,t_oth(1,:),t_oth(2,:),t_oth(3,:),1,'m');
        if RayTracingMode
            max_oth_scale = max(vecnorm(t_oth));
            max_oth_mr_scale = max(vecnorm(t_oth_mr));
            mr_quiver_ratio = max_oth_mr_scale / max_oth_scale;
%%            quiver3(t_cx,t_cy,zreliefmag*t_cz,t_oth_mr(1,:),t_oth_mr(2,:),t_oth_mr(3,:),mr_quiver_ratio,'r');
        end
    end
    if tt < 0  %if (max(max(n_af)) > 0)  if doing spin-up, show spin-up force vectors (easier to test for tt<0 rather than test for presence of n_af values)
%%        quiver3(n_x,n_y,zreliefmag*n_z,n_af(1,:),n_af(2,:),n_af(3,:),'r');
    end
    
    plot3(n_x(xradidx), n_y(xradidx), zreliefmag*n_z(xradidx), 'k.', 'MarkerSize',12); %one dot at the edge to help us track rotation...
    xlabel('x (mm)');
    ylabel('y (mm)');
    if zreliefmag~= 1
        zlabel(['z \times ' num2str(zreliefmag) ' (mm)']); 
    else
        zlabel('z (mm)');
    end

    if (colormode==0) title({'Strain (normalized to tensile limit)'}); end%, ['t=' num2str(tt) ' s' ] });
    if (colormode==2) title({'Temperature'}); end
    set(gca,'XLim',double([comx-xbox comx+xbox]));
    set(gca,'YLim',[comy-xbox comy+xbox]);
    %set(gca,'XLim',[comx comx+xbox]);
    %set(gca,'YLim',[comy comy+xbox]);
    %set(gca,'ZLim',[min(zminforplot,comz-zbox) max(zmaxforplot,comz+zbox/2)]);
    if zreliefmag == 1
        axis equal
        set(gca,'ZLim',[min(zminforplot,comz-zbox) max(zmaxforplot,comz+zbox)]);
    end
    plot3(comx, comy, zreliefmag*(ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
    if (tt>t_ramp_delay)
        plot3([1 1].*((0)), [1 1].*((0)), zreliefmag*[comz-zbox comz+zbox], 'm');
    end
    
   % view(0,90);
    
    
    %% Side plots setup
    if (sqrt((comx - ((0))).^2 + (comy - ((0))).^2) > xbox/2) || (abs(comy - 0) > xbox/2)
        sideY0 = comy;
        sideX0 = comx;
    else
        sideY0 = ((0));
        sideX0 = ((0));
    end
    
    
    
    %% Side plot 1:  Top down view, temperature
    ha2 = subplot(2,4,3) ;
    hold off;
    cdata = t_t;
    if (radius > showMeshCutoff )
        trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'LineStyle','none');
        hold on;
        %plot3([n_x(constrainedEdges(:,1)) n_x(constrainedEdges(end,2))],[n_y(constrainedEdges(:,1)) n_y(constrainedEdges(end,2))],[n_z(constrainedEdges(:,1)) n_z(constrainedEdges(end,2))],'k')
        plotMeshOutline;
        plot3(n_x(1), n_y(1), n_z(1), 'k.', 'MarkerSize',4);
    else
        trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs));
    end
    %if (nt < 2)
    colormap(ha2, 'jet');
    xlabel('x (mm)');
    ylabel('y (mm)');
    zlabel('z (mm)');
    colorbar;
    axis equal;

    %end
    hold on
    tscale_min = min(tscale_min,min(t_t));
    tscale_max = max(tscale_max,max(t_t));
    set(gca,'CLim',[tscale_min tscale_max]);
    view(0,90);
    title({'Temperature (\circK)'});%, ['t=' num2str(tt) ' s' ] });
    % set(gca,'ZLim',[comz-zbox comz+zbox]);
    %plot3(comx, comy, (ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
    %plot3([1 1].*((0)), [0 0], [comz-zbox comz+zbox], 'm');
    if (tt > t_ramp_delay) 
        %plot3( ((0)), ((0)), comz+zbox, 'mo');
        plot3([1 1].*((0)), [1 1].*((0)), [comz-zbox comz+zbox], 'm');
    end
    set(gca,'XLim',[sideX0-xbox-((0)) sideX0+xbox-((0))]);
    set(gca,'YLim',[sideY0-xbox-((0)) sideY0+xbox-((0))]);
    
    
    
    
    
    %% Side plot 2:  Top down view, absorbed power per unit area
    ha3 = subplot(2,4,4) ;
    hold off;
    if RayTracingMode
        cdata = t_abs+t_abs_mr;
    else
        cdata = t_abs;
    end
    cdata = cdata ./ t_a;
    if exist('max_abs_pwr_dens')
        max_abs_pwr_dens = max(max_abs_pwr_dens, max(cdata));
    else
        max_abs_pwr_dens = max(cdata);
    end
    max_abs_pwr_dens = max(max_abs_pwr_dens,1e-8);  %this has to be here, otherwise setting clim will fail at startup
    if (radius > showMeshCutoff)
        trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'LineStyle','none');
        hold on;
        %plot3([n_x(constrainedEdges(:,1)) n_x(constrainedEdges(end,2))],[n_y(constrainedEdges(:,1)) n_y(constrainedEdges(end,2))],[n_z(constrainedEdges(:,1)) n_z(constrainedEdges(end,2))],'k')
        plotMeshOutline;
        plot3(n_x(1), n_y(1), n_z(1), 'k.', 'MarkerSize',4);
    else
        trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs));
    end
   
    %if (nt < 2)
    colormap(ha3, 'jet');
    set(gca,'CLim',[0 max_abs_pwr_dens]);
    xlabel('x (mm)');
    ylabel('y (mm)');
    zlabel('z (mm)');
    axis equal
    hcb = colorbar;
    %end
    hold on
    %set(gca,'CLim',[0 I0]);
    view(0,90);
    title({'Absorbed power, W/mm^2'});%, ['t=' num2str(tt) ' s' ] });
    % set(gca,'ZLim',[comz-zbox comz+zbox]);
    %plot3(comx, comy, (ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
    %plot3([1 1].*((0)), [0 0], [comz-zbox comz+zbox], 'm');
    if (tt > t_ramp_delay)
        plot3( ((0)), ((0)), comz+zbox, 'mo');
    end
    set(gca,'XLim',[sideX0-xbox-((0)) sideX0+xbox-((0))]);
    set(gca,'YLim',[sideY0-xbox-((0)) sideY0+xbox-((0))]);
    
    
    
    
    
    
    
    %% Side plot 3:  Side view, mesh distortion
    ha4 = subplot(2,4,7) ;
    hold off;
    %cdata = sqrt( t_ix.^2 + t_iy.^2 + t_iz.^2);
    if (radius > 7)
        trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'FaceColor', [.5 .5 .5], 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.6,'SpecularStrength',0.3,'DiffuseStrength',0.6);
        lightangle(60,70);
    else
        trimesh(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,'EdgeColor',[0 0 0]);
    end
    %if (nt < 2)
    % colormap(ha3, 'jet');
    xlabel('x (mm)');
    ylabel('y (mm)');
    zlabel('z (mm)');
    axis equal
    %colorbar;
    %end
    hold on
    %set(gca,'CLim',[0 I0]);
    view(0,0);
    title({'XZ'});%, ['t=' num2str(tt) ' s' ] });
    
    if zreliefmag == 1
        set(gca,'ZLim',[min(zminforplot,comz-zbox) max(zmaxforplot,comz+zbox)]);
    end
    set(gca,'XLim',[sideX0-xbox-abs(((0))/2) sideX0+xbox+abs(((0))/2)]);
    if (tt > t_ramp_delay)
        plot3([1 1].*((0)), [1 1].*((0)), get(gca,'Zlim'), 'm');
    end
    par=get(gca,'PlotBoxAspectRatio');
    
    
    
    
    
    
    
    % ha5 = subplot(2,4,8) ;
    % hold off;
    % %cdata = sqrt( t_ix.^2 + t_iy.^2 + t_iz.^2);
    % trimesh(TRI(goodTs',:),n_x,n_y,n_z,'EdgeColor',[0 0 0]);
    % %if (nt < 2)
    %    % colormap(ha3, 'jet');
    %     xlabel('x (mm)');
    %     ylabel('y (mm)');
    %     zlabel('z (mm)');
    %     axis equal
    %     %colorbar;
    % %end
    % hold on
    % %set(gca,'CLim',[0 I0]);
    % view(90,0);
    % title({'YZ', ['t=' num2str(tt) ' s' ] });
    % set(gca,'ZLim',[comz-zbox comz+zbox]);
    % set(gca,'YLim',[sideY0-xbox sideY0+xbox]);
    % plot3([1 1].*((0)), [0 0], get(gca,'Zlim'), 'm');
    
    
    
    %graphtspan = 100000;
    
%     % side plot 4:  X/Y position vs time plot
%     if (nto > 1)
%         ha5 = subplot(2,4,8) ;
%         hold off;
%         hl1 = plot(graphts(1:nto),(comxs(1:nto)-((0))),'LineWidth',1);
%         hold on
%         hl2 = plot(graphts(1:nto),comys(1:nto),'LineWidth',1);
%         hl3 = plot(graphts(1:nto),vertxs(1:nto)-((0)),':','LineWidth',1);
%         hl4 = plot(graphts(1:nto),vertys(1:nto),':','LineWidth',1);
%         hl3.Color = hl1.Color;
%         hl4.Color = hl2.Color;
%         ylabel('(mm)');    xlabel('t (s)');
%         lgnd = legend({ 'X (COM)' , 'Y (COM)', 'X (vertex)', 'Y (vertex)' });
%         lgnd.NumColumns = 2;  %  R2018 only
%         title('Position rel. to beam center');
%         set(gca,'PlotBoxAspectRatio',[par(1) par(3) 1]);
%         alold = get(gca,'YLim');
%         set(gca,'YLim',[alold(1) alold(1)+1.5*(alold(2)-alold(1))]);
%         
%     end
    


%   %% side plot 4:  stress amplitude cross section
%   if (nto > 1)
%       ha5 = subplot(2,4,8) ;
%          hold off;
%          hl1 = plot(n_y(mycs),n_mf(2,mycs),'LineWidth',1);
% %         hold on
% %         hl2 = plot(graphts(1:nto),comys(1:nto),'LineWidth',1);
% %         hl3 = plot(graphts(1:nto),vertxs(1:nto)-((0)),':','LineWidth',1);
% %         hl4 = plot(graphts(1:nto),vertys(1:nto),':','LineWidth',1);
% %         hl3.Color = hl1.Color;
% %         hl4.Color = hl2.Color;
%         ylabel('(N)');    xlabel('y (mm)');
% %         lgnd = legend({ 'X (COM)' , 'Y (COM)', 'X (vertex)', 'Y (vertex)' });
% %         lgnd.NumColumns = 2;  %  R2018 only
%         title('Cross section along Y');
%          set(gca,'PlotBoxAspectRatio',[par(1) par(3) 1]);
% %         alold = get(gca,'YLim');
% %         set(gca,'YLim',[alold(1) alold(1)+1.5*(alold(2)-alold(1))]); 
%   end
  
  
  %% alt plot 4:  Plot radii of outer edge elements to look for elongation...
  if (nto > 1)
      ha5 = subplot(2,4,8) ;
         hold off;
         hl1 = plot(ringradii,'LineWidth',1);
%         hold on
%         hl2 = plot(graphts(1:nto),comys(1:nto),'LineWidth',1);
%         hl3 = plot(graphts(1:nto),vertxs(1:nto)-((0)),':','LineWidth',1);
%         hl4 = plot(graphts(1:nto),vertys(1:nto),':','LineWidth',1);
%         hl3.Color = hl1.Color;
%         hl4.Color = hl2.Color;
        ylabel('(mm)');    xlabel('index');
%         lgnd = legend({ 'X (COM)' , 'Y (COM)', 'X (vertex)', 'Y (vertex)' });
%         lgnd.NumColumns = 2;  %  R2018 only
        title('Hoop radii');
         set(ha5,'PlotBoxAspectRatio',[par(1) par(3) 1]);
%         alold = get(gca,'YLim');
%         set(gca,'YLim',[alold(1) alold(1)+1.5*(alold(2)-alold(1))]); 
  end
  
end



%%  Annotations
if (videoMode)
    plotAnoPos1 = [.02 .78 .2 .2];
    plotAnoPos2 = [.80 .78 .2 .2];
    plotAnoPos3 = [.02 .09 .2 .2];
    plotAnoPos4 = [.20 .78 .2 .2];
    plotAnoPos5 = [.85 .07 .2 .2];
    plotAnoPos6 = [.02 .65 .27 .25]; 
else
    plotAnoPos1 = [.02 .78 .2 .2];
    plotAnoPos2 = [.43 .78 .2 .2];
    plotAnoPos3 = [.02 .16 .2 .2];
    plotAnoPos4 = [.15 .78 .2 .2];
    plotAnoPos5 = [.43 .10  .2 .2];
    plotAnoPos6 = [.02 .65 .27 .2]; 
end



zunit = 'mm';
zmult = 1;
if (abs(comz) > 1e6 )
    zunit = 'km';
    zmult = 1e-6;
elseif (abs(comz) > 1000 )
    zunit = 'm';
    zmult = 1e-3;
end
zvelunit = 'mm/s';
zvelmult = 1;
if (abs(meanvz) > 1e6 )
    zvelunit = 'km/s';
    zvelmult = 1e-6;
elseif (abs(meanvz) > 1000 )
    zvelunit = 'm/s';
    zvelmult = 1e-3;
end
areaunit = 'mm^2';
areamult = 1;
if (totalarea >= 1e5)
    areaunit = 'm^2';
    areamult = 1e-6;
end

if (nto<1)
    disparea = sum(t_a0);
    dispareaxy = sum(t_a0xy);
    dispttlpwr = 0;
else
    disparea = areas(nto);
    dispareaxy = areasxy(nto);
    dispttlpwr = ttlinpwr(nto)./1e9;
end

plotAnoStr1 = ...
    {['t = ' num2str(tt) ' s_{ }' ]   ...
    ['z_{COM} = ' num2str(comz.*zmult) ' ' zunit] ...
    ['V_z = ' num2str( meanvz .* zvelmult) ' ' zvelunit ] ...
    ['A_z = ' num2str( accz ./ 9800 ) ' Gs' ] };

if (tt > 0)
    spinSplChar = '0';
else
    spinSplChar = '';
end
plotAnoStr2 = ...
   {['Study ' simdesc ' _{ }' ] ...
    ['Area = ' num2str(disparea .* areamult) ' _{ }' areaunit] ...
    ['Area (proj.) = ' num2str(dispareaxy .* areamult) ' _{ }' areaunit] ...
    ['Mass = ' num2str(sum(t_m(t_notbroken))) ' g' ] ...
    ['Spin0' spinSplChar ' = ' num2str(rot0) ' Hz ^{ }' ] };

if usetether
    tetheredstring = [num2str(tether_payload_mass) 'g tethered z=' num2str(tether_payload_distance) 'mm (' num2str(nte) ')' ];
else
    tetheredstring = 'No tethered payload';
end

plotAnoStr3 = { ...
    ['Mass: ' sprintf('%6g ', sum(t_m(t_notbroken)))      'g + ' sprintf('%6g', totalmass_nodes-totalmass_sail) 'g'  ] ...
    ['Area:   ' sprintf('%6g ', sum(t_a(t_notbroken)) /1e6) 'm2' ] ...
    ['Tensile margin = ' num2str(100*(1-max(e_s(e_notbroken)./e_tensile_strain_limit(e_notbroken)))) '% [' num2str(tensilefailures) ']' ] ...
    ['Thermal margin = ' num2str(-max(e_t(e_notbroken) - e_failtemp(e_notbroken))) ' (K) ['  num2str(thermalfailures) ']' ] ...
    [zMode ' AR= ' num2str(zAspectRatio) '  ' num2str(radius) ' rings' ] ...
         tetheredstring ...
    ['Nom rad: ' num2str(radActualx0mm) ' mm' ] ...
    ['dT = ' num2str( dt ) ' s'] ...
    ['nt/nto/nf = ' num2str(nt) ' / ' num2str(nto) ' / ' num2str(nf) ] };


plotAnoStr4 = { ...
    ['PE = ' num2str( PE ) ' J^{ }'] ...
    ['KE = ' num2str( KE ) ' J^{ }'] ...
    ['KEb = ' num2str( KE_COM ) ' J^{ }'] ...
    ['KEo = ' num2str( KEo ) ' J^{ }'] };


plotAnoStr5 = { ...
    ['I0 = ' num2str( I0*pwrRampVal/1e3 ) ' GW/m^2'] ...
    ['IR = ' num2str( Irad ) ' mm^{ }'] ...
    ['Ix = ' num2str(((0))/(radiusmm)) ' (' num2str(((0))) ' mm) ^{ }' ] ...
    ['Iy = ' num2str(((0))/(radiusmm)) ' (' num2str(((0))) ' mm) ^{ }' ] ...
    ['Pin = ' num2str( dispttlpwr ) ' GW^{ }'] };

mySailNorm = mean(t_norms(:,t_centerTris),2);
mySailNorm = mySailNorm ./ norm(mySailNorm);

    plotAnoStr6 = {  ...
        ['Pos = [' num2str(comx) ', ' num2str(comy) ', ' num2str(comz) ' ]' ] ...
        ['Velocity = [ ' num2str(meanvx) ', ' num2str(meanvy) ', ' num2str(meanvz) ' ]' ] ...
        ['Ang vel        = [ ' num2str(angVel(1)) ', ' num2str(angVel(2)) ', ' num2str(angVel(3)) ' ]' ] ...
        ['Up = [ ' num2str(mySailNorm(1)) ', ' num2str(mySailNorm(2)) ', ' num2str(mySailNorm(3)) ' ]' ] ...
        [sprintf('%d',movieNum) ]};

if  exist('plotAno1','var') && (isvalid(plotAno1))
    plotAno1.String = plotAnoStr1;
else
    plotAno1 = annotation('textbox', plotAnoPos1, 'String', plotAnoStr1, 'EdgeColor', 'none' );
end

if  exist('plotAno2','var') && (isvalid(plotAno2))
    plotAno2.String = plotAnoStr2;
else
    plotAno2 = annotation('textbox', plotAnoPos2, 'String', plotAnoStr2, 'EdgeColor', 'none' );
end

if exist('plotAno3','var') && (isvalid(plotAno3))
    plotAno3.String = plotAnoStr3;
else
    plotAno3 = annotation('textbox', plotAnoPos3, 'String', plotAnoStr3, 'EdgeColor', 'none' );
end

if exist('plotAno4','var') && (isvalid(plotAno4))
    plotAno4.String = plotAnoStr4;
else
    plotAno4 = annotation('textbox', plotAnoPos4, 'String', plotAnoStr4, 'EdgeColor', 'none' );
end

if exist('plotAno5','var') && (isvalid(plotAno5))
    plotAno5.String = plotAnoStr5;
else
    plotAno5 = annotation('textbox', plotAnoPos5, 'String', plotAnoStr5, 'EdgeColor', 'none' );
end

if exist('plotAno6','var') && (isvalid(plotAno6))
    plotAno6.String = plotAnoStr6;
else
    plotAno6 = annotation('textbox', plotAnoPos6, 'String', plotAnoStr6, 'EdgeColor', 'none' );
end

%['Beam radius = ' num2str((Irad)) ' mm^{ }' ] ...
%    ['  offset = ' num2str(((0))/(radiusmm)) ' (' num2str(((0))) ' mm) ^{ }' ] },

