% plotMesh

%% general plot settings
showMeshCutoff = 7;  % value of "radialRings" , above which we don't render the mesh edges. 
  % for meshes with fewer than this number of radial Rings, we draw each mesh edge with black lines
plotMeshOutlineZ0=0;
colormode = 3;  %0 use strain to color the mesh
                %1 use tex value to color the mesh
                %2 use temperature to color the mesh
                
                
if (videoMode > 1)
    return;
elseif videoMode > 0
    plotMeshFast_regenerateVideo_RG
%     plotMeshFast_v17
else
    
    %% Main 3D view plot:
    
    
    goodTs = ~t_broken;

    if isempty('straincolormap')
        load StrainColormap.mat
    end

    % Scientific color map from https://www.fabiocrameri.ch/colourmaps/
    if isempty('lajolla')
        load('lajolla.mat');
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
%         trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  cdata(goodTs),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
        trisurf(TRI(goodTs',:),n_x,n_y,zreliefmag*n_z,  ones(size(cdata(goodTs))),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
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
    end

    if colormode==0
        colormap(ha1, straincolormap );
        colorbar('southoutside')

    elseif colormode==1
        colormap(ha1,'parula');
        colorbar('southoutside')

    elseif colormode == 2
        colormap(flipud(lajolla));
        colorbar('southoutside')

    else
        colormap(ha1,'jet');
        colorbar('southoutside');
        set(gca,'CLim',[tscale_min tscale_max] );

    end
    hold on;
    
    
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

    if zreliefmag == 1
        axis equal
        set(gca,'ZLim',[min(zminforplot,comz-zbox) max(zmaxforplot,comz+zbox)]);
    end
    plot3(comx, comy, zreliefmag*(ceil((comz-zbox)/(2*zbox))*(2*xbox)), 'mx', 'MarkerSize',15);
    if (tt>t_ramp_delay)
        plot3([1 1].*((0)), [1 1].*((0)), zreliefmag*[comz-zbox comz+zbox], 'm');
    end
    
end
   
if 1


%%  Annotations
if (videoMode)
    plotAnoPos1 = [.02 .7 .2 .2];
    plotAnoPos2 = [.80 .78 .2 .2];
    plotAnoPos3 = [.02 .13 .4 .2];
%     plotAnoPos4 = [.02 .55 .2 .2];
    plotAnoPos5 = [.85 .2 .2 .2];
    plotAnoPos6 = [.02 .73 .5 .25]; 
else
    plotAnoPos1 = [.02 .78 .2 .2];
    plotAnoPos2 = [.43 .78 .2 .2];
    plotAnoPos3 = [.02 .16 .2 .2];
%     plotAnoPos4 = [.15 .78 .2 .2];
    plotAnoPos5 = [.43 .10  .2 .2];
    plotAnoPos6 = [.02 .67 .25 .2]; 
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
    ['z_{COM} = ' num2str( round(comz.*zmult,3) ) ' ' zunit] ...
    ['v_z = ' num2str( round(meanvz .* zvelmult,3) ) ' ' zvelunit ] ...
    ['a_z = ' num2str( round(accz ./ 9800) ) ' Gs' ] };

if (tt > 0)
    spinSplChar = '0';
else
    spinSplChar = '';
end

plotAnoStr2 = ...
   {['Study ' simdesc ' _{ }' ] ...
    [nameIntegrationMethod] ...
    ['Area = ' num2str(disparea .* areamult) '_{ }' areaunit] ...
    ['Mass = ' num2str(sum(t_m(t_notbroken))) ' g' ] ...
    ['Spin frequency = ' num2str(rot0) ' Hz ^{ }' ] };

if usetether
    tetheredstring = [num2str(tether_payload_mass) 'g tethered z=' num2str(tether_payload_distance) 'mm (' num2str(nte) ')' ];
else
    tetheredstring = 'No tethered payload';
end


plotAnoStr3 = { ...
    ['Radius: ' num2str(radActualx0mm) ' mm' ] ...
    ['dt = ' num2str( dt ) ' s'] ...
    ['Tensile margin = ' num2str(round(100*(1-max_strain_norm),2)) '% [' num2str(tensilefailures) ']' ] ...
    ['Minimum temperature = ' num2str(min_temp,5) ' (K)' ] ...
    ['Maximum temperature = ' num2str(max_temp,5) ' (K)' ] ...
    ['Average temperature = ' num2str(avg_temp,5) ' (K)' ] ...
    };


% plotAnoStr4 = { ...
%     ['PE = ' num2str( PE ) ' J^{ }'] ...
%     ['KE = ' num2str( KE ) ' J^{ }'] ...
%     ['KEb = ' num2str( KE_COM ) ' J^{ }'] ...
%     ['KEo = ' num2str( KEo ) ' J^{ }'] };


plotAnoStr5 = { ...
    ['I0 = ' num2str( I0*pwrRampVal/1e3 ) ' GW/m^2'] ...
    ['IR = ' num2str( round(Irad,3) ) ' mm^{ }'] ...
    ['Pin = ' num2str( dispttlpwr ) ' GW^{ }'] ...
    ['x0 = ' num2str(round(t0t_x0,3)) 'mm'] ...
    ['y0 = ' num2str(round(t0t_y0,3)) 'mm'] ...
    ['\theta_0 = ' num2str(round(rad2deg(t0t_t0_x),3)) ' deg'] ...
    ['\phi_0 = ' num2str(round(rad2deg(t0t_t0_y),3)) ' deg'] ...
    ['\psi_0 = ' num2str(round(rad2deg(t0t_t0_z),3)) ' deg'] ...
    };

mySailNorm = mean(t_norms(:,n_centerTris),2);
mySailNorm = mySailNorm ./ norm(mySailNorm);

if rigid_mode

    plotAnoStr6 = {  ...
        ['Pos = [' num2str(comx) ', ' num2str(comy) ', ' num2str(comz) ' ]' ] ...
        ['Velocity  = [ ' num2str(COMVel(1)) ', ' num2str(COMVel(2)) ', ' num2str(COMVel(3)) ' ]' ] ... 
        ['Ang vel   = [ ' num2str(angVel(1)) ', ' num2str(angVel(2)) ', ' num2str(angVel(3)) ' ]' ] ...
        ...['Ang vel Rigid = [ ' num2str(angVelRigid(1)) ', ' num2str(angVelRigid(2)) ', ' num2str(angVelRigid(3)) ' ]' ] ...
        ['Ang vel Rigid = [ ' num2str(rb_angVel(1,1)) ', ' num2str(rb_angVel(2,1)) ', ' num2str(rb_angVel(3,1)) ' ]' ] ...
        ['Up = [ ' num2str(mySailNorm(1)) ', ' num2str(mySailNorm(2)) ', ' num2str(mySailNorm(3)) ' ]' ] ...
        [ sprintf('%d',movieNum)  '   Mode = rigid'] };
        
else

      plotAnoStr6 = {  ...
        ['Pos = [' num2str(round(comx,3)) ', ' num2str(round(comy,3)) ', ' num2str(round(comz,3)) ' ] mm' ] ...
        ['Velocity = [ ' num2str(round(meanvx,3)) ', ' num2str(round(meanvy,3)) ', ' num2str(round(meanvz,3)) ' ] mm/s' ] };
end

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

% if exist('plotAno4','var') && (isvalid(plotAno4))
%     plotAno4.String = plotAnoStr4;
% else
%     plotAno4 = annotation('textbox', plotAnoPos4, 'String', plotAnoStr4, 'EdgeColor', 'none' );
% end

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

end
