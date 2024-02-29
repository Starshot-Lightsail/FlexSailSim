figure('Position',[48        500        1500         1500], 'Color','w')
set(gcf,'Renderer','painters');
load colormaps.mat

% 1600 1280 for other view

fs=14;  % font size
dolight=1;  % lighting settings, set to 0 to do non-shaded display
la1=-120;  % light angle
la2=30;    % light angle
intens=0.5;  % specular intensity
as=1;
ds=0.6;
ss=0.3; 
va1=-37.5000;  % view angle az
va2=30;   % view angle el

%alt view 1:
va1 = -190;
va2 = 45;
la1=-270;
la2=0;


load colormaps.mat

myidxs = 1:savemesh_num;  % which frames to plot (can be a single frame)
%myidxs=1000;
if length(myidxs)>1
    VW = VideoWriter([ filebasename '_RHT' ],'MPEG-4');  %alternate format for Ramon:  'Archival'
    VW.Quality = 100;
    open(VW);
end

frameno = 0;
for myidx=myidxs
    if graphts(savemesh_ntos(myidx)) < 0
        continue
    end
    frameno=frameno+1;
    myta = savemesht_a(:,myidx)';
    mynx = savemeshn_x(:,myidx)';
    myny = savemeshn_y(:,myidx)';
    mynz = savemeshn_z(:,myidx)';
    
    
    
    subplot(2,3,1);
    %trisurf(TRI,mynx,myny,mynz, savemesht_t(:,myidx)', 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds); 
    asdf=trisurf(TRI,mynx,myny,mynz,savemeshn_t(:,myidx)', 'FaceColor', 'interp', 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds); 
    set(gca,'Colormap',jet);
    axis equal
    axis off;
    hc = colorbar('southoutside'); 
    hc.Label.String = '\circK';
    set(gca,'fontsize',fs)
    view(va1, va2);
    title('Temperature')
    set(hc,'fontsize',fs);
    if dolight mylt=lightangle(la1,la2); mylt.Color=[intens intens intens]; end
    hold on; plot3([mynx(constrainedEdges(:,1)) mynx(constrainedEdges(end,2))],[myny(constrainedEdges(:,1)) myny(constrainedEdges(end,2))],zreliefmag*([mynz(constrainedEdges(:,1)) mynz(constrainedEdges(end,2))]),'k'); hold off 
    myclim=get(gca,'CLim');
    if frameno==1 
        ax1clim = myclim;
        axxspan = diff(get(gca,'XLim'))./1.98;
        axyspan = diff(get(gca,'YLim'))./1.98;
        axzspan = diff(get(gca,'ZLim'))./1.2;
    end
    set(gca,'XLim', mean(get(gca,'XLim')) + [-axxspan axxspan]);
    set(gca,'YLim', mean(get(gca,'YLim')) + [-axyspan axyspan]);
    set(gca,'ZLim', mean(get(gca,'ZLim')) + [-axzspan axzspan]);
    if frameno>1
        ax1clim(1) = min(myclim(1), ax1clim(1));
        ax1clim(2) = max(myclim(2), ax1clim(2));
        set(gca,'CLim', ax1clim);
        hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    end
    %hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    
    %myhb = savemesht_abs(:,myidx) + savemesht_abs_mr(:,myidx) + savemesht_rht(:,myidx) - ...
    %    savemesht_cnd(:,myidx) - savemesht_ems(:,myidx);
    
    subplot(2,3,4);
    %trisurf(TRI,mynx,myny,mynz, savemesht_cnd(:,myidx)'./myta*1e6, 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds); 
    trisurf(TRI,mynx,myny,mynz, savemeshn_hf(:,myidx)'./n_a*1e6, 'FaceColor', 'interp', 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds);
    set(gca,'Colormap', cmaphotcold);
    %mylims = get(gca','CLim');
    %set(gca,'CLim', [ -max(mylims) max(mylims) ]);
    axis equal
    axis off
    hc = colorbar('southoutside'); 
    hc.Label.String = 'W/m^2';
    set(gca,'fontsize',fs)
    view(va1, va2);
    title('Thermal conduction')
    set(hc,'fontsize',fs);
    if dolight mylt=lightangle(la1,la2); mylt.Color=[intens intens intens]; end
    hold on; plot3([mynx(constrainedEdges(:,1)) mynx(constrainedEdges(end,2))],[myny(constrainedEdges(:,1)) myny(constrainedEdges(end,2))],zreliefmag*([mynz(constrainedEdges(:,1)) mynz(constrainedEdges(end,2))]),'k'); hold off
    myclim=get(gca,'CLim');
    if frameno==1
        ax2clim = max(abs(myclim));
        set(gca,'CLim', [-ax2clim ax2clim]);
    end
    set(gca,'XLim', mean(get(gca,'XLim')) + [-axxspan axxspan]);
    set(gca,'YLim', mean(get(gca,'YLim')) + [-axyspan axyspan]);
    set(gca,'ZLim', mean(get(gca,'ZLim')) + [-axzspan axzspan]);
    if frameno>1
        ax2clim(1) = max( ax2clim, max(abs(myclim)));
        set(gca,'CLim', [-ax2clim ax2clim]);
        hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    end
    
    subplot(2,3,2);
    trisurf(TRI,mynx,myny,mynz, savemesht_abs(:,myidx)'./myta*1e6, 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds); 
    set(gca,'Colormap', cmaphot);
    axis equal
    axis off
    hc = colorbar('southoutside'); 
    hc.Label.String = 'W/m^2';
    set(gca,'fontsize',fs)
    view(va1, va2);
    title('Absorption (initial beam)');
    set(hc,'fontsize',fs);
    if dolight mylt=lightangle(la1,la2); mylt.Color=[intens intens intens]; end
    hold on; plot3([mynx(constrainedEdges(:,1)) mynx(constrainedEdges(end,2))],[myny(constrainedEdges(:,1)) myny(constrainedEdges(end,2))],zreliefmag*([mynz(constrainedEdges(:,1)) mynz(constrainedEdges(end,2))]),'k'); hold off
    myclim=get(gca,'CLim');
    if frameno==1
        ax3clim = myclim;
    end
    set(gca,'XLim', mean(get(gca,'XLim')) + [-axxspan axxspan]);
    set(gca,'YLim', mean(get(gca,'YLim')) + [-axyspan axyspan]);
    set(gca,'ZLim', mean(get(gca,'ZLim')) + [-axzspan axzspan]);
    if frameno>1
        ax3clim(1) = min(myclim(1), ax3clim(1));
        ax3clim(2) = max(myclim(2), ax3clim(2));
        set(gca,'CLim', ax3clim);
        hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    end
    
    subplot(2,3,3);
    trisurf(TRI,mynx,myny,mynz, savemesht_abs_mr(:,myidx)'./myta*1e6, 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds); 
    set(gca,'Colormap', cmaphot);
    axis equal
    axis off
    hc = colorbar('southoutside'); 
    hc.Label.String = 'W/m^2';
    set(gca,'fontsize',fs)
    view(va1, va2);
    title('Absorption (secondary reflections)')
    set(hc,'fontsize',fs);
    if dolight mylt=lightangle(la1,la2); mylt.Color=[intens intens intens]; end
    hold on; plot3([mynx(constrainedEdges(:,1)) mynx(constrainedEdges(end,2))],[myny(constrainedEdges(:,1)) myny(constrainedEdges(end,2))],zreliefmag*([mynz(constrainedEdges(:,1)) mynz(constrainedEdges(end,2))]),'k'); hold off
    myclim=get(gca,'CLim');
    if frameno==1
        ax4clim = myclim;
    end
    set(gca,'XLim', mean(get(gca,'XLim')) + [-axxspan axxspan]);
    set(gca,'YLim', mean(get(gca,'YLim')) + [-axyspan axyspan]);
    set(gca,'ZLim', mean(get(gca,'ZLim')) + [-axzspan axzspan]);
    if frameno>1
        ax4clim(1) = min(myclim(1), ax4clim(1));
        ax4clim(2) = max(myclim(2), ax4clim(2));
        set(gca,'CLim', ax4clim);
        hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    end
    
    
    subplot(2,3,6);
    trisurf(TRI,mynx,myny,mynz, savemesht_ems(:,myidx)'./myta*1e6, 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds);
    set(gca,'Colormap', cmapcold);
    axis equal
    axis off
    hc = colorbar('southoutside'); 
    hc.Label.String = 'W/m^2';
    set(gca,'fontsize',fs)
    view(va1, va2);
    title('Emission')
    set(hc,'fontsize',fs);
    hold on; plot3([mynx(constrainedEdges(:,1)) mynx(constrainedEdges(end,2))],[myny(constrainedEdges(:,1)) myny(constrainedEdges(end,2))],zreliefmag*([mynz(constrainedEdges(:,1)) mynz(constrainedEdges(end,2))]),'k'); hold off
    if dolight mylt=lightangle(la1,la2); mylt.Color=[intens intens intens]; end
    myclim=get(gca,'CLim');
    if frameno==1
        ax5clim = myclim;
    end
    set(gca,'XLim', mean(get(gca,'XLim')) + [-axxspan axxspan]);
    set(gca,'YLim', mean(get(gca,'YLim')) + [-axyspan axyspan]);
    set(gca,'ZLim', mean(get(gca,'ZLim')) + [-axzspan axzspan]);
    if frameno>1
        ax5clim(1) = min(myclim(1), ax5clim(1));
        ax5clim(2) = max(myclim(2), ax5clim(2));
        set(gca,'CLim', ax5clim);
        hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    end
    
    subplot(2,3,5);
    trisurf(TRI,mynx,myny,mynz, savemesht_rht(:,myidx)'./myta*1e6, 'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',as,'SpecularStrength',ss,'DiffuseStrength',ds); 
    %trisurf(TRI,mynx,myny,mynz, myhb'./myta*1e6, 'LineStyle','none');
    set(gca,'Colormap', cmaphot);
    axis equal
    axis off
    hc = colorbar('southoutside'); 
    hc.Label.String = 'W/m^2';
    set(gca,'fontsize',fs);
    view(va1, va2);
    title('Radiative heat transfer')
    set(hc,'fontsize',fs);
    if dolight mylt=lightangle(la1,la2); mylt.Color=[intens intens intens]; end
    hold on; plot3([mynx(constrainedEdges(:,1)) mynx(constrainedEdges(end,2))],[myny(constrainedEdges(:,1)) myny(constrainedEdges(end,2))],zreliefmag*([mynz(constrainedEdges(:,1)) mynz(constrainedEdges(end,2))]),'k'); hold off
    myclim=get(gca,'CLim');
    if frameno==1
        ax6clim = myclim;
    end
    set(gca,'XLim', mean(get(gca,'XLim')) + [-axxspan axxspan]);
    set(gca,'YLim', mean(get(gca,'YLim')) + [-axyspan axyspan]);
    set(gca,'ZLim', mean(get(gca,'ZLim')) + [-axzspan axzspan]);
    if frameno>1
        ax6clim(1) = min(myclim(1), ax6clim(1));
        ax6clim(2) = max(myclim(2), ax6clim(2));
        set(gca,'CLim', ax6clim);
        hold on; plot3( [0 0], [0 0], get(gca,'ZLim'), 'm'); hold off;
    end
    
    if length(myidxs)>1
        writeVideo(VW,getframe(gcf));
    end

end

if length(myidxs)>1
close(VW);
end
