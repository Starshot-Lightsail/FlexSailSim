%if (plotFast < 2)
    if (radius > 7)
        trisurf(TRI(~t_broken,:),n_x,n_y,n_z-comz+(ns-1)*xbox*2, t_t(~t_broken),'LineStyle','none','FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.5,'DiffuseStrength',0.5);
        plotMeshOutlineZ0=comz-(ns-1)*xbox*2;
        plotMeshOutline;
    else
        trisurf(TRI(~t_broken,:),n_x,n_y,n_z-comz+(ns-1)*xbox*2, t_t(~t_broken));%,'LineStyle','none');
        hold on;
        plot3(n_x(1), n_y(1), n_z(1)-comz+(ns-1)*xbox*2, 'k.');
    end
    set(gca,'CLim',[tscale_min tscale_max] );
    hold on;
%     text(-0,0,max(n_z)-comz+(ns-0.8)*xbox*2,['t=' num2str(tt)]);
    text(-0,0,max(n_z)-comz+(ns-0.8)*xbox*3,['t = ' num2str(tt), ' s']);
    
%end