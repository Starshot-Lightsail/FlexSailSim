%figure;
        set(gcf,'color','w');
        hline = plot3(comxs(1:nto_firstBroken-1),comys(1:nto_firstBroken-1),plot_time(1:nto_firstBroken-1),'DisplayName',sprintf('%d %.3g',movieNum, rot0));
        if (nto_firstBroken < nto)
            hold on;
            plot3(comxs(nto_firstBroken:nto),comys(nto_firstBroken:nto),plot_time(nto_firstBroken:nto),'Color',get(hline,'Color'),'LineStyle',':');
            hline2=plot3(comxs(nto_firstBroken),comys(nto_firstBroken),plot_time(nto_firstBroken),'*');
            set(hline2,'Color',get(hline,'Color'));
        end
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('t (s)');
        box on;
        view(0,90);
        
        annotation('textbox', [.02 .95 .95 .04], 'String', ...
                                sprintf(['Trajectory for %s       \n' ...
                                         'pD=%g\t      radius=%g mm\t     rot0=%d Hz\t   \n' ...
                                         'I0=%d   \t     x0=%g mm\t     y0=%g mm\t     beamRadius=%g mm\n' ], ... 
                                         filebasename,  pD, radiusmm, rot0, I0, t0t_x0, t0t_y0, Irad), ...
                                         'EdgeColor', 'none', 'Interpreter', 'none' );
        
        set(gca,'Position',[0.1300    0.1100    0.7750    0.7150]);
        axis equal; 
        axis square;