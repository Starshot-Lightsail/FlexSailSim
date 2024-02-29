figure
%imagesc([plot_time(1) plot_time(end)], [1 length(meshcs)], monitorcs)
imagesc(monitorcs)
title('sail cross section, z height');
xlabel('Time (s)')
ylabel('Position across cross section (A.U.)')
hcolorbar = colorbar;
hcolorbar.Label.String = ('Z (mm)');

figure  %118479
plot([1:length(meshcs)], monitorcs(:,126254));
