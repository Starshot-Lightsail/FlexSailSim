% RRT bin histogram thing

figure('pos',[100        100        780         1024]);
set(gcf,'color','w'); 

myvideofilename = [ filebasename '_RRT'] ;
if exist( [ myvideofilename '.mp4'], 'file' )
   for n=1:100
       if ~exist( [myvideofilename '_' num2str(n) '.mp4'], 'file' )
           myvideofilename = [myvideofilename '_' num2str(n)];
           break;
       end
   end
end

V1 = VideoWriter(myvideofilename,'MPEG-4');  %alternate format for Ramon:  'Archival'
            V1.Quality = 100;
            open(V1);
numBins = RRT_ang_bins;
angRange = 30; %degrees

bins = linspace(-angRange, angRange, numBins+1);
binCtrs = ( bins(1:end-1) + bins(2:end) ) / 2;

downsampler = 4;  % set to 1 for all frames.  I typically use 4.  
bax = zeros(1,floor(nto/downsampler));
bay = bax;
bhx = bax;
bhy = bax;
blx = bax;
bly = bax;
bax0 = bax;
bay0 = bax;
bts = bax;
bI0s = bax;
bwfraction = 0.5; % fraction of reflected rays used to define beam width.  0.5 would give you 50% of rays, or the center two quartiles.  Must be >0 and <1.
bwquantiles = [.5-bwfraction/2 0.5 0.5+bwfraction/2];

t_fom_start0 = 0;
t_fom_end0 = 0.1;

t_fom_start1 = 0.3;
t_fom_end1 = 0.4;

n_fom0 = 0;
x_fom0 = 0;
y_fom0 = 0;

n_fom1 = 0;
x_fom1 = 0;
y_fom1 = 0;


no = 0;  %output index for beam plots

for n=1:downsampler:nto

no = no+1;
N = histcounts2(RRT_ax(:,n), RRT_ay(:,n),bins, bins);

bax0(no) = median(RRT_ax(:,n));
bay0(no) = median(RRT_ay(:,n));

xqs = quantile(RRT_ax(:,n), bwquantiles);
yqs = quantile(RRT_ay(:,n), bwquantiles);

bhx(no) = xqs(3);
bhy(no) = yqs(3);

bax(no) = xqs(2);
bay(no) = yqs(2);

blx(no) = xqs(1);
bly(no) = yqs(1);

bI0s(no) = I0s(n);

bts(no) = plot_time(n);

if ( (plot_time(n) > t_fom_start0) && (plot_time(n) < t_fom_end0) )
    n_fom0 = n_fom0+1;
    x_fom0 = x_fom0 + ( xqs(3) - xqs(1) );
    y_fom0 = y_fom0 + ( yqs(3) - yqs(1) );
elseif ( (plot_time(n) > t_fom_start1) && (plot_time(n) < t_fom_end1) )
    n_fom1 = n_fom1+1;
    x_fom1 = x_fom1 + ( xqs(3) - xqs(1) );
    y_fom1 = y_fom1 + ( yqs(3) - yqs(1) );
end

subplot(8,1,1:6);
imagesc([-angRange angRange], [-angRange angRange], N);
%axis square
%if (n == 1)
    xlabel('\theta_x (\circ)', 'FontSize', 20);
    ylabel('\theta_y (\circ)', 'FontSize', 20);
    set(gca,'FontSize',16);
%end
upperpos = get(gca,'Position');
subplot(8,1,8)
plot(plot_time, I0s);
%plot( [-.01 .01 .011 .1 .11 .21], [0 0 1 1 0 0]);
hold on;
plot(plot_time(n), I0s(n), 'ko');
set(gca,'YLim', [-0.2 1.2]);
xlabel('\it{t} (s)');
set(gca,'FontSize',16);
lowerpos = get(gca,'Position');
set(gca,'Position', [upperpos(1) lowerpos(2) upperpos(3) lowerpos(4) ]);
hold off;
writeVideo(V1,getframe(gcf));

if mod(100,n) == 0
    disp(num2str(n/nto));
    pause(0.01);
end

end

close(V1);
disp(['Wrote RTT video ' myvideofilename ]);

if ( n_fom0 > 0)
    x_fom0 = x_fom0 / n_fom0;
    y_fom0 = y_fom0 / n_fom0;
end

if ( n_fom1 > 0)
    x_fom1 = x_fom1 / n_fom1;
    y_fom1 = y_fom1 / n_fom1;
end


figure;
set(gcf,'color','w'); 


subplot(3,1,1);
plot(bts, bax, 'k', 'LineWidth', 2);
hold on;
plot(bts, bhx, 'k--' );
plot(bts, blx, 'k--' );
ylabel('\theta_x (\circ)');
title(['Reflected beam angle and ' num2str(bwfraction*100) '% beam width' ]);
set(gca,'FontSize',14);

subplot(3,1,2);
plot(bts, bay, 'LineWidth', 2);
hold on;
plot(bts, bhy, '--' );
plot(bts, bly, '--' );
ylabel('\theta_y (\circ)');
title(myvideofilename);
set(gca,'FontSize',14);

subplot(3,1,3);
plot(bts, bI0s);
xlabel('\it{t} (s)');
set(gca,'YLim',[-0.1 1.1]);
title(['Light gate / FOM0=' num2str(sqrt(x_fom0.^2+y_fom0.^2)) ' FOM1=' num2str(sqrt(x_fom1.^2+y_fom1.^2))]);
set(gca,'FontSize',14);

saveas(gcf,[myvideofilename '_RRT.fig']);
saveas(gcf,[myvideofilename '_RRT.png']);
        
