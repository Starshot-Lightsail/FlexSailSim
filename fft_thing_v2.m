

fft_fund_rots = 1;
fft_stepback_ratio = 2;

fft_t_start = -spinup_time;
fft_t_stop = fft_fund_rots/rot0 + fft_t_start;
fft_t_start0 = fft_t_start;

fft_nto_start = round((fft_t_start + spinup_time) / (dnto ));
fft_nto_stop = round((fft_t_stop + spinup_time) / (dnto ));

if (fft_nto_start == 0)
    fft_nto_start = 1;
    fft_nto_stop = fft_nto_stop + 1;
end

fft_numpts = fft_nto_stop-fft_nto_start;
fft_nto_step = ceil(fft_numpts / fft_stepback_ratio );


Fs = 1/dnto;
f = Fs*(0:(fft_numpts/2))/fft_numpts;

nffts = floor((length(monitorrfs)-fft_nto_start)/fft_nto_step)-1-fft_stepback_ratio+1;
fftzarr = zeros(length(f),nffts);
fftrarr = zeros(length(f),nffts);
ffttarr = zeros(length(f),nffts);
ffttimes = (1:nffts) * (fft_t_stop - fft_t_start)/fft_stepback_ratio + fft_t_start; 
fftduration0 = (fft_t_stop - fft_t_start);

fft_nto_start0 = fft_nto_start;
fft_nto_stop0 = fft_nto_stop;

maxfftmag = 0;

for nfft = 1:nffts

myrfs = monitorrfs(:,fft_nto_start:fft_nto_stop);
mytfs = monitortfs(:,fft_nto_start:fft_nto_stop);
%myxfs = monitorfxs(:,fft_nto_start:fft_nto_stop);
%myyfs = monitorfys(:,fft_nto_start:fft_nto_stop);
myzfs = monitorfzs(:,fft_nto_start:fft_nto_stop);

%monitorfmxy = sqrt(myxfs.^2 + myyfs.^2);%+monitorzs.^2);
%monitorfmrt = sqrt(myrfs.^2 + mytfs.^2);

%xFFT = fft(sqrt(myrfs.^2+mytfs.*2)');
aFFT = fft(myzfs');

P2 = abs(aFFT)';
%P3 = angle(aFFT)';

P1 = P2(:,1:fft_numpts/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);

fftzarr(:,nfft)=sum(P1);
maxfftmag = max(maxfftmag,max(max(P1)));

%plot(f,sum(P1));

%set(gca,'YScale','log');
%set(gca,'XScale','log');
%hold on;

aFFT = fft(myrfs');

P2 = abs(aFFT)';
%P3 = angle(aFFT)';

P1 = P2(:,1:fft_numpts/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);

maxfftmag = max(maxfftmag,max(max(P1)));
fftrarr(:,nfft)=sum(P1);

aFFT = fft(mytfs');

P2 = abs(aFFT)';
%P3 = angle(aFFT)';

P1 = P2(:,1:fft_numpts/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);

maxfftmag = max(maxfftmag,max(max(P1)));
ffttarr(:,nfft)=sum(P1);

%figure;

% for n=1:50:length(n_x)
%     plot(plot_time,monitortfs(n,:));
%     hold on
% end
% 
% freqi = 90;
% 
% freqn = freqi / (Fs/nto);
% 
% cdata = P1(:,round(freqn));

%fft_nto_start = fft_nto_stop+1;
%fft_nto_stop = fft_nto_start + fft_numpts;
fft_nto_start = fft_nto_start+fft_nto_step;
fft_nto_stop = fft_nto_stop + fft_nto_step;


end

fftzarr_plot = log(fftzarr+1e-12);
fftrarr_plot = log(fftrarr+1e-12);
ffttarr_plot = log(ffttarr+1e-12);

redImage1 = fftzarr_plot';
greenImage1 = fftrarr_plot';
blueImage1 = ffttarr_plot';

minfftlogmag = min( min(  min(min(fftzarr_plot)), min(min(fftrarr_plot))), min(min(ffttarr_plot)));
fftzarr_plot = fftzarr_plot - minfftlogmag + 1;
fftrarr_plot = fftrarr_plot - minfftlogmag + 1;
ffttarr_plot = ffttarr_plot - minfftlogmag + 1;
maxfftlogmag = max( max(  max(max(fftzarr_plot)), max(max(fftrarr_plot))), max(max(ffttarr_plot)));


hWFF = figure('Position',[156         162        2566         1268]);
%hWFF = figure('Position',[2        1122        1918         428]);
hold on;

hafft1 = subplot(2,3,1);
surface(f,ffttimes,fftzarr_plot','LineStyle','none');
axis tight;
set(gca,'XScale','log');
box on;
set(gca, 'Layer', 'top');
xlabel('Frequency (Hz');
ylabel('Time (s)');
title('FFT Z');
set(gca,'CLim',[1 maxfftlogmag]);
%colorbar('eastoutside')
hold on;

hafft2 = subplot(2,3,2);
surface(f,ffttimes,fftrarr_plot','LineStyle','none');
axis tight;
set(gca,'XScale','log');
box on;
set(gca, 'Layer', 'top');
xlabel('Frequency (Hz');
ylabel('Time (s)');
set(gca,'CLim',[1 maxfftlogmag]);
%colorbar('eastoutside')
title('FFT R');
hold on;

hafft3 = subplot(2,3,3);
surface(f,ffttimes,ffttarr_plot','LineStyle','none');
axis tight;
set(gca,'XScale','log');
box on;
set(gca, 'Layer', 'top');
xlabel('Frequency (Hz');
ylabel('Time (s)');
title('FFT T');
set(gca,'CLim',[1 maxfftlogmag]);
%colorbar('eastoutside')
hold on;



the_real_maxest_max_of_all_the_max = max( max(  max(max(redImage1)), max(max(greenImage1)) ), max(max(blueImage1)) );

rgbImage1 = cat(3, redImage1, greenImage1, blueImage1 )./the_real_maxest_max_of_all_the_max;


hafft4 = subplot(2,3,[4 5 6]);
hold off;
imshow( rgbImage1 );
hafft4 = gca;
axis xy
axis normal
axis tight

saveas(gcf,['sim' num2str(movieNum) '_fft_' num2str(fft_fund_rots) '_' num2str(fft_stepback_ratio) '.fig']);
saveas(gcf,['sim' num2str(movieNum) '_fft_' num2str(fft_fund_rots) '_' num2str(fft_stepback_ratio) '.png']);

%set(gca,'Position',get(gca,'OuterPosition'));
pause;

%hfftmpf = figure('Position',[155         822        2566        1194]);
hfftmpf = figure('Position',[27        1452        1843         624]);

%% OK now it's time to pick the picker!!!

for numFFTclicks = 1:1000
    
    figure(hWFF);
    
    [csrx, csry, button] = ginput(1);
    
    if (isempty(button))
        break;
    elseif button > 3
        disp('ok you can zoom and stuff.')
        pause();
        continue
    end
    
    if isempty(csrx)
        break
    end
    
   % try
    % check which axis it's from.    For now I don't care/
    
    disp([' CSR ' num2str(csrx) '  ' num2str(csry) ] );
    
    if ( gca == hafft4 ) % image panel
        disp('found in image panel');
        freqi = f(max(floor(csrx),1));
        ncycle = csry;
        timei = ffttimes(round(csry));
    else
        freqi = csrx;
        timei = csry;
    end
    
   
    
    fft_t_start = timei;
    fft_t_stop = timei + fft_fund_rots/rot0;
    
    fft_nto_start = round((fft_t_start + spinup_time) / (dnto ));
    fft_nto_stop = round((fft_t_stop + spinup_time) / (dnto ));
    
    fft_numpts = fft_nto_stop-fft_nto_start;
    
    Fs = 1/dnto;
    f = Fs*(0:(fft_numpts/2))/fft_numpts;
    
    %nffts = floor((length(monitorrfs)-fft_nto_start)/fft_numpts)-1;
    fftzsgl = zeros(length(f),1);
    fftrsgl = zeros(length(f),1);
    ffttsgl = zeros(length(f),1);
    ffttimes_zoom = (1:nffts) * (fft_t_stop - fft_t_start) + fft_t_start;
    fftduration0 = (fft_t_stop - fft_t_start);
    
   % fft_nto_start0 = fft_nto_start;
   % fft_nto_stop0 = fft_nto_stop;
    
    myrfs = monitorrfs(:,fft_nto_start:fft_nto_stop);
    mytfs = monitortfs(:,fft_nto_start:fft_nto_stop);
    %myxfs = monitorfxs(:,fft_nto_start:fft_nto_stop);
    %myyfs = monitorfys(:,fft_nto_start:fft_nto_stop);
    myzfs = monitorfzs(:,fft_nto_start:fft_nto_stop);

    zFFT = fft(myzfs');
    rFFT = fft(myrfs');
    tFFT = fft(mytfs');
    
    P2z = abs(zFFT)';
    P3z = angle(zFFT)';
    P2r = abs(rFFT)';
    P3r = angle(rFFT)';
    P2t = abs(tFFT)';
    P3t = angle(tFFT)';

    P1z = P2z(:,1:fft_numpts/2+1);
    P1z(:,2:end-1) = 2*P1z(:,2:end-1);
    P1r = P2r(:,1:fft_numpts/2+1);
    P1r(:,2:end-1) = 2*P1r(:,2:end-1);
    P1t = P2t(:,1:fft_numpts/2+1);
    P1t(:,2:end-1) = 2*P1t(:,2:end-1);
    
    freqn = max(freqi / (Fs/fft_numpts)+1, 1) ;

    disp([' FREQ ' num2str( f(floor(freqn)) )  '    TIME ' num2str(graphts(fft_nto_start)) ' s' ] );
    
    myf_fzm = P2z(:,floor(freqn));
    myf_fzp = P3z(:,floor(freqn));
    myf_frm = P2r(:,floor(freqn));
    myf_frp = P3r(:,floor(freqn));
    myf_ftm = P2t(:,floor(freqn));
    myf_ftp = P3t(:,floor(freqn));
    
    
    if (button == 3)
        V2 = VideoWriter([ 'sim' num2str(movieNum) '_fft_' num2str(fft_fund_rots) '_' num2str(fft_stepback_ratio) '_F' num2str(round(freqn)) '_T' num2str(round(fft_nto_start)) ],'MPEG-4');
        V2.Quality = 100;
        open(V2);

        stupidcolormapgen;
        mymaxfftmag = max( [max(myf_fzm) max(myf_frm) max(myf_ftm) ]);
        hdeleteme = figure('Position',[84        1571         560         420]);
        colormap(mystupidmap);
        myf_fzm = myf_fzm./mymaxfftmag;
        myf_frm = myf_frm./mymaxfftmag;
        myf_ftm = myf_ftm./mymaxfftmag;
        
        mythetas = 0:pi/32:2*pi*fft_fund_rots;
        
        num_nto_per_cycle = 1/f(floor(freqn))/dnto;
        
        
        
        for n = 1:length(mythetas)
            mytheta = mythetas(n);
            if mytheta > 0
                oldviewangle = get(gca,'View');
            end
            mynto = fft_nto_start + round((n/length(mythetas))*num_nto_per_cycle*fft_fund_rots);
            
            %cdata3d = [ cmapsteps*round(myf_fzm.*(0.5+0.5*cos(mytheta - myf_fzp)));
            cdata3d = [ floor((0.5+0.5.*myf_fzm.*sin(myf_fzp+mytheta)).*cmapsteps) ...%+ ...  cmb*cmapsteps*
                floor((0.5+0.5.*myf_frm.*sin(myf_frp+mytheta)).*cmapsteps)  ...%+ ...  cmb*
                floor((0.5+0.5.*myf_ftm.*sin(myf_ftp+mytheta)).*cmapsteps) ];%+ 1;
            cdata = (cdata3d(:,1)*cmapsteps^2 + cdata3d(:,2)*cmapsteps + cdata3d(:,3)   )/(cmapsteps^3);
            
            
            
            %disp(['Step ' num2str(mynto) '   ' num2str(plot_time(mynto)) ]);
            
            %cdata = sqrt( myf_fzm.^2 + myf_frm.^2 + myf_ftm.^2);
            if (gcf ~= hdeleteme)
                figure(hdeleteme);
            end
            hold off;
            trimesh(TRI,monitorxs(:,mynto),monitorys(:,mynto),monitorzs(:,mynto),  cdata , ...
                'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5, ...
                'FaceColor','flat','EdgeColor','None' );
            set(gca,'CLim',[0 1-1/length(mystupidmap)]);
            myzlim=get(gca,'ZLim');
            hold on;
            plot3(monitorxs(end,mynto), monitorys(end,mynto), monitorzs(end,mynto), 'k*');
            plot3([Ictrx Ictrx], [0 0], myzlim, 'm-o');
            title(['\it{f} = \rm' num2str(f(floor(freqn))) ' Hz,   t = ' num2str(plot_time(mynto)) ' s,  \theta = ' num2str(mytheta./pi) '\pi rad']);
            xlabel('x (mm)');
            ylabel('y (mm)');
            %lightangle(20,0);
            %lightangle(-100,0);
            if mytheta>0
                view(oldviewangle)
            else
                view(0,90);
            end
            writeVideo(V2,getframe(gcf));
            pause(0.1)
        end
        %pause();
        close(hdeleteme);
        close(V2);
        
    else
        
figure(hfftmpf)

    cdata1 = log(myf_fzm)+17;
    cdata2 = myf_fzp;
    cdata3 = log(myf_frm)+17;
    cdata4 = myf_frp;
    cdata5 = log(myf_ftm)+17;
    cdata6 = myf_ftp;
        
   
asdf1=subplot(2,3,1);
oldview = get(gca,'View');
hold off;
trimesh(TRI,monitorxs(:,fft_nto_stop),monitorys(:,fft_nto_stop),monitorzs(:,fft_nto_stop),  cdata1, ...
        'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5, ...
        'FaceColor','interp','EdgeColor','None' );
hold on
colorbar('southoutside');
%set(gca,'CLim',[0 maxfftclim]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal
title({['ModeView, \it{F_Z}, amplitude, f=' num2str(freqi) ' Hz'] });%, ['t=' num2str(tt) ' s' ] });
set(gca,'View',oldview);


asdf2=subplot(2,3,4);
hold off
    hAwesome = trimesh(TRI,monitorxs(:,fft_nto_stop),monitorys(:,fft_nto_stop),monitorzs(:,fft_nto_stop),  cdata2,'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
    hold on
    set(hAwesome,'FaceColor','interp');
    set(hAwesome,'EdgeColor','None');
    %lightangle(-120,60);
    %plotMeshOutline;
hold on;
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal
% colorbar('southoutside')
title({['ModeView, \it{F_Z}, phase, \it{f}=' num2str(freqi) 'Hz'] });%, ['t=' num2str(tt) ' s' ] });
set(gca,'View',oldview);



asdf3=subplot(2,3,2);
hold off;
trimesh(TRI,monitorxs(:,fft_nto_stop),monitorys(:,fft_nto_stop),monitorzs(:,fft_nto_stop),  cdata3, ...
        'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5, ...
        'FaceColor','interp','EdgeColor','None' );
hold on
colorbar('southoutside');
%set(gca,'CLim',[0 maxfftclim]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal
title({['F_R , amplitude, f=' num2str(freqi) ' Hz'] });%, ['t=' num2str(tt) ' s' ] });
set(gca,'View',oldview);


asdf4=subplot(2,3,5);
hold off;
trimesh(TRI,monitorxs(:,fft_nto_stop),monitorys(:,fft_nto_stop),monitorzs(:,fft_nto_stop),  cdata4, ...
        'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5, ...
        'FaceColor','interp','EdgeColor','None' );
hold on
%set(gca,'CLim',[0 maxfftclim]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal
title({['F_R, phase, f=' num2str(freqi) ' Hz'] });%, ['t=' num2str(tt) ' s' ] });
set(gca,'View',oldview);


asdf5=subplot(2,3,3);
hold off;
trimesh(TRI,monitorxs(:,fft_nto_stop),monitorys(:,fft_nto_stop),monitorzs(:,fft_nto_stop),  cdata5, ...
        'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5, ...
        'FaceColor','interp','EdgeColor','None' );
hold on
colorbar('southoutside');
%set(gca,'CLim',[0 maxfftclim]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal
title({['F_T, amplitude, f=' num2str(freqi) ' Hz'] });%, ['t=' num2str(tt) ' s' ] });
set(gca,'View',oldview);



asdf6=subplot(2,3,6);
hold off;
trimesh(TRI,monitorxs(:,fft_nto_stop),monitorys(:,fft_nto_stop),monitorzs(:,fft_nto_stop),  cdata6, ...
        'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5, ...
        'FaceColor','interp','EdgeColor','None' );
hold on
%set(gca,'CLim',[0 maxfftclim]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('z (mm)');
axis equal
title({['F_T, phase, f=' num2str(freqi) ' Hz'] });%, ['t=' num2str(tt) ' s' ] });
set(gca,'View',oldview);
linkprop([asdf1 asdf2 asdf3 asdf4 asdf5 asdf6],'View');

saveas(gcf,['sim' num2str(movieNum) '_fft_' num2str(fft_fund_rots) '_' num2str(fft_stepback_ratio) '_F' num2str(round(freqn)) '_T' num2str(round(fft_nto_start)) '.fig']);
saveas(gcf,['sim' num2str(movieNum) '_fft_' num2str(fft_fund_rots) '_' num2str(fft_stepback_ratio) '_F' num2str(round(freqn)) '_T' num2str(round(fft_nto_start)) '.png']);


    end

end