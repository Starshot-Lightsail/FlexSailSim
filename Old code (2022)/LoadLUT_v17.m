%% LoadLUT -- loads lookup-table data of optical properties for lightsail simulations
% (c) May 2021 Michael Kelzenberg, California Institute of Technology
%
% Look up tables (LUTs) are used to determine the optical properties of lightsail surfaces in a computationally
% efficient way during simulations.  While the optical response of simple surfaces, e.g., homogonous isotripic
% dielecric stacks can be calculated quickly using closed-form equations (e.g., TMM), more generally, it is 
% inefficient or downright intractable to calculate the optical response of complex nanophotonic surfaces during 
% lightsail mesh simulations.  Instead, we assume that such surfaces have been characterized by the appropriate methods
% in advance (perhaps *gasp* even experimentally), and simply look up the correct pre-calculated value from a table
% during the mesh simulations.  Many such tables can be loaded into memory, allowing us to include many different
% materials or nanophotonic textures in our lightsail designs.
%
% This function loads a single set of data from a tabulated data source and constructs a LUT structure suitable for use
% in lightsail mesh simulations.  Specifically, it loads data from Ramon's simulations, which come as multi-page tables 
% with the following column format:
%
% Psi(rad)  Theta(rad)  PressForce_trans   PressForce_normal PressForce_parallel absorption
%
% (The last two columns I added for future expandability and are not present in Ramon's data as of 4/11/2021.)
%   - If only four data columns are present, the function adds a fifth column of zeros for PressForce_parallel.
%   - Then, if only five data columns are present, the funciton adds a sixth column, all with the value of 
%     'Iabs' from the global workspace, for the 'absorption' LUT.
%
% I'm not exactly sure what the units are for the pressure columns, but the values are such that multiplying each 
%  triangle's incident power (in watts), by that triangle's pressure values from the LUT, via the appropriate vector
%  math, produces the correct optical thrust force vector (in ... I'm pretty sure newtons--or whatever works correctly 
%  in the simulation code.)
%
% The units for the 'absorption' table are that of absolute absorption (0 to 1), such that multiplying the power
% (in watts) incident upon each triangle, by that triangle's absorption LUT value, produces the absorbed power (in
% watts) for that triangle.  NOTE:  Absorption values are used for thermal calculations only.  The simulator does not
% consider the momentum of absorbed light in its calculation of optical forces, as this force should already be present
% in the LUT values for pressure forces!
%
% Data tables can be split over arbitrary pages for convenience; they will
% be concatinated in order to produce a single listing of all simulation
% results with the above column structure.  Rows with all zeros will be
% skipped, i.e., each page can be zero-padded.  The spacing of Psi and Theta
% values can be varied to provide greater resolution in regions of greater
% interest.  However, the following ordering is required for the concatinated table:
%   - Column 1 (Psi) is the outer sweep variable.  It must remain constant throughout each inner variable (Theta) 
%     sweep, and must increase monotonically (by a nonzero amount) between each inner variable sweep.  
%
%     Note:  Psi is the *parallel tilt* direction.  The parallel plane is spanned by the surface normal (up) and the axis
%     of (along) the optical grating.  It is the "out of plane" direction for 2D grating simulations.  
%
%     Note 2:  Because some of Ramon's simulations occur in a rotated coordinate system with respect to the grating
%     axis, the function input variable "temp_dorotate" can be used to swap between Psi and Theta for some input files.
%     (See code below for implementation.)  However, I recommend that all future simulations be performed in a 
%     coordinate space consistent with the definitions of Psi and Theta described here, such that "temp_dorotate" can 
%     be zero.  In the future, we'll figure out how to deal with polarization in a consistent way, which will make this
%     convention even easier to follow.  
%
%   - Column 2 (Theta) is the inner sweep variable.  It must increase monotonically, without repeated values, within 
%     each block of Psi values.  The span and spacing of the Theta sweep need not be identical within each block of Psi 
%     values, and the Theta step can vary arbitrarily between each row, so long as it is positive and nonzero.    
%
%     Note:  Theta is the *transverse tilt* direction.  The transverse plane is spanned by the the surface normal (up) 
%     and the axis perpendicular to the grating along the surface.  In other words, the transverse plane is the one in 
%     which we draw the cross-section of a grating (for linear gratings).  Or in other other words, this is the angle
%     that is varied in a grating monochromator.  
%
%     (See also Note 2 above)
%
%% INPUTS:  
%
% filename:  The name/path of the .mat file contianing the data.  It must be properly specified for the current working
% path.  Its contents will be loaded in entirety into the base workspace.  
%
% varname:  The name of the data table as it will appear in the workspace after loading the above .mat file.
%
% LUTPsiStepDeg, LUTThetaStepDeg:  The step size for Psi and Theta to be used in constructing the LUT.  
%   Note: 2D linear interpolation will be used to construc the LUTs from the source data.
%   Note 2:  As of this writing, the main simulation script requries that all LUTs use the same value for Psi and Theta
%   steps.  This makes the LUT code easier and faster.  Thus, all LUTs must use the same values for angle steps.
%
% showplots:  Whether or not to display plots of the LUTs.  
%  Set to a nonzero value if you appreciate true art!
%  1: show colorful images of the LUTn and LUTt data (1 figure)
%  2: show the above plus a source data density plot (2 figures)
%  3: Overlay the source data density plot on the LUT images instead (1 figure)
%  0: Do not plot data because you're a mindless robot who doesn't appreciate the beauty of the source data.
%
% temp_dorotate:  a flag (of value 0 or 1) used to rotate Ramon's "TM" simulation tables to match the coordinate
% convention of the v12 simulator.  Use a value of 1 for "TM" data files, and 0 for "TE" data files.  See code below for
% implementation of the transformation.
% 
% 
%% OUTPUTS:
%
% This function produces a LUT structure with the following components:
% 
% LUT.n, LUT.t, LUT.p, LUT.a:  2D data tables for normal, transverse, and
% parallel forces, and for absorption.  The data tables span Psi and Theta
% ranges from -90 to +90 degrees.  If the input data tables do no span
% all the way to -90 or +90 degrees, values of zero will be assigned for
% angles exceeding available input data.  You should try to supply reasonable
% values for the entire incidence hemisphere. 
%
% LUT.pstep, LUT.tstep: The angle step size (radians) for Psi and Theta in the look up tables.
%  Note:  as mentioned above, the simulation code does not reference these values, because all LUTs are assumed to have
%  the same pstep and tstep values.  These outputs are for convience/forensic purposes only!
%
% LUT.file, LUT.name:  The file path (string) and workspace variable name from which the data was loaded (echoed from
%   the input parameters, again for convenience/forensic purposes only.)




%% Limitations of LUT implementation: 
%
%% Polariation:  
% As of v17, we still haven't developed support for polarization in the simulator.  So for now, you must design your 
% lightsail and LUTs in a way that is explicitly consistent with your assumed polarization scenario.  This is a pretty 
% substantial shortcoming, but so far we've been able to sleep at night by considering only horizontally and vertically
% oriented MEPS textures in our designs, and using TM and TE simulation data for each, respectively.  For spin-stabilized
% sails, this implies that the beam polarization must also rotate in perfect step with the sail.  
%
% I plan to add polarization support in the future, in which case the number of tables required for each LUT structure
% will be doubled, one for each of two orthogonal polarization states.  That will let us properly handle a broad set of 
% polarization scenarios for which incident light can be described by a combination of two orthogonal polarization 
% states, such as linear polarization (TE/TM) or circular polarization (CW/CCW).  
%
%% Temperature dependence
% As of v17, we haven't yet implemented support for temperature-dependent optical properties.  I plan to add this in the
% future, in which case multiple LUT data sets will be required, one for each temperature interval.
%
%% Reverse illumination
% Our simulator is designed to simulate flat or topoligically flattenable lightsail shapes only.  An implicit assumption
% is that all sail surfaces will be illuminated from the 'bottom' side only, at all times during a successful launch.  
% This precludes simulating such lightsail shapes as complete spheres or other closed surfaces.  Changing this
% assumption would require a substantial overhaul of the entire simulator.  We have no intention of ever doing this, 
% thus, LUTs only need to span the incidence hemisphere.  Simply expanding the LUTs to include reverse illumation angles
% won't make the simulator work correctly for spheres.    

function LUT = LoadLUT_v17(filepath, varname, LUTPsiStepDeg, LUTThetaStepDeg, showplots, temp_dorotate)

fprintf('\n');
disp('LoadLUT()  v12   MK 12 Apr 2021');
disp(['File source: ' filepath] );
disp(['Var name: ' varname] );

evalin('base',['load(''' filepath ''');']);

srctable = evalin('base', varname );

bigtable = (squeeze(srctable(:,:,1)));

for np = 2:(size(srctable,3))
    
    bigtable = [bigtable; (squeeze(srctable(:,:,np))) ];
    %table2 = (squeeze(pressure_TE_small_tables(:,:,2)));
%table3 = (squeeze(pressure_TE_small_tables(:,:,3)));
%table4 = (squeeze(pressure_TE_small_tables(:,:,4)));
%table5 = (squeeze(pressure_TE_small_tables(:,:,4)));
%table6 = (squeeze(pressure_TE_small_tables(:,:,4)));
%table7 = (squeeze(pressure_TE_small_tables(:,:,4)));
end
%table1 = table1(any(table1,2),:); %
%table2 = table2(any(table2,2),:); 
%table3 = table3(any(table3,2),:); 
%table4 = table4(any(table4,2),:); 
bigtable = bigtable(any(bigtable,2),:);

if temp_dorotate
    temp_var = LUTPsiStepDeg;
    LUTPsiStepDeg = LUTThetaStepDeg;
    LUTThetaStepDeg = temp_var;
end

S = whos('srctable');
disp(['Source table size: ' num2str(size(srctable)) ' (' num2str(S.bytes/1024) ' kB)' ]);

if (size(bigtable,2)) == 4
    bigtable(:,5) = 0;
end
if (size(bigtable,2)) == 5
    bigtable(:,6) = evalin('base','Iabs');  % temporary fix for our simulation tables as of v15 which lack absorption data
end
if size(bigtable,2) ~= 6
    error('Error, data table doesn''t have the correct number of columns');
end

S = whos('bigtable');
disp(['Reformed table: ' num2str(size(bigtable)) ' (' num2str(S.bytes/1024) ' kB)' ]);

if showplots==2
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/2-10 scrsz(4)/2-100 scrsz(3)/2 scrsz(4)/2]);
    plot(rad2deg(bigtable(:,1)), rad2deg(bigtable(:,2)), '.');
    xlabel('\psi (\circ)');
    ylabel('\theta (\circ)');
    set(gco,'MarkerSize',2)
    title([varname '  (' filepath ')' ],'Interpreter' , 'none');
end

disp(['Data table info:']);
disp([' Size:  ' num2str(size(bigtable))]);
minPsi = min(bigtable(:,1));
maxPsi = max(bigtable(:,1));
disp([' Psi range:  ' num2str(minPsi) ' to ' num2str(maxPsi) '   ( ' num2str(rad2deg(minPsi)) ' to ' num2str(rad2deg(maxPsi)) ' degrees )']);
B = abs(diff(bigtable(:,1)));
B(B==0) = NaN;
disp(['     steps:                  ' num2str(rad2deg(min(B))) ' to ' num2str(rad2deg(max((diff(bigtable(:,1)))))) ]);
minTheta = min(bigtable(:,2));
maxTheta = max(bigtable(:,2));
disp([' Theta range:  ' num2str(minTheta) ' to ' num2str(maxTheta) '   ( ' num2str(rad2deg(minTheta)) ' to ' num2str(rad2deg(maxTheta)) ' degrees )']);
disp(['     steps:                  ' num2str(rad2deg(min(abs(diff(bigtable(:,2)))))) ' to ' num2str(rad2deg(max((diff(bigtable(:,2)))))) ]);



%construct big table
%LUTPsiStepDeg = 0.125; %degrees
%LUTThetaStepDeg = 0.125;

%I've decided that theta and phi must always span +- 90 degrees, so the
%following variables are no longer used.  The values are hard coded.  
%Deal with it.
%LUTPsiStartDeg = -90;
%LUTPsiStopDeg = 90;
%LUTThetaStartDeg = -90;
%LUTThetaStopDeg = 90;

LUTPsis = deg2rad(-90:LUTPsiStepDeg:90);
LUTThetas = deg2rad(-90:LUTThetaStepDeg:90);

psiSteps = unique(bigtable(:,1));
psiStepsI = [];
nn = 1;
for n=1:length(bigtable(:,1))
    if (psiSteps(nn) == bigtable(n,1)) %found first instance of next value
        psiStepsI(nn) = n;
        if nn==length(psiSteps);
            break;
        else
        nn = nn+1;
        end
    else if psiSteps(nn) < bigtable(n,1)
        error('Error in ingest processing of MEPS data table: non-monotonic ordering of psi values');
        end
    end
end
if ( (nn) ~= length(psiSteps) )
    error('Error in ingest processing of MEPS data table, maybe non-monotonic ordering of psi values?');
end
psiStepsI(end+1) = length(bigtable(:,1))+1;
psiStepsI = psiStepsI';  % turn into row vector for no reason at all

% checking for sweep issues in source data...
% figure
% hold on;
% for n=1:5:length(psiSteps)
%     plot(bigtable(psiStepsI(n):psiStepsI(n+1),2),bigtable(psiStepsI(n):psiStepsI(n+1),3),'.-');
% end
    
    

LUTn = zeros(length(LUTPsis),length(LUTThetas));
LUTt = zeros(length(LUTPsis),length(LUTThetas));
LUTp = zeros(length(LUTPsis),length(LUTThetas));
LUTa = zeros(length(LUTPsis),length(LUTThetas));

psiIndexLower = 0;
psiIndexUpper = 0;
psiIndexF = 0;
myPsi = -1234;

tic;
LUT_t_start = tic;
for nx = 1:length(LUTPsis)
    if (toc) > 5
        disp(['Building LUT row ' num2str(nx) '/' num2str(length(LUTPsis)) ]);
        tic;
    end
    myPsi = LUTPsis(nx);
    if (myPsi < min(psiSteps)) || (myPsi > max(psiSteps))  % %  these 6 lines
        LUTn(nx,:) = zeros(1,length(LUTThetas));           % %  cause us to use zeros
        LUTt(nx,:) = zeros(1,length(LUTThetas));           % %  for all rows where psi
        LUTp(nx,:) = zeros(1,length(LUTThetas));           % %  exceeds the limits of 
        LUTa(nx,:) = zeros(1,length(LUTThetas));           % %  avaialable data
        continue;                                          % %  ...
%     if (myPsi < min(psiSteps))                   % these four lines
%         psiIndexF = 1;                           % cause us to stretch perimiter psi 
%     elseif (myPsi > max(psiSteps))               % values to the +-90d limits, instead 
%         psiIndexF = length(psiSteps);            % of using zeros (above)
    else
        psiIndexF = interp1(psiSteps,1:length(psiSteps),myPsi);
    end
    psiIndexLower = floor(psiIndexF);
    psiIndexUpper = ceil(psiIndexF);
   % if (psiIndexF ~= psiIndexLower)
        psiIndexF = psiIndexF - psiIndexLower;
   % end
   
   myThetasLower = bigtable(psiStepsI(psiIndexLower):(psiStepsI(psiIndexLower+1)-1),2);
   myThetasUpper = bigtable(psiStepsI(psiIndexUpper):(psiStepsI(psiIndexUpper+1)-1),2);
   
   mynLower = bigtable(psiStepsI(psiIndexLower):(psiStepsI(psiIndexLower+1)-1),4);
   mynUpper = bigtable(psiStepsI(psiIndexUpper):(psiStepsI(psiIndexUpper+1)-1),4);
  
   mytLower = bigtable(psiStepsI(psiIndexLower):(psiStepsI(psiIndexLower+1)-1),3);
   mytUpper = bigtable(psiStepsI(psiIndexUpper):(psiStepsI(psiIndexUpper+1)-1),3);
   
   mypLower = bigtable(psiStepsI(psiIndexLower):(psiStepsI(psiIndexLower+1)-1),5);
   mypUpper = bigtable(psiStepsI(psiIndexUpper):(psiStepsI(psiIndexUpper+1)-1),5);
      
   myaLower = bigtable(psiStepsI(psiIndexLower):(psiStepsI(psiIndexLower+1)-1),6);
   myaUpper = bigtable(psiStepsI(psiIndexUpper):(psiStepsI(psiIndexUpper+1)-1),6);
  
   LUTn(nx,:) = (1-psiIndexF)    *  interp1( myThetasLower, mynLower, LUTThetas )  + ...
                (0+psiIndexF)    *  interp1( myThetasUpper, mynUpper, LUTThetas );
   LUTt(nx,:) = (1-psiIndexF)    *  interp1( myThetasLower, mytLower, LUTThetas )  + ...
                (0+psiIndexF)    *  interp1( myThetasUpper, mytUpper, LUTThetas );
   LUTp(nx,:) = (1-psiIndexF)    *  interp1( myThetasLower, mypLower, LUTThetas )  + ...
                (0+psiIndexF)    *  interp1( myThetasUpper, mypUpper, LUTThetas );
   LUTa(nx,:) = (1-psiIndexF)    *  interp1( myThetasLower, myaLower, LUTThetas )  + ...
                (0+psiIndexF)    *  interp1( myThetasUpper, myaUpper, LUTThetas );            
   
  %wow, switching in the above four lines of vectorized code, in place of the below
  %inner for loop, reduced run time from ~25 seconds to .16 seconds!
            
%    for ny=1:length(LUTThetas)
%       if (LUTThetas(ny) < minTheta )
%           % fuck
%       elseif (LUTThetas(ny) > maxTheta) 
%           % dammit
%       else
%                
%         LUTn(nx,ny) = (1-psiIndexF)    *  interp1( myThetasLower, mynLower, LUTThetas(ny) )  + ...
%                       (0+psiIndexF) *  interp1( myThetasUpper, mynUpper, LUTThetas(ny) );
%                   
%         LUTt(nx,ny) = (1-psiIndexF)    *  interp1( myThetasLower, mytLower, LUTThetas(ny) )  + ...
%                       (0+psiIndexF) *  interp1( myThetasUpper, mytUpper, LUTThetas(ny) );
%       end               
%    end
end

LUTn(isnan(LUTn)) = 0;
LUTt(isnan(LUTt)) = 0;
LUTp(isnan(LUTp)) = 0;
LUTa(isnan(LUTa)) = 0;

 
if temp_dorotate
    LUTn = fliplr(LUTn');
    LUTt = fliplr(LUTt');
    LUTp = fliplr(LUTp');
    LUTa = fliplr(LUTa');
    temp_var = LUTPsiStepDeg;
    LUTPsiStepDeg = LUTThetaStepDeg;
    LUTThetaStepDeg = temp_var;
end


disp(['Elapsed time for LUT generation: ' num2str(toc(LUT_t_start)) ' sec.']);

if (~exist('forceColormap'))
        load forceColormap_v17.mat;
end

% test the stupid image plot...
% uncomment these lines to mark the (90p, 8?t) pixel on the data imagesc
% plots to confirm the axes are labelled correctly.
%LUTn(end,end-10) = -max(max(abs(LUTn)));
%LUTt(end,end-10) = max(max(abs(LUTt)));



LUT = [];
LUT.n = (-fliplr(LUTn));  %negative to correct for simulation orientation 
LUT.t = (-fliplr(LUTt));
LUT.p = (-fliplr(LUTp));
LUT.a = fliplr(LUTa);
LUT.pstep = deg2rad(LUTPsiStepDeg);
LUT.tstep = deg2rad(LUTThetaStepDeg);
LUT.file = filepath;
LUT.varname = varname;


disp(['LUT dimensions:  ' num2str(size(LUTn,1)) ' x ' num2str(size(LUTn,2)) ]);
disp(['LUT resolution: ' num2str(LUTPsiStepDeg) ' x ' num2str(LUTThetaStepDeg) ]);
S = whos('LUTn');
disp(['LUT page size: ' num2str(S.bytes/1024) 'kB' ]);
S = whos('LUT');
disp(['Total LUT memory footprint:  ' num2str(S.bytes/1024) ' kB']);
fprintf('\n');

    if showplots
    scrsz = get(0,'ScreenSize');
    
    figure('Position',[scrsz(3)/2-10 scrsz(4)/2-100 scrsz(3)/2 scrsz(4)/2]);
    colormap parula;

    subplot(1,2,1);
    imagesc(([-90 90]), ([ -90 90]), LUT.t);
    title('Transverse force');
    xlabel('\Theta (\circ)');
    ylabel('\Psi (\circ)');
    climmax = max(abs(get(gca,'CLim')));
    set(gca,'CLim',[-climmax climmax]);
    set(gca,'Colormap',forceColormap2)
    colorbar;
    
    if showplots==3
        hold on;
        plot(rad2deg(bigtable(:,2)), rad2deg(bigtable(:,1)), 'k.', 'MarkerSize', 1)
    end

    subplot(1,2,2);
    imagesc(([-90 90]), ([-90 90]), LUT.n);
    title('(-) Normal force');
    xlabel('\Theta (\circ)');
    ylabel('\Psi (\circ)');
    colorbar;
    
    if showplots==3
        hold on;
        plot(rad2deg(bigtable(:,2)), rad2deg(bigtable(:,1)), 'k.', 'MarkerSize', 1)
    end
    
    annotation('textbox', [0, 1, 0, 0], 'string', [varname '(' filepath ')' ],'Interpreter' , 'none');
    
end

