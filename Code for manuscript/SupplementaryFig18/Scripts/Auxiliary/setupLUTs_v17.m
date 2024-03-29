%% setupLUTs
% STARSHOT LIGHTSAIL SIMULATOR
% (c) 2018-2021 Michael Kelzenberg, Ramon Gao -- California Institute of Technology
% subroutine script:  setupLUTs
%
% Loads LUTs used in mesh simulations.  This is separated as a subroutine script now, so we can call it from multiple
% versions of the simulator code without dealing with inconsistencies between LUT conditions.  If you are working with
% multiple versions of LUT/textures for lightsail designs, it is generally easiest to load them all here, and assign
% whichever textures are used in the mesher code.  Aside from load-time, there is little computational penalty for 
% loading more LUTs than are actually used in your lightsail.  
%
% See LoadLUT_v17() for thorough document of LUT format and behavior.  


 
    LUTPsiStepDeg = 0.20;
    LUTThetaStepDeg = 0.20;
    showLUTPlots = 0;
    
    %clear('LUTs');  %LUTs is going to be a structure array for all look-up
    % tables in all regions.  That way we can easily add or remove patterned
    % regions in future designs.
    
    disp('Loading look-up talbles');
    
    if strcmpi(material,'silicon')
        % I have not yet tested silicon loading, but this block should work:
        file_name_collection_small_tables = ...
            '03212021_MEPS_SOI_Oggy_ReducedPitchRoll_Pressures_CollectionTables_';
        
        LUTs(1) = LoadLUT_v17(  ...
            [fullfile('Data',file_name_collection_small_tables) 'TE.mat'], ...
            'pressure_TE_small_tables', LUTPsiStepDeg, LUTThetaStepDeg, showLUTPlots, 0);
        
        LUTs(2) = LoadLUT_v17(  ...
            [fullfile('Data',file_name_collection_small_tables) 'TM.mat'], ...
            'pressure_TM_small_tables', LUTPsiStepDeg, LUTThetaStepDeg, showLUTPlots, 1);
        
    elseif strcmpi(material,'nitride')
        
        file_name_collection_small_tables = ...
          '03222021_MEPS_SiNx_Mark6e1_ReducedPitchRoll_Pressures_CollectionTables_';
        
        LUTs(1) = LoadLUT_v17(  ...
            [fullfile('Data',file_name_collection_small_tables) 'TE.mat'], ...
            'pressure_TE_small_tables', LUTPsiStepDeg, LUTThetaStepDeg, showLUTPlots,0);
        if (showLUTPlots)
            saveas(gcf,[filebasename '_LUT1.fig']);
            saveas(gcf,[filebasename '_LUT1.png']);
        end
        
        %warning('I am artificially increasing LUT(1).t!!!!');
        %LUTs(1).t = 10 * LUTs(1).t;
        
        LUTs(2) = LoadLUT_v17( ...
            [fullfile('Data',file_name_collection_small_tables) 'TM.mat'], ...
            'pressure_TM_small_tables', LUTPsiStepDeg, LUTThetaStepDeg, showLUTPlots,1);
        if showLUTPlots
            saveas(gcf,[filebasename '_LUT2.fig']);
            saveas(gcf,[filebasename '_LUT2.png']);
        end
        
        file_name_collection_small_tables = ...
             '08232021_MEPS_Si3N4_Mark1e_1064_PitchRoll_Pressures_CollectionTables_';
        %    '08222021_MEPS_Si3N4_Mark1c_1064_PitchRoll_Pressures_CollectionTables_';
        %    '04052021_MEPS_SiNx_Mark6e2_ReducedPitchRoll_Pressures_CollectionTables_';
        %    '03222021_MEPS_SiNx_Mark6e1_ReducedPitchRoll_Pressures_CollectionTables_';
        
        LUTs(3) = LoadLUT_v17(  ...
            [fullfile('Data',file_name_collection_small_tables) 'TE.mat'], ...
            'pressure_TE_small_tables', LUTPsiStepDeg, LUTThetaStepDeg, showLUTPlots,0);
        if (showLUTPlots)
            saveas(gcf,[filebasename '_LUT1.fig']);
            saveas(gcf,[filebasename '_LUT1.png']);
        end
        
        %warning('I am artificially increasing LUT(1).t!!!!');
        %LUTs(1).t = 10 * LUTs(1).t;
        
        LUTs(4) = LoadLUT_v17( ...
            [fullfile('Data',file_name_collection_small_tables) 'TM.mat'], ...
            'pressure_TM_small_tables', LUTPsiStepDeg, LUTThetaStepDeg, showLUTPlots,1);
        if showLUTPlots
            saveas(gcf,[filebasename '_LUT2.fig']);
            saveas(gcf,[filebasename '_LUT2.png']);
        end
        
    end