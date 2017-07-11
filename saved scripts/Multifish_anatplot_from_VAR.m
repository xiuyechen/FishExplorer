clear all; close all; clc;

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% params
for caseflag = 6
    switch caseflag
        case 1 % seeds
            stimrange = [];
            M_clus_range = {[11,1],[12,1]};
            M_clus_name = {'motorseed','ABN'};
            M_fishrange = {[1:3,5:18],[1:12,14:18]};
            n_reg = length(M_clus_name);
            
        case 2 % testing
            stimrange = [];
            M_clus_range = {[11,2],[11,3],[11,4]};
            M_clus_name = {'motorseed_v2','motorseed_LRmanual','motorseed_Fwmanual'};
            M_fishrange = {[5:13],[5:13],[5:13]};
            n_reg = length(M_clus_name);

        case 3 % conserved clusters
            stimrange = 1;
            M_clus_range = {[10,1]};
            M_clus_name = {'conservedclus_SLranked'};
            M_fishrange = {1:18};
            n_reg = length(M_clus_name);
            isStimLockRanking = 1;
            
        case 4 % all Auto0.7 clusters
            stimrange = [];
            M_clus_range = {[6,1]};
            M_clus_name = {'A0.7_def_stimlockrank'};
            M_fishrange = {1:18};
            n_reg = length(M_clus_name);
            isStimLockRanking = 1;
            
        case 5 % all Auto0.7 clusters
            stimrange = [];
            M_clus_range = {[6,1]};
            M_clus_name = {'A0.7_def_motorrank'};
            M_fishrange = {[1:3,5:18]};
            n_reg = length(M_clus_name);
            isStimLockRanking = 0;
            
        case 6 % HBO stripes
            stimrange = 2;
            M_clus_range = {[10,4]};
            M_clus_name = {'HBO_OMR_stimlockrank'};
            M_fishrange = {8:18};%{1:18};%{[1:3,5:18]};
            n_reg = length(M_clus_name);
            isStimLockRanking = 1;
            
            P_stimscore = zeros(18,4);
            
%                     case 6 % HBO stripes
%             stimrange = 1;
%             M_clus_range = {[10,4]};
%             M_clus_name = {'HBO_PT_stimlockrank'};
%             M_fishrange = {1:18};%{1:18};%{[1:3,5:18]};
%             n_reg = length(M_clus_name);
%             isStimLockRanking = 1;
%             
%             P_stimscore = zeros(18,4);
%             
%         case 3 % Looming
%             stimrange = 5;
%             M_stimmotorflag = [1,0]; % 1 for stim and 0 for motor
%             M_reg_name = {'Loom_LR','Loom_SwimLR'};
%             M_reg_range =  {[11,12],[1,3]};
%             M_fishrange = {[9:15,17:18],[9:15,17:18]};
%             n_reg = length(M_reg_name);
%             
%         case 4 % Dark Flash (black-white)
%             stimrange = 3;
%             M_stimmotorflag = [1,0]; % 1 for stim and 0 for motor
%             M_reg_name = {'DF_BW','DF_SwimLR'};
%             M_reg_range =  {[1,4],[1,3]};
%             M_fishrange = {[12:15,17:18],[12:15,17:18]};
%             n_reg = length(M_reg_name);
    end
    
    M_fishrange_im = M_fishrange;
    
    %% Load fish
    range_fish = 1:18;
    IM_full = cell(n_reg,max(range_fish));
    %%
    for i_fish = range_fish
        %% load fish
        if ismember(i_fish,cell2mat(M_fishrange))
            ClusterIDs = [1,1];
            if isempty(stimrange)
                LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
            else
                LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
            end
            %         [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
        end
        
        %% load clusters
        for i_set = 1:n_reg
            range_set = M_fishrange{i_set};
            if ismember(i_fish,range_set)
                ClusterIDs = M_clus_range{i_set};%[12,1];% GetClusterIDs('all');
                [cIX,gIX] = LoadCluster_Direct(i_fish,ClusterIDs(1),ClusterIDs(2));
                
%                 % override:
%                 if caseflag==6
%                     [cIX1,gIX1] = LoadCluster_Direct(i_fish,10,2);
%                     [cIX2,gIX2] = LoadCluster_Direct(i_fish,10,3);
%                     cIX = [cIX1;cIX2];
%                     gIX = [gIX1;5-gIX2];
                    
%                     absIX = getappdata(hfig,'absIX');
%                     name = 'HBO_4stripes';
%                     clusgroupID = 10;
%                     clusIDoverride = 4;
%                     SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
%                 end
                
                if isempty(cIX)
                    M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish);
                    continue;
                end
                
                %% [option: ranking]
                if caseflag>=3
                    if isStimLockRanking
                        M = UpdateIndices_Manual(hfig,cIX,gIX);
                        C = FindClustermeans(gIX,M);
                        [~,~,H] = GetTrialAvrLongTrace(hfig,C);
                        [gIX,rankscore] = SortGroupIXbyScore(H,gIX);
                        
                        if caseflag==6
                            P_stimscore(i_fish,:) = H;
                        end
                    else
                        M = UpdateIndices_Manual(hfig,cIX,gIX);
                        C = FindClustermeans(gIX,M);
                        numU = max(gIX);
                       [gIX,rankscore] = RankByMotorReg_Direct(hfig,gIX,numU,C,1); 
                    end
                end
                
                %% make figure
                
                if caseflag==3 || caseflag==4
                    setappdata(hfig,'clrmap_name','hsv_new');
                elseif caseflag==6
                    setappdata(hfig,'clrmap_name','jet');
                else
                    setappdata(hfig,'clrmap_name','hsv_old');
                end
                
                I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
                [h,im_full] = DrawCellsOnAnat(I);
                
                %% save figure
                close(h);
                IM_full{i_set,i_fish} = im_full;
            end
        end
    end
    
%     if caseflag==6
%         % save VAR
%         SaveVARtoMat(hfig);
%     end
    
    %% save as tiff stack
    for i_set = 1:n_reg
        range_im = M_fishrange_im{i_set};
        tiffdir = fullfile(outputDir,[M_clus_name{i_set},'_allfish.tiff']);
        IM = IM_full(i_set,range_im);
        
        SaveImToTiffStack(IM,tiffdir);
    end
    
    %% [for later] plot from tiff stack
    isPlotfromtiffstack = 0;
    
    if isPlotfromtiffstack
        IM_full = cell(n_reg,18);
        for i_set = 1:n_reg
            %% get tiff-stack path
            tiffdir = fullfile(outputDir,[M_clus_name{i_set},'_allfish.tiff']);
            
            %% load
            
            for i_fish = 1:18
                im = double(imread(tiffdir,i_fish))./255;
                IM_full{i_set,i_fish} = im(317:1236,1:621,:);
            end
        end
    end
    
    %% save average tiff image
    for i_set = 1:n_reg;
        range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
        cellarray = IM_full(i_set,range_im);
        
        % adjust params for visualization
        k_scale = 0.5;
        k_contrast = 1;% 1.2;
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_clus_name{i_set},'_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
    
end