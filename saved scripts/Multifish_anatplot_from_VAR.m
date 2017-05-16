clear all; close all; clc;

outputDir = GetOutputDataDir;

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% params
for caseflag = 2
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
                
                if isempty(cIX)
                    M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish);
                    continue;
                end
                
                %% make figure
                setappdata(hfig,'clrmap_name','hsv_old');
                I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
                [h,im_full] = DrawCellsOnAnat(I);
                
                %% save figure
                close(h);
                IM_full{i_set,i_fish} = im_full;
            end
        end
    end
    
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
    for i_set = 1:nreg;
        range_im = M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
        cellarray = IM_full(i_set,range_im);
        
        % adjust params for visualization
        k_scale = 1;
        k_contrast = 1;% 1.2;
        
        [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
        
        tiffdir = fullfile(outputDir,[M_clus_name{i_set},'_avr.tiff']);
        imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    end
    
end