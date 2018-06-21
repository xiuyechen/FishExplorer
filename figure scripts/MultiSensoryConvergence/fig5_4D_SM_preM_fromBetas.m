% key function: [C_trialAvr,C_trialRes,C_score,C_d2var_perstim] = GetTrialAvrLongTrace(hfig,C)

% right now: alternate i_lr = 1 or 2 manually for left/right motor side
% options: compute on single cell or AutoClus



clear all; close all; clc

%% folder setup
outputDir = GetOutputDataDir;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];%

%% set params!
% isCellbased = true;
ClusterIDs = [2,1];
% ClusterIDs = [7,1];

%%
tscriptstart = tic;
nSets = 1;
IM_1 = cell(nSets,18);
% IM_2 = cell(nSets,18);

M_reg_name = {'preM_0.02_0.1_0.1'};%{'SMT_vs_MO','SMT_vs_SO','MO_vs_SMT','SO_vs_SMT'};

%%
for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);
    
    %% load data
    [cIX_load,gIX_load] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,[],0);

    
    load(fullfile(outputDir,'4D_SM_betas.mat'));
    %%
    CIX = cell(1,2);
    GIX = cell(1,2);
    
    for i_lr = 1:2
        betas = Betas{i_lr,i_fish};
        b1 = betas(:,1);
        b2 = betas(:,2);
        b3 = betas(:,3);
        
        cIX_in = cIX_load;
        gIX_in = ones(size(cIX_in))*i_lr; %(1:length(cIX_load))'; % i.e. for loading all cells, gIX = cIX
        
        %% setdiff plots for 1D cutoff
        
        switch 1
            case 1 % works for preM fish8R
                [cIX_b2,gIX_b2,IX_b2] = thresY(b2,cIX_in,gIX_in,0.02);
                [cIX_b1,gIX_b1,IX_b1] = thresY(b1,cIX_in,gIX_in,0.1);
                [cIX_b3,gIX_b3,IX_b3] = thresY(b3,cIX_in,gIX_in,0.1);
                
                IX_pass = setdiff(IX_b2,union(IX_b1,IX_b3));
                
            case 2
                IX_pass = intersect(IX_b1,IX_b2);
            case 3
                [cIX_b1,gIX_b1,IX_b1] = thresY(b1,cIX_in,gIX_in,0.02);
                IX_pass = IX_b1;
            case 4
                [cIX_b3,gIX_b3,IX_b3] = thresY(b3,cIX_in,gIX_in,0.02);
                IX_pass = IX_b3;
            case 5
                [cIX_b2,gIX_b2,IX_b2] = thresY(b2,cIX_in,gIX_in,0.02);
                IX_pass = IX_b2;
        end
        
        CIX{i_lr} = cIX_in(IX_pass);
        GIX{i_lr} = gIX_in(IX_pass);
        
    end
    
    cIX_out = vertcat(CIX{1},CIX{2});
    gIX_out = vertcat(GIX{1},GIX{2});
    
    % get map color
%     clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres,1],64);
%     clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres,1],64);
%     %     clrIX1 = MapXto1Dcolormap(corr_max(cIX1),[reg_thres1,1],64);
%     %     clrIX2 = MapXto1Dcolormap(corr_max(cIX2),[reg_thres2,1],64);
%     
%     cIX = [cIX1;cIX2];
%     %     if isempty(cIX)
%     %         M_fishrange_im{i_set} = setdiff(M_fishrange_im{i_set} ,i_fish); %#ok<SAGROW>
%     %         continue;
%     %     end
%     
%     clrIX = [clrIX1;clrIX2];
%     gIX_offset = [ones(size(cIX1));2*ones(size(cIX2))];
%     gIX = clrIX+(gIX_offset-1)*64;
    
    % make double colormap
%     clr1 = [1,0,0];
%     clr1_ = [0.7,0.5,0.5];
%     clr2 = [0,1,1];
%     clr2_ = [0.5,0.7,0.7];
%     numC = 64;
%     clrmap1 = Make1DColormap([clr1_;clr1],numC);
%     clrmap2 = Make1DColormap([clr2_;clr2],numC);
%     clrmap = [clrmap1;clrmap2];
    
    %%
    clrmap = [1,0,0; 0,1,1];
    I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap);    
    [h,~,im] = DrawCellsOnAnat(I);
    %%
    close(h)
    IM_1{1,i_fish} = im;
        
    %%
    
%     figure('Position',[50,100,400,500]);
%     % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
%     setappdata(hfig,'isPlotBehavior',1);
%     setappdata(hfig,'isStimAvr',1);
%     %         setappdata(hfig,'stimrange',[1:2]);
%     UpdateTimeIndex(hfig);
%     DrawTimeSeries(hfig,cIX_out,ones(size(cIX_out)));
    
end
toc(tscriptstart)

% M_reg_name = {'bubbleplot_R_all_passy','multimotor_R_anat_all','multimotor_R_anat_all_passy'};
% M_reg_name = {'bubbleplot_L_all_passy','multimotor_L_anat_all','multimotor_L_anat_all_passy'};
% M_reg_name = {'bubbleplot_R_A0.5_passy','multimotor_R_anat_A0.5','multimotor_R_anat_A0.5_passy'};
% M_reg_name = {'bubbleplot_R_all','multimotor_R_anat_all','multimotor_R_anat_all_passangle'};
% M_reg_name = {'bubbleplot_R','multimotor_R_anat_allA0.5','multimotor_R_anat_passangle'};
%% save as tiff stack
% for i_set = 1:nSets
i_set = 1;
range_im = range_fish;%M_fishrange_im{i_set};

tiffdir = fullfile(outputDir,['4D_',M_reg_name{i_set},'_allfish.tiff']);
IM = IM_1(i_set,range_im);
SaveImToTiffStack(IM,tiffdir);

% end

%% Average Plot
% M_k_scale = {1,1.5,1};
% M_k_contrast = {1.2,1.5,1.2};

% optional override::::
% M_fishrange_im{1} = [1,3,5:17];
% M_fishrange_im{1} = 8:17;% for OMR

i_set = 1;
% for i_set = 1:nSets
    range_im = [1:3,5:10];%range_fish;% M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_1(i_set,range_im);
    
    %% adjust params for visualization
    k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.1;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    %%
    tiffdir = fullfile(outputDir,['4D_',M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    
    %%
%     cellarray = IM_2(i_set,range_im);
%     
%     %% adjust params for visualization
%     k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
%     k_contrast = 1.1;%M_k_contrast{i_set};
%     
%     [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
%     %%
%     tiffdir = fullfile(outputDir,['anat_',M_reg_name{i_set},'_avr.tiff']);
%     imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
% end

