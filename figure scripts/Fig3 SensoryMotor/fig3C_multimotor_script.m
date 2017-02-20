clear all; close all; clc

%% folder setup
saveFigFlag = 1;

outputDir = GetOutputDataDir;
saveDir1 = fullfile(outputDir,'multimotor_0216_motorseed_fixedscale_stimrangeP');
% saveDir2 = fullfile(outputDir,'motor_map_seed');
if ~exist(saveDir1, 'dir'), mkdir(saveDir1), end;
% if ~exist(saveDir2, 'dir'), mkdir(saveDir2), end;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
% i_fish = 8;
% setappdata(hfig,'isMotorseed',0);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];
% M_thres_reg = zeros(3,18);
% M_numTopCorr = zeros(1,18);
% M_motorseedRegs = cell(1,18);
% M_compareMotorCellNumber = zeros(2,18);
        
%%
stimflag = 'P';
ClusterIDs = [6,2];

% stimflag = []; % for default set
% ClusterIDs = [6,1];
M_stimrange = GetStimRange(stimflag);


tscriptstart = tic;
for i_fish = range_fish
    %% define stim range
    
    stimrange = M_stimrange{i_fish};
    if isempty(stimrange)
        continue;
    else
        [cIX,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    end
    
    %% Method 1: stim+motor regression
    fishset = getappdata(hfig,'fishset');
    [~,~,reg_sens] = GetStimRegressor(stim,fishset,i_fish);
    [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);

    % cell based
    gIX = (1:length(cIX))';
    tic
    [stimcorr,motorcorr] = MotorSourceCorrelation(M,reg_sens,reg_motor);
    toc
    [fig1,fig2] = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX,gIX);

    saveDir = fullfile(saveDir1,'motor rank - cell based - 2D');
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, saveDir, figName,fig1);
    
    saveDir = fullfile(saveDir1,'motor rank - cell based - anat');
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, saveDir, figName,fig2);
    
    %% cluster based
    % plot 1
%     setappdata(hfig,'isMotorseed',0);
%     %     setappdata(hfig,'stimrange',1:2);
%     
%     UpdateTimeIndex(hfig);
%     behavior = getappdata(hfig,'behavior');
%     [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);

    gIX = gIX_load;
    C = FindClustermeans(gIX,M);
    tic
    [stimcorr,motorcorr] = MotorSourceCorrelation(C,reg_sens,reg_motor);
    toc
    [fig1,fig2] = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX,gIX);
    
%     % plot 2
%     setappdata(hfig,'isMotorseed',1);
%     %     setappdata(hfig,'stimrange',1:2);
%     
%     UpdateTimeIndex(hfig);
%     behavior = getappdata(hfig,'behavior');
%     [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);
%     
%     gIX = gIX_load;
%     C = FindClustermeans(gIX,M);
%     tic
%     [stimcorr,motorcorr] = MotorSourceCorrelation(C,reg_sens,reg_motor);
%     toc
%     [fig3,fig4] = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX,gIX);
    %%
    saveDir = fullfile(saveDir1,'motor rank - cluster based - 2D');
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, saveDir, figName,fig1);
    
    saveDir = fullfile(saveDir1,'motor rank - cluster based - anat');
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, saveDir, figName,fig2);
    %%
    
    
% %  saveas(gcf, fn, 'png');    
% % 
% % %%
% %     [rawmotorcorr,IX] = max(corr(reg_motor',M'),[],1);
% %     thres_reg = 0.5;
% %     ix = find(rawmotorcorr>thres_reg);
    %% Method 2: stimAvr+motor regression
    reg_sens = GetStimAvrClusmean(hfig,gIX,M);
    [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);

    % cluster based
    gIX = gIX_load;
    C = FindClustermeans(gIX,M);
    tic
    [stimcorr,motorcorr] = MotorSourceCorrelation(C,reg_sens,reg_motor);
    toc
    [fig1,fig2] = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX,gIX);
    
    saveDir = fullfile(saveDir1,'stimAvr rank - cluster based - 2D');
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, saveDir, figName,fig1);
    
    saveDir = fullfile(saveDir1,'stimAvr rank - cluster based - anat');
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, saveDir, figName,fig2);
        
end
toc(tscriptstart)