clear all; close all; clc

%% folder setup
saveFigFlag = 1;

outputDir = GetOutputDataDir;
saveDir0 = fullfile(outputDir,'multimotor_allcells_0228');
if ~exist(saveDir0, 'dir'), mkdir(saveDir0), end;

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

%% set params!
isCellbased = true;
stimflag = 'P';

if isCellbased
    ClusterIDs = [2,1];
else
    ClusterIDs = [6,2];
end

% stimflag = []; % for default set
% ClusterIDs = [6,1];
M_stimrange = GetStimRange(stimflag);

range_fish_valid = [];
for i_fish = range_fish
    if ~isempty(M_stimrange{i_fish})
        range_fish_valid = [range_fish_valid,i_fish]; %#ok<AGROW>
    end
end

tscriptstart = tic;
for i_fish = range_fish_valid
    disp(['i_fish = ',num2str(i_fish)]);
    
    %% load data for chosen stim range
    stimrange = M_stimrange{i_fish};
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    
    %% get motor regressors
%     setappdata(hfig,'isMotorseed',0);
%     [~,~,behavior] = UpdateTimeIndex(hfig);
    [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);
    
    %% Method 1: stim regs + motor regs
    setID = 1;
    
    fishset = getappdata(hfig,'fishset');
    [~,~,reg_sens] = GetStimRegressor(stim,fishset,i_fish);
%     [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);
    
    if isCellbased
        gIX = (1:length(cIX_load))';
        Data = M;
        topDir = fullfile(saveDir0,'stimregs - cell based (new cmap)');
    else % cluster based
        gIX = gIX_load;
        C = FindClustermeans(gIX,M);
        Data = C;
        topDir = fullfile(saveDir0,'stimregs - cluster based');
    end

    [stimcorr,motorcorr] = MotorSourceCorrelation(Data,reg_sens,reg_motor);
    
    % make figures
    M_figs = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX_load,gIX,[1:5],setID);
    f = combineFiguresLR([M_figs{:}]);

    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, topDir, figName,f);    
    
    %% Method 2: stimAvr + motor regs
    setID = 2;
    
    if isCellbased
        gIX = (1:length(cIX_load))';
        Data = M;
        bottomDir = fullfile(saveDir0,'stimAvr - cell based (new cmap)');
    else % cluster based
        gIX = gIX_load;
        C = FindClustermeans(gIX,M);
        Data = C;
        bottomDir = fullfile(saveDir0,'stimAvr - cluster based');
    end

    reg_sens = GetTrialAvrLongTrace(hfig,Data);
%     [~,~,reg_motor] = GetMotorRegressor(behavior,i_fish);
    
    nUnit = size(Data,1);
    stimcorr = zeros(nUnit,1);
    motorcorr = zeros(nUnit,1);
    for i = 1:nUnit
        [stimcorr(i),motorcorr(i)] = MotorSourceCorrelation(Data(i,:),reg_sens(i,:),reg_motor);
    end    
    
    % make figures
    M_figs = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX_load,gIX,[1:5],setID);
    f = combineFiguresLR([M_figs{:}]);
    
    figName = ['Fish' num2str(i_fish)];
    SaveFigureHelper(saveFigFlag, bottomDir, figName,f);        
    
end
toc(tscriptstart)

%% compare params
% outputDir = GetOutputDataDir;
% saveDir0 = fullfile(outputDir,'multimotor_0228');

% topDir = fullfile(saveDir0,'stimregs - cell based');
% bottomDir  = fullfile(saveDir0,'stimAvr - cell based');

newDir = fullfile(saveDir0,'stimregs vs stimAvr (cell based, new cmap)');

compareFoldersTB(topDir,bottomDir,newDir);

