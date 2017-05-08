% LoadGUI

% manually check that directories are correct in 'FishExplorer2\dir setup'

clear all;close all;clc
% path(pathdef);

%%
% scriptName = mfilename('fullpath');
% [currentpath, filename, fileextension]= fileparts(scriptName);
% code_dir = currentpath;

%% setup dir
% file_masterdir = GetMasterFileDir;
file_masterdir = 'C:\Users\Xiu\Dropbox (Personal)';
if ~exist(file_masterdir,'dir')
    file_masterdir = 'C:\Users\fish\Dropbox';
end

code_dir = fullfile(file_masterdir,'FishExplorer2');
addpath(genpath(code_dir));

cd(code_dir)

%% setup
global VAR; %#ok<NUSED>
load_dir = GetOutputDataDir();
load(fullfile(load_dir,'VAR_new.mat'),'VAR'); % stores all clustering indices
VAR = CompileVARnames(VAR);

set(0, 'defaultUicontrolFontName', 'Arial');
set(0, 'defaultUitableFontName', 'Arial');
set(0, 'defaultUipanelFontName', 'Arial');
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');

%% Start GUI
% Run this, or press F5 in 'GUI_FishExplorer.m'
[hGUI,fcns] = GUI_FishExplorer;

% to push new functions to GUI, called with function 'push_cIX_gIX':
% hfig = GUI_FishExplorer(data_dir,name_CONSTs,name_MASKs,VAR,flag_script,var_script)

%%
% to get handle of local functions from GUI to use in workspace
% (call by using 'f.FunctionName')
f = [];
for i = 1:length(fcns),
    eval(['f.' char(fcns{i}) ' = fcns{i};']);
end

%%
% to recover after closing figure:
% global EXPORT_autorecover;

%% demo: load data without visualizaing in GUI (only using same infrastructure)
if false,
    %%
    h0 = figure;
    InitializeAppData(h0);
    i_fish = 3;
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,h0);
    
    figure('Position',[50,100,800,1000]);
    I = LoadCurrentFishForAnatPlot(h0);
    DrawCellsOnAnat(I);
    
    % filename = fullfile(saveDir, num2str(i_fish));
    % saveas(gcf, filename, 'png');
end

%%
% MASKs = getappdata(hGUI,'MASKs');
% Masks = MASKs.MaskDatabaseNames';