% LoadGUI

% change current folder to the folder containing this m-file

clear all;close all;clc
%%
% scriptName = mfilename('fullpath');
% [currentpath, filename, fileextension]= fileparts(scriptName);
% code_dir = currentpath;
code_dir = 'C:\Users\Xiu\Dropbox (Personal)\FishExplorer2';
addpath(genpath(code_dir));

% check that main data directory is correct in function'GetCurrentDataDir'
data_masterdir = GetCurrentDataDir();%'C:\Janelia2014'; %[pwd '\example data']; % 'F:\Janelia2014';%
disp(data_masterdir);
global VAR; %#ok<NUSED>
load(fullfile(data_masterdir,'VAR_new.mat'),'VAR'); % stores all clustering indices

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
