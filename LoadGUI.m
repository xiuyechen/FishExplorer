% LoadGUI

% change current folder to the folder containing this m-file

clear all;close all;clc
%%
% scriptName = mfilename('fullpath');
% [currentpath, filename, fileextension]= fileparts(scriptName);
% code_dir = currentpath;
code_dir = 'C:\Users\Xiu\Dropbox\FishExplorer2';
addpath(genpath(code_dir));

global data_masterdir name_MASKs name_ReferenceBrain VAR; %#ok<NUSED>

data_masterdir = GetCurrentDataDir();%'C:\Janelia2014'; %[pwd '\example data']; % 'F:\Janelia2014';%

load(fullfile(data_masterdir,'VAR_new.mat'),'VAR'); % stores all clustering indices

% load(fullfile(data_masterdir,'VAR_current.mat'),'VAR'); % stores all clustering indices

% save(fullfile(data_dir,'VAR_current.mat'),'VAR','-v6');

% name_CONSTs = 'CONSTs_current2.mat'; % stores selected data from all fish

% Example usage of name_CONSTs, used in 'Quickload' in GUI:
% matObj = matfile(fullfile(data_dir,'CONSTs_current2.mat'));
% i_fish = 1;
% eval(['CONST_s = matObj.CONST',num2str(i_fish),';']); 

% MASKs = load(fullfile(data_dir,'MaskDatabase.mat')); % loading this in GUI
name_MASKs = 'MaskDatabase.mat';

name_ReferenceBrain = 'ReferenceBrain.mat';

set(0, 'defaultUicontrolFontName', 'Arial');
set(0, 'defaultUitableFontName', 'Arial');
set(0, 'defaultUipanelFontName', 'Arial');
set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');

%% Start GUI
% Run this, or press F5 in 'GUI_FishExplorer.m'
[hfig,fcns] = GUI_FishExplorer;

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
