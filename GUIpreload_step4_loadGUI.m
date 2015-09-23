% Instant Load!

% change current folder to the folder containing this m-file

clear all;close all;clc

data_dir = GetCurrentDataDir();%'C:\Janelia2014'; %[pwd '\example data']; % 'F:\Janelia2014';%

global VAR;

load(fullfile(data_dir,'VAR_current.mat'),'VAR'); % stores all clustering indices
% save(fullfile(data_dir,'VAR_current.mat'),'VAR','-v6');

name_CONSTs = 'CONSTs_current2.mat'; % stores selected data from all fish

% Example usage of name_CONSTs, used in 'Quickload' in GUI:
% matObj = matfile(fullfile(data_dir,'CONSTs_current2.mat'));
% i_fish = 1;
% eval(['CONST_s = matObj.CONST',num2str(i_fish),';']); 

% MASKs = load(fullfile(data_dir,'MaskDatabase.mat')); % loading this in GUI
name_MASKs = 'MaskDatabase.mat'; % stores selected data from all fish

%% Start GUI

[hfig,fcns] = GUI_FishExplorer(data_dir,name_CONSTs,name_MASKs,VAR);

% to push new functions to GUI, called with function 'push_cIX_gIX':
% hfig = GUI_FishExplorer(data_dir,CONSTs,VAR,flag_script,var_script)

% to get handle of local functions from GUI to use in workspace 
% (call by using 'f.FunctionName')
f = [];
for i = 1:length(fcns),
    eval(['f.' char(fcns{i}) ' = fcns{i};']);
end

%%
% to recover after closing figure:
% global EXPORT_autorecover;
