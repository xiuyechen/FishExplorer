% Instant Load!

% change current folder to the folder containing this m-file

clear all;close all;clc

global VAR;
load('VAR_current.mat','VAR');

load('CONSTs_current.mat','CONSTs');

data_dir = [pwd '\example data']; % 'F:\Janelia2014';%

%% Start GUI

[hfig,fcns] = GUI_FishExplorer(data_dir,CONSTs,VAR);

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
