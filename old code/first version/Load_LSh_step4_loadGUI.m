% Instant Load!

clear all;close all;clc

global VAR;
load('VAR_current.mat','VAR');

load('CONSTs_current.mat','CONSTs');

data_dir = 'F:\Janelia2014';

%% Start GUI

[hfig,fcns] = GUI_FishExplorer(data_dir,CONSTs,VAR);

% hfig = GUI_FishExplorer(data_dir,CONSTs,VAR,flag_script,var_script)

f = [];
for i = 1:length(fcns),
    eval(['f.' char(fcns{i}) ' = fcns{i};']);
end

%%
% to recover after closing figure:
% global EXPORT_autorecover;
