%%
clear all; close all; clc

%% set MATLAB path
% set the Explore2p code folder as current directory
if true % set the
    scriptName = mfilename('fullpath');
    [code_dir, filename, fileextension]= fileparts(scriptName);
    addpath(genpath(code_dir));
    cd(code_dir)
    
else
    % (hard-coded, customize here)
    folder = 'C:\Users\xiuye\Dropbox\!Research\2Pcode\Explore2p'; %#ok<UNRCH>
    addpath(genpath(folder));
    cd(folder)
end

%% Automatically load demo data
% isDemo = true;
%     
% %% Start the GUI
% h = Explore2p(isDemo);

% to update the object handle in the workspace, go to the GUI Menu and press 
% File\Export to workspace. See demo\demo_workspace_interactive.m

%% For Scripts: see demo\demo_script.m for using Explore2p in custom scripts

%% To load the GUI with non-demo data, simply call -
% Explore2p;