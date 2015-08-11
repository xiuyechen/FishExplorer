%%% Interactive app for exploratory analysis of calcium imaging data 
% (with stimulus, behavior, and anatomy)

% Input calcium data: 1 trace per cell/ROI, ~50,000 cells per fish
% load collection of cells from multiple fish, or load full data of single fish individually
% main outputs: GUI plots, clusters saved into .mat, export variables to MATLAB workspace

% Tip: to see the structure of this code, use 'Ctrl' + 'm' + '=' to collapse all cells.
% UI controls are organized by tabs and then by rows, instructions and 
% comments are where they are constructed ('User Interface:' -> function hfig... ->)
% General internal functions are at the end, some specialized .m functions are outside.

% Written in Matlab R2014b running on Windows 7.

% - Xiuye Chen (xiuyechen@gmail.com), Engert Lab, Spring 2015


% To start, run the following in separate script or command window:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all;close all;clc
% global VAR;
% load('VAR_current.mat','VAR');
% load('CONSTs_current.mat','CONSTs');
% data_dir = [pwd '\example data'];
% [hfig,fcns] = GUI_FishExplorer(data_dir,CONSTs,VAR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full datasets are stored as 'CONST_F?_fast.mat' (? = fish number), to load from the GUI

%% User Interface:
% function hfig = GUI_FishExplorer(CONSTs,VAR,CONST,i_fish)
function [hfig,fcns] = GUI_FishExplorer(data_dir,CONSTs,VAR,flag_script,var_script)
if exist('flag_script','var')
    if ~exist('var_script','var')
        var_script={};
    end
    runscript(flag_script,var_script);
    return
end

%% Make figure
scrn = get(0,'Screensize');
hfig = figure('Position',[scrn(3)*0.1 scrn(4)*0.05 scrn(3)*0.85 scrn(4)*0.86],...% [50 100 1700 900]
    'Name','GUI_LSh','DeleteFcn',@closefigure_Callback);
hold off; axis off

%% Folder setup
% directory for full fish data (.mat)
setappdata(hfig,'data_dir',data_dir);

% copy of VAR files will be saved into this subfolder:
currentfolder = pwd;
arcmatfolder = [currentfolder '\arc mat'];
if ~exist(arcmatfolder, 'dir')
  mkdir(arcmatfolder);
end
setappdata(hfig,'arcmatfolder',arcmatfolder);

%% Pass external variables into appdata (stored with main figure handle)
setappdata(hfig,'CONSTs',CONSTs);
setappdata(hfig,'VAR',VAR);
nFish = length(CONSTs) + 1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fish protocol sets (different sets have different parameters)
M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 3, 3]; % M = Matrix
setappdata(hfig,'M_fish_set',M_fish_set);

% parameters / constants
setappdata(hfig,'z_res',19.7); % resoltion in z, um per slice
% fpsec = 1.97; % hard-coded in ext function 'GetStimRegressor.m'
% approx fpsec of 2 used in ext function 'DrawCluster.m'

% cache 
bC = []; % Cache for going b-ack (bad abbr)
fC = []; % Cache for going f-orward
bC.cIX = cell(1,1);
bC.gIX = cell(1,1);
bC.numK = cell(1,1);
fC.cIX = cell(1,1);
fC.gIX = cell(1,1);
fC.numK = cell(1,1);
setappdata(hfig,'bCache',bC);
setappdata(hfig,'fCache',fC);

% initialization
i_fish = 1;
QuickUpdateFish(hfig,i_fish,'init');

%% Initialize internal params into appdata

% thresholds
thres_merge = 0.9;
thres_split = 0.7;
thres_reg = 0.7;
thres_size = 10;
setappdata(hfig,'thres_merge',thres_merge);
setappdata(hfig,'thres_split',thres_split); % for function 'pushbutton_iter_split'
setappdata(hfig,'thres_reg',thres_reg); % regression threshold, ~correlation coeff
setappdata(hfig,'thres_size',thres_size); % min size for clusters

% variables 
% (not sure all these need to be initialized, probably not complete either)
setappdata(hfig,'clrmap','hsv');
setappdata(hfig,'opID',0);
setappdata(hfig,'rankID',0); 
setappdata(hfig,'isPlotLines',0); 
setappdata(hfig,'isPlotFictive',1); 
setappdata(hfig,'rankscore',[]);
setappdata(hfig,'isCentroid',0);
setappdata(hfig,'isWkmeans',1); % in autoclustering, with/without kmeans
setappdata(hfig,'isflipstim',0); % flag for flip updown
setappdata(hfig,'regchoice',{1,1}); % regressor choice; first cell,1=stim,2=motor,3=centroid
setappdata(hfig,'isfullfish',0); % no if QuickUpdateFish, yes if LoadFullFish
setappdata(hfig,'isPlotCorrHist',0); % option for regression
setappdata(hfig,'dataFR',1); % left plot data type, 1 for fluo trace and 0 for regression result
setappdata(hfig,'isPlotReg',1); % plot regressor when selecting it
setappdata(hfig,'hierinplace',1); % hier. partitioning, no reordering


Class = getappdata(hfig,'Class');
classID = numel(Class);
setappdata(hfig,'classID',classID);
setappdata(hfig,'classheader','Test:');
setappdata(hfig,'newclassname','');

Cluster = getappdata(hfig,'Cluster');
clusID = 1; %numel(Cluster);
setappdata(hfig,'clusID',clusID);
setappdata(hfig,'clusheader','Test:');
setappdata(hfig,'newclusname','');

setappdata(hfig,'cIX',Cluster(clusID).cIX); % cell-IndeX
setappdata(hfig,'gIX',Cluster(clusID).gIX); % group-IndeX

%% Create UI controls
set(gcf,'DefaultUicontrolUnits','normalized');
set(gcf,'defaultUicontrolBackgroundColor',[1 1 1]);

% tab group setup
tgroup = uitabgroup('Parent', hfig, 'Position', [0.05,0.88,0.91,0.12]);
numtabs = 5;
tab = cell(1,numtabs);
M_names = {'General','Operations','Regression','Clustering etc.','Saved Clusters'};
for i = 1:numtabs,
    tab{i} = uitab('Parent', tgroup, 'BackgroundColor', [1,1,1], 'Title', M_names{i});
end

% grid setup, to help align display elements
rheight = 0.2;
yrow = 0.7:-0.33:0;%0.97:-0.03:0.88;
bwidth = 0.03;
grid = 0:bwidth+0.001:1;

% various UI element handles ('global' easier than passing around..)
global hback hfwd hclassname hclassmenu hclusgroupmenu hclusmenu hclusname...
    hdatamenu hopID hdataFR hloadfish hfishnum hstimreg hmotorreg...
    hcentroidreg;

%% UI ----- tab one ----- (General)
i_tab = 1;

%% UI row 1: File 
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2; % saves 'VAR' to workspace. 'VAR' contains clustering indices (Class and Clusgroup) for all fish
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Quick save',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_save_Callback);

i=i+n;
n=2; % saves 'VAR' both to workspace and to 'VAR_current.mat' and to arc folder
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Save .mat',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_savemat_Callback);

i=i+n;
n=2; % plots selected cells on anatomy z-stack, display and save tiff stack in current directory
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Save Zstack',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_writeZstack_Callback);

i=i+n;
n=2; % plots selected cells on anatomy z-stack, tiled display
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Tile Zstack',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_tileZstack_Callback);

i=i+n;
n=2; % make new figure without the GUI components, can save manually from default menu
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Popup plot',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_popupplot_Callback);

i=i+n;
n=2; % popupplot option: whether to plot cluster mean lines instead of all raw traces
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot lines',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_isPlotLines_Callback);

i=i+n;
n=2; % popupplot option: whether to plot fictive behavior bar
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot fictive','Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_isPlotFictive_Callback);

i=i+n;
n=3; % export main working variables to workspace, can customize!
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Export to workspace',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_exporttoworkspace_Callback);

%% UI row 2: Load 
i_row = 2;
i = 1;n = 0;

i=i+n;
n=2; % this design is underused now... Quick-load only depends on CONSTs,
% which is a minimum collection of clusters from all fish, so you can load
% the program without full single-fish data. eventually can use this  
% platform to do things across fish (like after anatomical alignment).
uicontrol('Parent',tab{i_tab},'Style','text','String','Quick-load fish:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1; % loads 'CONSTs_current.mat' from current directory
temp = {}; for j = 1:nFish, temp = [temp,{num2str(j)}];end
hfishnum = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_fishmenu_Callback);

i=i+n;
n=2; % loads full single-fish data from CONST_F?.mat
uicontrol('Parent',tab{i_tab},'Style','text','String','Load full fish:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2; % these CONST_F?.mat are now saved as uncompressed version 6 .mat, loads faster
temp = {}; for j = 1:nFish, temp = [temp,{num2str(j)}];end
temp = [{'(choose)'},temp];
hloadfish = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_loadfullfishmenu_Callback);

i=i+n;
n=2; % load different presentations (average, all reps etc), pre-stored in CONST(s)
uicontrol('Parent',tab{i_tab},'Style','text','String','Data type:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=3; 
hdatamenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',CONSTs{i_fish}.datanames,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_datamenu_Callback);

i=i+n;
n=4; % only centroids (~mean) of clusters shown on left-side plot, the rest is unchanged
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Display centroids of clusters',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_showcentroids_Callback);

%% UI ----- tab two ----- (Operations)
i_tab = 2;

%% UI row 1: Range
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2; % saves up to 20 steps backwards (datatype/datamenu change does not count)
hback = uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Back',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_back_Callback);

i=i+n;
n=2; % same, 20 steps forward if applicable
hfwd = uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Forward',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_forward_Callback);

i=i+n;
n=3; % Choose range of clusters to keep. format: e.g. '1:2,4-6,8:end'
uicontrol('Parent',tab{i_tab},'Style','text','String','Choose cluster range:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1; 
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_choose_range_Callback);

i=i+n;
n=2; % Choose range of clusters to exclude. format: e.g. '1:2,4-6,8:end'
uicontrol('Parent',tab{i_tab},'Style','text','String','Exclude:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_exclude_range_Callback);

i=i+n;
n=1; % Choose range of clusters to fuse/combine into single cluster. format: e.g. '1:2,4-6,8:end'
uicontrol('Parent',tab{i_tab},'Style','text','String','Fuse:',... % (eg 1:2,3-5)
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_fuse_range_Callback);

i=i+n;
n=2; % can choose range in time, but can't see time axes to decide... resets when choosing datamenu
uicontrol('Parent',tab{i_tab},'Style','text','String','Time-range:',... % (eg 1:2,3-5)
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_t_range_Callback);

%% UI row 2: Operations
i_row = 2;
i = 1;n = 0;

i=i+n;
n=2; % operates between the current cell selection and the next (in this order). 
uicontrol('Parent',tab{i_tab},'Style','text','String','Set operations:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n; % 'setdiff' is current minus next, 'rev setdiff' is next minus current. 
n=2; % smartUnion = SmartUnique, cells belonging to 2 clusters goes to the more correlated one
menu = {'(choose)','union','intersect','setdiff','rev setdiff','setxor','smartUnion'};
hopID = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@popup_operations_Callback});

i=i+n;
n=2; % rank clusters based on various criteria (to choose)
uicontrol('Parent',tab{i_tab},'Style','text','String','Rank by:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n; % 'hier' is the same as default (used after every k-means);'stim-lock' uses std across reps;
n=2; % motor stuff uses the best alignment (by cross-correlation) with the fictive trace;
% L+R is average of L & R; stim-motor is combines 'stim-lock' w 'motor' with arbituary weighting.
menu = {'(choose)','hier.','size','stim-lock','corr',...
    'motor','L motor','R motor','L+R motor','stim-motor'}; 
uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@popup_ranking_Callback});

i=i+n;
n=3; % cluster indices will rank from 1 to number of clusters
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Sqeeze clusters',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_sqeeze_Callback);

i=i+n;
n=3; % flip the sequenc of clusters, first becomes last
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Flip up-down',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_flipud_Callback);

i=i+n;
n=3; % switch between 2 colormaps now, jet and a cropped version of hsv (so not all circular)
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Switch colormap',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_clrmap_Callback});

%% UI row 3: Anatomy
i_row = 3;
i = 1;n = 0;

i=i+n;
n=4; % Draw a polygon on anatomy maps to select the cells within those boundaries
uicontrol('Parent',tab{i_tab},'Style','text','String','Draw on anatomy map to crop:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n; % Draw on the yx-view (main view) 
n=2; % Click to make new vertex, double click to connect to first vertex, 
% then optionally drag vertices to reposition, and finally double click again to set
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Draw yx',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_polygon_yx_Callback);

i=i+n;
n=2; % Draw on yz-view (side projection), same as above
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Draw yz',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_polygon_yz_Callback);

i=i+n;
n=2; % Draw on zx-view (front projection), same as above
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Draw zx',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_polygon_zx_Callback);

i=i+n;
n=5; 
uicontrol('Parent',tab{i_tab},'Style','text','String','Select all cells within boundaries:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2; % selects ALL cells contained in dataset that are within the convex shape defined by current cells
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Convex hull',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_withinConvexHull_Callback);

%% UI ----- tab three ----- (Regression)
i_tab = 3;

%% UI row 1: regressor
i_row = 1; % Step 1:
i = 1;n = 0; % Choose one type of regressor here, choice highlighted in yellow

i=i+n;
n=2; % stimulus regressors, go to 'GetStimRegressor.m' to add/update
uicontrol('Parent',tab{i_tab},'Style','text','String','Stim reg.:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=3; % (updated when loading fish)
menu = {'(choose)',''};
hstimreg = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_getstimreg_Callback);

i=i+n; % motor regressors from fictive, not yet convolved/adjusted for time lag
n=2; % go to 'GetMotorRegressor.m' to add/update
uicontrol('Parent',tab{i_tab},'Style','text','String','Motor reg.:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2; % (unlike stim regressors, names hardcoded, not importet from regressor...)
menu = {'(choose)','left swims','right swims','forward swims','raw left','raw right','raw average'};
hmotorreg = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_getmotorreg_Callback);

i=i+n;
n=2; % if checked, plot regressor during selection, together with stim and motor
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot regressor',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',1,...
    'Callback',@checkbox_popupplotreg_Callback);

i=i+n+1;
n=4; % Get centroid (~mean) of selected cluster as regressor
uicontrol('Parent',tab{i_tab},'Style','text','String','Regressor from centroid #:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
hcentroidreg = uicontrol('Parent',tab{i_tab},'Style','edit','String',num2str(1),...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_ctrdID_as_reg_Callback);

i=i+n;
n=2; % too complicated... chomp up centroid, then order as if stimulus were left/right inverted
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Flip stim',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_flipstim_Callback);

%% UI row 2: regression
i_row = 2; % Step 2:
i = 1;n = 0; % Choose regression, using the regressor chosen above, search in full dataset

i=i+n;
n=3;
uicontrol('Parent',tab{i_tab},'Style','text','String','Choose regression ->',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2; % do regression, show all cells with correlation coeff (with regressor) above threshold
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Corr. threshold:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_thres_regression_Callback);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit','String',num2str(thres_reg),...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_regthres_Callback);

i=i+n;
n=2; % optionally plot histogram of correlation values for all cells in dataset, visualize cut-off
hdataFR = uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot corr. hist',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_plotcorrhist_Callback);

i=i+n;
n=3; % probably not so useful. Show top n cells of highest correlation coeff with regressor
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','top corr., number limit:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_topnum_regression_Callback});

i=i+n;
n=1; % specify n for above
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_topCorrNumber_Callback);

i=i+n+1; % more automatic, do a regression with every centroid, then combine (with 'SmartUnique'
n=4; % i.e. overlapping cells are assigned to the cluster with which the correlation is the highest)
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Regression with all centroids',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_AllCentroidRegression_Callback});

i=i+n; % this is a remnant button from a failed experiment, idea was to iterate the regression process
n=2; % until the cluster converges, but most of the time it doesn't...
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','iter.reg','Enable','off',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_IterCentroidRegression_Callback});

%% UI row 3: (TBD)
% i_row = 3; % old idea: display correlation values of multiple regressors (instead of fluo. trace)
% i = 1;n = 0; % and then can cluster that and further manipulate... 
% 
% i=i+n;
% n=4; % did not implement combination of regressors in this version... but display option is coded (dataFR==0)
% uicontrol('Parent',tab{i_tab},'Style','text','String','regressor combos??(TBD)',... 
%     'Position',[grid(i) yrow(i_row) bwidth*n rheight]);
% 
% i=i+n;
% n=4; % also data storage is coded, 'M' could be loaded from M_0_fluo or M_0_reg (storing reg corr values)
% hdataFR = uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Display fluo./regression',...
%     'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
%     'Callback',@checkbox_dataFR_Callback);

%% UI ----- tab four ----- (Clustering etc.)
i_tab = 4;

%% UI row 1: k-means
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2; % k-means clustering
uicontrol('Parent',tab{i_tab},'Style','text','String','k-means, k =',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_kmeans_Callback);

i=i+n;
n=4; % anatomy is added to the fluo trace as new dimensions, and (arbituarily) weighted strongly
uicontrol('Parent',tab{i_tab},'Style','text','String','k-means with anatomy, k =',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_kmeans2_Callback);

i=i+n; % trying to use Silhouette to evaluate cluster quality, find peak to determine optimal k,
n=3; % then display results with that k. But have not set k-means to replicate (speed concern), can be very noisy
uicontrol('Parent',tab{i_tab},'Style','text','String','Find best k in range:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight])

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_kmeans_elbow_Callback);

%% UI row 2: Auto-clustering
i_row = 2;
i = 1;n = 0;

i=i+n; % Adjacent clusters (arranged by hier.) will be merged
n=4;  % if correlation between centroids is above merging threshold
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Merge thres. (corr. based)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_merge_Callback});

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit','String',num2str(thres_merge),...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_mergethres_Callback});

i=i+n; % further split initial clusters so that the average within-cluster corr coeff is above thres
n=2; % (not so iterative anymore, but could easily restore)
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Iter. split',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_iter_split});

i=i+n;
n=1; 
uicontrol('Parent',tab{i_tab},'Style','edit','String',num2str(thres_split),...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_splitthres_Callback});

i=i+n;
n=3; % minimal size of cluster, otherwise delete at the end
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Cluster size thres.',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_thressize_Callback});

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit','String',num2str(thres_size),...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_sizethres_Callback});

i=i+n+1; % longest script here. Splits clusters and prunes them, to yield only very tight clusters.
n=3; % really pretty results, but takes a while when regressing with every centroid. Read code for details.
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Full Auto-Clustering',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_autoclus_Callback});

i=i+n;
n=4; % by default it starts with a k-mean of 20 of the current cells. Could skip that if already clustered.
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','(starting with k-mean)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',1,...
    'Callback',@checkbox_wkmeans_Callback);

%% UI row 3: misc plots
i_row = 3;
i = 1;n = 0;

i=i+n; % k-means are always ranked like in hierachical clustering ~ optimal leaf order
n=2; % here you can rank them again and plot the dengrogram.
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Hier. plot',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_hierplot_Callback);

i=i+n; % partitioning based on hier. clustering
n=2; % choose between max cluster numbers...
uicontrol('Parent',tab{i_tab},'Style','text','String','Hier.cut, max n:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_hierpartn_Callback);

i=i+n; % partitioning based on hier. clustering
n=3; % ...or set correlation value threshold
uicontrol('Parent',tab{i_tab},'Style','text','String','Hier.cut, corr thres:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_hierpartthres_Callback);

i=i+n; % hier. partition in place, i.e. without rearranging order clusters
n=3; 
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Hier.cut in place',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',1,...
    'Callback',@checkbox_hierinplace_Callback);

i=i+n;
n=2; % Plots the correlation between all current clusters as a matrix
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Corr. plot',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_corrplot_Callback});

i=i+n;
n=2; % random plots not really worth consolidating, only for 1-click convenience
uicontrol('Parent',tab{i_tab},'Style','text','String','temp. plots:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight])

i=i+n;
n=3; % e.g. this one looks at the history for the fish presented with 4 B/W stimuli
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','plot 4x4(stim)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_plot4x4_Callback});

%% UI ----- tab five ----- (Saved Clusters)
i_tab = 5;

%% UI row 1: Class
i_row = 1; % There is not really a difference between Class and Cluster
i = 1;n = 0; % I just like to have 2 layers, so I save bigger sets as Class, smaller things as Cluster
% also now Clusters come in cluster-groups ('Clusgroup'), i.e. multiple groups of Clusters...
% the saving is sort of a pain though, several Update___ internal functions dedicated to this.

i=i+n;
n=1; % 
uicontrol('Parent',tab{i_tab},'Style','text','String','Class:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=4;
temp = MakeMenu({Class.name});
hclassmenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_classmenu_Callback);

i=i+n;
n=1; % just edit and press enter to edit current name
uicontrol('Parent',tab{i_tab},'Style','text','String','Edit:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=3;
hclassname = uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_editclassname_Callback});
set(hclassname,'String',Class(classID).name);

i=i+n;
n=2; % save the current cells, with current cluster divisions, as 1 class
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Save class',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_saveclass_Callback);

i=i+n;
n=1; % not really useful, just to help me organize the Class name
uicontrol('Parent',tab{i_tab},'Style','text','String','Header:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1; 
uicontrol('Parent',tab{i_tab},'Style','edit','String','Test',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_class_header_Callback);

i=i+n;
n=1; % write down new name for new cluster
uicontrol('Parent',tab{i_tab},'Style','text','String','New:',... 
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','edit','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_newclassname_Callback);

i=i+n;
n=2; % new class is made when pressing this button
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Make class',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_makeclass_Callback);

i=i+n;
n=2; % will prompt you to confirm
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Delete Class!',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_delclass_Callback});

%% UI row 2: Cluster
i_row = 2; % Virtually the same as Class, except added a ranking option
i = 1;n = 0;

i=i+n;
i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','Cluster (group):',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
num = length(VAR(i_fish).ClusGroup);
temp = {}; for j = 1:num, temp = [temp,{num2str(j)}];end
hclusgroupmenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_clusgroupmenu_Callback);

i=i+n;
n=3;
temp = MakeMenu({Cluster.name});
hclusmenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_clusmenu_Callback);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','Edit:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=3;
hclusname = uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_editclusname_Callback});
set(hclusname,'String',Cluster(clusID).name);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Save cluster',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_saveclus_Callback);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','Header:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit','String','Test',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_clus_header_Callback);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','New:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','edit','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_newclusname_Callback);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Make cluster',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_makeclus_Callback);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','Set rank:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_setrank_Callback});

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','Notes:',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_notes_Callback});

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Delete Cluster!',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_delclus_Callback});

%% UI row 3: misc
i_row = 3; % continuation, dealing with the Cluster groups 
i = 1;n = 0; % (Cluster-group number: number menu before the Cluster-name menu)

i=i+n;
n=2; % Combines the chosen clusters into one view (can save as 1 class or cluster then)
uicontrol('Parent',tab{i_tab},'Style','text','String','Union(cluster):',... % (eg 1,3-5)
    'Position',[grid(i) yrow(i_row) bwidth*n rheight]);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_clusUnion_Callback);

i=i+n; % just adds a new number to the Clustergroup-number menu, 
n=3; % and saves current view as the first cluster in the new Clustergroup
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','new Clustergroup',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_newclusgroup_Callback});

i=i+n;
n=3; % delete current
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','del Clustergroup',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_delclusgroup_Callback});

%% Load figure

UpdateClusID(hfig,clusID);

%% get local function handles

fcns = localfunctions;

end

%% Callback functions for UI elements: 

%% ----- tab one ----- (General)

%% row 1: File

function pushbutton_save_Callback(hObject,~)
hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');
global VAR;
VAR(i_fish).Class = getappdata(hfig,'Class');
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);
disp('saved to workspace ''VAR''');
end

function pushbutton_savemat_Callback(hObject,~)
disp('saving...');
hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');
arcmatfolder = getappdata(hfig,'arcmatfolder');
global VAR;
VAR(i_fish).Class = getappdata(hfig,'Class');
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);
timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
save(fullfile(arcmatfolder,matname),'VAR','-v6');
save('VAR_current.mat','VAR','-v6');
disp('saved both to workspace and .mat');
end

function pushbutton_writeZstack_Callback(hObject,~) 
disp('preparing z-stack...');
isMarkAllCells = 1;

hfig = getParentFigure(hObject);
isfullfish = getappdata(hfig,'isfullfish');
if ~isfullfish,
    errordlg('Load full fish first!');
    return;
end

cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
CInfo = getappdata(hfig,'CInfo');
ave_stack = getappdata(hfig,'ave_stack');

timestamp = datestr(now,'mmddyy_HHMMSS');
tiffName = ['stack_' timestamp '.tif'];

U = unique(gIX);
numK = length(U);

% left half: stack with cells marked; 
% right half: original anatomy, or mark all cells

ave_stack2=zeros(size(ave_stack,1), size(ave_stack,2)*2, size(ave_stack,3) ,3);
nPlanes=size(ave_stack,3);
dimv_yxz=size(ave_stack);
stacklen=numel(ave_stack);

circle=makeDisk2(10,21);
[r, v]=find(circle);
r=r-11;v=v-11;
circle_inds  = r*dimv_yxz(1)+v;
cmap = hsv(numK);
weight = 0.3;

for i=1:nPlanes,
    ave_stack2(:,:,i,:)=repmat(imNormalize99(ave_stack(:,:,i)),[1 2 1 3]);
end

for j=1:length(cIX)
    cinds=(CInfo(cIX(j)).center(2)-1)*dimv_yxz(1)+CInfo(cIX(j)).center(1);
    labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yxz(1)*dimv_yxz(2));
    zinds=dimv_yxz(1)*dimv_yxz(2)*2*(CInfo(cIX(j)).slice-1);
    ix = find(U==gIX(j));
    ixs = cinds+circle_inds(labelinds)+zinds;
    ave_stack2(ixs)=cmap(ix,1)*weight + ave_stack2(ixs)*(1-weight);
    ixs = cinds+circle_inds(labelinds)+zinds+stacklen*2;
    ave_stack2(ixs)=cmap(ix,2)*weight + ave_stack2(ixs)*(1-weight);
    ixs = cinds+circle_inds(labelinds)+zinds+stacklen*4;
    ave_stack2(ixs)=cmap(ix,3)*weight + ave_stack2(ixs)*(1-weight);
end
if isMarkAllCells,
    clr = [0,1,0];
    numcell = getappdata(hfig,'numcell');
    for j=1:numcell,
        shift = dimv_yxz(1)*dimv_yxz(2);
        cinds=(CInfo(j).center(2)-1)*dimv_yxz(1)+CInfo(j).center(1);
        labelinds=find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yxz(1)*dimv_yxz(2));
        zinds=dimv_yxz(1)*dimv_yxz(2)*2*(CInfo(j).slice-1);
        ixs = cinds+circle_inds(labelinds)+zinds + shift;
        ave_stack2(ixs)=clr(1)*weight + ave_stack2(ixs)*(1-weight);
        ixs = cinds+circle_inds(labelinds)+zinds+stacklen*2 + shift;
        ave_stack2(ixs)=clr(2)*weight + ave_stack2(ixs)*(1-weight);
        ixs = cinds+circle_inds(labelinds)+zinds+stacklen*4 + shift;
        ave_stack2(ixs)=clr(3)*weight + ave_stack2(ixs)*(1-weight);
    end
end
h = figure;
for i_plane = 1:nPlanes,
    im = squeeze(ave_stack2(:,:,i_plane,:));
    image(im);
    % save tiff
    if (i_plane == 1)
        imwrite(im, tiffName, 'compression','none','writemode','overwrite')
    else
        imwrite(im, tiffName, 'compression','none','writemode','append')
    end
    pause(0.2)
end
close(h)
end

function pushbutton_tileZstack_Callback(hObject,~)
hfig = getParentFigure(hObject);
isfullfish = getappdata(hfig,'isfullfish');
if isfullfish,
    disp('Rendering...')
    DrawTiledPics(hfig);
else
    errordlg('Load full fish first!');
end
end

function out = makeDisk2(radius, dim)
center=floor(dim/2)+1;
out=zeros(dim);
for x=1:dim
    for y=1:dim
        if norm([x,y]-[center,center])<=radius
            out(x,y)=1;
        end
    end
end
end

function pushbutton_popupplot_Callback(hObject,~)
hfig = getParentFigure(hObject);

% same as function RefreshFigure(hfig), except new axes not global

% figure('Position',[50,100,853,512]); % right plot size: 2048x1706
% h1 = axes('Position',[0, 0, 0.5, 1]); % left ~subplot
% h2 = axes('Position',[0.50, 0, 0.5, 1]); % right ~subplot
figure('Position',[50,100,1600,800],'color',[1 1 1]);
h1 = axes('Position',[0.05, 0.03, 0.53, 0.94]); % left ~subplot
h2 = axes('Position',[0.61, 0.03, 0.35, 0.94]); % right ~subplot

CInfo = getappdata(hfig,'CInfo');
M = getappdata(hfig,'M');
dataFR = getappdata(hfig,'dataFR');
fictive = getappdata(hfig,'fictive');
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');

numK = getappdata(hfig,'numK');
stim = getappdata(hfig,'stim');
anat_yx = getappdata(hfig,'anat_yx');
anat_yz = getappdata(hfig,'anat_yz');
anat_zx = getappdata(hfig,'anat_zx');
isCentroid = getappdata(hfig,'isCentroid');
clrmap = getappdata(hfig,'clrmap');
rankscore = getappdata(hfig,'rankscore');
rankID = getappdata(hfig,'rankID');
iswrite = (rankID>=2);
isPlotLines = getappdata(hfig,'isPlotLines');
isPlotFictive = getappdata(hfig,'isPlotFictive');

isPopout = 1;

% left subplot
axes(h1);
if isCentroid,
    [C,~] = FindCentroid(gIX,M);
    DrawClusters(h1,C,unique(gIX),dataFR,numK,stim,fictive,clrmap,rankscore,...
        iswrite,isPopout,isPlotLines,isPlotFictive);
%     DrawClusters(h1,C,unique(gIX),dataFR,numK,stim,fictive,clrmap,rankscore,iswrite,ispopout);
else
    DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,...
        iswrite,isPopout,isPlotLines,isPlotFictive);
%     DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,iswrite,ispopout);
end

% right subplot
axes(h2);
DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,clrmap,'full');

end

function checkbox_isPlotLines_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotLines',get(hObject,'Value'));
end

function checkbox_isPlotFictive_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotFictive',get(hObject,'Value'));
end

function pushbutton_exporttoworkspace_Callback(hObject,~)
hfig = getParentFigure(hObject);
assignin('base', 'M', getappdata(hfig,'M'));
assignin('base', 'cIX', getappdata(hfig,'cIX'));
assignin('base', 'gIX', getappdata(hfig,'gIX'));
assignin('base', 'numK', getappdata(hfig,'numK'));
assignin('base', 'stim', getappdata(hfig,'stim'));
assignin('base', 'fictive', getappdata(hfig,'fictive'));
assignin('base', 'M_0_fluo', getappdata(hfig,'M_0_fluo'));
end

%% row 2: Load

function popup_fishmenu_Callback(hObject,~)
new_i_fish = get(hObject,'Value');
hfig = getParentFigure(hObject);
QuickUpdateFish(hfig,new_i_fish);
end

function QuickUpdateFish(hfig,new_i_fish,init) %#ok<INUSD>
M_fish_set = getappdata(hfig,'M_fish_set');
fishset = M_fish_set(new_i_fish);
setappdata(hfig,'fishset',fishset);

CONSTs = getappdata(hfig,'CONSTs');
setappdata(hfig,'isfullfish',0);

% load all fields from CONSTs, with names preserved
names = fieldnames(CONSTs{new_i_fish}); % cell of strings
for i = 1:length(names),
    setappdata(hfig,names{i},eval(['CONSTs{new_i_fish}.',names{i}]));
end

% recontruct the 3 big matrices from CIX to original size (numcell)
CIX = CONSTs{new_i_fish}.CIX;
numcell = CONSTs{new_i_fish}.numcell;

temp = CONSTs{new_i_fish}.CRAZ;
CellRespAvr = zeros(numcell,size(temp,2));
CellRespAvr(CIX,:) = temp;

temp = CONSTs{new_i_fish}.CRZt;
CellResp = zeros(numcell,size(temp,2));
CellResp(CIX,:) = temp;

temp = CONSTs{new_i_fish}.CInfo;
CInfo(numcell).center = '';
CInfo(numcell).area = '';
CInfo(numcell).slice = '';
CInfo(numcell).inds = '';
CInfo(numcell).y_minmax = '';
CInfo(numcell).x_minmax = '';
CInfo(CIX) = temp;

setappdata(hfig,'CellRespAvr',CellRespAvr); 
setappdata(hfig,'CellResp',CellResp); 
setappdata(hfig,'CInfo',CInfo); 

UpdateFishData(hfig,new_i_fish);
if ~exist('init','var'),
    UpdateFishDisplay(hfig);
end
end

function UpdateFishData(hfig,new_i_fish) % loading steps shared by both quick-load and full-load
UpdateDatamenu_Direct(hfig,1);

% save ClusGroup before updating, if applicable
clusgroupID = getappdata(hfig,'clusgroupID');
if ~isempty(clusgroupID),
    new_clusgroupID = 1;
    UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID,'norefresh'); % to save ClusGroup
end
% set new i_fish after UpdateClusGroupID, which saves into old 'i_fish'
setappdata(hfig,'i_fish',new_i_fish);

% load from VAR
global VAR;
Class = VAR(new_i_fish).Class; % actually collection of classes, each class can have many groups
setappdata(hfig,'Class',Class);

ClusGroup = VAR(new_i_fish).ClusGroup; % group of 'Cluster', i.e. 2nd level collections
setappdata(hfig,'ClusGroup',ClusGroup);

clusgroupID = 1;
setappdata(hfig,'clusgroupID',clusgroupID);

Cluster = ClusGroup{clusgroupID};% collection of clusters, 'Cluster' structurally same as 'Class'
setappdata(hfig,'Cluster',Cluster);
end

function UpdateFishDisplay(hfig) % loading steps that are not performed at very first initialization
stim = getappdata(hfig,'stim');
fishset = getappdata(hfig,'fishset');

[~, names] = GetStimRegressor(stim,fishset);

global hstimreg;
set(hstimreg,'String',['(choose)',(names)]);

% update display
UpdateClassID(hfig,1,'norefresh');
clusgroupID = 1;
new_clusgroupID = 1;
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID); % to display new menu
end

function popup_loadfullfishmenu_Callback(hObject,~)
new_i_fish = get(hObject,'Value')-1;
if new_i_fish>0,
    hfig = getParentFigure(hObject);
    LoadFullFish(hfig,new_i_fish);
end
global hfishnum;
set(hfishnum,'Value',new_i_fish);
end

function LoadFullFish(hfig,new_i_fish)
data_dir=getappdata(hfig,'data_dir');
disp(['loading fish #' num2str(new_i_fish) '...']);

M_fish_set = getappdata(hfig,'M_fish_set');
fishset = M_fish_set(new_i_fish);
setappdata(hfig,'fishset',fishset);

%%
if fishset == 1 || fishset == 2,
    tic
    load(fullfile(data_dir,['CONST_F' num2str(new_i_fish) '.mat']),'CONST');
    toc
    setappdata(hfig,'CellResp',CONST.CRZt);
    const = CONST;
    
else % load from parts
    tic
    fishdir = fullfile(data_dir,['CONST_F' num2str(new_i_fish) '_fast.mat']);
    load(fishdir,'const');
    load(fishdir,'dimCR');
    CellResp = zeros(dimCR);
    num = 0;
    nParts = round(dimCR(1)*dimCR(2)/(10^8));
    disp(['in ' num2str(nParts) ' parts']);
    for i = 1:nParts,
        disp(num2str(i));
        load(fishdir,['CellResp_' num2str(i)]);
        eval(['len = size(CellResp_' num2str(i) ',1);']);
        eval(['CellResp(num+1:num+len,:) = CellResp_' num2str(i) ';']);
        eval(['CellResp(num+1:num+len,:) = CellResp_' num2str(i) ';']);
        num = num+len;
    end
    % for i = 1:nParts,
    %     load(fishdir,['CellResp_' num2str(i)']);
    %     eval(['len = size(CellResp_' num2str(i) ',1);']);
    %     eval(['CellResp(num+1:num+len,:) = CellResp_' num2str(i) ';']);
    %     num = num+len;
    % end
    toc
    setappdata(hfig,'CellResp',CellResp);
end

    
%% load all fields from CONST, with names preserved
% names = fieldnames(const); % cell of strings
names = fieldnames(const); % cell of strings
for i = 1:length(names),
    
    
    
    % renaming exception!!!!!!!!!!!!should be obsolete!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if strcmp(names{i},'CRAZ'),
        setappdata(hfig,'CellRespAvr',const.CRAZ);
    elseif strcmp(names{i},'photostate'),
        setappdata(hfig,'stim_full',const.photostate);
        
        
        
    else
        setappdata(hfig,names{i},eval(['const.',names{i}]));
    end
end
setappdata(hfig,'numcell',length(const.CInfo));

UpdateFishData(hfig,new_i_fish);
UpdateFishDisplay(hfig);
setappdata(hfig,'isfullfish',1);
end

function popup_datamenu_Callback(hObject,~)
ID = get(hObject,'Value');
hfig = getParentFigure(hObject);
UpdateDatamenu_Direct(hfig,ID);
SetDataFR(hfig,1); % also set M (without updating cIX gIX)
RefreshFigure(hfig);
end

function UpdateDatamenu_Direct(hfig,ID)
CellResp = getappdata(hfig,'CellResp');
Fc = getappdata(hfig,'Fc');
datanames = getappdata(hfig,'datanames');
tlists = getappdata(hfig,'tlists');
stim_full = getappdata(hfig,'stim_full');
periods = getappdata(hfig,'periods');
shift = getappdata(hfig,'shift');
numcell = getappdata(hfig,'numcell');

fishset = getappdata(hfig,'fishset');

FcAvr = getappdata(hfig,'FcAvr');

if fishset == 1,
    CellRespAvr = getappdata(hfig,'CellRespAvr');    
else
    % find averages of different stimuli, and splice together
    % in future should allow user to chose stimtypelist to compile average
    % of selected stimuli groups
    CellRespAvr = [];   
    stimAvr = [];
    for i = 1:length(tlists)-1,
        M = circshift(CellResp(:,tlists{i}),shift,2);
        CellRespAvr = horzcat(CellRespAvr,mean(reshape(M,numcell,periods(i),[]),3));
        % stim from stim_full
        m = stim_full(tlists{i});
        stimAvr = horzcat(stimAvr,m(1:periods(i)));
    end
    setappdata(hfig,'CellRespAvr',CellRespAvr);
end


%% SWITCH LEFT/RIGHT DIRECTION AGAIN
FcAvr = vertcat(FcAvr(2,:),FcAvr(1,:),FcAvr(3:end,:));
Fc = vertcat(Fc(2,:),Fc(1,:),Fc(3:end,:));
% now 1=left, 2=right, 3=forward, 4 = raw left, 5 = raw right

%%
if ID == 1,
    M_0 = CellRespAvr;
    fictive = FcAvr;
elseif ID == 2,
    M_0 = CellResp;
    fictive = Fc;
else    
    if fishset == 1,
        IX = tlists{ID};
    else
        IX = tlists{ID-2};
    end
    M_0 = CellResp(:,IX);
    temp = Fc(:,IX);
    % normalize fictive channels, in 2 sets
    for k = 1:3,%size(F,1),
        m = temp(k,:);
        temp(k,:) = (m-min(m))/(max(m)-min(m));
    end
    m = temp(4:5,:);
    temp(4:5,:) = (m-min(min(m)))/(max(max(m))-min(min(m)));
    fictive = temp;
end

% set stim
if fishset == 1,
    stim = stim_full(tlists{ID});
else
    if ID == 1,
        stim = stimAvr;
    elseif ID == 2,
        stim = stim_full;
    else
        IX = tlists{ID-2};
        stim = stim_full(IX);
    end
end

setappdata(hfig,'M_0_fluo',M_0);
setappdata(hfig,'fictive',fictive);
setappdata(hfig,'fictive_0',fictive);
setappdata(hfig,'dataname',datanames(ID));
setappdata(hfig,'stim',stim);

global hdatamenu;
if ~isempty(hdatamenu), % before GUI initialization
    % 'global' remnant from previous run: hdatamenu is not empty but is invalid object
    if isvalid(hdatamenu), 
        set(hdatamenu,'String',datanames);
        set(hdatamenu,'Value',ID);
    end
end
end

function checkbox_showcentroids_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isCentroid',get(hObject,'Value'));
RefreshFigure(hfig);
end

%% ----- tab two ----- (Operations)

%% row 1: Range

function pushbutton_back_Callback(hObject,~)
global hback hfwd;
hfig = getParentFigure(hObject);
bC = getappdata(hfig,'bCache');
fC = getappdata(hfig,'fCache');

if ~isempty(bC.cIX{1}),
    % set last into forward cache
    fC.cIX = [getappdata(hfig,'cIX'),fC.cIX];
    fC.gIX = [getappdata(hfig,'gIX'),fC.gIX];
    fC.numK = [getappdata(hfig,'numK'),fC.numK];
    % retrieve
    cIX = bC.cIX{1};
    gIX = bC.gIX{1};
    numK = bC.numK{1};
    bC.cIX(1) = [];
    bC.gIX(1) = [];
    bC.numK(1) = [];
    % set M
    M_0 = GetM_0(hfig);
    M = M_0(cIX,:);
    % save
    setappdata(hfig,'M',M);
    setappdata(hfig,'bCache',bC);
    setappdata(hfig,'fCache',fC);
    setappdata(hfig,'cIX',cIX);
    setappdata(hfig,'gIX',gIX);
    setappdata(hfig,'numK',numK);
    % handle rankID: >=2 means write numbers as text next to colorbar
    setappdata(hfig,'rankID',0);
    % finish
    disp('back (from cache)')
    RefreshFigure(hfig);
    set(hfwd,'enable','on');
else % nothing to retrieve
    set(hback,'enable','off');
end
end

function pushbutton_forward_Callback(hObject,~)
global hback hfwd;
hfig = getParentFigure(hObject);
bC = getappdata(hfig,'bCache');
fC = getappdata(hfig,'fCache');

if ~isempty(fC.cIX{1}),
    % set last into (backward) cache
    bC.cIX = [getappdata(hfig,'cIX'),bC.cIX];
    bC.gIX = [getappdata(hfig,'gIX'),bC.gIX];
    bC.numK = [getappdata(hfig,'numK'),bC.numK];
    % retrieve
    cIX = fC.cIX{1};
    gIX = fC.gIX{1};
    numK = fC.numK{1};
    fC.cIX(1) = [];
    fC.gIX(1) = [];
    fC.numK(1) = [];
    % set M
    M_0 = GetM_0(hfig);
    M = M_0(cIX,:);
    % save
    setappdata(hfig,'M',M);
    setappdata(hfig,'bCache',bC);
    setappdata(hfig,'fCache',fC);
    setappdata(hfig,'cIX',cIX);
    setappdata(hfig,'gIX',gIX);
    setappdata(hfig,'numK',numK);
    % handle rankID: >=2 means write numbers as text next to colorbar
    setappdata(hfig,'rankID',0);
    % finish
    disp('forward (from cache)')
    RefreshFigure(hfig);
    set(hback,'enable','on');
else
    set(hfwd,'enable','off');
end
end

function edit_choose_range_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
% get/format range
str = get(hObject,'String'); 
if ~isempty(str),
    str = strrep(str,'end',num2str(max(gIX)));
    range = ParseRange(str);
    % update indices
    tempI = [];
    for i = range,
        tempI = [tempI;find(gIX==i)];
    end
    cIX = cIX(tempI);
    gIX = gIX(tempI);
    UpdateIndices(hfig,cIX,gIX);
    RefreshFigure(hfig);
end
end

function range = ParseRange(str)
str = strrep(str,':','-'); % e.g. str= '1,3,5:8';
C = textscan(str,'%d','delimiter',',');
m = C{:};
range = [];
for i = 1:length(m),
    if m(i)>0,
        range = [range,m(i)]; %#ok<AGROW>
    else % have '-'sign,
        range = [range,m(i-1)+1:-m(i)]; %#ok<AGROW>
    end
end
end

function edit_exclude_range_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
% get/format range
str = get(hObject,'String');
if ~isempty(str),
    str = strrep(str,'end',num2str(max(gIX)));
    range = ParseRange(str);
    range = setdiff(unique(gIX),range)';
    % update indices
    tempI = [];
    for i = range,
        tempI = [tempI;find(gIX==i)];
    end
    cIX = cIX(tempI);
    gIX = gIX(tempI);
    UpdateIndices(hfig,cIX,gIX);
    RefreshFigure(hfig);
end
end

function edit_fuse_range_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
% get/format range
str = get(hObject,'String');
if ~isempty(str),
    str = strrep(str,'end',num2str(max(gIX)));
    range = ParseRange(str);
    % update indices
    tempI = [];
    for i = range,
        tempI = [tempI;find(gIX==i)];
    end
    gIX(tempI) = gIX(tempI(1));
    [gIX, numK] = SqueezeGroupIX(gIX);
    UpdateIndices(hfig,cIX,gIX,numK);
    RefreshFigure(hfig);
end
end

function edit_t_range_Callback(hObject,~)
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');

% get/format range
str = get(hObject,'String'); 
if ~isempty(str),
    str = strrep(str,'end',num2str(size(M,2)));
    range = ParseRange(str);
    
    UpdateTRange(hfig,range)
    RefreshFigure(hfig);
end
end

function UpdateTRange(hfig,range) % this function is incomplete... delete?
% what about stimulus???

M_0 = GetM_0(hfig);
cIX = getappdata(hfig,'cIX');
fictive_0 = getappdata(hfig,'fictive_0');

% set M
M = M_0(cIX,range);
setappdata(hfig,'M',M);
% set fictive
fictive = fictive_0(:,range);
setappdata(hfig,'fictive',fictive);
end

%% row 2: Operations

function popup_operations_Callback(hObject,~)
opID = get(hObject,'Value') - 1;
hfig = getParentFigure(hObject);
setappdata(hfig,'opID',opID);
% highlight UI
global hopID;
set(hopID,'BackgroundColor',[1,1,0.8]);
end

function popup_ranking_Callback(hObject,~)
% menu = {'(ranking)','hier.','size','stim-lock','corr','noise'};
rankID = get(hObject,'Value') - 1;
hfig = getParentFigure(hObject);
setappdata(hfig,'rankID',rankID);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');

[gIX, numU] = SqueezeGroupIX(gIX);
switch rankID,
    case 1,
        disp('hier. (default)');
        [gIX, numU] = HierClus(M,gIX);
    case 2,
        disp('size');
        H = zeros(numU,1);
        for i = 1:numU,
            H(i) = length(find(gIX==i));
        end
        [gIX,rankscore] = SortH(H,gIX,numU,'descend');
    case 3,
        disp('stim-lock');
        [gIX,rankscore] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU);
    case 4,
        disp('corr');
        [~,D] = FindCentroid(gIX,M);
        [gIX,rankscore] = SortH(D,gIX,numU);
        rankscore = 1-rankscore;
    case 5,
        disp('motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,M,numU,1);
    case 6,
        disp('motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,M,numU,2);
    case 7,
        disp('motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,M,numU,3);
    case 8,
        disp('motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,M,numU,4);
    case 9,
        disp('stim-motor');
        [~,rankscore1] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU);
        [~,rankscore2] = RankByMotorStim_Direct(hfig,gIX,M,numU,1);
        rankscore = rankscore1 - rankscore2;
end
if rankID>1,
    setappdata(hfig,'rankscore',round(rankscore*100)/100);
    setappdata(hfig,'clrmap','jet');
else
    setappdata(hfig,'clrmap','hsv');
end
UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
disp('ranking complete');
end

function [gIX,rankscore] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU)
periods = getappdata(hfig,'periods');

if length(periods)==1,
    period = periods{1};
    [C,~] = FindCentroid(gIX,M);
else %if i_fish==8,
    tlists = getappdata(hfig,'tlists');
    CellResp = getappdata(hfig,'CellResp');
    IX = tlists{6}; % ptomr_circ
    M_ = CellResp(cIX,IX);
    [C,~] = FindCentroid(gIX,M_);
    periods = getappdata(hfig,'periods');
    period = periods{1}+periods{2};
% else ?

end
C_3D_0 = reshape(C,size(C,1),period,[]);
C_3D = zscore(C_3D_0,0,2);

H = nanmean(nanstd(C_3D,0,3),2);
[gIX,rankscore] = SortH(H,gIX,numU);
end

function [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,M,numU,option)
C = FindCentroid(gIX,M);
fictive = getappdata(hfig,'fictive');

reg = zeros(3,length(fictive));
reg(1,:) = fictive(5,:); % Left
reg(2,:) = fictive(4,:); % Right
reg(3,:) = mean(fictive(4:5,:));
H = zeros(numU,1);
shift = zeros(numU,1);
a = zeros(1,3);
I = zeros(1,3);
for i = 1:numU,
    switch option,
        case 1,
            for j = 1:3,
                [a(j),I(j)] = max(abs(xcorr(C(i,:),reg(j,:),'coeff')));
            end
            [H(i),jmax] = max(a);
            shift(i) = I(jmax) - length(fictive);
        case 2,
            [H(i),I] = max(xcorr(C(i,:),reg(1,:),'coeff'));
            shift(i) = I - length(fictive);
        case 3,
            [H(i),I] = max(xcorr(C(i,:),reg(2,:),'coeff'));
            shift(i) = I - length(fictive);
        case 4,
            [H(i),I] = max(xcorr(C(i,:),reg(3,:),'coeff'));
            shift(i) = I - length(fictive);
    end
end
[gIX,rankscore] = SortH(H,gIX,numU,'descend');
assignin('base', 'shift', shift);
end

function [gIX,B] = SortH(H,gIX,numU,descend) % new gIX is sorted based on H, size(H)=[numU,1];
if exist('descend','var'),
    [B,I] = sort(H,'descend');
else
    [B,I] = sort(H);
end
gIX_last = gIX;
for i = 1:numU,
    gIX(gIX_last==I(i)) = i;
end
end

function pushbutton_sqeeze_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
[gIX, numU] = SqueezeGroupIX(gIX);
UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
end

function pushbutton_flipud_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');

U = unique(gIX);
U_ = flipud(U);

gIX_ = gIX;
for i = 1:length(U),
    gIX_(gIX==U(i)) = U_(i);
end

UpdateIndices(hfig,cIX,gIX_);
RefreshFigure(hfig);
end

function pushbutton_clrmap_Callback(hObject,~)
hfig = getParentFigure(hObject);
clrmap = getappdata(hfig,'clrmap');
if strcmp(clrmap,'jet'),
    setappdata(hfig,'clrmap','hsv');
else
    setappdata(hfig,'clrmap','jet');
end
RefreshFigure(hfig);
end

%% row 3: Anatomy

function pushbutton_polygon_yx_Callback(hObject,~)
k_zres = 20;

hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
% M = getappdata(hfig,'M');
CInfo = getappdata(hfig,'CInfo');
numcell = getappdata(hfig,'numcell');
anat_yx = getappdata(hfig,'anat_yx');
anat_zx = getappdata(hfig,'anat_zx');
dimv_yx = size(anat_yx);
dimv_zx = size(anat_zx);
isfullfish = getappdata(hfig,'isfullfish');

h_poly_yx = impoly;
wait(h_poly_yx); % double click to finalize position!
% update finalized polygon in bright color
setColor(h_poly_yx,[0 1 1]);

if isfullfish,
    IJs = reshape([CInfo.center],2,[])';
else % Matlab can't handle CInfo with all the empty entries -> manual padding
    IJs_ = reshape([CInfo(cIX).center],2,[])';
    IJs = ones(numcell,2);
    IJs(cIX,:) = IJs_;
end
A = sub2ind(dimv_yx(1:2),IJs(:,1),IJs(:,2));
MaskArray = createMask(h_poly_yx);
MaskArray(1:dimv_zx*k_zres+10,:) = [];
B = find(MaskArray); % find indices of pixels within ROI

temp = find(ismember(A,B));
[IX,ia,~] = intersect(cIX,temp);

cIX = IX;
gIX = gIX(ia);

UpdateIndices(hfig,cIX,gIX);
RefreshFigure(hfig);
end

function pushbutton_polygon_yz_Callback(hObject,~)
k_zres = 20;

hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
% M = getappdata(hfig,'M');
CInfo = getappdata(hfig,'CInfo');
numcell = getappdata(hfig,'numcell');
anat_yx = getappdata(hfig,'anat_yx');
anat_yz = getappdata(hfig,'anat_yz');
anat_zx = getappdata(hfig,'anat_zx');
dimv_yx = size(anat_yx);
dimv_yz = size(anat_yz);
dimv_zx = size(anat_zx);
isfullfish = getappdata(hfig,'isfullfish');

h_poly_z = impoly;
wait(h_poly_z); % double click to finalize position!
% update finalized polygon in bright color
setColor(h_poly_z,[0 1 1]);

if isfullfish,
        IJs = reshape([CInfo.center],2,[])';
        Zs = [CInfo.slice]';
else % Matlab can't handle CInfo with all the empty entries -> manual padding
    IJs_ = reshape([CInfo(cIX).center],2,[])';
    IJs = ones(numcell,2);
    IJs(cIX,:) = IJs_;
    Zs_ = [CInfo(cIX).slice]';
    Zs = ones(numcell,1);
    Zs(cIX) = Zs_;
end
A = sub2ind(dimv_yz(1:2),IJs(:,1),Zs);

MaskArray = createMask(h_poly_z);
MaskArray(1:dimv_zx*k_zres+10,:) = [];
MaskArray(:,1:dimv_yx(2)+10) = [];
zMaskArray = imresize(MaskArray,[dimv_yz(1),dimv_yz(2)]);
B = find(zMaskArray); % find indices of pixels within ROI

temp = find(ismember(A,B));
[IX,ia,~] = intersect(cIX,temp);

cIX = IX;
gIX = gIX(ia);

UpdateIndices(hfig,cIX,gIX);
RefreshFigure(hfig);
end

function pushbutton_polygon_zx_Callback(hObject,~)
k_zres = 20;

hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
% M = getappdata(hfig,'M');
CInfo = getappdata(hfig,'CInfo');
numcell = getappdata(hfig,'numcell');
% anat_yx = getappdata(hfig,'anat_yx');
anat_zx = getappdata(hfig,'anat_zx');
% dimv_yx = size(anat_yx);
dimv_zx = size(anat_zx);
isfullfish = getappdata(hfig,'isfullfish');

h_poly_z = impoly;
wait(h_poly_z); % double click to finalize position!
% update finalized polygon in bright color
setColor(h_poly_z,[0 1 1]);

if isfullfish,
        IJs = reshape([CInfo.center],2,[])';
        Zs = [CInfo.slice]';
else % Matlab can't handle CInfo with all the empty entries -> manual padding
    IJs_ = reshape([CInfo(cIX).center],2,[])';
    IJs = ones(numcell,2);
    IJs(cIX,:) = IJs_;
    Zs_ = [CInfo(cIX).slice]';
    Zs = ones(numcell,1);
    Zs(cIX) = Zs_;
end
A = sub2ind(dimv_zx(1:2),Zs,IJs(:,2));

MaskArray = createMask(h_poly_z);

MaskArray(dimv_zx(1)*k_zres+1:end,:) = [];
MaskArray(:,dimv_zx(2)+1:end) = [];
zMaskArray = imresize(MaskArray,[dimv_zx(1),dimv_zx(2)]);
zMaskArray = flipud(zMaskArray); % careful!
B = find(zMaskArray); % find indices of pixels within ROI

temp = find(ismember(A,B));
[IX,ia,~] = intersect(cIX,temp);

cIX = IX;
gIX = gIX(ia);

UpdateIndices(hfig,cIX,gIX);
RefreshFigure(hfig);
end

function pushbutton_withinConvexHull_Callback(hObject,~) % only works with full load
disp('find points within boundary...');
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
CInfo = getappdata(hfig,'CInfo');
IJs = reshape([CInfo.center],2,[])';
Zs = [CInfo.slice]';
% point sets
xyz_all = horzcat(IJs, Zs);
xyz_cIX = horzcat(IJs(cIX,:), Zs(cIX));

tri = delaunayn(xyz_cIX); % Generate delaunay triangulization
t = tsearchn(xyz_cIX, tri, xyz_all); % Determine which triangle point is within
I_inside = ~isnan(t);

cIX_1 = cIX;
gIX_1 = gIX;
cIX_2 = setdiff(find(I_inside),cIX_1);
gIX_2 = ones(length(cIX_2),1)*(max(gIX_1)+1);
cIX = [cIX_1;cIX_2];
gIX = [gIX_1;gIX_2];

UpdateIndices(hfig,cIX,gIX,length(unique(gIX)));
RefreshFigure(hfig);
disp('search complete');
end

%% ----- tab three ----- (Regression)

%% row 1: regressor

function popup_getstimreg_Callback(hObject,~)
i_reg = get(hObject,'Value')-1;
hfig = getParentFigure(hObject);
if i_reg==0,
    return;
end
setappdata(hfig,'regchoice',{1,i_reg});
% highlight the choice (yellow)
global hstimreg hmotorreg hcentroidreg;
set(hstimreg,'BackgroundColor',[1,1,0.8]); % yellow
set(hmotorreg,'BackgroundColor',[1,1,1]);
set(hcentroidreg,'BackgroundColor',[1,1,1]);

isPlotReg = getappdata(hfig,'isPlotReg');
if isPlotReg,
    PlotRegWithStimMotor(hfig);
end
end

function PlotRegWithStimMotor(hfig)
stim = getappdata(hfig,'stim');
fictive = getappdata(hfig,'fictive');
regressor = GetRegressor(hfig);
% plot regressor
figure;
halfbarheight = 1;
stimbar = GetStimBar(halfbarheight,stim);
reg = imNormalize99(regressor);
regbar = repmat(reg,[1,1,3]);

h = subplot(3,1,1);image(stimbar);set(h,'box','on','xtick',[],'ytick',[]); title('stimulus');
h = subplot(3,1,2);image(regbar); set(h,'box','on','xtick',[],'ytick',[]); title('regressor');
subplot(3,1,3);

temp = fictive;
if 1, % plot all 5 lines
    fc = vertcat(temp(1,:),temp(3,:),temp(2,:),temp(4,:),temp(5,:));
    
    imagesc(fc);colormap gray
    set(gca,'YTick',[],'XTick',[]);
    set(gcf,'color',[1 1 1]);
    set(gca, 'box', 'off')
    hold on;axis ij;
    
    % plot division lines
%     for i = 0:3,
%         y = i+0.5;
%         plot([0.5,length(fictive)+0.5],[y,y],'w','Linewidth',0.5);
%     end
%     % labels
%     names = {'Left','Right','Forward','Raw L','Raw R'};
%     x = -s2*0.05;
%     for i = 1:5,
%         y = i;
%         text(x,y,names{i},'Fontsize',7);
%     end
    
else % only plot top 3 lines
    fc = vertcat(temp(1,:),temp(3,:),temp(2,:));
    imagesc(fc);colormap gray
    set(gca,'YTick',[],'XTick',[]);
    set(gcf,'color',[1 1 1]);
    set(gca, 'box', 'off')
    hold on;axis ij;
    
    % plot division lines
    for i = 0:2,
        y = i+0.5;
        plot([0.5,length(fictive)+0.5],[y,y],'w','Linewidth',0.5);
    end
    % labels
    names = {'Left','Forward','Right'};
    x = -s2*0.07;
    for i = 1:3,
        y = i;
        text(x,y,names{i},'Fontsize',10);
    end
    
end

imagesc(fictive);colormap hot; axis off; title('motor');
end

function out=imNormalize99(im)
im=double(im);
temp=sort(im(:),'descend');
th1=temp(round(length(im(:))/100));
th2=min(im(:));

out=(im-th2)/(th1-th2);
out(out>1)=1;
end

function popup_getmotorreg_Callback(hObject,~)
hfig = getParentFigure(hObject);
i_reg = get(hObject,'Value')-1;
% {'left swims','right swims','forward swims','raw left','raw right','raw average'};
if i_reg>0,
    setappdata(hfig,'regchoice',{2,i_reg});
end
% highlight the choice
global hstimreg hmotorreg hcentroidreg;
set(hstimreg,'BackgroundColor',[1,1,1]);
set(hmotorreg,'BackgroundColor',[1,1,0.8]); % yellow
set(hcentroidreg,'BackgroundColor',[1,1,1]);

isPlotReg = getappdata(hfig,'isPlotReg');
if isPlotReg,
    PlotRegWithStimMotor(hfig);
end
end

function checkbox_popupplotreg_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotReg',get(hObject,'Value'));
end

function edit_ctrdID_as_reg_Callback(hObject,~)
str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d');
    ctrdID = temp{:};
end
hfig = getParentFigure(hObject);
setappdata(hfig,'regchoice',{3,ctrdID});
% highlight the choice
global hstimreg hmotorreg hcentroidreg;
set(hstimreg,'BackgroundColor',[1,1,1]);
set(hmotorreg,'BackgroundColor',[1,1,1]);
set(hcentroidreg,'BackgroundColor',[1,1,0.8]); % yellow
end

function checkbox_flipstim_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isflipstim',get(hObject,'Value'));
end

%% row 2: regression

function regressor = GetRegressor(hObject)
hfig = getParentFigure(hObject);
regchoice = getappdata(hfig,'regchoice');
stim = getappdata(hfig,'stim');
    
if regchoice{1}==1, % stim Regressor
    fishset = getappdata(hfig,'fishset');
    [regressors, names] = GetStimRegressor(stim,fishset);    
    regressor = regressors(regchoice{2}).im;
    
elseif regchoice{1}==2, % motor Regressor
    fictive = getappdata(hfig,'fictive');
    regressors = GetMotorRegressor(fictive);
    regressor = regressors(regchoice{2}).im;
    
else % regchoice{1}==3, from Centroid
    ctrdID = regchoice{2};
    isflipstim = getappdata(hfig,'isflipstim');
    M = getappdata(hfig,'M');
    gIX = getappdata(hfig,'gIX');
    
    i = find(unique(gIX)==ctrdID);
    if isempty(i),
        disp('input is empty!');beep;
        regressor = [];
        return;
    end
    [C,~] = FindCentroid(gIX,M);
    regressor = C(i,:);
    
    if isflipstim, % swap left/right stim
        % 0 = all black; 1 = black/white; 2 = white/black; 3 = all white; 4 = all gray;
        % 10 = forward grating (very slow, more for calibration)
        % 11 = rightward grating
        % 12 = leftward grating
        % swap: 1~2,11~12
        
        %     if i_fish == 6 || i_fish == 7,
        %         halflength = length(stim)/2;
        %         stim = stim(1:halflength); % length is always even
        %     end
        flipstim = stim;
        flipstim(stim==1)=2;
        flipstim(stim==2)=1;
        flipstim(stim==11)=12;
        flipstim(stim==12)=11;
        
        phototrans = GetPhotoTrans(stim);
        fliptrans = GetPhotoTrans(flipstim);
        
        % transition e.g.:
        % stimulus    : 2     3     1     3     3     0     0     1     1     0     3     2     0     2     2     1
        % transitionID: 6    11    13     7    15    12     0     1     5     4     3    14     8     2    10     9
        
        IX = [];
        targetstates = unique(fliptrans,'stable');
        for i = 1:length(targetstates),
            IX = [IX, find(phototrans==targetstates(i))];
        end
        %     if i_fish == 6 || i_fish == 7,
        %         reg = regressor(1:halflength);
        %         reg = reg(IX);
        %         regressor = repmat(reg,1,2);
        %     else
        regressor = regressor(IX);
        %     end
    end    
end

end

function checkbox_plotcorrhist_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotCorrHist',get(hObject,'Value'));
end

function pushbutton_topnum_regression_Callback(hObject,~)
disp('regression...');
hfig = getParentFigure(hObject);
M_0 = getappdata(hfig,'M_0_fluo');
numTopCorr = getappdata(hfig,'numTopCorr');
regressor = GetRegressor(hfig);

%% for each cell, find correlation coeff
M_corr = corr(regressor',M_0');%cell_resp_ave');
% M_corr_ctrl = corr(regressor.ctrl',M_0');

[~,I] = sort(M_corr,'descend');
cIX = I(1:numTopCorr)';
gIX = ceil((1:length(cIX))'/length(cIX)*min(20,length(cIX)));

score = M_corr(cIX);
rankscore = zeros(1,20);
for i = 1:20,
    ix = find(gIX==i);
    rankscore(i) = mean(score(ix));
end
setappdata(hfig,'rankscore',round(rankscore*100)/100);
setappdata(hfig,'rankID',99);

UpdateIndices(hfig,cIX,gIX,20);
RefreshFigure(hfig);
end

function edit_topCorrNumber_Callback(hObject,~)
str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d');
    numTopCorr = temp{:};
    hfig = getParentFigure(hObject);
    setappdata(hfig,'numTopCorr',numTopCorr);
end
end

function pushbutton_thres_regression_Callback(hObject,~) 
disp('regression...');
hfig = getParentFigure(hObject);
M_0 = getappdata(hfig,'M_0_fluo');
thres_reg = getappdata(hfig,'thres_reg');
regressor = GetRegressor(hfig);
if isempty(regressor),
    return;
end

isPlotCorrHist = getappdata(hfig,'isPlotCorrHist');

[cIX,gIX] = Regression_direct(M_0,thres_reg,regressor,isPlotCorrHist);

if ~isempty(gIX),
    gIX = ceil((1:length(cIX))'/length(cIX)*min(20,length(cIX)));
    UpdateIndices(hfig,cIX,gIX,40);
    RefreshFigure(hfig);
end
end

function [cIX_,gIX_,R] = Regression_direct(M_0,thres_reg,regressor,isPlotCorrHist) % gIX_ is just ones
R = corr(regressor',M_0');
if thres_reg>0,
    cIX_ = (find(R>thres_reg))';
elseif thres_reg<0,
    cIX_ = (find(R<thres_reg))';
end
if isempty(cIX_),
    disp('result is empty!');beep;
    gIX_ = [];
    return;
end
[~,I] = sort(R(cIX_),'descend');
cIX_ = cIX_(I);
gIX_ = ones(size(cIX_));
%%
if exist('isPlotCorrHist','var'),
    if isPlotCorrHist,
        figure('Position',[500,200,500,150]);
        hold on;
        bins = -1:0.025:1;
        [N,~] = histcounts(R,bins);
        histogram(R,bins,'FaceColor',[0.4 0.4 0.4]);%,'EdgeColor','none'
        plot([thres_reg,thres_reg],[0,max(N)],'r--');
        xlim([-1,1]);ylim([0,max(N)]);
    end
end
end

function edit_regthres_Callback(hObject,~)
str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%f');
    thres_reg = temp{:};
    hfig = getParentFigure(hObject);
    setappdata(hfig,'thres_reg',thres_reg);
end
end

function pushbutton_AllCentroidRegression_Callback(hObject,~)
hfig = getParentFigure(hObject);
M_0 = getappdata(hfig,'M_0_fluo');
thres_reg = getappdata(hfig,'thres_reg');
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');

[cIX,gIX,numK] = AllCentroidRegression_Direct(M_0,thres_reg,cIX,gIX);

UpdateIndices(hfig,cIX,gIX,numK);
RefreshFigure(hfig);

disp('all regression complete');
end

function [cIX,gIX,numK] = AllCentroidRegression_Direct(M_0,thres_reg,cIX,gIX)
U = unique(gIX);
clusters(length(U)).cIX = [];
for i = 1:length(U),
    disp(['i = ' num2str(i)]);
    IX = find(gIX == U(i));    
    [~,regressor] = kmeans(M_0(cIX(IX),:),1,'distance','correlation'); %#ok<FNDSB>
    cIX_ = Regression_direct(M_0,thres_reg,regressor);
    clusters(i).cIX = cIX_;
end
disp('union...');
CIX = vertcat(clusters(1:end).cIX);
GIX = [];
for i = 1:numel(clusters),
    GIX = [GIX; ones(length(clusters(i).cIX),1)*i]; %#ok<AGROW>
end
[cIX,gIX,numK] = SmartUnique(CIX,GIX,M_0(CIX,:));
end

function pushbutton_IterCentroidRegression_Callback(hObject,~)
hfig = getParentFigure(hObject);
M_0 = getappdata(hfig,'M_0_fluo');
thres_reg = getappdata(hfig,'thres_reg');
M = getappdata(hfig,'M');
gIX = getappdata(hfig,'gIX');
regchoice = getappdata(hfig,'regchoice');
if regchoice{1}==3, % from Centroid
    ctrdID = regchoice{2};
end

i = find(unique(gIX)==ctrdID);
if isempty(i),
    disp('input is empty!');beep;
    return;
end
[C,~] = FindCentroid(gIX,M);
regressor = C(i,:);
[cIX_,gIX_] = Regression_direct(M_0,thres_reg,regressor);
numC = length(cIX_);
regressor_last = regressor;

if ~isempty(cIX_),
    for itr = 1:20,
        [cIX_,gIX_] = Regression_direct(M_0,thres_reg,regressor_last);
        numC = length(cIX_);
        disp(num2str(numC));
        M = M_0(cIX_,:);
        regressor = FindCentroid(gIX_,M);
        if itr>1 && abs(numC-numC_last)<5,%corr(regressor_last',regressor')>thres_iter,
            disp(['converged at itr = ' num2str(itr)]);
            break;
        end
        regressor_last = regressor;
        numC_last = numC;
    end
    gIX_ = ceil((1:length(cIX_))'/length(cIX_)*min(20,length(cIX_)));
    UpdateIndices(hfig,cIX_,gIX_,40);
    RefreshFigure(hfig);
end
end

%% row 3: (TBD)

% function checkbox_dataFR_Callback(hObject,~) % obsolete !!!
% hfig = getParentFigure(hObject);
% SetDataFR(hfig,get(hObject,'Value'));
% RefreshFigure(hfig);
% end

%% ----- tab four ----- (Clustering etc.)

%% row 1: k-means

function edit_kmeans_Callback(hObject,~)
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');

str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d');
    numK = temp{:};
    disp(['k-means k=' num2str(numK) '...']);
    tic
    rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK < 10^7 && numK~=1,
        disp('Replicates = 5');
        [gIX,C,sumd,D] = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        [gIX,C,sumd,D] = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        [gIX,C,sumd,D] = kmeans(M,numK,'distance','correlation');%,'Replicates',3);
    end
    toc
    beep
    
    if numK>1,
        gIX = HierClusDirect(C,gIX,numK);
    end
    
    cIX = getappdata(hfig,'cIX');
    UpdateIndices(hfig,cIX,gIX,numK);
    RefreshFigure(hfig);
end
end

function edit_kmeans2_Callback(hObject,~) % based on anatomical distance
hfig = getParentFigure(hObject);
z_res = getappdata(hfig,'z_res');
cIX = getappdata(hfig,'cIX');
CInfo = getappdata(hfig,'CInfo');
M = getappdata(hfig,'M');

XY = reshape([CInfo(cIX).center],2,[])';
Z = [CInfo(cIX).slice]'.*z_res;
XYZ = horzcat(XY,Z);
Combo = horzcat(XYZ/50,M);% how to weight????????????????

str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d');
    numK = temp{:};
    disp(['anat. weighted k-means k=' num2str(numK) '...']);
    tic
    rng('default');% default = 0, but can try different seeds if doesn't converge
    [gIX,C,sumd,D] = kmeans(Combo,numK,'distance','correlation');%,'Replicates',3);
    toc
    beep
    
    if numK>1,
        gIX = HierClusDirect(C,gIX,numK);
    end

    cIX = getappdata(hfig,'cIX');
    UpdateIndices(hfig,cIX,gIX,numK);
    RefreshFigure(hfig);
end
end

function edit_kmeans_elbow_Callback(hObject,~)
hfig = getParentFigure(hObject);
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');

str = get(hObject,'String');
if ~isempty(str),
    str = strrep(str,'end',num2str(max(gIX)));
    range = ParseRange(str);

    Silh_mean = zeros(1,length(range));
    for i = 1:length(range),
        disp(['k-means k=' num2str(range(i))]);
        gIX = kmeans(M,range(i),'distance','correlation');        
        silh = silhouette(M,gIX,'correlation');
        Silh_mean(i) = mean(silh);
    end
    [~,ix] = max(Silh_mean);    
    numK = range(ix);
    
    figure;hold on
    plot(range,Silh_mean,'o-','color',[0.5 0.5 0.5]);
    
    if length(range)~=1,
        plot([numK,numK],[min(Silh_mean),max(Silh_mean)],'r--');
        xlim([range(1),range(end)]);
        ylim([min(Silh_mean),max(Silh_mean)])
    end
    xlabel('k = # of clusters')
    ylabel('silhouette mean')
    
    %% perform regular k-means
    disp(['k-means k=' num2str(numK)]);
    rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK < 10^7 && numK~=1,
        disp('Replicates = 5');
        [gIX,C,sumd,D] = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        [gIX,C,sumd,D] = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        [gIX,C,sumd,D] = kmeans(M,numK,'distance','correlation');%,'Replicates',3);
    end       
    if numK>1,
        gIX = HierClusDirect(C,gIX,numK);
    end    
    cIX = getappdata(hfig,'cIX');    
    UpdateIndices(hfig,cIX,gIX,numK);
    RefreshFigure(hfig);
    beep;
end
end

%% row 2: Auto-clustering

function pushbutton_merge_Callback(hObject,~)
% disp('merging...');
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');
U = unique(gIX);
numU = length(U);
[C,D] = FindCentroid(gIX,M);

thres_merge = getappdata(hfig,'thres_merge');

i = 1;
while i<numU,    
    c = corr(C(i,:)',C(i+1,:)');
    if c > thres_merge,
        IX = find(gIX == U(i+1));
        gIX(IX)=U(i);
        U = unique(gIX);
        numU = length(U);
        
        IX = find(gIX == U(i));
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
        C(i+1,:) = [];
        D(i+1) = [];
    else
        i = i+1;
    end
end

if numU>1,
    [gIX, numU] = HierClus(M,gIX);
end

UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
disp('merging complete');
end

function edit_mergethres_Callback(hObject,~)
str = get(hObject,'String');
temp = textscan(str,'%f',1);
hfig = getParentFigure(hObject);
setappdata(hfig,'thres_merge',temp{:});
end

function pushbutton_iter_split(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
numU = getappdata(hfig,'numK');
thres_split = getappdata(hfig,'thres_split');
M_0 = getappdata(hfig,'M_0_fluo');
classheader = getappdata(hfig,'classheader');

disp('iter. split all, beep when done...');
thres_size = 10;
thres_H = thres_split; 
% thres_H = [0.2;0.15;0.1;0.05]; % could have more rounds...

% initialization
I_rest = [];
% loop
tic
for round = 1:length(thres_H),
    disp(['round ' num2str(round) ', numU = ' num2str(numU)]);
    dthres = 1-thres_H(round);
    gIX_last = gIX;
    I_clean_last = cIX;
    cIX = [];
    gIX = [];
    for i = 1:numU,
        disp(['i = ' num2str(i)]);
        ix = find(gIX_last == i);
%         IX = ix;
        IX = I_clean_last(ix);
        M_s = M_0(IX,:);
        [I_rest,cIX,gIX,numU] = CleanClus(M_s,IX,I_rest,cIX,gIX,numU,dthres,thres_size);
%         cIX = I_clean_last(I_clean);
    end

    [gIX, numU] = SqueezeGroupIX(gIX);
    SaveClass(hfig,cIX,gIX,classheader,['clean_round' num2str(round)]);

    SaveClass(hfig,I_rest,ones(length(I_rest),1),classheader,...
       ['rest_round' num2str(round)]);
end
toc
beep
end

function [I_rest,I_clean,gIX_clean,numU] = CleanClus(M_s,IX,I_rest,I_clean,gIX_clean,numU,dthres,thres_size)
I_clean_s = [];

% find numK_s for kmeans
kmax = min(round(size(M_s,1)/thres_size),30);
% try numK_s = 1
numK_s = 1;
rng('default');
[gIX_s,~,~,D] = kmeans(M_s,numK_s,'distance','correlation');
Dist = min(D,[],2);
if mean(Dist)>dthres,
    % try numK_s = kmax
    numK_s = kmax;
    rng('default');
    [gIX_s,~,~,D] = kmeans(M_s,numK_s,'distance','correlation');
    Dist = min(D,[],2);
    if mean(Dist)<dthres, % find in between value for numK_s
        numK_s = 2;
        while 1,
            % kmeans-cluster by numK_s
            rng('default');
            [gIX_s,~,~,D] = kmeans(M_s,numK_s,'distance','correlation');
            Dist = min(D,[],2);
            if mean(Dist)<dthres,
                break;
            end
            
            if numK_s < kmax,
                numK_s = numK_s+1;
                disp(['numK_s = ' num2str(numK_s)]);
            else % numK_s = kmax;
                break;
            end
        end
    else disp(['numK_s = ' num2str(numK_s)]);
    end
end
% have numK_s that makes mean(Dist) < thres, or numK_s = kmax

for i = 1:numK_s,
    IX_s = find(gIX_s == i);
    if length(IX_s)>thres_size, % cluster big enough to start
        dst = Dist(IX_s);
        if mean(dst) < dthres,
            I_clean_s = [I_clean_s; IX_s];
            gIX_clean = [gIX_clean; gIX_s(IX_s)+double(numU)];
        else
            ix = find(dst<dthres); % clean
            if length(ix)>=thres_size, % clean cluster still big enough
                I = IX_s(ix);
                I_clean_s = [I_clean_s; I];
                gIX_clean = [gIX_clean; gIX_s(I)+double(numU)];
            end
        end
    end
end

numU = numU + numK_s;
I_rest_s = setdiff(1:size(M_s,1),I_clean_s);
I_rest = [I_rest; IX(I_rest_s)];
I_clean = [I_clean; IX(I_clean_s)];
end

function edit_splitthres_Callback(hObject,~)
str = get(hObject,'String');
temp = textscan(str,'%f',1);
hfig = getParentFigure(hObject);
setappdata(hfig,'thres_split',temp{:});
end

function pushbutton_autoclus_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');
M_0 = getappdata(hfig,'M_0_fluo');
classheader = getappdata(hfig,'classheader');
thres_size = 10;
thres_split = getappdata(hfig,'thres_split');
thres_stimlock = 1.0;
thres_merge = getappdata(hfig,'thres_merge');
thres_silh = 0.4;

isWkmeans = getappdata(hfig,'isWkmeans');

%% kmeans
if isWkmeans,
    numK = 20;
    disp(['kmeans k = ' num2str(numK)]);
    tic
    rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK < 10^7 && numK~=1,
        disp('Replicates = 5');
        gIX = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        gIX = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        gIX = kmeans(M,numK,'distance','correlation');
    end
    toc
end
[gIX, numU] = SqueezeGroupIX(gIX);

%% pushbutton_iter_split(hObject,~);
disp('iter. split all...');
I_rest = [];
iter = 1;
gIX_last = gIX;
I_clean_last = cIX;
cIX = [];
gIX = [];
for i = 1:numU,
    disp(['i = ' num2str(i)]);
    ix = gIX_last == i;
    IX = I_clean_last(ix);
    M_s = M_0(IX,:);
    [I_rest,cIX,gIX,numU] = CleanClus(M_s,IX,I_rest,cIX,gIX,numU,1-thres_split,thres_size);
end
[gIX, ~] = SqueezeGroupIX(gIX);
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
SaveClass(hfig,cIX,gIX,classheader,['clean_round' num2str(iter)]);
SaveClass(hfig,I_rest,ones(length(I_rest),1),classheader,...
    ['rest_round' num2str(iter)]);

[gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX);

%% rank by stim-lock
disp('stim-lock');
M = M_0(cIX,:);
[gIX,rankscore] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU);
disp('ranking complete');
% and threshold
IX = find(rankscore<thres_stimlock);
ix = ismember(gIX,IX);
gIX = gIX(ix);
cIX = cIX(ix);

%% ~ pushbutton_autoregclus_Callback(hObject,~);
thres_reg = getappdata(hfig,'thres_reg');
[cIX,gIX,~] = AllCentroidRegression_Direct(M_0,thres_reg,cIX,gIX);
disp('auto-reg-clus complete');

[gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX);

SaveClass(hfig,cIX,gIX,classheader,'clean_round2');

%% Silhouette
disp('silhouette analysis');
gIX_last = gIX;
for i = 1:numU,
    disp(['i = ' num2str(i)]);
    IX = find(gIX_last == i);
    cIX_2 = cIX(IX);
    M_s = M_0(cIX_2,:);
    % try k-means with k=2, see whether to keep
    gIX_ = kmeans(M_s,2,'distance','correlation');
    silh = silhouette(M_s,gIX_,'correlation');
    if mean(silh)>thres_silh,
        % keep the k-means k=2 subsplit
        disp('split');
        gIX(IX) = gIX_ + numU; % reassign (much larger) gIX
    end
end
[gIX, ~] = SqueezeGroupIX(gIX);

SaveClass(hfig,cIX,gIX,classheader,'clean_round3');

%% rank by stim-lock ?? bug?
% disp('stim-lock');
% M = M_0(cIX,:);
% [gIX,rankscore] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU);
% disp('ranking complete');
% % and threshold
% IX = find(rankscore<thres_stimlock);
% ix = ismember(gIX,IX);
% gIX = gIX(ix);
% cIX = cIX(ix);
% 
% [gIX, ~] = Merge_direct(thres_merge,M_0,cIX,gIX);

% size threshold
thres_size = getappdata(hfig,'thres_size');
[cIX, gIX, numU] = ThresSize(cIX,gIX,thres_size);

%% update GUI
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
beep;

end

function pushbutton_thressize_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
thres_size = getappdata(hfig,'thres_size');
[cIX, gIX, numU] = ThresSize(cIX,gIX,thres_size);
UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
end

function edit_sizethres_Callback(hObject,~)
str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%f');
    thres_size = temp{:};
end
hfig = getParentFigure(hObject);
setappdata(hfig,'thres_size',thres_size);
end

function [cIX, gIX, numU] = ThresSize(cIX,gIX,thres_size)
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_size,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX, numU] = SqueezeGroupIX(gIX);
end

function [gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX)
M = M_0(cIX,:);
[gIX, numU] = HierClus(M,gIX);
U = unique(gIX);
M = M_0(cIX,:);
[C,D] = FindCentroid(gIX,M);
i = 1;
while i<numU,    
    c = corr(C(i,:)',C(i+1,:)');
    if c > thres_merge,
        IX = find(gIX == U(i+1));
        gIX(IX)=U(i);
        U = unique(gIX);
        numU = length(U);
        
        IX = find(gIX == U(i));
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
        C(i+1,:) = [];
        D(i+1) = [];
    else
        i = i+1;
    end
end
[gIX, numU] = HierClus(M,gIX);
disp('merging complete');
end

function checkbox_wkmeans_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isWkmeans',get(hObject,'Value'));
end

%% row 3: misc plots

function pushbutton_hierplot_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');

[gIX2, numU] = HierClus(M,gIX,'isplotfig');
if ~isequal(gIX,gIX2),
    UpdateIndices(hfig,cIX,gIX2,numU);
    RefreshFigure(hfig);
end
end

function edit_hierpartn_Callback(hObject,~)
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');
gIX = getappdata(hfig,'gIX');
hierinplace = getappdata(hfig,'hierinplace');

[gIX, numU] = SqueezeGroupIX(gIX);
[C,~] = FindCentroid(gIX,M);

str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d',1);
    numCuts = temp{:};
    
%     tree = linkage(C,'average','correlation');
%     IX_tree = cluster(tree,'maxclust',numCuts);
    IX_tree = clusterdata(C,'criterion','distance','distance','correlation','maxclust',numCuts)
    
    if hierinplace,
        % sort to keep clusters in place
        U = unique(IX_tree,'stable');
        temp = zeros(size(IX_tree));
        for i = 1:length(U),
            temp(IX_tree==U(i)) = i;
        end
        IX_tree = temp;
    end
    
    % update gIX
    temp = zeros(size(gIX));
    for i = 1:numU,
        temp(gIX==i) = IX_tree(i);
    end
    gIX = temp;

    cIX = getappdata(hfig,'cIX');
    UpdateIndices(hfig,cIX,gIX,length(unique(IX_tree)));
    RefreshFigure(hfig);
end
end

function edit_hierpartthres_Callback(hObject,~)
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');
gIX = getappdata(hfig,'gIX');
hierinplace = getappdata(hfig,'hierinplace');

[gIX, numU] = SqueezeGroupIX(gIX);
[C,~] = FindCentroid(gIX,M);

str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%f',1);
    thres= temp{:};
    
    IX_tree = clusterdata(C,'criterion','distance',...
        'distance','correlation','cutoff',thres);
        
    if hierinplace,
        % sort to keep clusters in place
        U = unique(IX_tree,'stable');
        temp = zeros(size(IX_tree));
        for i = 1:length(U),
            temp(IX_tree==U(i)) = i;
        end
        IX_tree = temp;
    end
    
    % update gIX
    temp = zeros(size(gIX));
    for i = 1:numU,
        temp(gIX==i) = IX_tree(i);
    end
    gIX = temp;

    cIX = getappdata(hfig,'cIX');
    UpdateIndices(hfig,cIX,gIX,length(unique(IX_tree)));
    RefreshFigure(hfig);
end
end

function checkbox_hierinplace_Callback(hObject,~)
hierinplace = get(hObject,'Value');
hfig = getParentFigure(hObject);
setappdata(hfig,'hierinplace',hierinplace);
end

function pushbutton_corrplot_Callback(hObject,~)
hfig = getParentFigure(hObject);
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');
[C,~] = FindCentroid(gIX,M);
coeffs = corr(C');%corr(C(1,:)',C(2,:)')
figure('Position',[1000,200,500,500]);
CorrPlot(coeffs);
end

function CorrPlot(coeffs)
im = coeffs;

% red-white-blue colormap
cmap = zeros(64,3);
cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
minlim = -1; %min(min(im));
maxlim = 1; %max(max(im));

RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap

image(RGB); axis equal; axis tight;
set(gca,'XTick',1:length(im),'YTick',1:length(im));

for i = 1:size(im,2), % horizontal.. because of image axis
    for j = 1:size(im,1),
        text(i-0.3, j, num2str(round(im(i,j)*100)/100));%, 'Units', 'data')
    end
end
end

function RGB = ImageToRGB(im,cmap,minlim,maxlim)
L = size(cmap,1);
ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
end

function pushbutton_plot4x4_Callback(hObject,~)
hfig = getParentFigure(hObject);
PlotBy16StimsFromGUI(hfig);
end

%% ----- tab five ----- (Saved Clusters)

%% row 1: Class

function popup_classmenu_Callback(hObject,~)
classID = get(hObject,'Value') - 1;
if classID>0,
    hfig = getParentFigure(hObject);
    UpdateClassID(hfig,classID); % update popupmenu
end
end

function edit_editclassname_Callback(hObject,~)
str = get(hObject,'String');

hfig = getParentFigure(hObject);
Class = getappdata(hfig,'Class');
classID = getappdata(hfig,'classID');
Class(classID).name = str;
setappdata(hfig,'Class',Class);
UpdateClassID(hfig,classID); % update popupmenu
end

function pushbutton_saveclass_Callback(hObject,~)
hfig = getParentFigure(hObject);
SaveClass(hfig);
end

function edit_class_header_Callback(hObject,~)
% get/format range
str = get(hObject,'String');

hfig = getParentFigure(hObject);
setappdata(hfig,'classheader',[str ':']);
end

function edit_newclassname_Callback(hObject,~)
str = get(hObject,'String');
hfig = getParentFigure(hObject);
setappdata(hfig,'newclassname',str);
end

function pushbutton_makeclass_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
classheader = getappdata(hfig,'classheader');
newclassname = getappdata(hfig,'newclassname');
SaveClass(hfig,cIX,gIX,classheader,newclassname);
end

function pushbutton_delclass_Callback(hObject,~)
choice = questdlg('Delete current class?','','Cancel','Yes','Yes');
if strcmp(choice,'Yes'), 
    hfig = getParentFigure(hObject);
    Class = getappdata(hfig,'Class');
    classID = getappdata(hfig,'classID');
    disp(['delete ' num2str(classID)]);
    Class(classID) = [];
    setappdata(hfig,'Class',Class);
        
    classID = max(1,classID-1);
    disp(['new classID: ' num2str(classID)]);
    UpdateClassID(hfig,classID);
end
end

%% row 2: Cluster

function popup_clusgroupmenu_Callback(hObject,~)
hfig = getParentFigure(hObject);
clusgroupID = getappdata(hfig,'clusgroupID');
new_clusgroupID = get(hObject,'Value');
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID);
end

function popup_clusmenu_Callback(hObject,~)
clusID = get(hObject,'Value') - 1;
if clusID>0,
    hfig = getParentFigure(hObject);
    UpdateClusID(hfig,clusID);
end
end

function edit_editclusname_Callback(hObject,~)
str = get(hObject,'String');

hfig = getParentFigure(hObject);
Cluster = getappdata(hfig,'Cluster');
clusID = getappdata(hfig,'clusID');
Cluster(clusID).name = str;
setappdata(hfig,'Cluster',Cluster);
UpdateClusID(hfig,clusID); % update popupmenu
end

function pushbutton_saveclus_Callback(hObject,~)
hfig = getParentFigure(hObject);
SaveCluster(hfig,'current');
end

function edit_clus_header_Callback(hObject,~)
% get/format range
str = get(hObject,'String');

hfig = getParentFigure(hObject);
setappdata(hfig,'clusheader',[str ':']);
end

function edit_newclusname_Callback(hObject,~)
str = get(hObject,'String');
hfig = getParentFigure(hObject);
setappdata(hfig,'newclusname',str);
end

function pushbutton_makeclus_Callback(hObject,~)
hfig = getParentFigure(hObject);
SaveCluster(hfig,'new');
end

function edit_setrank_Callback(hObject,~)
str = get(hObject,'String');
C = textscan(str,'%d');
rank = C{:};

hfig = getParentFigure(hObject);
Cluster = getappdata(hfig,'Cluster');
clusID = getappdata(hfig,'clusID');
% insert current cluster into new position = 'rank'
if rank > 1,
    temp = Cluster(clusID);
    Cluster(clusID) = [];
    Cluster = [Cluster(1:rank-1),temp,Cluster(rank:end)];
else 
    temp = Cluster(clusID);
    Cluster(clusID) = [];
    Cluster = [temp,Cluster(rank:end)];
end
setappdata(hfig,'Cluster',Cluster);
% update pointer = clusID
clusID = rank;
UpdateClusID(hfig,clusID);
end

function edit_notes_Callback(hObject,~)
str = get(hObject,'String');
hfig = getParentFigure(hObject);
Cluster = getappdata(hfig,'Cluster');
clusID = getappdata(hfig,'clusID');
Cluster(clusID).notes = str;
setappdata(hfig,'Cluster',Cluster);
end

function pushbutton_delclus_Callback(hObject,~)
choice = questdlg('Delete current cluster?','','Cancel','Yes','Yes');
if strcmp(choice,'Yes'), 
    hfig = getParentFigure(hObject);
    Cluster = getappdata(hfig,'Cluster');
    clusID = getappdata(hfig,'clusID');
    disp(['delete ' num2str(clusID)]);
    Cluster(clusID) = [];
    setappdata(hfig,'Cluster',Cluster);
        
    clusID = max(1,clusID-1);
    disp(['new clusID: ' num2str(clusID)]);
    UpdateClusID(hfig,clusID);
end
end

%% row 3: misc

function edit_clusUnion_Callback(hObject,~)
disp('union processing...');
hfig = getParentFigure(hObject);
M_0 = GetM_0(hfig);
Cluster = getappdata(hfig,'Cluster');
str = get(hObject,'String');
if ~isempty(str),
    % get/format range
    str = strrep(str,':','-'); % e.g. str= '1,3,5:8';
    C = textscan(str,'%d','delimiter',',');
    m = C{:};
    range = [];
    for i = 1:length(m),
        if m(i)>0,
            range = [range,m(i)]; %#ok<AGROW>
        else % have '-'sign,
            range = [range,m(i-1)+1:-m(i)]; %#ok<AGROW>
        end
    end

    % combine
    CIX = vertcat(Cluster(range).cIX);
    GIX = []; % gIX to match A
    for i = 1:length(range),
        GIX = [GIX; ones(length(Cluster(range(i)).gIX),1)*i];
    end    
    [cIX,gIX,numK] = SmartUnique(CIX,GIX,M_0(CIX,:));    
    UpdateIndices(hfig,cIX,gIX,numK);    
    RefreshFigure(hfig);
end
disp('union complete');
end

function [cIX,gIX,numK] = SmartUnique(CIX,GIX,M) % input: simply concatenated groups
disp('unique based on corr coeff...');
CTRD = FindCentroid(GIX,M); 

[uCIX,ia,~] = unique(CIX);
uGIX = GIX(ia);
if length(uCIX) == 1, 
    counts = length(x);
else counts = hist(CIX,uCIX); 
end
IX1 = find(counts==1);
CIX1 = uCIX(IX1); % single occurence
GIX1 = uGIX(IX1);
ia = find(~ismember(CIX,CIX1));
C = CIX(ia); % multiple copies/occurence, keep all copies
G = GIX(ia);
M_s = M(ia,:);
Coeff = zeros(length(ia),1);
for i = 1:length(ia),
    reg = CTRD(G(i),:);
    Coeff(i) = corr(M_s(i,:)',reg');
end
[~,I] = sort(Coeff,'descend'); % rank by coeff, so unique ('stable') gets highest
C = C(I); % finish sorting
G = G(I);
[CIX2,ia,~] = unique(C,'stable'); % keep the order!
GIX2 = G(ia); % the chosen copy for those with multiple copies
cIX = vertcat(CIX1,CIX2);
gIX = vertcat(GIX1,GIX2);
[gIX, numK] = SqueezeGroupIX(gIX);
disp('smart union complete');
end

function pushbutton_newclusgroup_Callback(hObject,~)
hfig = getParentFigure(hObject);
ClusGroup = getappdata(hfig,'ClusGroup');
ClusGroup = [ClusGroup,cell(1,1)];
setappdata(hfig,'ClusGroup',ClusGroup);

clusgroupID = getappdata(hfig,'clusgroupID');
new_clusgroupID = length(ClusGroup);
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID);
end

function pushbutton_delclusgroup_Callback(hObject,~)
choice = questdlg('Delete current Cluster-group?','','Cancel','Yes','Yes');
if strcmp(choice,'Yes'),
    hfig = getParentFigure(hObject);
    ClusGroup = getappdata(hfig,'ClusGroup');
    clusgroupID = getappdata(hfig,'clusgroupID');
    ClusGroup(clusgroupID) = [];
    setappdata(hfig,'ClusGroup',ClusGroup);
    UpdateClusGroupID(hfig,[],max(1,clusgroupID-1));
end
end

%% Internal functions

function SetDataFR(hfig,dataFR)
cIX = getappdata(hfig,'cIX');
setappdata(hfig,'dataFR',dataFR);
if dataFR,
    M_0_fluo = getappdata(hfig,'M_0_fluo');
    setappdata(hfig,'M',M_0_fluo(cIX,:));
else % Unchecked
    M_0_reg = getappdata(hfig,'M_0_reg'); % obsolete !!!
    setappdata(hfig,'M',M_0_reg(cIX,:));
end
global hdataFR;
set(hdataFR,'Value',dataFR);
end

function SaveClass(hfig,cIX,gIX,classheader,newclassname)
if ~exist('cIX','var'),
    cIX = getappdata(hfig,'cIX');
end
if ~exist('gIX','var'),
    gIX = getappdata(hfig,'gIX');
end

Class = getappdata(hfig,'Class');
dataname = getappdata(hfig,'dataname');

if ~exist('newclassname','var'),
    classID = getappdata(hfig,'classID');
else
    classID = numel(Class)+1;
    Class(classID).name = [classheader newclassname];
end

Class(classID).cIX = cIX;
Class(classID).gIX = gIX;
Class(classID).numel = length(cIX);
U = unique(gIX);
numU = length(U);
Class(classID).numK = numU;
Class(classID).dataname = dataname;

setappdata(hfig,'Class',Class);
UpdateClassID(hfig,classID); 
disp('class saved');
end

function UpdateClassID(hfig,classID,norefresh)
Class = getappdata(hfig,'Class');
numK = Class(classID).numK;
i_fish = getappdata(hfig,'i_fish');
% save
setappdata(hfig,'classID',classID);
% update GUI
global hclassname hclassmenu;
if ishandle(hclassname),
    set(hclassname,'String',Class(classID).name);
    temp = MakeMenu({Class.name});
    set(hclassmenu, 'String', temp);
    set(hclassmenu, 'Value', classID+1);
end

if ~exist('norefresh','var'),
    UpdateIndices(hfig,Class(classID).cIX,Class(classID).gIX,numK);
    RefreshFigure(hfig);
end
global VAR;
VAR(i_fish).Class = getappdata(hfig,'Class');
end

function UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID,norefresh)
% save/update old Cluster into ClusGroup before exiting, 
% as Cluster is the variable handled in hfig but not saved elsewhere
ClusGroup = getappdata(hfig,'ClusGroup');
Cluster = getappdata(hfig,'Cluster');
i_fish = getappdata(hfig,'i_fish');

ClusGroup{clusgroupID} = Cluster;
setappdata(hfig,'ClusGroup',ClusGroup);
global VAR;
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);

Cluster = ClusGroup{new_clusgroupID};
setappdata(hfig,'Cluster',Cluster);
setappdata(hfig,'clusgroupID',new_clusgroupID);

global hclusgroupmenu;
if ishandle(hclusgroupmenu),
    num = length(ClusGroup);
    menu = {}; for i = 1:num, menu = [menu,{num2str(i)}];end
    set(hclusgroupmenu,'String', menu,'Value',new_clusgroupID);
end

if ~exist('norefresh','var'),    
    if numel(Cluster) == 0,
        SaveCluster(hfig,'new');
    else
        clusID = 1;
        UpdateClusID(hfig,clusID);
    end
end
end

function ClusGroup = CurrentClusGroup(hfig)
ClusGroup = getappdata(hfig,'ClusGroup');
Cluster = getappdata(hfig,'Cluster');
clusgroupID = getappdata(hfig,'clusgroupID');
ClusGroup{clusgroupID} = Cluster;
setappdata(hfig,'ClusGroup',ClusGroup);
end

function [Cluster,clusID] = SaveCluster(hfig,state,clusheader,name)
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');

Cluster = getappdata(hfig,'Cluster');
dataname = getappdata(hfig,'dataname');

if strcmp(state,'current'),
    clusID = getappdata(hfig,'clusID');
else %if strcmp(state,'new'),
    clusID = numel(Cluster)+1;
    if ~exist('name','var'),
        name = getappdata(hfig,'newclusname');
    end
    if ~exist('clusheader','var'),
        clusheader = getappdata(hfig,'clusheader');
    end
    Cluster(clusID).name = [clusheader name];
end

Cluster(clusID).cIX = cIX;
Cluster(clusID).gIX = gIX;
Cluster(clusID).numel = length(cIX);
U = unique(gIX);
numU = length(U);
Cluster(clusID).numK = numU;
Cluster(clusID).dataname = dataname;

setappdata(hfig,'Cluster',Cluster);
UpdateClusID(hfig,clusID); 
disp('cluster saved');
end

function UpdateClusID(hfig,clusID)
Cluster = getappdata(hfig,'Cluster');
    
% save
setappdata(hfig,'clusID',clusID);
% update GUI
global hclusname hclusmenu;
set(hclusname,'String',Cluster(clusID).name);
menu = MakeMenu({Cluster.name});
set(hclusmenu,'String', menu,'Value',clusID+1);
numK = Cluster(clusID).numK;

UpdateIndices(hfig,Cluster(clusID).cIX,Cluster(clusID).gIX,numK);
RefreshFigure(hfig);
end

function menu = MakeMenu(name) % careful: input {Class.name}, otherwise only get 1 somehow
menu = [{'(choose)'},name];
for j=2:length(menu),menu(j)={[num2str(j-1) ': ' menu{j}]};end
end

function [gIX, numU] = HierClus(M,gIX,isplotfig) %#ok<INUSD>
[gIX, numU] = SqueezeGroupIX(gIX);
[C,~] = FindCentroid(gIX,M);
D = pdist(C,'correlation');
tree = linkage(C,'average','correlation');
leafOrder = optimalleaforder(tree,D);

if numU>1,
    if exist('isplotfig','var'),
            figure('Position',[100 400 1300 400]);
            subplot(1,3,1);
            CORR = corr(C');
            CorrPlot(CORR);
            
            subplot(1,3,2);
            dendrogram(tree,numU,'orientation','right','reorder',leafOrder);
            set(gca,'YDir','reverse');
            set(gca,'XTick',[]);
%             colormap(jet);

            subplot(1,3,3);
            C2 = C(leafOrder,:);
            CORR2 = corr(C2');
            CorrPlot(CORR2);
%             imagesc(CORR2);axis equal;axis tight
    end
    % sort for uniform colorscale
    temp = zeros(size(gIX));
    for i = 1:numU,
        temp(gIX==leafOrder(i)) = i; % = T(i) for clusters segmented from tree
    end
    gIX = temp;
end
end

function gIX = HierClusDirect(C,gIX,numU)
D = pdist(C,'correlation');
tree = linkage(C,'average','correlation');
leafOrder = optimalleaforder(tree,D);

% sort for uniform colorscale
temp = zeros(size(gIX));
for i = 1:numU,
    temp(gIX==leafOrder(i)) = i; % = T(i) for clusters segmented from tree
end
gIX = temp;
end

function [gIX, numK] = SqueezeGroupIX(gIX)
U = unique(gIX);
numK = length(U);
for i = 1:numK,
    old = U(i);
    gIX(gIX==old) = i;
end
end

% frequently used, updates cell-index,group-index,cluster-number. set-operations included in here.
function UpdateIndices(hfig,cIX,gIX,numK)
global hback hopID;
M_0 = GetM_0(hfig);
if ~exist('gIX','var'),
    gIX = getappdata(hfig,'gIX');
end
if ~exist('cIX','var'),
    cIX = getappdata(hfig,'cIX');
end

% update cache
bC = getappdata(hfig,'bCache');
cIX_last = getappdata(hfig,'cIX');
gIX_last = getappdata(hfig,'gIX');
if ~(isequal(cIX_last,cIX) && isequal(gIX_last,gIX)),
    bC.cIX = [cIX_last,bC.cIX];
    bC.gIX = [gIX_last,bC.gIX];
    bC.numK = [getappdata(hfig,'numK'),bC.numK];
    set(hback,'enable','on');
    if length(bC.cIX)>20,
        bC.cIX(end) = [];
        bC.gIX(end) = [];
        bC.numK(end) = [];
    end
end

% set operations, if applicable
opID = getappdata(hfig,'opID');
if opID ~= 0,
    switch opID,
        case 1,
            disp('union');
            [~,ia,ib] = union(cIX_last,cIX,'stable');
            IX = vertcat(cIX_last(ia),cIX(ib));% IX;
        case 2,
            disp('intersect');
            [IX,ia,~] = intersect(cIX_last,cIX);
            ib = [];
        case 3,
            disp('setdiff');
            [IX,ia] = setdiff(cIX_last,cIX);
            ib = [];
        case 4,
            disp('rev setdiff');
            % swap sequence, then same as opID==3
            temp = cIX;
            cIX = cIX_last;
            cIX_last = temp;
            temp = gIX;
            gIX = gIX_last;
            gIX_last = temp;
            [IX,ia] = setdiff(cIX_last,cIX);
            ib = [];
        case 5,
            disp('setxor');
            [IX,ia,ib] = setxor(cIX_last,cIX);
        case 6,
            disp('smartUnion');
            CIX = vertcat(cIX_last,cIX);
            GIX = [gIX_last;gIX+max(gIX_last)]; % gIX to match
            [cIX,gIX,numK] = SmartUnique(CIX,GIX,M_0(CIX,:));              
    end
    if opID<6,
        if ~isempty(IX),
            cIX = IX;
            gIX = vertcat(gIX_last(ia),gIX(ib)+max(gIX_last(ia)));
            numK = length(unique(gIX));
            %         [gIX, numK] = SqueezeGroupIX(gIX);
        else
            errordlg('operation result is empty set!')
            waitforbuttonpress;
        end
    end
    set(hopID,'Value',1,'BackgroundColor',[1,1,1]); % reset
    setappdata(hfig,'opID',0);
end

if length(cIX)*size(M_0,2)<5*10^8,
    % set M
    M = M_0(cIX,:);
    
    setappdata(hfig,'M',M);
    setappdata(hfig,'bCache',bC);
    setappdata(hfig,'cIX',cIX);
    setappdata(hfig,'gIX',gIX);
    if exist('numK','var'),
        setappdata(hfig,'numK',double(numK));
    end
else
    errordlg('dataset too large!')
    waitforbuttonpress;
end

% handle rankID: >=2 means write numbers as text next to colorbar
% first UpdateIndices sets rankID to 100, second sets back to 0 
rankID = getappdata(hfig,'rankID');
if rankID>=2,
    if rankID==100,
        setappdata(hfig,'rankID',0);
    else
        setappdata(hfig,'rankID',100);
    end
end
end

% frequently used, 2 plotting functions are outside ('DrawClusters.m' and 'DrawClustersOnMap_LSh.m')
function RefreshFigure(hfig)
CInfo = getappdata(hfig,'CInfo');
M = getappdata(hfig,'M');
dataFR = getappdata(hfig,'dataFR');
fictive = getappdata(hfig,'fictive');
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');

numK = getappdata(hfig,'numK');
stim = getappdata(hfig,'stim');
anat_yx = getappdata(hfig,'anat_yx');
anat_yz = getappdata(hfig,'anat_yz');
anat_zx = getappdata(hfig,'anat_zx');
isCentroid = getappdata(hfig,'isCentroid');
clrmap = getappdata(hfig,'clrmap');
rankscore = getappdata(hfig,'rankscore');
rankID = getappdata(hfig,'rankID');
iswrite = (rankID>=2);
isPlotLines = 0; %getappdata(hfig,'isPlotLines');
isPlotFictive = 1; %getappdata(hfig,'isPlotFictive');

isPopout = 0;

if isempty(cIX),
    errordlg('empty set! go back!');
    return;
end

allAxesInFigure = findall(hfig,'type','axes');
if ~isempty(allAxesInFigure)
    delete(allAxesInFigure);
end

figure(hfig);
% left subplot
h1 = axes('Position',[0.05, 0.04, 0.55, 0.83]);
% if isCentroid,
%     [C,~] = FindCentroid(gIX,M);
%     DrawClusters(h1,C,unique(gIX),dataFR,numK,stim,fictive,clrmap,rankscore,iswrite);
% else
%     DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,iswrite);
% end
if isCentroid,
    [C,~] = FindCentroid(gIX,M);
    DrawClusters(h1,C,unique(gIX),dataFR,numK,stim,fictive,clrmap,rankscore,...
        iswrite,isPopout,isPlotLines,isPlotFictive);
%     DrawClusters(h1,C,unique(gIX),dataFR,numK,stim,fictive,clrmap,rankscore,iswrite,ispopout);
else
    DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,...
        iswrite,isPopout,isPlotLines,isPlotFictive);
%     DrawClusters(h1,M,gIX,dataFR,numK,stim,fictive,clrmap,rankscore,iswrite,ispopout);
end

% right subplot
h2 = axes('Position',[0.61, 0.04, 0.35, 0.83]);
DrawClustersOnMap_LSh(CInfo,cIX,gIX,numK,anat_yx,anat_yz,anat_zx,clrmap);
end

function M_0 = GetM_0(hfig)
dataFR = getappdata(hfig,'dataFR');
if dataFR,
    M_0 = getappdata(hfig,'M_0_fluo');
else
    M_0 = getappdata(hfig,'M_0_reg');
end
end

function [C,D] = FindCentroid(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
    end
end
end

function closefigure_Callback(hfig,~)
global EXPORT_autorecover;
EXPORT_autorecover = getappdata(hfig);
end

function fig = getParentFigure(fig)
% if the object is a figure or figure descendent, return the figure. Otherwise return [].
while ~isempty(fig) && ~strcmp('figure', get(fig,'type'))
  fig = get(fig,'parent');
end
end

function runscript(flag_script,var_script)
switch flag_script
    case 'push_cIX_gIX'
        UpdateIndices(var_script{:});
        RefreshFigure(var_script{1});        
end
end
% end
