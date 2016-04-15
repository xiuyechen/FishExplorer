%%% Interactive app for exploratory analysis of calcium imaging data
% (with stimulus, behavior, and anatomy)
%{
Input calcium data: 1 trace per cell/ROI, ~50,000 cells per fish
load collection of cells from multiple fish, or load full data of single fish individually
main outputs: GUI plots, clusters saved into .mat, export variables to MATLAB workspace

Tip: to see the structure of this code, use '(Right click -> )Code Folding\Fold all'
(or hotkey) to collapse all cells.
UI controls are organized by tabs and then by rows, instructions and
comments are where they are constructed ('User Interface:' -> function hfig... ->)
General internal functions are at the end, some specialized .m functions are outside.

Written in Matlab R2014b/R2015a running on Windows 7.

- Xiuye Chen (xiuyechen@gmail.com), Engert Lab, Spring 2015


To start, run the following in separate script or command window:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;clc
global VAR;
load('VAR_current.mat','VAR');
load('CONSTs_current.mat','CONSTs');
data_dir = [pwd '\example data'];
[hfig,fcns] = GUI_FishExplorer(data_dir,CONSTs,VAR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Full datasets are stored as 'CONST_F?.mat' or 'CONST_F?_fast.mat' (? = fish number), to load from the GUI
%}

%% User Interface:
function [hfig,fcns] = GUI_FishExplorer(flag_script,var_script)
global data_masterdir name_MASKs name_ReferenceBrain VAR;

if exist('flag_script','var')
    if ~exist('var_script','var')
        var_script={};
    end
    runscript(flag_script,var_script);
    return
end

%% Make figure
scrn = get(0,'Screensize');
hfig = figure('Position',[scrn(3)*0.2 scrn(4)*0.05 scrn(3)*0.75 scrn(4)*0.86],...% [50 100 1700 900]
    'Name','GUI_LSh','DeleteFcn',@closefigure_Callback,...
    'ToolBar', 'none'); % 'MenuBar', 'none'
hold off; axis off

%% Make menu
global hm1;
hm1 = uimenu(hfig,'Label','My File');
hm1_1 = uimenu(hm1,'Label','Quick save to workspace');
hm1_2 = uimenu(hm1,'Label','Save to file (default path)');

%% Folder setup
% directory for full fish data (.mat)
setappdata(hfig,'data_masterdir',data_masterdir);

%% Pass external variables into appdata (stored with main figure handle)
setappdata(hfig,'VAR',VAR);
nFish = length(VAR);

%% Load ZBrain Atlas  
if exist(fullfile(data_masterdir,name_MASKs),'file') ... 
    && exist(fullfile(data_masterdir,name_ReferenceBrain),'file'),

    % load masks for ZBrain Atlas
    MASKs = load(fullfile(data_masterdir,name_MASKs));
    setappdata(hfig,'MASKs',MASKs);
    % load reference brain image stack with the 3 projections
    load(fullfile(data_masterdir,name_ReferenceBrain));
    setappdata(hfig,'anat_stack_norm',anat_stack_norm);
    setappdata(hfig,'anat_yx_norm',anat_yx_norm);
    setappdata(hfig,'anat_yz_norm',anat_yz_norm);
    setappdata(hfig,'anat_zx_norm',anat_zx_norm);
end

%% fish protocol sets (different sets have different parameters)
M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2]; % M = Matrix
setappdata(hfig,'M_fish_set',M_fish_set);

%% parameters / constants
setappdata(hfig,'z_res',19.7); % resoltion in z, um per slice
% fpsec = 1.97; % hard-coded in ext function 'GetStimRegressor.m'
% approx fpsec of 2 used in ext function 'DrawCluster.m'

% cache
bCache = []; % Cache for going b-ack (bad abbr)
fCache = []; % Cache for going f-orward
bCache.cIX = cell(1,1);
bCache.gIX = cell(1,1);
bCache.numK = cell(1,1);
fCache.cIX = cell(1,1);
fCache.gIX = cell(1,1);
fCache.numK = cell(1,1);
setappdata(hfig,'bCache',bCache);
setappdata(hfig,'fCache',fCache);

%% Initialize internal params into appdata

% thresholds
thres_merge = 0.9;
thres_split = 0.7;
thres_reg = 0.7;
thres_size = 10;
thres_ttest = 0.001;
setappdata(hfig,'thres_merge',thres_merge);
setappdata(hfig,'thres_split',thres_split); % for function 'pushbutton_iter_split'
setappdata(hfig,'thres_reg',thres_reg); % regression threshold, ~correlation coeff
setappdata(hfig,'thres_size',thres_size); % min size for clusters
setappdata(hfig,'thres_ttest',thres_ttest); % min size for clusters

% variables
% (not sure all these need to be initialized, probably not complete either)
setappdata(hfig,'clrmap','hsv');
setappdata(hfig,'opID',0);
setappdata(hfig,'rankID',0);
setappdata(hfig,'isPlotLines',0);
setappdata(hfig,'isPlotBehavior',1);
setappdata(hfig,'rankscore',[]);
setappdata(hfig,'isCentroid',0);
setappdata(hfig,'isWkmeans',1); % in autoclustering, with/without kmeans
setappdata(hfig,'isflipstim',0); % flag for flip updown %%%%%%obsolete ????????????????????????????????????
setappdata(hfig,'regchoice',{1,1}); % regressor choice; first cell,1=stim,2=motor,3=centroid
% setappdata(hfig,'isfullfish',0); % no if QuickUpdateFish, yes if LoadFullFish
setappdata(hfig,'isPlotCorrHist',0); % option for regression
setappdata(hfig,'isPlotReg',1); % plot regressor when selecting it
setappdata(hfig,'hierinplace',1); % hier. partitioning, no reordering
setappdata(hfig,'isAvr',1); % show average/full stimulus
setappdata(hfig,'isRawtime',0); % show stimulus in original order or sorted
setappdata(hfig,'isZscore',0); % show normalized (z-score) version of fluorescent
setappdata(hfig,'isShowMasks',0);
setappdata(hfig,'isShowMskOutline',0);
setappdata(hfig,'isWeighAlpha',0);
setappdata(hfig,'isPlotAnatomyOnly',0);
setappdata(hfig,'isRefAnat',0);

%% Create UI controls
set(gcf,'DefaultUicontrolUnits','normalized');
set(gcf,'defaultUicontrolBackgroundColor',[1 1 1]);

% tab group setup
tgroup = uitabgroup('Parent', hfig, 'Position', [0.05,0.88,0.91,0.12]);
numtabs = 6;
tab = cell(1,numtabs);
M_names = {'General','Operations','Regression','Clustering etc.','Saved Clusters','Atlas'};
for i = 1:numtabs,
    tab{i} = uitab('Parent', tgroup, 'BackgroundColor', [1,1,1], 'Title', M_names{i});
end

% grid setup, to help align display elements
rheight = 0.2;
yrow = 0.7:-0.33:0;%0.97:-0.03:0.88;
dTextHt = 0.05; % dTextHt = manual adjustment for 'text' controls:
% (vertical alignment is top instead of center like for all other controls)
bwidth = 0.03;
grid = 0:bwidth+0.001:1;

%% global variables: various UI element handles
global hback hfwd hclusgroupmenu hclusgroupname hclusmenu hclusname...
    hstimrangemenu hopID hloadfish hstimreg hmotorreg...
    hcentroidreg hcentroid hstimrange hmasklistbox h_isavr h_israwtime h_iszscore; % hfishnum

%% UI ----- tab one ----- (General)
i_tab = 1;

%% UI row 1: File
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2; % saves clusters to workspace (global) variable 'VAR'
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Quick save',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_save_Callback);

i=i+n;
n=2; % saves both to workspace and to 'VAR_current.mat' and to arc folder
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
n=2; % plots selected cells on anatomy z-stack, tiled display
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Tile Clusters',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_tileClusters_Callback);

i=i+n;
n=2; % export main working variables to workspace, can customize!
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Export(workspace)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_exporttoworkspace_Callback);

i=i+n;
n=2; % 
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Import(VAR)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_loadVARfromworkspace_Callback);

i=i+n;
n=2; % 
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Import(current)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_loadCurrentClustersfromworkspace_Callback);

i=i+n;
n=2; % create popup figure without the GUI components, can save manually from default menu
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Popup plot',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_popupplot_Callback);

i=i+n;
n=2; % popupplot option: whether to plot cluster mean lines instead of all raw traces
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot lines',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_isPlotLines_Callback);

i=i+n;
n=2; % popupplot option: whether to plot behavior bar
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot behavior','Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_isPlotBehavior_Callback);

i=i+n;
n=2; % popupplot option: whether to only plot anatomy map (right half)
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot anatomy only',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_isPlotAnatomyOnly_Callback);

%% UI row 2: Load
i_row = 2;
i = 1;n = 0;

% i=i+n;
% n=2; % this design is underused now... Quick-load only depends on CONSTs,
% % which is a minimum collection of clusters from all fish, so you can load
% % the program without full single-fish data. eventually can use this
% % platform to do things across fish (like after anatomical alignment).
% uicontrol('Parent',tab{i_tab},'Style','text','String','Quick-load fish:',...
%     'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');
%
% i=i+n;
% n=1; % loads 'CONSTs_current.mat' from current directory
% temp = {}; for j = 1:nFish, temp = [temp,{num2str(j)}];end
% hfishnum = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
%     'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
%     'Callback',@popup_quickloadfishmenu_Callback);

i=i+n;
n=2; % loads full single-fish data from CONST_F?.mat
uicontrol('Parent',tab{i_tab},'Style','text','String','Load fish #:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=2; %
temp = {}; for j = 1:nFish, temp = [temp,{num2str(j)}];end
temp = [{'(choose)'},temp];
hloadfish = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',temp,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_loadfullfishmenu_Callback);

i=i+n;
n=2; %
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','or choose files',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_choosefilestoload_Callback);

i=i+n;
n=2; % options to load different stimulus types (if applicable for this fish)
uicontrol('Parent',tab{i_tab},'Style','text','String','Stim type:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=3;
hstimrangemenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String','(empty)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_stimrangemenu_Callback);

%% UI row 3: Display
i_row = 3;
i = 1;n = 0;

i=i+n;
n=2; % only centroids (~mean) of clusters shown on left-side plot, the rest is unchanged
h_isavr = uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show stim-avr',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',1,...
    'Callback',@checkbox_isStimAvr_Callback);

i=i+n;
n=2; % only centroids (~mean) of clusters shown on left-side plot, the rest is unchanged
h_israwtime = uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show raw-time',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',0,...
    'Callback',@checkbox_isRawtime_Callback);

i=i+n;
n=2; % showing z-scored version (each cell normalized to mean=0, std=1) on left-side plot
h_iszscore = uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show z-score',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',0,...
    'Callback',@checkbox_isZscore_Callback);

i=i+n;
n=2; % choose stimulus range - use numbers indicated in stimrangemenu % (eg 1:2,3-5)
uicontrol('Parent',tab{i_tab},'Style','text','String','Stim range:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
hstimrange = uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'String','1:end',...
    'Callback',@edit_stimrange_Callback);

i=i+n;
n=3; % only centroids (~mean) of clusters shown on left-side plot, the rest is unchanged
hcentroid = uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show cluster-centroids',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_showcentroids_Callback);

i=i+n;
n=3; % only centroids (~mean) of clusters shown on left-side plot, the rest is unchanged
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show normalized anatomy',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@checkbox_showrefanat_Callback);

%% UI ----- tab two ----- (Operations)
i_tab = 2;

%% UI row 1: Range
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2; % saves up to 20 steps backwards (datatype/stimrangemenu change does not count)
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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_choose_range_Callback);

i=i+n;
n=2; % Choose range of clusters to exclude. format: e.g. '1:2,4-6,8:end'
uicontrol('Parent',tab{i_tab},'Style','text','String','Exclude:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_exclude_range_Callback);

i=i+n;
n=1; % Choose range of clusters to fuse/combine into single cluster. format: e.g. '1:2,4-6,8:end'
uicontrol('Parent',tab{i_tab},'Style','text','String','Fuse:',... % (eg 1:2,3-5)
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_fuse_range_Callback);

%% UI row 2: Operations
i_row = 2;
i = 1;n = 0;

i=i+n;
n=2; % operates between the current cell selection and the next (in this order).
uicontrol('Parent',tab{i_tab},'Style','text','String','Set operations:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n; % 'setdiff' is current minus next, 'rev setdiff' is next minus current.
n=2; % smartUnion = SmartUnique, cells belonging to 2 clusters goes to the more correlated one
menu = {'(choose)','union','intersect','setdiff','rev setdiff','setxor','smartUnion'};
hopID = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@popup_operations_Callback});

i=i+n;
n=2; % rank clusters based on various criteria (to choose)
uicontrol('Parent',tab{i_tab},'Style','text','String','Rank by:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n; % 'hier' is the same as default (used after every k-means);'stim-lock' uses std across reps;
n=2; % motor stuff uses the best alignment (by cross-correlation) with the behavior trace;
% L+R is average of L & R; stim-motor is combines 'stim-lock' w 'motor' with arbituary weighting.
menu = {'(choose)','hier.','size','stim-lock','corr','motor','L motor','R motor','L+R motor',...
    'multi-motor','multi-motor w/ stim-avr','multi-stim w/ stim-avr'};
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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=3; % (updated when loading fish)
menu = {'(choose)',''};
hstimreg = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_getstimreg_Callback);

i=i+n;
n=2; % stimulus regressors, type range of stim-reg ID (e.g. '1-3,5',but can't use 'end')
uicontrol('Parent',tab{i_tab},'Style','text','String','stim combo:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_getstimregcombo_Callback);

i=i+n; % motor regressors from behavior, not yet convolved/adjusted for time lag
n=2; % go to 'GetMotorRegressor.m' to add/update
uicontrol('Parent',tab{i_tab},'Style','text','String','Motor reg.:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Plot corr. hist',...
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
    'Callback',{@pushbutton_allCentroidRegression_Callback});

i=i+n; % this is a remnant button from a failed experiment, idea was to iterate the regression process
n=2; % until the cluster converges, but most of the time it doesn't...
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','iter.reg','Enable','off',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_IterCentroidRegression_Callback});

%% UI row 3: t-tests
i_row = 3;
i = 1;n = 0;

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','Choose stim pair:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_tteststimrange_Callback); % e.g. '1-3,5', but can't use 'end'

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','t-test thres:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit','String',num2str(thres_ttest),...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_ttestthres_Callback);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','t-test',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_ttest_Callback);

% i=i+n;
% n=2;
% s = 'Choose single motor-regressor';
% uicontrol('Parent',tab{i_tab},'Style','text','String','Motor reg:',...
%     'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');
%
% i=i+n;
% n=2; % (unlike stim regressors, names hardcoded, not importet from regressor...)
% menu = {'(choose)','left swims','right swims','forward swims','raw left','raw right','raw average'};
% uicontrol('Parent',tab{i_tab},'Style','popupmenu','String',menu,'Value',1,...
%     'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
%     'Callback',@popup_getmotorreg_Callback);

%% UI ----- tab four ----- (Clustering etc.)
i_tab = 4;

%% UI row 1: k-means
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2; % k-means clustering
uicontrol('Parent',tab{i_tab},'Style','text','String','k-means, k =',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_kmeans_Callback);

i=i+n;
n=4; % anatomy is added to the fluo trace as new dimensions, and (arbituarily) weighted strongly
uicontrol('Parent',tab{i_tab},'Style','text','String','k-means with anatomy, k =',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_kmeans2_Callback);

i=i+n; % trying to use Silhouette to evaluate cluster quality, find peak to determine optimal k,
n=3; % then display results with that k. But have not set k-means to replicate (speed concern), can be very noisy
uicontrol('Parent',tab{i_tab},'Style','text','String','Find best k in range:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_hierpartn_Callback);

i=i+n; % partitioning based on hier. clustering
n=3; % ...or set correlation value threshold
uicontrol('Parent',tab{i_tab},'Style','text','String','Hier.cut, corr thres:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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

%% UI ----- tab five ----- (Saved Clusters)
i_tab = 5;

%% UI row 1: Cluster-Group (groups of cluster, one level above 'Clusters')
i_row = 1;
i = 1;n = 0;

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','ClusterGroup:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=2;
% menu = MakeNumberedMenu(VAR(i_fish).ClusGroupName);
hclusgroupmenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_clusgroupmenu_Callback);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','Edit name:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=2;
hclusgroupname = uicontrol('Parent',tab{i_tab},'Style','edit','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_clusgroupname_Callback);

i=i+n; % just adds a new number to the ClusterGroup-number menu,
n=3; % and saves current view as the first cluster in the new Clustergroup
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','New ClusterGroup',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_newclusgroup_Callback});

i=i+n;
n=3; % delete current ClusterGroup
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Delete ClusterGroup',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_delclusgroup_Callback});

%% UI row 2: Clusters
i_row = 2;
i = 1;n = 0;

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','Cluster:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=3;
% menu = MakeNumberedMenu({Cluster.name});
hclusmenu = uicontrol('Parent',tab{i_tab},'Style','popupmenu','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@popup_clusmenu_Callback);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','Edit name:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=3;
hclusname = uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_editclusname_Callback});
% set(hclusname,'String',Cluster(clusID).name);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Save cluster',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_saveclus_Callback);

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','New name:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','edit','String','(blank)',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_newclusname_Callback);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Create cluster',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@pushbutton_makeclus_Callback);

i=i+n;
n=2;
uicontrol('Parent',tab{i_tab},'Style','text','String','Set rank:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@edit_setrank_Callback});

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','text','String','Notes:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

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
n=2; % Combines the chosen clusters into one view
uicontrol('Parent',tab{i_tab},'Style','text','String','Union(cluster):',... % (eg 1,3-5)
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1;
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_clusUnion_Callback);

%% UI ----- tab six ----- (Atlas)
i_tab = 6;

%% UI row 1: find masks
i_row = 1;
i = 1;n = 0;

% row-height exception! listboxes are tall
i=i+n+5;
n=2;
s = 'plot histogram of all masks, also printing thresholded mask-names';
uicontrol('Parent',tab{i_tab},'Style','pushbutton','String','Find masks',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',{@pushbutton_findmasks_Callback},'TooltipString',s);

i=i+n;
n=2; % if checked, show thresholded masks on right-side plot
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show masks',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',0,...
    'Callback',@checkbox_showmasks_Callback);

%% UI row 2:
i_row = 2;
i = 1;n = 0;

i=i+n+5;
n=2; % display selected mask(s)
uicontrol('Parent',tab{i_tab},'Style','text','String','Choose mask ID:',...
    'Position',[grid(i) yrow(i_row)-dTextHt bwidth*n rheight],'HorizontalAlignment','right');

i=i+n;
n=1; % manually input (can choose from histogram recommendation)
uicontrol('Parent',tab{i_tab},'Style','edit',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],...
    'Callback',@edit_chooseMskID_Callback);

i=i+n;
n=3; % if checked, show thresholded masks on right-side plot
uicontrol('Parent',tab{i_tab},'Style','checkbox','String','Show outline only',...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight],'Value',1,...
    'Callback',@checkbox_showmskoutline_Callback);

%% UI row 3: find masks
i_row = 3;
i = 1;n = 0;

% row-height exception! listboxes are tall
i=i+n;
n=5; % manually input (can choose from histogram recommendation)
hmasklistbox = uicontrol('Parent',tab{i_tab},'Style','listbox','Max',10,'Min',0,...
    'Position',[grid(i) yrow(i_row) bwidth*n rheight*5],...
    'Callback',@listbox_chooseMskID_Callback);

%% Load figure

% UpdateClusID(hfig,clusID);

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
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);
disp('saved to workspace ''VAR''');
end

function pushbutton_savemat_Callback(hObject,~)
disp('saving...');
hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');
data_masterdir = getappdata(hfig,'data_masterdir');

% copy of VAR files will be saved into this subfolder:
arcmatfolder = fullfile(data_masterdir, 'arc mat');
if ~exist(arcmatfolder, 'dir')
    mkdir(arcmatfolder);
end

global VAR;
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);
timestamp  = datestr(now,'mmddyy_HHMM');
matname = [timestamp '.mat'];
save(fullfile(arcmatfolder,matname),'VAR','-v6');

% also save the current VAR file
save(fullfile(data_masterdir,'VAR_current.mat'),'VAR','-v6');
disp('saved both to workspace and .mat');
end

function pushbutton_writeZstack_Callback(hObject,~)
disp('preparing z-stack...');
isMarkAllCells = 1;

hfig = getParentFigure(hObject);
% isfullfish = getappdata(hfig,'isfullfish');
% if ~isfullfish,
%     errordlg('Load full fish first!');
%     return;
% end

cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
CellXYZ = getappdata(hfig,'CellXYZ');
absIX = getappdata(hfig,'absIX');
anat_stack = getappdata(hfig,'anat_stack');

timestamp = datestr(now,'mmddyy_HHMMSS');
tiffName = ['stack_' timestamp '.tif'];

U = unique(gIX);
numK = length(U);

% left half: stack with cells marked;
% right half: original anatomy, or mark all cells

anat_stack2 = zeros(size(anat_stack,1), size(anat_stack,2)*2, size(anat_stack,3) ,3);
nPlanes = size(anat_stack,3);
dimv_yxz = size(anat_stack);
stacklen = numel(anat_stack);

circle = makeDisk2(10,21);
[r, v] = find(circle);
r = r-11; v = v-11;
circle_inds  = r*dimv_yxz(1)+v;
cmap = hsv(numK);
weight = 0.3;

for i = 1:nPlanes,
    anat_stack2(:,:,i,:) = repmat(imNormalize99(anat_stack(:,:,i)),[1 2 1 3]);
end

cIX_abs = absIX(cIX);
for j = 1:length(cIX),
    cinds = (CellXYZ(cIX_abs(j),2)-1)*dimv_yxz(1)+CellXYZ(cIX_abs(j),1);
    labelinds = find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yxz(1)*dimv_yxz(2));
    zinds = dimv_yxz(1)*dimv_yxz(2)*2*(CellXYZ(cIX_abs(j),3)-1);
    ix = find(U==gIX(j));
    ixs = cinds+circle_inds(labelinds)+zinds;
    anat_stack2(ixs) = cmap(ix,1)*weight + anat_stack2(ixs)*(1-weight);
    ixs = cinds+circle_inds(labelinds)+zinds+stacklen*2;
    anat_stack2(ixs) = cmap(ix,2)*weight + anat_stack2(ixs)*(1-weight);
    ixs = cinds+circle_inds(labelinds)+zinds+stacklen*4;
    anat_stack2(ixs) = cmap(ix,3)*weight + anat_stack2(ixs)*(1-weight);
end
if isMarkAllCells,
    clr = [0,1,0];
    for j = 1:size(CellXYZ,1),
        shift = dimv_yxz(1)*dimv_yxz(2);
        cinds = (CellXYZ(j,2)-1)*dimv_yxz(1)+CellXYZ(j,1);
        labelinds = find((cinds+circle_inds)>0 & (cinds+circle_inds)<=dimv_yxz(1)*dimv_yxz(2));
        zinds = dimv_yxz(1)*dimv_yxz(2)*2*(CellXYZ(j,3)-1);
        ixs = cinds+circle_inds(labelinds)+zinds + shift;
        anat_stack2(ixs) = clr(1)*weight + anat_stack2(ixs)*(1-weight);
        ixs = cinds+circle_inds(labelinds)+zinds+stacklen*2 + shift;
        anat_stack2(ixs) = clr(2)*weight + anat_stack2(ixs)*(1-weight);
        ixs = cinds+circle_inds(labelinds)+zinds+stacklen*4 + shift;
        anat_stack2(ixs) = clr(3)*weight + anat_stack2(ixs)*(1-weight);
    end
end
h = figure;
for i_plane = 1:nPlanes,
    im = squeeze(anat_stack2(:,:,i_plane,:));
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
% isfullfish = getappdata(hfig,'isfullfish');
% if isfullfish,
    disp('Rendering...')
    DrawTiledPics(hfig);
% else
%     errordlg('Load full fish first!');
% end
end

function pushbutton_tileClusters_Callback(hObject,~)
hfig = getParentFigure(hObject);
% isfullfish = getappdata(hfig,'isfullfish');
% if isfullfish,
    disp('Rendering...')
    DrawTiledPics_clus(hfig);
% else
%     errordlg('Load full fish first!');
% end
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

function pushbutton_exporttoworkspace_Callback(hObject,~)
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');
M_0_fluo = GetTimeIndexedData(hfig,'isAllCells');
assignin('base', 'M', M);
assignin('base', 'M_0_fluo', M_0_fluo);
assignin('base', 'behavior', getappdata(hfig,'behavior'));
assignin('base', 'stim', getappdata(hfig,'stim'));

assignin('base', 'cIX', getappdata(hfig,'cIX'));
assignin('base', 'gIX', getappdata(hfig,'gIX'));
assignin('base', 'tIX', getappdata(hfig,'tIX'));
assignin('base', 'numK', getappdata(hfig,'numK'));
end

function pushbutton_loadVARfromworkspace_Callback(hObject,~)
% (code copied from UpdateFishData and UpdateFishDisplay)
disp('import VAR from workspace');
hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');

global VAR;
ClusGroup = VAR(i_fish).ClusGroup; % group of 'Cluster', i.e. 2nd level collections
setappdata(hfig,'ClusGroup',ClusGroup);

clusgroupID = 1;
setappdata(hfig,'clusgroupID',clusgroupID);

Cluster = ClusGroup{clusgroupID};% collection of clusters
setappdata(hfig,'Cluster',Cluster);

% update display
clusgroupID = 1;
new_clusgroupID = 1;
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID); % to display new menu
end

function pushbutton_loadCurrentClustersfromworkspace_Callback(hObject,~)
disp('import current clusters from workspace');
hfig = getParentFigure(hObject);

cIX = evalin('base','cIX');
gIX = evalin('base','gIX');
numK = max(gIX);
UpdateIndices(hfig,cIX,gIX,numK);
RefreshFigure(hfig);
end

function pushbutton_popupplot_Callback(hObject,~)
% very similar as function RefreshFigure(hfig)
isPopout = 1; % no down-sampling in plots
% load
hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');
isCentroid = getappdata(hfig,'isCentroid');
isRefAnat = getappdata(hfig,'isRefAnat');
isPlotLines = getappdata(hfig,'isPlotLines');
isPlotBehavior = getappdata(hfig,'isPlotBehavior');
isPlotAnatomyOnly = getappdata(hfig,'isPlotAnatomyOnly');

if ~isPlotAnatomyOnly,
    figure('Position',[50,100,1600,800],'color',[1 1 1],...
        'Name',['Fish#' num2str(i_fish)]);
    h1 = axes('Position',[0.05, 0.03, 0.53, 0.94]); % left ~subplot
    h2 = axes('Position',[0.61, 0.03, 0.35, 0.94]); % right ~subplot
    
    % left subplot
    axes(h1);
    DrawTimeSeries(hfig,h1,isPopout,isCentroid,isPlotLines,isPlotBehavior);
    
    % right subplot
    axes(h2);
    DrawCellsOnAnatProj(hfig,isRefAnat,isPopout);
    
else
    figure('Position',[600,300,600,900],'color',[1 1 1],...
        'Name',['Fish#' num2str(i_fish)]);
    axes('Position',[0.03, 0.03, 0.94, 0.94]); % right ~subplot
    % right subplot
    DrawCellsOnAnatProj(hfig,isRefAnat,isPopout);
end
end

function checkbox_isPlotLines_Callback(hObject,~)
hfig = getParentFigure(hObject);
value = get(hObject,'Value');
setappdata(hfig,'isPlotLines',value);
end

function checkbox_isPlotBehavior_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotBehavior',get(hObject,'Value'));
end

function checkbox_isPlotAnatomyOnly_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotAnatomyOnly',get(hObject,'Value'));
end

%% row 2: Load

function popup_loadfullfishmenu_Callback(hObject,~)
i_fish = get(hObject,'Value')-1;
hfig = getParentFigure(hObject);
if i_fish>0,
    tic    
    watchon; drawnow;
    LoadFullFish(hfig,i_fish);
    UpdateFishData(hfig,i_fish);
    UpdateFishDisplay(hfig);
    toc
    watchoff;
end
end

function pushbutton_choosefilestoload_Callback(hObject,~)
hfig = getParentFigure(hObject);
data_masterdir = getappdata(hfig,'data_masterdir');

% get fish number
prompt={'Enter fish number:'};
answer = inputdlg(prompt);
if ~isempty(answer),
    new_i_fish = str2num(answer{:});
    
    [FileName1,PathName] = uigetfile('*.h5','Select the hdf5(.h5) file for TimeSeries data',data_masterdir);
    hdf5_dir = fullfile(PathName,FileName1);
    [FileName2,PathName] = uigetfile('*.mat','Select the .mat file for other data',PathName);
    mat_dir = fullfile(PathName,FileName2);
    
    if isequal(FileName1,0) || isequal(FileName2,0),
        disp('User selected Cancel')
    else
        % display fish-number in hloadfish
        global hloadfish; %#ok<TLEV>
        set(hloadfish,'Value',new_i_fish+1);
        
        LoadFullFish(hfig,new_i_fish,hdf5_dir,mat_dir);
    end
end
end

function LoadFullFish(hfig,i_fish,hdf5_dir,mat_dir)
disp(['loading fish #' num2str(i_fish) '...']);

M_fish_set = getappdata(hfig,'M_fish_set');
fishset = M_fish_set(i_fish);
setappdata(hfig,'fishset',fishset);

%% load data from file

disp(['load fish ' num2str(i_fish) '...']);

if ~exist('hdf5_dir','var') || ~exist('mat_dir','var'),
    data_masterdir = getappdata(hfig,'data_masterdir');
    data_dir = fullfile(data_masterdir,['subject_' num2str(i_fish)]);
    hdf5_dir = fullfile(data_dir,'TimeSeries.h5');
    mat_dir = fullfile(data_dir,'data_full.mat');
end

tic

% load 'data'
load(mat_dir,'data'); % struct with many fields
names = fieldnames(data); % cell of strings
for i = 1:length(names),
    setappdata(hfig,names{i},eval(['data.',names{i}]));
end

% load time series (hdf5 file)
CellResp = h5read(hdf5_dir,'/CellResp');
CellRespZ = h5read(hdf5_dir,'/CellRespZ');
CellRespAvr = h5read(hdf5_dir,'/CellRespAvr');
CellRespAvrZ = h5read(hdf5_dir,'/CellRespAvrZ');
absIX = h5read(hdf5_dir,'/absIX');

setappdata(hfig,'CellResp',CellResp);
setappdata(hfig,'CellRespZ',CellRespZ);
setappdata(hfig,'CellRespAvr',CellRespAvr);
setappdata(hfig,'CellRespAvrZ',CellRespAvrZ);
setappdata(hfig,'absIX',absIX);

toc
% setappdata(hfig,'isfullfish',1);
end

function UpdateFishData(hfig,new_i_fish)
%% set Cluster display
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
ClusGroup = VAR(new_i_fish).ClusGroup; % group of 'Cluster', i.e. 2nd level collections
setappdata(hfig,'ClusGroup',ClusGroup);

clusgroupID = 1;
setappdata(hfig,'clusgroupID',clusgroupID);

Cluster = ClusGroup{clusgroupID};% collection of clusters
setappdata(hfig,'Cluster',Cluster);

%% update GUI timelists selection, if applicable
global hstimrangemenu;
if ~isempty(hstimrangemenu), % before GUI initialization,
    % 'global' remnant from previous run: hstimrangemenu is not empty but is invalid object
    if isvalid(hstimrangemenu),
        timelists_names = getappdata(hfig,'timelists_names');
        s = cell(size(timelists_names));
        for i = 1:length(timelists_names),
            s{i} = [num2str(i),': ',timelists_names{i}];
        end
        set(hstimrangemenu,'String',s);
        set(hstimrangemenu,'Value',1);
    end
end

%% set current timelists
periods = getappdata(hfig,'periods');
setappdata(hfig,'stimrange',1:length(periods));
UpdateTimeIndex(hfig,'isSkipcIX'); % doesn't include Refresh

end

function UpdateFishDisplay(hfig) % loading steps that are not performed at very first initialization
stim = getappdata(hfig,'stim');
fishset = getappdata(hfig,'fishset');

[~, names] = GetStimRegressor(stim,fishset);

s = cell(size(names));
for i = 1:length(names),
    s{i} = [num2str(i),': ',names{i}];
end

global hstimreg;
set(hstimreg,'String',['(choose)',s]);

% update display
clusgroupID = 1;
new_clusgroupID = 1;
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID); % to display new menu
end

function popup_stimrangemenu_Callback(hObject,~)
ID = get(hObject,'Value');
hfig = getParentFigure(hObject);
setappdata(hfig,'stimrange',ID);
hfig = getParentFigure(hObject);

UpdateTimeIndex(hfig);
RefreshFigure(hfig);
end

%% row 3: Display

function checkbox_isStimAvr_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isAvr',get(hObject,'Value'));
UpdateTimeIndex(hfig);
RefreshFigure(hfig);
end

function checkbox_isRawtime_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isRawtime',get(hObject,'Value'));
disp('updating display...');
UpdateTimeIndex(hfig);
RefreshFigure(hfig);
disp('complete');
end

function checkbox_isZscore_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isZscore',get(hObject,'Value'));
disp('updating display...');
UpdateTimeIndex(hfig);
RefreshFigure(hfig);
disp('complete');
end

function edit_stimrange_Callback(hObject,~)
hfig = getParentFigure(hObject);
periods = getappdata(hfig,'periods');

% get/format range
str = get(hObject,'String');
if ~isempty(str),
    str = strrep(str,'end',num2str(length(periods)));
    stimrange = ParseRange(str);
    setappdata(hfig,'stimrange',stimrange);
    UpdateTimeIndex(hfig);
    RefreshFigure(hfig);
end
end

function checkbox_showcentroids_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isCentroid',get(hObject,'Value'));
RefreshFigure(hfig);
end

function checkbox_showrefanat_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isRefAnat',get(hObject,'Value'));
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
    % save
    setappdata(hfig,'bCache',bC);
    setappdata(hfig,'fCache',fC);
    setappdata(hfig,'cIX',cIX);
    setappdata(hfig,'gIX',gIX);
    setappdata(hfig,'numK',numK);
    % handle rankID: >=2 means write numbers as text next to colorbar
    setappdata(hfig,'rankID',0);
    setappdata(hfig,'isWeighAlpha',0);
    
    M = GetTimeIndexedData(hfig);
    setappdata(hfig,'M',M);
    
    % FindCentroid reset:
    setappdata(hfig,'Centroids',[]);
    setappdata(hfig,'D_ctrd',[]);
    
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
    % save
    setappdata(hfig,'bCache',bC);
    setappdata(hfig,'fCache',fC);
    setappdata(hfig,'cIX',cIX);
    setappdata(hfig,'gIX',gIX);
    setappdata(hfig,'numK',numK);
    % handle rankID: >=2 means write numbers as text next to colorbar
    setappdata(hfig,'rankID',0);
    setappdata(hfig,'isWeighAlpha',0);
    
    M = GetTimeIndexedData(hfig);
    setappdata(hfig,'M',M);
    
    % FindCentroid reset:
    setappdata(hfig,'Centroids',[]);
    setappdata(hfig,'D_ctrd',[]);
    
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
        range = [range,m(i)];
    else % have '-'sign,
        range = [range,m(i-1)+1:-m(i)];
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

[gIX, numU] = SqueezeGroupIX(gIX);
switch rankID,
    case 1,
        disp('hier. (default)');
        M = getappdata(hfig,'M');
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
        [gIX,rankscore] = RankByStimLock_Direct(hfig,gIX,numU);
    case 4,
        disp('corr');
        [~,D] = FindCentroid(hfig);
        [gIX,rankscore] = SortH(D,gIX,numU);
        rankscore = 1-rankscore;
    case 5,
        disp('motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,numU,1);
    case 6,
        disp('L motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,numU,2);
    case 7,
        disp('R motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,numU,3);
    case 8,
        disp('L+R motor');
        [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,numU,4);
    case 9,
        %         disp('stim-motor');
        %         [~,rankscore1] = RankByStimLock_Direct(hfig,cIX,gIX,numU);
        %         [~,rankscore2] = RankByMotorStim_Direct(hfig,gIX,numU,1);
        %         rankscore = rankscore1 - rankscore2;
        disp('multi-motor');
        [gIX,rankscore] = RankByMultiRegression_Direct(hfig,gIX,numU,1);
    case 10,
        disp('multi-motor w/ stim-avr');
        [gIX,rankscore] = RankByMultiRegression_Direct(hfig,gIX,numU,3);
    case 11,
        disp('multi-stim w/ stim-avr');
        [gIX,rankscore] = RankByMultiRegression_Direct(hfig,gIX,numU,4);
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

function [gIX,rankscore] = RankByStimLock_Direct(hfig,gIX,numU)
periods = getappdata(hfig,'periods');
fishset = getappdata(hfig,'fishset');

isAvr = getappdata(hfig,'isAvr');
isRawtime = getappdata(hfig,'isRawtime');
if isAvr == 1 || isRawtime == 1,
    setappdata(hfig,'isAvr',0);
    setappdata(hfig,'isRawtime',0);
    global h_isavr h_israwtime; %#ok<TLEV>
    h_isavr.Value = 0;
    h_israwtime.Value = 0;
    UpdateTimeIndex(hfig);
end

C = FindCentroid(hfig);

if fishset == 1,
    period = periods;
    C_3D_0 = reshape(C,size(C,1),period,[]);
    C_3D = zscore(C_3D_0,0,2);
    
    H = nanmean(nanstd(C_3D,0,3),2);
else
    stimset = getappdata(hfig,'stimset');
    stimrange = getappdata(hfig,'stimrange');
    periods = getappdata(hfig,'periods');
    timelists = getappdata(hfig,'timelists');
    
    range = [];
    H_raw = [];
    for i = 1:length(stimrange),
        i_stim = stimrange(i);
        if sum(stimset(i_stim).nReps)>3,
            range = [range,i_stim];
            offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
            tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
            period = periods(i_stim);
            C_3D_0 = reshape(C(:,tIX_),size(C,1),period,[]);
            C_3D = zscore(C_3D_0,0,2);
            H_raw = horzcat(H_raw,nanmean(nanstd(C_3D,0,3),2));
        end
    end
    if isempty(range),
        errordlg('chosen stimulus range not suitable for stim-lock analysis');
        rankscore = 1:numU;
        return;
    end
    H = zeros(size(H_raw,1),1);
    for i = 1:size(H_raw,1),
        H(i) = min(H_raw(i,:));
    end
    
    %     % instead of algebraic average along 2nd dimension, use
    %     % inverse of geometric average... large value~low variation. geometric
    %     % mean biases towards large values, i.e. good stim-lock of any stimulus
    %     % is emphasized.
    %     H = zeros(size(H_raw,1),1);
    %     for i = 1:size(H_raw,1),
    %         temp = sum((1./H_raw(i,:)).^2);
    %         H(i) = 1./sqrt(temp);
    %     end
    
end

[gIX,rankscore] = SortH(H,gIX,numU);
end

function [gIX,rankscore] = RankByMotorStim_Direct(hfig,gIX,numU,option)
C = FindCentroid(hfig);
behavior = getappdata(hfig,'behavior');
regressors = GetMotorRegressor(behavior);
% % behavior(1,:);   %weighted: right turns
% % behavior(2,:);   %weighted: left turns
% % behavior(3,:);   %weighted: forward swims
% % behavior(4,:);  %analog: right channel
% % behavior(5,:);  %analog: left channel
% % behavior(4,:)+behavior(5,:);   %analog: average

H = zeros(numU,1);
% shift = zeros(numU,1);
% a = zeros(1,3);
% I = zeros(1,3);
for i = 1:numU,
    switch option,
        %         case 1, % 'motor'
        %             for j = 1:3,
        %                 [a(j),I(j)] = max(abs(xcorr(C(i,:),reg(j,:),'coeff')));
        %             end
        %             [H(i),jmax] = max(a);
        %             shift(i) = I(jmax) - length(behavior);
        %         case 2, % 'L motor'
        %             [H(i),I] = max(xcorr(C(i,:),reg(1,:),'coeff'));
        %             shift(i) = I - length(behavior);
        %         case 3, % 'R motor'
        %             [H(i),I] = max(xcorr(C(i,:),reg(2,:),'coeff'));
        %             shift(i) = I - length(behavior);
        %         case 4, % 'L+R motor'
        %             [H(i),I] = max(xcorr(C(i,:),reg(3,:),'coeff'));
        %             shift(i) = I - length(behavior);
        case 1, % 'motor'
            R = zeros(1,4);
            for j = 1:3,
                regressor = regressors(j).im;
                R(j) = corr(regressor',C(i,:)');
            end
            regressor = 0.5*(regressors(4).im+regressors(5).im);
            H(4) = corr(regressor',C(i,:)');
            H(i) = max(R);
        case 2, % 'L motor'
            regressor = regressors(1).im;
            H(i) = corr(regressor',C(i,:)');
        case 3, % 'R motor'
            regressor = regressors(2).im;
            H(i) = corr(regressor',C(i,:)');
        case 4, % 'L+R motor'
            regressor = 0.5*(regressors(4).im+regressors(5).im);
            H(i) = corr(regressor',C(i,:)');
    end
end
[gIX,rankscore] = SortH(H,gIX,numU,'descend');
% assignin('base', 'shift', shift);
end

function [gIX,rankscore] = RankByMultiRegression_Direct(hfig,gIX,numU,option)
%% get cluster centroids (means) from GUI current selection
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
C = FindCentroid(hfig);
nClus = size(C,1);

%% get stimulus regressor
switch option
    case 1; % use all pre-defined stim-regressors
        regressors = GetStimRegressor(stim,fishset);
        M_regressor = zeros(length(regressors),length(regressors(1).im));
        for i = 1:length(regressors),
            M_regressor(i,:) = regressors(i).im;
        end
        regressor_s = M_regressor;
        
    case 2; % get stim regressors and find best match
        stimregset = 1:8;
        regressors = GetStimRegressor(stim,fishset);
        M_regressor = zeros(length(stimregset),length(regressors(1).im));
        for i = 1:length(stimregset),
            M_regressor(i,:) = regressors(stimregset(i)).im;
        end
        
        R = corr(M_regressor',C'); % row of R: each regressor
        [~,IX] = max(R,[],1);
        regressor_s_allclus = M_regressor(IX,:);
        
    case {3,4}; % alternative: using 'stim-lock' means
        periods = getappdata(hfig,'periods');
        
        if fishset == 1,
            period = periods;
            C_3D_0 = reshape(C,size(C,1),period,[]);
            C_period = mean(C_3D_0,3);
            nPeriods = round(size(C,2)/period);
            C_mean = repmat(C_period,1,nPeriods);
        else
            stimset = getappdata(hfig,'stimset');
            stimrange = getappdata(hfig,'stimrange');
            periods = getappdata(hfig,'periods');
            timelists = getappdata(hfig,'timelists');
            
            C_mean = [];
            for i = 1:length(stimrange),
                i_stim = stimrange(i);
                if sum(stimset(i_stim).nReps)>3,
                    offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
                    tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
                    period = periods(i_stim);
                    C_3D_0 = reshape(C(:,tIX_),size(C,1),period,[]);
                    C_period = mean(C_3D_0,3);
                    nPeriods = length(tIX_)/period;
                    C_mean = horzcat(C_mean,repmat(C_period,1,nPeriods));
                    %             C_3D = zscore(C_3D_0,0,2);
                    %             H_raw = horzcat(H_raw,nanmean(nanstd(C_3D,0,3),2));
                else
                    offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
                    tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
                    C_mean = horzcat(C_mean,zeros(size(C,1),length(tIX_)));
                end
            end
        end
        
        regressor_s_allclus = C_mean;
end

%% get motor regressor
regressors = GetMotorRegressor(behavior);
motorregset = 1:3;
M_regressor = zeros(length(motorregset),length(regressors(1).im));
for i = 1:length(motorregset),
    M_regressor(i,:) = regressors(motorregset(i)).im;
end

% switch option
%     case 1;
regressor_m = M_regressor;
%     case {2,3};
%         R = corr(M_regressor',C'); % row of R: each regressor
%         [~,IX_m] = max(R,[],1);
%         regressor_m_allclus = M_regressor(IX_m,:);
% end

%% multi-regression
switch option
    case 1; % regression with all regs, stim-regs before motor regs
        regs = vertcat(regressor_s,regressor_m);
        orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?
        betas = C * orthonormal_basis; % to reconstitute: betas(i,:)*orthonormal_basis'
        % get ranking score: combined of motor coeffs
        H = sqrt(sum((betas(:,end-2:end)).^2,2));
    case {2,3}; % regression with stim-avr reg before motor regs
        betas = zeros(nClus,1+length(motorregset)); %size(regressor_s,1)+1);
        for i_clus = 1:nClus,
            regs = vertcat(regressor_s_allclus(i_clus,:),regressor_m);
            %             regs = vertcat(regressor_s_allclus(i_clus,:),regressor_m_allclus(i_clus,:));
            
            orthonormal_basis = Gram_Schmidt_Process(regs');
            
            %% check orthonormal_basis
            % figure;imagesc(orthonormal_basis)
            % norm(orthonormal_basis(:,1))
            % norm(orthonormal_basis(:,2))
            % dot(orthonormal_basis(:,1),orthonormal_basis(:,2))
            % % compare to raw regressors
            % figure;
            % subplot(311);hold on;
            % i = 1;
            % plot(regs(i,:)/norm(regs(i,:)),'r');plot(orthonormal_basis(:,1)','k:')
            % subplot(312);hold on;
            % i = 2;
            % plot(regs(i,:)/norm(regs(i,:)),'r');plot(orthonormal_basis(:,i)','k:')
            % subplot(313);hold on;
            % i = 3;
            % plot(regs(i,:)/norm(regs(i,:)),'r');plot(orthonormal_basis(:,i)','k:')
            
            %% Get multi-regression coefficients
            % FunctionalActivity = beta1*reg1 + beita2*reg2_orth
            betas(i_clus,:) = C(i_clus,:) * orthonormal_basis;
        end
        % get ranking score: combined motor coeffs
        H = sqrt(sum((betas(:,end-2:end)).^2,2));
    case 4; % regression with motor regs before stim-avr reg
        betas = zeros(nClus,1+length(motorregset)); %size(regressor_s,1)+1);
        for i_clus = 1:nClus,
            regs = vertcat(regressor_m,regressor_s_allclus(i_clus,:));
            orthonormal_basis = Gram_Schmidt_Process(regs');
            betas(i_clus,:) = C(i_clus,:) * orthonormal_basis;
        end
        % get ranking score: coeffs of (the single) stim-avr reg
        H = betas(:,end);
end

%% rank and plot
[gIX,rankscore] = SortH(H,gIX,numU,'descend');

figure;
h1 = subplot(2,1,1);
imagesc(orthonormal_basis);
title('All orthonormal bases');
xlabel('orthonormal basis #');
% make x-lables
nBases = size(orthonormal_basis,2);
s = [];
for i_basis = 1:nBases-3,
    s = [s,{['stim.' num2str(i_basis)]}];
end
s = [s,{'motor.R','motor.L','motor.F'}];
h1.XTick = 1:nBases;
h1.TickLength = [0,0];
h1.XTickLabel = s;
h1.XTickLabelRotation = 45;
ylabel('time ~ frames');

h2 = subplot(2,1,2);
imagesc(betas)
title('multi-regression coefficients');
xlabel('orthonormal basis #');
h2.XTick = 1:nBases;
h2.TickLength = [0,0];
h2.XTickLabel = s;
h2.XTickLabelRotation = 45;
ylabel('cluster ID');

end

function [gIX,B] = SortH(H,gIX,numU,descend) %#ok<INUSD> % new gIX is sorted based on H, size(H)=[numU,1];
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
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
absIX = getappdata(hfig,'absIX');
isRefAnat = getappdata(hfig,'isRefAnat');
if ~isRefAnat,
    CellXYZ = getappdata(hfig,'CellXYZ');
    anat_yx = getappdata(hfig,'anat_yx');
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    k_zres = 20;
else
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    anat_yx = getappdata(hfig,'anat_yx_norm');
    anat_yz = getappdata(hfig,'anat_yz_norm');
    anat_zx = getappdata(hfig,'anat_zx_norm');
    k_zres = 2.5;
end
dimv_yx = size(anat_yx);
dimv_zx = size(anat_zx);

h_poly_yx = impoly;
wait(h_poly_yx); % double click to finalize position!
% update finalized polygon in bright color
setColor(h_poly_yx,[0 1 1]);

A = sub2ind(dimv_yx(1:2),CellXYZ(absIX,1),CellXYZ(absIX,2));

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
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
absIX = getappdata(hfig,'absIX');
isRefAnat = getappdata(hfig,'isRefAnat');
if ~isRefAnat,
    CellXYZ = getappdata(hfig,'CellXYZ');
    anat_yx = getappdata(hfig,'anat_yx');
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    k_zres = 20;
else
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    anat_yx = getappdata(hfig,'anat_yx_norm');
    anat_yz = getappdata(hfig,'anat_yz_norm');
    anat_zx = getappdata(hfig,'anat_zx_norm');
    k_zres = 2.5;
end
dimv_yx = size(anat_yx);
dimv_yz = size(anat_yz);
dimv_zx = size(anat_zx);

h_poly_z = impoly;
wait(h_poly_z); % double click to finalize position!
% update finalized polygon in bright color
setColor(h_poly_z,[0 1 1]);

A = sub2ind(dimv_yz(1:2),CellXYZ(absIX,1),CellXYZ(absIX,3));

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
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
absIX = getappdata(hfig,'absIX');
isRefAnat = getappdata(hfig,'isRefAnat');
if ~isRefAnat,
    CellXYZ = getappdata(hfig,'CellXYZ');
    anat_yx = getappdata(hfig,'anat_yx');
    anat_yz = getappdata(hfig,'anat_yz');
    anat_zx = getappdata(hfig,'anat_zx');
    k_zres = 20;
else
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    anat_yx = getappdata(hfig,'anat_yx_norm');
    anat_yz = getappdata(hfig,'anat_yz_norm');
    anat_zx = getappdata(hfig,'anat_zx_norm');
    k_zres = 2.5;
end
dimv_zx = size(anat_zx);

h_poly_z = impoly;
wait(h_poly_z); % double click to finalize position!
% update finalized polygon in bright color
setColor(h_poly_z,[0 1 1]);

A = sub2ind(dimv_zx(1:2),CellXYZ(absIX,3),CellXYZ(absIX,2));

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
CellXYZ = getappdata(hfig,'CellXYZ');
absIX = getappdata(hfig,'absIX');
cIX_abs = absIX(cIX);

tri = delaunayn(CellXYZ(cIX_abs,:)); % Generate delaunay triangulization
t = tsearchn(CellXYZ(cIX_abs,:), tri, CellXYZ(absIX,:)); % Determine which triangle point is within
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

function edit_getstimregcombo_Callback(hObject,~)
hfig = getParentFigure(hObject);
% get/format range
str = get(hObject,'String');
if ~isempty(str),
    %     str = strrep(str,'end',num2str(nMasks));
    range = ParseRange(str);
    
    setappdata(hfig,'regchoice',{1,range});
    % highlight the choice (yellow)
    global hstimreg hmotorreg hcentroidreg; %#ok<TLEV>
    set(hstimreg,'BackgroundColor',[1,1,1]);
    %         set(hstimreg,'BackgroundColor',[1,1,0.8]); % yellow
    set(hmotorreg,'BackgroundColor',[1,1,1]);
    set(hcentroidreg,'BackgroundColor',[1,1,1]);
end
end

function PlotRegWithStimMotor(hfig)
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
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

temp = behavior;
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
    %         plot([0.5,length(behavior)+0.5],[y,y],'w','Linewidth',0.5);
    %     end
    %     % labels
    %     names = {'Left','Right','Forward','Raw L','Raw R'};
    %     x = -s2*0.05;
    %     for i = 1:5,
    %         y = i;
    %         text(x,y,names{i},'Fontsize',7);
    %     end
    
else % only plot top 3 lines
    fc = vertcat(temp(1,:),temp(3,:),temp(2,:)); %#ok<UNRCH>
    imagesc(fc);colormap gray
    set(gca,'YTick',[],'XTick',[]);
    set(gcf,'color',[1 1 1]);
    set(gca, 'box', 'off')
    hold on;axis ij;
    
    % plot division lines
    for i = 0:2,
        y = i+0.5;
        plot([0.5,length(behavior)+0.5],[y,y],'w','Linewidth',0.5);
    end
    % labels
    names = {'Left','Forward','Right'};
    x = -s2*0.07;
    for i = 1:3,
        y = i;
        text(x,y,names{i},'Fontsize',10);
    end
    
end

imagesc(behavior);colormap hot; axis off; title('motor');
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
    regressors = GetStimRegressor(stim,fishset);
    if length(regchoice{2})>1,
        regressor = zeros(length(regchoice{2}),length(regressors(1).im));
        for i = 1:length(regchoice{2}),
            regressor(i,:) = regressors(regchoice{2}(i)).im;
        end
    else
        regressor = regressors(regchoice{2}).im;
    end
    
elseif regchoice{1}==2, % motor Regressor
    behavior = getappdata(hfig,'behavior');
    regressors = GetMotorRegressor(behavior);
    regressor = regressors(regchoice{2}).im;
    
else % regchoice{1}==3, from Centroid
    ctrdID = regchoice{2};
    isflipstim = getappdata(hfig,'isflipstim');
    gIX = getappdata(hfig,'gIX');
    
    i = find(unique(gIX)==ctrdID);
    if isempty(i),
        disp('input is empty!');beep;
        regressor = [];
        return;
    end
    C = FindCentroid(hfig);
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

function pushbutton_thres_regression_Callback(hObject,~)
disp('regression...');
hfig = getParentFigure(hObject);
thres_reg = getappdata(hfig,'thres_reg');
regressor = GetRegressor(hfig);
if isempty(regressor),
    return;
end

isCentroid = getappdata(hfig,'isCentroid');

tic
isPlotCorrHist = getappdata(hfig,'isPlotCorrHist');
if size(regressor,1) == 1,
    [cIX,gIX,wIX] = Regression_Direct(hfig,thres_reg,regressor,isCentroid,isPlotCorrHist);
else
    CIX = []; GIX = []; WIX = [];
    for i = 1:size(regressor,1),
        [cIX_,gIX_,wIX_] = Regression_Direct(hfig,thres_reg,regressor(i,:),isCentroid,isPlotCorrHist);
        gIX_ = i*ones(size(gIX_));
        CIX = [CIX;cIX_]; GIX = [GIX;gIX_]; WIX = [WIX;wIX_];
    end
    [cIX,gIX,wIX] = SmartUnique_weighted(CIX,GIX,WIX);
end

if ~isempty(gIX),
    %     gIX = ceil((1:length(cIX))'/length(cIX)*min(20,length(cIX)));
    setappdata(hfig,'wIX',wIX);
    setappdata(hfig,'isWeighAlpha',1);
    UpdateIndices(hfig,cIX,gIX,length(unique(gIX)));
    %     UpdateIndices(hfig,cIX,gIX,40);
    RefreshFigure(hfig);
end
toc
beep

end

function [cIX,gIX,wIX] = SmartUnique_weighted(CIX,GIX,WIX)
[uCIX,ia,~] = unique(CIX);
uGIX = GIX(ia);
uWIX = WIX(ia);
if length(uCIX) == 1,
    counts = length(x);
else counts = hist(CIX,uCIX);
end
IX1 = find(counts==1);
CIX1 = uCIX(IX1); % single occurence
GIX1 = uGIX(IX1);
WIX1 = uWIX(IX1);
ia = find(~ismember(CIX,CIX1));
C = CIX(ia); % multiple copies/occurence, keep all copies
G = GIX(ia);
W = WIX(ia);
[~,I] = sort(W,'descend'); % rank by coeff, so unique ('stable') gets highest
C = C(I); % finish sorting
G = G(I);
[CIX2,ia,~] = unique(C,'stable'); % keep the order!
GIX2 = G(ia); % the chosen copy for those with multiple copies
WIX2 = W(ia);
cIX = vertcat(CIX1,CIX2);
gIX = vertcat(GIX1,GIX2);
wIX = vertcat(WIX1,WIX2);
disp('smart union (weighted) complete');
end

function [cIX_,gIX_,wIX_] = Regression_Direct(hfig,thres_reg,regressor,isCentroid,isPlotCorrHist) % gIX_ is just ones
if ~isCentroid,
    M_ = getappdata(hfig,'M_0');
else
    M_ = FindCentroid(hfig);
end
R = corr(regressor',M_');

if thres_reg>=0,
    cIX_ = (find(R>thres_reg))';
elseif thres_reg<0,
    cIX_ = (find(R<thres_reg))';
end
if isempty(cIX_),
    disp('result is empty!');beep;
    gIX_ = [];
    wIX_ = [];
    return;
end
% re-order cIX_ based on corr coeff
wIX = R(cIX_); % weight, here set to equal corr coeff


if isCentroid,
    % cIX_ is actually the selection of clusters
    clusterIX = cIX_;
    weightIX = wIX;
    
%     M = getappdata(hfig,'M');
    cIX = getappdata(hfig,'cIX');
    gIX = getappdata(hfig,'gIX');
    cIX_ = [];
    gIX_ = [];
    wIX_ = [];
    for i=1:length(clusterIX),
        IX = find(gIX == clusterIX(i));
        cIX_ = [cIX_; cIX(IX)];
        gIX_ = [gIX_; i*ones(length(IX),1)];
        wIX_ = [wIX_; weightIX(i)*ones(length(IX),1)];
    end
else
    [~,I] = sort(wIX,'descend');
    cIX_ = cIX_(I);
    gIX_ = ones(size(cIX_));
    wIX_ = wIX(I)';
end

% option: plot histogram
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

function checkbox_plotcorrhist_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isPlotCorrHist',get(hObject,'Value'));
end

function pushbutton_topnum_regression_Callback(hObject,~)
disp('regression...');
hfig = getParentFigure(hObject);
M_0 = getappdata(hfig,'M_0');
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

UpdateIndices(hfig,cIX,gIX,40);
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

function pushbutton_allCentroidRegression_Callback(hObject,~)
watchon; drawnow;
hfig = getParentFigure(hObject);
[cIX,gIX,numK] = AllCentroidRegression_direct(hfig);

UpdateIndices(hfig,cIX,gIX,numK);
RefreshFigure(hfig);
disp('all regression complete');
watchoff;
end

function pushbutton_IterCentroidRegression_Callback(hObject,~)
hfig = getParentFigure(hObject);
M_0 = getappdata(hfig,'M_0');
thres_reg = getappdata(hfig,'thres_reg');
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
C = FindCentroid(hfig);
regressor = C(i,:);
cIX_= Regression_Direct(hfig,thres_reg,regressor,0);
regressor_last = regressor;

if ~isempty(cIX_),
    for itr = 1:20,
        [cIX_,gIX_] = Regression_Direct(hfig,thres_reg,regressor_last,0);
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

%% row 3: t-tests

function edit_tteststimrange_Callback(hObject,~)
hfig = getParentFigure(hObject);
str = get(hObject,'String');
if ~isempty(str),
    % e.g. '1-3,5', but can't use 'end'
    tteststimrange = ParseRange(str);
    setappdata(hfig,'tteststimrange',tteststimrange);
    %     UpdateTimeIndex(hfig);
    %     RefreshFigure(hfig);
end
end

function edit_ttestthres_Callback(hObject,~)
str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%f');
    thres_ttest = temp{:};
    hfig = getParentFigure(hObject);
    setappdata(hfig,'thres_ttest',thres_ttest);
end
end

function pushbutton_ttest_Callback(hObject,~)
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');
cIX = getappdata(hfig,'cIX');
stim = getappdata(hfig,'stim');
fishset = getappdata(hfig,'fishset');
thres_ttest = getappdata(hfig,'thres_ttest');
tteststimrange = getappdata(hfig,'tteststimrange');

[stimStateBinary_all, names] = GetStimStateBinary(stim,fishset);
stimStateBinary = stimStateBinary_all(tteststimrange,:);
if isempty(stimStateBinary),
    return;
end

%% Misha's method of averaging over 5 frames;
if 1,
    samples = cell(1,length(tteststimrange));
    for i = 1:length(tteststimrange),
        IX = find(stimStateBinary(i,:));
        len = 5*floor(length(IX)/5);
        M_ = M(:,IX(1:len));
        temp = reshape(M_,size(M,1),5,[]);
        samples{i} = squeeze(mean(temp,2));
    end
    
    nCells = size(M,1);
    pvals = zeros(1,nCells);
    signs = zeros(1,nCells);
    for i_cell = 1:nCells,
        [~, p] = ttest2(samples{1}(i_cell,:),samples{2}(i_cell,:));
        pvals(i_cell) = p;
        signs(i_cell) = sign(mean(samples{1}(i_cell,:)) -  mean(samples{2}(i_cell,:)));
    end
else
    %% Linear Discrimination Analysis with HotellingT2 test
    disp('t-test with HotellingT2...');
    % reshape
    sampling_interval = 5; % (down-sampling the time-dimension)
    M_3D = cell(1,length(tteststimrange));
    nCells = size(M,1);
    for i = 1:length(tteststimrange),
        diffStim = horzcat(stimStateBinary(i,1)==1,diff(stimStateBinary(i,:)));
        startIX = find(diffStim==1);
        stopIX = find(diffStim==-1);
        nBlocks = min(length(startIX),length(stopIX));
        startIX = startIX(1:nBlocks);
        stopIX = stopIX(1:nBlocks);
        
        intervals = stopIX-startIX;
        period = mode(intervals);
        IX = find(intervals==period);
        if length(IX)<length(startIX),
            startIX = startIX(IX);
            stopIX = stopIX(IX);
        end
        period_effective = length(1:sampling_interval:period);
        M_3D{i} = zeros(nCells,period_effective,nBlocks);
        for i_block = 1:nBlocks,
            M_3D{i}(:,:,i_block) = M(:,startIX(i_block):sampling_interval:stopIX(i_block)-1);
        end
    end
    
    % Test with 'HotellingT2' function
    pvals = zeros(1,nCells);
    tic
    for i_cell = 1:nCells,
        M_1 = squeeze(M_3D{1}(i_cell,:,:))'; % rows ~ observations; col ~ dimensions
        M_2 = squeeze(M_3D{2}(i_cell,:,:))'; % rows ~ observations; col ~ dimensions
        
        %         Miu1 = mean(M_1,1)';
        %         Covar1 = cov(M_1);
        %         Miu2 = mean(M_2,1)';
        %         Covar2 = cov(M_2);
        %         Sigma_inv = inv(0.5*(Covar1+Covar2));
        %         w = Sigma_inv * (Miu1 - Miu2);
        
        % format input matrix for HotellingT2
        X1 = horzcat(ones(size(M_1,1),1),M_1);
        X2 = horzcat(2*ones(size(M_2,1),1),M_2);
        X = vertcat(X1,X2);
        % direct call to branch-function of 'HotellingT2' package
        pvals(i_cell) = T2Hot2ihe_Direct(X,thres_ttest);
    end
    toc
end

% plot
mean(pvals)
std(pvals)
% set lower bound
pvals_ = pvals;
pvals_(find(pvals==0)) = 10^-100;
figure;hold on;
xv = -100:1:0;
counts = hist(log(pvals_),xv);
bar(-100:1:0,counts);
plot([log(thres_ttest),log(thres_ttest)],[0,max(counts)],'r--');

% rank cells, and threshold if desired
[B,I] = sort(pvals_);
cIX = cIX(I);
if thres_ttest>0,
    ix = find(B<thres_ttest,1,'last');
    cIX = cIX(1:ix);
end
if length(cIX)>100,
    interval = round(length(cIX)/20);
    gIX = ceil((1:length(cIX))/interval)';
    gIX(gIX>20) = 20;
else
    gIX = ones(size(cIX));
end
UpdateIndices(hfig,cIX,gIX);
RefreshFigure(hfig);


end

%% ----- tab four ----- (Clustering etc.)

%% row 1: k-means

function edit_kmeans_Callback(hObject,~)
watchon; drawnow;
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
        [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        [gIX,C] = kmeans(M,numK,'distance','correlation');%,'Replicates',3);
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
watchoff
end

function edit_kmeans2_Callback(hObject,~) % based on anatomical distance
hfig = getParentFigure(hObject);
M = getappdata(hfig,'M');
% z_res = getappdata(hfig,'z_res');
cIX = getappdata(hfig,'cIX');
CellXYZ = getappdata(hfig,'CellXYZ');
absIX = getappdata(hfig,'absIX');
cIX_abs = absIX(cIX);
Combo = horzcat(CellXYZ(cIX_abs,:)/50,M);% how to weight????????????????

str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d');
    numK = temp{:};
    disp(['anat. weighted k-means k=' num2str(numK) '...']);
    tic
    rng('default');% default = 0, but can try different seeds if doesn't converge
    [gIX,C] = kmeans(Combo,numK,'distance','correlation');%,'Replicates',3);
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
        [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        [gIX,C] = kmeans(M,numK,'distance','correlation');%,'Replicates',3);
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
[C,D] = FindCentroid(hfig);

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
M_0 = getappdata(hfig,'M_0');

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
    SaveCluster_Direct(hfig,cIX,gIX,['clean_round' num2str(round)]);
    
    SaveCluster_Direct(hfig,I_rest,ones(length(I_rest),1),['rest_round' num2str(round)]);
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
            I_clean_s = [I_clean_s; IX_s]; %#ok<*AGROW>
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
M_0 = getappdata(hfig,'M_0');
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
    SaveCluster_Direct(hfig,cIX,gIX,'k=20');
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
SaveCluster_Direct(hfig,cIX,gIX,['clean_round' num2str(iter)]);
SaveCluster_Direct(hfig,I_rest,ones(length(I_rest),1),['rest_round' num2str(iter)]);

[gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX);

%% rank by stim-lock
disp('stim-lock');
UpdateIndices(hfig,cIX,gIX,numU);
[gIX,rankscore] = RankByStimLock_Direct(hfig,gIX,numU);
disp('ranking complete');
% and threshold
IX = find(rankscore<thres_stimlock);
ix = ismember(gIX,IX);
gIX = gIX(ix);
cIX = cIX(ix);
UpdateIndices(hfig,cIX,gIX);

%% Regression with the centroid of each cluster
[cIX,gIX,~] = AllCentroidRegression_direct(hfig);
disp('auto-reg-clus complete');

[gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX);
SaveCluster_Direct(hfig,cIX,gIX,'clean_round2');

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

SaveCluster_Direct(hfig,cIX,gIX,'clean_round3');
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
[C,D] = FindCentroid_Direct(gIX,M);
i = 1;
while i<numU,
    c = corr(C(i,:)',C(i+1,:)');
    if c > thres_merge,
        IX = find(gIX == U(i+1));
        gIX(IX)=U(i); %#ok<*FNDSB>
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
gIX = getappdata(hfig,'gIX');
hierinplace = getappdata(hfig,'hierinplace');

[gIX, numU] = SqueezeGroupIX(gIX);
C = FindCentroid(hfig);

str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%d',1);
    numCuts = temp{:};
    
    %     tree = linkage(C,'average','correlation');
    %     IX_tree = cluster(tree,'maxclust',numCuts);
    IX_tree = clusterdata(C,'criterion','distance','distance','correlation','maxclust',numCuts);
    
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
gIX = getappdata(hfig,'gIX');
hierinplace = getappdata(hfig,'hierinplace');

[gIX, numU] = SqueezeGroupIX(gIX);
C = FindCentroid(hfig);

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
C = FindCentroid(hfig);
coeffs = corr(C');%corr(C(1,:)',C(2,:)')
figure('Position',[1000,200,500,500]);
CorrPlot(coeffs);
end

% function CorrPlot(coeffs)
% im = coeffs;
% 
% % red-white-blue colormap
% cmap = zeros(64,3);
% cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
% cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
% cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];
% minlim = -1; %min(min(im));
% maxlim = 1; %max(max(im));
% 
% RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap
% 
% image(RGB); axis equal; axis tight;
% set(gca,'XTick',1:length(im),'YTick',1:length(im));
% 
% for i = 1:size(im,2), % horizontal.. because of image axis
%     for j = 1:size(im,1),
%         text(i-0.3, j, num2str(round(im(i,j)*100)/100));%, 'Units', 'data')
%     end
% end
% end

% function RGB = ImageToRGB(im,cmap,minlim,maxlim)
% L = size(cmap,1);
% ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
% RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
% end

%% ----- tab five ----- (Saved Clusters)

%% row 1: Cluster-Group

function popup_clusgroupmenu_Callback(hObject,~)
hfig = getParentFigure(hObject);
clusgroupID = getappdata(hfig,'clusgroupID');
new_clusgroupID = get(hObject,'Value')-1;
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID);
end

function edit_clusgroupname_Callback(hObject,~)
str = get(hObject,'String');

hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');
clusgroupID = getappdata(hfig,'clusgroupID');
global VAR;
VAR(i_fish).ClusGroupName(clusgroupID) = str; %{}???

menu = MakeNumberedMenu(VAR(i_fish).ClusGroupName);
global hclusgroupmenu;
set(hclusgroupmenu,'String', menu);
end

function pushbutton_newclusgroup_Callback(hObject,~)
hfig = getParentFigure(hObject);
i_fish = getappdata(hfig,'i_fish');
ClusGroup = getappdata(hfig,'ClusGroup');
clusgroupID = getappdata(hfig,'clusgroupID');

% create new Group
ClusGroup = [ClusGroup,cell(1,1)];
setappdata(hfig,'ClusGroup',ClusGroup);
new_clusgroupID = length(ClusGroup); % = oldlength + 1
global VAR;
VAR(i_fish).ClusGroupName(new_clusgroupID) = {'(blank)'};

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

%% row 2: Clusters

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
M_0 = getappdata(hfig,'M_0');
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
            range = [range,m(i)];
        else % have '-'sign,
            range = [range,m(i-1)+1:-m(i)];
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

function [cIX,gIX,numK] = SmartUnique(CIX,GIX,M)
% input: simply concatenated groups
% output: every cell only appears once - in group with the highest
% correlation to group mean/centroid
disp('unique based on corr coeff...');
CTRD = FindCentroid_Direct(GIX,M);

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

%% ----- tab five ----- (Atlas)

function pushbutton_findmasks_Callback(hObject,~)
hfig = getParentFigure(hObject);
MASKs = getappdata(hfig,'MASKs');
% load to shorter variable names
height = MASKs.height; % 1406;
width = MASKs.width; % 621;
Zs = MASKs.Zs; % 138;

cIX = getappdata(hfig,'cIX');
absIX = getappdata(hfig,'absIX');
% isRefAnat = getappdata(hfig,'isRefAnat');
setappdata(hfig,'isRefAnat',1);

%% load normalized data:
% if ~isRefAnat, % make fake data...
%     errordlg('not using normalized cell coordinates!')
%     CellXYZ = getappdata(hfig,'CellXYZ');
%     anat_stack = getappdata(hfig,'anat_stack');
%     [s1,s2,s3] = size(anat_stack);
%     
%     % fake
%     X_raw = CellXYZ(absIX(cIX),1);
%     Y_raw = CellXYZ(absIX(cIX),2);
%     Z_raw =  CellXYZ(absIX(cIX),3);
% 
%     X = ceil(X_raw.*((height-1)/(s1-1)));
%     Y = ceil(Y_raw.*((width-1)/(s2-1)));
%     Z = ceil(1+(Z_raw-1).*((Zs-1)/(s3-1)));
% else
    CellXYZ = getappdata(hfig,'CellXYZ_norm');
    X = CellXYZ(absIX(cIX),1);
    Y = CellXYZ(absIX(cIX),2);
    Z = CellXYZ(absIX(cIX),3);
% end

%%
pxID = sub2ind([height,width,Zs],X',Y',Z');
Msk_hist = full(sum(MASKs.MaskDatabase(pxID,:),1));

% choose top hits
thres_2std = mean(Msk_hist)+2*std(Msk_hist);
Msk_IDs = find(Msk_hist>thres_2std);

% plot histogram
figure('Position',[300,600,1000,300]);
hold on;
nMasks = length(Msk_hist);
plot([1,nMasks],[thres_2std/sum(Msk_hist),thres_2std/sum(Msk_hist)],'r:');
yspace = max(Msk_hist(Msk_IDs)/sum(Msk_hist))/50;
plot(Msk_IDs,Msk_hist(Msk_IDs)/sum(Msk_hist)+yspace,'r.')
bar(Msk_hist/sum(Msk_hist),'facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5]);
xlim([-1,nMasks+2]);
xlabel('mask ID');
ylabel('probability')

% display names in command window
names = MASKs.MaskDatabaseNames(Msk_IDs);
names_numbered = cell(size(names));
disp('Regions:');
for i = 1:length(names),
    names_numbered{i} = [num2str(Msk_IDs(i)), ': ', names{i}];
    disp(names_numbered{i});
end

setappdata(hfig,'Msk_hist',Msk_hist);
setappdata(hfig,'Msk_IDs',Msk_IDs);
setappdata(hfig,'Msk_IDlist',Msk_IDs);

global hmasklistbox;
set(hmasklistbox,'String',names_numbered);

RefreshFigure(hfig);
end

function checkbox_showmasks_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isShowMasks',get(hObject,'Value'));
RefreshFigure(hfig);
end

function edit_chooseMskID_Callback(hObject,~)
hfig = getParentFigure(hObject);
MASKs = getappdata(hfig,'MASKs');
nMasks = size(MASKs.MaskDatabase,2);
% get/format range
str = get(hObject,'String');
if ~isempty(str),
    str = strrep(str,'end',num2str(nMasks));
    range = ParseRange(str);
    
    setappdata(hfig,'Msk_IDs',range);
    RefreshFigure(hfig);
end
end

function listbox_chooseMskID_Callback(hObject,~)
hfig = getParentFigure(hObject);
Msk_IDlist = getappdata(hfig,'Msk_IDlist');
% names = get(hObject,'String');
IX = get(hObject,'Value');
setappdata(hfig,'Msk_IDs',Msk_IDlist(IX));
RefreshFigure(hfig);
end

function checkbox_showmskoutline_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isShowMskOutline',get(hObject,'Value'));
RefreshFigure(hfig);
end

%% Internal functions

function UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID,norefresh) %#ok<INUSD>
% save/update old Cluster into ClusGroup before exiting,
% as Cluster is the variable handled in hfig but not saved elsewhere
ClusGroup = getappdata(hfig,'ClusGroup');
Cluster = getappdata(hfig,'Cluster');
i_fish = getappdata(hfig,'i_fish');

% update into workspace
ClusGroup{clusgroupID} = Cluster;
setappdata(hfig,'ClusGroup',ClusGroup);
global VAR;
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);

% load new 'Cluster'
Cluster = ClusGroup{new_clusgroupID};
setappdata(hfig,'Cluster',Cluster);
setappdata(hfig,'clusgroupID',new_clusgroupID);

% update GUI: hclusgroupmenu
global hclusgroupmenu hclusgroupname;
if ishandle(hclusgroupmenu),
    menu = MakeNumberedMenu(VAR(i_fish).ClusGroupName);
    set(hclusgroupmenu,'String',menu,'Value',new_clusgroupID+1);
    set(hclusgroupname,'String',VAR(i_fish).ClusGroupName(new_clusgroupID));
end

if ~exist('norefresh','var'),
    if numel(Cluster) == 0, % i.e. for newly created ClusGroup
        SaveCluster(hfig,'new');
    else % load this ClusGroup
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
absIX = getappdata(hfig,'absIX');
Cluster = getappdata(hfig,'Cluster');

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

cIX_abs = absIX(cIX);

Cluster(clusID).cIX_abs = cIX_abs;
Cluster(clusID).gIX = gIX;
Cluster(clusID).numel = length(cIX_abs);
U = unique(gIX);
numU = length(U);
Cluster(clusID).numK = numU;

setappdata(hfig,'Cluster',Cluster);
UpdateClusID(hfig,clusID);
disp('cluster saved');
end

function [Cluster,clusID] = SaveCluster_Direct(hfig,cIX,gIX,name) %,clusheader,name)
new_clusgroupID = 1;
clusgroupID = getappdata(hfig,'clusgroupID');
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID);

if ~exist('cIX','var'),
    cIX = getappdata(hfig,'cIX');
end
if ~exist('gIX','var'),
    gIX = getappdata(hfig,'gIX');
end

Cluster = getappdata(hfig,'Cluster');
absIX = getappdata(hfig,'absIX');
cIX_abs = absIX(cIX);

clusID = numel(Cluster)+1;
% if ~exist('name','var'),
%     name = getappdata(hfig,'newclusname');
% end
% if ~exist('clusheader','var'),
%     clusheader = getappdata(hfig,'clusheader');
% end
Cluster(clusID).name = name; %[clusheader name];
Cluster(clusID).cIX_abs = cIX_abs;
Cluster(clusID).gIX = gIX;
Cluster(clusID).numel = length(cIX);
Cluster(clusID).numK = length(unique(gIX));

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
menu = MakeNumberedMenu({Cluster.name});
set(hclusmenu,'String', menu,'Value',clusID+1);
numK = Cluster(clusID).numK;
gIX = Cluster(clusID).gIX;

% convert absolute index to index used for this dataset
absIX = getappdata(hfig,'absIX');
cIX_abs = Cluster(clusID).cIX_abs;
[~,cIX] = ismember(cIX_abs,absIX);

if ~isempty(find(cIX==0,1)),
    errordlg('cell index out of bound for currently loaded dataset');
    IX = cIX==0;
    cIX(IX) = [];
    gIX(IX) = [];
end

UpdateIndices(hfig,cIX,gIX,numK);
RefreshFigure(hfig);
end

function menu = MakeNumberedMenu(name) % e.g. name = {Cluster.name} (note {})
menu = [{'(choose)'},name];
for j=2:length(menu),menu(j)={[num2str(j-1) ': ' menu{j}]};end
end

function [gIX, numU] = HierClus(M,gIX,isplotfig) %#ok<INUSD>
[gIX, numU] = SqueezeGroupIX(gIX);
[C,~] = FindCentroid_Direct(gIX,M);
D = pdist(C,'correlation');
tree = linkage(C,'average','correlation');
leafOrder = optimalleaforder(tree,D);

if numU>1,
    if exist('isplotfig','var'),
        figure('Position',[100 100 600 600]);
        %             subplot(1,3,1);
        %             CORR = corr(C');
        %             CorrPlot(CORR);
        %
        %             subplot(1,3,2);
        dendrogram(tree,numU,'orientation','right','reorder',leafOrder);
        set(gca,'YDir','reverse');
        set(gca,'XTick',[]);
        
        %             subplot(1,3,3);
        %             C2 = C(leafOrder,:);
        %             CORR2 = corr(C2');
        %             CorrPlot(CORR2);
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
            M_0 = getappdata(hfig,'M_0');
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

setappdata(hfig,'bCache',bC);
setappdata(hfig,'cIX',cIX);
setappdata(hfig,'gIX',gIX);

M = GetTimeIndexedData(hfig);
setappdata(hfig,'M',M);

if exist('numK','var'),
    setappdata(hfig,'numK',double(numK));
end

%% Resets: reset flags the NEXT time this function is called (so they only apply to this particular plot)
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

% toggle 'isWeighAlpha'
isWeighAlpha = getappdata(hfig,'isWeighAlpha');
if isWeighAlpha == 1,
    setappdata(hfig,'isWeighAlpha',100);
elseif isWeighAlpha == 100,
    setappdata(hfig,'isWeighAlpha',0);
end

% FindCentroid reset:
setappdata(hfig,'Centroids',[]);
setappdata(hfig,'D_ctrd',[]);
end

% frequently used, 2 plotting functions are outside ('DrawTimeSeries.m' and 'DrawCellsOnAnatProj.m')
function RefreshFigure(hfig)
watchon; drawnow;
isPopout = 0; % with down-sampling in plots

% clean-up canvas
allAxesInFigure = findall(hfig,'type','axes');
if ~isempty(allAxesInFigure)
    delete(allAxesInFigure);
end

figure(hfig);
h1 = axes('Position',[0.05, 0.04, 0.55, 0.83]);
h2 = axes('Position',[0.63, 0.04, 0.35, 0.83]);

isCentroid = getappdata(hfig,'isCentroid');
isRefAnat = getappdata(hfig,'isRefAnat');

isPlotLines = 0; %getappdata(hfig,'isPlotLines');
isPlotBehavior = 1; %getappdata(hfig,'isPlotBehavior');

% double-check if cIX is valid
cIX = getappdata(hfig,'cIX');
if isempty(cIX),
    errordlg('empty set!');
    % GO BACK to the last step (presumably not empty)
    pushbutton_back_Callback(h1); % using h1 instaed of the usual 'hObject'
    return;
end

% left subplot
axes(h1);
DrawTimeSeries(hfig,h1,isPopout,isCentroid,isPlotLines,isPlotBehavior);

% right subplot
axes(h2);
DrawCellsOnAnatProj(hfig,isRefAnat,isPopout);
watchoff;
end

function [C,D] = FindCentroid_Direct(gIX,M)
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

function UpdateTimeIndex(hfig,isSkipcIX) %#ok<INUSD>
% input params
isAvr = getappdata(hfig,'isAvr');
isRawtime = getappdata(hfig,'isRawtime');
stimrange = getappdata(hfig,'stimrange');
% load
timelists = getappdata(hfig,'timelists');
periods = getappdata(hfig,'periods');
fishset = getappdata(hfig,'fishset');

if fishset == 1,
    if isAvr,
        tIX = 1:periods;
    else
        tIX = timelists{1};
    end
    
else % fishset>1,
    if isAvr,
        tIX = [];
        for i = 1:length(stimrange),
            ix = stimrange(i);
            i_start = sum(periods(1:ix-1)); % if ix-1<1, sum = 0
            tIX = horzcat(tIX,(i_start+1:i_start+periods(ix)));
            %             tIX = vertcat(tIX,(i_start+1:i_start+periods(ix))');
        end
    else % full range
        if ~isRawtime,
            tIX = cat(2, timelists{stimrange});
        else
            tIX = sort(cat(2, timelists{stimrange}));
        end
    end
end

setappdata(hfig,'tIX',tIX);

% set Matrices to hold time-series
M_0 = GetTimeIndexedData(hfig,'isAllCells');
setappdata(hfig,'M_0',M_0);
if ~exist('isSkipcIX','var'),
    cIX = getappdata(hfig,'cIX');
    setappdata(hfig,'M',M_0(cIX,:));
end
end

function [M,behavior,stim] = GetTimeIndexedData(hfig,isAllCells) %#ok<INUSD>
%{
% naming convention used:
M = GetTimeIndexedData(hfig);
M_0 = GetTimeIndexedData(hfig,'isAllCells');
%}

isZscore = getappdata(hfig,'isZscore');
% main data input
if ~isZscore,
    cellResp = getappdata(hfig,'CellResp');
    cellRespAvr = getappdata(hfig,'CellRespAvr');
else
    cellResp = getappdata(hfig,'CellRespZ');
    cellRespAvr = getappdata(hfig,'CellRespAvrZ');
end
Behavior_full = getappdata(hfig,'Behavior_full');
BehaviorAvr = getappdata(hfig,'BehaviorAvr');
stim_full = getappdata(hfig,'stim_full');
stimAvr = getappdata(hfig,'stimAvr');
% other params
isAvr = getappdata(hfig,'isAvr');
cIX = getappdata(hfig,'cIX');
tIX = getappdata(hfig,'tIX');

%% set data
if isAvr,
    if exist('isAllCells','var'),
        M = cellRespAvr(:,tIX);
    else
        M = cellRespAvr(cIX,tIX);
    end
    behavior = BehaviorAvr(:,tIX);
    stim = stimAvr(:,tIX);
else
    if exist('isAllCells','var'),
        M = cellResp(:,tIX);
    else
        M = cellResp(cIX,tIX);
    end
    behavior = Behavior_full(:,tIX);
    stim = stim_full(:,tIX);
end

setappdata(hfig,'behavior',behavior);
setappdata(hfig,'stim',stim);
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
