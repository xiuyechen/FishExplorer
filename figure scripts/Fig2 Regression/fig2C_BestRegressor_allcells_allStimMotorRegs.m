% fig2C: best clusters/cells? for an array of stim/motor regressors, GUI function
clear all;close all;clc;

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

setappdata(hfig,'isMotorseed',1);
i_fish = 8;
LoadSingleFishDefault(i_fish,hfig);

%% regression with all regs
% ResetDisplayParams(hfig,i_fish);
setappdata(hfig,'stimrange',1:2);
UpdateTimeIndex(hfig);

isRegIndividualCells = 1;
isRegCurrentCells = 0;
setappdata(hfig,'thres_reg',0.4);
[cIX_reg,gIX_reg,numK,IX_regtype,corr_max] = AllRegsRegression(hfig,isRegIndividualCells,isRegCurrentCells);

regtypes_plot = [2,3,8,9,16,17];%[2,3,4,7,8,9,16,17]; % manual input % PT, OMR, joint
[cIX,gIX] = SelectClusterRange(cIX_reg,gIX_reg,regtypes_plot);
[gIX, numU] = SqueezeGroupIX(gIX);
% manually switch order for the last two groups (motor)!
m1 = regtypes_plot(end-1);
m2 = regtypes_plot(end);
gIX_old = gIX;
gIX(gIX_old==m1) = m2;
gIX(gIX_old==m2) = m1;

UpdateIndices_Manual(hfig,cIX,gIX,numU);

%% left plot
figure('Position',[50,100,400,500]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',0);
setappdata(hfig,'isStimAvr',1);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX);

%% right plot
% figure('Position',[50,100,800,1000]);
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
DrawCellsOnAnat(I);

%% right plot of PT/OMR joint cells only
IX = union(find(gIX==length(regtypes_plot)-1),find(gIX==length(regtypes_plot)));
cIX2 = cIX(IX);
gIX2 = gIX(IX);
% UpdateIndices_Manual(hfig,cIX2,gIX2,numU);
I = LoadCurrentFishForAnatPlot(hfig,cIX2,gIX2);
DrawCellsOnAnat(I);

%% left plot for the motor traces, not plotted in anat
regtypes_motor = [18,19];
[cIX,gIX] = SelectClusterRange(cIX_reg,gIX_reg,regtypes_motor);
UpdateIndices_Manual(hfig,cIX,gIX,numU);
clrmap = [0.5,0.5,0.5;0.5,0.5,0.5];
opts = [];
opts.clrmap = clrmap;

figure('Position',[50,100,400,500]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',0);
setappdata(hfig,'isStimAvr',1);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX,opts);
% %% left-right combined plot
% figure('Position',[50,100,1400,800]);
% % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
% subplot(121)
% setappdata(hfig,'isPlotBehavior',0);
% setappdata(hfig,'isStimAvr',1);
% UpdateTimeIndex(hfig);
% DrawTimeSeries(hfig,cIX,gIX);
% 
% % right plot
% subplot(122)
% I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
% DrawCellsOnAnat(I);

%%
% dataDir = GetCurrentDataDir;
% saveDir = fullfile(dataDir,'motorsourceplot');
% if ~exist(saveDir, 'dir'), mkdir(saveDir), end;
% filename = fullfile(saveDir, num2str(i_fish));
% saveas(gcf, filename, 'png');
% close(gcf)