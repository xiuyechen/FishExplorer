% fig2C: best clusters/cells? for an array of stim/motor regressors, GUI function

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
i_fish = 8;
LoadSingleFishDefault(i_fish,hfig);

%% regression with all regs
% ResetDisplayParams(hfig,i_fish);
setappdata(hfig,'stimrange',1:2);
UpdateTimeIndex(hfig);

isRegIndividualCells = 1;
isRegCurrentCells = 0;
setappdata(hfig,'thres_reg',0.4);
[cIX,gIX,numK,IX_regtype,corr_max] = AllRegsRegression(hfig,isRegIndividualCells,isRegCurrentCells);

regtypes_keep = [2,3,4,7,8,9,16,17]; % manual input
[cIX,gIX] = SelectClusterRange(cIX,gIX,regtypes_keep);
[gIX, numU] = SqueezeGroupIX(gIX);
UpdateIndices_Manual(hfig,cIX,gIX,numU);

%% left plot
figure('Position',[50,100,800,1000]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',0);
setappdata(hfig,'isStimAvr',1);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX);

%% right plot
figure('Position',[50,100,800,1000]);
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
DrawCellsOnAnat(I);

%%
% dataDir = GetCurrentDataDir;
% saveDir = fullfile(dataDir,'motorsourceplot');
% if ~exist(saveDir, 'dir'), mkdir(saveDir), end;
% filename = fullfile(saveDir, num2str(i_fish));
% saveas(gcf, filename, 'png');
% close(gcf)