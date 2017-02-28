% fig2A: correlation to phototaxis stimulus
% for single regressor, detail
% option set to do correlation on trial-averaged traces
% plot corresponding histogram of corr.coeffs, 
% and plot functional as well as anat plot

% compare to fig3A_sensory_motor, which for now plots a pair of regressors
% in a double colormap, and batch saves figures.



hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
i_fish = 6;
ClusterIDs = [2,1];
[cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);

%% corr with PT regressor, stim-avr

setappdata(hfig,'isStimAvr',1);
[M_0,M,behavior,stim] = UpdateTimeIndex(hfig);

fishset = getappdata(hfig,'fishset');
[~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
reg_range = 2;        
Reg = regressors(reg_range,:);
R = corr(Reg',M_0');       

thres_reg = 0.5;
cIX = (find(R>thres_reg))';
wIX = R(cIX); % weight, i.e. corr.coeff
[~,I] = sort(wIX,'descend');
cIX = cIX(I);
gIX = ones(size(cIX));
wIX = wIX(I)';

UpdateIndices_Manual(hfig,cIX,gIX);

%% histogram of corr.coeff, for all cells
figure('Position',[500,200,500,150]);
hold on;
bins = -1:0.05:1;
[N,~] = histcounts(R,bins);
histogram(R,bins,'FaceColor',[0.4 0.4 0.4]);%,'EdgeColor','none'
plot([thres_reg,thres_reg],[0,max(N)],'r--');
xlim([-1,1]);ylim([0,max(N)]);

%% histogram with rainbow colors
clrbins = 0.5:0.05:1;
cmap = hot(length(clrbins)-1);
n_blank = length(bins)-length(clrbins);
cmap = vertcat(ones(n_blank,3),cmap);
DrawRainbowHistogram(N,bins,cmap);

%% left plot
figure('Position',[50,500,300,200])%[50,100,800,1000]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
setappdata(hfig,'isPlotBehavior',0);
setappdata(hfig,'isStimAvr',1);
UpdateTimeIndex(hfig);

% set color
clrmap = [0,0,0];
numK = max(gIX);
assert(size(clrmap,1)==numK);

opts = [];
opts.clrmap = clrmap;
opts.isAxTight = 1;
DrawTimeSeries(hfig,cIX,gIX,opts);
       
%% set up custom colormap
clrbins = 0.5:0.05:1;
clrmap = hot(length(clrbins));
gIX = interp1(clrbins,1:length(clrbins),wIX,'nearest');

%% right plot
figure('Position',[50,100,800,1000]);
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
I.clrmap = clrmap;
DrawCellsOnAnat(I);

%%
DrawCustomColorbar(clrmap,bins,numTicks);

%%




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

%%
%% motor regressor
%             [~,names,regressors] = GetMotorRegressor(behavior,i_fish);
        
