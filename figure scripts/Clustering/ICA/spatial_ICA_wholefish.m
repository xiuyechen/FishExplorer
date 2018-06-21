tic

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
i_fish = 8;
[cIX_0,gIX_0,M] = LoadSingleFishDefault(i_fish,hfig,[2,1]);

% [cIX_0,gIX_0] = LoadCluster_Direct(i_fish,2,1);

cIX = cIX_0;
gIX = gIX_0;

%% PCA
[coeff,score,latent,tsquared,explained,mu] = pca(M);
% figure;subplot(121);imagesc(coeff);subplot(122);imagesc(score);title('M');

% data to pass on from PCA:
k = 100;
Z0 = score(:,1:k);
% figure;imagesc(Z0);
t1 = toc

%% ICA
tic
numK = 139;%20;
Z = Z0';
[Zica, W, T, mu] = fastICA(Z,numK);

t2 = toc
%%
n = size(M',2);
Zr = T \ W' * Zica + repmat(mu,1,n);

figure;
subplot(131)
imagesc(Zica);
subplot(132)
imagesc(W);
subplot(133)
imagesc(Zr);

%% clustering
gIX = kmeans(Zr',numK,'distance','correlation');

UpdateIndices_Manual(hfig,cIX,gIX);

figure('Position',[50,100,1400,800]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
subplot(121)
setappdata(hfig,'isPlotBehavior',1);
setappdata(hfig,'isStimAvr',0);
setappdata(hfig,'isPlotLines',0);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX);

% right plot
subplot(122)
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
DrawCellsOnAnat(I);

% kmeansPlot(Zr,gIX)
%% bad? pick top 100 cells for each component
% gIX = zeros(size(gIX));
% for i = 1:numK
%     [B,IX] = sort(Zica(i,:),'descend');
%     IX_keep = IX(1:100); % this is too crude ~~~~~~
%     
%     gIX(IX_keep) = i;
% end
% IX = (gIX~=0);
% cIX = cIX(IX);
% gIX = gIX(IX);

%% pick top ? cells for each component % still crude, no 
cutoff_std = 3;

Z2 = Zica;Z2(Z2<cutoff_std)=nan;
% manual code to not assign values to NaN columns (facepalm)
[~,gIX_1] = max(vertcat(nan(1,size(Z2,2)),Z2));
% gIX_1 = gIX_1';
gIX = gIX_1-1;
IX = (gIX~=0);
cIX = cIX_0(IX);
gIX = gIX_1(IX);

%%
% [cIX,gIX] = SelectClusterRange(cIX,gIX,[3,4]);

UpdateIndices_Manual(hfig,cIX,gIX);
%%
figure('Position',[50,100,1400,800]);
% isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
subplot(121)
setappdata(hfig,'isPlotBehavior',1);
setappdata(hfig,'isStimAvr',0);
setappdata(hfig,'isPlotLines',0);
UpdateTimeIndex(hfig);
DrawTimeSeries(hfig,cIX,gIX);

% right plot
subplot(122)
I = LoadCurrentFishForAnatPlot(hfig,cIX,gIX);
DrawCellsOnAnat(I);
