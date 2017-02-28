% [figC: 'AllRegsRegression.m' in 'GUI functions']

% figE: histogram of regressions, different regressors, for all cells

%%
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
i_fish = 8;
ClusterIDs = [2,1];
[cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);

%%
% get stim/motor regressors
fishset = getappdata(hfig,'fishset');
[~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);

isMotorseed = 0;
setappdata(hfig,'isMotorseed',isMotorseed);
[~,~,behavior] = UpdateTimeIndex(hfig);

[~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
Reg = vertcat(regressor_s,regressor_m);
regnames = [names_s,names_m]';

% regression
Corr = corr(Reg',M_0');

% make control distribution
IX = randperm(size(M_0,2));
Data_shf = M_0(:,IX);
Corr_shf = corr(Reg',Data_shf');

%% plot distribution 
figure('Position',[500,200,400,850]);
if i_fish==8
    rangeReg = [2,3,7,8,9,18,19,20];
else % i_fish--6
    
    if isMotorseed
        rangeReg = [2,3,5,6,9,10];%[2,3,6,7,8,9,38,39,40];
    else
        rangeReg = [2,3,5,6,9,11];%10,11];
    end
end
numReg = length(rangeReg);
for i = 1:numReg
    i_reg = rangeReg(i);
    
    subplot(numReg,1,i);
    R = Corr(i_reg,:);
    R_shf = Corr_shf(i_reg,:);
    
    hold on;
    bins = -1:0.05:1;%-1:0.025:1;
    [N,~] = histcounts(R,bins);
    histogram(R,bins,'FaceColor',[0.4 0.4 0.4]);%,'EdgeColor','none'
    %     plot([thres_reg,thres_reg],[0,max(N)],'r--');
    [N_shf,~] = histcounts(R_shf,bins);
    histogram(R_shf,bins,'FaceColor',[1 0.4 0.4],'FaceAlpha',0.4);
    ymax = max(max(N),max(N_shf));
    plot([0,0],[0,ymax],'r--');
    xlim([-1,1]);
    ylim([0,ymax]);
    set(gca,'YTick',[])
    ylabel('a.u.')
    title(regnames{rangeReg(i)})
end

%


%% fig2F: representational similarity
data = Corr(rangeReg,:);

D = pdist(data,'correlation');
tree = linkage(data,'average','correlation');
leafOrder = optimalleaforder(tree,D);
    
data2 = data(leafOrder,:);
RS = data2*data2';
names = regnames(rangeReg);
names2 = names(leafOrder);

figure;
set(gcf,'color','w');
imagesc(RS)
set(gca,'YTick',1:length(names2),'YTickLabel',names2,'TickLength',[0,0]);
set(gca,'XTick',1:length(names2),'XTickLabel',names2,'XTickLabelRotation',90);
axis equal;axis tight
