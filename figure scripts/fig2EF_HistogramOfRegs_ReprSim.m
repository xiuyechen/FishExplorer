% [figC: 'AllRegsRegression.m' in 'GUI functions']

% figE: histogram of regressions, different regressors, for all cells
% load a fish with the appropriate stim-range

Data = getappdata(hfig,'M_0');
%         Data = getappdata(hfig,'M');
thres_reg = getappdata(hfig,'thres_reg');
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
% cIX_in = getappdata(hfig,'cIX');
% gIX_in = getappdata(hfig,'gIX');
% numK = getappdata(hfig,'numK');
i_fish = getappdata(hfig,'i_fish');

% get stim/motor regressors
[~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
[~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
Reg = vertcat(regressor_s,regressor_m);
regnames = [names_s,names_m]';

% regression
Corr = corr(Reg',Data');

% make control distribution
IX = randperm(size(Data,2));
Data_shf = Data(:,IX);
Corr_shf = corr(Reg',Data_shf');

%% plot distribution 
figure('Position',[500,200,500,850]);
rangeReg = [2,3,6,7,8,9,38,39,40];
numReg = length(rangeReg);
for i = 1:numReg,
    subplot(numReg,1,i);
    R = Corr(i,:);
    R_shf = Corr_shf(i,:);
    
    hold on;
    bins = -1:0.01:1;%-1:0.025:1;
    [N,~] = histcounts(R,bins);
    histogram(R,bins,'FaceColor',[0.4 0.4 0.4]);%,'EdgeColor','none'
    %     plot([thres_reg,thres_reg],[0,max(N)],'r--');
    [N_shf,~] = histcounts(R_shf,bins);
    histogram(R_shf,bins,'FaceColor',[1 0.4 0.4]);
    ymax = max(max(N),max(N_shf));
    plot([0,0],[0,ymax],'r--');
    xlim([-0.5,1]);ylim([0,ymax]);
    title(regnames{rangeReg(i)})
end

%


%% fig2F: representational similarity
RS = Corr*Corr';
figure;imagesc(RS)
set(gca,'YTick',1:length(regnames),'YTickLabel',regnames,'TickLength',[0,0]);
set(gca,'XTick',1:length(regnames),'XTickLabel',regnames,'XTickLabelRotation',90);
axis equal;axis tight
