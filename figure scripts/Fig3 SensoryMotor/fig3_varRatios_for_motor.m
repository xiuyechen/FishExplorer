
clear all;close all;clc
%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
range_fish = [1:3,5:18];
P_ratio = zeros(18,5,2);
P_ratio_tAvr = zeros(18,5,2);
P_ratio_tRes = zeros(18,5,2);
P_tAvrcorr= zeros(18,5,2);
P_tRescorr= zeros(18,5,2);

for i_fish = range_fish
    setappdata(hfig,'isMotorseed',0);
    ClusterIDs = [1,1];
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    [bh_tAvr,bh_tRes] = GetTrialAvrLongTrace(hfig,behavior);
    for i_lfr = 1:5
        vAvr_bh = var(bh_tAvr(i_lfr,:));
        vTot_bh = var(behavior(i_lfr,:));
        P_ratio(i_fish,i_lfr,1) = vAvr_bh/vTot_bh;
        P_tAvrcorr(i_fish,i_lfr,1) = corr(behavior(i_lfr,:)',bh_tAvr(i_lfr,:)');
        P_tRescorr(i_fish,i_lfr,1) = corr(behavior(i_lfr,:)',bh_tRes(i_lfr,:)');
    end
    
    %%
    setappdata(hfig,'isMotorseed',1);
    [~,~,behavior] = UpdateTimeIndex(hfig);
    [bh_tAvr,bh_tRes] = GetTrialAvrLongTrace(hfig,behavior);
    for i_lfr = 1:2
        vAvr_bh = var(bh_tAvr(i_lfr,:));
        vRes_bh = var(bh_tRes(i_lfr,:));
        vTot_bh = var(behavior(i_lfr,:));
        P_ratio(i_fish,i_lfr,2) = vAvr_bh/vTot_bh;
%         P_ratio_tAvr(i_fish,i_lfr,2) = vAvr_bh/vTot_bh;
%         P_ratio_tRes(i_fish,i_lfr,2) = vRes_bh/vTot_bh; 
% these two add to one
        P_tAvrcorr(i_fish,i_lfr,2) = corr(behavior(i_lfr,:)',bh_tAvr(i_lfr,:)');
        P_tRescorr(i_fish,i_lfr,2) = corr(behavior(i_lfr,:)',bh_tRes(i_lfr,:)');
    end
end

%%
save('varRatios.mat','P_ratio','P_tAvrcorr','P_tRescorr');

%% bar plot by fish (l/r in different colors)
% figure;
% bar(P_ratio(:,:,2));
% ylim([0,1])
% ylabel('var.expl by t.Avr')
% xlabel('fish ID')

%% 2D plot left against right, motorseed/fictive/difference
switch 1
    case 1
        data = P_ratio(range_fish,1:2,2);
        %         data = P_tAvrcorr(range_fish,1:2,2);
    case 2
        %          data = P_ratio(:,[1,3],1);
        data = P_tAvrcorr(:,[1,3],1);
    case 3
        %         data = P_tAvrcorr(:,1:2,2)-P_tAvrcorr(:,[1,3],1);
end

figure('Position',[100,100,150,150]);
hold on;
scatter(data(:,1),data(:,2),5,[0.5,0.5,0.5])
axis equal
xlim([0,1])
xlabel('left')
ylim([0,1])
ylabel('right')
plot([0;1],[0;1],'r:')

%% compare motorseed with raw fictive LR
switch 2
    case 1 
        P = P_ratio;
    case 2 
        P = P_tAvrcorr
end

figure('Position',[100,100,450,150]);
for i = 1:3
    switch i
        case 1
             data = P(range_fish,[4,5],1);
        case 2
            data = P(range_fish,[1,3],1);
        case 3
            data = P(range_fish,1:2,2);
%             data = P(range_fish,1:2,2)-P(range_fish,[4,5],1);
    end
    
    subplot(1,3,i)
    hold on;
    scatter(data(:,1),data(:,2),5,[0.5,0.5,0.5])
    axis equal
    xlim([0,1])
    xlabel('left')
    ylim([0,1])
    ylabel('right')
    plot([0;1],[0;1],'r:')
end

%% single histogram (pooling left and right)
switch 4
    case 1 
        P = P_ratio;
        xlabelname = '%var expl. by stim';
    case 2
        P = 1-P_ratio;
        xlabelname = '%var expl. by residual';
    case 3
        P = P_tAvrcorr;
        xlabelname = 'corr(motor, motor.tAvr)';
    case 4
        P = P_tRescorr;
        xlabelname = 'corr(motor, motor.tRes)';
end

figure('Position',[100,100,550,150]);
for i = 1:3
    switch i
        case 1
            data = P(range_fish,4:5,1);
        case 2
            data = P(range_fish,[1,3],1);
        case 3
            data = P(range_fish,1:2,2);
    end
    
    subplot(1,3,i)
    hold on;
    histogram(data(:),0:0.1:1,'facecolor',[0.5,0.5,0.5],'edgecolor','w')
    xlabel(xlabelname)
    ylabel('count')
    box off
    set(gca,'Layer','top')

end

% data = P(range_fish,1:3,1);
% data = P(range_fish,1:2,2);
% figure('Position',[200,200,150,120]);
% histogram(data(:),0:0.1:1,'facecolor',[0.5,0.5,0.5],'edgecolor','w')
% xlabel('%var expl. by stim')
% ylabel('count')
% box off

%%
% P1 = P_tAvrcorr
%  P2 = P_tRescorr;
range_fish = [1:3,5:18];
  data1 = P_tAvrcorr(range_fish,1:2,2);
 data2 = P_tRescorr(range_fish,1:2,2);

 figure('Position',[400,400,100,100]);
 scatter(data1(:),data2(:),20,[1,0.5,0.5])
 axis equal
 xlim([0,1])
 ylim([0,1])
%% raw data summaries
% figure
% subplot(2,2,1)
% bar(P_tAvrcorr(:,:,1));
% subplot(2,2,2)
% bar(P_tAvrcorr(:,:,2));
% subplot(2,2,3)
% bar(P_tRescorr(:,:,1));
% subplot(2,2,4)
% bar(P_tRescorr(:,:,2));
% %%
% figure
% bar(P_tAvrcorr(:,1:2,2)-P_tAvrcorr(:,[1,3],1))
