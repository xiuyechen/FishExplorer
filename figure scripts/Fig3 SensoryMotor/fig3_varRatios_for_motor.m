
clear all;close all;clc
%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
range_fish = 1:18;
P_ratio = zeros(18,3,2);
P_tAvrcorr= zeros(18,3,2);
P_tRescorr= zeros(18,3,2);
for i_fish = range_fish
    setappdata(hfig,'isMotorseed',0);
    ClusterIDs = [1,1];
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    [bh_tAvr,bh_tRes] = GetTrialAvrLongTrace(hfig,behavior);
    for i_lfr = 1:3
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
        vTot_bh = var(behavior(i_lfr,:));
        P_ratio(i_fish,i_lfr,2) = vAvr_bh/vTot_bh;
        P_tAvrcorr(i_fish,i_lfr,2) = corr(behavior(i_lfr,:)',bh_tAvr(i_lfr,:)');
        P_tRescorr(i_fish,i_lfr,2) = corr(behavior(i_lfr,:)',bh_tRes(i_lfr,:)');
    end    
end
%%
figure;
bar(P_ratio(:,:,2));
ylim([0,1])
ylabel('var.expl by t.Avr')
xlabel('fish ID')
%%
% figure
% subplot(2,2,1)
% bar(P_tAvrcorr(:,:,1));
% subplot(2,2,2)
% bar(P_tAvrcorr(:,:,2));
% subplot(2,2,3)
% bar(P_tRescorr(:,:,1));
% subplot(2,2,4)
% bar(P_tRescorr(:,:,2));