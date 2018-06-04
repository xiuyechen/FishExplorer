% this script is directly adapted from fig3_varRatios_for_motor

clear all; close all; clc

%% folder setup
outputDir = GetOutputDataDir;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];%

M_stim = {'5','P','O','L','D','Y'};
M_stimname = {'phT-DF','phT','OMR','Loom','DF','Dot'};

nStim = length(M_stim);
% for fictive channels (5 channels)
P_tAvrcorr_5= nan(18,5,nStim);
P_tRescorr_5= nan(18,5,nStim);

% for motor-seeds (L/R)
P_tAvrcorr_2= nan(18,2,nStim);
P_tRescorr_2= nan(18,2,nStim);

ClusterIDs = [1,1];%GetClusterIDs('all');
% M_stimname = {'OMR'};
for i_stim = 1:length(M_stim)
    [M_stimrange,stimrange] = GetStimRange(M_stim{i_stim});
    
    for i_fish = range_fish
        stimrange = M_stimrange{i_fish};
        if ~isempty(stimrange)
            setappdata(hfig,'isMotorseed',0);
            LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange,0);
            stim = getappdata(hfig,'stim');
            behavior = getappdata(hfig,'behavior');
            
            
            
            [bh_tAvr,bh_tRes] = GetTrialAvrLongTrace(hfig,behavior);
            for i_lfr = 1:5
                vAvr_bh = var(bh_tAvr(i_lfr,:));
                vTot_bh = var(behavior(i_lfr,:));
                P_ratio(i_fish,i_lfr,1) = vAvr_bh/vTot_bh;
                P_tAvrcorr_5(i_fish,i_lfr,i_stim) = corr(behavior(i_lfr,:)',bh_tAvr(i_lfr,:)');
                P_tRescorr_5(i_fish,i_lfr,i_stim) = corr(behavior(i_lfr,:)',bh_tRes(i_lfr,:)');
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
                P_tAvrcorr_2(i_fish,i_lfr,i_stim) = corr(behavior(i_lfr,:)',bh_tAvr(i_lfr,:)');
                P_tRescorr_2(i_fish,i_lfr,i_stim) = corr(behavior(i_lfr,:)',bh_tRes(i_lfr,:)');
            end
        end
    end
end


%%
save(fullfile(outputDir,'varRatios_stimspecific.mat'),'P_tAvrcorr','P_tRescorr');

%% bar plot, pooling left and right sides
switch 1 % 1
    case 1
        P = P_tAvrcorr_2;
        xlabelname = 'frac. var. m.tAvr'; % i.e. 'corr(m,m.tAvr).^2';
    case 2
        P = P_tRescorr_2;
        xlabelname = 'frac. var. m.tRes'; % i.e. 'corr(m,m.tRes).^2';
end

figure('Position',[100,100,1050,150]);
for i_stim = 1:nStim

    data = P(range_fish,1:2,i_stim).^2;
    
    subplot(1,nStim,i_stim)
    hold on;
    histogram(data(:),0:0.1:1,'facecolor',[0.5,0.5,0.5],'edgecolor','w')
    xlabel(xlabelname)
    ylabel('count')
    box off
    set(gca,'Layer','top')
    title(M_stimname{i_stim})
end
