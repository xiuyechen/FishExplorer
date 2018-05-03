
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

M_pair_range = {[2,1],[2,1],[12,11],[14,15],[0,3],[13,0]}; % [3,2] for PT L vs R; [9,8] for OMR L vs R
%     leftstates_full = [2,12,14,22];
%     rightstates_full = [1,11,15,21];
%     leftstates_full = [2,12,14,22];
%     rightstates_full = [1,11,15,21];
%     
% if fishset == 1,
%     States = [0,1,2,3];
%     singleNames = {'black','phototaxis R','phototaxis L','white'};%,...
% %         'L on','R on','L off','R off'};
% % elseif fishset == 2,
% %     States = [0,1,2,3,4,10,11,12];
% %     names = {'black','phototaxis left','phototaxis right','white','grey',...
% %         'OMR forward','OMR left','OMR right',...
% %         'left PT&OMR','right PT&OMR'};
% else%if fishset == 2,
%        States = [0,1,2,3,4,9,10,11,12,13,14,15,21,22,23];
%     singleNames = {'black','phototaxis R','phototaxis L','white','grey',...% 0-4
%         'OMR backward', 'OMR forward','OMR right','OMR left',... % 9-12
%         'Dot','looming L','looming R',... % 13-15
%         'red/blue R','red/blue L','red/red',...
%         };
% end
nStim = length(M_stim);
% for fictive channels (5 channels)
P_pval5= nan(18,5,nStim);
% for motor-seeds (L/R)
P_pval2= nan(18,2,nStim);

ClusterIDs = [1,1];%GetClusterIDs('all');
% M_stimname = {'OMR'};
for i_stim = 1:length(M_stim)
    [M_stimrange,stimrange] = GetStimRange(M_stim{i_stim});
    
    for i_fish = range_fish
        stimrange = M_stimrange{i_fish};
        
        tteststimrange = M_pair_range{i_stim};
        
        if ~isempty(stimrange)
            setappdata(hfig,'isMotorseed',0);
            LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange,0);
            stim = getappdata(hfig,'stim');
            behavior = getappdata(hfig,'behavior');

            %%
            samples = cell(1,length(tteststimrange));
            for i = 1:length(tteststimrange),
                IX = find(stim==tteststimrange(i));
                samples{i} = behavior(:,IX)';
            end
            
            [~, p] = ttest2(samples{1},samples{2});
            
            P_pval5(i_fish,:,i_stim) = p;
            
            %%
            setappdata(hfig,'isMotorseed',1);
            [~,~,behavior] = UpdateTimeIndex(hfig);
            
            samples = cell(1,length(tteststimrange));
            for i = 1:length(tteststimrange),
                IX = find(stim==tteststimrange(i));
                samples{i} = behavior(:,IX)';
            end
            
            [~, p] = ttest2(samples{1},samples{2});
            
            P_pval2(i_fish,:,i_stim) = p;
                        
        end
    end
end

%%
figure('Position',[550,100,800,120]);
for i_stim = 2:6 
    subplot(1,5,i_stim-1); hold on;
    y = P_pval2(:,:,i_stim);
    h = histogram(y(:),[0:0.05:1]);
    h.FaceColor = [0.5,0.5,0.5];
    xlim([0,1]);ylim([0,30])
    plot([0.05,0.05],[0,30],'r--')
    title(M_stimname{i_stim})
    length(find(y(:)>1))
end