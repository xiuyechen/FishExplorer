% load motor seeds from VAR, save into file and save into GUIdata for each
% fish % careful!!! overwrite!!! 
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
range_fish = [1:12,14:18];%GetFishRange;%[1:3,5:18];
for i_fish = range_fish,
    %% load
    save_masterdir = GetCurrentDataDir();
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    % load 'data'
    load(fullfile(save_dir,'data_full.mat'),'data'); % struct with many fields
%     save(fullfile(save_dir,'data_full_050617.mat'),'data');
    
    %% get full-length motorseed regressor
%     ClusterIDs = [11,1];
%     [~,stimrange] = GetStimRange('F',i_fish); % full range
%     setappdata(hfig,'isMotorseed',0);
%     [cIX,gIX,M] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);      
%     Behavior_full_motorseed = FindClustermeans(gIX,M);
    % need to do this in raw time instead! (note that full range is not the
    % same sequence as CellResp for fishset=2)
    LoadFullFish(hfig,i_fish);
    ClusterIDs = [12,1];
    [cIX,gIX] = LoadCluster_Direct(i_fish,ClusterIDs(1),ClusterIDs(2));
    
    CellResp = getappdata(hfig,'CellResp');
    CellRespAvr = getappdata(hfig,'CellRespAvr');
    
    Behavior_full_motorseed = FindClustermeans(gIX,CellResp(cIX,:));
    BehaviorAvr_motorseed = FindClustermeans(gIX,CellRespAvr(cIX,:));
    
    %% compute stim-rep averaged motor regressor
    
%     if length(data.timelists_names)==1, % fish with single stimulus set
%         period = data.periods;
%         BehaviorAvr_motorseed = mean(reshape(Behavior_full_motorseed,size(Behavior_full_motorseed,1),period,[]),3);
%         
%     else % fish with multiple stimulus sets
%         nTypes = length(data.timelists_names);
%         periods = data.periods;
%         BehaviorAvr_motorseed = [];
%         for i = 1:nTypes,
%             M_behav = Behavior_full_motorseed(:,data.timelists{i});
%             if mod(numel(M_behav),size(M_behav,1)*periods(i))==0,            
%                 avr = mean(reshape(M_behav,size(M_behav,1),periods(i),[]),3);
%             else % for patterns with dummy periods where the stimulus is constant, like 'spont'
%                 nrep = floor(numel(M_behav)/(size(M_behav,1)*periods(i)));
%                 M2_behav = M_behav(:,1:nrep*periods(i));
%                 avr = mean(reshape(M2_behav,size(M_behav,1),periods(i),[]),3);
%             end
%             BehaviorAvr_motorseed = horzcat(BehaviorAvr_motorseed,avr); %#ok<AGROW>
%         end
%     end
    
    %% save into GUIdata % careful!!! overwrite!!!   
    data.Eye_full = Behavior_full_motorseed;
    data.Eye_avr = BehaviorAvr_motorseed;
    save(fullfile(save_dir,'data_full.mat'),'data');
end