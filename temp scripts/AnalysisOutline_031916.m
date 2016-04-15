%
range_fish = 8;

M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2];

%%
% ---- multi-fish analysis ----

range_stim = 1:12;


%% prep
% need to interpolate data to account for different scanning speeds???!!!!!!!

%% load 'data'
for i_fish = range_fish;

    f.LoadFullFish(hfig,i_fish);

%     global hm1;
%     hObject = hm1;
    
    periods = getappdata(hfig,'periods');
    setappdata(hfig,'stimrange',1:length(periods));

    setappdata(hfig,'cIX',1);
    setappdata(hfig,'gIX',1);
    UpdateTimeIndex(hfig);
    
    %% 
    stim = getappdata(hfig,'stim');
    fishset = getappdata(hfig,'fishset');
    regressors = GetStimRegressor(stim,fishset);
    Reg = zeros(length(regressors),length(stim));
    for i = range_stim,
        Reg(i,:) = regressors(i).im;
    end
    
    M_0 = getappdata(hfig,'M_0');
    Corr = corr(Reg',M_0');
    
    thres_reg = 0.3;
    [corr_max,IX] = max(Corr,[],1);
    cIX = find(corr_max>thres_reg)';
    gIX = IX(cIX)';
    numK = length(range_stim);
    
    f.UpdateIndices(hfig,cIX,gIX,numK);
    
    %%
    figure(hfig)
    % change fish number???
    f.pushbutton_loadCurrentClustersfromworkspace_Callback(hObject);
    
    %%
    M = getappdata(hfig,'M');
    CellXYZ = getappdata(hfig,'CellXYZ');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    absIX = getappdata(hfig,'absIX');
    cIX_abs = absIX(cIX);
    
    
    xyz = CellXYZ(cIX_abs,:);
    xyz_norm = CellXYZ_norm(cIX_abs,:);
    numStim = length(range_stim);
    xyz_avr = zeros(length(U),3);
    xyz_norm_avr = zeros(numStim,3);
    for i = 1:numStim,
        IX = find(gIX==range_stim(i));
        if ~isempty(IX),
            xyz_avr(i,:) = mean(xyz(IX,:),1);
            xyz_norm_avr(i,:) = mean(xyz_norm(IX,:),1);
        else
            xyz_avr(i,:) = NaN;
            xyz_norm_avr(i,:) = NaN;
        end
    end
    %%
    % Fish = cell(1,13);
    Fish{i_fish} = [];
    Fish{i_fish}.cIX = cIX;
    Fish{i_fish}.gIX = gIX;
    Fish{i_fish}.xyz = xyz;
    Fish{i_fish}.xyz_avr = xyz_avr;
    Fish{i_fish}.xyz_norm = xyz_norm;
    Fish{i_fish}.xyz_norm_avr = xyz_norm_avr;

end
%%
%% regression

% % get stim and motor regressors
% % motor: L/R/F, full length data
% % stim: use regressor and data for given stimlus range only
% % (use stimset or manually code matching reg across fish??)
% 
% % e.g.:
% % stim_range = [2,3];
% % regchoice = {1,stim_range};
% 
% 
% % stim = getappdata(hfig,'stim');
% 
% regchoice = {2,1};
% 
% if regchoice{1}==1, % stim Regressor
%     regressors = GetStimRegressor(stim,fishset);
%     if length(regchoice{2})>1,
%         regressor = zeros(length(regchoice{2}),length(regressors(1).im));
%         for i = 1:length(regchoice{2}),
%             regressor(i,:) = regressors(regchoice{2}(i)).im;
%         end
%     else
%         regressor = regressors(regchoice{2}).im;
%     end
%     
% elseif regchoice{1}==2, % motor Regressor
%     behavior = getappdata(hfig,'behavior');
%     regressors = GetMotorRegressor(behavior);
%     regressor = regressors(regchoice{2}).im;
% end
% 
% % regression
% thres_reg = 0.5;
% isCentroid = 1;
% [cIX,gIX,wIX] = f.Regression_Direct(hfig,thres_reg,regressor,isCentroid);

%%

% screen clusters for given regressor
% need to pass 2 thres: corr coeff and top percentage?? (need to hand-tune)

% -- between 2 fish --
% compare anatomical location:
% for given screened cluster, is there another screened cluster anatomically close
% (anatomically close defined as cluster centroids distance below
% threshold, and difference in distributedness below threshold)
% save selected pairs/multiple clusters
% output percentage of "conserved" clusters
% -- get average percentage for all pairs of fish --

% -- multi-fish clustering --
% or, plot clusters for all fish for given regressor onto standard brain,
% and cluster centroid location and distributedness; this "conserved"
% cluster needs to span multiple fish, ideally all fish

%% other clusters ((spontaneous))

% like above, find anatomically similar candidate clusters, and
% plot out their activity traces
% select by hand
% save selected pairs/multiple clusters

