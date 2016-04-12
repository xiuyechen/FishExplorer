range_fish = 8:11;

M_fish_set = [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2];

%%
% ---- multi-fish analysis ----

range_stimreg = 1:12;
range_motorreg = 1:3;
numReg = length(range_stimreg+range_motorreg);

FC = cell(14,numReg); % cell(numFish,numReg)

i_ClusGroup = 3;
i_Cluster = 1;

%% prep
% need to interpolate data to account for different scanning speeds???!!!!!!!

%% load 'data'
for i_fish = range_fish;
    
    f.LoadFullFish(hfig,i_fish);
    
    %     global hm1;
    %     hObject = hm1;
    
    periods = getappdata(hfig,'periods');
    if periods>1,
        setappdata(hfig,'stimrange',1:length(periods));
    else
        setappdata(hfig,'stimrange',1);
    end
    
    setappdata(hfig,'cIX',1);
    setappdata(hfig,'gIX',1);
    UpdateTimeIndex(hfig);
    
    %% Cluster indexing
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};

    absIX = getappdata(hfig,'absIX');
    cIX_abs = Cluster(i_Cluster).cIX_abs; % convert absolute index to index used for this dataset
    [~,cIX] = ismember(cIX_abs,absIX);
%     cIX = cIX(find(cIX));
    gIX = Cluster(i_Cluster).gIX;
    numK = Cluster(i_Cluster).numK;
%     setappdata(hfig,'cIX',cIX);
%     setappdata(hfig,'gIX',gIX);
%     setappdata(hfig,'numK',numK);
    f.UpdateIndices(hfig,cIX,gIX,numK);
    C_0 = FindCentroid(hfig);
    
    %% regression
    stim = getappdata(hfig,'stim');
    stim_regressors = GetStimRegressor(stim,fishset);
    
    behavior = getappdata(hfig,'behavior');
    motor_regressors = GetMotorRegressor(behavior);

    Reg = zeros(numReg,length(stim));
    
    for i = 1:length(range_stimreg),
        Reg(i,:) = stim_regressors(range_stimreg(i)).im;
    end
    for i = 1:length(range_motorreg),
        Reg(i+length(range_stimreg),:) = motor_regressors(range_motorreg(i)).im;
    end

    Corr = corr(Reg',C_0');
    
    for i_reg = 1:numReg,
        % set fixed threshold
        thres_reg = 0.3;
        [corr_max,IX] = max(Corr,[],1);
        clusIX_allreg = find(corr_max>thres_reg)';
        gIX_clus = IX(clusIX_allreg);
        
        % or find top match?
        
        
        clusIX = clusIX_allreg(gIX_clus==regID);
        %%
        [cIX,gIX,numK] = FindCellsFromCentroids(clusIX,cIX,gIX);
        f.UpdateIndices(hfig,cIX,gIX,numK);
        
        
        %     figure(hfig)
        %     % change fish number???
        %     f.pushbutton_loadCurrentClustersfromworkspace_Callback(hObject);
        
        %%
        M = getappdata(hfig,'M');
        CellXYZ = getappdata(hfig,'CellXYZ');
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        absIX = getappdata(hfig,'absIX');
        cIX_abs = absIX(cIX);
        
        
        xyz = CellXYZ(cIX_abs,:);
        xyz_norm = CellXYZ_norm(cIX_abs,:);
        numStim = length(range_motorreg);
        xyz_avr = zeros(length(U),3);
        xyz_norm_avr = zeros(numStim,3);
        for i = 1:numStim,
            IX = find(gIX==range_motorreg(i));
            if ~isempty(IX),
                xyz_avr(i,:) = mean(xyz(IX,:),1);
                xyz_norm_avr(i,:) = mean(xyz_norm(IX,:),1);
            else
                xyz_avr(i,:) = NaN;
                xyz_norm_avr(i,:) = NaN;
            end
        end
        %%
        FC{i_fish,i_reg} = [];
        FC{i_fish,i_reg}.cIX = cIX;
        FC{i_fish,i_reg}.gIX = gIX;
        FC{i_fish,i_reg}.xyz = xyz;
        FC{i_fish,i_reg}.xyz_avr = xyz_avr;
        FC{i_fish,i_reg}.xyz_norm = xyz_norm;
        FC{i_fish,i_reg}.xyz_norm_avr = xyz_norm_avr;
        
    end
end

datadir = GetCurrentDataDir();
save(fullfile(datadir,'FC.mat'),'FC');

% TtestTrial_loopall;
