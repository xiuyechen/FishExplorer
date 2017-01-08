%% Regression screen all clus all fish

range_fish = 8:18;
Coeff_Allfish = cell(1,length(range_fish));
Coeff_names = cell(1,9);
%%
for i_fishnum = 1:length(range_fish),
    i_fish = range_fish(i_fishnum);
    % i_fish = 8;
    disp(i_fish);
    
    isFullData = true;
    LoadFullFish(hfig,i_fish,isFullData);
    timelists = getappdata(hfig,'timelists');
    %%
    i_ClusGroup = 6;
    i_Cluster = 1;
    
    [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster);
    
    numK = length(unique(gIX));
    Coeff = zeros(numK,9);
    %%
    for i_stim = 1:3,
        if i_stim == 1,
            M_stimrange = GetStimRange('P'); % Phototaxis
            regrange = 1:3;
        elseif i_stim == 2,
            M_stimrange = GetStimRange('O'); % OMR
            regrange = 4:6;
            %     elseif i_stim == 3,
            %         M_stimrange = GetStimRange('S'); % Spont
        elseif i_stim == 3,
            M_stimrange = GetStimRange('3'); % the above three
            regrange = 7:9;
        end
        
        stimrange = M_stimrange{i_fish};
        
        fishset = getappdata(hfig,'fishset');
        
        % Load cluster data
        tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
        [M,behavior,stim] = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
        
        Data = FindCentroid_Direct(gIX,M);
        
        % get stim/motor regressors
        if i_stim == 1,
            [~,stimregnames,Reg_raw] = GetStimRegressor(stim,fishset,i_fish);

            IX = find(sum(Reg_raw,2));
            Reg = Reg_raw(IX,:);
            Coeff_names(regrange) = stimregnames(IX);
        elseif i_stim == 2,
            [~,stimregnames,Reg_raw] = GetStimRegressor(stim,fishset,i_fish);            
            IX = [7;8;9]; % hardcode!! for OMR f/r/l, ignore baseline and backward grating
%             IX = find(sum(Reg_raw,2));
            Reg = Reg_raw(IX,:);
            Coeff_names(regrange) = stimregnames(IX);
        elseif i_stim == 3,
            [~,~,Reg,motorregnames] = GetMotorRegressor(behavior,i_fish);
            Coeff_names(regrange) = motorregnames;
        end
        
        % regression
        Coeff(:,regrange) = corr(Reg',Data')';
    end
    Coeff_Allfish{i_fishnum} = Coeff;
end
data_dir = GetCurrentDataDir();
save(fullfile(data_dir,'RegsCoeff.mat'),'Coeff','Coeff_names');

%% Multi-fish clustering
nTotalClus = 0;
clus2fish = [];
M_coeffs = [];
for i_fishnum = 1:length(range_fish),
    i_fish = range_fish(i_fishnum);
    Coeff = Coeff_Allfish{i_fishnum};
    numK = size(Coeff,1);
    nTotalClus = nTotalClus+numK;
    clus2fish = vertcat(clus2fish,ones(numK,1)*i_fish);
    M_coeffs = vertcat(M_coeffs,Coeff);
end
%%
numK = 30;
[gIX,C] = kmeans(M_coeffs,numK);
% figure;
% [~,IX] = sort(gIX);
% imagesc(M_coeffs(IX,:))
%%
kmeansPlot(M_coeffs,gIX,C)
%%
clus_multifishcount = zeros(numK,1);
for i = 1:numK,
    temp = clus2fish(gIX==i);
    clus_multifishcount(i) = length(unique(temp));
end
figure;hist(clus_multifishcount,1:11)

