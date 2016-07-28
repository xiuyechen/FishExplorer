function MakefMasksForAllFish(hfig)

% batch input data setup:
isFullData = 1;
data_masterdir = GetCurrentDataDir();

% M_fishset = GetFishStimset();
% M_stimrange = GetStimRange();

range_fish =  11:18; % range_fish = GetFishRange();

%% custom params here:
thres_minsize = 10000;
thres_maxsize = 6*10^5;

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    
        option
    [cIX,gIX,M_xyz_norm] = GetDefaultClustersFromLoad(hfig,i_fish,option);
    
    % Load fish
    LoadFullFish(hfig,i_fish,-1);    
    absIX = getappdata(hfig,'absIX');
    ClusterIDs = GetClusterIDs();
    [cIX,gIX] = LoadCluster_Direct(i_fish,ClusterIDs(1),ClusterIDs(2),absIX);
    anat_stack_norm = getappdata(hfig,'anat_stack_norm');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
    
    % make fMask for each auto-cluster
    fMASKs = [];
    MaskDatabase = sparse(0,0);
    mask_clusID = [];
    gIX = SqueezeGroupIX(gIX);
    for i_clus = 1:length(unique(gIX)),
        IX = find(gIX==i_clus);
        newMask = MakeFunctionalMask(anat_stack_norm,CellXYZ_norm,absIX,cIX(IX),'isBatch');
        sz = length(find(newMask));
        if sz>thres_minsize && sz<thres_maxsize,
            MaskDatabase(:, end+1) = newMask;
            mask_clusID = [mask_clusID;i_clus];
        end
    end
    fMASKs.MaskDatabase = MaskDatabase;
    fMASKs.mask_clusID = mask_clusID;
    fMASKs.i_fish = i_fish;
    fMASKs.note = 'Auto_defS_woA_Master0.5';
    % save:
    name = ['fMASKs_',num2str(i_fish),'.mat'];
    save(fullfile(data_masterdir,'fMASKs',name),'fMASKs');
end