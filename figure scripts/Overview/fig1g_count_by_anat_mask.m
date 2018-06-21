% fig2j. Anatomy screen: report cell number
% and find matching clusters

%% batch input data setup:
isFullData = 1;
data_masterdir = GetCurrentDataDir();

ClusterIDs = GetClusterIDs(); % Auto_M0.7_woA
const_ClusGroup = ClusterIDs(1);
const_Cluster = ClusterIDs(2); 
M_fishset = GetFishStimset();
M_stimrange = GetStimRange();

range_fish = 1:18; % range_fish = GetFishRange();

%% custom params here:
% load:
MASKs = load(fullfile(data_masterdir,'MaskDatabase.mat'));

% manually selected mask lists
anat_list1 = {275,1,94,114,77:93}; % gross divisions
anat_list1_names = {'Telencephalon','Diencephalon','Mesencephalon','Rhombencephalon','Ganglia'};
    
anat_list2 = {[5,35,74],[13,76],15,60,66,97,108,110,131,175,176,...
217,218,279,283,291}; % secondary anat features
anat_list2_names = {'Hypothalamus','Thalamus','Habenula','Preoptic Area','Pretectum','NucMLF','Tegmentum',...
    'Torus Semicircularis','Cerebellum','Inferior Olive','Interpeduncular Nucleus',...
    'Raphe (inferior)','Raphe (superior)','Olfactory Bulb','Pallium','Subpallium'};

anat_list3 = {[13,76],15,60,66,108,110,131,175,...
218,279,283,291}; % secondary anat features
anat_list3_names = {'Thalamus','Habenula','Preoptic Area','Pretectum','Tegmentum',...
    'Torus Semicircularis','Cerebellum','Inferior Olive',...
   'Raphe (superior)','Olfactory Bulb','Pallium','Subpallium'};

anat_list4 = {275,1,94,114,77:93,[5,35,74],[13,76],15,60,66,97,108,110,131,175,176,...
217,218,279,283,291}; % secondary anat features
anat_list4_names = {'Telencephalon','Diencephalon','Mesencephalon','Rhombencephalon','Ganglia','Hypothalamus','Thalamus','Habenula','Preoptic Area','Pretectum','NucMLF','Tegmentum',...
    'Torus Semicirc.','Cerebellum','Inferior Olive','IPN',...
    'Raphe (inferior)','Raphe (superior)','Olfactory Bulb','Pallium','Subpallium'}; % IPN = Interpeduncular Nucleus; 'Torus Semicircularis'

anat_listm = {201,117}; % collection of small masks with known matching clusters
anat_list_full = num2cell(1:294);

AnatScreen = cell(1,length(range_fish));

%% chose:
anat_list = anat_list4;
anat_list_names = anat_list4_names;

%%
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    disp(i_fish);
    
    % Load fish
    LoadFullFish(hfig,i_fish,isFullData);

    %% 1. 
    % setup
    absIX = getappdata(hfig,'absIX');
    fishset = M_fishset(i_fish);
    i_ClusGroup = const_ClusGroup;
    i_Cluster = const_Cluster;
    stimrange = M_stimrange{i_fish};
%     timelists = getappdata(hfig,'timelists');
    CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');    
    
    % Load cluster data
    [cIX,gIX] = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);
%     tIX = GetTimeIndex_Direct(stimrange,timelists,fishset);
%     M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);
    %M_0 = GetTimeIndexedData_Default_Direct(hfig,[],tIX,'isAllCells');
        
    %% ------custom code here---------          
    nClus = length(unique(gIX));
    nMasks = length(anat_list);
    Grid = zeros(nClus,length(anat_list),3); % matrix of all clusters x all masks
    for i_anat = 1:nMasks,
        Msk_IDs = anat_list{i_anat};
        [~,gIX_out] = ScreenCellsWithMasks(Msk_IDs,cIX,gIX,MASKs,CellXYZ_norm,absIX);
        % count number of cells per clus given mask
        Grid(:,i_anat,1) = histcounts(gIX_out,1:nClus+1);        
    end

    % normalize by mask: % for given map, cell% out of all cells in this mask
    temp = sum(Grid(:,:,1),1);
    knorm1 = repmat(temp,[nClus,1]);
    Grid(:,:,2) = Grid(:,:,1)./knorm1;
    
    % normalize by cluster: % for given cluster, cell% out of all cells in
    % this cluster
    temp = sum(Grid(:,:,1),2);
    knorm2 = repmat(temp,[1,nMasks]);
    Grid(:,:,3) = Grid(:,:,1)./knorm2;
        
    AnatScreen{i_fish} = Grid;
    
%     %%
% %     clusscore = max(Grid(:,:,2),[],2);    
% %     figure;imagesc(mean(Grid(:,:,2:3),3))
% 
%     figure;
%     subplot(1,2,1);imagesc(Grid(:,:,2));
%     title('non-prmsc mask')
%     subplot(1,2,2);imagesc(Grid(:,:,3))
%     title('non-prmsc cluster')
%     
%     %% save cluster
%     name = '?';
%     clusgroupID = 2;
%     clusIDoverride = 3;
%     SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
end

%%
data = zeros(length(range_fish),length(anat_list));
for i = 1:length(range_fish)
    i_fish = range_fish(i);
    Grid = AnatScreen{i_fish};
    data(i,:) = sum(Grid(:,:,1),1);    
end
%%
% figure;
% % imagesc(Grid(:,:,1));
% % Count = zeros(1,length(anat_list));
% % for i = 1:length(anat_list),
% %     Count(i) = sum(Grid(:,anat_list{i},1),1);
% % end
% bar(sum(Grid(:,:,1),1),'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1])
% title('# of cells by anatomical region')
% set(gca,'XTick',1:length(anat_list),'XTickLabel',anat_list_names,'XTickLabelRotation',45);
% xlim([0.5,length(anat_list)+0.5]);

figure;
% plot(mean(data,1));
errorbar(mean(data,1),std(data,1)/sqrt(length(range_fish)));

% bar(sum(Grid(:,:,1),1),'FaceColor',[0.5,0.5,0.5],'EdgeColor',[1,1,1])
title('# of cells by anatomical region')
set(gca,'XTick',1:length(anat_list),'XTickLabel',anat_list_names,'XTickLabelRotation',45);
xlim([0.5,length(anat_list)+0.5]);

%%


mask_range = 1:4;
D = data(:,mask_range);
[~,IX] = sort(mean(D,1),'descend');
D_sorted = D(:,IX);
anat_list_names_sorted = anat_list_names(mask_range(IX));

figure('Position',[500,500,200,250]); hold on

for i = 1:length(mask_range)
    scatter(ones(18,1)*i,D_sorted(:,i),10,[0.8,0.8,0.8],'filled');
    avr = mean(D_sorted(:,i));
    scatter(i,avr,10,[0.2,0.2,0.2],'filled');
    SEM = std(D_sorted(:,i))/sqrt(size(D_sorted,1));
    plot([i,i],[avr-SEM, avr+SEM],'color',[0.2,0.2,0.2])
end
% title('# of cells by anatomical region')
set(gca,'XTick',1:length(mask_range),'XTickLabel',anat_list_names_sorted,'XTickLabelRotation',45);
xlim([0.5,length(mask_range)+0.5]);
% ylim([0,600])
% semilogy


%%
mask_range = 6:21;
D = data(:,mask_range);
[~,IX] = sort(mean(D,1),'descend');
D_sorted = D(:,IX);
anat_list_names_sorted = anat_list_names(mask_range(IX));

figure('Position',[500,500,250,250]); hold on

plot_range = 1:9;

for i = plot_range
    scatter(ones(18,1)*i,D_sorted(:,i),10,[0.8,0.8,0.8],'filled');
    avr = mean(D_sorted(:,i));
    scatter(i,avr,10,[0.2,0.2,0.2],'filled');
    SEM = std(D_sorted(:,i))/sqrt(size(D_sorted,1));
    plot([i,i],[avr-SEM, avr+SEM],'color',[0.2,0.2,0.2])
end
% title('# of cells by anatomical region')
set(gca,'XTick',1:length(plot_range),'XTickLabel',anat_list_names_sorted(plot_range),'XTickLabelRotation',45);
xlim([0.5,length(plot_range)+0.5]);
% ylim([0,600])

figure('Position',[500,500,200,250]); hold on
plot_range = 10:16;

for i = 1:length(plot_range)
    i_mask = plot_range(i);
    scatter(ones(18,1)*i,D_sorted(:,i_mask),10,[0.8,0.8,0.8],'filled');
    avr = mean(D_sorted(:,i_mask));
    scatter(i,avr,10,[0.2,0.2,0.2],'filled');
    SEM = std(D_sorted(:,i_mask))/sqrt(size(D_sorted,1));
    plot([i,i],[avr-SEM, avr+SEM],'color',[0.2,0.2,0.2])
end
% title('# of cells by anatomical region')
set(gca,'XTick',1:length(plot_range),'XTickLabel',anat_list_names_sorted(plot_range),'XTickLabelRotation',45);
xlim([0.5,length(plot_range)+0.5]);
% ylim([0,600])


%%

% V-cell list:
% RoV3, MiV1 and MiV2: 232,193,194,

% Rh4-6: 222-224

% Mask 111: hypothalamic dopaminergic cluster (putative MLR)

% RS list: [126,127,128,184,187:194,226:232]

% 124: Rhombencephalon - Anterior Cluster of nV Trigeminal Motorneurons[Fish6]
% jaw 213: Rhombencephalon - Posterior Cluster of nV Trigeminal Motorneurons [Fish 13]
% jaw MiD3: 188???

% 187: Rhombencephalon - MiD2
% 188: Rhombencephalon - MiD3

% RS ref list: (Mauthner & homologues, and V-cells, 6 pairs total) 184,187,188,232,193,194

% 183: Rhombencephalon - Locus Coreuleus (fw motor seed)