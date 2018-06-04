function [fMaskList_matchZ,fMask_matchZ_ZList] = ScreenAllfMasksWithzMasks(hfig)

% fMASKs.MaskDatabase = MaskDatabase;
%     fMASKs.mask_clusID = mask_clusID;
%     fMASKs.i_fish = i_fish;
%     fMASKs.note = 'Auto_defS_woA_Master0.5';

data_masterdir = GetCurrentDataDir();

range_fish =  1:4; % range_fish = GetFishRange();

%% custom params here:
MASKs = getappdata(hfig,'MASKs');
nzMasks = size(MASKs.MaskDatabase,2);
thres_ratio = 0.33;
thres_maxsize = 6*10^5; % ~ size of subpallium

%%
fMaskList_matchZ = cell(1,18);
fMask_matchZ_ZList = cell(1,18);
for i_fishnum = 1:length(range_fish),
    i_fish = range_fish(i_fishnum);
    disp(i_fish);
    
    %% Load all fMasks from this fish
    name = ['fMASKs_',num2str(i_fish),'.mat'];
    load(fullfile(data_masterdir,'fMASKs',name),'fMASKs');
    nfMasks = size(fMASKs.MaskDatabase,2);
    mask_clusID = fMASKs.mask_clusID;
    %%
    Grid_ratio = zeros(nzMasks,nfMasks);
    for i_zMask = 1:nzMasks,
       zMsk = MASKs.MaskDatabase(:,i_zMask); 
       if length(find(zMsk))<thres_maxsize,
           % compare this Z-Brain mask with each fMask
           for i_fMask = 1:nfMasks,
               [tf,ix] = ismember(i_fMask,mask_clusID);
               if tf,
                   fMsk = fMASKs.MaskDatabase(:,ix);
                   if length(find(fMsk))<thres_maxsize,
                       IX_intersect = find(zMsk.*fMsk);
                       r1 = length(IX_intersect)/length(find(zMsk));
                       r2 = length(IX_intersect)/length(find(fMsk));
                       Grid_ratio(i_zMask,i_fMask) = mean([r1,r2]);%max(r1,r2);
                   end
               end
               
           end
       end
    end
    
    % Pool from Grid
    [A,IX_zMask] = max(Grid_ratio,[],1);
    IX = find(A>thres_ratio);
    fMaskList_matchZ{i_fish} = IX;
    fMask_matchZ_ZList{i_fish} = IX_zMask(IX);
    
    %% save    
    save(fullfile(data_masterdir,'fMASKs','fMASK_matchZ_Output.mat'),...
        'fMaskList_matchZ','fMask_matchZ_ZList');
    timestamp  = datestr(now,'mmddyy_HHMM');
    backupname = ['fMASK_matchZ_Output_',timestamp,'.mat'];
    save(fullfile(data_masterdir,'fMASKs',backupname),...
        'fMaskList_matchZ','fMask_matchZ_ZList');    
end

end