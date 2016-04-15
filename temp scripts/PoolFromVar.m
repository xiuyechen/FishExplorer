
save_masterdir = GetCurrentDataDir();
%%
C = [];

range_fish = [8,9,11];
M_ClusGroup = [1,2,2];
M_Cluster = [7,5,5];

for i = 1:length(range_fish),
    i_fish = range_fish(i)
    i_ClusGroup = M_ClusGroup(i);
    i_Cluster = M_Cluster(i);
    cIX_abs = VAR(i_fish).ClusGroup{i_ClusGroup}(i_Cluster).cIX_abs;
    
    %%
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    %%
    filename = fullfile(save_dir,'TimeSeries_full.h5');
    tic
    CellResp = h5read(filename,'/CellResp');
    CellRespZ = h5read(filename,'/CellRespZ');
    CellRespAvr = h5read(filename,'/CellRespAvr');
    CellRespAvrZ = h5read(filename,'/CellRespAvrZ');
    toc
    %%
    C(i_fish).cIX_abs = cIX_abs;
    C(i_fish).CellResp = CellResp(cIX_abs,:);
    C(i_fish).CellRespZ = CellRespZ(cIX_abs,:);
    C(i_fish).CellRespAvr = CellRespAvr(cIX_abs,:);
    C(i_fish).CellRespAvrZ = CellRespAvrZ(cIX_abs,:);
    
    %%
    load(fullfile(save_dir,'data_full.mat'));
    names = {'stim_full','stimAvr','BehaviorAvr','Behavior_full',...
        'CellXYZ_norm','timelists','timelists_names'};
    for i = 1:length(names),
        eval(['C(i_fish).',names{i},' = data.',names{i},';']);
    end
   
end
%%
tic
timestamp  = datestr(now,'mmddyy_HHMM');
save(fullfile(save_masterdir,['C_' timestamp '.mat']),'C');
toc
