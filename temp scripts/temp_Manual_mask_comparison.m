
i_fish = 1;
LoadFullFish(hfig,i_fish,-1);

absIX = getappdata(hfig,'absIX');
fishset = M_fishset(i_fish);
i_ClusGroup = 1;
i_Cluster = 2;
stimrange = M_stimrange{i_fish};
timelists = getappdata(hfig,'timelists');

% Load cluster data
cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);

anat_stack_norm = getappdata(hfig,'anat_stack_norm');
CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
m1 = MakeFunctionalMask(anat_stack_norm,CellXYZ_norm,absIX,cIX);

%%
i_fish = 11;
LoadFullFish(hfig,i_fish,-1);

absIX = getappdata(hfig,'absIX');
fishset = M_fishset(i_fish);
i_ClusGroup = 1;
i_Cluster = 2;
stimrange = M_stimrange{i_fish};
timelists = getappdata(hfig,'timelists');

% Load cluster data
cIX = LoadCluster_Direct(i_fish,i_ClusGroup,i_Cluster,absIX);

anat_stack_norm = getappdata(hfig,'anat_stack_norm');
CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
m2 = MakeFunctionalMask(anat_stack_norm,CellXYZ_norm,absIX,cIX);

%% 
% m1 = MASKs.MaskDatabase(:,219);
% m2 = MASKs.MaskDatabase(:,247);
tic
m0 = m1.*m2;
mIX = find(m0);
r1 = length(mIX)/length(find(m1))
r2 = length(mIX)/length(find(m2))
toc
