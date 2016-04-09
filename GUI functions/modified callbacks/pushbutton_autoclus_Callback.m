function pushbutton_autoclus_Callback(hObject,f,i_fish)
%% test
% hfig = f.getParentFigure(hObject);
% cIX = getappdata(hfig,'cIX');
% gIX = getappdata(hfig,'gIX');
% SaveCluster_Direct(hfig,cIX(1:100),gIX(1:100),'test');

%%
tic
hfig = f.getParentFigure(hObject);
setappdata(hfig,'i_fish',i_fish);

cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');
M_0 = getappdata(hfig,'M_0');
thres_size = 10;
thres_split = getappdata(hfig,'thres_split');
thres_stimlock = 1.0;
thres_merge = getappdata(hfig,'thres_merge');
thres_silh = 0.4;

isWkmeans = getappdata(hfig,'isWkmeans');

%% kmeans
if isWkmeans,
    numK = 20;
    disp(['kmeans k = ' num2str(numK)]);
    tic
    rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK < 10^7 && numK~=1,
        disp('Replicates = 5');
        gIX = kmeans(M,numK,'distance','correlation','Replicates',5);
    elseif numel(M)*numK < 10^8 && numK~=1,
        disp('Replicates = 3');
        gIX = kmeans(M,numK,'distance','correlation','Replicates',3);
    else
        gIX = kmeans(M,numK,'distance','correlation');
    end
    toc
    SaveCluster_Direct(hfig,cIX,gIX,'k=20');
end
[gIX, numU] = f.SqueezeGroupIX(gIX);

%% pushbutton_iter_split(hObject,~);
disp('iter. split all...');
I_rest = [];
iter = 1;
gIX_last = gIX;
I_clean_last = cIX;
cIX = [];
gIX = [];
for i = 1:numU,
    disp(['i = ' num2str(i)]);
    ix = gIX_last == i;
    IX = I_clean_last(ix);
    M_s = M_0(IX,:);
    [I_rest,cIX,gIX,numU] = f.CleanClus(M_s,IX,I_rest,cIX,gIX,numU,1-thres_split,thres_size);
end
[gIX, ~] = f.SqueezeGroupIX(gIX);
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
% SaveCluster_Direct(hfig,cIX,gIX,['clean_round' num2str(iter)]);
% SaveCluster_Direct(hfig,I_rest,ones(length(I_rest),1),['rest_round' num2str(iter)]);

[gIX, numU] = f.Merge_direct(thres_merge,M_0,cIX,gIX);

%% rank by stim-lock
disp('stim-lock');
f.UpdateIndices(hfig,cIX,gIX,numU);
[gIX,rankscore] = f.RankByStimLock_Direct(hfig,gIX,numU);
disp('ranking complete');
% and threshold
IX = find(rankscore<thres_stimlock);
ix = ismember(gIX,IX);
gIX = gIX(ix);
cIX = cIX(ix);
f.UpdateIndices(hfig,cIX,gIX);

%% Regression with the centroid of each cluster
[cIX,gIX,~] = AllCentroidRegression_direct(hfig);
disp('auto-reg-clus complete');

[gIX, numU] = f.Merge_direct(thres_merge,M_0,cIX,gIX);
% SaveCluster_Direct(hfig,cIX,gIX,'clean_round2');

%% Silhouette
disp('silhouette analysis');
gIX_last = gIX;
for i = 1:numU,
    disp(['i = ' num2str(i)]);
    IX = find(gIX_last == i);
    cIX_2 = cIX(IX);
    M_s = M_0(cIX_2,:);
    % try k-means with k=2, see whether to keep
    gIX_ = kmeans(M_s,2,'distance','correlation');
    silh = silhouette(M_s,gIX_,'correlation');
    if mean(silh)>thres_silh,
        % keep the k-means k=2 subsplit
        disp('split');
        gIX(IX) = gIX_ + numU; % reassign (much larger) gIX
    end
end
[gIX, ~] = f.SqueezeGroupIX(gIX);

%% rank by stim-lock ?? bug?
% disp('stim-lock');
% M = M_0(cIX,:);
% [gIX,rankscore] = RankByStimLock_Direct(hfig,cIX,gIX,M,numU);
% disp('ranking complete');
% % and threshold
% IX = find(rankscore<thres_stimlock);
% ix = ismember(gIX,IX);
% gIX = gIX(ix);
% cIX = cIX(ix);
%
% [gIX, ~] = Merge_direct(thres_merge,M_0,cIX,gIX);

% size threshold
thres_size = getappdata(hfig,'thres_size');
[cIX, gIX, numU] = f.ThresSize(cIX,gIX,thres_size);

%% update GUI
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
f.UpdateIndices(hfig,cIX,gIX,numU);
% RefreshFigure(hfig);

SaveCluster_Direct(hfig,cIX,gIX,'Full_autoclus');
beep;
toc
end

function [Cluster,clusID] = SaveCluster_Direct(hfig,cIX,gIX,name) %,clusheader,name)
% new_clusgroupID = 1;
% UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID);

if ~exist('cIX','var'),
    cIX = getappdata(hfig,'cIX');
end
if ~exist('gIX','var'),
    gIX = getappdata(hfig,'gIX');
end

i_fish = getappdata(hfig,'i_fish');
absIX = getappdata(hfig,'absIX');
cIX_abs = absIX(cIX);

i_ClusGroup = 3;
global VAR;
if length(VAR(i_fish).ClusGroup)<i_ClusGroup,
    VAR(i_fish).ClusGroup = [VAR(i_fish).ClusGroup,cell(1,1)];
end
Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};

clusID = numel(Cluster)+1;
Cluster(clusID).name = name; %[clusheader name];
Cluster(clusID).cIX_abs = cIX_abs;
Cluster(clusID).gIX = gIX;
Cluster(clusID).numK = length(unique(gIX));

setappdata(hfig,'Cluster',Cluster);
% UpdateClusID(hfig,clusID);
VAR(i_fish).ClusGroup{i_ClusGroup} = Cluster;
disp('cluster saved');
end