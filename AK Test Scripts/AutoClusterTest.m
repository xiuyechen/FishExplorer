function AutoClusterTest(hObject,~)
%hfig = getParentFigure(hObject);
hfig = hObject;% pass it hFig for now

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
[gIX, numU] = SqueezeGroupIX(gIX);

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
    [I_rest,cIX,gIX,numU] = CleanClus(M_s,IX,I_rest,cIX,gIX,numU,1-thres_split,thres_size);
end
[gIX, ~] = SqueezeGroupIX(gIX);
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
SaveCluster_Direct(hfig,cIX,gIX,['clean_round' num2str(iter)]);
SaveCluster_Direct(hfig,I_rest,ones(length(I_rest),1),['rest_round' num2str(iter)]);

[gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX);

%% rank by stim-lock
disp('stim-lock');
UpdateIndices(hfig,cIX,gIX,numU);
[gIX,rankscore] = RankByStimLock_Direct(hfig,gIX,numU);
disp('ranking complete');
% and threshold
IX = find(rankscore<thres_stimlock);
ix = ismember(gIX,IX);
gIX = gIX(ix);
cIX = cIX(ix);
UpdateIndices(hfig,cIX,gIX);

%% Regression with the centroid of each cluster
[cIX,gIX,~] = AllCentroidRegression_direct(hfig);
disp('auto-reg-clus complete');

[gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX);
SaveCluster_Direct(hfig,cIX,gIX,'clean_round2');

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
[gIX, ~] = SqueezeGroupIX(gIX);

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
[cIX, gIX, numU] = ThresSize(cIX,gIX,thres_size);

%% update GUI
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);

SaveCluster_Direct(hfig,cIX,gIX,'clean_round3');
beep;

end

%% Functions
function pushbutton_thressize_Callback(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
thres_size = getappdata(hfig,'thres_size');
[cIX, gIX, numU] = ThresSize(cIX,gIX,thres_size);
UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
end

function edit_sizethres_Callback(hObject,~)
str = get(hObject,'String');
if ~isempty(str),
    temp = textscan(str,'%f');
    thres_size = temp{:};
end
hfig = getParentFigure(hObject);
setappdata(hfig,'thres_size',thres_size);
end

function [cIX, gIX, numU] = ThresSize(cIX,gIX,thres_size)
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_size,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX, numU] = SqueezeGroupIX(gIX);
end

function [gIX, numU] = Merge_direct(thres_merge,M_0,cIX,gIX)
M = M_0(cIX,:);
[gIX, numU] = HierClus(M,gIX);
U = unique(gIX);
M = M_0(cIX,:);
[C,D] = FindCentroid_Direct(gIX,M);
i = 1;
while i<numU,
    c = corr(C(i,:)',C(i+1,:)');
    if c > thres_merge,
        IX = find(gIX == U(i+1));
        gIX(IX)=U(i); %#ok<*FNDSB>
        U = unique(gIX);
        numU = length(U);
        
        IX = find(gIX == U(i));
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
        C(i+1,:) = [];
        D(i+1) = [];
    else
        i = i+1;
    end
end
[gIX, numU] = HierClus(M,gIX);
disp('merging complete');
end

function checkbox_wkmeans_Callback(hObject,~)
hfig = getParentFigure(hObject);
setappdata(hfig,'isWkmeans',get(hObject,'Value'));
end
function pushbutton_merge_Callback(hObject,~)
% disp('merging...');
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');
U = unique(gIX);
numU = length(U);
[C,D] = FindCentroid(hfig);

thres_merge = getappdata(hfig,'thres_merge');

i = 1;
while i<numU,
    c = corr(C(i,:)',C(i+1,:)');
    if c > thres_merge,
        IX = find(gIX == U(i+1));
        gIX(IX)=U(i);
        U = unique(gIX);
        numU = length(U);
        
        IX = find(gIX == U(i));
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
        C(i+1,:) = [];
        D(i+1) = [];
    else
        i = i+1;
    end
end

if numU>1,
    [gIX, numU] = HierClus(M,gIX);
end

UpdateIndices(hfig,cIX,gIX,numU);
RefreshFigure(hfig);
disp('merging complete');
end

function edit_mergethres_Callback(hObject,~)
str = get(hObject,'String');
temp = textscan(str,'%f',1);
hfig = getParentFigure(hObject);
setappdata(hfig,'thres_merge',temp{:});
end

function pushbutton_iter_split(hObject,~)
hfig = getParentFigure(hObject);
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
numU = getappdata(hfig,'numK');
thres_split = getappdata(hfig,'thres_split');
M_0 = getappdata(hfig,'M_0');

disp('iter. split all, beep when done...');
thres_size = 10;
thres_H = thres_split;
% thres_H = [0.2;0.15;0.1;0.05]; % could have more rounds...

% initialization
I_rest = [];
% loop
tic
for round = 1:length(thres_H),
    disp(['round ' num2str(round) ', numU = ' num2str(numU)]);
    dthres = 1-thres_H(round);
    gIX_last = gIX;
    I_clean_last = cIX;
    cIX = [];
    gIX = [];
    for i = 1:numU,
        disp(['i = ' num2str(i)]);
        ix = find(gIX_last == i);
        %         IX = ix;
        IX = I_clean_last(ix);
        M_s = M_0(IX,:);
        [I_rest,cIX,gIX,numU] = CleanClus(M_s,IX,I_rest,cIX,gIX,numU,dthres,thres_size);
        %         cIX = I_clean_last(I_clean);
    end
    
    [gIX, numU] = SqueezeGroupIX(gIX);
    SaveCluster_Direct(hfig,cIX,gIX,['clean_round' num2str(round)]);
    
    SaveCluster_Direct(hfig,I_rest,ones(length(I_rest),1),['rest_round' num2str(round)]);
end
toc
beep
end

function [I_rest,I_clean,gIX_clean,numU] = CleanClus(M_s,IX,I_rest,I_clean,gIX_clean,numU,dthres,thres_size)
I_clean_s = [];

% find numK_s for kmeans
kmax = min(round(size(M_s,1)/thres_size),30);
% try numK_s = 1
numK_s = 1;
rng('default');
[gIX_s,~,~,D] = kmeans(M_s,numK_s,'distance','correlation');
Dist = min(D,[],2);
if mean(Dist)>dthres,
    % try numK_s = kmax
    numK_s = kmax;
    rng('default');
    [gIX_s,~,~,D] = kmeans(M_s,numK_s,'distance','correlation');
    Dist = min(D,[],2);
    if mean(Dist)<dthres, % find in between value for numK_s
        numK_s = 2;
        while 1,
            % kmeans-cluster by numK_s
            rng('default');
            [gIX_s,~,~,D] = kmeans(M_s,numK_s,'distance','correlation');
            Dist = min(D,[],2);
            if mean(Dist)<dthres,
                break;
            end
            
            if numK_s < kmax,
                numK_s = numK_s+1;
                disp(['numK_s = ' num2str(numK_s)]);
            else % numK_s = kmax;
                break;
            end
        end
    else disp(['numK_s = ' num2str(numK_s)]);
    end
end
% have numK_s that makes mean(Dist) < thres, or numK_s = kmax

for i = 1:numK_s,
    IX_s = find(gIX_s == i);
    if length(IX_s)>thres_size, % cluster big enough to start
        dst = Dist(IX_s);
        if mean(dst) < dthres,
            I_clean_s = [I_clean_s; IX_s]; %#ok<*AGROW>
            gIX_clean = [gIX_clean; gIX_s(IX_s)+double(numU)];
        else
            ix = find(dst<dthres); % clean
            if length(ix)>=thres_size, % clean cluster still big enough
                I = IX_s(ix);
                I_clean_s = [I_clean_s; I];
                gIX_clean = [gIX_clean; gIX_s(I)+double(numU)];
            end
        end
    end
end

numU = numU + numK_s;
I_rest_s = setdiff(1:size(M_s,1),I_clean_s);
I_rest = [I_rest; IX(I_rest_s)];
I_clean = [I_clean; IX(I_clean_s)];
end

function edit_splitthres_Callback(hObject,~)
str = get(hObject,'String');
temp = textscan(str,'%f',1);
hfig = getParentFigure(hObject);
setappdata(hfig,'thres_split',temp{:});
end


%% Internal functions

function UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID,norefresh) %#ok<INUSD>
% save/update old Cluster into ClusGroup before exiting,
% as Cluster is the variable handled in hfig but not saved elsewhere
ClusGroup = getappdata(hfig,'ClusGroup');
Cluster = getappdata(hfig,'Cluster');
i_fish = getappdata(hfig,'i_fish');

% update into workspace
ClusGroup{clusgroupID} = Cluster;
setappdata(hfig,'ClusGroup',ClusGroup);
global VAR;
VAR(i_fish).ClusGroup = CurrentClusGroup(hfig);

% load new 'Cluster'
Cluster = ClusGroup{new_clusgroupID};
setappdata(hfig,'Cluster',Cluster);
setappdata(hfig,'clusgroupID',new_clusgroupID);

% update GUI: hclusgroupmenu
global hclusgroupmenu hclusgroupname;
if ishandle(hclusgroupmenu),
    menu = MakeNumberedMenu(VAR(i_fish).ClusGroupName);
    set(hclusgroupmenu,'String',menu,'Value',new_clusgroupID+1);
    set(hclusgroupname,'String',VAR(i_fish).ClusGroupName(new_clusgroupID));
end

if ~exist('norefresh','var'),
    if numel(Cluster) == 0, % i.e. for newly created ClusGroup
        SaveCluster(hfig,'new');
    else % load this ClusGroup
        clusID = 1;
        UpdateClusID(hfig,clusID);
    end
end
end

function ClusGroup = CurrentClusGroup(hfig)
ClusGroup = getappdata(hfig,'ClusGroup');
Cluster = getappdata(hfig,'Cluster');
clusgroupID = getappdata(hfig,'clusgroupID');
ClusGroup{clusgroupID} = Cluster;
setappdata(hfig,'ClusGroup',ClusGroup);
end

function [Cluster,clusID] = SaveCluster(hfig,state,clusheader,name)
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
absIX = getappdata(hfig,'absIX');
Cluster = getappdata(hfig,'Cluster');

if strcmp(state,'current'),
    clusID = getappdata(hfig,'clusID');
else %if strcmp(state,'new'),
    clusID = numel(Cluster)+1;
    if ~exist('name','var'),
        name = getappdata(hfig,'newclusname');
    end
    if ~exist('clusheader','var'),
        clusheader = getappdata(hfig,'clusheader');
    end
    Cluster(clusID).name = [clusheader name];
end

cIX_abs = absIX(cIX);

Cluster(clusID).cIX_abs = cIX_abs;
Cluster(clusID).gIX = gIX;
Cluster(clusID).numel = length(cIX_abs);
U = unique(gIX);
numU = length(U);
Cluster(clusID).numK = numU;

setappdata(hfig,'Cluster',Cluster);
UpdateClusID(hfig,clusID);
disp('cluster saved');
end

function [Cluster,clusID] = SaveCluster_Direct(hfig,cIX,gIX,name) %,clusheader,name)
new_clusgroupID = 1;
clusgroupID = getappdata(hfig,'clusgroupID');
UpdateClusGroupID(hfig,clusgroupID,new_clusgroupID);

if ~exist('cIX','var'),
    cIX = getappdata(hfig,'cIX');
end
if ~exist('gIX','var'),
    gIX = getappdata(hfig,'gIX');
end

Cluster = getappdata(hfig,'Cluster');
absIX = getappdata(hfig,'absIX');
cIX_abs = absIX(cIX);

clusID = numel(Cluster)+1;
% if ~exist('name','var'),
%     name = getappdata(hfig,'newclusname');
% end
% if ~exist('clusheader','var'),
%     clusheader = getappdata(hfig,'clusheader');
% end
Cluster(clusID).name = name; %[clusheader name];
Cluster(clusID).cIX_abs = cIX_abs;
Cluster(clusID).gIX = gIX;
Cluster(clusID).numel = length(cIX);
Cluster(clusID).numK = length(unique(gIX));

setappdata(hfig,'Cluster',Cluster);
UpdateClusID(hfig,clusID);
disp('cluster saved');
end

function UpdateClusID(hfig,clusID)
Cluster = getappdata(hfig,'Cluster');

% save
setappdata(hfig,'clusID',clusID);
% update GUI
global hclusname hclusmenu;
set(hclusname,'String',Cluster(clusID).name);
menu = MakeNumberedMenu({Cluster.name});
set(hclusmenu,'String', menu,'Value',clusID+1);
numK = Cluster(clusID).numK;
gIX = Cluster(clusID).gIX;

% convert absolute index to index used for this dataset
absIX = getappdata(hfig,'absIX');
cIX_abs = Cluster(clusID).cIX_abs;
[~,cIX] = ismember(cIX_abs,absIX);

if ~isempty(find(cIX==0,1)),
    errordlg('cell index out of bound for currently loaded dataset');
    IX = cIX==0;
    cIX(IX) = [];
    gIX(IX) = [];
end

UpdateIndices(hfig,cIX,gIX,numK);
RefreshFigure(hfig);
end

function menu = MakeNumberedMenu(name) % e.g. name = {Cluster.name} (note {})
menu = [{'(choose)'},name];
for j=2:length(menu),menu(j)={[num2str(j-1) ': ' menu{j}]};end
end

function [gIX, numU] = HierClus(M,gIX,isplotfig) %#ok<INUSD>
[gIX, numU] = SqueezeGroupIX(gIX);
[C,~] = FindCentroid_Direct(gIX,M);
D = pdist(C,'correlation');
tree = linkage(C,'average','correlation');
leafOrder = optimalleaforder(tree,D);

if numU>1,
    if exist('isplotfig','var'),
        figure('Position',[100 100 600 600]);
        %             subplot(1,3,1);
        %             CORR = corr(C');
        %             CorrPlot(CORR);
        %
        %             subplot(1,3,2);
        dendrogram(tree,numU,'orientation','right','reorder',leafOrder);
        set(gca,'YDir','reverse');
        set(gca,'XTick',[]);
        
        %             subplot(1,3,3);
        %             C2 = C(leafOrder,:);
        %             CORR2 = corr(C2');
        %             CorrPlot(CORR2);
    end
    % sort for uniform colorscale
    temp = zeros(size(gIX));
    for i = 1:numU,
        temp(gIX==leafOrder(i)) = i; % = T(i) for clusters segmented from tree
    end
    gIX = temp;
end
end

function gIX = HierClusDirect(C,gIX,numU)
D = pdist(C,'correlation');
tree = linkage(C,'average','correlation');
leafOrder = optimalleaforder(tree,D);

% sort for uniform colorscale
temp = zeros(size(gIX));
for i = 1:numU,
    temp(gIX==leafOrder(i)) = i; % = T(i) for clusters segmented from tree
end
gIX = temp;
end

function [gIX, numK] = SqueezeGroupIX(gIX)
U = unique(gIX);
numK = length(U);
for i = 1:numK,
    old = U(i);
    gIX(gIX==old) = i;
end
end

% frequently used, updates cell-index,group-index,cluster-number. set-operations included in here.
function UpdateIndices(hfig,cIX,gIX,numK)
global hback hopID;
if ~exist('gIX','var'),
    gIX = getappdata(hfig,'gIX');
end
if ~exist('cIX','var'),
    cIX = getappdata(hfig,'cIX');
end

% update cache
bC = getappdata(hfig,'bCache');
cIX_last = getappdata(hfig,'cIX');
gIX_last = getappdata(hfig,'gIX');
if ~(isequal(cIX_last,cIX) && isequal(gIX_last,gIX)),
    bC.cIX = [cIX_last,bC.cIX];
    bC.gIX = [gIX_last,bC.gIX];
    bC.numK = [getappdata(hfig,'numK'),bC.numK];
    set(hback,'enable','on');
    if length(bC.cIX)>20,
        bC.cIX(end) = [];
        bC.gIX(end) = [];
        bC.numK(end) = [];
    end
end

% set operations, if applicable
opID = getappdata(hfig,'opID');
if opID ~= 0,
    switch opID,
        case 1,
            disp('union');
            [~,ia,ib] = union(cIX_last,cIX,'stable');
            IX = vertcat(cIX_last(ia),cIX(ib));% IX;
        case 2,
            disp('intersect');
            [IX,ia,~] = intersect(cIX_last,cIX);
            ib = [];
        case 3,
            disp('setdiff');
            [IX,ia] = setdiff(cIX_last,cIX);
            ib = [];
        case 4,
            disp('rev setdiff');
            % swap sequence, then same as opID==3
            temp = cIX;
            cIX = cIX_last;
            cIX_last = temp;
            temp = gIX;
            gIX = gIX_last;
            gIX_last = temp;
            [IX,ia] = setdiff(cIX_last,cIX);
            ib = [];
        case 5,
            disp('setxor');
            [IX,ia,ib] = setxor(cIX_last,cIX);
        case 6,
            disp('smartUnion');
            CIX = vertcat(cIX_last,cIX);
            GIX = [gIX_last;gIX+max(gIX_last)]; % gIX to match
            M_0 = getappdata(hfig,'M_0');
            [cIX,gIX,numK] = SmartUnique(CIX,GIX,M_0(CIX,:));
    end
    if opID<6,
        if ~isempty(IX),
            cIX = IX;
            gIX = vertcat(gIX_last(ia),gIX(ib)+max(gIX_last(ia)));
            numK = length(unique(gIX));
            %         [gIX, numK] = SqueezeGroupIX(gIX);
        else
            errordlg('operation result is empty set!')
            waitforbuttonpress;
        end
    end
    set(hopID,'Value',1,'BackgroundColor',[1,1,1]); % reset
    setappdata(hfig,'opID',0);
end

setappdata(hfig,'bCache',bC);
setappdata(hfig,'cIX',cIX);
setappdata(hfig,'gIX',gIX);

M = GetTimeIndexedData(hfig);
setappdata(hfig,'M',M);

if exist('numK','var'),
    setappdata(hfig,'numK',double(numK));
end

%% Resets: reset flags the NEXT time this function is called (so they only apply to this particular plot)
% handle rankID: >=2 means write numbers as text next to colorbar
% first UpdateIndices sets rankID to 100, second sets back to 0
rankID = getappdata(hfig,'rankID');
if rankID>=2,
    if rankID==100,
        setappdata(hfig,'rankID',0);
    else
        setappdata(hfig,'rankID',100);
    end
end

% toggle 'isWeighAlpha'
isWeighAlpha = getappdata(hfig,'isWeighAlpha');
if isWeighAlpha == 1,
    setappdata(hfig,'isWeighAlpha',100);
elseif isWeighAlpha == 100,
    setappdata(hfig,'isWeighAlpha',0);
end

% FindCentroid reset:
setappdata(hfig,'Centroids',[]);
setappdata(hfig,'D_ctrd',[]);
end

% frequently used, 2 plotting functions are outside ('DrawTimeSeries.m' and 'DrawCellsOnAnatProj.m')
function RefreshFigure(hfig)
watchon; drawnow;
isPopout = 0; % with down-sampling in plots

% clean-up canvas
allAxesInFigure = findall(hfig,'type','axes');
if ~isempty(allAxesInFigure)
    delete(allAxesInFigure);
end

figure(hfig);
h1 = axes('Position',[0.05, 0.04, 0.55, 0.83]);
h2 = axes('Position',[0.63, 0.04, 0.35, 0.83]);

isCentroid = getappdata(hfig,'isCentroid');
isRefAnat = getappdata(hfig,'isRefAnat');

isPlotLines = 0; %getappdata(hfig,'isPlotLines');
isPlotBehavior = 1; %getappdata(hfig,'isPlotBehavior');

% double-check if cIX is valid
cIX = getappdata(hfig,'cIX');
if isempty(cIX),
    errordlg('empty set!');
    % GO BACK to the last step (presumably not empty)
    pushbutton_back_Callback(h1); % using h1 instaed of the usual 'hObject'
    return;
end

% left subplot
axes(h1);
DrawTimeSeries(hfig,h1,isPopout,isCentroid,isPlotLines,isPlotBehavior);

% right subplot
axes(h2);
DrawCellsOnAnatProj(hfig,isRefAnat,isPopout);
watchoff;
end

function [C,D] = FindCentroid_Direct(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
    end
end
end

function UpdateTimeIndex(hfig,isSkipcIX) %#ok<INUSD>
% input params
isAvr = getappdata(hfig,'isAvr');
isRawtime = getappdata(hfig,'isRawtime');
stimrange = getappdata(hfig,'stimrange');
% load
timelists = getappdata(hfig,'timelists');
periods = getappdata(hfig,'periods');
fishset = getappdata(hfig,'fishset');

if fishset == 1,
    if isAvr,
        tIX = 1:periods;
    else
        tIX = timelists{1};
    end
    
else % fishset>1,
    if isAvr,
        tIX = [];
        for i = 1:length(stimrange),
            ix = stimrange(i);
            i_start = sum(periods(1:ix-1)); % if ix-1<1, sum = 0
            tIX = horzcat(tIX,(i_start+1:i_start+periods(ix)));
            %             tIX = vertcat(tIX,(i_start+1:i_start+periods(ix))');
        end
    else % full range
        if ~isRawtime,
            tIX = cat(2, timelists{stimrange});
        else
            tIX = sort(cat(2, timelists{stimrange}));
        end
    end
end

setappdata(hfig,'tIX',tIX);

% set Matrices to hold time-series
M_0 = GetTimeIndexedData(hfig,'isAllCells');
setappdata(hfig,'M_0',M_0);
if ~exist('isSkipcIX','var'),
    cIX = getappdata(hfig,'cIX');
    setappdata(hfig,'M',M_0(cIX,:));
end
end

function [M,behavior,stim] = GetTimeIndexedData(hfig,isAllCells) %#ok<INUSD>
%{
% naming convention used:
M = GetTimeIndexedData(hfig);
M_0 = GetTimeIndexedData(hfig,'isAllCells');
%}

isZscore = getappdata(hfig,'isZscore');
% main data input
if ~isZscore,
    cellResp = getappdata(hfig,'CellResp');
    cellRespAvr = getappdata(hfig,'CellRespAvr');
else
    cellResp = getappdata(hfig,'CellRespZ');
    cellRespAvr = getappdata(hfig,'CellRespAvrZ');
end
Behavior_full = getappdata(hfig,'Behavior_full');
BehaviorAvr = getappdata(hfig,'BehaviorAvr');
stim_full = getappdata(hfig,'stim_full');
stimAvr = getappdata(hfig,'stimAvr');
% other params
isAvr = getappdata(hfig,'isAvr');
cIX = getappdata(hfig,'cIX');
tIX = getappdata(hfig,'tIX');

%% set data
if isAvr,
    if exist('isAllCells','var'),
        M = cellRespAvr(:,tIX);
    else
        M = cellRespAvr(cIX,tIX);
    end
    behavior = BehaviorAvr(:,tIX);
    stim = stimAvr(:,tIX);
else
    if exist('isAllCells','var'),
        M = cellResp(:,tIX);
    else
        M = cellResp(cIX,tIX);
    end
    behavior = Behavior_full(:,tIX);
    stim = stim_full(:,tIX);
end

setappdata(hfig,'behavior',behavior);
setappdata(hfig,'stim',stim);
end

function closefigure_Callback(hfig,~)
global EXPORT_autorecover;
EXPORT_autorecover = getappdata(hfig);
end

function fig = getParentFigure(fig)
% if the object is a figure or figure descendent, return the figure. Otherwise return [].
while ~isempty(fig) && ~strcmp('figure', get(fig,'type'))
    fig = get(fig,'parent');
end
end

function runscript(flag_script,var_script)
switch flag_script
    case 'push_cIX_gIX'
        UpdateIndices(var_script{:});
        RefreshFigure(var_script{1});
end
end