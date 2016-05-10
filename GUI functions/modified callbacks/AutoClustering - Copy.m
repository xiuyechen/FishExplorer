function AutoClustering(hfig,f,i_fish,savename,isWkmeans)
tic
setappdata(hfig,'i_fish',i_fish);

cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
M = getappdata(hfig,'M');
M_0 = getappdata(hfig,'M_0');
thres_size = 10;
thres_split = getappdata(hfig,'thres_split');
thres_stimlock = 1.0;
thres_merge = 0.9;%getappdata(hfig,'thres_merge');
% thres_silh = 0.4;
% isWkmeans = getappdata(hfig,'isWkmeans');

%% kmeans
if isWkmeans,
    numK1 = 20;
    disp(['kmeans k = ' num2str(numK1)]);
    tic
    rng('default');% default = 0, but can try different seeds if doesn't converge
    if numel(M)*numK1 < 10^7 && numK1~=1,
        disp('Replicates = 5');
        gIX = kmeans(M,numK1,'distance','correlation','Replicates',5);
    elseif numel(M)*numK1 < 10^8 && numK1~=1,
        disp('Replicates = 3');
        gIX = kmeans(M,numK1,'distance','correlation','Replicates',3);
    else
        gIX = kmeans(M,numK1,'distance','correlation');
    end
    toc
%     f.SaveCluster_Direct(hfig,cIX,gIX,'k=20');
else
    [gIX, numK1] = SqueezeGroupIX(gIX);
end

%% pushbutton_iter_split(hObject,~);
numK2 = 20;
disp(['2nd tier kmeans k = ' num2str(numK2)]);

gIX_old = gIX;
for i = 1:numK1,
%     disp(['i = ' num2str(i)]);
    IX = find(gIX_old == i);
    M_sub = M_0(IX,:);
    
    if numK2<length(IX),
        [gIX_sub,C] = kmeans(M_sub,numK2,'distance','correlation');
    else
        [gIX_sub,C] = kmeans(M_sub,length(IX),'distance','correlation');
    end
    gIX(IX) = (i-1)*numK2+gIX_sub;
end

% f.SaveCluster_Direct(hfig,cIX,gIX,'k20x20',1);

%% Regression with the centroid of each cluster
disp('auto-reg-clus');
f.UpdateIndices(hfig,cIX,gIX);
[cIX,gIX,~] = AllCentroidRegression_direct(hfig);
f.UpdateIndices(hfig,cIX,gIX);

SaveCluster_Direct(hfig,cIX,gIX,'k20x20_reg',1);

%% Merge
disp('merge');
[gIX, numU] = Merge_direct(f,thres_merge,M_0,cIX,gIX);

% [gIX, numU] = SqueezeGroupIX(gIX);
% C = FindCentroid(hfig);
% IX_tree = clusterdata(C,'criterion','distance',...
%     'distance','correlation','cutoff',0.1);
% 
% % update gIX
% temp = zeros(size(gIX));
% for i = 1:numU,
%     temp(gIX==i) = IX_tree(i);
% end
% gIX = temp;
% cIX = getappdata(hfig,'cIX');
% f.UpdateIndices(hfig,cIX,gIX,length(unique(IX_tree)));
% f.RefreshFigure(hfig);
%% Regression with the centroid of each cluster, round 2
disp('auto-reg-clus, round 2');
f.UpdateIndices(hfig,cIX,gIX);
[cIX,gIX,~] = AllCentroidRegression_direct(hfig);
f.UpdateIndices(hfig,cIX,gIX);

SaveCluster_Direct(hfig,cIX,gIX,'k20x20_reg_round2',1);

%% size threshold
thres_size = getappdata(hfig,'thres_size');
U = unique(gIX);
numU = length(U);
for i=1:numU,
    if length(find(gIX==U(i)))<thres_size,
        cIX(gIX==U(i)) = [];
        gIX(gIX==U(i)) = [];
    end
end
[gIX, numU] = SqueezeGroupIX(gIX);

%% update GUI
if isempty(gIX),
    errordlg('nothing to display!');
    return;
end
f.UpdateIndices(hfig,cIX,gIX,numK1);
SaveCluster_Direct(hfig,cIX,gIX,savename,1);%'Full_autoclus');
f.RefreshFigure(hfig);
UpdateClustersGUI(hfig);
beep;
toc
end
