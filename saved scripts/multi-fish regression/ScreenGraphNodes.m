% summarize conserved clusters as graph

% P1 = vertcat(Pairs_AllClusAllFish{:,1});
% P2 = vertcat(Pairs_AllClusAllFish{:,2});
% Nodes = unique([P1;P2],'rows');
% numNodes = length(Nodes);
% 
% [~,IX1] = ismember(P1,Nodes,'rows');
% [~,IX2] = ismember(P2,Nodes,'rows');
% 
% G = graph(IX1,IX2); figure;plot(G)
% bins = conncomp(G);
% %%
% thres_fishcount = round(length(range_fish)/2);
% range_bins = unique(bins);
% fishCountForG = zeros(length(range_bins),1);
% for i = range_bins,
%     IX = find(bins==range_bins(i));
%     fishCountForG(i) = length(unique(Nodes(IX,1)));
% end
%     
% 
% %% Plot anat
% anat_yx_norm = getappdata(hfig,'anat_yx_norm');
% 
% figure;
% hold on;
% image(anat_yx_norm)
% axis equal
% axis ij
% axis off
% 
% clrmap = hsv(round(max(bins)*1.1));
% 
% for i_node = 1:length(Nodes),
%     i_fish = Nodes(i_node,1);
%     clusIX = Nodes(i_node,2);
%     
%     GIX = VAR(i_fish).ClusGroup{3}.gIX;
%     IX = find(GIX == clusIX);
%     CIX_abs = VAR(i_fish).ClusGroup{3}.cIX_abs;
%     cIX_abs = CIX_abs(IX);
%     LoadFishDataWithoutTS(hfig,i_fish);
%     CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
%     xyz_norm = CellXYZ_norm(cIX_abs,:);
%     %    cmap = repmat(clrmap(bins(i_node),:)',1,length(cIX_abs));
%     scatter(xyz_norm(:,2),xyz_norm(:,1),2,clrmap(bins(i_node),:));%,20,clrmap(IX,:))
% end
% 
% %%
% anat_yx_norm = getappdata(hfig,'anat_yx_norm');
% 
% figure;
% hold on;
% image(anat_yx_norm)
% axis equal
% axis ij
% axis off
% 
% clrmap = hsv(round(max(bins)*1.1));
% range_fish_plot = unique(Nodes(:,1))';
% for i_fish = range_fish_plot,
%     [cIX,gIX,M_xyz_norm] = GetDefaultClustersFromLoad(hfig,i_fish);
%     IX_node = find(Nodes(:,1)==i_fish);
%     U = Nodes(IX_node,2);
%     clr = bins(IX_node);
%     for i = 1:length(U),
%         IX = find(gIX==U(i));
%         scatter(M_xyz_norm(IX,2),M_xyz_norm(IX,1),2,clrmap(clr(i),:));%,20,clrmap(IX,:))
%     end    
% end


%% JUST SCREEN NORMALLY


% k_consrv_ratio = 1/3;% or just set int number of fish......
k_consrv = 5; % at least in 5 more fish within this set


range_fish = 1:18;
for i_fish = range_fish,
    
    m = TF_fishrange{i_fish};
    m1 = sum(logical(m),2);
%     m1 = sum(m,2);
    U = find(m1>=k_consrv);
    [cIX_in,gIX_in,~,~,~,absIX] = GetDefaultClustersFromLoad(hfig,i_fish);
    
    cIX = [];gIX = [];
    for i = 1:length(U),
        IX = find(gIX_in==U(i));
        cIX = [cIX;cIX_in(IX)];
        gIX = [gIX;gIX_in(IX)];
    end
    
    % save cluster
    name = 'D12screen_6/18fish';
    clusgroupID = 10;
    clusIDoverride = 1;
    SaveCluster_Direct(cIX,gIX,absIX,i_fish,name,clusgroupID,clusIDoverride);
end
%%
% SaveVARwithBackup();
% %%
% 
% % fish6_clus47:right abducence
% i_fish = 6;
% ix = 47;
% 
% m = TF_fishrange{i_fish};
% fIX = find(m(ix,:));
% 
