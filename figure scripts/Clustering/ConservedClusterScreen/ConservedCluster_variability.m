% load('C:\Users\xiuye\Dropbox\!Proj FishExplorer\output\D12screen_allfish_062018.mat')

i_fish = 8;

k_consrv = 6;
m = TF_fishrange{i_fish};
m1 = sum(logical(m),2);
U = find(m1>=k_consrv);

% ID_fishrange{1,i_fish};
%
P_sizes = cell(1,length(U));
for i_cluscount = 1:length(U)
    i_clus = U(i_cluscount);
    m2 = ID_fishrange{1,i_fish}(i_clus,:);
    %%
    P = [];
    for ii = 1:18
        IX = m2{ii};
        P = [P;ClusSize_fishrange{ii}(IX)];
    end
    P_sizes{i_cluscount} = P;
% m3 = cell2mat(m2');
end

%%
figure;hold on;
for i_cluscount = 1:length(U)
    y = P_sizes{i_cluscount};
   xv = ones(size(y))*i_cluscount;
    plot(xv,y,'o')
end
set(gca, 'YScale', 'log')
ylim([0,10000])

%% rank based on std
M_std = nan(1,length(U));
for i_cluscount = 1:length(U)
    y = P_sizes{i_cluscount};
    if length(y)>=k_consrv
        M_std(i_cluscount) = std(y)/mean(y);
    end
end

[~,IX_sort] = sort(M_std,'ascend');

[cIX,gIX,M_xyz_norm] = GetDefaultClustersFromLoad(hfig,i_fish);

[cIX,gIX] = SelectClusterRange(cIX,gIX,U);
gIX_last = gIX;
% numU = length(unique(gIX));

for i = 1:length(U)
    ix = IX_sort(i);    
    gIX(gIX_last==U(ix)) = i;
end


