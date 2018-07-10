i_fish = 8;

k_consrv = 6;
m = TF_fishrange{i_fish};
m1 = sum(logical(m),2);
U = find(m1>=k_consrv);

m = ID_fishrange{1};
%
P_sizes = cell(1,length(U));
for i_cluscount = 1:length(U)
    i_clus = U(i_cluscount);
    m2 = m(i_clus,:);
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