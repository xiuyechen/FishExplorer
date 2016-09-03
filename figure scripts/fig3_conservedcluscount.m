
Ratio = zeros(1,18);
for i_fish = 1:18,
    L1 = length(VAR(i_fish).ClusGroup{10}(1).gIX);
    L2 = length(VAR(i_fish).ClusGroup{6}(1).gIX);
    Ratio(i_fish) = L1/L2;
end

Ratio_clus = zeros(1,18);
for i_fish = 1:18,
    L1 = length(unique(VAR(i_fish).ClusGroup{10}(1).gIX));
    L2 = length(unique(VAR(i_fish).ClusGroup{6}(1).gIX));
    Ratio_clus(i_fish) = L1/L2;
end

%%
figure('Position',[500,500,300,200]);
bar(Ratio,'EdgeColor','w','FaceColor',[0.5,0.5,0.5])
ylim([0,1])
xlim([0.5,18.5])
xlabel('fish #')
ylabel('% cells conserved')
%%
figure('Position',[500,500,300,200]);
bar(Ratio_clus,'EdgeColor','w','FaceColor',[0.5,0.5,0.5])
ylim([0,1])
xlim([0.5,18.5])
xlabel('fish #')
ylabel('% clusters conserved')