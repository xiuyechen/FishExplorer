% Pairwise distance for each auto-cluster
data_dir = GetCurrentDataDir();
load(fullfile(data_dir,'AllCentroids.mat'));

range_fish = 1:13;

xv = 1:20:800;
numFish = length(range_fish);
N = zeros(numFish,length(xv)-1);
N_sim = N;
Avr = zeros(numFish,1);
Avr_sim = Avr;
for i_fish = range_fish,
    LoadFishDataWithoutTS(hfig,i_fish);
    
    gIX = VAR(i_fish).ClusGroup{3}(1).gIX;
    
    numClus = length(AllCentroids{i_fish}.XYZn);
    sd_tot = zeros(numClus,1);
    for i_clus = 1:numClus,
        clus_XYZ = AllCentroids{i_fish}.XYZn{i_clus};
        
        sd = std(clus_XYZ);
        sd_tot(i_clus) = sqrt(sd(1)^2 + sd(2)^2 + sd(3)^2);       
    end
    
    % Simulate clusters
    cIX_abs = VAR(i_fish).ClusGroup{3}(1).cIX_abs;    
    sim_sd_tot = zeros(numClus,1);
    for i_clus = 1:numClus,
        clus_XYZ = AllCentroids{i_fish}.XYZn{i_clus};
        clussize = size(clus_XYZ,1);
        bootXYZ = zeros(clussize,3);
        
        IX = randi(length(cIX_abs),clussize,1);
        sim_cIX_abs = cIX_abs(IX);
        CellXYZ = getappdata(hfig,'CellXYZ');
        sim_clus_XYZ = CellXYZ(sim_cIX_abs,:);

        sd = std(sim_clus_XYZ);
        sim_sd_tot(i_clus) = sqrt(sd(1)^2 + sd(2)^2 + sd(3)^2);        
    end
    
    %%
% figure;
%     hold on;
    h = histogram(sd_tot,xv);
    N(i_fish,:) = h.Values';
    h = histogram(sim_sd_tot,xv);
    N_sim(i_fish,:) = h.Values';
    Avr(i_fish) = mean(sd_tot);
    Avr_sim(i_fish) = mean(sim_sd_tot);    
end

%% One-time conversion
k = 0.406; % um per pixel

Avr = k*Avr;
Avr_sim = k*Avr_sim;
xv = k*xv;

%% Plot figure
figure('Position',[100,400,100,300]);hold on;
scatter(ones(13,1),Avr,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
scatter(2*ones(13,1),Avr_sim,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
inc = 0.3;
scatter(1+inc,mean(Avr),10,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',1,'MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',1);
scatter(2+inc,mean(Avr_sim),10,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',1,'MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',1);

% plot([1-inc,1+inc],[mean(Avr),mean(Avr)],'k','Linewidth',1.5)
% plot([2-inc,2+inc],[mean(Avr_sim),mean(Avr_sim)],'k','Linewidth',1.5)
errorbar(1+inc,mean(Avr),std(Avr)/sqrt(length(Avr)),'k','Linewidth',1.5)
errorbar(2+inc,mean(Avr_sim),std(Avr_sim)/sqrt(length(Avr)),'k','Linewidth',1.5)
xlim([0.5,2.5])
ylim([0,200])
set(gca,'XTickLabels',{'Data','Shuffled'},'XTickLabelRotation',45);
ylabel('average spread (micron)')

%%
% figure;hold on;
% plot(mean(N,1),'b')
% plot(mean(N_sim,1),'r')