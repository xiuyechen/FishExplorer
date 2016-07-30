% Pairwise distance for each auto-cluster

range_fish = 1:18;

xv = 1:20:800;
numFish = length(range_fish);
N = zeros(numFish,length(xv)-1);
N_sim = N;
Avr = zeros(numFish,1);
Avr_sim = Avr;

 k_xy_um = 0.798; % 0.406 um per pixel for XY res FOR NON-REGISTERED FISH
 k_z_um = 2;
%   ZBrain: x/y/z = 0.798/0.798/2um

% k_um = 0.406; % um per pixel
% k_zres_toXYratio = 2.5;

%%
for i_fish = range_fish,
    [cIX,gIX,M_xyz_norm,~,cIX_abs] = GetDefaultClustersFromLoad(hfig,i_fish,'0.7');
    
    U = unique(gIX);
    numClus = length(U);
    sd_tot = zeros(numClus,1);
    for i_clus = 1:numClus,
        IX = find(gIX==U(i_clus));
        clus_XYZ = M_xyz_norm(IX,:);
        
        % scale X/Y/Z distance from coordinates (Z res is lower)
        clus_XYZ_scaled = zeros(size(clus_XYZ));
        clus_XYZ_scaled(:,1) = clus_XYZ(:,1)*k_xy_um;
        clus_XYZ_scaled(:,2) = clus_XYZ(:,2)*k_xy_um;
        clus_XYZ_scaled(:,3) = clus_XYZ(:,3)*k_z_um;
        
        sd = std(clus_XYZ_scaled);
        sd_tot(i_clus) = sqrt(sd(1)^2 + sd(2)^2 + sd(3)^2);
    end
    
    % Simulate clusters
    sim_sd_tot = zeros(numClus,1);
    for i_clus = 1:numClus,
        IX = find(gIX==U(i_clus));
        clussize = length(IX);
        bootXYZ = zeros(clussize,3);
        
        IX_rand = randi(length(cIX_abs),clussize,1);
        sim_cIX_abs = cIX_abs(IX_rand);
        CellXYZ_norm = getappdata(hfig,'CellXYZ_norm');
        sim_clus_XYZ = CellXYZ_norm(sim_cIX_abs,:);
        
        % scale X/Y/Z distance from coordinates (Z res is lower)
        sim_clus_XYZ_scaled = zeros(size(sim_clus_XYZ));
        sim_clus_XYZ_scaled(:,1) = sim_clus_XYZ(:,1)*k_xy_um;
        sim_clus_XYZ_scaled(:,2) = sim_clus_XYZ(:,2)*k_xy_um;
        sim_clus_XYZ_scaled(:,3) = sim_clus_XYZ(:,3)*k_z_um;
        
        sd = std(sim_clus_XYZ_scaled);
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

%% Plot figure
figure('Position',[100,400,150,160]);hold on;
scatter(ones(numFish,1),Avr,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
scatter(2*ones(numFish,1),Avr_sim,20,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',0.5,'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
inc1 = 0.2;
inc = 0.15;
x = 1;
Y = Avr;
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);

x = 2;
Y = Avr_sim;
plot([x-inc1,x+inc1],[mean(Y),mean(Y)],'k');
err = std(Y)/sqrt(length(Y));
plot([x-inc,x+inc],[mean(Y)+err,mean(Y)+err],'color',[1,0.5,0.5]);
plot([x-inc,x+inc],[mean(Y)-err,mean(Y)-err],'color',[1,0.5,0.5]);

% scatter(1+inc,mean(Avr),10,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',1,'MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',1);
% scatter(2+inc,mean(Avr_sim),10,'MarkerEdgeColor',[0,0,0],'MarkerEdgeAlpha',1,'MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',1);

% plot([1-inc,1+inc],[mean(Avr),mean(Avr)],'k','Linewidth',1.5)
% plot([2-inc,2+inc],[mean(Avr_sim),mean(Avr_sim)],'k','Linewidth',1.5)
% errorbar(1+inc,mean(Avr),std(Avr)/sqrt(length(Avr)),'k','Linewidth',1.5)
% errorbar(2+inc,mean(Avr_sim),std(Avr_sim)/sqrt(length(Avr)),'k','Linewidth',1.5)
xlim([0.5,2.5])
ylim([0,220])
set(gca,'XTickLabels',{'Data','Shuffled'},'XTickLabelRotation',45);
ylabel('avr. spread (um)')

%%
% figure;hold on;
% plot(mean(N,1),'b')
% plot(mean(N_sim,1),'r')