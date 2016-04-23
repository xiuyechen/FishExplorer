% Pairwise distance for each auto-cluster

i_fish = 8;
% LoadFishDataWithoutTS(hfig,i_fish);
LoadFullFish(hfig,i_fish);

numClus = length(AllCentroids{i_fish}.XYZn);
sd_tot = zeros(numClus,3);
H = zeros(numClus,4);
H_corr = zeros(numClus,1);
for i_clus = 1:numClus,
    clus_XYZ = AllCentroids{i_fish}.XYZn{i_clus};
    D = pdist(clus_XYZ);
    D_srt = sort(D);
    H(i_clus,1) = mean(D_srt);
    H(i_clus,2) = mean(D_srt(1:round(end/2)));
    H(i_clus,3) = median(D);
    H(i_clus,4) = std(D);
    H_corr(i_clus) = ;
    
    sd = std(clus_XYZ);
    sd_tot(i_clus,1) = sqrt(sd(1)^2 + sd(2)^2 + sd(3)^2);            
    sd_tot(i_clus,2) = size(clus_XYZ,1);
    sd_tot(i_clus,3) = sd_tot(i_clus,1)/size(clus_XYZ,1);
    %%
    cIX_abs = VAR(i_fish).ClusGroup{3}(1).cIX_abs;
    gIX = VAR(i_fish).ClusGroup{3}(1).gIX;
end

%%
figure;
xv = 1:10:400;
subplot(141);histogram(H(:,1),xv);
subplot(142);histogram(H(:,2),xv);
subplot(143);histogram(H(:,3),xv);
subplot(144);histogram(H(:,4),xv);
%%
figure; histogram(sd_tot(:,3));
%%
figure;scatter(sd_tot(:,2),sd_tot(:,1));

%% Simulate clusters
anat_stack = getappdata(hfig,'anat_stack');
cIX_abs = VAR(i_fish).ClusGroup{3}(1).cIX_abs;
[s1,s2,s3] = size(anat_stack);

sim_sd_tot = zeros(numClus,3);
sim_H = zeros(numClus,4);
for i_clus = 1:numClus,
    clus_XYZ = AllCentroids{i_fish}.XYZn{i_clus};
    clussize = size(clus_XYZ,1);
    bootXYZ = zeros(clussize,3);
    
    
    IX = randi(length(cIX_abs),clussize,1);
    sim_cIX_abs = cIX_abs(IX);
    CellXYZ = getappdata(hfig,'CellXYZ');
    sim_clus_XYZ = CellXYZ(sim_cIX_abs,:);
    
    D = pdist(sim_clus_XYZ);
    D_srt = sort(D);
    sim_H(i_clus,1) = mean(D_srt);
    sim_H(i_clus,2) = mean(D_srt(1:round(end/2)));
    sim_H(i_clus,3) = median(D);
    sim_H(i_clus,4) = std(D);
    
    
    sd = std(sim_clus_XYZ);
    sim_sd_tot(i_clus,1) = sqrt(sd(1)^2 + sd(2)^2 + sd(3)^2);            
    sim_sd_tot(i_clus,2) = size(sim_clus_XYZ,1);
    sim_sd_tot(i_clus,3) = sim_sd_tot(i_clus,1)/size(sim_clus_XYZ,1);

end
figure;
xv = 1:10:800;
subplot(141);
nout = histogram(H(:,1),xv);
subplot(142);histogram(H(:,2),xv);
subplot(143);histogram(H(:,3),xv);
subplot(144);histogram(H(:,4),xv);
figure; histogram(sim_sd_tot(:,1));