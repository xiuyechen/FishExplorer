
% GLM for all cells, using all normal regressors
% load a fish with the appropriate stim-range

Data = getappdata(hfig,'M_0');
%         Data = getappdata(hfig,'M');
thres_reg = getappdata(hfig,'thres_reg');
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
% cIX_in = getappdata(hfig,'cIX');
% gIX_in = getappdata(hfig,'gIX');
% numK = getappdata(hfig,'numK');
i_fish = getappdata(hfig,'i_fish');

% get stim/motor regressors
[~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
[~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
Reg = vertcat(regressor_s,regressor_m)';
regnames = [names_s,names_m]';

%%
tic
GLM_reg = zeros(size(Data,1),4);
for i = 1:size(Data,1),
    X = [ones(size(Reg,1),1) Reg];
    y = Data(i,:)';
    [B,~,~,~,stat] = regress(y,X);
   GLM_reg(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
end
t1=toc
% run-time: ~20min for all ~80k cells * 5k frames

% GLM with Auto-cluster centroids as regressors
Ctrd = Centroids_Auto07_multi{8}';
tic
GLM_ctrd = zeros(size(Data,1),4);
for i = 1:size(Data,1),
    X = [ones(size(Ctrd,1),1) Ctrd];
    y = Data(i,:)';
    [~,~,~,~,stat] = regress(y,X);
   GLM_ctrd(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
end
t2=toc

%% k-means
numK = 139;
% M_0 = getappdata(hfig,'M_0');
[~,Ctrd_k] = kmeans(Data,numK,'distance','correlation');
tic
GLM_ctrd_k = zeros(size(Data,1),4);
for i = 1:size(Data,1),
    X = [ones(size(Ctrd_k,1),1) Ctrd_k];
    y = Data(i,:)';
    [~,~,~,~,stat] = regress(y,X);
   GLM_ctrd_k(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
end
t3=toc

save('GLM_fish8.mat','GLM_ctrd_k','GLM_reg','GLM_ctrd');

%% prediction

GLM_CV_allfish_script; % using GLM_allcells_CV.m
% GLMfit_CV_comparison; % obsl


%%
figure;hist(GLM_reg)
xlabel('r^2');
ylabel('cell count');
