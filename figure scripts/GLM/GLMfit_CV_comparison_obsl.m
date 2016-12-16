function [GLM_reg,GLM_ctrd,GLM_ctrd_k,gIX] = GLMfit_CV_comparison(M_0,Ctrd,behavior,stim,fishset,i_fish) % obsl

Data = M_0; %getappdata(hfig,'M_0');
%         Data = getappdata(hfig,'M');
% fishset = getappdata(hfig,'fishset');
% stim = getappdata(hfig,'stim');
% behavior = getappdata(hfig,'behavior');
% i_fish = getappdata(hfig,'i_fish');

% get stim/motor regressors
[~,names_s,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
[~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
Reg = vertcat(regressor_s,regressor_m)';
regnames = [names_s,names_m]';

%% Main 
myparpool = gcp;
% GLM with regressors as feature-matrix
tic
GLM_reg = zeros(size(Data,1),4);
parfor i = 1:size(Data,1),
    X = [ones(size(Reg,1),1) Reg];
    y = Data(i,:)';
    [~,~,~,~,stat] = regress(y,X);
   GLM_reg(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
end
t1=toc

% GLM with Auto-cluster centroids as regressors
tic
GLM_ctrd = zeros(size(Data,1),4);
parfor i = 1:size(Data,1),
    X = [ones(size(Ctrd,1),1) Ctrd];
    y = Data(i,:)';
    [~,~,~,~,stat] = regress(y,X);
   GLM_ctrd(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
end
t2=toc

% k-means
numK = size(Ctrd,2);
% M_0 = getappdata(hfig,'M_0');
[gIX,Ctrd_k] = kmeans(Data,numK,'distance','correlation');
tic
GLM_ctrd_k = zeros(size(Data,1),4);
parfor i = 1:size(Data,1),
    X = [ones(size(Ctrd_k,1),1) Ctrd_k];
    y = Data(i,:)';
    [~,~,~,~,stat] = regress(y,X);
   GLM_ctrd_k(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
end
t3=toc

filename = ['GLM_fish' num2str(i_fish) '.mat'];
save(filename,'GLM_reg','GLM_ctrd','GLM_ctrd_k','i_fish');

end

% function output = GLM_allcells(Data,Reg)
% myparpool = gcp;
% % GLM with regressors as feature-matrix
% tic
% GLM_reg = zeros(size(Data,1),4);
% parfor i = 1:size(Data,1),
%     X = [ones(size(Reg,1),1) Reg];
%     y = Data(i,:)';
%     [~,~,~,~,stat] = regress(y,X);
%    GLM_reg(i,:) = stat; % the R2 statistic % The R2 statistic can be negative for models without a constant, indicating that the model is not appropriate for the data.
% end
% t1=toc
% end
