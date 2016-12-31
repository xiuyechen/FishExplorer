function [betas,stimcorr,motorcorr] = MultiMotorRegression(i_fish,M,stim,behavior)

% stim = getappdata(hfig,'stim');
% behavior = getappdata(hfig,'behavior');
% i_fish = getappdata(hfig,'i_fish');
M_fish_set = GetFishStimset();
fishset = M_fish_set(i_fish); % one fewer param to load
% fishset = getappdata(hfig,'fishset'); 

% (this code works the same way for clusters or individual cells)
nClus = size(M,1);

[~,~,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
[~,~,regressor_m] = GetMotorRegressor(behavior,i_fish);

regs = vertcat(regressor_s,regressor_m);
orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?

betas = zeros(nClus,size(orthonormal_basis,2)+1);
X = [ones(size(orthonormal_basis,1),1),orthonormal_basis];
for i_clus = 1:nClus,
    y = M(i_clus,:)';
    betas(i_clus,:) = regress(y,X)';
end

stimcorr = max(betas(:,1:end-3),[],2);
motorcorr = max(betas(:,end-2:end),[],2);
