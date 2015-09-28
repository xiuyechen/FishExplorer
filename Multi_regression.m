
addpath('C:\Users\engertlab\Dropbox\Github\FishExplorer\Gram-Schmidt Process')

%% get cluster centroids (means) from GUI current selection
C = f.FindCentroid(hfig);
nClus = size(C,1);

%% get stim regressors and find best match
% stimregset = [8,9,16];

fishset = getappdata(hfig,'fishset');
[regressors, names] = GetStimRegressor(stim,fishset);
M_regressor_s = zeros(length(stimregset),size(regressors,2));
for i = 1:length(stimregset),
    M_regressor_s(i,:) = regressors(stimregset(i)).im;
end

R = corr(M_regressor_s',C'); % row of R: each regressor
[~,IX] = max(R,[],1);
regressor_s_allclus = M_regressor_s(IX,:);

%% get motor regressor
regressors = GetMotorRegressor(fictive);
regressor_m = regressors(1).im;

%%
betas = zeros(nClus,2);
for i_clus = 1:nClus,
    regs = vertcat(regressor_s_allclus(i_clus,:),regressor_m);
    
    % get orthonormal_basis    
    orthonormal_basis = Gram_Schmidt_Process(regs');
    
    %% check orthonormal_basis
    % figure;imagesc(orthonormal_basis)
    % norm(orthonormal_basis(:,1))
    % norm(orthonormal_basis(:,2))
    % dot(orthonormal_basis(:,1),orthonormal_basis(:,2))
    % % compare to raw regressors
    % figure;
    % subplot(311);hold on;
    % i = 1;
    % plot(regs(i,:)/norm(regs(i,:)),'r');plot(orthonormal_basis(:,1)','k:')
    % subplot(312);hold on;
    % i = 2;
    % plot(regs(i,:)/norm(regs(i,:)),'r');plot(orthonormal_basis(:,i)','k:')
    % subplot(313);hold on;
    % i = 3;
    % plot(regs(i,:)/norm(regs(i,:)),'r');plot(orthonormal_basis(:,i)','k:')
    
    %% Get multi-regression coefficients
    % FunctionalActivity = beta1*reg1 + beita2*reg2_orth
    betas(i_clus,:) = C(i_clus,:) * orthonormal_basis;
end

%% or:
% betas = C * orthonormal_basis;% each row correspond to [beta1, beta2] of a cluster
motorScore = betas(:,end);% for each cluster

