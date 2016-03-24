
% addpath('C:\Users\engertlab\Dropbox\Github\FishExplorer\Gram-Schmidt Process')
% http://www.mathworks.com/matlabcentral/fileexchange/18843-gram-schmidt-process

%% get cluster centroids (means) from GUI current selection
C = f.FindCentroid(hfig);
nClus = size(C,1);
fishset = getappdata(hfig,'fishset');

%%
choice = 1;
switch choice
    case 1;
        [regressors, names] = GetStimRegressor(stim,fishset);
        M_regressor = zeros(length(regressors),length(regressors(1).im));
        for i = 1:length(regressors),
            M_regressor(i,:) = regressors(i).im;
        end
        regressor_s = M_regressor;
        
    case 2;
        %% get stim regressors and find best match
        stimregset = [1:8];
        [regressors, names] = GetStimRegressor(stim,fishset);
        M_regressor = zeros(length(stimregset),length(regressors(1).im));
        for i = 1:length(stimregset),
            M_regressor(i,:) = regressors(stimregset(i)).im;
        end
        
        R = corr(M_regressor',C'); % row of R: each regressor
        [~,IX] = max(R,[],1);
        regressor_s_allclus = M_regressor(IX,:);
        
    case 3;
        %% alternative: using 'stim-lock' means
        periods = getappdata(hfig,'periods');
        
        if fishset == 1,
            period = periods;
            C_3D_0 = reshape(C,size(C,1),period,[]);
            C_period = mean(C_3D_0,3);
            nPeriods = round(size(C,2)/period);
            C_mean = repmat(C_period,1,nPeriods);
        else
            stimset = getappdata(hfig,'stimset');
            stimrange = getappdata(hfig,'stimrange');
            periods = getappdata(hfig,'periods');
            tlists = getappdata(hfig,'tlists');
            
            C_mean = [];
            for i = 1:length(stimrange),
                i_stim = stimrange(i);
                if sum(stimset(i_stim).nReps)>3,
                    offset = length(horzcat(tlists{stimrange(1:i-1)}));% works for i=0 too
                    tIX_ = 1+offset:length(tlists{stimrange(i)})+offset;
                    period = periods(i_stim);
                    C_3D_0 = reshape(C(:,tIX_),size(C,1),period,[]);
                    C_period = mean(C_3D_0,3);
                    nPeriods = length(tIX_)/period;
                    C_mean = horzcat(C_mean,repmat(C_period,1,nPeriods)); %#ok<AGROW>
                    %             C_3D = zscore(C_3D_0,0,2);
                    %             H_raw = horzcat(H_raw,nanmean(nanstd(C_3D,0,3),2));
                else
                    offset = length(horzcat(tlists{stimrange(1:i-1)}));% works for i=0 too
                    tIX_ = 1+offset:length(tlists{stimrange(i)})+offset;
                    C_mean = horzcat(C_mean,zeros(size(C,1),length(tIX_))); %#ok<AGROW>
                end
            end
        end
        
        regressor_s_allclus = C_mean;
end

%% get motor regressor
regressors = GetMotorRegressor(fictive);
motorregset = 1:3;

switch choice
    case 1;
        M_regressor = zeros(length(motorregset),length(regressors(1).im));
        for i = 1:length(motorregset),
            M_regressor(i,:) = regressors(motorregset(i)).im;
        end
        regressor_m = M_regressor;
    case {2,3};                
        M_regressor = zeros(length(motorregset),length(regressors(1).im));
        for i = 1:length(motorregset),
            M_regressor(i,:) = regressors(motorregset(i)).im;
        end
        
        R = corr(M_regressor',C'); % row of R: each regressor
        [~,IX_m] = max(R,[],1);
        regressor_m_allclus = M_regressor(IX_m,:);
end

%%
switch choice
    case 1;
        regs = vertcat(regressor_s,regressor_m);
        orthonormal_basis = Gram_Schmidt_Process(regs');
        betas = C * orthonormal_basis;
    case {2,3};
        betas = zeros(nClus,2); %size(regressor_s,1)+1);
        for i_clus = 1:nClus,
            regs = vertcat(regressor_s_allclus(i_clus,:),regressor_m_allclus(i_clus,:));
            
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
        motorScore = horzcat(betas(:,end),IX_m');% for each cluster
end

figure;
subplot(2,1,1)
imagesc(orthonormal_basis)
subplot(2,1,2)
imagesc(betas)

%% or:
% betas = C * orthonormal_basis;% each row correspond to [beta1, beta2] of a cluster


