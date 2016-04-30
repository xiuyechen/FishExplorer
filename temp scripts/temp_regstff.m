%% get cluster centroids (means) from GUI current selection
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
C = FindCentroid(hfig);
nClus = size(C,1);

stim = stim(2000:4000);
behavior = behavior(2000:4000);
%% get stimulus regressor
switch option
    case 1; % use all pre-defined stim-regressors
        regressors = GetStimRegressor(stim,fishset);
        M_regressor = zeros(length(regressors),length(regressors(1).im));
        for i = 1:length(regressors),
            M_regressor(i,:) = regressors(i).im;
        end
        regressor_s = M_regressor;
        
    case 2; % get stim regressors and find best match
        stimregset = 1:8;
        regressors = GetStimRegressor(stim,fishset);
        M_regressor = zeros(length(stimregset),length(regressors(1).im));
        for i = 1:length(stimregset),
            M_regressor(i,:) = regressors(stimregset(i)).im;
        end
        
        R = corr(M_regressor',C'); % row of R: each regressor
        [~,IX] = max(R,[],1);
        regressor_s_allclus = M_regressor(IX,:);
        
    case {3,4}; % alternative: using 'stim-lock' means
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
            timelists = getappdata(hfig,'timelists');
            
            C_mean = [];
            for i = 1:length(stimrange),
                i_stim = stimrange(i);
                if sum(stimset(i_stim).nReps)>3,
                    offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
                    tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
                    period = periods(i_stim);
                    C_3D_0 = reshape(C(:,tIX_),size(C,1),period,[]);
                    C_period = mean(C_3D_0,3);
                    nPeriods = length(tIX_)/period;
                    C_mean = horzcat(C_mean,repmat(C_period,1,nPeriods));
                    %             C_3D = zscore(C_3D_0,0,2);
                    %             H_raw = horzcat(H_raw,nanmean(nanstd(C_3D,0,3),2));
                else
                    offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
                    tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
                    C_mean = horzcat(C_mean,zeros(size(C,1),length(tIX_)));
                end
            end
        end
        
        regressor_s_allclus = C_mean;
end

%% get motor regressor
regressors = GetMotorRegressor(behavior);
motorregset = 1:3;
M_regressor = zeros(length(motorregset),length(regressors(1).im));
for i = 1:length(motorregset),
    M_regressor(i,:) = regressors(motorregset(i)).im;
end

% switch option
%     case 1;
regressor_m = M_regressor;
%     case {2,3};
%         R = corr(M_regressor',C'); % row of R: each regressor
%         [~,IX_m] = max(R,[],1);
%         regressor_m_allclus = M_regressor(IX_m,:);
% end

%% multi-regression
switch option
    case 1; % regression with all regs, stim-regs before motor regs
        regs = vertcat(regressor_s,regressor_m);
        orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?
        betas = C * orthonormal_basis; % to reconstitute: betas(i,:)*orthonormal_basis'
        % get ranking score: combined of motor coeffs
        H = sqrt(sum((betas(:,end-2:end)).^2,2));
    case {2,3}; % regression with stim-avr reg before motor regs
        betas = zeros(nClus,1+length(motorregset)); %size(regressor_s,1)+1);
        for i_clus = 1:nClus,
            regs = vertcat(regressor_s_allclus(i_clus,:),regressor_m);
            %             regs = vertcat(regressor_s_allclus(i_clus,:),regressor_m_allclus(i_clus,:));
            
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
        % get ranking score: combined motor coeffs
        H = sqrt(sum((betas(:,end-2:end)).^2,2));
    case 4; % regression with motor regs before stim-avr reg
        betas = zeros(nClus,1+length(motorregset)); %size(regressor_s,1)+1);
        for i_clus = 1:nClus,
            regs = vertcat(regressor_m,regressor_s_allclus(i_clus,:));
            orthonormal_basis = Gram_Schmidt_Process(regs');
            betas(i_clus,:) = C(i_clus,:) * orthonormal_basis;
        end
        % get ranking score: coeffs of (the single) stim-avr reg
        H = betas(:,end);
end

%% rank and plot
[gIX,rankscore] = SortH(H,gIX,numU,'descend');

figure;
h1 = subplot(2,1,1);
imagesc(orthonormal_basis);
title('All orthonormal bases');
xlabel('orthonormal basis #');
% make x-lables
nBases = size(orthonormal_basis,2);
s = [];
for i_basis = 1:nBases-3,
    s = [s,{['stim.' num2str(i_basis)]}];
end
s = [s,{'motor.R','motor.L','motor.F'}];
h1.XTick = 1:nBases;
h1.TickLength = [0,0];
h1.XTickLabel = s;
h1.XTickLabelRotation = 45;
ylabel('time ~ frames');

h2 = subplot(2,1,2);
imagesc(betas)
title('multi-regression coefficients');
xlabel('orthonormal basis #');
h2.XTick = 1:nBases;
h2.TickLength = [0,0];
h2.XTickLabel = s;
h2.XTickLabelRotation = 45;
ylabel('cluster ID');