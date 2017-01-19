
h0 = figure;
InitializeAppData(h0);
% range_fish = [1:3,5:18];
for i_fish = range_fish
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,h0);
    
    % tic
    % [betas,stimcorr,motorcorr] = MultiMotorRegression(i_fish,M,stim,behavior);
    % toc
    
    M_fish_set = GetFishStimset();
    fishset = M_fish_set(i_fish); % one fewer param to load
    [~,~,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
    [~,~,regressor_m] = GetMotorRegressor(behavior,i_fish);
    regs = vertcat(regressor_s,regressor_m);
    orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?
    
    stimregs = orthonormal_basis(:,1:end-3);
    orthmotorregs = orthonormal_basis(:,end-2:end);
    stimcorr = max(corr(stimregs,M'),[],1)';
    orthmotorcorr = max(corr(orthmotorregs,M'),[],1)';
    rawmotorcorr = max(corr(regressor_m',M'),[],1)';
    
    [score,IX] = max([stimcorr,orthmotorcorr,rawmotorcorr],[],2);
    
    ix = find(IX==2);
    cIX_plot = cIX(ix);
    gIX_plot = IX(ix);
    %
    figure('Position',[50,100,800,1000]);
    I = LoadCurrentFishForAnatPlot(h0,cIX_plot,gIX_plot);
    DrawCellsOnAnat(I);
    
    dataDir = GetCurrentDataDir;
    saveDir = fullfile(dataDir,'motorsourceplot');
    if ~exist(saveDir, 'dir'), mkdir(saveDir), end;
    filename = fullfile(saveDir, num2str(i_fish));
    saveas(gcf, filename, 'png');
    close(gcf)
    % figure;DrawTimeSeries(h0,cIX_plot,gIX_plot);
end

%% this doesn't work well because hard thresholds work quite differently on different fish

% setappdata(h0,'isRefAnat',1);
% range_fish = [1:3,5:18];
for i_fish = range_fish;
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,h0);
    
    % tic
    % [betas,stimcorr,motorcorr] = MultiMotorRegression(i_fish,M,stim,behavior);
    % toc
    
    M_fish_set = GetFishStimset();
    fishset = M_fish_set(i_fish); % one fewer param to load
    [~,~,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
    [~,~,regressor_m] = GetMotorRegressor(behavior,i_fish);
%     regs = vertcat(regressor_s,regressor_m);
%     orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?
%     
%     stimregs = orthonormal_basis(:,1:end-3);
%     orthmotorregs = orthonormal_basis(:,end-2:end);
%     stimcorr = max(corr(stimregs,M'),[],1)';
%     orthmotorcorr = max(corr(orthmotorregs,M'),[],1)';
    [rawmotorcorr,IX] = max(corr(regressor_m',M'),[],1);
    
%     [score,IX] = max([stimcorr,orthmotorcorr,rawmotorcorr],[],2);
%     
    ix = find(rawmotorcorr>0.4);
    cIX_plot = cIX(ix);
    gIX_plot = IX(ix);
    %
    figure('Position',[50,100,800,1000]);
    I = LoadCurrentFishForAnatPlot(h0,cIX_plot,gIX_plot);
    DrawCellsOnAnat(I);
    
    dataDir = GetCurrentDataDir;
    saveDir = fullfile(dataDir,'simple_regression_plot_motor');
    if ~exist(saveDir, 'dir'), mkdir(saveDir), end;
    fn = fullfile(saveDir, num2str(i_fish));
    saveas(gcf, fn, 'png');
    close(gcf)
    % figure;DrawTimeSeries(h0,cIX_plot,gIX_plot);
end

%% need to develop a script that finds a good cutoff, to generate nice simple regression plots (e.g. raw motor correlation)


