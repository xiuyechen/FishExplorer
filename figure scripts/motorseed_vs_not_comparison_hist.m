% motorseed comparison


%%
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
binEdges = -1:0.05:1;
nbins = length(binEdges)-1;
N_hist_notseed = zeros(18,nbins,2);
N_hist_seed = zeros(18,nbins,2);
for i_fish = 1:18
    ClusterIDs = [2,1];
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    %%
    isMotorseed = 0;
    setappdata(hfig,'isMotorseed',isMotorseed);
    [~,~,behavior] = UpdateTimeIndex(hfig);
    
    [~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
    if isMotorseed
        Reg = regressor_m;
    else
        Reg = regressor_m([1,3],:);
    end
    % regression
    Corr = corr(Reg',M_0');
    for i_reg = 1:2
        R = Corr(i_reg,:);
        [N,~] = histcounts(R,binEdges);
        N_hist_notseed(i_fish,:,i_reg) = N;
    end
    %%
    isMotorseed = 1;
    setappdata(hfig,'isMotorseed',isMotorseed);
    [~,~,behavior] = UpdateTimeIndex(hfig);
    
    [~,~,regressor_m,names_m] = GetMotorRegressor(behavior,i_fish);
    if isMotorseed
        Reg = regressor_m;
    else
        Reg = regressor_m([1,3],:);
    end
    % regression
    Corr = corr(Reg',M_0');
    for i_reg = 1:2
        R = Corr(i_reg,:);
        [N,~] = histcounts(R,binEdges);
        N_hist_seed(i_fish,:,i_reg) = N;
    end
end

%%
figure;
cmap = [1 0.4 0.4;0.4 0.4 0.4];
range_fish = 1:18;
m = length(range_fish);
N_type = cell(1,2);
N_type{1} = N_hist_seed;
N_type{2} = N_hist_notseed;
% i_plot = 0;

% left half
i_reg = 1;
for i_fish = range_fish    
    for i_type = 1:2
        i_plot = 1+2*(i_fish-1);%i_plot+1;
        subplot(m,2,i_plot);
        hold on;

        counts = squeeze(N_type{i_type}(i_fish,:,i_reg));
        histogram('BinEdges',binEdges,'BinCounts',counts,'FaceColor',cmap(i_type,:))
        
        %     ymax = max(max(N),max(N_shf));
        %     plot([0,0],[0,ymax],'r--');
        ymax = counts(24)+eps;
        xlim([-1,1]);
        ylim([0,ymax]);
        set(gca,'YTick',[])
        ylabel('a.u.')
        title(['Fish ' num2str(i_fish)])
    end
    end

% right half
i_reg = 2;
for i_fish = range_fish    
    for i_type = 1:2
        i_plot = 2*i_fish;%i_plot+1;
        subplot(m,2,i_plot);
        hold on;

        counts = squeeze(N_type{i_type}(i_fish,:,i_reg));
        histogram('BinEdges',binEdges,'BinCounts',counts,'FaceColor',cmap(i_type,:))
        
        %     ymax = max(max(N),max(N_shf));
        %     plot([0,0],[0,ymax],'r--');
        ymax = counts(24)+eps;
        xlim([-1,1]);
        ylim([0,ymax]);
        set(gca,'YTick',[])
        ylabel('a.u.')
        title(['Fish ' num2str(i_fish)])
    end
end