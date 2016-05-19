

data_masterdir = GetCurrentDataDir();
range_fish = [2,3,4,5,6];
%% a) Histogram of Pairwise Distances for One Fish

%
clear T N
figure;
hold on;
for i = 1:length(range_fish),
    i_fish = range_fish(i);
    [f_cells,M,N(i),T(i),hfig] = getFullFishData(hfig,i_fish,VAR);
    
    num_select = 1000;%num_cells;
    idx_select = randperm(N(i), num_select);
    
    f_cells_select = f_cells(idx_select,:);
 
    
    %% Correlation 
    tic
    Dpairs = 1-pdist(f_cells_select,'correlation');
    % Dpairs is in correlation, not corr dist
    toc

    %% Plot Histogram
    binwidth = .005;
    edges = -1:binwidth:1;
    

    [counts,edges] = histcounts(Dpairs,edges);
    allBins = edges(1:end-1) + binwidth/2;
    
    allCorrVals = counts;
    leg{i} = ['Fish ',num2str(i_fish)];
    
    plot(allBins,allCorrVals/trapz(allBins,allCorrVals));
    xlabel('Pairwise Correlations');
    ylabel('Normalized Frequency');
    xlim([-0.2,0.2]);
    

end

% Compare to Random Data
T_control = T(end);
f_cells_control = normrnd(0,1,num_select,T_control);
Dpairs_control = 1-pdist(f_cells_control,'correlation');

[ctrlCounts,edgesCtrl] = histcounts(Dpairs_control,edges);
control_allBins = edgesCtrl(1:end-1) + binwidth/2;
control_allCorrVals = ctrlCounts;

plot(control_allBins,control_allCorrVals/trapz(control_allBins,control_allCorrVals)...
    ,'k');

leg{i+1} = 'Control';
legend(leg)

%% b) Percentile Threshold Plot

% Percentile Thresh
threshPerc = .99;
figure;
xlabel('Number of Frames');
ylabel('99th Percentile Correlation');

hold on;
for num_fish = 1:length(range_fish),
    
    i_fish = range_fish(num_fish);
    [f_cells,M,N(num_fish),T(num_fish),hfig] = getFullFishData(hfig,i_fish,VAR);
    
    num_select = 1000;%num_cells;
    idx_select = randperm(N(num_fish), num_select);
    
    f_cells_select = f_cells(idx_select,:);

   
    T_trunc = T(num_fish)/10:T(num_fish)/10:T(num_fish); %hardcode for now
    clear threshCorr;
    for i = 1:length(T_trunc)
        f_trunc = f_cells_select(:,1:T_trunc(i));             
        
        % Correlation Distance
        tic
        Dpairs = pdist(f_trunc,'correlation');
        % 138 sec for 86k cells
        % 1.8 sec for 10k cells
        toc
        %figure;
        %h = histogram(1-Dpairs);%Plot corr, not corr dist

        sortCorr = sort(1-Dpairs);
        threshIDX = round(length(sortCorr)*threshPerc);
        threshCorr(i) = sortCorr(threshIDX);
    end
    plot(T_trunc,threshCorr,'o-')
    leg2{num_fish} = ['Fish ',num2str(i_fish)];
end

% Compare to Random Data
T_control = T(end);
f_cells_control = normrnd(0,1,num_select,T_control);
Dpairs_control = 1-pdist(f_cells_control,'correlation');


T_trunc = T_control/10:T_control/10:T_control; %hardcode for now
clear threshCorr;
clear threshCorr;
for i = 1:length(T_trunc)
    f_trunc = f_cells_control(:,1:T_trunc(i));
    
    % Correlation Distance
    tic
    Dpairs = pdist(f_trunc,'correlation');
    % 138 sec for 86k cells
    % 1.8 sec for 10k cells
    toc
    %figure;
    %h = histogram(1-Dpairs);%Plot corr, not corr dist
    sortCorr = sort(1-Dpairs);
    threshIDX = round(length(sortCorr)*threshPerc);
    threshCorr(i) = sortCorr(threshIDX);
end
plot(T_trunc,threshCorr,'ko-')
leg2{num_fish+1} = 'control';
ylim([0 .5]);
legend(leg2);
grid on;