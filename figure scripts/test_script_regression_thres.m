% test_script_normalized_corrcoeff

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);
%%
t1 = tic;
fishrange = 1;%GetFishRange;
Ys = zeros(length(fishrange),5);
Ym = zeros(length(fishrange),5);
Y_size = zeros(length(fishrange),2);
Y_shuffstd = zeros(length(fishrange),2);

for i_fishcount = 1:length(fishrange)
    i_fish = fishrange(i_fishcount);
    ClusterIDs = [2,1];
%     stimrange = 1;
%     [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs,stimrange);
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);
    
    %% corr with PT regressor
    
    % setappdata(hfig,'isStimAvr',1);
    % [M_0,M,behavior,stim] = UpdateTimeIndex(hfig);
    
    fishset = getappdata(hfig,'fishset');
    [~,names,regressors] = GetStimRegressor(stim,fishset,i_fish);
    reg_range = 2; % PT
    Reg = regressors(reg_range,:);
    R = corr(Reg',M_0');
    
    Ys(i_fishcount,1) = prctile(R,96);
    Ys(i_fishcount,2) = prctile(R,98);
    Ys(i_fishcount,3) = prctile(R,99);
    Ys(i_fishcount,4) = prctile(R,99.5);
    Ys(i_fishcount,5) = prctile(R,99.9);
    
    % bootstrap
    nFrames = size(M_0,2);
    M_shuff = M_0(:,randperm(nFrames));
    R = corr(Reg',M_shuff');
    Y_shuffstd(i_fishcount,1) = std(R);
    
    %% motor
    [~,~,regressor_m] = GetMotorRegressor(behavior,i_fish);
    Reg = regressor_m(1,:);
    R = corr(Reg',M_0');
    Ym(i_fishcount,1) = prctile(R,96);
    Ym(i_fishcount,2) = prctile(R,98);
    Ym(i_fishcount,3) = prctile(R,99);
    Ym(i_fishcount,4) = prctile(R,99.5);
    Ym(i_fishcount,5) = prctile(R,99.9);

    Y_size(i_fishcount,:) = size(M_0);
    
    % bootstrap
    nFrames = size(M_0,2);
    M_shuff = M_0(:,randperm(nFrames));
    R = corr(Reg',M_shuff');
    Y_shuffstd(i_fishcount,2) = std(R);
end
toc(t1)

%%
figure;
plot(fishrange,Ys)
xlabel('fish ID')
ylabel('corr coeff, cutoff value')
legend('96%','98%','99%','99.5%','99.9%','location','eastoutside')
title('stim')

figure;
plot(fishrange,Ym)
xlabel('fish ID')
ylabel('corr coeff, cutoff value')
legend('96%','98%','99%','99.5%','99.9%','location','eastoutside')
title('motor')

figure;
subplot(2,1,1);plot(fishrange,Y_size(:,1))
subplot(2,1,2);plot(fishrange,Y_size(:,2))

figure;
plot(fishrange,Y_shuffstd)

savedir = GetOutputDataDir;
M_prct = [96,98,99,99.5,99.9];
save('reg_cutoff_defS.mat','Ys','Ym','Y_size','Y_shuffstd','fishrange','M_prct');
%%
thres_reg = 0.5;
cIX = (find(R>thres_reg))';
wIX = R(cIX); % weight, i.e. corr.coeff
[~,I] = sort(wIX,'descend');
cIX = cIX(I);
gIX = ones(size(cIX));
wIX = wIX(I)';

UpdateIndices_Manual(hfig,cIX,gIX);

%% bootstrap
nFrames = size(M_0,2);
M_shuff = M_0(:,randperm(nFrames));
R = corr(Reg',M_shuff'); 
max(R) % ~0.07
prctile(R,99.9) % 0.1738; 0.0348
%%
nFrames = size(M_0,2);
step = floor(nFrames/10);
Ys = zeros(1,10);
X = step:step:nFrames;
for i = 1:length(X)
    xv = 1:X(i);
    M_shuff = M_0(:,randperm(length(xv)));
    R = corr(Reg(xv)',M_shuff');
    Ys(i) = prctile(R,99.9);
end
figure;
plot(X,Ys)
xlabel('# frames');
ylabel('99.9% cutoff for regression w PT')
% ans =    0.1186 % for 95%
    
%% histogram of corr.coeff, for all cells
figure('Position',[500,200,500,150]);
hold on;
bins = -1:0.05:1;
[N,~] = histcounts(R,bins);
histogram(R,bins,'FaceColor',[0.4 0.4 0.4]);%,'EdgeColor','none'
plot([thres_reg,thres_reg],[0,max(N)],'r--');
xlim([-1,1]);ylim([0,max(N)]);