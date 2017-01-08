%%
i_fish = 8;
[cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,h0);

%% 
% % run this again, use M_0 instead of CellResp (stimrange)
% M = getappdata(hfig,'M'); % load A0.7
% gIX = getappdata(hfig,'gIX');
% cIX = getappdata(hfig,'cIX');
% numClus = length(unique(gIX));
% 
% M_0 = getappdata(hfig,'M_0');
tic;
[coeff,score,latent,tsquared,explained,mu] = pca(M_0); % ~52 sec on linux
toc

im = coeff;
im = im.*repmat(explained',size(coeff,1),1);
figure;
imagesc(im(:,1:numClus))
topPCs = coeff(:,1:numClus);

%% GLM fit

clusmeans = FindClustermeans(gIX,M);
C_betas = cell(1,2);
C_rsq = cell(1,2);
C_bases = {clusmeans,topPCs'};

M_this = M_0;%M;

for i_model = 1:2
    bases = C_bases{i_model};
    range_cell = 1:size(M_this);%1:100:size(M_0,1);
    M_betas = zeros(140,length(range_cell));
    M_rsq = zeros(1,length(range_cell));
    tic
    parfor i_count = 1:length(range_cell),
        i_cell = range_cell(i_count);
        X = [ones(size(bases,2),1),bases'];
        y = M_this(i_cell,:)';
        
        [b,~,~,~,stat] = regress(y,X);
        M_betas(:,i_count) = b;
        M_rsq(:,i_count) = stat(1);
    end
    toc
    C_betas{i_model} = M_betas;
    C_rsq{i_model} = M_rsq;
end

%% plot betas, given cluster rank
% im = M;
[~,I_clus] = sort(gIX);
% im = im(I_clus,:);

figure;
m = [];

for i = 1:2
    M_betas = C_betas{i};
%     M_rsq = C_rsq{i};
%     m{i} = M_betas(:,find(M_rsq > thres_rsq));
%     [~,IX] = sort(max(m{i},[],1),'descend');
%     m_srt{i} = m{i}(:,IX);
%     [~,IX] = sort(max(m_srt{i},[],2),'descend');
%     m_srt{i} = m_srt{i}(IX,:);
    
    subplot(1,2,i)
    imagesc(M_betas(:,I_clus))
    colormap jet
end

%% plot max beta trace, given cluster rank

figure;

for i = 1:2
    M_betas = C_betas{i};
    A = M_betas(:,I_clus);
    [~,IX] = max(abs(A),[],1);
    subplot(1,2,i)
    plot(IX,'.')
end

%% rank cells based on max beta only! (not knowing cluster assignment)

figure;
m = [];

for i = 1:2
    M_betas = C_betas{i};
%     M_rsq = C_rsq{i};
    %     m{i} = M_betas(:,find(M_rsq > thres_rsq));
    [~,IX_raw] = max(abs(M_betas),[],1);
    [a,b]=sort(IX_raw);
    
    subplot(1,2,i)
    imagesc(M_betas(:,b))
    colormap jet
end


%% for M_0 run: rank cells based on max beta only! (not knowing cluster assignment)

figure;
m = [];

for i = 1,%:2
    M_betas = C_betas{i}(:,cIX);
%     M_rsq = C_rsq{i};
    %     m{i} = M_betas(:,find(M_rsq > thres_rsq));
    [maxvalues,IX_raw] = max(abs(M_betas),[],1);
    [a,b]=sort(IX_raw);
    
    subplot(1,2,i)
    imagesc(M_betas(:,b))
    colormap jet
end

%% find threshold to screen cells from M_0 that are well spanned

thres_maxbeta = 0.5;
M_rsq(cIX);

figure;hist(M_rsq(cIX),0:0.01:1)
counts = hist(maxvalues,0:0.1:6);

thres_rsq = 0.5;
%%

figure;
m = [];

for i = 1,%:2
    M_betas = C_betas{i};
    M_rsq = C_rsq{i};
    %     m{i} = M_betas(:,find(M_rsq > thres_rsq));
    [maxvalues,IX_raw] = max(abs(M_betas),[],1);
    IX_cell = find(M_rsq>thres_rsq);
    IX_thres = IX_raw(IX_cell);
    [a,b]=sort(IX_thres);
    
    subplot(1,2,i)
    imagesc(M_betas(:,IX_cell(b)))
    colormap hot
end

%% percentage of overlap of well-spanned cells
thres_rsq = 0.5;
i = 1;
M_rsq = C_rsq{i};
[maxvalues,IX_raw] = max(abs(M_betas),[],1);
IX_cell_clus = find(M_rsq>thres_rsq);
i = 2;
M_rsq = C_rsq{i};
[maxvalues,IX_raw] = max(abs(M_betas),[],1);
IX_cell_pca = find(M_rsq>thres_rsq);
r_clus_pca = length(intersect(IX_cell_clus,IX_cell_pca))/length(union(IX_cell_clus,IX_cell_pca))
r_clus_cIX = length(intersect(IX_cell_clus,cIX))/length(union(IX_cell_clus,cIX))
r_pca_cIX = length(intersect(IX_cell_pca,cIX))/length(union(IX_cell_pca,cIX))
% ans =
% 
%     0.7860
    
%%
% use first _n_ PC's to reconstruct functional data
M_1 = score*coeff';

% Autoclus
cIX_in = (1:size(M_1,1))';
gIX_in = ones(size(cIX_in));
cIX_reg = cIX_in; 
tic;
[cIX_pca,gIX_pca] = AutoClustering(cIX_in,gIX_in,M_1,cIX_reg);
toc;

% CV between this and A0.7
figure;
% I = LoadCurrentFishForAnatPlot(h0);
I = LoadCurrentFishForAnatPlot(h0,cIX_pca,gIX_pca);
DrawCellsOnAnat(I);
figure;DrawTimeSeries(h0,cIX_pca,gIX_pca);

%% save time-consuming step, ~30min on linux
datadir = GetCurrentDataDir;
save(fullfile(datadir,'PCA_GLM_F8_Auto7.mat'),'C_betas','C_rsq','coeff');