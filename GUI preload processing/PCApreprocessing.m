% PCA on Fish data
tic
[coeff,score,latent,tsquared,explained,mu] = pca(M');
toc
code_dir = GetCurrentCodeDir();
save(fullfile(code_dir,'\saved mats\Fish3_full100_PCA.mat'),'coeff','score','latent','tsquared','explained','mu');

cIX_full = cIX;
gIX_full = gIX;

%%
cIX = cIX_full;
gIX = gIX_full;

%%
nPCs = 200;
figure;
im = coeff; % % rows reordered with optimal-leaf-order from dendrogram
% im = im.*repmat(explained',size(coeff,1),1); % PC's weighted by percent-explained
imagesc(im(:,1:nPCs));
% axis equal;
% set(gca,'YTick',1:size(coeff,1),'YTickLabel',ylabels);
xlim([0.5,nPCs+0.5])
xlabel('PC''s')

%% screen for cells that have high PC scores


coeff_crp = coeff(:,1:nPCs); % cropped
cellscore = abs(sum(coeff_crp,2));
[A,IX] = sort(cellscore,'descend');

percCell_PC = 0.1;
numCell = size(M,1);
nCell_PC = round(numCell*percCell_PC);
cIX = cIX(IX(1:nCell_PC));
gIX = ceil((1:length(cIX))'/length(cIX)*min(20,length(cIX)));
