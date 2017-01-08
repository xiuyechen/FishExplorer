cIX = getappdata(hfig,'cIX');
% gIX = (1:length(cIX))';%getappdata(hfig,'gIX');

[gIX, numU] = SqueezeGroupIX(gIX);

fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
i_fish = getappdata(hfig,'i_fish');

C = getappdata(hfig,'M');
gIX = (1:size(C,1))';

% C = FindCentroid(hfig);
nClus = size(C,1);

setappdata(hfig,'gIX_betas',gIX);

[~,~,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
[~,~,regressor_m] = GetMotorRegressor(behavior,i_fish);

regs = vertcat(regressor_s,regressor_m);
orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?

betas = zeros(nClus,size(orthonormal_basis,2)+1);
X = [ones(size(orthonormal_basis,1),1),orthonormal_basis];
for i_clus = 1:nClus,
    y = C(i_clus,:)';
    betas(i_clus,:) = regress(y,X)';
end
%         betas = C * orthonormal_basis; % to reconstitute: betas(i,:)*orthonormal_basis'
% get ranking score: combined of motor coeffs
%         H = -sqrt(sum((betas(:,1:end-3)).^2,2)); % least stim-dependent
H = sqrt(sum((betas(:,end-2:end)).^2,2));
%%
% [gIX,rankscore] = SortH(H,gIX,numU,'descend');
% 
% function [gIX,B] = SortH(H,gIX,numU,descend) %#ok<INUSD> % new gIX is sorted based on H, size(H)=[numU,1];
descend = true;
if exist('descend','var'),
    [B,I] = sort(H,'descend');
else
    [B,I] = sort(H);
end
gIX_last = gIX;
for i = 1:numU,
    gIX(gIX_last==I(i)) = i;
end
% end
% setappdata(hfig,'gIX',gIX);
%% 1D colormap to visualize H-score ...
thres_multimotor = 0;
IX = find(H>thres_multimotor);
cIX_plot = cIX(IX);
H_select = H(IX);
numC = 64;
% clrmap = hot(numC);
clrmap = makeColormap([0 0 0],[1 0 0],numC);
gIX_plot = GetColorIndexFromScore(H_select,numC);

figure('Position',[600,50,600,900],'color',[1 1 1],...
    'Name',['Fish#' num2str(i_fish)]);
axes('Position',[0.03, 0.03, 0.94, 0.94]); % right ~subplot
% right subplot
I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap);
% I = LoadCurrentFishForAnatPlot(hfig);
DrawCellsOnAnat(I);

colormap(clrmap);
caxis([thres_multimotor,max(H)]);
% caxis([min(H),max(H)]);
colorbar('Location','manual','Position',[0.9,0.8,0.05,0.15])
