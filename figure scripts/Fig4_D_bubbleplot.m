 
 %%
stimcorr = max(betas(:,1:end-3),[],2);
motorcorr = max(betas(:,end-2:end),[],2);
% figure;scatter(stimcorr,motorcorr)%scatter(motorcorr,stimcorr)

%%
% multi-motor

% make colormap
res = 100;
grad = linspace(0,1,res);
rev_grad = linspace(1,0,res);

grid = ones(res,res,3);

grid(:,:,3) = repmat(grad,res,1)';
grid(:,:,1) = 0.5*repmat(rev_grad',1,res)'+0.5*repmat(rev_grad,res,1)';
grid(:,:,2) = repmat(grad',1,res)';

clrmap_2D = reshape(grid,res*res,3);

figure;imagesc(grid)
axis xy
axis off
axis equal
%%
% get new gIX to use with this custom colormap
gIX_x = round((stimcorr-min(stimcorr))/(max(stimcorr)-min(stimcorr))*(res-1))+1;
gIX_y = round((motorcorr-min(motorcorr))/(max(motorcorr)-min(motorcorr))*(res-1))+1;

gIX_old = gIX;
U = unique(gIX_old);
% % gIX = zeros(size(gIX_old));
% % for i = 1:length(U);
% %     ix = gIX_old == U(i);
% %     gIX(ix) = sub2ind([res,res],gIX_x(U(i))',gIX_y(U(i))');
% % end

%%
U_size = zeros(size(gIX_x));
cmap_U = zeros(length(gIX_x),3);
for i = 1:length(U);
    ix = find(gIX_old == U(i));
    U_size(i) = length(ix);
    ix = sub2ind([res,res],gIX_y(U(i))',gIX_x(U(i))');
    cmap_U(i,:) = clrmap_2D(ix,:);
end
%%
% figure('Position',[500,500,250,200]);scatter(stimcorr,motorcorr,U_size)%,cmap)
% xlabel('stimulus corr.');ylabel('motor corr.');
%%
figure('Position',[500,500,250,200]);scatter(stimcorr,motorcorr,U_size,cmap_U)
xlabel('stimulus corr.');ylabel('motor corr.');

%%
radius = sqrt(stimcorr.^2 + motorcorr.^2);
[A,U_sorted] = sort(radius,'descend');
i_end = find(A>0.5,1,'last')
cIX_plot = []; gIX_plot = [];
cIX = getappdata(hfig,'cIX');
gIX = getappdata(hfig,'gIX');
for i = 1:i_end,
    IX = find(gIX == U_sorted(i));
    cIX_plot = [cIX_plot; cIX(IX)];
    gIX_plot = [gIX_plot; gIX(IX)] ;   
end