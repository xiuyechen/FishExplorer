figure('Position',[100 0 1300 900]);hold on;

% subplot(1,2,1)
cIX = find(temp==0);
% gIX = round(CellXYZ_norm(cIX,3)/max(CellXYZ_norm(cIX,3))*10)+1;
gIX = CellXYZ(cIX,3);
numK = max(gIX);
% cIX = cIX(CellXYZ_ref(cIX,3)>=68 & CellXYZ_ref(cIX,3)<=72);
% cIX = cIX(CellXYZ_ref(cIX,3)==70);
% gIX = ones(size(cIX));
% numK = 1;
% BasicDrawCellsOnAnatProj(data.CellXYZ,cIX,gIX,numK,data.anat_yx,data.anat_yz);
BasicDrawCellsOnAnatProj(CellXYZ_norm,cIX,gIX,numK,anat_yx,anat_yz);

%%

% OTec_left./OTec_total:   0.5991    0.4163    0.2874    0.4619    0.5020    0.4518
% AH_left./AH_total:          0.7892    0.1315    0.1610    0.2826    0.3827    0.1498
% abs((AH_left./AH_total-0.5)./(OTec_left./OTec_total - 0.5)): 2.9194    4.4049    1.5945    5.7093   58.8847    7.2593


 m1 = [0.5991,0.4163,0.2874,0.4619,0.5020,0.4518];
 m2 = [0.7892,0.1315,0.1610,0.2826,0.3827,0.1498];
 m3 = [2.9194,4.4049,1.5945,5.7093,58.8847,7.2593];
 y = vertcat(m1,m2)';
 y2 = vertcat(abs(m1-0.5)+0.5,abs(m2-0.5)+0.5)';
 
 figure;
 bar(y,'BaseValue',0.5);xlim([0,length(m1)+1]);ylim([0,1]);
 ylabel('left/total')
 xlabel('fish #')
 legend('upstream','downstream')
%  bar(y2,'BaseValue',0.5);ylim([0.5,1])
 
 %%
stimcorr = max(betas(:,1:end-3),[],2);
motorcorr = max(betas(:,end-2:end),[],2);
figure;scatter(stimcorr,motorcorr)%scatter(motorcorr,stimcorr)

%%
% multi-motor

% make colormap
res = 100;
grad = linspace(0,1,res);
rev_grad = linspace(1,0,res);

grid = ones(res,res,3);
grid(:,:,3) = repmat(grad,res,1);
grid(:,:,1) = 0.5*repmat(rev_grad',1,res)+0.5*repmat(rev_grad,res,1);
grid(:,:,2) = repmat(grad',1,res);

figure;imagesc(grid)
clrmap_2D = reshape(grid,res*res,3);
axis xy
%%
% get new gIX to use with this custom colormap
gIX_x = round((stimcorr-min(stimcorr))/(max(stimcorr)-min(stimcorr))*(res-1))+1;
gIX_y = round((motorcorr-min(motorcorr))/(max(motorcorr)-min(motorcorr))*(res-1))+1;

gIX_old = gIX;
U = unique(gIX_old);
gIX = zeros(size(gIX_old));
for i = 1:length(U);
    ix = gIX_old == U(i);
    gIX(ix) = sub2ind([res,res],gIX_y(U(i))',gIX_x(U(i))');
end
%%
U_size = zeros(size(gIX_x));
cmap = zeros(length(gIX_x),3);
for i = 1:length(U);
    ix = find(gIX_old == U(i));
    U_size(i) = length(ix);
    ix = sub2ind([res,res],gIX_x(U(i))',gIX_y(U(i))');
    cmap(i,:) = clrmap_2D(ix,:);
end
%%
figure('Position',[500,500,250,200]);scatter(stimcorr,motorcorr,U_size)%,cmap)
xlabel('stimulus corr.');ylabel('motor corr.');
