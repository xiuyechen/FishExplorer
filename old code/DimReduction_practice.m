%%
hfig = gcf;
M = getappdata(hfig,'M');
W_nnmf = getappdata(hfig,'W_nnmf');
H_nnmf = getappdata(hfig,'H_nnmf');
%%
tic
[coeff,score,latent,tsquared,explained,mu] = pca(M);
toc
beep



%% PCA face exercise

%% plot first faces of dataset
figure;
n = 3;
for i = 1:n^2,
    subplot(n,n,i);
    imagesc(reshape(fea(i,:),[faceH,faceW]));colormap gray;axis image
end

%% PCA
[coeff,score,latent,tsquared,explained]=pca(fea);

%% visualize first 9 bases
figure
for i = 1:9,
    subplot(3,3,i);
    imagesc(reshape(coeff(:,i),[faceH,faceW]));colormap gray;axis image
end

%% reconstruct! with pc's for 90% variance
temp = cumsum(explained);
k = find(temp>90,1,'first')

mu = mean(mean(fea));
score_ = score';
im = (coeff(:,1:k)*score_(1:k,:))'+mu;
figure;
for i = 1:9,
    subplot(3,3,i);
    imagesc(reshape(im(i,:),[faceH,faceW]),[0,255]);colormap gray;axis image
end

%% NMF
k = 43;
[W,H]=nnmf(fea,k);
%% visualize first 9 bases
figure;
for i = 1:9,
    subplot(3,3,i);
    imagesc(reshape(H(i,:),[faceH,faceW]));colormap gray;axis image
end

%% auto-contrast?
figure;
for i = 1:9,
    subplot(3,3,i);
    imagesc(reshape(H(i,:),[faceH,faceW]));colormap gray;axis image
end

%% reconstruct!
im = W*H;
figure;
for i = 1:9,
    subplot(3,3,i);
    imagesc(reshape(im(i,:),[faceH,faceW]));colormap gray;axis image
%     imagesc(reshape(im(i,:),[faceH,faceW]),[0,255]);colormap gray;axis image
end

%%
j = 3;
im = H(j,:)'*W(:,j)';
figure;
for i = 1:9,
    subplot(3,3,i);
    imagesc(reshape(im(:,i),[faceH,faceW]));colormap gray;axis image
end
%%
j = 2;

figure;
for i = 1:9,
    subplot(3,3,i);
    im = H(j,:)'*W(:,j)';
    imagesc(reshape(im(:,i),[faceH,faceW]));colormap gray;axis image
end



