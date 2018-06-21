addpath(genpath('C:\Users\Xiu\Dropbox\t-SNE'));

tic;[mappedA, mapping] = compute_mapping(M,'tSNE');toc

%%
            % Compute t-SNE mapping
% 			if isempty(varargin), mappedA = tsne(A, [], no_dims);
%             else mappedA = tsne(A, [], no_dims, varargin{1}); end
%             mapping.name = 't-SNE';

mappedA = tsne(M, [], 2); % no_dims = 2
%%
figure;
n = round(numK*1.1);
cmap = hsv(max(1,n));
gscatter(mappedA(:,1),mappedA(:,2),gIX,cmap,'.',[],'off');

%%
figure;hold on;%image([1,2],[3,4],pic)

k_scale = 800;
for i = 1:size(M,1),
    name = Astar.itemName{i};
    pic = imread([fullfile(picdir,name),'.jpg']);    
    dx = size(pic,2)/k_scale;
    dy = size(pic,1)/k_scale;
    image([mappedA(i,1),mappedA(i,1)+dx],[mappedA(i,2),mappedA(i,2)+dy],pic);
end
axis ij
axis equal

%%
save('C:\Janelia2014\tSNE_F10_autorightmotor.mat','mappedA','mapping','M','cIX','gIX','tIX','numK','stim','fictive');

%%
figure('Position',[50,100,1200,1200])
numK = 14;
rng('default');
[gIX,C,~,D] = kmeans(mappedA,numK,'Replicates',5);
Dist = min(D,[],2);

%
dthres = 80;
I_clean = [];
gIX_clean = [];
thres_size = 10;
for i = 1:numK,
    IX = find(gIX == i);
    dst = Dist(IX);
%     if mean(dst) < dthres,
%         I_clean_s = [I_clean_s; IX_s];
%         gIX_clean = [gIX_clean; gIX(IX)+double(numU)];
%     else
        ix = find(dst<dthres); % clean
%         if length(ix)>=thres_size, % clean cluster still big enough
            I = IX(ix);
            I_clean = [I_clean; I];
            gIX_clean = [gIX_clean; gIX(I)];
%         end
%     end

end
%
linewidth = 5;
subplot(2,2,1);
scatter(mappedA(:,1),mappedA(:,2),'.','MarkerEdgeColor',[0.5 0.5 0.5]);
xlim([-100,100]);ylim([-100,100]);
subplot(2,2,2)
n = round(numK*1.1);
cmap = hsv(max(1,n));
gscatter(mappedA(:,1),mappedA(:,2),gIX,cmap,'.',linewidth,'off');
xlim([-100,100]);ylim([-100,100]);
%
subplot(2,2,3)
n = round(numK*1.1);
cmap = hsv(max(1,n));
gscatter(mappedA(I_clean,1),mappedA(I_clean,2),gIX_clean,cmap,'.',linewidth,'off');
xlim([-100,100]);ylim([-100,100]);
%
subplot(2,2,4);
scatter(mappedA(:,1),mappedA(:,2),'.','MarkerEdgeColor',[0.5 0.5 0.5]);
hold on;
n = round(numK*1.1);
cmap = hsv(max(1,n));
gscatter(mappedA(I_clean,1),mappedA(I_clean,2),gIX_clean,cmap,'.',linewidth,'off');
xlim([-100,100]);ylim([-100,100]);

%%
figure
X = mappedA;

scatter(X(:,1),X(:,2),10,'ko')
%%
rng default;
numK = 14;
options = statset('Display','final');
gm = fitgmdist(X,numK,'Options',options);
%%
hold on
ezcontour(@(x,y)pdf(gm,[x y]));
hold off
%%
idx = cluster(gm,X);
Cluster = cell(1,numK);
for i = 1:numK,
    Cluster{i} = (idx == i);
end

n = round(numK*1.1);
cmap = hsv(max(1,n));

figure
hold on
for i = 1:numK,
scatter(X(Cluster{i},1),X(Cluster{i},2),10,cmap(i,:));
end
% hold off
% legend('Cluster 1','Cluster 2','Location','NW')
%%
thres = 0.3;
P = posterior(gm,X);
figure;
hold on
for i = 1:numK,
    ix = find(P(:,i)>thres);
scatter(X(ix,1),X(ix,2),10,P(ix,i),'o');
end
clrmap = jet(80); colormap(clrmap(9:72,:))
ylabel(colorbar,'Component 1 Posterior Probability')
xlim([-100,100]);ylim([-100,100]);
%%
P = posterior(gm,X);

scatter(X(cluster1,1),X(cluster1,2),10,P(cluster1,1),'+')
hold on
scatter(X(cluster2,1),X(cluster2,2),10,P(cluster2,1),'o')
hold off
legend('Cluster 1','Cluster 2','Location','NW')
clrmap = jet(80); colormap(clrmap(9:72,:))
ylabel(colorbar,'Component 1 Posterior Probability')

