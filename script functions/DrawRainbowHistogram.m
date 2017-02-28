function DrawRainbowHistogram(N,bins,cmap)
figure; hold on
if ~exist('cmap') || size(cmap,1)~=length(N)
    cmap = jet(length(N));
end
for i = 1:length(N)
    x1 = bins(i);
    x2 = bins(i+1);
    y = N(i);    
    patch([x1,x2,x2,x1],[y,y,0,0],cmap(i,:));
end
end