function hclrbar = Add2DcolorbarToAnat(grid,cmins,cmaxs)
axes('Position',[0.8,0.8,0.05,0.15],'Units','normalized');
%%
res = size(grid,1);
% figure
imagesc(grid)
axis equal; axis tight
set(gca,'YDir','normal')
set(gca,'XTick',[1,res])
set(gca,'YTick',[1,res])

if exist('cmins','var') && exist('cmaxs','var')
    assert(length(cmins)==2);
    assert(length(cmaxs)==2);
    
    set(gca,'XTickLabels',[cmins(1),cmaxs(1)])
    set(gca,'YTickLabels',[cmins(2),cmaxs(2)])
end

if nargout>0
    hclrbar = gca;
end
end
