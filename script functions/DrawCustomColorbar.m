function DrawCustomColorbar(clrmap,bins,numTicks,ax)
%%
if exist('ax','var')
    axes(ax);
else
    figure('Position',[500,200,10,300]);
    set(gcf,'color','w');
    set(gca,'Position',[0.4,0.1,0.2,0.8])
end
n = length(clrmap);
im = zeros(1,n,3);
im(1,:,:) = clrmap;
im2 = permute(im,[2 1 3]);
image(im2);

set(gca,'YDir','normal');

if numTicks == 2
    set(gca,'YTick',[1,n]);
    set(gca,'YTickLabels',[bins(1),bins(end)]);
    set(gca,'TickLength',[0,0]);
elseif numTicks == 3
        set(gca,'YTick',[1,n/2,n]);
    set(gca,'YTickLabels',[bins(1),mean([bins(1),bins(end)]),bins(end)]);
    set(gca,'TickLength',[0,0]);
    
else
set(gca,'YTick',1:length(bins));
set(gca,'YTickLabels',bins);
end
set(gca,'XTick',[]);

end