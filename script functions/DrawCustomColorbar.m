function DrawCustomColorbar(clrmap,bins,numTicks,h)
%%
if exist('h','var')
    if ishandle(h) && strcmp(get(h,'type'),'figure')
        % this is customized for default figure from 'DrawCellsOnAnat'
        axes('Position',[0.8,0.8,0.05,0.15],'Units','normalized');
        axis off;
    elseif ishandle(h) && strcmp(get(h,'type'),'axes')
        axes(h);
        axis off;
    end     
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
% box off
% set(gca,'XColor','w');
% set(gca,'YColor','w');

set(gca,'YDir','normal');

if numTicks == 2
    set(gca,'YTick',[1,n]);
    set(gca,'YTickLabels',{num2str(bins(1),1),num2str(bins(end),1)});
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