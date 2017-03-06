function f = combineFiguresTB(hTop, hBottom)

% get frame:
figure(hTop) 
pos1 = get(hTop,'Position');
cdata1 = print('-RGBImage','-r0');

figure(hBottom)
pos2 = get(hBottom,'Position');
cdata2 = print('-RGBImage','-r0');

% gaa edited to get size of 2nd dimension
h = pos1(4) + pos2(4);
w = max([pos1(3) pos2(3)]);
figPos = [50 50 w h];
r = pos1(4)/h;

% open new figure,  put in subplot:
f = figure();
set(f, 'Position', figPos);
set(f, 'Color', [0 0 0]);
set(f, 'InvertHardCopy', 'off');
set(f, 'PaperPositionMode', 'auto');
subplot('Position',[0 (1-r) 1 r]), imagesc(cdata1)
axis('equal'), axis('tight'), axis('off')
subplot('Position',[0 0 1 (1-r)]), imagesc(cdata2)
axis('equal'), axis('tight'), axis('off')
drawnow();

% close singleton figures
close(hTop);
close(hBottom);

