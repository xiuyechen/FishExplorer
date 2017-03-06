function f = combineFiguresLR(hLeft, hRight)
% 2 usages: 
% 1) 2 input (hLeft, hRight)
% 2) input array of figure handles as first parameter ('hRight' not used)

if length(hLeft)==1
    f = combine2FiguresLR(hLeft,hRight);
else % multiple inputs, combine recursively
    M = hLeft;
    f = M(1);
    for i = 2:length(M)
        hRight = M(i);
        f = combine2FiguresLR(f,hRight);
    end
end
end

function f = combine2FiguresLR(hLeft, hRight)        
% get frame:
figure(hLeft) 
pos1 = get(hLeft,'Position');
cdata1 = print('-RGBImage','-r0');

figure(hRight)
pos2 = get(hRight,'Position');
cdata2 = print('-RGBImage','-r0');

% gaa edited to get size of 2nd dimension
w = pos1(3) + pos2(3);
h = max([pos1(4) pos2(4)]);
figPos = [50 50 w h];
r = pos1(3)/w;

% open new figure,  put in subplot:
f = figure();
set(f, 'Position', figPos)
% set(f, 'Color', [0 0 0]);
set(f, 'InvertHardCopy', 'off');
set(f, 'PaperPositionMode', 'auto');
subplot('Position',[0  0  r 1]), imagesc(cdata1)
axis('equal'), axis('tight'), axis('off')
subplot('Position',[r 0 (1-r) 1]), imagesc(cdata2)
axis('equal'), axis('tight'), axis('off')
drawnow();

% close singleton figures
close(hLeft);
close(hRight);
end