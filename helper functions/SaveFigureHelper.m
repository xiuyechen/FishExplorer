function SaveFigureHelper(saveFigFlag, saveDir, figName, figHandle)
 
if saveFigFlag
    if nargin>3
        figure(figHandle);
    end
    set(gcf, 'PaperPositionMode', 'auto');
    if ~exist(saveDir, 'dir'), mkdir(saveDir), end;
    fn = fullfile(saveDir, figName);    
    saveas(gcf, fn, 'png');    
    disp('figure saved')
    close(gcf)
end