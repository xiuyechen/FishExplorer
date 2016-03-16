% f.pushbutton_autoclus_Callback

global hm1;
hObject = hm1;

save_masterdir = GetCurrentDataDir();
%%
range_fish = 1:11;
global hclusmenu;
for i_fish = range_fish,
    i_fish
    i_ClusGroup = 2;
    clusID = 2;
    Cluster = VAR(i_fish).ClusGroup{i_ClusGroup};
    
    %%
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    %%
    filename = fullfile(save_dir,'TimeSeries.h5');
    tic
    CellResp = h5read(filename,'/CellResp');
    CellRespZ = h5read(filename,'/CellRespZ');
    CellRespAvr = h5read(filename,'/CellRespAvr');
    CellRespAvrZ = h5read(filename,'/CellRespAvrZ');
    absIX = h5read(filename,'/absIX');
    
    setappdata(hfig,'CellResp',CellResp);
    setappdata(hfig,'CellRespZ',CellRespZ);
    setappdata(hfig,'CellRespAvr',CellRespAvr);
    setappdata(hfig,'CellRespAvrZ',CellRespAvrZ);
    setappdata(hfig,'absIX',absIX);
    toc
    %%
    
    numK = Cluster(clusID).numK;
    gIX = Cluster(clusID).gIX;
    
    % convert absolute index to index used for this dataset
    % absIX = getappdata(hfig,'absIX');
    
    cIX_abs = Cluster(clusID).cIX_abs;
    [~,cIX] = ismember(cIX_abs,absIX);
    
    % if ~isempty(find(cIX==0,1)),
    %     errordlg('cell index out of bound for currently loaded dataset');
    %     IX = cIX==0;
    %     cIX(IX) = [];
    %     gIX(IX) = [];
    % end
    f.UpdateIndices(hfig,cIX,gIX,numK);
    
    %%
    tic
    pushbutton_autoclus_Callback(hObject,f,i_fish);
    toc
end