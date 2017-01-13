function ResetDisplayParams(hfig,i_fish)
% for general display
setappdata(hfig,'isPopout',1);
setappdata(hfig,'isPlotAnatomyOnly',0);
setappdata(hfig,'clrmap_name','hsv_new');

% for functional display: need to UpdateTimeIndex if changed
setappdata(hfig,'isStimAvr',0); % show average/full stimulus
setappdata(hfig,'isRawtime',0); % show stimulus in original order or sorted
setappdata(hfig,'isZscore',1); % show normalized (z-score) version of fluorescent

% for functional display: display only
setappdata(hfig,'isCentroid',0);
setappdata(hfig,'isPlotLines',1);
setappdata(hfig,'isPlotBehavior',1);
setappdata(hfig,'isPlotRegWithTS',0);
setappdata(hfig,'rankscore',[]); % 
setappdata(hfig,'rankID',0); % whether to plot scores or not, toggled away

% for anat display:
setappdata(hfig,'isRefAnat',1);
setappdata(hfig,'isShowFishOutline',1);
setappdata(hfig,'isShowMasks',1);
setappdata(hfig,'isShowMskOutline',0);
setappdata(hfig,'isWeighAlpha',0);

%% for specific fish
if exist('i_fish','var'),
    [~,stimrange] = GetStimRange([],i_fish);
    setappdata(hfig,'stimrange',stimrange); % need to UpdateTimeIndex if changed
    UpdateTimeIndex(hfig);
end
end